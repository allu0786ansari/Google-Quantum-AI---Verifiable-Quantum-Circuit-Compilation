"""Augmented-Lagrangian refit of (w, W) at frozen angle words; optional copy removal.

    python refit.py IN_NPZ OUT_NPZ [--drop K] [--rounds 6] [--iters 400]

--drop K removes K copies per group (C -> C - K) before refitting, choosing in each group
the copies whose removal raises error_2 least (greedy, re-scored after each removal).
"""
import argparse
import time

import numpy as np
from scipy.optimize import minimize

import fastfit as F
import jax
import jax.numpy as jnp

ap = argparse.ArgumentParser()
ap.add_argument('inp'); ap.add_argument('out')
ap.add_argument('--drop', type=int, default=0)
ap.add_argument('--add', type=int, default=0)
ap.add_argument('--jitter', type=float, default=1e-3)
ap.add_argument('--seed', type=int, default=0)
ap.add_argument('--rounds', type=int, default=6)
ap.add_argument('--iters', type=int, default=400)
ap.add_argument('--e2max', type=float, default=0.019990)
ap.add_argument('--gapmargin', type=float, default=2e-6)
args = ap.parse_args()

t0 = time.time()
P = F.Problem(args.inp)
parts = P.jax_parts()
R, B, C = P.R, P.B, P.C
w, W = P.w0.copy(), P.W0.copy()

if args.drop:
    e2f = jax.jit(lambda w_, W_: parts(w_, W_)['e2'])
    for _ in range(args.drop):
        keep = np.ones((R, w.shape[2]), bool)
        drop_idx = []
        for r in range(R):
            best = None
            for c in range(w.shape[2]):
                w2 = w.copy(); w2[r, :, c] = 0
                e = float(e2f(jnp.asarray(w2), jnp.asarray(W)))
                if best is None or e < best[0]:
                    best = (e, c)
            drop_idx.append(best[1])
        w = np.stack([np.delete(w[r], drop_idx[r], axis=1) for r in range(R)])
        W = np.stack([np.delete(W[r], drop_idx[r]) for r in range(R)])
    print(f'dropped to C={w.shape[2]}: e2 {float(e2f(jnp.asarray(w), jnp.asarray(W))) * 1e3:.3f} mHa '
          f'({time.time() - t0:.0f}s)', flush=True)

if args.add:
    # exact split: copy c -> two copies scaled 1/sqrt(2) leaves X(x)X sums, sum W X, sum W^2 and the
    # one-norm unchanged; perturb slightly so the optimizer can separate them.
    rng = np.random.default_rng(args.seed)
    ws, Ws = [], []
    for r in range(R):
        strength = np.abs(W[r]) * np.abs(w[r]).sum(0)
        top = np.argsort(-strength)[:args.add]
        wr, Wr = w[r].copy(), W[r].copy()
        wr[:, top] /= np.sqrt(2); Wr[top] /= np.sqrt(2)
        new_w = wr[:, top] * (1 + args.jitter * rng.standard_normal(wr[:, top].shape))
        ws.append(np.concatenate([wr, new_w], axis=1)); Ws.append(np.concatenate([Wr, Wr[top]]))
    w, W = np.stack(ws), np.stack(Ws)
    print(f'split to C={w.shape[2]}', flush=True)

shape_w, n_w = w.shape, w.size


def unpack(x):
    return x[:n_w].reshape(shape_w), x[n_w:].reshape(W.shape)


floor = F.GAP_FLOOR + args.gapmargin
SC = 1e-3                                       # constraint scale: 1 mHa


def psi(c, nu, rho):
    return jnp.where(c >= -nu / rho, nu * c + 0.5 * rho * c * c, -nu * nu / (2 * rho))


def make_obj(nu1, nu2, rho):
    def L(x):
        w_, W_ = unpack(x)
        p = parts(w_, W_)
        c1 = (p['e2'] - args.e2max) / SC
        c2 = (floor - p['gap']) / SC
        return p['lam'] + psi(c1, nu1, rho) + psi(c2, nu2, rho)
    return jax.jit(jax.value_and_grad(L))


x = np.concatenate([w.ravel(), W.ravel()])
evalp = jax.jit(lambda x_: parts(*unpack(x_)))
p = evalp(jnp.asarray(x))
print('start', {k: round(float(v), 7) for k, v in p.items()}, flush=True)
nu1 = nu2 = 0.0
rho = float(__import__("os").environ.get("RHO0", "1000"))
best = None
for rd in range(args.rounds):
    vg = make_obj(nu1, nu2, rho)
    res = minimize(lambda z: tuple(np.asarray(v, float) for v in vg(jnp.asarray(z))), x, jac=True,
                   method='L-BFGS-B', options=dict(maxiter=args.iters, maxcor=30))
    x = res.x
    p = {k: float(v) for k, v in evalp(jnp.asarray(x)).items()}
    c1, c2 = (p['e2'] - args.e2max) / SC, (floor - p['gap']) / SC
    nu1, nu2 = max(0.0, nu1 + rho * c1), max(0.0, nu2 + rho * c2)
    rho = min(rho * 2, 1e6)
    feas = p['e2'] <= 0.019999 and p['gap'] > F.GAP_FLOOR
    print(f'round {rd}: lam {p["lam"]:.6f} e2 {p["e2"] * 1e3:.4f} gap {p["gap"]:.7f} Lam {p["Lam"]:.3f} '
          f'feasible {feas} nu ({nu1:.3g},{nu2:.3g}) ({time.time() - t0:.0f}s)', flush=True)
    if feas and (best is None or p['lam'] < best[0]):
        best = (p['lam'], x.copy())
if best is None:
    raise SystemExit('no feasible point')
w1, W1 = unpack(best[1])
e1 = P.error_1(w1, W1)
print(f'best lam {best[0]:.6f}, error_1 {e1 * 1e3:.3f} mHa, C={w1.shape[2]}')
P.save(args.out, w1, W1)
