"""Probe: are the frontier's 8-bit two-body angle words locally optimal on the grid?

For random single-word moves (vector b, tree angle t, +-1 step) compute the exact change of
lambda_eff and error_2 at fixed (w, W), and the Lagrangian change dL = d lambda + nu * d error_2[mHa]
with nu = 0.36 lambda/mHa (the KKT multiplier measured at the frontier).

    python gridmoves.py [N_MOVES]
"""
import os
import sys
import time

import numpy as np

import fastfit as F
from dfthcsos.scoreinfra import dfthc_math as M

NU = 0.36
P = F.Problem('/home/claude/work/lt/output/submission/dfthcsos/submission.npz')
z = P.z
brot, edges, words = int(z['brot']), z['tree_edges'], z['angle_words']
N, R, B, C = P.N, P.R, P.B, P.C
two = words[N:].copy()                                   # (R*B, N-1)
tei = np.asarray(P.inp.ref['tei'], float)
t0 = time.time()
ev, V = np.linalg.eigh(tei.reshape(N * N, N * N))
keep = ev > ev.max() * 1e-13
Lf = V[:, keep] * np.sqrt(ev[keep])                      # tei = Lf Lf^T  (N^2, r)
Lr = Lf.T.reshape(-1, N, N)                              # (r, N, N)
print(f'tei factor rank {keep.sum()} ({time.time() - t0:.0f}s)', flush=True)
S = P.S
U = P.U.reshape(R * B, N)
w, W = P.w0, P.W0


def A_of(u):                                             # (n, r): u^T L_r u
    return np.einsum('np,rpq,nq->nr', u, Lr, u, optimize=True)


A = A_of(U)
sb = np.einsum('np,pq,nq->n', U, S, U)


def metrics(U_, A_, sb_):
    Ug = U_.reshape(R, B, N)
    G = np.einsum('rbp,sdp->rsbd', Ug, Ug) ** 2
    gram = np.einsum('rbc,rsbd,sde->rcse', w, G, w)
    Ag = A_.reshape(R, B, -1); sg = sb_.reshape(R, B)
    Mb = np.einsum('rbk,rdk->rbd', Ag, Ag) + 0.5 * (sg[:, :, None] + sg[:, None, :])
    e2sq = np.sum(gram ** 2) - 2 * np.einsum('rbc,rbd,rdc->', w, Mb, w) + P.T2
    coef = np.einsum('rc,rbc->rb', W, w)
    hp = P.A - np.einsum('rb,rbp,rbq->pq', coef, Ug, Ug)
    e = np.linalg.eigvalsh(hp)
    lh1 = np.abs(e).sum()
    t = W + w.sum(1)
    lh2 = 0.25 * np.sum((np.abs(w).sum(1) + np.abs(t)) ** 2)
    Lam = lh1 + lh2
    Esos = P.const - (lh1 - e.sum()) - 0.5 * np.sum(W ** 2)
    gap = P.Eup - Esos
    return np.sqrt(gap * (2 * Lam - gap)), np.sqrt(max(e2sq, 0)), gap


lam0, e20, gap0 = metrics(U, A, sb)
print(f'base: lambda {lam0:.6f} e2 {e20 * 1e3:.5f} mHa gap {gap0:.7f} ({time.time() - t0:.0f}s)', flush=True)
rng = np.random.default_rng(0)
n_moves = int(sys.argv[1]) if len(sys.argv) > 1 else 400
good = []
for it in range(n_moves):
    i, t, s = rng.integers(R * B), rng.integers(N - 1), rng.choice([-1, 1])
    row = two[i].copy(); row[t] = (row[t] + s) % (1 << brot)
    u = M.reconstruct_u_vectors(row[None] / (1 << (brot - 1)), edges)[0]
    u /= np.linalg.norm(u)
    U2 = U.copy(); U2[i] = u
    A2 = A.copy(); A2[i] = A_of(u[None])[0]
    sb2 = sb.copy(); sb2[i] = u @ S @ u
    lam, e2, gap = metrics(U2, A2, sb2)
    dL = (lam - lam0) + NU * (e2 - e20) * 1e3
    if gap > F.GAP_FLOOR and dL < 0:
        good.append((dL, i, t, s, lam - lam0, (e2 - e20) * 1e3))
dLs = sorted(good)
print(f'{len(good)} of {n_moves} random moves lower the Lagrangian; best 5:', [tuple(round(float(x), 5) for x in g) for g in dLs[:5]],
      f'({time.time() - t0:.0f}s)')
