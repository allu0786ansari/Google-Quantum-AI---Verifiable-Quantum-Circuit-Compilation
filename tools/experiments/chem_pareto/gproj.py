"""Gradient projection with restoration: minimize lambda_eff over (U, w, W), keeping
error_2 = target and E_gap = floor + margin (both active at the frontier), continuous U.

    python gproj.py ITERS [WHICH]   WHICH in {all, ww} (ww: only weights move, as a control)
"""
import os, sys, time
import numpy as np
sys.path.insert(0, '/home/claude/frontier/chem_w'); sys.path.insert(0, '/home/claude/frontier/pareto')
os.environ.setdefault('VERIFIERS', '/home/claude/data/comp/verifiers')
import jax, jax.numpy as jnp
jax.config.update('jax_enable_x64', True)
import importlib.util
spec = importlib.util.spec_from_file_location('pm', '/home/claude/frontier/pareto/pareto_metrics.py')
pm = importlib.util.module_from_spec(spec); spec.loader.exec_module(pm)

iters = int(sys.argv[1]); which = sys.argv[2] if len(sys.argv) > 2 else 'all'
U0, w0, W0 = pm.start()
shapes = [U0.shape, w0.shape, W0.shape]; sizes = [int(np.prod(s)) for s in shapes]
mask = np.concatenate([np.full(sizes[0], which == 'all'), np.ones(sizes[1] + sizes[2], bool)]).astype(float)


def unpack(x):
    a, b = sizes[0], sizes[0] + sizes[1]
    return x[:a].reshape(shapes[0]), x[a:b].reshape(shapes[1]), x[b:].reshape(shapes[2])


f_all = jax.jit(lambda x: jnp.stack(pm.metrics(*unpack(x))[:3]))           # lam, e2, gap
J_all = jax.jit(jax.jacrev(lambda x: jnp.stack(pm.metrics(*unpack(x))[:3])))
x = jnp.concatenate([U0.ravel(), w0.ravel(), W0.ravel()])
lam0, E2T, gap0 = [float(v) for v in f_all(x)]
GAPT = pm.GAP_FLOOR + 2e-6
print(f'start lambda {lam0:.6f} e2 {E2T*1e3:.5f} gap {gap0:.7f} (which={which})', flush=True)
mk = jnp.asarray(mask)


def restore(x, n=6):
    for _ in range(n):
        v = f_all(x); c = jnp.array([v[1] - E2T, v[2] - GAPT])
        if float(jnp.max(jnp.abs(c / jnp.array([1e-6, 1e-7])))) < 1:
            return x, True
        J = J_all(x)[1:] * mk
        x = x - J.T @ jnp.linalg.solve(J @ J.T, c)
    v = f_all(x)
    return x, bool(abs(float(v[1]) - E2T) < 1e-6 and float(v[2]) > pm.GAP_FLOOR)


x, ok = restore(x)
alpha = 1e-4; t0 = time.time(); lam = float(f_all(x)[0])
for it in range(1, iters + 1):
    J = J_all(x) * mk
    gl, C = J[0], J[1:]
    d = -(gl - C.T @ jnp.linalg.solve(C @ C.T, C @ gl))
    dn = float(jnp.linalg.norm(d))
    accepted = False
    for _ in range(12):
        xn, ok = restore(x + alpha * d / dn)
        ln = float(f_all(xn)[0])
        if ok and ln < lam - 1e-9:
            x, lam, accepted = xn, ln, True
            alpha *= 1.6
            break
        alpha *= 0.3
    if not accepted:
        print(f'iter {it}: no descent (alpha {alpha:.2e}); stop', flush=True)
        break
    if it <= 3 or it % 10 == 0:
        v = f_all(x)
        print(f'iter {it}: lambda {lam:.6f} ({100*(lam0-lam)/lam0:.3f}% lower) e2 {float(v[1])*1e3:.5f} '
              f'gap {float(v[2]):.7f} step {alpha:.2e} |d| {dn:.3g} ({time.time()-t0:.0f}s)', flush=True)
np.save(f'gproj_{which}.npy', np.asarray(x))
print(f'END lambda {lam:.6f} vs start {lam0:.6f}: {100*(lam0-lam)/lam0:.3f}% lower')
