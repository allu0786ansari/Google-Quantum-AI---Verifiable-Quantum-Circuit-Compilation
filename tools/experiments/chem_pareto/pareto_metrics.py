"""Shared metrics for the continuous chemistry tests (truncated tei, rank 786 at 1e-5)."""
import os, sys
import numpy as np
sys.path.insert(0, '/home/claude/frontier/chem_w')
os.environ.setdefault('VERIFIERS', '/home/claude/data/comp/verifiers')
import fastfit as F
import jax, jax.numpy as jnp
jax.config.update('jax_enable_x64', True)
GAP_FLOOR = F.GAP_FLOOR
P = F.Problem('/home/claude/work/lt/output/submission/dfthcsos/submission.npz')
N = P.N
eg = np.load('/home/claude/frontier/pareto/tei_eig.npy', allow_pickle=True).item()
keep = eg['ev'] > float(os.environ.get('TEI_TOL', '1e-5'))
Lr = jnp.asarray((eg['V'][:, keep] * np.sqrt(eg['ev'][keep])).T.reshape(-1, N, N))
S, A0, T2, const, Eup = jnp.asarray(P.S), jnp.asarray(P.A), P.T2, P.const, P.Eup


def start():
    return jnp.asarray(P.U), jnp.asarray(P.w0), jnp.asarray(P.W0)


def metrics(U, w, W):
    U = U / jnp.linalg.norm(U, axis=2, keepdims=True)
    G = jnp.einsum('rbp,sdp->rsbd', U, U) ** 2
    gram = jnp.einsum('rbc,rsbd,sde->rcse', w, G, w)
    A = jnp.einsum('rbp,kpq,rbq->rbk', U, Lr, U)
    sb = jnp.einsum('rbp,pq,rbq->rb', U, S, U)
    Mb = jnp.einsum('rbk,rdk->rbd', A, A) + 0.5 * (sb[:, :, None] + sb[:, None, :])
    e2sq = jnp.sum(gram ** 2) - 2 * jnp.einsum('rbc,rbd,rdc->', w, Mb, w) + T2
    coef = jnp.einsum('rc,rbc->rb', W, w)
    hp = A0 - jnp.einsum('rb,rbp,rbq->pq', coef, U, U)
    e = jnp.linalg.eigvalsh(hp)
    lh1 = jnp.sum(jnp.abs(e))
    t = W + jnp.sum(w, axis=1)
    lh2 = 0.25 * jnp.sum((jnp.sum(jnp.abs(w), axis=1) + jnp.abs(t)) ** 2)
    Lam = lh1 + lh2
    Esos = const - (lh1 - jnp.sum(e)) - 0.5 * jnp.sum(W ** 2)
    gap = Eup - Esos
    return jnp.sqrt(gap * (2 * Lam - gap)), jnp.sqrt(jnp.maximum(e2sq, 0)), gap, Lam
