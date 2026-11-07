"""Exact, fast refit of the continuous DF-THC weights (w, W) with the angle words frozen.

Frozen words => frozen vectors u_b => the ANGLE circuits are unchanged. Everything the
verifier scores then depends on (w, W) through small matrices:

  <X_k, X_l>        = w_k^T K_(rho rho') w_l,   K = (U_rho U_rho'^T) ** 2        (143 x 143)
  <X_k (x) X_k, T'> = w_k^T M_rho w_k,          M_bb' = (bb|b'b') + (s_b + s_b')/2
  error_2^2         = sum_kl (w_k^T K w_l)^2 - 2 sum_k w_k^T M w_k + ||T'||^2

with T'_pqrs = tei_pqrs + S_pq d_rs / 2 + d_pq S_rs / 2 (the verifier's shifted target),
s_b = u_b^T S u_b. h1_prime = A - sum_rho U_rho^T diag(sum_c W_rc w_rbc) U_rho.
error_1 uses the frozen one-body rows, matched to h1_prime's eigenvectors as the verifier does.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.environ['VERIFIERS'])
from dfthcsos.verify import load_inputs                       # noqa: E402
from dfthcsos.scoreinfra import dfthc_math as M               # noqa: E402

import jax                                                     # noqa: E402
jax.config.update('jax_enable_x64', True)
import jax.numpy as jnp                                        # noqa: E402

GAP_FLOOR = 3.914898


class Problem:
    def __init__(self, npz):
        inp = load_inputs(npz)
        z = dict(np.load(npz))
        self.z, self.inp = z, inp
        brot, edges, words = int(z['brot']), z['tree_edges'], z['angle_words']
        u = M.reconstruct_u_vectors(words / (1 << (brot - 1)), edges)
        u /= np.linalg.norm(u, axis=1, keepdims=True)
        w = np.asarray(inp.w, float)
        self.R, self.B, self.C = w.shape
        N = self.N = inp.num_orbitals
        self.one, ub = M.split_u_vectors(u, w)
        self.U = np.asarray(ub, float).reshape(self.R, self.B, N)
        ref = inp.ref
        self.eta = ref['electrons']
        S = self.S = np.asarray(inp.S, float)
        self.b1 = float(inp.beta_1)
        tei = np.asarray(ref['tei'], float)
        self.h1 = ref['oei'] - 0.5 * np.einsum('prrq->pq', tei)
        self.A = self.h1 - 0.5 * self.eta * S + self.b1 * np.eye(N)
        self.const = ref['ecore'] - self.eta * self.b1
        self.Eup = float(inp.energy_upper)
        # K blocks and M blocks
        G = np.einsum('rbp,sdp->rsbd', self.U, self.U)
        self.K = G ** 2                                        # (R, R, B, B)
        Mb = np.zeros((self.R, self.B, self.B))
        sb = np.einsum('rbp,pq,rbq->rb', self.U, S, self.U)
        for r in range(self.R):
            V = self.U[r]
            Y = np.einsum('bp,bq,pqrs->brs', V, V, tei, optimize=True)
            Mb[r] = np.einsum('brs,cr,cs->bc', Y, V, V, optimize=True) + 0.5 * (sb[r][:, None] + sb[r][None, :])
        self.M = Mb
        d = np.eye(N)
        T2 = (np.sum(tei ** 2) + 2 * 0.5 * np.einsum('pqrs,pq,rs->', tei, S, d) * 2 / 2
              + 2 * 0.5 * np.einsum('pqrs,pq,rs->', tei, d, S) * 2 / 2)
        # ||tei + S(x)d/2 + d(x)S/2||^2 expanded exactly
        a = 0.5
        T2 = (np.sum(tei ** 2)
              + 2 * a * np.einsum('pqrr,pq->', tei, S) + 2 * a * np.einsum('pprs,rs->', tei, S)
              + a * a * N * np.sum(S ** 2) + a * a * N * np.sum(S ** 2)
              + 2 * a * a * np.trace(S) ** 2)
        self.T2 = T2
        self.w0 = w
        self.W0 = np.asarray(inp.W, float)

    # ---- metrics ------------------------------------------------------------
    def jax_parts(self):
        K, Mb, U = jnp.asarray(self.K), jnp.asarray(self.M), jnp.asarray(self.U)
        A, T2, const, Eup = jnp.asarray(self.A), self.T2, self.const, self.Eup

        def parts(w, W):
            # w (R, B, C), W (R, C)
            gram = jnp.einsum('rbc,rsbd,sde->rcse', w, K, w)    # <X_rc, X_se>
            e2sq = jnp.sum(gram ** 2) - 2 * jnp.einsum('rbc,rbd,rdc->', w, Mb, w) + T2
            coef = jnp.einsum('rc,rbc->rb', W, w)
            hp = A - jnp.einsum('rb,rbp,rbq->pq', coef, U, U)
            e = jnp.linalg.eigvalsh(hp)
            lh1 = jnp.sum(jnp.abs(e))
            t = W + jnp.sum(w, axis=1)
            lh2 = 0.25 * jnp.sum((jnp.sum(jnp.abs(w), axis=1) + jnp.abs(t)) ** 2)
            Lam = lh1 + lh2
            Esos = const - (lh1 - jnp.sum(e)) - 0.5 * jnp.sum(W ** 2)
            gap = Eup - Esos
            lam = jnp.sqrt(gap * (2 * Lam - gap))
            return dict(e2=jnp.sqrt(jnp.maximum(e2sq, 0.0)), lam=lam, gap=gap, Lam=Lam)
        return parts

    def error_1(self, w, W):
        X = np.einsum('rbc,rbp,rbq->rcpq', w, self.U, self.U)
        return M.h1_thc_reconstruction_error(self.one, self.h1, self.S, self.b1, self.eta, W, X)

    def save(self, path, w, W):
        z = dict(self.z)
        z['W'] = w.astype(self.z['W'].dtype)
        z['H2_I'] = W.astype(self.z['H2_I'].dtype)
        np.savez(path, **z)


if __name__ == '__main__':
    P = Problem(sys.argv[1])
    parts = jax.jit(P.jax_parts())
    p = parts(jnp.asarray(P.w0), jnp.asarray(P.W0))
    print({k: float(v) for k, v in p.items()}, 'error_1', P.error_1(P.w0, P.W0))
