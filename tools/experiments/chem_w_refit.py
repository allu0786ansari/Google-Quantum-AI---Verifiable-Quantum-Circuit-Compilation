"""Optimize W (= H2_I) only, with the loaded words, w and the fixed shifts unchanged.

error_2 does not depend on W, so only lambda_eff, E_gap and error_1 move.
    python wopt.py IN_NPZ OUT_NPZ
"""
import sys, os
import numpy as np
from scipy.optimize import minimize
sys.path.insert(0, os.environ['VERIFIERS'])
from dfthcsos.verify import load_inputs
from dfthcsos.scoreinfra import dfthc_math as M

GAP_FLOOR, E1_MAX = 3.914898, 0.020
inp = load_inputs(sys.argv[1])
z = dict(np.load(sys.argv[1]))
brot, edges, words = int(z['brot']), z['tree_edges'], z['angle_words']
u = M.reconstruct_u_vectors(words / (1 << (brot - 1)), edges)
u /= np.linalg.norm(u, axis=1, keepdims=True)
w = np.asarray(inp.w, float); R, B, C = w.shape; N = inp.num_orbitals
one, ub = M.split_u_vectors(u, w)
X = M.reconstruct_x(ub, w); X = np.tril(X) + np.swapaxes(np.tril(X, k=-1), -1, -2)
ref = inp.ref; eta = ref['electrons']; S = np.asarray(inp.S, float); b1 = float(inp.beta_1)
h1 = ref['oei'] - 0.5 * np.einsum('prrq->pq', ref['tei'])
A = h1 - 0.5 * eta * S + b1 * np.eye(N)
Xk = X.reshape(R * C, N, N)
a = np.abs(w).sum(1).reshape(-1)          # ||w_(rho c)||_1, (R*C,)
s = w.sum(1).reshape(-1)                  # sum_b w_b
trX = np.einsum('kpp->k', Xk)
const = ref['ecore'] - eta * b1
Eup = float(inp.energy_upper)


def parts(Wf, grad=True):
    hp = A - np.einsum('k,kpq->pq', Wf, Xk)
    e, V = np.linalg.eigh(hp)
    lh1 = np.abs(e).sum()
    t = Wf + s
    lh2 = 0.25 * np.sum((a + np.abs(t)) ** 2)
    Lam = lh1 + lh2
    Esos = const - (lh1 - e.sum()) - 0.5 * np.sum(Wf ** 2)
    g = Eup - Esos
    lam = np.sqrt(g * (2 * Lam - g))
    out = dict(lam=lam, gap=g, Lam=Lam, lh1=lh1, lh2=lh2)
    if grad:
        G = np.einsum('pi,kpq,qi->ik', V, Xk, V)        # v_i^T X_k v_i
        dlh1 = -(np.sign(e) @ G)
        dsum = -trX
        dLam = dlh1 + 0.5 * (a + np.abs(t)) * np.sign(t)
        dEsos = -(dlh1 - dsum) - Wf
        dg = -dEsos
        out['dlam'] = (dg * (2 * Lam - g) + g * (2 * dLam - dg)) / (2 * lam)
        out['dgap'] = dg
    return out


def error_1(Wf):
    return M.h1_thc_reconstruction_error(one, h1, S, b1, eta, Wf.reshape(R, C), X)


W0 = np.asarray(inp.W, float).reshape(-1)
p0 = parts(W0)
print('start: lam %.6f gap %.7f Lam %.4f (h1 %.3f, h2 %.3f) error_1 %.3f mHa' %
      (p0['lam'], p0['gap'], p0['Lam'], p0['lh1'], p0['lh2'], 1e3 * error_1(W0)))
margin = float(os.environ.get('GAP_MARGIN', '2e-6'))
res = minimize(lambda x: (parts(x)['lam'], parts(x)['dlam']), W0, jac=True, method='SLSQP',
               constraints=[dict(type='ineq', fun=lambda x: parts(x, False)['gap'] - GAP_FLOOR - margin,
                                 jac=lambda x: parts(x)['dgap'])],
               options=dict(maxiter=int(os.environ.get('ITERS', '300')), ftol=1e-12))
W1 = res.x
# keep error_1 within budget: back off toward W0 if needed
for t in np.linspace(1, 0, 41):
    Wt = W0 + t * (W1 - W0)
    pt = parts(Wt, False)
    if error_1(Wt) <= E1_MAX * 0.995 and pt['gap'] > GAP_FLOOR:
        break
p1 = parts(Wt, False)
print(res.message, 'step', t)
print('end:   lam %.6f gap %.7f Lam %.4f (h1 %.3f, h2 %.3f) error_1 %.3f mHa' %
      (p1['lam'], p1['gap'], p1['Lam'], p1['lh1'], p1['lh2'], 1e3 * error_1(Wt)))
print('lambda_eff ratio %.5f' % (p0['lam'] / p1['lam']))
z['H2_I'] = Wt.reshape(R, C).astype(z['H2_I'].dtype)
np.savez(sys.argv[2], **z)
