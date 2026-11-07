"""Our initialisation of a compact DF-THC factorization (R groups x B shared vectors x C copies).

1. Eigendecompose the symmetry-shifted two-body tensor in the weighted pair basis
   (off-diagonal pairs weighted by sqrt 2, so truncation is optimal in the full four-index norm).
2. Keep the top R*C positive components as symmetric 76x76 matrices L_k = sqrt(s_k) M_k and split
   them into R groups of C by decreasing eigenvalue.
3. Per group, pick B unit vectors greedily from the eigenvectors of the group's matrices (largest
   |eigenvalue| first, skipping near-duplicates), then fit W by least squares for every copy.
Writes an .npz the starter trainer can be seeded from (same keys as the baselines).
    python df_init.py R B C out.npz
"""
import sys
import numpy as np
KIT = '/home/claude/data/starter/starter_kit'
sys.path.insert(0, KIT)
from dfthc.train_df_thc_sga import Initialize
from starter.task1_factorization.train_df_thc_qaggle import PACKAGED_FCIDUMP, FIXED_SHIFTS_FILE, PACKAGED_BASELINE
from starter.task1_factorization.solution_io import read_factorization

R, B, C, out = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
ref = Initialize(PACKAGED_FCIDUMP, -22140.27335)
N = ref['orbitals']
ti, tj = np.asarray(ref['tril'][0]), np.asarray(ref['tril'][1])
diag = np.asarray(ref['tril_diag'])
H2 = np.asarray(ref['H2_lower'], dtype=np.float64)
S = np.load(FIXED_SHIFTS_FILE)['Hsym_lower'].astype(np.float64)
T = H2.copy()
T[diag, :] += 0.5 * S[None, :]
T[:, diag] += 0.5 * S[:, None]
P = len(ti)
w = np.where(ti == tj, 1.0, np.sqrt(2.0))
Tw = (w[:, None] * T) * w[None, :]
s, V = np.linalg.eigh(Tw)
order = np.argsort(-s)
s, V = s[order], V[:, order]
K = R * C
print(f'pairs {P}; top eigenvalues {s[:4]}; kept {K}; positive {int((s > 0).sum())}; '
      f'tail norm {np.sqrt((s[K:] ** 2).sum()) * 1e3:.1f} mHa-ish', flush=True)


def to_mat(v):
    M = np.zeros((N, N)); M[ti, tj] = v; M[tj, ti] = v
    return M


Ls = [np.sqrt(max(s[k], 0)) * to_mat(V[:, k] / w) for k in range(K)]
U = np.zeros((R, B, N)); W = np.zeros((R, B, C))
for g in range(R):
    group = Ls[g * C:(g + 1) * C]
    cands = []
    for L in group:
        lam, vec = np.linalg.eigh(L)
        for a in range(N):
            cands.append((abs(lam[a]), vec[:, a]))
    cands.sort(key=lambda c: -c[0])
    chosen = []
    for val, v in cands:
        if all(abs(v @ u) < 0.98 for u in chosen):
            chosen.append(v)
        if len(chosen) == B:
            break
    Ug = np.array(chosen)
    U[g] = Ug
    A = np.stack([np.outer(u, u)[ti, tj] * w for u in Ug], axis=1)       # P x B, weighted
    for l, L in enumerate(group):
        W[g, :, l] = np.linalg.lstsq(A, L[ti, tj] * w, rcond=None)[0]
# report the initial two-body error in the trainer's norm
H2s = np.einsum('tjp,tjq,tjl->tlpq', U, U, W)[:, :, ti, tj]
rec = np.einsum('tlp,tlq->pq', H2s, H2s)
res = rec - T
d = res ** 2
err = np.sqrt(4 * d.sum() - 4 * d[diag, :].sum() + d[:, diag][diag, :].sum())
print(f'initial two-body error {err * 1e3:.1f} mHa', flush=True)
seed_params, _ = read_factorization(PACKAGED_BASELINE)
np.savez(out, U=U.astype(np.float32), W=W.astype(np.float32), H2_I=np.zeros((R, C), np.float32),
         Hsym_lower=np.asarray(seed_params['Hsym_lower']), H1_I=np.asarray(seed_params['H1_I']))
print('saved', out)
