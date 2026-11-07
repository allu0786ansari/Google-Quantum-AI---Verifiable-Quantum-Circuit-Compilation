"""Exact integer model of layer I-2 (the arithmetic the circuit performs), over every x.

idx = top K bits of x_norm (22-bit), t = next M bits. Tables: C0[idx] (F0 fractional bits),
C1[idx] = |slope| * 2^(F1) (F1 fractional bits, < 2^W1). Circuit arithmetic:
    y0 = C0 - sum_j t_j * shift(C1, j + F0 - F1 - M)        (shift right truncates)
Coefficients: least-squares relative fit per interval, then C0 bias-corrected for the mean truncation.
"""
import numpy as np
from model_l1_l2 import xn, ytrue

K, M, F0, F1 = 7, 8, 18, 15


def tables(K=K, M=M, F0=F0, F1=F1):
    idx = xn >> (22 - K); t = (xn >> (22 - K - M)) & ((1 << M) - 1)
    C0 = np.zeros(1 << K, np.int64); C1 = np.zeros(1 << K, np.int64)
    for i in np.unique(idx):
        sel = idx == i; tf = t[sel] / float(1 << M)
        A = np.stack([np.ones(sel.sum()), -tf], 1)
        c0, c1 = np.linalg.lstsq(A / ytrue[sel, None], np.ones(sel.sum()), rcond=None)[0]
        C1[i] = int(round(c1 * (1 << F1)))
        C0[i] = int(round(c0 * (1 << F0)))
    return idx, t, C0, C1


def seed_int(idx, t, C0, C1, M=M, F0=F0, F1=F1):
    y = C0[idx].copy()
    for j in range(M):
        sh = j + F0 - F1 - M
        term = (C1[idx] << sh) if sh >= 0 else (C1[idx] >> (-sh))
        y -= ((t >> j) & 1) * term
    return y


idx, t, C0, C1 = tables()
y = seed_int(idx, t, C0, C1)
d0 = y / float(1 << F0) / ytrue - 1
# bias correction: shift C0 per entry so the max relative error is centred
for i in np.unique(idx):
    sel = idx == i
    lo, hi = d0[sel].min(), d0[sel].max()
    C0[i] -= int(round((lo + hi) / 2 * ytrue[sel].mean() * (1 << F0)))
y = seed_int(idx, t, C0, C1)
d0 = y / float(1 << F0) / ytrue - 1
d1 = -1.5 * d0 ** 2 - 0.5 * d0 ** 3
print(f'K={K} M={M} F0={F0} F1={F1}: entries {len(np.unique(idx))}, C1 bits {int(C1.max()).bit_length()}, '
      f'y0 bits {int(y.max()).bit_length()}')
print(f'seed worst relative error 2^{np.log2(np.abs(d0).max()):.2f}; after exact Newton 2^{np.log2(np.abs(d1).max()):.2f}')
np.savez('seed_tables.npz', C0=C0, C1=C1, K=K, M=M, F0=F0, F1=F1)
