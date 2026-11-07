"""Layer I-3, cheaper formulation: column-truncated partial products, band-limited residual, phase split.

  v = u * y0                 keep partial products of weight >= 2^-T1 (+ bias constant)
  w = v * y0  (mod 2^-10)    keep only partial products of weight in [2^-T2, 2^-11]; heavier ones are
                             multiples of 2^-10 and vanish mod 2^-10, and w = 1 - eps with |eps| < 2^-12
  E = signed window          -eps in units of 2^-T2 (T2 - 10 bits)
  c = y0 * eps / 2           keep partial products of weight >= 2^-T3 (+ bias)
  phase(y1) = phase(y0) + phase(c): y1 = y0 + c is never formed; c's bits get their own rotations.

Every kept partial product is one AND (free when one factor is a known-1 bit); a column-compression
(carry-save) tree turns N products into an n-bit result with ~N - n full adders, each 1 Toffoli.
Exact integer model over every x; reports worst relative error of y0 + c and the product counts.
"""
import sys
import numpy as np
from seed_model import seed_int, xn, ytrue

T = np.load('seed_tables.npz')
K, M, F0 = int(T['K']), int(T['M']), int(T['F0'])
idx = xn >> (22 - K); tt = (xn >> (22 - K - M)) & ((1 << M) - 1)
Y0 = seed_int(idx, tt, T['C0'], T['C1']).astype(np.int64)
XN = xn.astype(np.int64)
assert Y0.min() >= 1 << (F0 - 1) and Y0.max() < 1 << F0, 'y0 top bit must be a known 1'




def band_mul(A, na, a_exp, B, nb, b_exp, L, H=10**9, known1_b=()):
    """Sum of partial products A_i B_j 2^(a_exp+i+b_exp+j) with exponent in [L, H], row by row
    (row j = B_j times the slice of A whose products land in the band). Units 2^L.
    Returns (sum, #products, #ANDs)."""
    acc = np.zeros(len(A), np.int64); n = ands = 0
    for j in range(nb):
        ilo = max(0, L - a_exp - b_exp - j); ihi = min(na - 1, H - a_exp - b_exp - j)
        if ilo > ihi:
            continue
        sl = (A >> ilo) & ((1 << (ihi - ilo + 1)) - 1)
        acc += ((B >> j) & 1) * (sl << (a_exp + ilo + b_exp + j - L))
        k = ihi - ilo + 1; n += k; ands += 0 if j in known1_b else k
    return acc, n, ands


def dropped(na, a_exp, nb, b_exp, L):
    """Max total of the products below 2^L (all bits 1), in units 2^L."""
    return sum(2.0 ** (a_exp + i + b_exp + j - L) for i in range(na) for j in range(nb) if a_exp + i + b_exp + j < L)


def model(T1, T2, T3, bias=True, verbose=False):
    # A: v = u*y0 with u = XN 2^-20, y0 = Y0 2^-F0
    vA, nA, andA = band_mul(XN, 22, -20, Y0, F0, -F0, -T1, known1_b=(F0 - 1,))
    if bias:   # the dropped products average ~1/4 of their all-ones maximum: add it back as a constant
        vA = vA + int(round(dropped(22, -20, F0, -F0, -T1) / 4))
    nv = int(vA.max()).bit_length()
    # B: band [2^-T2, 2^-11] of v*y0, read mod 2^-10
    wB, nB, andB = band_mul(vA, nv, -T1, Y0, F0, -F0, -T2, -11, known1_b=(F0 - 1,))
    if bias:
        wB = wB + int(round(dropped(nv, -T1, F0, -F0, -T2) / 4))
    WEb = T2 - 10
    Wm = wB & ((1 << WEb) - 1)
    E = np.where(Wm >= 1 << (WEb - 1), Wm - (1 << WEb), Wm)      # E = -eps * 2^T2 (signed)
    assert np.abs(E).max() < 1 << (WEb - 2), 'window too narrow'
    # C: c = y0*eps/2 = -Y0 E 2^-(F0+T2+1); the top (sign) bit of E has weight -2^(WEb-1)
    Eu = E & ((1 << WEb) - 1)
    pos, nC1, andC1 = band_mul(Y0, F0, -F0, Eu & ((1 << (WEb - 1)) - 1), WEb - 1, -T2 - 1, -T3)
    neg, nC2, andC2 = band_mul(Y0, F0, -F0, Eu >> (WEb - 1), 1, -T2 - 1 + WEb - 1, -T3)
    acc = pos - neg
    nC, andC = nC1 + nC2, andC1 + andC2
    if bias:
        acc = acc + int(round(dropped(F0, -F0, WEb, -T2 - 1, -T3) / 4))
    Cv = -acc                                  # c in units 2^-T3
    y1 = Y0 / 2.0 ** F0 + Cv / 2.0 ** T3
    rel = np.abs(y1 / ytrue - 1).max()
    cbits = int(np.abs(Cv).max()).bit_length() + 1
    prods = nA + nB + nC; ands = andA + andB + andC
    tof = ands + (nA - nv) + (nB - WEb) + (nC - cbits)     # ANDs + full adders (forward)
    if verbose:
        print(f'   A: {nA} products ({andA} ANDs), v {nv} bits;  B: {nB} ({andB});  C: {nC} ({andC}); c {cbits} bits')
    return rel, prods, ands, tof, cbits


if __name__ == '__main__':
    print(f'seed: F0={F0}, Y0 in [2^{np.log2(Y0.min()):.2f}, 2^{np.log2(Y0.max()):.2f}]')
    print(' T1  T2  T3  bias   worst rel err   products  ANDs  est Toffolis (fwd)  c bits')
    best = []
    for T1 in range(25, 31):
        for T2 in range(25, 30):
            for T3 in range(26, 30):
                for b in (True,):
                    r, p, a, t, cb = model(T1, T2, T3, b)
                    best.append((t, T1, T2, T3, r, p, a, cb))
                    if np.log2(r) < -25.0:
                        print(f' {T1}  {T2}  {T3}  {b!s:5}   2^{np.log2(r):7.2f}     {p:5d}   {a:5d}   {t:6d}       {cb}', flush=True)
    ok = sorted(x for x in best if np.log2(x[4]) < -25.3)
    print('\ncheapest passing 2^-25.3:', ok[:5])
