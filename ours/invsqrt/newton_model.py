"""Exact integer model of layer I-3 (Newton step) as the circuit computes it, over every x.

Inputs: XN (22-bit, u = XN / 2^20 in [1, 4)), Y0 (18-bit, y0 = Y0 / 2^18).
  SQ   = sum_j Y0_j * floor(Y0 * 2^j / 2^D1)                    (truncated square, y0^2 ~ SQ / 2^(36-D1))
  BW   = sum_j XN_j * floor(SQ * 2^j / 2^D2)  mod 2^WE            (window of u*y0^2 below the known top bits)
  E    = signed(BW) ; u*y0^2 = 1 + E * 2^(D2 - 20 - FA) (mod the dropped top), so eps = -E * 2^(...)
  Y1   = Y0 * 2^(F1 - 18) - sum_k E_k * floor(Y0 * 2^k * 2^SH)   (y1 = y0 - y0*E*scale/2; signed top row adds)
Every floor() mirrors a controlled add of a sliced register in the circuit.
"""
import numpy as np
from seed_model import seed_int, xn, ytrue
T = np.load('seed_tables.npz')
K, M = int(T['K']), int(T['M'])
idx = xn >> (22 - K); tt = (xn >> (22 - K - M)) & ((1 << M) - 1)
Y0 = seed_int(idx, tt, T['C0'], T['C1']).astype(object)
XN = xn.astype(object)


def rows_mul(A, Bv, nb, shift_bits, width):
    """sum_j Bv_j * floor(A * 2^j / 2^shift_bits) mod 2^width (shift may be negative = left)."""
    acc = np.zeros(len(A), dtype=object)
    for j in range(nb):
        sh = shift_bits - j
        term = (A >> sh) if sh >= 0 else (A << (-sh))
        acc = acc + ((Bv >> j) & 1) * term
    return acc % (1 << width), nb


def newton(FA, WE, F1, g2=3):
    D1 = 36 - FA
    SQ, _ = rows_mul(Y0, Y0, 18, D1, FA + 1)                      # y0^2 < 1 so FA bits (+1 safety)
    # u*y0^2 * 2^(20+FA) = XN * SQ ; window bits [D2, D2 + WE) where WE covers |eps| < 2^-12.4 plus sign
    D2 = 20 + FA - (12 + WE - 2) - 1                              # window top sits just above |eps| bits
    D2 = 20 + FA - (WE - 1) - 11                                  # |E| < 2^(WE-2): eps up to 2^-12 needs top at 2^-11
    BW, _ = rows_mul(SQ, XN, 22, D2, WE)
    E = np.where(BW >= (1 << (WE - 1)), BW - (1 << WE), BW)       # signed
    # y1 = y0 * (1 - (u y0^2 - 1)/2) = y0 - y0 * E * 2^(D2 - 20 - FA) / 2
    SH = (D2 - 20 - FA - 1) - 18 + F1                             # Y0 * E * 2^SH in Y1 units (F1 frac bits)
    Y1 = Y0 * (1 << (F1 - 18))
    Eu = E % (1 << WE)
    corr = np.zeros(len(Y0), dtype=object)
    for k in range(WE):
        sh = -(SH + k)
        term = (Y0 >> sh) if sh >= 0 else (Y0 << (-sh))
        bit = (Eu >> k) & 1
        corr = corr + (bit * term if k < WE - 1 else -bit * term)  # top bit has weight -2^(WE-1)
    Y1 = Y1 - corr
    y1 = np.array([float(v) for v in Y1]) / 2.0 ** F1
    rel = y1 / ytrue - 1
    # rough cost: controlled-add Toffolis ~ 2.4 x target width per row (kickmix iadd)
    cost = 2.4 * (18 * (FA + 1) + 22 * WE + WE * (F1 - 16))
    return np.abs(rel).max(), cost, D2


print(' FA  WE  F1   worst rel err (log2)   est. Toffolis')
for FA in (28, 30, 32):
    for WE in (17, 19, 21):
        for F1 in (28, 30):
            e, c, D2 = newton(FA, WE, F1)
            print(f' {FA}  {WE}  {F1}      {np.log2(e):8.2f}            {c:6.0f}   (D2={D2})', flush=True)
