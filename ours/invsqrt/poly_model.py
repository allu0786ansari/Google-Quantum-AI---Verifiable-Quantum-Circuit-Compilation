"""Layer I-2/I-3 alternative: no Newton step; a direct piecewise polynomial in fixed point.

  idx = top K bits of x_norm, z = remaining 22-K bits as a fraction in [0, 1)
  y = c0 - z*(c1 - z*(c2 - z*c3))      all c_k >= 0 (1/sqrt alternates in sign), Horner with
                                       column-truncated products: stage k keeps products of weight
                                       >= 2^-T and uses only the top NZ[k] bits of z
c0 is fitted last, per interval, on the exact integer pipeline: it centres the worst error, absorbing
every truncation bias. Exact over every x. Reports worst relative error, products, ANDs, word widths.
"""
import sys
import numpy as np
from model_l1_l2 import xn, ytrue
from newton_band import band_mul


def fit(K, deg, P, NZ, T):
    """P[k] = fractional bits of c_k (k >= 1); NZ[k] = z bits used at stage k; T = product cut."""
    nzb = 22 - K
    idx = xn >> nzb
    Z = (xn & ((1 << nzb) - 1)).astype(np.int64)        # z = Z / 2^nzb
    zf = Z / float(1 << nzb)
    ent = np.unique(idx)
    C = {k: np.zeros(1 << K, np.int64) for k in range(deg + 1)}
    for i in ent:                                        # relative least squares, signs alternate
        sel = idx == i
        A = np.vander(zf[sel], deg + 1, increasing=True)
        coef = np.linalg.lstsq(A / ytrue[sel, None], np.ones(sel.sum()), rcond=None)[0]
        for k in range(1, deg + 1):
            v = coef[k] * (-1) ** k
            assert v >= 0, (i, k, v)
            C[k][i] = int(round(v * (1 << P[k])))
    # Horner from the top: acc (units 2^-P_acc), nonnegative
    acc = C[deg][idx]; pacc = P[deg]
    prods = ands = 0
    for k in range(deg, 0, -1):
        nz = NZ[k]
        Zk = Z >> (nzb - nz)                              # top nz bits of z, units 2^-nz
        na = int(acc.max()).bit_length()
        pr, n, a = band_mul(acc, na, -pacc, Zk, nz, -nz, -T)
        prods += n; ands += a
        if k > 1:
            nxt = (C[k - 1][idx] << (T - P[k - 1])) if T >= P[k - 1] else (C[k - 1][idx] >> (P[k - 1] - T))
            acc = nxt - pr; pacc = T
            assert acc.min() >= 0
        else:
            prod1 = pr                                    # z*acc1 in units 2^-T
    rest = -prod1 / 2.0 ** T
    for i in ent:                                         # c0 centres the error, per interval
        sel = idx == i
        d = ytrue[sel] - rest[sel]
        lo, hi = (d / ytrue[sel]).min(), (d / ytrue[sel]).max()
        C[0][i] = int(round(np.median(d) * 2 ** T))
        # minimax on relative error: shift so the relative error extremes balance
        best = None
        for c in range(C[0][i] - 6, C[0][i] + 7):
            e = np.abs((c / 2.0 ** T + rest[sel]) / ytrue[sel] - 1).max()
            if best is None or e < best[0]:
                best = (e, c)
        C[0][i] = best[1]
    y = C[0][idx] / 2.0 ** T + rest
    rel = np.abs(y / ytrue - 1).max()
    widths = {k: int(C[k][ent].max()).bit_length() for k in range(deg + 1)}
    return rel, prods, ands, widths, len(ent)


if __name__ == '__main__':
    print(' K deg  P            NZ          T   worst      prods  ANDs  word bits      entries')
    for K, deg, P, NZ, T in [
        (7, 3, {1: 28, 2: 28, 3: 28}, {1: 15, 2: 15, 3: 15}, 28),
        (7, 3, {1: 28, 2: 28, 3: 28}, {1: 15, 2: 13, 3: 10}, 28),
        (7, 3, {1: 28, 2: 28, 3: 27}, {1: 15, 2: 12, 3: 8}, 28),
        (7, 3, {1: 27, 2: 27, 3: 26}, {1: 15, 2: 12, 3: 8}, 27),
        (6, 3, {1: 28, 2: 28, 3: 28}, {1: 16, 2: 14, 3: 11}, 28),
        (6, 4, {1: 28, 2: 28, 3: 28, 4: 28}, {1: 16, 2: 14, 3: 11, 4: 8}, 28),
        (5, 4, {1: 28, 2: 28, 3: 28, 4: 28}, {1: 17, 2: 15, 3: 12, 4: 9}, 28),
        (8, 2, {1: 28, 2: 28}, {1: 14, 2: 12}, 28),
        (9, 2, {1: 28, 2: 28}, {1: 13, 2: 10}, 28),
    ]:
        rel, p, a, w, ent = fit(K, deg, P, NZ, T)
        print(f' {K}  {deg}  {list(P.values())}  {list(NZ.values())}  {T}  2^{np.log2(rel):6.2f}  {p:5d}  {a:5d}  '
              f'{list(w.values())}  {ent}', flush=True)
