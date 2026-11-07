"""Per-stage precision model for the Horner circuit (generalises poly_model).

acc_k is held in units 2^-TT[k] (k = 3 innermost ... 0 = Y). Stage k -> k-1 subtracts, for each of the
top NZ[k] bits z_j of z, the row  (acc_k >> d) << off,  d = NZ[k] - j + TT[k] - TT[k-1]  (off = -d if d < 0).
This is exactly what the circuit's row adders compute. Fitted coefficients are rounded per stage; c0 is
refitted last on the exact integer pipeline (minimax relative error per interval).

Returns error and an estimated cost: rows cost (n_row ANDs + target slice adder), Q estimate.
"""
import numpy as np
from model_l1_l2 import xn, ytrue

K = 7; NZB = 22 - K
IDX = xn >> NZB
Z = (xn & ((1 << NZB) - 1)).astype(np.int64)
ENT = np.unique(IDX)
_coef = {}


def lsq(deg):
    if deg not in _coef:
        zf = Z / float(1 << NZB)
        co = np.zeros((1 << K, deg + 1))
        order = np.argsort(IDX, kind='stable'); bounds = np.searchsorted(IDX[order], ENT)
        for n, i in enumerate(ENT):
            sel = order[bounds[n]: bounds[n + 1] if n + 1 < len(ENT) else len(order)]
            A = np.vander(zf[sel], deg + 1, increasing=True)
            co[i] = np.linalg.lstsq(A / ytrue[sel, None], np.ones(len(sel)), rcond=None)[0]
        _coef[deg] = (co, order, bounds)
    return _coef[deg]


def rows_sum(acc, wacc, nz, d0):
    """sum_j z_j * row_j with row_j = (acc >> d) << off, d = nz - j + d0. Returns (sum, rows meta)."""
    Zk = Z >> (NZB - nz)
    tot = np.zeros(len(acc), np.int64); meta = []
    for j in range(nz):
        d = nz - j + d0
        if d >= wacc:
            continue
        row = (acc >> d) if d >= 0 else (acc << (-d))
        tot += ((Zk >> j) & 1) * row
        meta.append((j, max(d, 0), max(-d, 0)))     # (z bit, src low bit, target offset)
    return tot, meta


def fit(TT, NZ, deg=3, tables=False):
    co, order, bounds = lsq(deg)
    C = {k: np.round(co[:, k] * (-1) ** k * 2.0 ** TT[k - 1 if k else 0]).astype(np.int64) for k in range(1, deg + 1)}
    # c_k (k >= 1) is loaded into acc_{k-1}'s register? No: c_k starts stage k's target, i.e. units TT[k]... see below
    # Convention: acc_deg = c_deg in units TT[deg]; acc_{k-1} = c_{k-1} (units TT[k-1]) - z * acc_k.
    C = {k: np.round(co[:, k] * (-1) ** k * 2.0 ** TT[k]).astype(np.int64) for k in range(1, deg + 1)}
    for k in C:
        assert C[k][ENT].min() >= 0
    acc = C[deg][IDX]; meta_all = {}; widths = {deg: int(C[deg][ENT].max()).bit_length()}
    for k in range(deg, 0, -1):
        wacc = int(acc.max()).bit_length()
        tot, meta = rows_sum(acc, wacc, NZ[k], TT[k] - TT[k - 1])
        meta_all[k] = meta
        if k > 1:
            acc = C[k - 1][IDX] - tot
            assert acc.min() >= 0
            widths[k - 1] = int(C[k - 1][ENT].max()).bit_length()
        else:
            prod1 = tot
    rest = -prod1 / 2.0 ** TT[0]
    C0 = np.zeros(1 << K, np.int64)
    for n, i in enumerate(ENT):
        sel = order[bounds[n]: bounds[n + 1] if n + 1 < len(ENT) else len(order)]
        yt = ytrue[sel]; r = rest[sel]
        # minimise max |c/2^T + r - yt| / yt over integer c: solve the 1-D minimax exactly
        lo = (yt - r) * 2.0 ** TT[0]
        cands = np.arange(int(np.floor(lo.min())) - 1, int(np.ceil(lo.max())) + 2)
        best = min(cands, key=lambda c: np.abs((c - lo) / yt).max())
        C0[i] = best
    Yint = C0[IDX] - prod1
    rel = np.abs(Yint / 2.0 ** TT[0] / ytrue - 1).max()
    widths[0] = int(C0[ENT].max()).bit_length()
    if tables:
        return dict(C={0: C0, **C}, Y=Yint, rel=rel, widths=widths, meta=meta_all, TT=TT, NZ=NZ)
    return rel, widths, meta_all


def est_cost(widths, meta, TT):
    """Forward rows: sum over rows of (row width ANDs + adder over target[off:] ~ width-1)."""
    tw = {k - 1: widths[k - 1] - (1 if k == 1 else 0) for k in meta}      # Y stores one bit less
    sw = {k: widths[k] for k in meta}
    tof = 0; peak_tmp = 0
    for k, rows in meta.items():
        for j, lo, off in rows:
            n = sw[k] - lo; w = tw[k - 1] - off
            tof += n + w - 1; peak_tmp = max(peak_tmp, n + w)
    return tof


if __name__ == '__main__':
    import itertools, time
    t0 = time.time()
    r = fit({0: 29, 1: 29, 2: 29, 3: 29}, {1: 15, 2: 14, 3: 11})
    print('check vs poly_model (2^-25.29):', np.log2(r[0]), r[1], f'{time.time() - t0:.0f}s')
    res = []
    for T0 in (28, 29):
        for T1 in range(T0 - 4, T0 + 1):
            for T2 in range(T1 - 6, T1 + 1, 2):
                for T3 in range(T2 - 6, T2 + 1, 2):
                    for nz1, nz2, nz3 in [(15, 14, 11), (15, 13, 10), (15, 12, 9), (15, 12, 8), (15, 11, 8), (14, 12, 9)]:
                        TT = {0: T0, 1: T1, 2: T2, 3: T3}; NZ = {1: nz1, 2: nz2, 3: nz3}
                        try:
                            rel, w, meta = fit(TT, NZ)
                        except AssertionError:
                            continue
                        if np.log2(rel) > -25.1:
                            continue
                        tof = est_cost(w, meta, TT)
                        live = 26 + w[2] + w[1] + max(w[0] - 1 + 30, 39)   # stage-1 peak ~ +tmp+ws vs shift
                        res.append((live * (2 * tof + 1200), np.log2(rel), TT, NZ, w, tof))
    res.sort(key=lambda r: r[0])
    for r in res[:12]:
        print(f'{r[0]:>10,.0f}  err 2^{r[1]:.2f}  T {list(r[2].values())}  NZ {list(r[3].values())}  widths {r[4]}  rows fwd {r[5]}')
    print(f'{len(res)} passing, {time.time() - t0:.0f}s')
