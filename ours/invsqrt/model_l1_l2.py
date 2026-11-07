"""Classical models for Benchmark I layers I-1 (normalisation) and I-2 (initial guess).

I-1: binary leading-zero-pair count s (4 bits) by stages of 8/4/2/1 pairs; x_norm = x << 2s in [2^20, 2^22).
I-2: piecewise seed for y0 ~ 1/sqrt(x_norm / 2^20) in (1/2, 1]; index = top k bits of x_norm,
     t = next m bits; degree 0 (table), 1 (c0 + c1 t) or 2. Coefficients minimax-fitted then rounded
     to w fractional bits. Reports worst relative seed error and the error after one exact Newton step.
"""
import numpy as np

MAXX = 3133452
x = np.arange(1, MAXX + 1, dtype=np.int64)

# ---- I-1 -------------------------------------------------------------------------------
def normalise(x):
    v = x.copy(); s = np.zeros_like(x)
    for pairs in (8, 4, 2, 1):
        d = 2 * pairs
        top_zero = (v >> (22 - d)) == 0          # top d bits of the 22-bit register are zero
        v = np.where(top_zero, v << d, v); s += np.where(top_zero, pairs, 0)
    return v, s

xn, s = normalise(x)
assert ((xn >> 20) > 0).all() and (xn < (1 << 22)).all(), 'x_norm out of [2^20, 2^22)'
assert (xn == x << (2 * s)).all() and s.max() <= 10
print(f'I-1 ok: s in [{s.min()}, {s.max()}] for all {MAXX:,} x; x_norm in [2^20, 2^22)')

# ---- I-2 -------------------------------------------------------------------------------
u = xn / float(1 << 20)                          # in [1, 4)
ytrue = 1 / np.sqrt(u)


def seed_error(k, m, deg, w):
    idx = xn >> (22 - k)                          # top k bits
    t = ((xn >> (22 - k - m)) & ((1 << m) - 1)) if m else np.zeros_like(xn)
    tf = t / float(1 << m) if m else np.zeros(len(xn))
    y0 = np.empty(len(xn))
    for i in np.unique(idx):
        sel = idx == i
        if deg == 0:
            lo, hi = ytrue[sel].min(), ytrue[sel].max()
            c = 2 * lo * hi / (lo + hi)           # minimax relative constant
            y0[sel] = np.round(c * (1 << w)) / (1 << w)
        else:
            A = np.vander(tf[sel], deg + 1, increasing=True)
            # minimise max relative error approx via weighted least squares then round
            coef = np.linalg.lstsq(A / ytrue[sel, None], np.ones(sel.sum()), rcond=None)[0]
            coef = np.round(coef * (1 << w)) / (1 << w)
            y0[sel] = A @ coef
    d0 = y0 / ytrue - 1
    d1 = -1.5 * d0 ** 2 - 0.5 * d0 ** 3              # exact one-step Newton relative error
    return np.abs(d0).max(), np.abs(d1).max(), len(np.unique(idx))


print('\nI-2 seed options (worst relative error over all x; budget after Newton: < 2^-24 = 5.96e-8, keep <= 2^-25.5 for arithmetic)')
print(' deg  k   m   w   entries   seed err (log2)   after Newton (log2)')
for deg, k, m, w in [(0, 10, 0, 14), (0, 12, 0, 15), (0, 13, 0, 16),
                     (1, 5, 8, 16), (1, 6, 8, 16), (1, 7, 7, 16), (1, 7, 8, 17), (1, 8, 7, 17),
                     (2, 3, 10, 18), (2, 4, 9, 18), (2, 5, 8, 18)]:
    e0, e1, ent = seed_error(k, m, deg, w)
    print(f' {deg}   {k:2d}  {m:2d}  {w:2d}  {ent:6d}      {np.log2(e0):7.2f}            {np.log2(e1):7.2f}')
