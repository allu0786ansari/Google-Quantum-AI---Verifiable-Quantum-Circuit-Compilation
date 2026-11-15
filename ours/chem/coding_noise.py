"""How much error does forcing a bit-plane into RM(1,4) per 16-address block cost, and how much of it is
avoided by (a) sensitivity-weighted decoding and (b) reordering the vectors inside each group?

Input: a public factorization table (U, W, angle_words, tree_edges, codes). Own code, numpy only.
Sensitivity of the two-body tensor to an angle error at (vector a, tree node t):
    s(a, t) = r(a, t) * || sum_c W[a, c] L_c ||_F,   r = amplitude entering the node.
Cost model of flipping bit k of a word: s^2 * 4^k (in squared angle steps); plain 8-bit rounding is s^2 / 12.
    python coding_noise.py table.npz
"""
import sys, itertools
import numpy as np

t = np.load(sys.argv[1])
U = t['U'].astype(np.float64); R, B, N = U.shape
U /= np.linalg.norm(U, axis=2, keepdims=True)
W = t['W'].astype(np.float64); edges = t['tree_edges']; words = t['angle_words'].astype(np.int64)
codes = t['codes'] if 'codes' in t else np.zeros((75, 8), np.int64)
L = np.einsum('rbp,rbq,rbc->rcpq', U, U, W)
g = np.linalg.norm(np.einsum('rbc,rcpq->rbpq', W, L).reshape(R, B, -1), axis=2).reshape(-1)
u = np.abs(U.reshape(-1, N)).copy(); r = np.zeros((R * B, 75))
for tt in range(74, -1, -1):
    p, c = edges[tt]; r[:, tt] = np.hypot(u[:, p], u[:, c]); u[:, p] = r[:, tt]; u[:, c] = 0
s2 = (r * g[:, None]) ** 2                               # (1144, 75)
D0 = s2.sum() / 12
print(f'vectors {R*B}, rounding variance D0 = {D0:.4e} (reference = 1.0 below)')
print('already coded (column, bit) pairs:', int((codes > 0).sum()), 'of', codes.size)

pts = np.array([[(x >> i) & 1 for i in range(4)] for x in range(16)])
RM1 = np.array([(c0 + pts @ np.array([c1, c2, c3, c4])) % 2 for c0, c1, c2, c3, c4 in itertools.product((0, 1), repeat=5)], np.uint8)  # 32 x 16


def plane_cost(order, k, weighted):
    """Total cost of forcing bit k of every uncoded column into RM(1,4) per block, for a vector order."""
    w_all = np.vstack([words[:76], words[76:][order]])
    s_all = np.vstack([np.full((76, 75), np.inf), s2[order]])          # one-body rows are never changed
    n = len(w_all); H = -(-n // 16)
    pad = H * 16 - n
    w_all = np.vstack([w_all, np.zeros((pad, 75), np.int64)]); s_all = np.vstack([s_all, np.zeros((pad, 75))])
    bits = ((w_all >> k) & 1).reshape(H, 16, 75)
    wt = s_all.reshape(H, 16, 75) * 4.0 ** k
    cols = [c for c in range(75) if codes[c, k] == 0]
    total = 0.0; flips = 0; infeasible = 0
    for c in cols:
        b = bits[:, :, c]; w = wt[:, :, c]
        diff = b[:, None, :] ^ RM1[None, :, :]                          # (H, 32, 16)
        if weighted:
            wf = np.where(np.isinf(w), 1e30, w)
            cost = (diff * wf[:, None, :]).sum(2)
            best = cost.argmin(1)
        else:
            best = diff.sum(2).argmin(1)
        d = diff[np.arange(H), best]
        wfin = np.where(np.isinf(w), 0, w)
        infeasible += int((d & np.isinf(w)).sum())
        total += (d * wfin).sum(); flips += int(d.sum())
    return total / D0, flips, len(cols), infeasible


ident = np.arange(R * B)
by_g = np.concatenate([rho * B + np.argsort(g[rho * B:(rho + 1) * B]) for rho in range(R)])   # importance-sorted inside each group
print('\nforcing one more bit-plane into RM(1,4) on all not-yet-coded columns: added variance / D0')
print(' bit   order        decoding    variance/D0    flips   columns')
for k in (0, 1, 2):
    for oname, order in (('as trained', ident), ('by importance', by_g)):
        for wname, wflag in (('plain', False), ('weighted', True)):
            v, fl, nc, inf = plane_cost(order, k, wflag)
            print(f'  {k}    {oname:13} {wname:9}  {v:10.3f}   {fl:7d}   {nc}' + (f'   (touches {inf} one-body bits)' if inf else ''))
