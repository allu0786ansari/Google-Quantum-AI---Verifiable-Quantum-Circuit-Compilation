"""Flip-graph search (Kauers & Moosbauer) on rank-1 pieces a (x) b (x) c of a loader tensor.

Pieces (a: 15-bit low form, b: 85-bit high form, c: 35-bit fan-out). Flips between two pieces
that share a factor keep the tensor and the piece count; a reduction merges two pieces that
share two factors (or drops a zero piece). Cost in Toffoli terms = sum over distinct c of
ceil(#pieces with that c / 2), since two pieces with one fan-out set form one paired term.

    python flipgraph.py SLOT [FLIPS] [SEED]
"""
import random
import sys
import time
from collections import defaultdict

import numpy as np

from load import slot

L, H, O = 15, 85, 35


def pieces_from_terms(terms):
    P = []
    for a1, b1, a2, b2, s in terms:
        P.append([a1, b1, s])
        if a2 and b2:
            P.append([a2, b2, s])
    return P


def tensor_of(P):
    T = np.zeros((L, H, O), np.uint8)
    for a, b, c in P:
        av = np.array([(a >> i) & 1 for i in range(L)], np.uint8)
        bv = np.array([(b >> i) & 1 for i in range(H)], np.uint8)
        cv = np.array([(c >> i) & 1 for i in range(O)], np.uint8)
        T ^= av[:, None, None] & bv[None, :, None] & cv[None, None, :]
    return T


def term_cost(P):
    cnt = defaultdict(int)
    for _, _, c in P:
        cnt[c] += 1
    return sum((n + 1) // 2 for n in cnt.values())


def reduce_all(P):
    """Merge pieces sharing two factors; drop zero pieces. Repeat until stable."""
    changed = True
    while changed:
        changed = False
        P[:] = [p for p in P if p[0] and p[1] and p[2]]
        for (x, y, z) in ((0, 1, 2), (0, 2, 1), (1, 2, 0)):
            seen = {}
            for i, p in enumerate(P):
                key = (p[x], p[y])
                if key in seen:
                    j = seen.pop(key)
                    P[j][z] ^= p[z]
                    p[0] = 0          # mark removed
                    changed = True
                else:
                    seen[key] = i
            if changed:
                break
    return P


def walk(P, flips, rng, report=500000, t0=None):
    best = (len(P), term_cost(P), [list(p) for p in P])
    for step in range(flips):
        axis = rng.randrange(3)
        buckets = defaultdict(list)
        # sampling a shared-factor pair: pick random piece, find partners sharing factor `axis`
        i = rng.randrange(len(P))
        key = P[i][axis]
        partners = [j for j in range(len(P)) if j != i and P[j][axis] == key]
        if not partners:
            continue
        j = rng.choice(partners)
        o1, o2 = [ax for ax in range(3) if ax != axis]
        if rng.random() < 0.5:
            o1, o2 = o2, o1
        # (k, x_i, y_i) + (k, x_j, y_j) -> (k, x_i + x_j, y_i) + (k, x_j, y_i + y_j)
        P[i][o1] ^= P[j][o1]
        P[j][o2] ^= P[i][o2]
        n0 = len(P)
        reduce_all(P)
        if len(P) < n0 or step % 1000 == 0:
            tc = term_cost(P)
            if (tc, len(P)) < (best[1], best[0]):
                best = (len(P), tc, [list(p) for p in P])
        if report and step % report == 0:
            print(f'  flip {step}: pieces {len(P)} terms {term_cost(P)} | best pieces {best[0]} terms {best[1]}'
                  f' ({time.time() - t0:.0f}s)', flush=True)
    return best


if __name__ == '__main__':
    k = int(sys.argv[1])
    flips = int(sys.argv[2]) if len(sys.argv) > 2 else 2000000
    seed = int(sys.argv[3]) if len(sys.argv) > 3 else 0
    plan, delta, C, terms = slot(k)
    T = C[1:, 1:]
    P = pieces_from_terms(terms)
    assert np.array_equal(tensor_of(P), T)
    print(f'slot {k}: frontier {len(terms)} terms, {len(P)} pieces, term cost {term_cost(P)}', flush=True)
    rng = random.Random(seed)
    best = walk(P, flips, rng, t0=time.time())
    assert np.array_equal(tensor_of(best[2]), T)
    print(f'slot {k}: best pieces {best[0]}, terms {best[1]} (frontier {len(terms)})')
