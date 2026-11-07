"""Search for cheaper paired-product decompositions of a loader's mixed ANF tensor.

T (15 low x 85 high x 35 out, GF(2)) = sum_g M_g (x) c_g; cost = sum_g ceil(rank(M_g) / 2), the
number of Toffoli terms (a rank-2 matrix a1 b1^T + a2 b2^T with one fan-out set is one term).
Flip moves keep T exact:  M_g c_g + M_h c_h = M_g (c_g + c_h) + (M_g + M_h) c_h.
Groups with equal c are merged. Low-cache scheduling is ignored here (checked later).

    python anneal.py SLOT [STEPS] [SEED]
"""
import math
import random
import sys
import time

import numpy as np

from load import slot

L, H, O = 15, 85, 35


def rank(rows):
    basis = []
    for r in rows:
        for b in basis:
            r = min(r, r ^ b)
        if r:
            basis.append(r)
    return len(basis)


def mat_from_term(a1, b1, a2, b2):
    return tuple(((b1 if (a1 >> i) & 1 else 0) ^ (b2 if (a2 >> i) & 1 else 0)) for i in range(L))


def add(M, N):
    return tuple(x ^ y for x, y in zip(M, N))


def tensor_of(groups):
    T = np.zeros((L, H, O), np.uint8)
    for M, c in groups:
        Mm = np.array([[(row >> j) & 1 for j in range(H)] for row in M], np.uint8)
        cv = np.array([(c >> o) & 1 for o in range(O)], np.uint8)
        T ^= Mm[:, :, None] & cv[None, None, :]
    return T


def cost_of(M):
    return (rank(M) + 1) // 2


def run(k, steps, seed):
    rng = random.Random(seed)
    plan, delta, C, terms = slot(k)
    T = C[1:, 1:]
    groups = {}
    for a1, b1, a2, b2, s in terms:
        M = mat_from_term(a1, b1, a2, b2)
        groups[s] = add(groups[s], M) if s in groups else M
    G = [(M, c) for c, M in groups.items() if any(M)]
    assert np.array_equal(tensor_of(G), T)
    costs = [cost_of(M) for M, _ in G]
    total = sum(costs)
    best = total
    best_G = list(G)
    print(f'slot {k}: frontier terms {len(terms)}, as merged groups cost {total} ({len(G)} groups)', flush=True)
    t0 = time.time()
    temp0, temp1 = 0.6, 0.05
    for step in range(steps):
        temp = temp0 * (temp1 / temp0) ** (step / steps)
        n = len(G)
        g, h = rng.randrange(n), rng.randrange(n)
        if g == h:
            continue
        (Mg, cg), (Mh, ch) = G[g], G[h]
        if rng.random() < 0.5:
            # M_h += M_g ; c_g += c_h
            newMh = add(Mh, Mg)
            newcg = cg ^ ch
            new_h_cost = cost_of(newMh) if any(newMh) else 0
            delta_c = new_h_cost - costs[h]
            if newcg == 0:
                delta_c -= costs[g]
            if delta_c <= 0 or rng.random() < math.exp(-delta_c / temp):
                G[h] = (newMh, ch); costs[h] = new_h_cost
                G[g] = (Mg, newcg)
                drop = [i for i in (g, h) if G[i][1] == 0 or not any(G[i][0])]
                for i in sorted(drop, reverse=True):
                    G.pop(i); costs.pop(i)
                total += delta_c
        else:
            # merge-and-resplit-free move: swap roles (dual flip) M_g += M_h ; c_h += c_g
            newMg = add(Mg, Mh)
            newch = ch ^ cg
            new_g_cost = cost_of(newMg) if any(newMg) else 0
            delta_c = new_g_cost - costs[g]
            if newch == 0:
                delta_c -= costs[h]
            if delta_c <= 0 or rng.random() < math.exp(-delta_c / temp):
                G[g] = (newMg, cg); costs[g] = new_g_cost
                G[h] = (Mh, newch)
                drop = [i for i in (g, h) if G[i][1] == 0 or not any(G[i][0])]
                for i in sorted(drop, reverse=True):
                    G.pop(i); costs.pop(i)
                total += delta_c
        if total < best:
            best, best_G = total, list(G)
        if step % 200000 == 0:
            print(f'  step {step} temp {temp:.3f} cost {total} best {best} groups {len(G)} ({time.time() - t0:.0f}s)', flush=True)
    assert np.array_equal(tensor_of(best_G), T), 'tensor changed'
    print(f'slot {k}: best cost {best} vs frontier {len(terms)}', flush=True)
    return best, best_G


if __name__ == '__main__':
    k = int(sys.argv[1])
    steps = int(sys.argv[2]) if len(sys.argv) > 2 else 1000000
    seed = int(sys.argv[3]) if len(sys.argv) > 3 else 0
    run(k, steps, seed)
