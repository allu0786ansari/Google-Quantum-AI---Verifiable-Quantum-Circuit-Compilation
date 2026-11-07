"""Scan a compiled circuit for windows (<= K wires) whose net action is Clifford: an affine permutation
with a phase polynomial of degree <= 2. Such windows can be rebuilt with X/CX/SWAP/CZ/Z only, so any
Toffolis inside are removable. Exact for all 2^k inputs of the window wires.

    python rewrite_scan.py OPS_PREFIX STATS.npy|- SHOTS [KMAX]
"""
import sys, json, itertools
import numpy as np

prefix, stats, shots = sys.argv[1], sys.argv[2], float(sys.argv[3])
KMAX = int(sys.argv[4]) if len(sys.argv) > 4 else 8
import os
full = np.memmap(prefix + '.bin', dtype=np.int32, mode='r').reshape(-1, 5)
LO, HI = int(os.environ.get('LO', 0)), int(os.environ.get('HI', len(full)))
ops = np.array(full[LO:HI])
n = len(ops)
act = np.ones(n)
if stats != '-':
    st = np.load(stats); act = np.zeros(n); g = st[:, 0].astype(np.int64); m = (g >= LO) & (g < HI); act[g[m] - LO] = st[m, 1] / shots
QUBIT_ARITY = {1: 2, 2: 3, 3: 1, 4: 2, 5: 3, 13: 1, 14: 2, 15: 3}   # pure gates we can absorb
TOF = {2, 5, 15}
# conditioned (by bit or by PUSH_CONDITION stack) gates are barriers
d = np.cumsum((ops[:, 0] == 11).astype(np.int64) - (ops[:, 0] == 12).astype(np.int64))
cond = (d > 0) | (ops[:, 0] == 12) | (ops[:, 4] >= 0)
wires_of = [tuple(int(x) for x in ops[i, 1:1 + QUBIT_ARITY[ops[i, 0]]]) if ops[i, 0] in QUBIT_ARITY
            else (tuple(int(x) for x in ops[i, 1:2]) if ops[i, 0] in (6, 7, 17) else ()) for i in range(n)]


def simulate(idx, W):
    k = len(W); pos = {w: j for j, w in enumerate(W)}
    x = np.arange(1 << k)
    s = ((x[None, :] >> np.arange(k)[:, None]) & 1).astype(bool)
    ph = np.zeros(1 << k, bool)
    for i in idx:
        op = ops[i, 0]; q = [pos[w] for w in wires_of[i]]
        if op == 1: s[q[1]] ^= s[q[0]]
        elif op == 2: s[q[2]] ^= s[q[0]] & s[q[1]]
        elif op == 3: s[q[0]] = ~s[q[0]]
        elif op == 4: s[[q[0], q[1]]] = s[[q[1], q[0]]]
        elif op == 5:
            c = s[q[0]]; a, b = s[q[1]].copy(), s[q[2]].copy(); s[q[1]] = np.where(c, b, a); s[q[2]] = np.where(c, a, b)
        elif op == 13: ph ^= s[q[0]]
        elif op == 14: ph ^= s[q[0]] & s[q[1]]
        elif op == 15: ph ^= s[q[0]] & s[q[1]] & s[q[2]]
    return s, ph


def anf_degree(tt, k):
    a = tt.astype(np.uint8).copy()
    for j in range(k):
        a = a.reshape(-1, 2, 1 << j); a[:, 1, :] ^= a[:, 0, :]; a = a.reshape(-1)
    nz = np.flatnonzero(a)
    return max((bin(int(m)).count('1') for m in nz), default=0)


def clifford(idx, W):
    s, ph = simulate(idx, W); k = len(W)
    return all(anf_degree(s[j], k) <= 1 for j in range(k)) and anf_degree(ph, k) <= 2


found = []
tof_idx = [i for i in range(n) if ops[i, 0] in TOF and not cond[i]]
for t in tof_idx:
    W = list(wires_of[t]); idx = [t]; blocked = set(); best = None
    for j in range(t + 1, min(n, t + 4000)):
        wj = wires_of[j]
        if not wj or not (set(wj) & set(W)) and not (set(wj) & blocked):
            if wj and set(wj) & blocked:
                pass
            continue
        if ops[j, 0] not in QUBIT_ARITY or cond[j] or set(wj) & blocked:
            blocked |= set(wj) & set(W)
            if blocked >= set(W):
                break
            continue
        newW = W + [w for w in wj if w not in W]
        if len(newW) > KMAX:
            blocked |= set(wj) & set(W)
            if blocked >= set(W):
                break
            continue
        W = newW; idx.append(j)
        if ops[j, 0] in TOF and clifford(idx, W):
            saved = sum(act[i] for i in idx if ops[i, 0] in TOF)
            best = (saved, list(idx), list(W))
    if best:
        found.append(best)
found.sort(key=lambda b: -b[0])
used = set(); chosen = []
for saved, idx, W in found:
    if used.isdisjoint(idx):
        chosen.append((saved, idx, W)); used |= set(idx)
total = sum(c[0] for c in chosen)
print(f'Toffolis scanned {len(tof_idx)}; Clifford windows found {len(found)}; disjoint {len(chosen)}; '
      f'Toffolis removable per shot {total:.3f}')
for saved, idx, W in chosen[:15]:
    print(f'  saves {saved:.3f}/shot: {len(idx)} gates on {len(W)} wires {W}, ops {idx[0]}..{idx[-1]}')
json.dump([dict(saved=s, idx=[int(i) for i in idx], wires=[int(w) for w in W]) for s, idx, W in chosen],
          open(prefix + '_clifford_windows.json', 'w'))
