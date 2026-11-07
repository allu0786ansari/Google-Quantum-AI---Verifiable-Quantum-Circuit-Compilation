"""GF(2^12) arithmetic matching km.GF2Field(12) (log/antilog tables) and test-distribution sampling."""
import numpy as np
import kickmix as km
F = km.GF2Field(12)
Q = 1 << 12
ALPHA = int(F.primitive_element)
EXP = np.zeros(2 * Q, np.int64); LOG = np.zeros(Q, np.int64)
v = 1
for i in range(Q - 1):
    EXP[i] = v; LOG[v] = i
    v = int(F.mul(v, ALPHA))
for i in range(Q - 1, 2 * Q): EXP[i] = EXP[i - (Q - 1)]
assert v == 1


def mul(a, b):
    return 0 if a == 0 or b == 0 else int(EXP[LOG[a] + LOG[b]])


def inv(a):
    return 0 if a == 0 else int(EXP[(Q - 1 - LOG[a]) % (Q - 1)])


def div(a, b):
    return mul(a, inv(b))


def apow(e):
    return int(EXP[e % (Q - 1)])


M_POS, N_SYN, L_MAX = 4095, 70, 35


def weight_probs(r=2016):
    m, l, q = M_POS, L_MAX, Q
    A = np.zeros((l + 1, l + 1))
    for k in range(l + 1): A[k, k] = k * (q - 2 * r) / np.sqrt(r * (q - r))
    for k in range(1, l + 1): A[k - 1, k] = A[k, k - 1] = np.sqrt(k * (m - k + 1))
    w, vv = np.linalg.eigh(A); p = vv[:, -1] ** 2
    return p / p.sum()


def sample(rng, p=None):
    p = weight_probs() if p is None else p
    k = rng.choice(len(p), p=p)
    pos = rng.choice(M_POS, size=k, replace=False)
    val = rng.integers(1, Q, size=k)
    return pos, val


def syndromes(pos, val):
    return [int(np.bitwise_xor.reduce([mul(int(e), apow((i + 1) * int(j))) for j, e in zip(pos, val)])) if len(pos) else 0
            for i in range(N_SYN)]
