"""How often is the EEA on (x^70, S) 'normal' (every round: degree-1 quotient, remainder degree -1) until collapse at round k?"""
import numpy as np
from gf import *


def deg(p):
    for i in range(len(p) - 1, -1, -1):
        if p[i]: return i
    return -1


def polymod_rounds(S):
    a = [0] * 70 + [1]; b = list(S)          # r_{-1} = x^70, r_0 = S
    degs = [deg(b)]
    while deg(b) >= 35:
        # a mod b
        a = a[:]; db = deg(b)
        qdeg = deg(a) - db
        for s in range(deg(a) - db, -1, -1):
            c = a[db + s]
            if c:
                f = div(c, b[db])
                for i in range(db + 1): a[i + s] ^= mul(f, b[i])
        a, b = b, a
        degs.append(deg(b))
    return degs


rng = np.random.default_rng(1)
N = 4000; normal = 0; bad = []
for n in range(N):
    pos, val = sample(rng); k = len(pos)
    S = syndromes(pos, val)
    d = polymod_rounds(S)
    expect = [69 - i for i in range(k)]
    ok = d[:k] == expect and len(d) == k + 1 and d[k] < k
    normal += ok
    if not ok: bad.append((k, d[:8], len(d)))
print(f'normal-case shots: {normal}/{N} = {normal / N:.5f}; failure {1 - normal / N:.5f}')
print(bad[:5])
