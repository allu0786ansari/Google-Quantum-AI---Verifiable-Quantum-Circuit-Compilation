"""Word-level model of our DQI key-equation layer: normal-schedule, register-shared Euclid.

Two pairs (A, B), each 71 words: A = (r_a, t_a), B = (r_b, t_b), deg r_b = deg r_a + 1, deg t_b = deg t_a - 1.
Start: A = (S, 1), B = (x^70, 0). Round i = 1..35 (Da = 70 - i, dt = i - 1):
  f_i = f_{i-1} or [A.r[Da] == 0]                      (sticky freeze flag; k < 35 errors collapse r_k)
  q1 = B.r[Da+1] / A.r[Da];  q1e = q1 and not f_i;   B.r[c] ^= q1e * A.r[c-1]   (c = 1..Da+1)
  q0 = B.r[Da]   / A.r[Da];  q0e = q0 and not f_i;   B.r[c] ^= q0e * A.r[c]     (c = 0..Da)
  the two top words of B.r (now zero) become B.t[dt], B.t[dt+1]
  B.t[c] ^= q1e * A.t[c-1] (c = 1..dt+1);  B.t[c] ^= q0e * A.t[c] (c = 0..dt)
  q1e, q0e are uncomputed from B.t[dt+1], B.t[dt] (gated by not f_i); roles swap.
After round 35: if the number of frozen rounds is odd, swap A and B (word matching below).
Then sigma = A.t (36 words), Omega = A.r (35 words); e_j = a^-j Omega(a^-j) / sigma_odd(a^-j).

Physical words are integers in a pool; the model tracks which word holds which coefficient, exactly as
the circuit builder will, and asserts every relabelled word is zero when it is supposed to be.
"""
import numpy as np
from gf import mul, div, inv, apow, sample, syndromes, weight_probs, N_SYN

T = 35


def run(S, check_zero=True):
    W = {}                                   # physical word id -> value
    def new(v=0):
        w = len(W); W[w] = v; return w
    Ar = [None] * 70; At = []
    for i in range(70):
        Ar[i] = new(S[i])
    At = [new(1)]
    Br = [new(0) for _ in range(70)] + [new(1)]; Bt = []
    flags = []
    f = 0
    for i in range(1, T + 1):
        Da, dt = 70 - i, i - 1
        assert len(Ar) == Da + 1 and len(At) == dt + 1 and len(Br) == Da + 2 and len(Bt) == dt
        z = int(W[Ar[Da]] == 0)
        f = f | z; flags.append(f)
        q1 = div(W[Br[Da + 1]], W[Ar[Da]]); q1e = 0 if f else q1
        for c in range(1, Da + 2):
            W[Br[c]] ^= mul(q1e, W[Ar[c - 1]])
        q0 = div(W[Br[Da]], W[Ar[Da]]); q0e = 0 if f else q0
        for c in range(0, Da + 1):
            W[Br[c]] ^= mul(q0e, W[Ar[c]])
        top1 = Br.pop(); top0 = Br.pop()
        if check_zero and not f:
            assert W[top1] == 0 and W[top0] == 0
        Bt = Bt + [top0, top1]
        for c in range(1, dt + 2):
            W[Bt[c]] ^= mul(q1e, W[At[c - 1]])
        for c in range(0, dt + 1):
            W[Bt[c]] ^= mul(q0e, W[At[c]])
        if not f:   # uncompute check: q's recoverable from the cofactor's top words
            q1r = div(W[Bt[dt + 1]], W[At[dt]])
            u = mul(q1r, W[At[dt - 1]]) if dt >= 1 else 0
            q0r = div(W[Bt[dt]] ^ u, W[At[dt]])
            assert (q1r, q0r) == (q1e, q0e), (i, q1r, q1e, q0r, q0e)
        Ar, At, Br, Bt = Br, Bt, Ar, At
    parity = sum(flags) & 1
    if parity:
        # A has (35 r, 36 t), B has (36 r, 35 t): match A.r[c]<->B.r[c] (c<35), A.t[c]<->B.t[c] (c<35),
        # A.t[35] <-> B.r[35]
        assert len(Ar) == 35 and len(At) == 36 and len(Br) == 36 and len(Bt) == 35
        pairs = [(Ar[c], Br[c]) for c in range(35)] + [(At[c], Bt[c]) for c in range(35)] + [(At[35], Br[35])]
        for a, b in pairs:
            W[a], W[b] = W[b], W[a]
    sigma = [W[w] for w in At]; omega = [W[w] for w in Ar]
    return sigma, omega, flags, W, (Ar, At, Br, Bt)


def phase_from(sigma, omega, M):
    tot = 0
    for j in range(4095):
        x = apow(-j)
        s = 0; so = 0; om = 0; xp = 1
        for c, v in enumerate(sigma):
            t = mul(v, xp); s ^= t
            if c & 1: so ^= t
            xp = mul(xp, x)
        xp = x
        for v in omega:
            om ^= mul(v, xp); xp = mul(xp, x)
        if s == 0:
            e = div(om, so)
            tot ^= bin(e & M[j]).count('1') & 1
    return tot


if __name__ == '__main__':
    import sys
    rng = np.random.default_rng(int(sys.argv[1]) if len(sys.argv) > 1 else 5)
    p = weight_probs()
    N = 300; ok = 0; fails = []
    for n in range(N):
        pos, val = sample(rng, p)
        S = syndromes(pos, val)
        try:
            sigma, omega, flags, W, _ = run(S)
        except AssertionError as ex:
            fails.append(('assert', len(pos), str(ex)[:60])); continue
        # compare with the true locator up to scale
        true_sig = [1]
        for j in pos:
            Xj = apow(int(j)); nxt = [0] * (len(true_sig) + 1)
            for c, v in enumerate(true_sig):
                nxt[c] ^= v; nxt[c + 1] ^= mul(v, Xj)
            true_sig = nxt
        lam = sigma[0]
        good = lam != 0 and all(sigma[c] == mul(lam, true_sig[c] if c < len(true_sig) else 0) for c in range(36))
        if n < 40:
            M = rng.integers(0, 4096, size=4095)
            want = 0
            for j, e in zip(pos, val):
                want ^= bin(int(e) & int(M[j])).count('1') & 1
            good = good and phase_from(sigma, omega, M) == want
        ok += good
        if not good: fails.append(('wrong', len(pos), sum(flags)))
    print(f'{ok}/{N} correct; failures: {fails[:10]}')
