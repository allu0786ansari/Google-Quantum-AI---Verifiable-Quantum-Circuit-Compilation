"""Our addition for dfthcsos: code bit 1 of the angle words into RM(1,4) per 16-address block, choosing for
every block the code word that disturbs the Hamiltonian least (sensitivity-weighted nearest code word).

Bit 0 of the public table is already coded this way; bit 1 was left free because plain nearest-code-word
rounding costs 11x the rounding variance. Weighted, it costs 1.7x (coding_noise.py).

    VERIFIERS=<dir> python code_bit1.py table.npz out.npz [first_batch] [last_batch]

Columns of batches first..last (5 columns each, default all 15) get bit 1 coded. The weights (w, W) are not
changed here; refit them afterwards (tools/experiments/refit.py) to restore the 20 mHa limit.
"""
import itertools, os, sys
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', 'tools', 'experiments'))
src, out = sys.argv[1], sys.argv[2]
b0, b1 = (int(sys.argv[3]), int(sys.argv[4])) if len(sys.argv) > 4 else (0, 14)
z = dict(np.load(src))
words = z['angle_words'].astype(np.int64).copy(); edges = z['tree_edges']; brot = int(z['brot'])
import fastfit as F                                              # our exact evaluator (decodes the words)
P = F.Problem(src)
U, w = P.U, P.w0                                                 # decoded two-body vectors, weights
R, B, N = U.shape
Lm = np.einsum('rbp,rbq,rbc->rcpq', U, U, w)
g = np.linalg.norm(np.einsum('rbc,rcpq->rbpq', w, Lm).reshape(R, B, -1), axis=2).reshape(-1)


def amplitudes(vecs):
    u = np.abs(vecs).copy(); r = np.zeros((len(vecs), N - 1))
    for t in range(N - 2, -1, -1):
        p, c = edges[t]
        r[:, t] = np.hypot(u[:, p], u[:, c]); u[:, p] = r[:, t]; u[:, c] = 0
    return r


s2_two = (amplitudes(U.reshape(-1, N)) * g[:, None]) ** 2
# one-body rows: eigenvectors of h1'; weight by |eigenvalue| (the one-body residual is linear in it)
coef = np.einsum('rc,rbc->rb', P.W0, w)
hp = P.A - np.einsum('rb,rbp,rbq->pq', coef, U, U)
ev = np.linalg.eigvalsh(hp)
one = np.asarray(P.one, float)
lam_one = np.abs(np.einsum('kp,pq,kq->k', one, hp, one))        # Rayleigh quotient of each loaded row
s2_one = (amplitudes(one) * lam_one[:, None]) ** 2 * float(os.environ.get('ONE_WEIGHT', '1'))
s2 = np.vstack([s2_one, s2_two])
n = len(words); H = -(-n // 16); pad = H * 16 - n
wp = np.vstack([words, np.zeros((pad, 75), np.int64)]); sp = np.vstack([s2, np.zeros((pad, 75))])
pts = np.array([[(x >> i) & 1 for i in range(4)] for x in range(16)])
RM1 = np.array([(c[0] + pts @ np.array(c[1:])) % 2 for c in itertools.product((0, 1), repeat=5)], np.uint8)
cols = [c for c in range(75) if b0 <= c // 5 <= b1]
flips = 0; cost = 0.0
for c in cols:
    bits = ((wp[:, c] >> 1) & 1).reshape(H, 16); wt = sp[:, c].reshape(H, 16)
    diff = bits[:, None, :] ^ RM1[None, :, :]
    best = (diff * wt[:, None, :]).sum(2).argmin(1)
    d = diff[np.arange(H), best].reshape(-1)
    flips += int(d[:n].sum()); cost += float((d * sp[:, c]).sum())
    wp[:, c] ^= (d.astype(np.int64) << 1)
new = wp[:n]
assert not np.any(((new >> (brot - 1)) ^ (new >> (brot - 2))) & 1), 'signed gauge broken'
z['angle_words'] = new.astype(z['angle_words'].dtype)
z['coded2'] = np.array(cols, np.int64)
np.savez(out, **z)
print(f'bit 1 coded on {len(cols)} columns (batches {b0}..{b1}): {flips} words changed of {n * len(cols)}, '
      f'weighted cost {cost / (s2_two.sum() / 12):.2f} x rounding variance (two-body scale)')
