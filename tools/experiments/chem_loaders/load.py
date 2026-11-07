import json, sys, numpy as np
sys.path.insert(0, '/home/claude/frontier/v2src')
import tensorqrom as tq
WORDS = np.load('/home/claude/data/kg_assets/qcc-reproducible-circuit-assets-v2/v2-c39-factor.npz')['angle_words']
def slot(k):
    plan = json.load(open(f'/home/claude/frontier/v2src/plans/slot-{k}.json'))
    cur = WORDS[:, k*5:k*5+5]; prev = WORDS[:, (k-1)*5:k*5] if k else np.zeros_like(cur)
    delta = np.vstack([cur ^ prev, np.asarray(plan['unused_address_completion'], dtype=np.int64)])
    C = tq.anf_tensor(delta, 8, 4)
    terms = [tuple(int(v, 16) for v in t) for t in plan['terms']]
    return plan, delta, C, terms
def gf2rank(M):
    M = M.copy() % 2; r = 0; rows, cols = M.shape
    for c in range(cols):
        p = next((i for i in range(r, rows) if M[i, c]), None)
        if p is None: continue
        M[[r, p]] = M[[p, r]]
        for i in range(rows):
            if i != r and M[i, c]: M[i] ^= M[r]
        r += 1
        if r == rows: break
    return r
