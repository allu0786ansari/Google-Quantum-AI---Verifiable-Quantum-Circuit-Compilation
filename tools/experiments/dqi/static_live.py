import sys, numpy as np, collections
path = sys.argv[1]
first, last = {}, {}
starts, ends = [], []
def close(q):
    if q in first:
        starts.append(first.pop(q)); ends.append(last.pop(q))
with open(path) as f:
    for i, line in enumerate(f):
        p = line.split()
        if not p or p[0] in ('INFO_TAG', 'REGISTER', 'APPEND_TO_REGISTER'): continue
        op = p[0]
        toks = p[1:p.index('if')] if 'if' in p else p[1:]
        qs = [int(t[1:]) for t in toks if t[0] == 'q']
        if op == 'R':
            for q in qs: close(q)
            continue
        for q in qs:
            if q not in first: first[q] = i
            last[q] = i
        if op == 'HMR':
            for q in qs: close(q)
for q in list(first): close(q)
s, e = np.array(starts), np.array(ends)
ev = np.concatenate([np.stack([s, np.ones_like(s)], 1), np.stack([e + 1, -np.ones_like(e)], 1)])
ev = ev[np.lexsort((ev[:, 1], ev[:, 0]))]
live = np.cumsum(ev[:, 1])
k = int(np.argmax(live))
print('intervals', len(s), 'conservative peak live', int(live.max()), 'at line', int(ev[k, 0]))
