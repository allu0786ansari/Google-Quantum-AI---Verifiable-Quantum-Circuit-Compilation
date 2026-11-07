"""kmx text -> ops.bin (int32 x5 per op: op, a, b, c, cond) + tags.json. Streams; low memory."""
import sys, json, numpy as np
OPS = {'CX': 1, 'CCX': 2, 'X': 3, 'SWAP': 4, 'CSWAP': 5, 'R': 6, 'HMR': 7, 'BIT_INVERT': 8,
       'BIT_STORE0': 9, 'BIT_STORE1': 10, 'PUSH_CONDITION': 11, 'POP_CONDITION': 12,
       'Z': 13, 'CZ': 14, 'CCZ': 15, 'NEG': 16, 'Z_POW': 17}
src, out = sys.argv[1], sys.argv[2]
buf = np.empty((1 << 20, 5), np.int32); k = 0; n_ops = 0; tags = []
fo = open(out + '.bin', 'wb')
with open(src) as f:
    for line in f:
        p = line.split()
        if not p or p[0] in ('REGISTER', 'APPEND_TO_REGISTER'):
            continue
        if p[0] == 'INFO_TAG':
            tags.append((n_ops, line.split(None, 1)[1].strip().strip('"'))); continue
        cond = -1
        if 'if' in p:
            j = p.index('if'); cond = int(p[j + 1][1:]); p = p[:j]
        args = [int(t[1:]) for t in (p[1:2] if p[0] == 'Z_POW' else p[1:])]
        assert len(args) <= 3, line
        buf[k, 0] = OPS[p[0]]; buf[k, 4] = cond
        buf[k, 1:4] = -1
        buf[k, 1:1 + len(args)] = args
        k += 1; n_ops += 1
        if k == len(buf):
            fo.write(buf.tobytes()); k = 0
fo.write(buf[:k].tobytes()); fo.close()
json.dump(dict(n_ops=n_ops, tags=tags), open(out + '.json', 'w'))
print(n_ops, 'ops', len(tags), 'tags')
