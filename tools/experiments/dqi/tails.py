"""Count Toffolis in XOR-only qubit tails (q only ever an XOR target until R q / end of circuit)."""
import sys, re, collections
path = sys.argv[1]
XOR_ARITY = {'CX': 2, 'CCX': 3, 'X': 1}
protected = int(sys.argv[2]) if len(sys.argv) > 2 else 840     # syndrome qubits q0..q839
tail = collections.defaultdict(list)      # q -> [(line, op)] candidate tail ops
saved = collections.Counter(); seg_ends = collections.Counter(); break_reasons = collections.Counter()
nline = 0
def qnum(t): return int(t[1:])
with open(path) as f:
    for nline, line in enumerate(f):
        parts = line.split()
        if not parts: continue
        op = parts[0]
        if op in ('INFO_TAG', 'REGISTER', 'APPEND_TO_REGISTER'): continue
        cond = 'if' in parts
        toks = parts[1:parts.index('if')] if cond else parts[1:]
        qs = [qnum(t) for t in toks if t[0] == 'q']
        if op == 'R':
            for q in qs:
                if tail[q]:
                    for _, o in tail[q]: saved[o] += 1
                    seg_ends['R'] += 1
                tail[q] = []
            continue
        if op == 'HMR':
            for q in qs:
                tail[q] = []; seg_ends['HMR'] += 1
            continue
        if op in XOR_ARITY and not cond:
            k = XOR_ARITY[op]
            for i in range(0, len(qs), k):
                grp = qs[i:i + k]
                *ctrl, t = grp
                for c in ctrl:
                    if tail[c]: break_reasons['read_as_control'] += 1
                    tail[c] = []
                if t >= protected:
                    tail[t].append((nline, op))
            continue
        # anything else touching q (CZ, Z, SWAP, CSWAP, CCZ, conditioned gates, ...) reads q
        for q in qs:
            if tail[q]: break_reasons[op + ('?' if cond else '')] += 1
            tail[q] = []
for q, lst in tail.items():
    if lst:
        for _, o in lst: saved[o] += 1
        seg_ends['end'] += 1
print('convertible tail ops:', dict(saved), 'segments ending:', dict(seg_ends))
print('tail breaks:', break_reasons.most_common(8), 'lines', nline)
