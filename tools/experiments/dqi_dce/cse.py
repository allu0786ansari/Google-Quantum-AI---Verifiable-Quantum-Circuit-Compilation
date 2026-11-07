"""Repeated AND products: CCX/CCZ on a control pair (a, b) whose values are unchanged since an earlier
CCX/CCZ on the same pair. Weighted by executions per shot from the simulator stats."""
import numpy as np, collections, sys
prefix, stats = sys.argv[1], sys.argv[2]
ops = np.fromfile(prefix + '.bin', dtype=np.int32).reshape(-1, 5)
st = np.load(stats); act = np.zeros(len(ops)); act[st[:, 0].astype(np.int64)] = st[:, 1] / 4096.0
ver = np.zeros(1200, np.int64); bump = 0
seen = {}; rep_ops = 0; rep_T = 0.0; gaps = []
WRITE1 = {1: [2], 2: [3], 3: [1], 6: [1], 7: [1]}           # op -> written arg positions (1-based)
for i in range(len(ops)):
    op, a, b, c, cd = ops[i]
    if op in (2, 15):
        key = (min(a, b), max(a, b), ver[min(a, b)], ver[max(a, b)], cd)
        if key in seen:
            rep_ops += 1; rep_T += act[i]; gaps.append(i - seen[key])
        seen[key] = i
    if op in WRITE1:
        for pos in WRITE1[op]:
            q = ops[i, pos]; bump += 1; ver[q] = bump
    elif op == 4:
        bump += 1; ver[a] = bump; bump += 1; ver[b] = bump
    elif op == 5:
        bump += 1; ver[b] = bump; bump += 1; ver[c] = bump
gaps = np.array(gaps)
print(f'repeated products: {rep_ops} ops, {rep_T:.1f} Toffolis/shot ({100 * rep_T / act[np.isin(ops[:, 0], [2, 5, 15])].sum():.2f}% of T)')
if len(gaps):
    print('distance (ops) between repeats: median', int(np.median(gaps)), 'p10', int(np.percentile(gaps, 10)), 'p90', int(np.percentile(gaps, 90)))
