"""Accumulate per-op firing over natural shots and per-weight strata.
    python sweep.py STATS_PREFIX NAT_RUNS"""
import os, subprocess, sys, numpy as np
prefix, nat_runs = sys.argv[1], int(sys.argv[2])
cand = np.load('cand_stats_nat_4096.npy', allow_pickle=True).item()
ids = np.unique(np.concatenate([cand['dead'], cand['a_one']]))
res = {}
def run(seed, weight=None):
    out = f'tmp_{seed}_{weight}.npy'
    cmd = ['python3', 'run_stats.py', 'ops893', '893', '4096', str(seed), out] + ([str(weight)] if weight is not None else [])
    r = subprocess.run(cmd, capture_output=True, text=True)
    assert 'syndromes restored: True' in r.stdout and 'ancilla bits at end: 0' in r.stdout, r.stdout + r.stderr
    st = np.load(out); os.remove(out)
    m = dict(zip(st[:, 0].tolist(), range(len(st))))
    rows = np.array([m.get(int(i), -1) for i in ids])
    fire = np.where(rows >= 0, st[np.maximum(rows, 0), 4], 0)
    a1 = np.where(rows >= 0, st[np.maximum(rows, 0), 2], 0)
    b1 = np.where(rows >= 0, st[np.maximum(rows, 0), 3], 0)
    act = np.where(rows >= 0, st[np.maximum(rows, 0), 1], 0)
    return fire, a1, b1, act
acc = {}
for s in range(nat_runs):
    f, a1, b1, act = run(100 + s)
    for k, v in (('fire', f), ('a1', a1), ('b1', b1), ('act', act)):
        acc[('nat', k)] = acc.get(('nat', k), 0) + v
    print('nat run', s, 'cumulative candidates fired:', int((acc[('nat', 'fire')] > 0).sum()), flush=True)
for w in range(14, 36):
    f, a1, b1, act = run(1000 + w, w)
    acc[(w, 'fire')], acc[(w, 'a1')], acc[(w, 'b1')], acc[(w, 'act')] = f, a1, b1, act
    print('weight', w, 'candidates fired:', int((f > 0).sum()), flush=True)
np.save(prefix + '_sweep.npy', dict(ids=ids, acc=acc), allow_pickle=True)
