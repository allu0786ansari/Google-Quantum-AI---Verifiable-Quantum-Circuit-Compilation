"""All x in [1, max_x], bit-sliced, several HMR seeds; per-op firing stats summed over seeds."""
import json, subprocess, sys, numpy as np
MAXX = 3133452; seeds = [int(s) for s in sys.argv[1:]] or [1, 2, 3, 4]
meta = json.load(open('opsinv.json')); nops = meta['n_ops']
ops = np.fromfile('opsinv.bin', dtype=np.int32).reshape(-1, 5)
nb = int(max(ops[np.isin(ops[:, 0], [8, 9, 10]), 1].max(initial=0), ops[ops[:, 0] == 7, 2].max(initial=0), ops[:, 4].max())) + 1
nq = 62; W = (MAXX + 63) // 64
xs = np.arange(1, W * 64 + 1); xs[xs > MAXX] = 1
bits = ((xs[None, :] >> np.arange(22)[:, None]) & 1).astype(np.uint8)
packed = np.packbits(bits.reshape(22, W, 64)[:, :, ::-1], axis=2).view('>u8').astype(np.uint64).reshape(22, W)
Q = np.zeros((nq, W), np.uint64); Q[:22] = packed; Q.tofile('q0.bin')
np.zeros((nb, W), np.uint64).tofile('b0.bin')
total = None
for s in seeds:
    r = subprocess.run(['/home/claude/frontier/dce/sim', 'opsinv.bin', str(nops), str(nq), str(nb), str(W), 'q0.bin', 'b0.bin', str(s), 'st.bin'],
                       capture_output=True, text=True)
    st = np.fromfile('st.bin', dtype=np.uint32).reshape(nops, 4).astype(np.int64)
    fin = np.fromfile('final_q.bin', dtype=np.uint64)  # first 840 rows dumped by sim; only 62 exist -> read what's there
    print('seed', s, r.stderr.strip(), flush=True)
    total = st if total is None else total + st
np.save('inv_stats.npy', total)
code = ops[:, 0]
tof = np.isin(code, [2, 5, 15]) & (total[:, 0] > 0)
shots = W * 64 * len(seeds)
print(f'mean Toffoli-like per shot (sim): {total[tof, 0].sum() / shots:.4f}')
dead = tof & (total[:, 3] == 0)
print(f'never fire over all x and {len(seeds)} seeds: {dead.sum()} ops, {total[dead, 0].sum() / shots:.4f} Toffolis/shot')
ccx = np.isin(code, [2, 15]) & tof & ~dead
one = ccx & ((total[:, 1] == total[:, 0]) | (total[:, 2] == total[:, 0]))
print(f'one control always 1 (-> CX/CZ): {one.sum()} ops, {total[one, 0].sum() / shots:.4f} Toffolis/shot')
np.save('inv_cands.npy', dict(dead=np.flatnonzero(dead), one=np.flatnonzero(one)), allow_pickle=True)
