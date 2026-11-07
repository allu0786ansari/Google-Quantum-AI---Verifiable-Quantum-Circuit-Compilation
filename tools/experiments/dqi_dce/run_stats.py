"""Draw verifier-style shots (optionally at one error weight), run sim, return per-op stats.
    python run_stats.py OPS_PREFIX NQ SHOTS SEED OUT.npy [WEIGHT]"""
import json, os, subprocess, sys, numpy as np
sys.path.insert(0, '/home/claude/google-quantum-ai---verifiable-quantum-circuit-compilation/tools')
os.environ.setdefault('VERIFIERS', '/home/claude/data/comp/verifiers')
import harness as H
import kickmix as km
prefix, nq, shots, seed, out = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), sys.argv[5]
weight = int(sys.argv[6]) if len(sys.argv) > 6 else None
assert shots % 64 == 0
W = shots // 64
meta = json.load(open(prefix + '.json')); nops = meta['n_ops']
nb = int(meta.get('nb', 0))
if not nb:
    ops = np.memmap(prefix + '.bin', dtype=np.int32, mode='r').reshape(-1, 5)
    bitcols = []
    for code, cols in ((7, [2]), (8, [1]), (9, [1]), (10, [1])):
        sel = ops[ops[:, 0] == code]
        bitcols.append(sel[:, cols].max())
    nb = int(max(max(bitcols), ops[:, 4].max())) + 1
    meta['nb'] = nb; json.dump(meta, open(prefix + '.json', 'w'))
V = H._load_dqi_verifier()
spec = V.ProblemSpecDQI.benchmark_instance(num_shots=shots)
probs = V.optimal_dqi_weights(spec.m, spec.max_errors, spec.degree, spec.r) ** 2
if weight is not None:
    probs = np.zeros_like(probs); probs[weight] = 1.0
tables = V.FieldTables(km.GF2Field(spec.degree))
rng = np.random.default_rng(seed)
syn, mbuc, phases = V.sample_shots(spec, tables, rng, shots, probs)
sb = V.to_bits(syn, spec.degree).astype(np.uint8)        # (840, shots)
mb = V.to_bits(mbuc, spec.degree).astype(np.uint8)       # (49140, shots)
def pack(bits):                                          # (n, shots) -> (n, W) uint64
    return np.packbits(bits.reshape(bits.shape[0], W, 64)[:, :, ::-1], axis=2).view('>u8').astype(np.uint64).reshape(bits.shape[0], W)
Q = np.zeros((nq, W), np.uint64); Q[:840] = pack(sb)
B = np.zeros((nb, W), np.uint64); B[:49140] = pack(mb)
tag = f'{os.getpid()}'
Q.tofile(f'q_{tag}.bin'); B.tofile(f'b_{tag}.bin')
r = subprocess.run(['./sim', prefix + '.bin', str(nops), str(nq), str(nb), str(W), f'q_{tag}.bin', f'b_{tag}.bin',
                    str(seed * 7919 + 1), f's_{tag}.bin'], capture_output=True, text=True)
print(r.stderr.strip())
fin = np.fromfile('final_q.bin', dtype=np.uint64).reshape(840, W)
print('syndromes restored:', bool(np.array_equal(fin, Q[:840])))
st = np.fromfile(f's_{tag}.bin', dtype=np.uint32).reshape(nops, 4)
idx = np.flatnonzero(st[:, 0])
np.save(out, np.concatenate([idx[:, None].astype(np.uint32), st[idx]], 1))
for f in (f'q_{tag}.bin', f'b_{tag}.bin', f's_{tag}.bin'):
    os.remove(f)
