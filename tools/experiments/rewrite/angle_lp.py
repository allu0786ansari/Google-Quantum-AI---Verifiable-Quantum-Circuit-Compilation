"""Retune invsqrt's non-Clifford Z_POW angles by a minimax LP over every x; price the result.

The phase is linear in each rotation exponent t_r (half turns): changing t_r by d_r adds
d_r / 2 * b_r(x) turns, b_r(x) = value of the rotated qubit at that gate. We minimise
max_x |e(x) + sum_r d_r b_r(x) / 2| / tol(x) by cutting planes, then recompute the verifier's
rotation cost from the new margin.
    python angle_lp.py circuit.kmb
"""
import math, os, sys, time
from fractions import Fraction
import numpy as np
from scipy.optimize import linprog
sys.path.insert(0, '/home/claude/data/comp/verifiers')
import kickmix as km
from invsqrt import verify_submission_invsqrt as V

path = sys.argv[1]
circ = km.Circuit.from_file(path, format='kmb')
text = str(circ).splitlines()
spec = V.ProblemSpecInvSqrt.for_n_p(10)
MAXX = spec.max_x
target, tol = V.reference_values(spec, np.arange(1, MAXX + 1, dtype=np.uint64))

# split points: non-Clifford Z_POW lines (exponent not a multiple of 1/4 half-turn ... i.e. not T/Clifford)
splits, rot = [], []
for i, line in enumerate(text):
    if line.startswith('Z_POW') and ' if ' not in line:
        q, t = line.split()[1], Fraction(line.split()[2])
        if (t * 4).denominator != 1:
            splits.append(i); rot.append((q, t, i))
print('non-Clifford rotations:', len(rot), flush=True)
meta = [l for l in text if l.startswith(('REGISTER', 'APPEND'))]
body = [l for l in text if not l.startswith(('REGISTER', 'APPEND'))]
# map split indices into body indices
body_idx = []; bi = 0
for i, l in enumerate(text):
    if not l.startswith(('REGISTER', 'APPEND')):
        if i in set(splits): body_idx.append(bi)
        bi += 1
chunks, prev = [], 0
for b in body_idx:
    chunks.append(km.Circuit('\n'.join(meta + body[prev:b]))); prev = b
chunks.append(km.Circuit('\n'.join(meta + body[prev:])))
rq = [km.q(int(r[0][1:])) for r in rot]
xq = km.array([km.q(i) for i in range(22)])
BATCH = 1 << 16
nrot = len(rot)
B = np.zeros((nrot, MAXX), dtype=bool)
err = np.zeros(MAXX, dtype=np.int64)
t0 = time.time()
sim = km.Simulator(batch_size=BATCH, debug_print=False, seed=12345)
for start in range(0, MAXX, BATCH):
    xs = np.arange(start + 1, start + BATCH + 1, dtype=np.uint64); valid = xs <= MAXX
    xs = np.where(valid, xs, 1)
    sim.clear_for_shot()
    sim.write_across_shots(xq, V.to_bits(xs, 22))
    for r, ch in enumerate(chunks):
        sim.do(ch)
        if r < nrot:
            B[r, start:start + BATCH][valid[:min(BATCH, MAXX - start)] if False else slice(None)] = \
                sim.read_across_shots(rq[r])[:min(BATCH, MAXX - start)]
    theta = np.array([V.half_turns_to_units(sim.read_shot_phase(i)) for i in range(BATCH)], dtype=np.uint64)
    e = (theta - target[xs - np.uint64(1)]).view(np.int64)
    n = min(BATCH, MAXX - start)
    err[start:start + n] = e[:n]
    if start % (BATCH * 8) == 0:
        print(f'  x up to {start + n:,} ({time.time() - t0:.0f}s)', flush=True)
tolf = tol.astype(np.float64)
ratio0 = np.abs(err) / tolf
print(f'current worst |error|/tol {ratio0.max():.6f} at x={int(ratio0.argmax()) + 1}', flush=True)

# LP in units of 2^-64 turns: maximise m with |e + B d| <= tol - m for all x
Bf = B.astype(np.float64)
slack0 = tolf - np.abs(err)
act = set(np.argsort(slack0)[:6000].tolist())
for it in range(40):
    rows = np.array(sorted(act))
    A1 = np.hstack([Bf[:, rows].T, np.ones((len(rows), 1))]); b1 = tolf[rows] - err[rows]
    A2 = np.hstack([-Bf[:, rows].T, np.ones((len(rows), 1))]); b2 = tolf[rows] + err[rows]
    c = np.zeros(nrot + 1); c[-1] = -1
    res = linprog(c, A_ub=np.vstack([A1, A2]), b_ub=np.concatenate([b1, b2]),
                  bounds=[(None, None)] * nrot + [(None, None)], method='highs')
    d = res.x[:nrot]; m = res.x[-1]
    slack = tolf - np.abs(err + d @ Bf)
    viol = np.flatnonzero(slack < m * (1 - 1e-9))
    print(f'  LP round {it}: margin {m / 2**64:.4e} turns, rows {len(rows)}, full min slack {slack.min() / 2**64:.4e}, violated {len(viol)}', flush=True)
    if len(viol) == 0:
        break
    act |= set(viol[np.argsort(slack[viol])[:6000]].tolist())


def cost(ratio_worst_units):
    return None


# price: margin = min_x (tol - |e'|), eps = pi * margin / 2^64 / nrot, rotation T = 0.56 log2(1/eps) + 5.3
for name, e_new in (('current', err.astype(np.float64)), ('retuned', err + d @ Bf)):
    margin = float(np.min(tolf - np.abs(e_new)))
    eps = math.pi * margin / 2**64 / nrot
    tr = 0.56 * math.log2(1 / eps) + 5.3
    print(f'{name}: worst ratio {np.max(np.abs(e_new) / tolf):.6f}, margin {margin / 2**64:.3e} turns, '
          f'eps {eps:.3e}, T per rotation {tr:.3f}, rotation Toffoli-eq {nrot * tr / 4:.2f}')
new_t = [r[1] + Fraction(int(round(di)), 1 << 63) for r, di in zip(rot, d)]   # d in 2^-64 turns = 2^-63 half turns
np.save('angle_lp.npy', dict(delta_units=d, rot=[(r[0], str(r[1]), r[2]) for r in rot]), allow_pickle=True)
