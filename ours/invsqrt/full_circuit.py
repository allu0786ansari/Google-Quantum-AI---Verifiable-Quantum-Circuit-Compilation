"""Our own Benchmark I circuit, layers I-1..I-5 (polynomial route), plus an all-x test.

I-1 normalise x -> x_norm, s (binary, 4 flags)              layer1.normalise
I-2/3 Y = 2^29 * y ~ 2^29 / sqrt(u) by truncated Horner      layer_poly.horner
I-4 R = Y << s in a 39-bit register (CSWAP stages 1,2,4,8; small first, only data-carrying pairs);
    phase Y * 2^(s - 39) turns = sum_k R_k 2^(k - 39): one Z_POW per bit (k >= 37 Clifford, 36 = T)
I-5 unshift, unhorner (row reversal + measurement-based lookup deletes), denormalise (MBU flags)

    python full_circuit.py [out.kmb]
"""
import sys, time
from fractions import Fraction
import numpy as np
import kickmix as km

sys.path.insert(0, '/home/claude/google-quantum-ai---verifiable-quantum-circuit-compilation/ours/invsqrt')
sys.path.insert(0, '/home/claude/data/comp/verifiers')
from layer1 import normalise, denormalise
from layer_poly import horner, unhorner, WY

RW = 39
# (flag index j, shift d, highest data index before this stage)
SHIFT_STAGES = [(0, 1, 28), (1, 2, 29), (2, 4, 31), (3, 8, 30)]


def shift(b, s, R, inverse=False):
    stages = SHIFT_STAGES if not inverse else list(reversed(SHIFT_STAGES))
    for j, d, top in stages:
        idxs = range(top, -1, -1) if not inverse else range(0, top + 1)
        for i in idxs:
            b.cswap(s[j], R[i], R[i + d])


def phase(b, R, delta=None):
    for k in range(RW):
        t = Fraction(1, 1 << (38 - k)) if k <= 38 else 0
        if delta is not None:
            t += delta[k]
        b.z_pow(R[k], t)


def build(max_qubits=None, delta=None):
    b = km.CircuitBuilder()
    x = b.create_quantum_register(22, 'x')
    if max_qubits:
        b.set_max_qubits(max_qubits)
    s = normalise(b, x)
    y, a1, a2 = horner(b, x)
    top = b.alloc_qubits(1)[0]; b.x(top)                    # the constant bit 28 of Y
    ext = b.alloc_qubits(RW - WY - 1)
    R = list(y) + [top] + list(ext)
    shift(b, s, R)
    phase(b, R, delta)
    shift(b, s, R, inverse=True)
    b.free(ext); b.x(top); b.free(top)
    unhorner(b, x, y, a1, a2)
    denormalise(b, x, s)
    return b.finish_circuit()


def test_all_x(circ, seed=7):
    from invsqrt import verify_submission_invsqrt as V
    spec = V.ProblemSpecInvSqrt.for_n_p(10)
    MAXX = spec.max_x
    target, tol = V.reference_values(spec, np.arange(1, MAXX + 1, dtype=np.uint64))
    B = 1 << 16
    sim = km.Simulator(batch_size=B, debug_print=False, seed=seed, counting_operations=True)
    nq = circ.num_qubits
    xq = km.array([km.q(i) for i in range(22)]); sq = km.array([km.q(i) for i in range(22, nq)])
    worst = 0.0; min_margin = 1 << 64; fails = 0; tof = 0; n = 0
    err_all = np.zeros(MAXX, np.int64)
    for start in range(0, MAXX, B):
        xs = np.arange(start + 1, start + B + 1, dtype=np.uint64); nn = min(B, MAXX - start)
        xs[nn:] = 1
        sim.clear_for_shot(); sim.clear_operation_statistics()
        sim.write_across_shots(xq, V.to_bits(xs, 22))
        sim.do(circ)
        ok = V.from_bits(sim.read_across_shots(xq)) == xs
        ok &= ~sim.read_across_shots(sq).any(axis=0)
        theta = np.array([V.half_turns_to_units(sim.read_shot_phase(i)) for i in range(B)], dtype=np.uint64)
        e = (theta - target[xs - 1]).view(np.int64)
        used = np.abs(e).view(np.uint64); allowed = tol[xs - 1]
        ok &= used <= allowed
        fails += int((~ok[:nn]).sum())
        worst = max(worst, float((used[:nn] / allowed[:nn]).max()))
        min_margin = min(min_margin, int((allowed[:nn].astype(np.int64) - used[:nn].astype(np.int64)).min()))
        err_all[start:start + nn] = e[:nn]
        tof += sim.operation_statistics().num_toffoli_like; n += B
    return dict(fails=fails, worst=worst, min_margin=min_margin, toffoli=tof / n, qubits=nq, err=err_all)


def cost(circ, res, nrot=36):
    import math
    eps = math.pi * res['min_margin'] / 2 ** 64 / nrot     # verifier: n_rot * 2 eps <= 2 pi margin
    tr = 0.56 * math.log2(1 / eps) + 5.3
    tbase = 1                                              # one T gate (k = 36)
    C = res['toffoli'] + (tbase + nrot * tr) / 4
    return C, res['qubits'] * C, 2076774 / (res['qubits'] * C)


if __name__ == '__main__':
    out = sys.argv[1] if len(sys.argv) > 1 else 'ours_invsqrt.kmb'
    t0 = time.time()
    circ = build()
    print(f'built: qubits {circ.num_qubits}, static Toffoli-like {circ.operation_statistics().num_toffoli_like} ({time.time() - t0:.0f}s)', flush=True)
    res = test_all_x(circ)
    C, QC, ratio = cost(circ, res)
    print(f"all x: failures {res['fails']}, worst |err|/tol {res['worst']:.4f}, mean Toffolis {res['toffoli']:.1f}")
    print(f'estimated C {C:.1f}, Q*C {QC:,.0f}, ratio {ratio:.3f}  ({time.time() - t0:.0f}s)')
    circ.to_file(out, format='kmb')
