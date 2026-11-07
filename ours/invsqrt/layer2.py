"""Layer I-2: seed y0 = C0[idx] - C1[idx] * t / 2^M from two table lookups and M controlled subtractions.

x holds x_norm (after layer I-1). idx = x[22-K:22], t = x[22-K-M:22-K]. y0 has F0 fractional bits.
The C1 register is temporary: cleared by the measurement-based del_lookup right after use.
"""
import numpy as np
import kickmix as km

T = np.load('/home/claude/google-quantum-ai---verifiable-quantum-circuit-compilation/ours/invsqrt/seed_tables.npz')
C0, C1 = [int(v) for v in T['C0']], [int(v) for v in T['C1']]
K, M, F0, F1 = (int(T[k]) for k in ('K', 'M', 'F0', 'F1'))
W0 = 18
W1 = max(v.bit_length() for v in C1)


def seed(b, x):
    addr = km.array([x[i] for i in range(22 - K, 22)])
    y0 = b.alloc_qubits(W0); c1 = b.alloc_qubits(W1)
    b.init_lookup(table=[a | (c << W0) for a, c in zip(C0, C1)], address=addr,
                  target=km.array(list(y0) + list(c1)))          # one lookup loads both coefficients
    _apply_slope(b, x, y0, c1, sign=-1)
    b.del_lookup(table=C1, address=addr, target=c1); b.free(c1)
    return y0


def _apply_slope(b, x, y0, c1, sign):
    js = range(M) if sign < 0 else reversed(range(M))
    for j in js:
        sh = j + F0 - F1 - M
        ctrl = x[22 - K - M + j]
        low = list(y0)[:W0 - 1]          # y0 stays in (1/2, 1): its top bit is always 1 and never borrowed from
        if sh >= 0:
            off, tgt = km.array(list(c1)), km.array(low[sh:])
        else:
            if -sh >= W1: continue
            off, tgt = km.array(list(c1)[-sh:]), km.array(low)
        (b.isub if sign < 0 else b.iadd)(off, target=tgt, control=ctrl)


def unseed(b, x, y0):
    addr = km.array([x[i] for i in range(22 - K, 22)])
    c1 = b.alloc_qubits(W1)
    b.init_lookup(table=C1, address=addr, target=c1)
    _apply_slope(b, x, y0, c1, sign=+1)
    b.del_lookup(table=C1, address=addr, target=c1); b.free(c1)
    b.del_lookup(table=C0, address=addr, target=y0); b.free(y0)


if __name__ == '__main__':
    import sys, time
    sys.path.insert(0, '/home/claude/data/comp/verifiers'); sys.path.insert(0, '.')
    from invsqrt import verify_submission_invsqrt as V
    from layer1 import normalise, denormalise
    from seed_model import seed_int, xn, ytrue
    MAXX = 3133452
    idx = xn >> (22 - K); tt = (xn >> (22 - K - M)) & ((1 << M) - 1)
    want = seed_int(idx, tt, np.array(C0), np.array(C1))
    xs_all = np.arange(1, MAXX + 1, dtype=np.uint64)

    def run(circ, check):
        sim = km.Simulator(batch_size=1 << 16, debug_print=False, seed=3, counting_operations=True)
        bad = tof = n = 0
        for start in range(0, MAXX, 1 << 16):
            xs = xs_all[start:start + (1 << 16)]; nn = len(xs)
            if nn < (1 << 16): xs = np.concatenate([xs, np.ones((1 << 16) - nn, np.uint64)])
            sim.clear_for_shot(); sim.clear_operation_statistics()
            sim.write_across_shots(km.array([km.q(i) for i in range(22)]), V.to_bits(xs, 22))
            sim.do(circ)
            bad += check(sim, start, nn)
            tof += sim.operation_statistics().num_toffoli_like; n += 1 << 16
        return bad, tof / n

    # forward: normalise + seed; compare y0 with the model
    b = km.CircuitBuilder(); xr = b.create_quantum_register(22, 'x')
    s = normalise(b, xr); y0 = seed(b, xr)
    fwd = b.finish_circuit()
    yq = km.array(list(y0))

    def chk(sim, start, nn):
        got = (sim.read_across_shots(yq).astype(np.int64) * (1 << np.arange(W0))[:, None]).sum(0)[:nn]
        return int((got != want[start:start + nn]).sum())
    bad, mt = run(fwd, chk)
    print(f'I-1 + I-2 forward: qubits {fwd.num_qubits}, mean Toffolis {mt:.2f}; y0 mismatches vs model: {bad}')
    # round trip: normalise, seed, unseed, denormalise
    b = km.CircuitBuilder(); xr = b.create_quantum_register(22, 'x')
    s = normalise(b, xr); y0 = seed(b, xr); unseed(b, xr, y0); denormalise(b, xr, s)
    rt = b.finish_circuit()
    nq = rt.num_qubits; allq = km.array([km.q(i) for i in range(nq)])

    def chk2(sim, start, nn):
        st = sim.read_across_shots(allq)
        ok = (V.from_bits(st[:22]) == np.concatenate([xs_all[start:start + nn], np.ones((1 << 16) - nn, np.uint64)])) \
            & ~st[22:].any(axis=0) & ~sim.read_phase_flipped_across_shots()
        return int((~ok[:nn]).sum())
    bad, mt = run(rt, chk2)
    print(f'round trip: qubits {nq}, mean Toffolis {mt:.2f}; failures {bad}')
