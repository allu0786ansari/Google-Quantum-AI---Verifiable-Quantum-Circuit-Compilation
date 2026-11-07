"""Layer I-1: in-place normalisation x -> x_norm = x << 2s with s in binary (4 qubits).

Stage p in (8, 4, 2, 1) pairs, d = 2p bits: flag_p = AND(NOT top d bits) via an AND tree whose
temporaries are cleared by measurement (del_and), then x <<= d controlled by flag_p, using only the
22 - d CSWAPs that move data (the vacated top bits are zero when the flag is set).
"""
import kickmix as km

STAGES = (8, 4, 2, 1)


def _and_tree(b, bits):
    """Return (result, temps) with result = AND(bits); temps to be cleared in reverse."""
    layer, temps = list(bits), []
    while len(layer) > 1:
        nxt = []
        for i in range(0, len(layer) - 1, 2):
            t = b.alloc_qubits(1)[0]; b.init_and(layer[i], layer[i + 1], t)
            temps.append((layer[i], layer[i + 1], t)); nxt.append(t)
        if len(layer) % 2: nxt.append(layer[-1])
        layer = nxt
    return layer[0], temps


def _clear_tree(b, temps, keep):
    for a, c, t in reversed(temps):
        if t is keep: continue
        b.del_and(a, c, t); b.free(t)


def normalise(b, x):
    """x: 22 qubits little-endian. Returns s (4 qubits, s[j] = stage with 2^j pairs)."""
    s = [None] * 4
    for j, pairs in zip((3, 2, 1, 0), STAGES):
        d = 2 * pairs
        top = [x[i] for i in range(22 - d, 22)]
        for q in top: b.x(q)
        flag, temps = _and_tree(b, top)                    # the tree root is the flag itself
        _clear_tree(b, temps, flag)
        for q in top: b.x(q)
        for i in range(21 - d, -1, -1):
            b.cswap(flag, x[i], x[i + d])
        s[j] = flag
    return s


def denormalise(b, x, s):
    """Exact inverse of normalise: unshift, then erase each flag by measurement with a phase fix-up."""
    for j, pairs in zip((0, 1, 2, 3), reversed(STAGES)):
        d = 2 * pairs; flag = s[j]
        for i in range(0, 22 - d):
            b.cswap(flag, x[i], x[i + d])
        top = [x[i] for i in range(22 - d, 22)]
        for q in top: b.x(q)
        m = b.alloc_bits(1)[0]; b.hmr(flag, m); b.free(flag)
        # phase fix-up (-1)^(m * AND(NOT top)) : AND tree, then CZ-like with the measured bit
        b.push_condition(m)                                # fix-up only runs when m = 1 (half the shots)
        res, temps = _and_tree(b, top)
        b.z(res)
        _clear_tree(b, temps, None)
        b.pop_condition()
        for q in top: b.x(q)
        b.free(m)


if __name__ == '__main__':
    import sys, time
    import numpy as np
    sys.path.insert(0, '/home/claude/data/comp/verifiers')
    from invsqrt import verify_submission_invsqrt as V
    MAXX = 3133452
    # forward-only test circuit: check x_norm and s on every x
    b = km.CircuitBuilder(); xr = b.create_quantum_register(22, 'x')
    s = normalise(b, xr)
    sq = list(s)
    fwd = b.finish_circuit()
    st = fwd.operation_statistics()
    print('forward: qubits', fwd.num_qubits, 'Toffoli-like', st.num_toffoli_like)
    sim = km.Simulator(batch_size=1 << 16, debug_print=False, seed=1)
    bad = 0; t0 = time.time()
    xs_all = np.arange(1, MAXX + 1, dtype=np.uint64)
    for start in range(0, MAXX, 1 << 16):
        xs = xs_all[start:start + (1 << 16)]
        if len(xs) < (1 << 16): xs = np.concatenate([xs, np.ones((1 << 16) - len(xs), np.uint64)])
        sim.clear_for_shot()
        sim.write_across_shots(km.array(list(xr)), V.to_bits(xs, 22))
        sim.do(fwd)
        got_x = V.from_bits(sim.read_across_shots(km.array(list(xr))))
        sbits = sim.read_across_shots(km.array(sq)).astype(np.int64)
        got_s = (sbits * (1 << np.arange(4))[:, None]).sum(0)
        v = xs.astype(np.int64); ss = np.zeros_like(v)
        for pairs in (8, 4, 2, 1):
            dd = 2 * pairs; z = (v >> (22 - dd)) == 0; v = np.where(z, v << dd, v); ss += np.where(z, pairs, 0)
        bad += int(((got_x.astype(np.int64) != v) | (got_s != ss)).sum())
    print(f'forward check on all x: {bad} mismatches ({time.time() - t0:.0f}s)')
    # round trip: normalise then denormalise must restore x and clear everything
    b = km.CircuitBuilder(); xr = b.create_quantum_register(22, 'x')
    s = normalise(b, xr); denormalise(b, xr, s)
    rt = b.finish_circuit()
    print('round trip: qubits', rt.num_qubits, 'Toffoli-like (static)', rt.operation_statistics().num_toffoli_like)
    rt.to_file('layer1_roundtrip.kmb', format='kmb')
    sim = km.Simulator(batch_size=1 << 16, debug_print=False, seed=2, counting_operations=True)
    allq = km.array([km.q(i) for i in range(rt.num_qubits)])
    bad = 0; tof = 0; n = 0
    for start in range(0, MAXX, 1 << 16):
        xs = xs_all[start:start + (1 << 16)]
        if len(xs) < (1 << 16): xs = np.concatenate([xs, np.ones((1 << 16) - len(xs), np.uint64)])
        sim.clear_for_shot(); sim.clear_operation_statistics()
        sim.write_across_shots(km.array(list(xr)), V.to_bits(xs, 22))
        sim.do(rt)
        state = sim.read_across_shots(allq)
        restored = V.from_bits(state[:22]) == xs
        clean = ~state[22:].any(axis=0)
        phase_ok = ~sim.read_phase_flipped_across_shots()
        bad += int((~(restored & clean & phase_ok)).sum())
        tof += sim.operation_statistics().num_toffoli_like; n += 1 << 16
    print(f'round trip on all x: {bad} failures; mean Toffolis per shot {tof / n:.2f}')
