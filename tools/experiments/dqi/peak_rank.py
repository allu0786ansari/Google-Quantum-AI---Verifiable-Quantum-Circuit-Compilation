"""GF(2) rank (with a constant row: affine dependencies) of the full qubit state at given kmx lines."""
import os, sys, time
import numpy as np
sys.path.insert(0, '/home/claude/google-quantum-ai---verifiable-quantum-circuit-compilation/tools')
import harness as H
import kickmix as km
kmb, kmx = sys.argv[1], sys.argv[2]; targets = sorted(int(x) for x in sys.argv[3:])
SH = 4096
V = H._load_dqi_verifier(); spec = V.ProblemSpecDQI.benchmark_instance(num_shots=SH)
full = km.Circuit.from_file(kmb); syn_q, mbuc_b, _ = V.validate_registers(full, spec); Q = full.num_qubits
w = V.optimal_dqi_weights(spec.m, spec.max_errors, spec.degree, spec.r) ** 2
syn, mbuc, _ = V.sample_shots(spec, V.FieldTables(km.GF2Field(spec.degree)), np.random.default_rng(7), SH, w)
sim = km.Simulator(batch_size=SH, debug_print=False, seed=99); sim.use_same_registers_as(full)
sim.write_across_shots(km.array(syn_q), V.to_bits(syn, spec.degree)); sim.write_across_shots(km.array(mbuc_b), V.to_bits(mbuc, spec.degree))
allq = km.array([km.q(i) for i in range(Q)])


def rank_gf2(M):
    rows = [int.from_bytes(np.packbits(r).tobytes(), 'big') for r in M]
    basis = {}
    for r in rows:
        while r:
            h = r.bit_length() - 1
            if h in basis: r ^= basis[h]
            else: basis[h] = r; break
    return len(basis)


buf, depth, n, ti = [], 0, 0, 0
with open(kmx) as f:
    for n, line in enumerate(f, 1):
        op = line.split(None, 1)[0] if line.strip() else ''
        if op not in ('REGISTER', 'APPEND_TO_REGISTER', ''):
            buf.append(line.rstrip('\n'))
            depth += (op == 'PUSH_CONDITION') - (op == 'POP_CONDITION')
        if len(buf) >= 200000 and depth == 0:
            sim.do(km.Circuit('\n'.join(buf))); buf = []
        if ti < len(targets) and n >= targets[ti] and depth == 0:
            if buf: sim.do(km.Circuit('\n'.join(buf))); buf = []
            S = sim.read_across_shots(allq).astype(np.uint8)          # (Q, SH)
            live = S.any(axis=1); const1 = S.all(axis=1)
            Mv = S[live & ~const1]
            r = rank_gf2(Mv); ra = rank_gf2(np.vstack([Mv, np.ones((1, SH), np.uint8)]))
            print(f'line {n}: live {int(live.sum())} (always-1 {int(const1.sum())}); GF(2) rank {r}, '
                  f'affine rank {ra - 1}; linear dependencies {Mv.shape[0] - r}', flush=True)
            ti += 1
            if ti == len(targets): break
