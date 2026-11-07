"""Exact (sampled) qubit occupancy of a dqi circuit: run it in chunks on verifier-style
inputs and record, after each chunk, which qubits are nonzero on any shot.

    VERIFIERS=... python chunk_live.py circuit.kmb circuit.kmx OUT.npz [CHUNK_LINES] [SHOTS]
"""
import os
import sys
import time

import numpy as np

sys.path.insert(0, '/home/claude/google-quantum-ai---verifiable-quantum-circuit-compilation/tools')
import harness as H
import kickmix as km

kmb, kmx, out = sys.argv[1:4]
chunk = int(sys.argv[4]) if len(sys.argv) > 4 else 20000
shots = int(sys.argv[5]) if len(sys.argv) > 5 else 256
V = H._load_dqi_verifier()
spec = V.ProblemSpecDQI.benchmark_instance(num_shots=shots)
full = km.Circuit.from_file(kmb)
syn_q, mbuc_b, anc_q = V.validate_registers(full, spec)
Q = full.num_qubits
weights = V.optimal_dqi_weights(spec.m, spec.max_errors, spec.degree, spec.r) ** 2
tables = V.FieldTables(km.GF2Field(spec.degree))
rng = np.random.default_rng(12345)
syndromes, mbuc, phases = V.sample_shots(spec, tables, rng, shots, weights)
sim = km.Simulator(batch_size=shots, debug_print=False, seed=777)
sim.use_same_registers_as(full)
sim.write_across_shots(km.array(syn_q), V.to_bits(syndromes, spec.degree))
sim.write_across_shots(km.array(mbuc_b), V.to_bits(mbuc, spec.degree))
allq = km.array([km.q(i) for i in range(Q)])

rows, line_at, stage_at = [], [], []
buf, depth, n = [], 0, 0
stack = []
t0 = time.time()


def flush():
    if buf:
        sim.do(km.Circuit('\n'.join(buf)))
        buf.clear()
    rows.append(np.packbits(sim.read_across_shots(allq).any(axis=1)))
    line_at.append(n)
    stage_at.append('/'.join(stack[:2]))


with open(kmx) as f:
    for n, line in enumerate(f, 1):
        op = line.split(None, 1)[0] if line.strip() else ''
        if op in ('REGISTER', 'APPEND_TO_REGISTER', ''):
            continue
        if op == 'INFO_TAG':
            tag = line.split(None, 1)[1].strip().strip('"')
            if tag.startswith('@enter:'): stack.append(tag[7:])
            elif tag.startswith('@exit:'): stack.pop()
        buf.append(line.rstrip('\n'))
        if op == 'PUSH_CONDITION':
            depth += 1
        elif op == 'POP_CONDITION':
            depth -= 1
        if len(buf) >= chunk and depth == 0:
            flush()
            if len(rows) % 200 == 0:
                print(len(rows), n, f'{time.time() - t0:.0f}s', flush=True)
flush()
nz = np.unpackbits(np.array(rows), axis=1)[:, :Q].astype(bool)
counts = nz.sum(1)
k = int(counts.argmax())
print('boundaries', len(rows), 'max nonzero qubits', int(counts.max()), 'at line', line_at[k],
      'final nonzero ancillas', int(nz[-1, 840:].sum()), f'{time.time() - t0:.0f}s')
np.savez_compressed(out, nz=np.array(rows), line=np.array(line_at), Q=Q, stage=np.array(stage_at))
