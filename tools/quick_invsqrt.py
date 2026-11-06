"""Quick invsqrt check on the first MAX_X inputs (default 65,536), 4 shots each.

Uses the official verifier. The score it prints is only indicative: the
rotation cost depends on the worst margin over *all* x, so run the full
verifier before submitting.

    VERIFIERS=/path/to/verifiers python -I tools/quick_invsqrt.py circuit.kmb [MAX_X]
"""
import os
import sys
import time

import kickmix as km

sys.path.insert(0, os.environ.get('VERIFIERS', 'verifiers'))
from invsqrt import verify_submission_invsqrt as V  # noqa: E402

circuit = km.Circuit.from_file(sys.argv[1])
max_x = int(sys.argv[2]) if len(sys.argv) > 2 else 65536
spec = V.ProblemSpecInvSqrt(width=22, max_x=max_x, beta=V.Fraction(1), precision_bits=24,
                            shots_per_x=4, system_qubits=0)
t = time.time()
r = V.verify(circuit, spec, batch_size=256, max_workers=os.cpu_count())
keys = ['score', 'toffoli_equivalent', 'toffoli_count', 't_count', 'rotation_count',
        'rotation_eps', 'num_qubits', 'worst_error_over_tolerance', 'worst_x', 'verified_shots']
print({k: r[k] for k in keys}, f'{time.time() - t:.1f}s')
