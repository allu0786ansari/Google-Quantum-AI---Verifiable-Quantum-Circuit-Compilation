"""Quick dqi check on a reduced shot count (default 512).

Runs the official verifier source, skipping only the `reaction_depth()`
report, which is not part of the score and runs out of memory on small
machines. The default log2_delta of -1 loosens the confidence bound so a
short run passes; the leaderboard uses -128 over 102,400 shots.

    VERIFIERS=/path/to/verifiers python -I tools/quick_dqi.py circuit.kmb [NUM_SHOTS] [LOG2_DELTA]
"""
import os
import sys
import time
import types

import kickmix as km

src_path = os.path.join(os.environ.get('VERIFIERS', 'verifiers'), 'dqi', 'verify_submission_dqi.py')
src = open(src_path).read().replace('circuit.reaction_depth()', '-1')
V = types.ModuleType('vdqi')
V.__file__ = src_path
sys.modules['vdqi'] = V
exec(compile(src, src_path, 'exec'), V.__dict__)

circuit = km.Circuit.from_file(sys.argv[1])
n = int(sys.argv[2]) if len(sys.argv) > 2 else 512
delta = float(sys.argv[3]) if len(sys.argv) > 3 else -1.0
spec = V.ProblemSpecDQI.benchmark_instance(num_shots=n)
t = time.time()
r = V.verify(circuit, spec, max_workers=1, log2_delta=delta)
ops = r.pop('shot_operation_counts')
print(r)
print({k: v / n for k, v in ops.items()})
print(f'{time.time() - t:.1f}s')
