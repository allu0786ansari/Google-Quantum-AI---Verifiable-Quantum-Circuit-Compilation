"""Turn a trained factorization into a verified dfthcsos submission (our pipeline).

    python finish.py trained.npz brot refit_epochs refit_lr outdir

1. Round the two-body vectors to brot-bit Givens angles (log-depth tree), freeze them and refit the
   continuous weights (W, H2_I) so they absorb the rounding error (starter-kit trainer, objective eq. 48).
2. One-body rows: eigenvectors of h1' rebuilt from the decoded two-body vectors.
3. Largest J under the 128-qubit cap; one plain QROM per batch.
4. Official verifier.
"""
import sys, time, json, types
from pathlib import Path
import numpy as np
REPO = Path('/home/claude/google-quantum-ai---verifiable-quantum-circuit-compilation')
VER = Path('/home/claude/data/comp/verifiers')
sys.path.insert(0, str(VER))
src, brot, epochs, lr, outdir = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), float(sys.argv[4]), Path(sys.argv[5])
outdir.mkdir(parents=True, exist_ok=True)
ns = dict(np=np, Path=Path, sys=sys, time=time, VERIFIERS=VER, OUT=outdir, SUB=outdir, SAMPLE=None,
          STARTER_KIT='/home/claude/data/starter/starter_kit')
code = open(REPO / 'notebooks/src/dfthcsos.py').read()
code = code.split('def run_dfthcsos')[0]
exec(compile(code, 'dfthcsos_src', 'exec'), ns)
t0 = time.time()
fact = Path(src)
if epochs > 0:
    fact = ns['dft_refit'](fact, brot, epochs, lr)
    print('refit done', fact, round(time.time() - t0), 's', flush=True)
from dfthcsos.verify import load_inputs, verify as verify_dfthcsos
sys.path.insert(0, str(VER / 'dfthcsos' / 'example_submission'))
from angle_k_reference_circuit import build_angle_k_family
inputs = load_inputs(str(fact))
edges = ns['dft_log_depth_edges'](inputs.num_orbitals)
words = ns['dft_encode'](inputs, edges, brot)
J = ns['dft_batch_size'](inputs.num_vectors, brot)
circuits = build_angle_k_family(words, brot, J)
sub = outdir / 'submission'
sub.mkdir(exist_ok=True)
for f in sub.glob('*'):
    f.unlink()
for k, c in enumerate(circuits):
    c.to_file(str(sub / f'angle_k_{k}.kmb'), format='kmb')
with np.load(fact) as z:
    np.savez(sub / 'submission.npz', **{k: z[k] for k in ('U', 'W', 'H2_I', 'H1_I', 'Hsym_lower')},
             tree_edges=edges, angle_words=words)
try:
    r = verify_dfthcsos(str(sub))
    m = r['metrics']
    res = dict(status=r['status'], score=r['score'], ratio=3355599588.3155375 / r['score'], lambda_eff=r['lambda_eff'],
               toffolis=r['toffoli_count'], qubits=r['num_qubits'], brot=brot, J=J, K=len(circuits),
               error_1_mHa=1e3 * m['error_1'], error_2_mHa=1e3 * m['error_2'])
except Exception as ex:
    from dfthcsos.verify import evaluate_words
    met = evaluate_words(inputs, words, edges, brot)
    res = dict(status='FAIL', error=str(ex)[:300], metrics={k: float(v) for k, v in met.items() if np.ndim(v) == 0})
print(json.dumps(res, indent=1), flush=True)
