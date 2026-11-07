import sys, os, time
import numpy as np
sys.path.insert(0, os.environ['VERIFIERS'])
from dfthcsos.verify import load_inputs, evaluate_words
from dfthcsos.scoreinfra import dfthc_math as M
src = sys.argv[1]; brots = [int(b) for b in sys.argv[2:]]
inp = load_inputs(src)
z = np.load(src); edges = z['tree_edges']; w8 = z['angle_words']; b0 = int(z['brot'])
N = inp.num_orbitals
print('reference', {k: round(float(v), 5) for k, v in evaluate_words(inp, w8, edges, b0).items() if k in ('error_1','error_2','lambda_eff','E_gap')})
theta = w8 / (1 << (b0 - 1))                     # theta/pi in [0, 2)
def words_at(t, b):
    return np.floor(np.mod(t, 2.0) * (1 << (b - 1)) + 0.5).astype(np.int64) % (1 << b)
for b in brots:
    t0 = time.time()
    two = words_at(theta[N:], b)
    vec = M.reconstruct_u_vectors(two / (1 << (b - 1)), edges); vec /= np.linalg.norm(vec, axis=1, keepdims=True)
    w = np.asarray(inp.w); R, B, C = w.shape
    X = M.reconstruct_x(vec.reshape(R, B, N), w); X = np.tril(X) + np.swapaxes(np.tril(X, k=-1), -1, -2)
    h1 = inp.ref['oei'] - 0.5 * np.einsum('prrq->pq', inp.ref['tei'])
    hp = M.build_h1_prime(h1, inp.S, inp.beta_1, inp.ref['electrons'], inp.W, X)
    v1 = np.linalg.eigh(hp)[1].T
    one = words_at(M.__dict__.get('_tree_angles', None) and 0 or 0, b) if False else None
    # tree angles for one-body eigenvectors
    u = v1.copy(); th = np.zeros((N, N - 1))
    for k in range(N - 2, -1, -1):
        p, c = edges[k]; th[:, k] = np.arctan2(u[:, c], u[:, p]); u[:, p] = np.hypot(u[:, p], u[:, c]); u[:, c] = 0
    one = words_at(np.mod(th / np.pi, 2.0), b)
    m = evaluate_words(inp, np.vstack([one, two]), edges, b)
    print(f'brot {b}:', {k: round(float(v), 5) for k, v in m.items() if k in ('error_1','error_2','lambda_eff','E_gap')}, f'{time.time()-t0:.0f}s')
