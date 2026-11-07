import numpy as np, jax, jax.numpy as jnp
import fastfit as F
P = F.Problem('/home/claude/work/lt/output/submission/dfthcsos/submission.npz')
parts = P.jax_parts(); w0, W0 = P.w0, P.W0; nw = w0.size
def unpack(x): return x[:nw].reshape(w0.shape), x[nw:].reshape(W0.shape)
x0 = jnp.asarray(np.concatenate([w0.ravel(), W0.ravel()]))
gl = np.asarray(jax.grad(lambda x: parts(*unpack(x))['lam'])(x0))
ge = np.asarray(jax.grad(lambda x: parts(*unpack(x))['e2'] / 1e-3)(x0))
gg = np.asarray(jax.grad(lambda x: parts(*unpack(x))['gap'])(x0))
A = np.stack([ge, -gg], 1)
nu, *_ = np.linalg.lstsq(A, -gl, rcond=None)
r = gl + A @ nu
print('nu (per mHa, gap):', nu, '|grad lam|', np.linalg.norm(gl), 'KKT residual', np.linalg.norm(r), 'rel', np.linalg.norm(r)/np.linalg.norm(gl))
# predicted first-order gain: move along -r (projected), keeping e2 and gap fixed to first order
d = -r / np.linalg.norm(r)
for s in (1e-6, 1e-5, 3e-5, 1e-4, 3e-4, 1e-3):
    p = {k: float(v) for k, v in parts(*unpack(x0 + s * d)).items()}
    print(f's {s:g}: lam {p["lam"]:.6f} e2 {p["e2"]*1e3:.5f} gap {p["gap"]:.7f}')
