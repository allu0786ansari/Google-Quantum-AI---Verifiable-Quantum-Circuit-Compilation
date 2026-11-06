# %% [markdown]
# ## 4. dfthcsos: chemistry block encoding
#
# **Contract.** Submit a factorization of the FeMoCo Hamiltonian (`submission.npz`)
# plus the circuits `angle_k_0 ... angle_k_{K-1}` that load each vector's Givens
# angles as `brot`-bit words, `J` words per circuit. The verifier simulates the
# circuits, rebuilds the Hamiltonian from the words they load, and requires both
# residuals to be at most 20 mHa.
# Cost = lambda_eff x Toffolis x qubits for the whole block encoding.
#
# **Baseline:** brot = 12, J = 8, K = 10 circuits, each a plain lookup over 3,976
# addresses. The angle circuits are 85% of the Toffolis, so fewer bits per angle
# means more angles per circuit and fewer circuits:
#
# | brot | J (128-qubit cap) | K | angle Toffolis |
# |---|---|---|---|
# | 12 | 8 | 10 | 79,480 |
# | 11 | 9 | 9 | 71,532 |
# | 10 | 10 | 8 | 63,584 |
#
# Lower precision raises the two-body error (11 bits: 25.0 mHa with plain rounding,
# over the limit). Our fix is to round the two-body vectors first, freeze them, then
# refit the continuous weights to those rounded vectors with the starter kit's
# trainer, so the weights absorb the rounding error (optional cell below; needs JAX).
#
# The pipeline below works for any factorization and precision. It rebuilds the
# one-body rows as eigenvectors after the two-body rounding (as the starter kit
# does), picks the largest J that fits the qubit cap, builds the circuits with
# the competition's reference lookup builder, verifies them, and falls back to
# the baseline if anything fails.

# %%
DFT_BROT = 12                  # bits per angle word
DFT_FACTORIZATION = None       # None = the baseline in sample_submission/; or a path to an .npz
DFT_REFIT_EPOCHS = 0           # > 0: refit the weights to the rounded vectors first (needs JAX + starter kit)
DFT_REFIT_LR = 1e-5
STARTER_KIT = None             # local starter kit folder (the one holding starter_kit/); found automatically on Kaggle


def dft_log_depth_edges(n):
    """Breadth-first edges of the log-depth binary Givens tree (arXiv:2605.30455)."""
    edges, level, depth = [], [0], 0
    while len(edges) < n - 1:
        nxt = []
        for label in level:
            nxt.append(label)
            right = label + (1 << depth)
            if right < n:
                edges.append((label, right))
                nxt.append(right)
        level, depth = nxt, depth + 1
    return np.asarray(edges, dtype=np.int64).reshape(-1, 2)


def dft_tree_angles(vectors, edges):
    """theta / pi in [0, 2) such that the tree maps e_0 to each (normalized) row."""
    u = np.array(vectors, dtype=np.float64, ndmin=2)
    u /= np.linalg.norm(u, axis=1, keepdims=True)
    theta = np.zeros((u.shape[0], u.shape[1] - 1))
    for t in range(u.shape[1] - 2, -1, -1):
        p, c = edges[t]
        theta[:, t] = np.arctan2(u[:, c], u[:, p])
        u[:, p] = np.hypot(u[:, p], u[:, c])
        u[:, c] = 0.0
    return np.mod(theta / np.pi, 2.0)


def dft_words(theta_over_pi, brot):
    return np.floor(np.mod(theta_over_pi, 2.0) * (1 << (brot - 1)) + 0.5).astype(np.int64) % (1 << brot)


def dft_encode(inputs, edges, brot):
    """(N + R*B, N - 1) angle words: two-body rows rounded, then one-body rows as the
    eigenvectors of h1_prime built from the decoded two-body vectors."""
    from dfthcsos.scoreinfra.dfthc_math import reconstruct_u_vectors, reconstruct_x, build_h1_prime
    U = np.asarray(inputs.params['U'], dtype=np.float64)
    R, B, N = U.shape
    two = dft_words(dft_tree_angles(U.reshape(-1, N), edges), brot)
    vec = reconstruct_u_vectors(two / (1 << (brot - 1)), edges)
    vec /= np.linalg.norm(vec, axis=1, keepdims=True)
    X = reconstruct_x(vec.reshape(R, B, N), inputs.w)
    X = np.tril(X) + np.swapaxes(np.tril(X, k=-1), -1, -2)
    h1 = inputs.ref['oei'] - 0.5 * np.einsum('prrq->pq', inputs.ref['tei'])
    hp = build_h1_prime(h1, inputs.S, inputs.beta_1, inputs.ref['electrons'], inputs.W, X)
    one = dft_words(dft_tree_angles(np.linalg.eigh(hp)[1].T, edges), brot)
    return np.vstack([one, two])


def dft_batch_size(num_vectors, brot, cap=128):
    """Largest J with selection + J*brot + (selection - 1) unary-iteration ancillas <= cap."""
    sel = (num_vectors - 1).bit_length()
    return (cap - sel - (sel - 1)) // brot


def find_starter_kit():
    tried = []
    if STARTER_KIT:
        tried.append(str(STARTER_KIT))
        hits = sorted(Path(STARTER_KIT).expanduser().glob('**/dfthc/df_thc_sga.py'))
        if hits:
            return hits[0].parents[1]
    try:
        import kagglehub
        p = kagglehub.dataset_download('sohier/verifiable-circuit-compilation-starter-kit')
        hits = sorted(Path(p).glob('**/dfthc/df_thc_sga.py'))
        if hits:
            return hits[0].parents[1]
        tried.append(f'kagglehub ({p})')
    except Exception as ex:
        tried.append(f'kagglehub failed: {ex}')
    hits = sorted(Path('/kaggle/input').glob('**/dfthc/df_thc_sga.py')) if Path('/kaggle/input').exists() else []
    if hits:
        return hits[0].parents[1]
    raise FileNotFoundError('starter kit not found: ' + '; '.join(tried))


def dft_refit(seed_npz, brot, epochs, lr):
    """Round the two-body vectors to brot bits, freeze them, refit the weights.
    Returns the path of the refitted factorization."""
    kit = find_starter_kit()
    sys.path.insert(0, str(kit))
    import jax.numpy as jnp
    from dfthc import df_thc_sga
    from dfthc.train_df_thc_sga import Initialize, Train
    from starter.task1_factorization.train_df_thc_qaggle import load_seed_params, frozen_symmetry_params, save_solution
    fcidump = VERIFIERS / 'dfthcsos' / 'FeMoCo76o113e.fcidump'
    ref = Initialize(str(fcidump), -22140.27335)
    params, _ = load_seed_params(str(seed_npz), orbitals=ref['orbitals'])
    U = np.asarray(params['U'], dtype=np.float64)
    R, B, N = U.shape
    edges = dft_log_depth_edges(N)
    from dfthcsos.scoreinfra.dfthc_math import reconstruct_u_vectors
    Uq = reconstruct_u_vectors(dft_words(dft_tree_angles(U.reshape(-1, N), edges), brot) / (1 << (brot - 1)), edges)
    Uq = (Uq / np.linalg.norm(Uq, axis=1, keepdims=True)).reshape(R, B, N)
    shifts = np.load(VERIFIERS / 'dfthcsos' / 'FeMoCo76o113e.fixed_shifts.npz')
    frozen = {'U': jnp.asarray(Uq, dtype=df_thc_sga.dtype),
              'Hsym_lower': jnp.asarray(shifts['Hsym_lower'], dtype=df_thc_sga.dtype),
              'H1_I': jnp.asarray(shifts['H1_I'], dtype=df_thc_sga.dtype)}
    sol = Train(generating_set='Qaggle', method='DF-THC-SOS', epochs=epochs, learning_rate=lr, rho=1e-4,
                params={**params, **frozen}, ref=ref, print_interval=max(1, epochs // 10), epsilon_dfthc=0.0,
                frozen_params=frozen, lambda_floor=49.0, sos_gap_floor=3.914898, learning_rate_final=lr / 5)
    out = OUT / f'refit_b{brot}_e{epochs}.npz'
    if out.exists():
        out.unlink()
    save_solution(str(out), sol['params_best'], sol['values'], ref=ref, metadata={'brot': brot, 'epochs': epochs})
    return out


def run_dfthcsos():
    from dfthcsos.verify import load_inputs, verify as verify_dfthcsos
    sys.path.insert(0, str(VERIFIERS / 'dfthcsos' / 'example_submission'))
    from angle_k_reference_circuit import build_angle_k_family
    t0 = time.time()
    src = Path(DFT_FACTORIZATION) if DFT_FACTORIZATION else SAMPLE / 'dfthcsos' / 'submission.npz'
    if DFT_REFIT_EPOCHS > 0:
        src = dft_refit(src, DFT_BROT, DFT_REFIT_EPOCHS, DFT_REFIT_LR)
    inputs = load_inputs(str(src))
    edges = dft_log_depth_edges(inputs.num_orbitals)
    words = dft_encode(inputs, edges, DFT_BROT)
    J = dft_batch_size(inputs.num_vectors, DFT_BROT)
    circuits = build_angle_k_family(words, DFT_BROT, J)
    out = SUB / 'dfthcsos'
    for f in out.glob('*'):
        f.unlink()
    for k, c in enumerate(circuits):
        c.to_file(str(out / f'angle_k_{k}.kmb'), format='kmb')
    with np.load(src) as z:
        np.savez(out / 'submission.npz', **{k: z[k] for k in ('U', 'W', 'H2_I', 'H1_I', 'Hsym_lower')},
                 tree_edges=edges, angle_words=words)
    r = verify_dfthcsos(str(out))          # raises if any gate fails
    m = r['metrics']
    return dict(brot=DFT_BROT, batch_size=J, circuits=len(circuits), score=r['score'], lambda_eff=r['lambda_eff'],
                toffolis=r['toffoli_count'], qubits=r['num_qubits'], error_1_mHa=round(1e3 * m['error_1'], 3),
                error_2_mHa=round(1e3 * m['error_2'], 3), seconds=round(time.time() - t0, 1))


def use_baseline_dfthcsos():
    out = SUB / 'dfthcsos'
    for f in out.glob('*'):
        f.unlink()
    for f in sorted((SAMPLE / 'dfthcsos').iterdir()):
        shutil.copy(f, out / f.name)


if USE_IMPROVED['dfthcsos']:
    try:
        REPORT['dfthcsos'] = dict(version='improved', **run_dfthcsos())
    except Exception as ex:
        print('dfthcsos: improved submission rejected, using baseline:', ex)
        use_baseline_dfthcsos()
        REPORT['dfthcsos'] = dict(version='baseline (fallback)', error=str(ex))
else:
    use_baseline_dfthcsos()
    REPORT['dfthcsos'] = dict(version='baseline')
    if VERIFY in ('quick', 'full'):
        from dfthcsos.verify import verify as verify_dfthcsos
        r = verify_dfthcsos(str(SUB / 'dfthcsos'))
        REPORT['dfthcsos'].update(verifier_status=r['status'], score=r['score'])
print(json.dumps(REPORT['dfthcsos'], indent=2, default=str))
