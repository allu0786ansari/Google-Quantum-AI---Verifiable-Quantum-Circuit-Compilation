# %% [markdown]
# ## 1b. Public frontier components (credited)
#
# The strongest public circuits for all three tasks are built by many Kaggle
# participants and shared under Apache-2.0. We start from them, verify every one
# with the official verifiers, and keep our own circuits (sections 2-4) as the
# fallback. They come from **ahoirg's [QCC Reproducible Circuit Assets](https://www.kaggle.com/datasets/ahoirg/qcc-reproducible-circuit-assets), release 2**
# (attach it to this notebook):
#
# | Task | Component | Credits (see the dataset's CREDITS files for the full lineage) |
# |---|---|---|
# | invsqrt | 62-qubit circuit with five proved rewrites | mrmeows / team selenium (architecture), busyaprime, yoyobar, Ahmed Berat Özer, Ehan Mallick, SomeJoint; rewrite matcher by Kirill Tymchishin; lineage heliosli, Phuc Thai, Dmitrii Gluzdov, marsyas |
# | dqi | 891-qubit decoder | SomeJoint / Agniv Banerjee and the team credited in their source (Ahmed Berat Özer, Alex Chilton, jday96314, mrmeows, manwithacat, Dmitrii Gluzdov, ViberTeam, Phuc Thai), yoyobar, ahoirg |
# | dfthcsos | C39 factorization + paired-product loaders | marculera (factorization pipeline, tensorqrom), Syed Asad Ali / SomeJoint (C39 loaders), mrmeows, heliosli, Ahmed Berat Özer, kozykappa, Michael Rustad and others |
#
# What this notebook runs from that dataset (all hash-checked against the
# dataset's manifest first):
#
# - `inverse_factor.py`: an exact five-window gate rewrite of the invsqrt circuit,
#   each edit proved on all 128 inputs of its 7 wires.
# - `scheduled_emit.py` and `tensorqrom.py`: emit the 15 chemistry loaders from the
#   C39 angle table and its algebraic plans.
# - `v2-dqi.kmb.xz`: the prebuilt dqi decoder, decompressed.

# %%
USE_FRONTIER = {'invsqrt': True, 'dqi': True, 'dfthcsos': True}
LOCAL_ASSETS = None          # local folder holding qcc-reproducible-circuit-assets(-v2); found automatically on Kaggle
ASSETS_DATASET = 'ahoirg/qcc-reproducible-circuit-assets'


def _sha256(path):
    import hashlib
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(1 << 24), b''):
            h.update(chunk)
    return h.hexdigest()


def find_assets():
    """Folder holding release 2 of the assets (v2-assets-manifest.json), or None."""
    roots = [Path(p) for p in (LOCAL_ASSETS, os.environ.get('VQCC_ASSETS')) if p]
    try:
        import kagglehub
        roots.append(Path(kagglehub.dataset_download(ASSETS_DATASET)))
    except Exception as ex:
        print('kagglehub assets lookup failed:', ex)
    roots.append(Path('/kaggle/input'))
    for root in roots:
        if root.exists():
            for hit in sorted(root.expanduser().glob('**/v2-assets-manifest.json')):
                return hit.parent
    return None


def load_frontier():
    """Authenticate the assets, extract the source bundle and build the frontier
    invsqrt and chemistry circuits. Returns a dict of paths, or raises."""
    assets = find_assets()
    if assets is None:
        raise FileNotFoundError(f'attach the {ASSETS_DATASET} dataset (release 2)')
    manifest = json.loads((assets / 'v2-assets-manifest.json').read_text())
    for name, pin in manifest['files'].items():
        if (assets / name).exists() and _sha256(assets / name) != pin['sha256']:
            raise ValueError(f'{name}: SHA-256 differs from the assets manifest')
    work = OUT / 'frontier_build'
    if work.exists():
        shutil.rmtree(work)
    src = work / 'source'
    src.mkdir(parents=True)
    import hashlib
    bundle = json.loads((assets / 'v2-source-bundle.json').read_text())
    for name, entry in bundle['files'].items():
        if hashlib.sha256(entry['text'].encode()).hexdigest() != entry['sha256']:
            raise ValueError(f'source bundle member {name} fails its hash')
        (src / name).parent.mkdir(parents=True, exist_ok=True)
        (src / name).write_text(entry['text'])
    sys.path.insert(0, str(src))
    from inverse_factor import factor_windows, prove_edits
    from scheduled_emit import build_circuit
    pins = json.loads((src / 'manifest.json').read_text())
    # invsqrt: base circuit -> five proved rewrites
    base = assets / 'v2-inverse-control.kmb'
    text = work / 'inverse_control.kmx'
    km.Circuit.from_file(str(base), format='kmb').to_file(str(text), format='kmx')
    lines, edits = factor_windows(text.read_text().splitlines())
    proof = prove_edits(edits)
    inv_path = work / 'invsqrt.kmb'
    tmp = work / 'inverse_transformed.kmx'
    tmp.write_text('\n'.join(lines) + '\n')
    km.Circuit.from_file(str(tmp), format='kmx').to_file(str(inv_path), format='kmb')
    # chemistry: 15 loaders from the C39 table and its plans
    chem = work / 'dfthcsos'
    chem.mkdir()
    factor = assets / 'v2-c39-factor.npz'
    with np.load(factor, allow_pickle=False) as z:
        words = z['angle_words']
    for slot in range(15):
        plan = json.loads((src / 'plans' / f'slot-{slot}.json').read_text())
        cur = words[:, slot * 5:slot * 5 + 5]
        prev = words[:, (slot - 1) * 5:slot * 5] if slot else np.zeros_like(cur)
        delta = np.vstack([cur ^ prev, np.asarray(plan['unused_address_completion'], dtype=np.int64)])
        terms = [tuple(int(v, 16) for v in t) for t in plan['terms']]
        build_circuit(delta, terms, plan['low_cache_schedule']).to_file(str(chem / f'angle_k_{slot}.kmb'), format='kmb')
    shutil.copyfile(factor, chem / 'submission.npz')
    pinned = {**{k: v['sha256'] for k, v in pins['output_members'].items()}}
    mism = [k for k, h in pinned.items() if (work / k).exists() and _sha256(work / k) != h]
    if mism:
        raise ValueError(f'rebuilt files differ from the pinned outputs: {mism}')
    # dqi: decompress the prebuilt decoder
    import lzma
    dqi_path = work / 'dqi.kmb'
    with lzma.open(assets / 'v2-dqi.kmb.xz') as fin, open(dqi_path, 'wb') as fout:
        shutil.copyfileobj(fin, fout, 1 << 24)
    if _sha256(dqi_path) != manifest['build_inputs']['dqi']['sha256']:
        raise ValueError('decompressed dqi.kmb fails its hash')
    print('frontier built:', dict(proof=proof, files_checked=len(pinned) + 1))
    return dict(invsqrt=inv_path, dqi=dqi_path, dfthcsos=chem)


FRONTIER = {}
if any(USE_FRONTIER.values()):
    try:
        FRONTIER = load_frontier()
    except Exception as ex:
        print('frontier components unavailable, using our own circuits:', ex)
