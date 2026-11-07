# %% [markdown]
# ## 1b. Public frontier components (credited)
#
# The strongest public circuits for all three tasks are built by many Kaggle
# participants and shared under Apache-2.0. We start from them, verify every one
# with the official verifiers, and keep our own circuits (sections 2-4) as the
# fallback. They come from **ahoirg's [QCC Reproducible Circuit Assets](https://www.kaggle.com/datasets/ahoirg/qcc-reproducible-circuit-assets), release 2**
# (downloaded with `kagglehub.dataset_download` in the setup cell):
#
# | Task | Component | Credits (see the dataset's CREDITS files for the full lineage) |
# |---|---|---|
# | invsqrt | 62-qubit circuit with five proved rewrites | mrmeows / team selenium (architecture), busyaprime, yoyobar, Ahmed Berat Özer, Ehan Mallick, SomeJoint; rewrite matcher by Kirill Tymchishin; lineage heliosli, Phuc Thai, Dmitrii Gluzdov, marsyas |
# | dqi | 891-qubit decoder | SomeJoint / Agniv Banerjee and the team credited in their source (Ahmed Berat Özer, Alex Chilton, jday96314, mrmeows, manwithacat, Dmitrii Gluzdov, ViberTeam, Phuc Thai), yoyobar, ahoirg |
# | dfthcsos | C39 factorization + paired-product loaders | marculera (factorization pipeline, tensorqrom), Syed Asad Ali / SomeJoint (C39 loaders), mrmeows, heliosli, Ahmed Berat Özer, kozykappa, Michael Rustad and others |
#
# **Newer public release (preferred when attached):** heliosli's
# [qcc-circuits-live](https://www.kaggle.com/code/heliosli/the-dataset-is-public-now-enjoy) dataset,
# portfolio `v812e4edf89e_y33_ay3` (local GM 24.94): yoyobar's version-33 62-qubit invsqrt; their
# update of yoyobar's version-38 892-qubit decoder (lineage kirilltym, daifanhao, syed3000, SomeJoint,
# ahoirg, Atakan Aldemir and the earlier decoder line); syed3000's (8, 143, 44) chemistry table and
# angle circuits with weights refined from maoxin123's refit (lineage marculera, Ahmed Berat Özer,
# heliosli and others). Each file is used only if its SHA-256 matches the pinned value below.
#
# What this notebook runs from the ahoirg dataset (all hash-checked against the
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
    roots = [Path(p) for p in (LOCAL_ASSETS, os.environ.get('VQCC_ASSETS'), ASSETS_PATH) if p]
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


# Newer public release (heliosli/qcc-circuits-live, portfolio v812e4edf89e_y33_ay3, local GM 24.94):
# prebuilt files only, each pinned here by SHA-256. Used per file when present and exact; otherwise the
# ahoirg build above is kept.
PUBLIC_RELEASE_DATASET = 'heliosli/qcc-circuits-live'
PUBLIC_RELEASE_PORTFOLIO = 'v812e4edf89e_y33_ay3'
PUBLIC_RELEASE_SHA256 = {
    'invsqrt.kmb': 'c17971aff984830dc68f4f1459cdc549fd2637dc656bc7d3d70aeabc7c8eb775',
    'dqi.kmb': '4de6b023fbf146d29cc608fd08c15ce12c298f5771fd0d8a678efc403929d212',
    'dfthcsos/submission.npz': '388b6a7ee7d882e9e8d165a79812c3ecb6c7e2486de59fd63733ad82cf61baac',
    'dfthcsos/angle_k_0.kmb': 'd3aab7479fe3b8510e3af3797a14ec8d9bfa226d7aa7717ae4e93f496d2d9538',
    'dfthcsos/angle_k_1.kmb': '4595e5dcc982295db2520c8497219e1768757edf3aa9e706fa44b1716f996ed2',
    'dfthcsos/angle_k_2.kmb': '77dec28c5ebffe723fb58717a31c017d66c44cbec220bc3cc5a15a2fff051aa3',
    'dfthcsos/angle_k_3.kmb': '8e137a014a1325450eaf9ec15bc824475b757b3a05bbc542eed9cccfd6b5ec53',
    'dfthcsos/angle_k_4.kmb': 'd9f2a0d739142a1d20789352c617cccc702a8f0b4b17f91971e4e4ef3e6f6c3f',
    'dfthcsos/angle_k_5.kmb': 'bd0886bb6c3d5c04f5e1fd29dbd2562d1c7e5b1a7406bc1abde6819093a0880e',
    'dfthcsos/angle_k_6.kmb': '12c0c6f92d708ff7d852f2c6f36a40ecf4369249fa1f0ab4527b08e0152e7558',
    'dfthcsos/angle_k_7.kmb': '93dc9742392a0dd6e754aebbbea1519014292a56b6908dfce5ae417f2cde73d6',
    'dfthcsos/angle_k_8.kmb': 'd43bd96b88e6ef38c7ea9e697856646af4259889b75d1bc67a94a1768670b5a6',
    'dfthcsos/angle_k_9.kmb': 'a5960a572c39a3996b3919eb8bc18066b8aa78d3fc83bc7435d7e5550de3a2a6',
    'dfthcsos/angle_k_10.kmb': 'f28f2ead2b7eb3d45aadf4e0bb9f9399c7d8f016247d7020ae39dbd9c062cace',
    'dfthcsos/angle_k_11.kmb': '7b9094ddc215c6c7db223fe319e6d2a6f67c1669f412e4d3281c280723df5b81',
    'dfthcsos/angle_k_12.kmb': 'd49dfbfdb8cff3c3d680e4ef65a3d1a01db7a2fc666f49a522d7f58f50305e78',
    'dfthcsos/angle_k_13.kmb': 'cb7e48a7029323a6ec11c0f6dede5f5e0443f67da92e0e6554ff9120b8279efd',
    'dfthcsos/angle_k_14.kmb': '07fcdcb75d3a9e7f2573e057376ec20b7a56bf9a53fe7d7f00fb351614c2003f',
}
LOCAL_RELEASE = None   # local folder holding QCC_RELEASE_MANIFEST.json; found automatically on Kaggle


def load_public_release():
    """Copy the pinned public-release files that are present and exact. Returns {task: path}."""
    import lzma
    roots = [Path(p) for p in (LOCAL_RELEASE, os.environ.get('QCC_DS_DIR')) if p]
    got = download_dataset(PUBLIC_RELEASE_DATASET)
    if got:
        roots.append(got)
    roots.append(Path('/kaggle/input'))
    manifest_path = next((m for r in roots if r.exists() for m in _find_files(r, 'QCC_RELEASE_MANIFEST.json')), None)
    if manifest_path is None:
        print(f'{PUBLIC_RELEASE_DATASET}: not attached; keeping the ahoirg build')
        return {}
    ds = manifest_path.parent
    manifest = json.loads(manifest_path.read_text())
    if manifest.get('portfolio') != PUBLIC_RELEASE_PORTFOLIO:
        print(f"{PUBLIC_RELEASE_DATASET}: portfolio {manifest.get('portfolio')} is not the pinned "
              f'{PUBLIC_RELEASE_PORTFOLIO}; only files with pinned hashes are used')
    out = OUT / 'public_release'
    if out.exists():
        shutil.rmtree(out)
    exact = set()
    for item in manifest['files']:
        name = item['file']
        if name not in PUBLIC_RELEASE_SHA256:
            continue
        src, dst = ds / item['stored_as'], out / name
        if not src.exists():
            continue
        dst.parent.mkdir(parents=True, exist_ok=True)
        if item['stored_as'].endswith('.xz'):
            with lzma.open(src) as fi, open(dst, 'wb') as fo:
                shutil.copyfileobj(fi, fo, 1 << 22)
        else:
            shutil.copyfile(src, dst)
        if _sha256(dst) == PUBLIC_RELEASE_SHA256[name]:
            exact.add(name)
        else:
            dst.unlink()
    found = {}
    if 'invsqrt.kmb' in exact:
        found['invsqrt'] = out / 'invsqrt.kmb'
    if 'dqi.kmb' in exact:
        found['dqi'] = out / 'dqi.kmb'
    chem = [n for n in PUBLIC_RELEASE_SHA256 if n.startswith('dfthcsos/')]
    if all(n in exact for n in chem):
        found['dfthcsos'] = out / 'dfthcsos'
    print(f'{PUBLIC_RELEASE_DATASET}: exact files for {sorted(found)}')
    return found


FRONTIER = {}
if any(USE_FRONTIER.values()):
    problems = []
    try:
        FRONTIER = load_frontier()
    except Exception as ex:
        problems.append(f'ahoirg build: {type(ex).__name__}: {ex}')
    try:
        FRONTIER.update(load_public_release())
    except Exception as ex:
        problems.append(f'{PUBLIC_RELEASE_DATASET}: {type(ex).__name__}: {ex}')
    missing = [t for t, on in USE_FRONTIER.items() if on and t not in FRONTIER]
    if missing:
        msg = (f'frontier circuits missing for {missing}, so the submission would score far lower. '
               f'Attach {ASSETS_DATASET} and/or {PUBLIC_RELEASE_DATASET}, or set REQUIRE_FRONTIER = False. '
               f'Causes: {problems}')
        if REQUIRE_FRONTIER:
            raise RuntimeError(msg)
        print(msg)
print('frontier sources:', {k: str(v) for k, v in FRONTIER.items()})
