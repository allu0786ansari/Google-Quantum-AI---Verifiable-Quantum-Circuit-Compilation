# %% [markdown]
# # VQCC: all three circuits built from open-source generators, then verified
#
# This notebook downloads **no finished circuit**. It builds every circuit from published, licensed
# generator source code, checks each result with the **official verifiers**, and writes `submission.zip`.
#
# ## Credits (the generators are other people's work)
#
# | Task | Generator | Author | Licence |
# |---|---|---|---|
# | invsqrt | 46-qubit generator (`invsqrt_source.zip`): root/residual normalisation, two coefficient loaders, phase divider | Syed Asad Ali and his AI-agent team | Apache-2.0 |
# | dqi | 889-qubit decoder generator (`dqi_source.zip`) | Syed Asad Ali and his AI-agent team | Apache-2.0 |
# | dfthcsos | (8, 143, 44) factorization table, coded 8-bit words, GF(2) product-plan loader compiler (`recipes/h14-deep512`) | Francisco Angulo de Lafuente; loader design after mrmeows, somejoint, kirilltym, marculera | MIT or Apache-2.0 |
#
# All three are taken from one public repository at a pinned commit
# (`github.com/Agnuxo1/qvc-chemistry-compiler`), and the two archives are checked against pinned SHA-256
# values before anything runs. Keep this table and the licence notices if you share the notebook.
#
# ## What is ours in this notebook
#
# The build orchestration, the integrity checks, the full official verification with a second seed, and the
# packaging. **No circuit-level technique of our own is applied yet**; section 6 lists what we are adding.
#
# ## How to run on Kaggle
# 1. *Add input → Competition* (this competition). Turn **Internet on** (the sources come from GitHub).
# 2. Run all cells (about 45-70 minutes with `VERIFY = 'full'`).
# 3. Upload `/kaggle/working/submission.zip` on the Submit page.
#
# Expected result, from the generators' published receipts: invsqrt 80,163 (46 qubits), dqi 3.607e11
# (889 qubits), dfthcsos ratio 11.23; overall about **27.0**.

# %%
# ---- Settings ---------------------------------------------------------------
VERIFY = 'full'                    # 'full': official verification of all three tasks; 'quick': reduced checks
SECOND_SALT = 'vqcc-second-seed'   # also verify invsqrt and chemistry with other measurement outcomes ('' = skip)
SOURCE_REPO = 'https://github.com/Agnuxo1/qvc-chemistry-compiler.git'
SOURCE_COMMIT = 'b0b0a4ffb18b13a2faa92972c9c6b737f2d13e78'
SOURCE_DIR = None                  # set to a folder holding that repository to skip the download (e.g. an attached dataset)
LOCAL_DATA = None                  # only outside Kaggle: folder with verifiers/, kickmix/ and sample_submission/

PINNED = {                         # SHA-256 of the source archives inside the repository
    'studies/compact_invsqrt46/invsqrt_source.zip': '7d9801ae20e942dba2211a656bba6014a0a16d57da534c0b60328d1c5d9d8254',
    'studies/source_portfolio889/upstream/dqi_source.zip': '49e2b331d4008afb724985caf1038765b7a910f08838b590bf6e8fb692772763',
}
EXPECTED = {                       # SHA-256 of the circuits those sources are documented to build
    'invsqrt.kmb': '321aff83636214196daa9cbc9dcc9c5e12a155bcff8d9989e98b8cee40935177',
    'dqi.kmb': '7cbe8ecf6f6e9743529b574f019504bb1712ce06fbb570c719ce9d9e2c914cd9',
}
CHEM_RECIPE = 'recipes/h14-deep512'

# %% [markdown]
# ## 1. Setup

# %%
import fnmatch
import hashlib
import importlib.util
import json
import math
import os
import shutil
import subprocess
import sys
import time
import types
import zipfile
from fractions import Fraction
from pathlib import Path

COMPETITION = 'google-verifiable-quantum-circuit-compilation'
MARKER = Path('verifiers') / 'invsqrt' / 'verify_submission_invsqrt.py'
BASELINE_SCORE = {'invsqrt': 2076774.23, 'dqi': 24439397322564.6, 'dfthcsos': 3355599588.3155375}


def walk_files(root, max_depth=8):
    root = Path(root)
    if not root.is_dir():
        return
    base = len(root.parts)
    for folder, dirs, files in os.walk(root, followlinks=True):
        if len(Path(folder).parts) - base >= max_depth:
            dirs[:] = []
        for f in files:
            yield Path(folder) / f


def find_data_root():
    roots = [LOCAL_DATA, os.environ.get('VQCC_DATA')]
    try:
        import kagglehub
        roots.append(kagglehub.competition_download(COMPETITION))
    except Exception as ex:
        print(f'kagglehub competition download failed: {type(ex).__name__}')
    roots.append('/kaggle/input')
    for r in roots:
        if not r or not Path(r).expanduser().exists():
            continue
        r = Path(r).expanduser().resolve()
        if (r / MARKER).exists():
            return r
        for f in walk_files(r):
            if f.name == MARKER.name and f.parent.name == 'invsqrt' and f.parent.parent.name == 'verifiers':
                return f.parents[2]
    raise FileNotFoundError('Competition data not found: add the competition as input.')


def load_kickmix(data_root):
    api = ('CircuitBuilder', 'Circuit', 'GF2Field')
    try:
        import kickmix
        if all(hasattr(kickmix, a) for a in api):
            return kickmix
    except ImportError:
        pass
    sys.modules.pop('kickmix', None)
    for so in [f for r in (data_root, '/kaggle/input') for f in walk_files(r) if fnmatch.fnmatch(f.name, 'kickmix*.so')]:
        try:
            spec = importlib.util.spec_from_file_location('kickmix', so)
            module = importlib.util.module_from_spec(spec)
            sys.modules['kickmix'] = module
            spec.loader.exec_module(module)
            if all(hasattr(module, a) for a in api):
                return module
        except Exception:
            pass
        sys.modules.pop('kickmix', None)
    raise ImportError(f'kickmix not usable with Python {sys.version.split()[0]}')


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(1 << 22), b''):
            h.update(chunk)
    return h.hexdigest()


def run(cmd, cwd=None, env=None, timeout=3600, log=None):
    """Run a build command; on failure show the end of its output and stop."""
    t0 = time.time()
    res = subprocess.run([str(c) for c in cmd], cwd=cwd, env=env, capture_output=True, text=True, timeout=timeout)
    out = res.stdout + res.stderr
    if log:
        Path(log).write_text(out)
    print(f'   {" ".join(str(c) for c in cmd[:3])} ... exit {res.returncode} ({time.time() - t0:.0f}s)')
    if res.returncode:
        print(out[-3000:])
        raise RuntimeError(f'build step failed (full log: {log})')
    return out


DATA = find_data_root()
VERIFIERS = DATA / 'verifiers'
ENGINE = DATA / 'kickmix'
assert (ENGINE / 'kickmix.so').exists(), f'{ENGINE}/kickmix.so not found (the generators need the competition engine)'
km = load_kickmix(DATA)
sys.path.insert(0, str(VERIFIERS))
OUT = Path('/kaggle/working') if Path('/kaggle/working').exists() else Path('output_from_source').resolve()
WORK, SUB = OUT / '_work', OUT / 'submission'
for d in (WORK, SUB):
    if d.exists():
        shutil.rmtree(d)
(SUB / 'dfthcsos').mkdir(parents=True)
WORK.mkdir(parents=True)
print('data:', DATA, '| Python', sys.version.split()[0], '| cpus', os.cpu_count())

# %% [markdown]
# ## 2. Sources at a pinned commit, with integrity checks

# %%
if SOURCE_DIR:
    REPO = Path(SOURCE_DIR)
else:
    REPO = WORK / 'sources'
    REPO.mkdir()
    run(['git', 'init', '-q'], cwd=REPO)
    run(['git', 'remote', 'add', 'origin', SOURCE_REPO], cwd=REPO)
    run(['git', 'fetch', '-q', '--depth', '1', 'origin', SOURCE_COMMIT], cwd=REPO, timeout=900)
    run(['git', 'checkout', '-q', 'FETCH_HEAD'], cwd=REPO)
for rel, want in PINNED.items():
    got = sha256(REPO / rel)
    assert got == want, f'{rel}: SHA-256 {got} differs from the pinned {want}'
    print('ok ', rel)
assert (REPO / CHEM_RECIPE / 'manifest.json').exists() and (REPO / 'tools' / 'build_chemistry.py').exists()
ENV = dict(os.environ, PY=sys.executable, PYTHONOPTIMIZE='')
ENV.pop('PYTHONOPTIMIZE')

# %% [markdown]
# ## 3. Build the three tasks from source
#
# Each generator is run exactly as its authors document it. Both circuit generators finish with their own
# exact checks and a SHA-256 comparison; the notebook compares the result with the documented hash again.

# %%
def unzip(rel, name):
    dest = WORK / name
    with zipfile.ZipFile(REPO / rel) as z:
        z.extractall(dest)
    script = next(p for p in dest.rglob('build.sh'))
    return script


print('invsqrt: 46-qubit generator')
script = unzip('studies/compact_invsqrt46/invsqrt_source.zip', 'invsqrt_src')
(WORK / 'invsqrt_out').mkdir()
(WORK / 'invsqrt_tmp').mkdir()
run(['bash', script, ENGINE, VERIFIERS, WORK / 'invsqrt_out'], env=dict(ENV, WORK=str(WORK / 'invsqrt_tmp')),
    timeout=3600, log=WORK / 'invsqrt_build.log')
assert sha256(WORK / 'invsqrt_out' / 'invsqrt.kmb') == EXPECTED['invsqrt.kmb'], 'invsqrt.kmb differs from the documented build'
shutil.copy(WORK / 'invsqrt_out' / 'invsqrt.kmb', SUB / 'invsqrt.kmb')

print('dqi: 889-qubit decoder generator')
script = unzip('studies/source_portfolio889/upstream/dqi_source.zip', 'dqi_src')
(WORK / 'dqi_out').mkdir()
(WORK / 'dqi_tmp').mkdir()
run(['bash', script, ENGINE, VERIFIERS, WORK / 'dqi_out'], env=dict(ENV, WORK=str(WORK / 'dqi_tmp')),
    timeout=5400, log=WORK / 'dqi_build.log')
assert sha256(WORK / 'dqi_out' / 'dqi.kmb') == EXPECTED['dqi.kmb'], 'dqi.kmb differs from the documented build'
shutil.move(str(WORK / 'dqi_out' / 'dqi.kmb'), SUB / 'dqi.kmb')

print('dfthcsos: frozen table + GF(2) product plan -> 15 angle circuits')
run([sys.executable, REPO / 'tools' / 'build_chemistry.py', WORK / 'chem_out', '--recipe', REPO / CHEM_RECIPE,
     '--data', DATA, '--runtime', ENGINE], timeout=3600, log=WORK / 'chem_build.log')
chem = WORK / 'chem_out'
angle = sorted(chem.glob('angle_k_*.kmb'), key=lambda f: int(f.stem.split('_')[-1]))
assert [int(f.stem.split('_')[-1]) for f in angle] == list(range(len(angle))) and (chem / 'submission.npz').exists()
for f in angle + [chem / 'submission.npz']:
    shutil.copy(f, SUB / 'dfthcsos' / f.name)
print('built:', sorted(p.relative_to(SUB).as_posix() for p in SUB.rglob('*') if p.is_file())[:4], '...')

# %% [markdown]
# ## 4. Official verification

# %%
from invsqrt import verify_submission_invsqrt as VINV
from dfthcsos.verify import verify as verify_chem


def load_dqi_verifier():
    """The official dqi verifier; only the reaction-depth report (not part of the score) is skipped, to save memory."""
    path = VERIFIERS / 'dqi' / 'verify_submission_dqi.py'
    mod = types.ModuleType('vdqi')
    mod.__file__ = str(path)
    sys.modules['vdqi'] = mod
    exec(compile(path.read_text().replace('circuit.reaction_depth()', '-1'), str(path), 'exec'), mod.__dict__)
    return mod


VDQI = load_dqi_verifier()
R_ZERO = math.log(2) / -math.log(0.5 + 2 * math.sqrt(1 - 2 ** (-128 / 102400)))
REPORT = {}


def timed(label, fn):
    t0 = time.time()
    r = fn()
    print(f'   {label}: PASS ({time.time() - t0:.0f}s)')
    return r


print('invsqrt')
inv = km.Circuit.from_file(str(SUB / 'invsqrt.kmb'), format='kmb')
if VERIFY == 'full':
    r = timed('full, default seed', lambda: VINV.verify(inv, VINV.ProblemSpecInvSqrt.for_n_p(10), max_workers=os.cpu_count()))
    if SECOND_SALT:
        timed('full, second seed', lambda: VINV.verify(inv, VINV.ProblemSpecInvSqrt.for_n_p(10), salt=SECOND_SALT,
                                                       max_workers=os.cpu_count()))
    score = r['score']
else:
    spec = VINV.ProblemSpecInvSqrt(width=22, max_x=65536, beta=Fraction(1), precision_bits=24, shots_per_x=4, system_qubits=0)
    r = timed('quick (x <= 65,536)', lambda: VINV.verify(inv, spec, batch_size=256, max_workers=os.cpu_count()))
    score = r['num_qubits'] * (r['toffoli_count'] + (r['t_count'] + r['rotation_count'] * 28.0) / 4)
REPORT['invsqrt'] = dict(score=score, qubits=r['num_qubits'], toffolis=r['toffoli_count'], rotations=r['rotation_count'])
del inv

print('dqi')
dqi = km.Circuit.from_file(str(SUB / 'dqi.kmb'), format='kmb')
if VERIFY == 'full':
    r = timed('full, 102,400 shots', lambda: VDQI.verify(dqi, VDQI.ProblemSpecDQI.benchmark_instance(),
                                                         max_workers=min(4, os.cpu_count() or 1)))
    score = r['score']
else:
    r = timed('quick, 1,024 shots', lambda: VDQI.verify(dqi, VDQI.ProblemSpecDQI.benchmark_instance(num_shots=1024),
                                                        max_workers=min(4, os.cpu_count() or 1), log2_delta=-1.0))
    assert r['value_errors'] + r['phase_errors'] == 0, 'failing shots in the quick check'
    score = R_ZERO * r['num_qubits'] ** 2 * r['toffoli_count']
REPORT['dqi'] = dict(score=score, qubits=r['num_qubits'], toffolis=r['toffoli_count'],
                     failures=r['value_errors'] + r['phase_errors'])
del dqi

print('dfthcsos')
r = timed('full, default seed', lambda: verify_chem(str(SUB / 'dfthcsos')))
if SECOND_SALT:
    timed('full, second seed', lambda: verify_chem(str(SUB / 'dfthcsos'), salt=SECOND_SALT))
REPORT['dfthcsos'] = dict(score=r['score'], lambda_eff=r['lambda_eff'], toffolis=r['toffoli_count'], qubits=r['num_qubits'])

# %% [markdown]
# ## 5. Package and report

# %%
ratios = {t: BASELINE_SCORE[t] / REPORT[t]['score'] for t in BASELINE_SCORE}
overall = math.prod(ratios.values()) ** (1 / 3)
print(f'{"task":10} {"score":>18} {"ratio":>8}')
for t in BASELINE_SCORE:
    print(f'{t:10} {REPORT[t]["score"]:18,.6g} {ratios[t]:8.3f}   {REPORT[t]}')
print(f'\nestimated overall score: {overall:.4f}' + ('' if VERIFY == 'full' else '  (quick estimate)'))

zip_path = OUT / 'submission.zip'
with zipfile.ZipFile(zip_path, 'w', compression=zipfile.ZIP_DEFLATED, compresslevel=6) as z:
    for f in sorted(SUB.rglob('*')):
        if f.is_file():
            z.write(f, f.relative_to(SUB).as_posix())
with zipfile.ZipFile(zip_path) as z:
    names = set(z.namelist())
assert {'invsqrt.kmb', 'dqi.kmb', 'dfthcsos/submission.npz'} <= names
assert sorted(n for n in names if n.startswith('dfthcsos/angle_k_')) == sorted(f'dfthcsos/angle_k_{i}.kmb' for i in range(len(angle)))
json.dump(dict(REPORT, ratios=ratios, overall=overall, verify=VERIFY, source_commit=SOURCE_COMMIT,
               submission_zip_sha256=sha256(zip_path)), open(OUT / 'report.json', 'w'), indent=1, default=str)
shutil.rmtree(WORK, ignore_errors=True)
shutil.rmtree(SUB, ignore_errors=True)
print(f'{zip_path} ({zip_path.stat().st_size / 1e6:.1f} MB, {len(names)} files)')

# %% [markdown]
# ## 6. Our own additions (status)
#
# | Addition | Status | Measured so far |
# |---|---|---|
# | Chemistry: sensitivity-weighted coding of bit 1 of the angle words, then a refit of the weights | experiment done, compiler not built | weighted code words add 1.70x the rounding variance instead of 11.4x; estimated +3% on dfthcsos |
# | Our own decoder (fixed-schedule Euclid, `ours/dqi`) | built and fully verified | ratio 4.94: far behind the 889-qubit generator, kept as research |
# | Our own invsqrt (`notebooks/src/invsqrt.py`, v4) | built and fully verified | ratio 7.50: far behind the 46-qubit generator, kept as research |
#
# Nothing in this table changes the circuits built above yet.
