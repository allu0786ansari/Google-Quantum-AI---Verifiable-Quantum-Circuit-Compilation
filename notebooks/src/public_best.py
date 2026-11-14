# %% [markdown]
# # VQCC: best-of-public portfolio (verified)
#
# Collects the best publicly shared circuits for the three tasks, checks every candidate with the
# **official verifiers**, keeps the cheapest verified one per task and writes `submission.zip`.
#
# Sources (all public, credited to their authors):
#
# | Task | Source | Author |
# |---|---|---|
# | invsqrt | notebook `leotang09/vqcc-t7421-submission` | leotang09 |
# | invsqrt, dqi | dataset `syed3000/vqcc-own-circuit-portfolio` and its notebook | Syed Asad Ali (Apache-2.0) |
# | dqi | notebook `atakanaldemir/qcc-888q-repair-and-native-chemistry-co-design` | atakanaldemir |
# | dfthcsos | notebook `maoxin123/qcc-v25-coded-word-chemistry` | maoxin123 |
# | fallback | dataset `heliosli/qcc-circuits-live` | heliosli |
#
# **How to run on Kaggle**
# 1. Add the competition as input (*Add input → Competition*). Turn **Internet on** in the notebook settings.
# 2. Run all cells. If a source cannot be downloaded, attach it by hand (*Add input → Notebooks / Datasets*):
#    everything under `/kaggle/input` is scanned as well.
# 3. Download `/kaggle/working/submission.zip` and upload it on the competition's Submit page.
#
# Nothing here trusts a file by its name: a circuit is used only after the official verifier passes it.

# %%
# ---- Settings ---------------------------------------------------------------
VERIFY = 'full'          # 'full': official verification of the chosen circuits (about 30-60 min on Kaggle)
                         # 'quick': reduced checks, scores are estimates (about 10 min)
INVSQRT_FULL_TOP = 2     # fully verify this many invsqrt candidates (best by the quick estimate) and keep the best
STOP_IF_BELOW = 26.0     # raise instead of writing a zip when the estimated overall score is lower than this
LOCAL_DATA = None        # only for running outside Kaggle: folder with verifiers/ and sample_submission/
EXTRA_INPUTS = []        # extra folders to scan for circuits (outside Kaggle)

SOURCES = [              # (kind, handle); kind is 'notebook' (its output files) or 'dataset'
    ('notebook', 'leotang09/vqcc-t7421-submission'),
    ('notebook', 'maoxin123/qcc-v25-coded-word-chemistry'),
    ('dataset', 'syed3000/vqcc-own-circuit-portfolio'),
    ('notebook', 'syed3000/vqcc-own-circuit-portfolio-notebook'),
    ('notebook', 'atakanaldemir/qcc-888q-repair-and-native-chemistry-co-design'),
    ('dataset', 'heliosli/qcc-circuits-live'),
]

# %% [markdown]
# ## 1. Setup

# %%
import fnmatch
import hashlib
import importlib
import importlib.util
import json
import lzma
import math
import os
import shutil
import sys
import time
import types
import zipfile
from fractions import Fraction
from pathlib import Path

import numpy as np

COMPETITION = 'google-verifiable-quantum-circuit-compilation'
MARKER = Path('verifiers') / 'invsqrt' / 'verify_submission_invsqrt.py'
BASELINE_SCORE = {'invsqrt': 2076774.23, 'dqi': 24439397322564.6, 'dfthcsos': 3355599588.3155375}


def walk_files(root, max_depth=8):
    """All files under root, following symlinked folders."""
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
    tried = []
    roots = [LOCAL_DATA, os.environ.get('VQCC_DATA')]
    try:
        import kagglehub
        roots.append(kagglehub.competition_download(COMPETITION))
    except Exception as ex:
        tried.append(f'kagglehub competition download failed: {type(ex).__name__}: {ex}')
    roots.append('/kaggle/input')
    for r in roots:
        if not r or not Path(r).expanduser().exists():
            continue
        r = Path(r).expanduser().resolve()
        tried.append(str(r))
        if (r / MARKER).exists():
            return r
        for f in walk_files(r):
            if f.name == MARKER.name and f.parent.name == 'invsqrt' and f.parent.parent.name == 'verifiers':
                return f.parents[2]
    raise FileNotFoundError('Competition data not found. Add the competition as input. Tried: ' + '; '.join(tried))


def load_kickmix(data_root):
    """The real kickmix module (a folder named kickmix in the data can shadow it)."""
    api = ('CircuitBuilder', 'Circuit', 'GF2Field')
    try:
        import kickmix
        if all(hasattr(kickmix, a) for a in api):
            return kickmix
    except ImportError:
        pass
    sys.modules.pop('kickmix', None)
    tried = []
    for so in [f for r in (data_root, '/kaggle/input') for f in walk_files(r) if fnmatch.fnmatch(f.name, 'kickmix*.so')]:
        try:
            spec = importlib.util.spec_from_file_location('kickmix', so)
            module = importlib.util.module_from_spec(spec)
            sys.modules['kickmix'] = module
            spec.loader.exec_module(module)
            if all(hasattr(module, a) for a in api):
                return module
        except Exception as ex:
            tried.append(f'{so}: {type(ex).__name__}: {ex}')
        sys.modules.pop('kickmix', None)
    raise ImportError(f'kickmix not usable (Python {sys.version.split()[0]}): ' + '; '.join(tried))


DATA = find_data_root()
VERIFIERS = DATA / 'verifiers'
SAMPLE = DATA / 'sample_submission'
km = load_kickmix(DATA)
sys.path.insert(0, str(VERIFIERS))
OUT = Path('/kaggle/working') if Path('/kaggle/working').exists() else Path('output_public_best').resolve()
OUT.mkdir(parents=True, exist_ok=True)
CACHE = OUT / '_cache'
SUB = OUT / 'submission'
for d in (CACHE, SUB):
    if d.exists():
        shutil.rmtree(d)
(SUB / 'dfthcsos').mkdir(parents=True)
CACHE.mkdir()
print('data:', DATA, '| kickmix:', km.__file__, '| Python', sys.version.split()[0], '| cpus', os.cpu_count())

# %% [markdown]
# ## 2. Collect candidates
#
# Downloads each source, then scans every input folder. Archives (`.zip`) and compressed circuits (`.xz`)
# are unpacked into a cache. A file is a candidate by its name only (`*invsqrt*`, `*dqi*`, a folder with
# `submission.npz` and `angle_k_0`); whether it is *used* is decided by the verifiers below.

# %%
def fetch(kind, handle):
    try:
        import kagglehub
        if kind == 'dataset':
            path = kagglehub.dataset_download(handle)
        elif hasattr(kagglehub, 'notebook_output_download'):
            path = kagglehub.notebook_output_download(handle)
        else:
            raise RuntimeError('this kagglehub has no notebook_output_download; attach the notebook output as input')
        print(f'  ok    {kind:8} {handle} -> {path}')
        return Path(path)
    except Exception as ex:
        print(f'  FAIL  {kind:8} {handle}: {type(ex).__name__}: {str(ex)[:150]}')
        return None


def sha256(path, limit=None):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(1 << 22), b''):
            h.update(chunk)
    return h.hexdigest()


WANTED = ('.kmb', '.kmx', '.npz', '.kmb.xz', '.kmx.xz')


def unpack(roots):
    """Extract circuit-bearing members of zips and .xz files below roots into CACHE; returns the new folders."""
    new = []
    for root in roots:
        for f in list(walk_files(root)):
            low = f.name.lower()
            try:
                if low.endswith('.zip') and f.stat().st_size < 4 << 30:
                    with zipfile.ZipFile(f) as z:
                        members = [m for m in z.namelist() if m.lower().endswith(WANTED)]
                        if members:
                            dest = CACHE / ('zip_' + hashlib.sha256(str(f).encode()).hexdigest()[:12])
                            if not dest.exists():
                                z.extractall(dest, members=members)
                                new.append(dest)
                elif low.endswith(('.kmb.xz', '.kmx.xz')):
                    dest = CACHE / ('xz_' + hashlib.sha256(str(f).encode()).hexdigest()[:12])
                    if not dest.exists():
                        dest.mkdir()
                        with lzma.open(f) as src, open(dest / f.name[:-3], 'wb') as dst:
                            shutil.copyfileobj(src, dst)
                        new.append(dest)
            except Exception as ex:
                print(f'  skip {f}: {type(ex).__name__}: {ex}')
    return new


print('downloading sources')
roots = [p for p in (fetch(k, h) for k, h in SOURCES) if p]
roots += [Path(p) for p in EXTRA_INPUTS]
if Path('/kaggle/input').exists():
    roots.append(Path('/kaggle/input'))
level = roots
for _ in range(3):                       # zips inside zips, at most three levels
    level = unpack(level)
    roots += level
    if not level:
        break

cands = {'invsqrt': {}, 'dqi': {}, 'dfthcsos': {}}
for root in roots:
    for f in walk_files(root):
        low = f.name.lower()
        if not low.endswith(('.kmb', '.kmx')):
            continue
        if SAMPLE in f.parents:
            continue
        if low.startswith('angle_k_0.'):
            if (f.parent / 'submission.npz').exists():
                key = sha256(f.parent / 'submission.npz')[:16] + sha256(f)[:16]
                cands['dfthcsos'].setdefault(key, f.parent)
        elif low.startswith('angle_k_'):
            continue
        elif 'invsqrt' in low:
            cands['invsqrt'].setdefault(sha256(f), f)
        elif 'dqi' in low:
            cands['dqi'].setdefault(sha256(f), f)
        elif f.stat().st_size < 2 << 30:       # unnamed circuit: classify by its registers
            try:
                regs = set(km.Circuit.from_file(str(f), format='kmx' if low.endswith('.kmx') else 'kmb').register_data)
            except Exception:
                continue
            if regs == {'x'}:
                cands['invsqrt'].setdefault(sha256(f), f)
            elif 'S0' in regs and 'E0' in regs:
                cands['dqi'].setdefault(sha256(f), f)
for task, d in cands.items():
    print(f'{task}: {len(d)} distinct candidate(s)')
    for key, p in d.items():
        print(f'   {key[:12]}  {p}')

# %% [markdown]
# ## 3. Official verifiers

# %%
from invsqrt import verify_submission_invsqrt as VINV
from dfthcsos.verify import verify as verify_chem


def load_dqi_verifier():
    """The official dqi verifier; only the reaction-depth report (not part of the score) is skipped, to save memory."""
    path = VERIFIERS / 'dqi' / 'verify_submission_dqi.py'
    src = path.read_text().replace('circuit.reaction_depth()', '-1')
    mod = types.ModuleType('vdqi')
    mod.__file__ = str(path)
    sys.modules['vdqi'] = mod
    exec(compile(src, str(path), 'exec'), mod.__dict__)
    return mod


VDQI = load_dqi_verifier()
R_ZERO = math.log(2) / -math.log(0.5 + 2 * math.sqrt(1 - 2 ** (-128 / 102400)))   # R for 0 failures in 102,400 shots


def load_circuit(path):
    return km.Circuit.from_file(str(path), format='kmx' if str(path).lower().endswith('.kmx') else 'kmb')


def invsqrt_quick(path):
    """Estimated full score from the official verifier on x <= 65,536 (4 shots each)."""
    c = load_circuit(path)
    spec = VINV.ProblemSpecInvSqrt(width=22, max_x=65536, beta=Fraction(1), precision_bits=24,
                                   shots_per_x=4, system_qubits=0)
    r = VINV.verify(c, spec, batch_size=256, max_workers=os.cpu_count())
    # r['t_count'] holds the T gates only; price each rotation at 28 T, about its cost on the full domain
    return r['num_qubits'] * (r['toffoli_count'] + (r['t_count'] + r['rotation_count'] * 28.0) / 4), r


def invsqrt_full(path):
    r = VINV.verify(load_circuit(path), VINV.ProblemSpecInvSqrt.for_n_p(10), max_workers=os.cpu_count())
    return r['score'], r


def dqi_static(path):
    c = load_circuit(path)
    q, t = c.num_qubits, c.operation_statistics().num_toffoli_like
    return R_ZERO * q * q * t, dict(num_qubits=q, static_toffolis=t)


def dqi_verify(path, full):
    c = load_circuit(path)
    spec = VDQI.ProblemSpecDQI.benchmark_instance() if full else VDQI.ProblemSpecDQI.benchmark_instance(num_shots=1024)
    r = VDQI.verify(c, spec, max_workers=min(4, os.cpu_count() or 1), **({} if full else {'log2_delta': -1.0}))
    r.pop('shot_operation_counts', None)
    fails = r.get('value_errors', 0) + r.get('phase_errors', 0)
    score = r['score'] if full else R_ZERO * r['num_qubits'] ** 2 * r['toffoli_count']
    if not full and fails:
        raise AssertionError(f'{fails} failing shots in the quick check')
    return score, r


def chem_verify(folder):
    r = verify_chem(str(folder))
    return r['score'], r


def attempt(label, fn, *args):
    t0 = time.time()
    try:
        score, info = fn(*args)
        print(f'   {label}: {score:,.6g}  ({time.time() - t0:.0f}s)')
        return score, info
    except Exception as ex:
        print(f'   {label}: REJECTED ({type(ex).__name__}: {str(ex)[:160]})  ({time.time() - t0:.0f}s)')
        return None, None

# %% [markdown]
# ## 4. Choose the best verified circuit per task

# %%
CHOSEN = {}

# --- dfthcsos: the official verifier is fast, so every candidate is fully verified
print('dfthcsos')
best = None
for key, folder in cands['dfthcsos'].items():
    score, info = attempt(f'{key[:12]} full', chem_verify, folder)
    if score and (best is None or score < best[0]):
        best = (score, folder, info, 'full')
if best:
    CHOSEN['dfthcsos'] = dict(score=best[0], path=best[1], check=best[3],
                              detail={k: best[2].get(k) for k in ('lambda_eff', 'toffoli_count', 'num_qubits', 'brot')})

# --- invsqrt: quick estimate for all, full verification for the best few
print('invsqrt')
ranked = []
for key, path in cands['invsqrt'].items():
    est, info = attempt(f'{key[:12]} quick', invsqrt_quick, path)
    if est:
        ranked.append((est, path, info))
ranked.sort(key=lambda r: r[0])
best = None
if VERIFY == 'full':
    for est, path, _ in ranked[:INVSQRT_FULL_TOP]:
        score, info = attempt(f'{path.name} full', invsqrt_full, path)
        if score and (best is None or score < best[0]):
            best = (score, path, info, 'full')
    for est, path, _ in ranked[INVSQRT_FULL_TOP:]:        # only if none of the top ones passed
        if best:
            break
        score, info = attempt(f'{path.name} full', invsqrt_full, path)
        if score:
            best = (score, path, info, 'full')
elif ranked:
    best = (ranked[0][0], ranked[0][1], ranked[0][2], 'quick estimate')
if best:
    CHOSEN['invsqrt'] = dict(score=best[0], path=best[1], check=best[3],
                             detail={k: best[2].get(k) for k in ('num_qubits', 'toffoli_count', 'rotation_count')})

# --- dqi: rank by qubits^2 x static Toffolis, then verify in that order until one passes
print('dqi')
ranked = []
for key, path in cands['dqi'].items():
    est, info = attempt(f'{key[:12]} static', dqi_static, path)
    if est:
        ranked.append((est, path, info))
ranked.sort(key=lambda r: r[0])
for est, path, _ in ranked:
    score, info = attempt(f'{path.name} {"full" if VERIFY == "full" else "quick"}', dqi_verify, path, VERIFY == 'full')
    if score:
        CHOSEN['dqi'] = dict(score=score, path=path, check='full' if VERIFY == 'full' else 'quick estimate',
                             detail={k: info.get(k) for k in ('num_qubits', 'toffoli_count', 'value_errors', 'phase_errors')})
        break

# %% [markdown]
# ## 5. Package and report

# %%
def put_circuit(src, dst):
    if str(src).lower().endswith('.kmb'):
        shutil.copy(src, dst)
    else:
        load_circuit(src).to_file(str(dst), format='kmb')


missing = [t for t in BASELINE_SCORE if t not in CHOSEN]
for t in missing:
    print(f'WARNING: no verified public circuit for {t}; the competition baseline is used (ratio 1.0)')
if 'invsqrt' in CHOSEN:
    put_circuit(CHOSEN['invsqrt']['path'], SUB / 'invsqrt.kmb')
else:
    shutil.copy(SAMPLE / 'invsqrt.kmb', SUB / 'invsqrt.kmb')
if 'dqi' in CHOSEN:
    put_circuit(CHOSEN['dqi']['path'], SUB / 'dqi.kmb')
else:
    shutil.copy(SAMPLE / 'dqi.kmb', SUB / 'dqi.kmb')
chem_src = CHOSEN['dfthcsos']['path'] if 'dfthcsos' in CHOSEN else SAMPLE / 'dfthcsos'
angle = sorted((f for f in Path(chem_src).iterdir() if f.name.startswith('angle_k_') and f.suffix in ('.kmb', '.kmx')),
               key=lambda f: int(f.stem.split('_')[-1]))
assert [int(f.stem.split('_')[-1]) for f in angle] == list(range(len(angle))), 'angle files are not contiguous'
for f in angle:
    put_circuit(f, SUB / 'dfthcsos' / (f.stem + '.kmb'))
shutil.copy(Path(chem_src) / 'submission.npz', SUB / 'dfthcsos' / 'submission.npz')

ratios = {t: (BASELINE_SCORE[t] / CHOSEN[t]['score'] if t in CHOSEN else 1.0) for t in BASELINE_SCORE}
overall = math.prod(ratios.values()) ** (1 / 3)
print(f'\n{"task":10} {"score":>18} {"ratio":>8}  {"check":15} source')
for t in BASELINE_SCORE:
    c = CHOSEN.get(t)
    print(f'{t:10} {(c["score"] if c else BASELINE_SCORE[t]):18,.6g} {ratios[t]:8.3f}  {(c["check"] if c else "baseline"):15} '
          f'{c["path"] if c else "sample_submission"}')
    if c:
        print(f'{"":10} {c["detail"]}')
print(f'\nestimated overall score: {overall:.4f}')
if overall < STOP_IF_BELOW:
    raise RuntimeError(f'estimated score {overall:.3f} is below STOP_IF_BELOW = {STOP_IF_BELOW}: a source is missing. '
                       'Attach the missing notebook outputs or datasets as inputs and run again.')

zip_path = OUT / 'submission.zip'
with zipfile.ZipFile(zip_path, 'w', compression=zipfile.ZIP_DEFLATED, compresslevel=6) as z:
    for f in sorted(SUB.rglob('*')):
        if f.is_file():
            z.write(f, f.relative_to(SUB).as_posix())
with zipfile.ZipFile(zip_path) as z:
    names = set(z.namelist())
assert {'invsqrt.kmb', 'dqi.kmb', 'dfthcsos/submission.npz'} <= names
assert sorted(n for n in names if n.startswith('dfthcsos/angle_k_')) == sorted(f'dfthcsos/angle_k_{i}.kmb' for i in range(len(angle)))
shutil.rmtree(CACHE, ignore_errors=True)
shutil.rmtree(SUB, ignore_errors=True)
print(f'{zip_path} ({zip_path.stat().st_size / 1e6:.1f} MB, {len(names)} files): upload this file on the Submit page')
json.dump({t: dict(score=c['score'], ratio=ratios[t], check=c['check'], source=str(c['path']), **c['detail'])
           for t, c in CHOSEN.items()} | {'overall': overall}, open(OUT / 'report.json', 'w'), indent=1, default=str)
