# %%
# ---- Settings ---------------------------------------------------------------
VERIFY = 'quick'      # 'none' | 'quick' (minutes) | 'full' (official verifiers, ~30-60 min)
USE_IMPROVED = {'invsqrt': True, 'dqi': True, 'dfthcsos': False}
REQUIRE_FRONTIER = True   # stop with an error if the frontier circuits can't be built (instead of a much lower score)

# %% [markdown]
# ### Optional: run locally
#
# On Kaggle, leave this as `None`: the data is found with `kagglehub`.
# To run on your own machine, set `LOCAL_DATA` to the folder where you unpacked the
# competition data (the folder that contains `verifiers/` and `sample_submission/`,
# or any folder above it).

# %%
LOCAL_DATA = None     # e.g. 'data/'

# %% [markdown]
# ## 1. Setup
#
# Finds the competition data, imports kickmix and the official verifiers, and
# prepares the output folder. Data is looked up in this order:
#
# 1. `LOCAL_DATA` (the optional cell above), or the `VQCC_DATA` environment variable;
# 2. `kagglehub.competition_download('google-verifiable-quantum-circuit-compilation')`;
# 3. any copy mounted under `/kaggle/input`.
#
# kickmix is imported and checked: if the import picks up the data's `kickmix/` folder
# instead of the real module (the cause of the "no attribute CircuitBuilder" errors),
# the compiled `kickmix*.so` from the data is loaded directly.
#
# Two public datasets are downloaded with `kagglehub.dataset_download`:
# `ahoirg/qcc-reproducible-circuit-assets` (frontier circuits, section 1b) and
# `marculera/dfthc-compact-pipeline` (chemistry factorization pipeline).

# %%
import glob
import io
import itertools
import json
import math
import os
import shutil
import sys
import time
import zipfile
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path

import numpy as np


COMPETITION = 'google-verifiable-quantum-circuit-compilation'
MARKER = Path('verifiers') / 'invsqrt' / 'verify_submission_invsqrt.py'


def _data_root_under(folder) -> Path | None:
    """The folder holding verifiers/ (the given folder or one below it), if any."""
    folder = Path(folder).expanduser().resolve()
    if (folder / MARKER).exists():
        return folder
    for hit in sorted(folder.glob('**/' + MARKER.as_posix())):
        return hit.parents[2]
    return None


def find_data_root() -> Path:
    tried = []
    for folder in (LOCAL_DATA, os.environ.get('VQCC_DATA')):
        if folder:
            tried.append(str(folder))
            root = _data_root_under(folder) if Path(folder).expanduser().exists() else None
            if root:
                return root
    try:
        import kagglehub
        path = kagglehub.competition_download(COMPETITION)
        print('kagglehub:', path)
        tried.append(f'kagglehub ({path})')
        root = _data_root_under(path)
        if root:
            return root
    except Exception as ex:
        tried.append(f'kagglehub failed: {ex}')
    if Path('/kaggle/input').exists():
        tried.append('/kaggle/input')
        root = _data_root_under('/kaggle/input')
        if root:
            return root
    raise FileNotFoundError('Competition data not found. Tried: ' + '; '.join(tried) +
                            '. On Kaggle, attach the competition data; locally, set LOCAL_DATA.')


KICKMIX_API = ('CircuitBuilder', 'Circuit', 'GF2Field')


def _is_kickmix(module):
    return all(hasattr(module, name) for name in KICKMIX_API)


def _find_files(root, pattern, max_depth=6):
    """Files under root matching pattern, following symlinked folders (Path.glob('**') does not)."""
    import fnmatch
    root, hits = Path(root), []
    if not root.is_dir():
        return hits
    base = len(root.parts)
    for folder, dirs, files in os.walk(root, followlinks=True):
        if len(Path(folder).parts) - base >= max_depth:
            dirs[:] = []
        hits += [Path(folder) / f for f in files if fnmatch.fnmatch(f, pattern)]
    return sorted(hits)


def load_kickmix(data_root):
    """Return the real kickmix module.

    A plain `import kickmix` can pick up a *folder* named kickmix (such as the one in
    the competition data) as an empty namespace package, which then has no
    CircuitBuilder. So the imported module is checked, and if it is not the real one,
    the compiled kickmix*.so shipped with the data is loaded from its file path.
    """
    import importlib
    import importlib.util
    tried = []
    try:
        import kickmix
        if _is_kickmix(kickmix):
            return kickmix
        where = getattr(kickmix, '__file__', None) or list(getattr(kickmix, '__path__', []))
        tried.append(f'import kickmix gave {where} without {KICKMIX_API}')
    except ImportError as ex:
        tried.append(f'import kickmix: {ex}')
    sys.modules.pop('kickmix', None)
    roots = [Path(data_root)] + ([Path('/kaggle/input')] if Path('/kaggle/input').exists() else [])
    for so in [p for r in roots for p in _find_files(r, 'kickmix*.so')]:
        try:
            spec = importlib.util.spec_from_file_location('kickmix', so)
            module = importlib.util.module_from_spec(spec)
            sys.modules['kickmix'] = module
            spec.loader.exec_module(module)
            if _is_kickmix(module):
                return module
            tried.append(f'{so}: loaded but lacks {KICKMIX_API}')
        except Exception as ex:
            tried.append(f'{so}: {type(ex).__name__}: {ex}')
        sys.modules.pop('kickmix', None)
    for wheel in [p for r in roots for p in _find_files(r, 'kickmix*.whl')]:
        import subprocess
        res = subprocess.run([sys.executable, '-m', 'pip', 'install', '-q', '--no-deps', str(wheel)],
                             capture_output=True, text=True)
        importlib.invalidate_caches()
        try:
            import kickmix
            if _is_kickmix(kickmix):
                return kickmix
        except ImportError:
            pass
        sys.modules.pop('kickmix', None)
        tried.append(f'{wheel}: pip install -> {res.returncode} {res.stderr.strip()[-200:]}')
    folder = Path(data_root) / 'kickmix'
    listing = sorted(p.name for p in folder.iterdir()) if folder.is_dir() else 'missing'
    raise ImportError(f'kickmix not usable (Python {sys.version.split()[0]}). Tried: ' + '; '.join(tried) +
                      f'. Contents of {folder}: {listing}')


def download_dataset(handle):
    """kagglehub.dataset_download(handle), or None (with the reason printed)."""
    try:
        import kagglehub
        path = kagglehub.dataset_download(handle)
        print(f'{handle}: {path}')
        return Path(path)
    except Exception as ex:
        print(f'{handle}: kagglehub download failed ({type(ex).__name__}: {ex})')
        return None


DATA = find_data_root()
VERIFIERS = DATA / 'verifiers'
SAMPLE = DATA / 'sample_submission'
km = load_kickmix(DATA)
sys.path.insert(0, str(VERIFIERS))

# Public datasets: frontier circuit assets (section 1b) and the chemistry pipeline.
ASSETS_PATH = download_dataset('ahoirg/qcc-reproducible-circuit-assets')
CHEM_PIPELINE_PATH = download_dataset('marculera/dfthc-compact-pipeline')

OUT = Path('/kaggle/working') if Path('/kaggle/working').exists() else Path('output').resolve()
SUB = OUT / 'submission'
if SUB.exists():
    shutil.rmtree(SUB)
(SUB / 'dfthcsos').mkdir(parents=True)
REPORT = {}
print('data:', DATA)
print('kickmix:', km.__file__, '| version', getattr(km, '__version__', '?'), '| Python', sys.version.split()[0])
print('output:', SUB)
