# %%
# ---- Settings ---------------------------------------------------------------
VERIFY = 'quick'      # 'none' | 'quick' (minutes) | 'full' (official verifiers, ~30-60 min)
USE_IMPROVED = {'invsqrt': True, 'dqi': False, 'dfthcsos': False}

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

# %%
import glob
import io
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


DATA = find_data_root()
VERIFIERS = DATA / 'verifiers'
SAMPLE = DATA / 'sample_submission'
try:
    import kickmix as km
except ImportError:
    sys.path.insert(0, str(DATA / 'kickmix'))
    import kickmix as km
sys.path.insert(0, str(VERIFIERS))

OUT = Path('/kaggle/working') if Path('/kaggle/working').exists() else Path('output').resolve()
SUB = OUT / 'submission'
if SUB.exists():
    shutil.rmtree(SUB)
(SUB / 'dfthcsos').mkdir(parents=True)
REPORT = {}
print('data:', DATA)
print('kickmix:', km.__file__)
print('output:', SUB)
