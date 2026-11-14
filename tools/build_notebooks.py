"""Generate all notebooks from the cell-marked sources in notebooks/src/.

    python tools/build_notebooks.py

Sources (plain Python, '# %%' code cells and '# %% [markdown]' text cells):

    notebooks/src/intro.py      title of the merged notebook
    notebooks/src/setup.py      settings, optional LOCAL_DATA, data lookup, output folder
    notebooks/src/frontier.py   public frontier components (credited), built from ahoirg's assets
    notebooks/src/invsqrt.py    task 1
    notebooks/src/dqi.py        task 2
    notebooks/src/dfthcsos.py   task 3
    notebooks/src/package.py    zip + report (merged notebook only)

Outputs:

    notebooks/vqcc_submission.ipynb   setup + all three tasks + package (upload this one)
    notebooks/invsqrt.ipynb           setup + invsqrt only
    notebooks/dqi.ipynb               setup + dqi only
    notebooks/dfthcsos.ipynb          setup + dfthcsos only
    notebooks/vqcc_public_best.ipynb  best-of-public portfolio, verified (from notebooks/src/public_best.py)
"""
from pathlib import Path

import nbformat

ROOT = Path(__file__).resolve().parents[1] / 'notebooks'
SRC = ROOT / 'src'
TASKS = {
    'invsqrt': 'invsqrt (fusion): phasing by 1/sqrt(x)',
    'dqi': 'dqi (optimization): Reed-Solomon decoder',
    'dfthcsos': 'dfthcsos (chemistry): block encoding',
}


def cells_from(text):
    cells, kind, buf = [], None, []

    def flush():
        body = '\n'.join(buf).strip('\n')
        if kind == 'markdown':
            body = '\n'.join(l[2:] if l.startswith('# ') else l.lstrip('#') for l in body.splitlines())
            cells.append(nbformat.v4.new_markdown_cell(body))
        elif kind == 'code' and body:
            cells.append(nbformat.v4.new_code_cell(body))

    for line in text.splitlines():
        if line.startswith('# %%'):
            flush()
            kind, buf = ('markdown' if '[markdown]' in line else 'code'), []
        else:
            buf.append(line)
    flush()
    return cells


def src(name):
    return (SRC / f'{name}.py').read_text(encoding='utf-8')


def write(path, parts):
    nb = nbformat.v4.new_notebook(cells=[c for p in parts for c in cells_from(p)])
    nb.metadata['kernelspec'] = {'name': 'python3', 'display_name': 'Python 3', 'language': 'python'}
    nbformat.write(nb, path)
    print(f'{path.relative_to(ROOT.parent)}: {len(nb.cells)} cells')


write(ROOT / 'vqcc_submission.ipynb',
      [src('intro'), src('setup'), src('frontier'), *(src(t) for t in TASKS), src('package')])

write(ROOT / 'vqcc_public_best.ipynb', [src('public_best')])

for task, title in TASKS.items():
    intro = (f'# %% [markdown]\n# # {title}\n#\n'
             f'# Builds and checks only this task, writing it to the `submission/` folder.\n'
             f'# The same code runs inside the merged notebook `vqcc_submission.ipynb`,\n'
             f'# which builds all three tasks and writes `submission.zip`.\n')
    summary = ('# %% [markdown]\n# ## Summary\n\n# %%\n'
               'print(json.dumps(REPORT, indent=2, default=str))\n'
               "for f in sorted(SUB.rglob('*')):\n"
               '    if f.is_file():\n'
               "        print(f'{f.relative_to(SUB)}  {f.stat().st_size / 1e6:.1f} MB')\n")
    write(ROOT / f'{task}.ipynb', [intro, src('setup'), src('frontier'), src(task), summary])
