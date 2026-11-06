"""Convert a '# %%' cell-marked Python file into a Jupyter notebook.

    python tools/py_to_ipynb.py notebooks/vqcc_submission.py notebooks/vqcc_submission.ipynb
"""
import sys

import nbformat

src, dst = sys.argv[1], sys.argv[2]
cells, kind, buf = [], None, []


def flush():
    text = '\n'.join(buf).strip('\n')
    if kind == 'markdown':
        text = '\n'.join(line[2:] if line.startswith('# ') else line.lstrip('#') for line in text.splitlines())
        cells.append(nbformat.v4.new_markdown_cell(text))
    elif kind == 'code' and text:
        cells.append(nbformat.v4.new_code_cell(text))


for line in open(src, encoding='utf-8').read().splitlines():
    if line.startswith('# %%'):
        flush()
        kind, buf = ('markdown' if '[markdown]' in line else 'code'), []
    else:
        buf.append(line)
flush()
nb = nbformat.v4.new_notebook(cells=cells)
nb.metadata['kernelspec'] = {'name': 'python3', 'display_name': 'Python 3', 'language': 'python'}
nbformat.write(nb, dst)
print(f'{dst}: {len(cells)} cells')
