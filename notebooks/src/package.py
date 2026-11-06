# %% [markdown]
# ## 5. Package and report
#
# Writes `submission.zip` with the three tasks at the archive root and checks
# its layout. The ratio column is baseline cost / our cost (higher is better);
# the overall score is the geometric mean of the three ratios.

# %%
# Baseline costs: invsqrt from the reference model's exact margin; dqi = R(0 failures) x 1893^2 x mean Toffolis.
BASELINE_SCORE = {'invsqrt': 2076774.23, 'dqi': 24439397322564.6, 'dfthcsos': 3355599588.3155375}

zip_path = OUT / 'submission.zip'
with zipfile.ZipFile(zip_path, 'w', compression=zipfile.ZIP_DEFLATED, compresslevel=6) as z:
    for f in sorted(SUB.rglob('*')):
        if f.is_file():
            z.write(f, f.relative_to(SUB).as_posix())
with zipfile.ZipFile(zip_path) as z:
    names = set(z.namelist())
need = {'invsqrt.kmb', 'dqi.kmb', 'dfthcsos/submission.npz'}
angles = sorted(n for n in names if n.startswith('dfthcsos/angle_k_'))
assert need <= names, f'missing {need - names}'
assert angles == [f'dfthcsos/angle_k_{i}.kmb' for i in range(len(angles))], 'angle files not contiguous'

ratios = {}
for task, info in REPORT.items():
    if info.get('version', '').startswith('baseline'):
        ratios[task] = 1.0
    elif info.get('score') and BASELINE_SCORE.get(task):
        ratios[task] = BASELINE_SCORE[task] / info['score']
print(f'{"task":10} {"version":22} {"score":>16} {"ratio":>8}')
for task, info in REPORT.items():
    s = info.get('score')
    print(f'{task:10} {info["version"]:22} {(f"{s:,.0f}" if s else "-"):>16} {ratios.get(task, float("nan")):8.2f}')
overall = math.prod(ratios.values()) ** (1 / len(ratios))
print(f'\nestimated overall score: {overall:.3f}')
print(f'{zip_path} ({zip_path.stat().st_size / 1e6:.1f} MB), {len(names)} files')
