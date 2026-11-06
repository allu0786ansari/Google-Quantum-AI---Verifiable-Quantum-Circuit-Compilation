# %% [markdown]
# ## 4. dfthcsos: chemistry block encoding
#
# Baseline for now: the sample submission's 10 angle circuits and factorization.

# %%
for f in sorted((SAMPLE / 'dfthcsos').iterdir()):
    shutil.copy(f, SUB / 'dfthcsos' / f.name)
REPORT['dfthcsos'] = dict(version='baseline', files=len(list((SUB / 'dfthcsos').iterdir())))
if VERIFY in ('quick', 'full'):
    from dfthcsos.verify import verify as verify_dfthcsos
    r = verify_dfthcsos(str(SUB / 'dfthcsos'))
    REPORT['dfthcsos'].update(verifier_status=r['status'], score=r['score'])
print(REPORT['dfthcsos'])
