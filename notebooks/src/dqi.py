# %% [markdown]
# ## 3. dqi: Reed-Solomon decoder
#
# Baseline for now: the reference decoder (1,893 qubits, about 5.73M Toffolis).
# It is copied from the sample submission when available, or rebuilt from the
# reference script (about 7 s; the file is about 890 MB before zipping).

# %%
def use_baseline_dqi():
    src = SAMPLE / 'dqi.kmb'
    if src.exists():
        shutil.copy(src, SUB / 'dqi.kmb')
    else:
        from dqi.dqi_reference_circuit import build_reference_circuit
        build_reference_circuit(4095, 70, 12).to_file(str(SUB / 'dqi.kmb'), format='kmb')


t0 = time.time()
use_baseline_dqi()
REPORT['dqi'] = dict(version='baseline', seconds=round(time.time() - t0, 1))
print(REPORT['dqi'])
