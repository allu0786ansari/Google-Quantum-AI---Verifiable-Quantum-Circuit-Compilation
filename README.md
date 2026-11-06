# Google Quantum AI – Verifiable Quantum Circuit Compilation

Working repo for the Kaggle competition. The goal is to beat three reference quantum
circuits (scored by geometric mean of baseline cost / our cost) while passing Google's
verifiers.

## Docs

- [docs/01_paper_analysis.md](docs/01_paper_analysis.md): the companion paper
  (arXiv:2610.23832), task contracts and main levers.
- [docs/02_verifiers_findings.md](docs/02_verifiers_findings.md): what the verifier code
  actually checks and charges, plus baseline cost breakdowns.
- [docs/03_project_proposal.md](docs/03_project_proposal.md): problem statement, what makes
  each task hard, and the solution approach.
- [docs/04_implementation_plan.md](docs/04_implementation_plan.md): workstreams, checks,
  notebook design and schedule to the 20 Nov milestone.
- [docs/05_results_log.md](docs/05_results_log.md): verified results for each version.

## Baselines (all reproduced locally)

| Task | Cost formula | Baseline | Score |
|---|---|---|---|
| invsqrt | Q × (Toffolis + T/4) | Q 260, 7,735 Toffolis, 1 T, 37 rotations | 2,076,774 |
| dqi | R(ε₊) × Q² × T | Q 1,893, T 5,725,093/shot, R 1.19 | ≈ 2.44×10¹³ |
| dfthcsos | λ_eff × T × Q | λ_eff 49.30, T 93,237, Q 730 | 3.356×10⁹ |

`results/dfthcsos_baseline.json` holds the full dfthcsos verifier output.

## Status

| Task | Best score ratio | Notes |
|---|---|---|
| invsqrt | **7.50** | v4: 128 qubits, 1,914 Toffolis, score 276,794 (full verify PASS) |
| dqi | 1.00 | baseline |
| dfthcsos | 1.00 | baseline |

Estimated overall score with the current notebook: **1.958** (= 7.50^(1/3)).

## Submission notebook

`notebooks/vqcc_submission.ipynb` builds every circuit from source, checks it and writes
`submission.zip`.

- **On Kaggle:** add the competition data to the notebook and run all cells. The data is
  found with `kagglehub.competition_download('google-verifiable-quantum-circuit-compilation')`
  (falling back to anything mounted under `/kaggle/input`). The zip lands in
  `/kaggle/working`.
- **Locally:** set `LOCAL_DATA = 'data/'` in the optional cell, pointing at the folder where
  you unpacked the competition data (it can be the folder holding `verifiers/` or any folder
  above it). Install kickmix first with `pip install kickmix`. The zip lands in
  `output/` next to the notebook.

Settings in the second cell:

- `VERIFY`: `'none'`, `'quick'` (a few minutes) or `'full'` (official verifiers on
  everything, about 30–60 minutes).
- `USE_IMPROVED`: per-task switch between our circuit and the baseline.

Any task whose checks fail falls back to its baseline circuit, so the zip is always valid.

The notebook is generated from `notebooks/vqcc_submission.py` (plain Python with `# %%`
cell markers). Edit that file, then regenerate:

```bash
python tools/py_to_ipynb.py notebooks/vqcc_submission.py notebooks/vqcc_submission.ipynb
```

The `VQCC_DATA` environment variable works as an alternative to `LOCAL_DATA`.

## Setup

```bash
pip install kickmix numpy scipy
```

The competition data (`verifiers/`, `sample_submission/`) is not committed here; download
it from Kaggle. Point the tools at it with `VERIFIERS=/path/to/verifiers`.

- `python -I tools/quick_invsqrt.py circuit.kmb [MAX_X]`: fast partial invsqrt check.
- `python -I tools/quick_dqi.py circuit.kmb [NUM_SHOTS]`: fast partial dqi check.
- `python -I tools/profile_marks.py circuit.kmb [DEPTH]`: Toffolis per named stage.

Regenerate the reference circuits:

```bash
python verifiers/invsqrt/nr_reference_circuit.py invsqrt.kmb
python verifiers/dqi/dqi_reference_circuit.py dqi.kmb
```

## Key dates

- 20 Nov 2026: Milestone prize (top public leaderboard, public notebook required).
- 4 Jan 2027: Longest Leader prize cutoff.
- 13 Jan 2027: entry and team merger deadline.
- 20 Jan 2027: final submission deadline.
