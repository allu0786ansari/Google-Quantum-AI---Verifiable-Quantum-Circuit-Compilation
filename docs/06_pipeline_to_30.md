# Pipeline: from 24.5 toward 30+

Status: **draft for review, nothing implemented yet.** Each phase ends with a measured
result and a go/no-go decision before the next one starts.

## 1. Where we stand and what the targets require

Score = geometric mean of three ratios (baseline cost / our cost). Every task counts the
same, so a 20% gain on any task lifts the overall score by the same amount.

| Task | Current ratio | Cost formula | Current cost terms |
|---|---|---|---|
| invsqrt | 21.10 | Q × C | 62 qubits × 1,587.6 Toffoli-eq (1,289 Toffolis + 41 rotations ≈ 19%) |
| dqi | 65.11 | R(ε₊) × Q² × T | 1.19126 (floor, 0 failures) × 891² × 396,922 |
| dfthcsos | 10.70 | λ_eff × T × Q | 47.03 × 15,585 × 428 |
| **Overall** | **24.49** | product = 14,697 | |

What each target needs, multiplying the three ratios together:

| Target | Product of ratios | Gain over today | Example split |
|---|---|---|---|
| 25 (today's top) | 15,625 | ×1.06 | dqi −6% Toffolis |
| **30** | 27,000 | **×1.84** | each task ×1.23 |
| 40 | 64,000 | ×4.35 | each task ×1.63 |
| 50 | 125,000 | ×8.5 | each task ×2.04 |

**Honest assessment.** 30 is a stretch target that needs real gains on all three tasks.
40–50 would mean cutting the cost of every task by about 40–50% on top of a frontier that
many teams have already squeezed. dqi alone can't carry it: 840 of its 891 qubits hold the
syndrome input and can't be removed, so the Q² term can improve by at most ~11%. Unless
Phase 1 finds a much larger opening than the profiles below show, I would plan for 27–32
and treat anything above that as a bonus.

## 2. Phase 0: make the Kaggle notebook run correctly (blocker, do first)

### 2a. The kickmix errors: one cause

In the Kaggle run, `import kickmix` "succeeded" but printed `kickmix: None`. Python found a
**folder** named `kickmix` (the competition data's `kickmix/` directory) on the import path
and loaded it as an empty namespace package. That shadowed the real compiled module. Every
later error follows from this:

- `module 'kickmix' has no attribute 'Circuit'`: the frontier build failed and silently fell back;
- `... no attribute 'CircuitBuilder'`: our invsqrt fell back to the baseline;
- `... no attribute 'GF2Field'`, then `name 'TW_TI' is not defined`: the dqi cells crashed.

Fix:
1. `load_kickmix()`: import normally, and **reject the module unless it has
   `CircuitBuilder`, `Circuit` and `GF2Field`**. Otherwise drop it from `sys.modules`, find
   `kickmix*.so` under the competition data (and `/kaggle/input`), and load that file
   explicitly with `importlib.util.spec_from_file_location`. If there's a wheel instead,
   `pip install` it. If nothing works, stop with a clear message listing what was found.
2. Print the kickmix file path, version and Python version in the setup cell.
3. **No silent fallbacks.** A `REQUIRE_FRONTIER = True` switch makes the notebook stop if
   the frontier build fails, instead of quietly submitting a ~2.3 score.

### 2b. Data access through kagglehub

```python
ASSETS = kagglehub.dataset_download("ahoirg/qcc-reproducible-circuit-assets")
CHEM_PIPELINE = kagglehub.dataset_download("marculera/dfthc-compact-pipeline")
```

Both are used directly, with `/kaggle/input` and `LOCAL_*` folders kept as fallbacks. The
release-2 manifest check stays: if the downloaded version isn't release 2, say so and stop.

**Gate:** the notebook runs on Kaggle end to end, prints the three frontier rows, and the
estimate is 24.49.

## 3. Phase 1: measurement harness (shared by all tasks)

Results can't be trusted without reproducible, fast measurement.

- One `score_<task>()` per task that returns the exact verifier score. Quick mode takes
  seconds; full mode is used only before submitting.
- Per-stage Toffoli/qubit profilers (stage tags in the circuit, conditioned gates counted
  by firing rate).
- An analytic chemistry cost model: call the verifier's own `heuristic_costs` for any
  shape (R, B, C, brot) without building circuits.
- A regression check: every candidate is compared with the frontier on the same inputs,
  and its SHA-256 and scores are logged in `results/`.

**Gate:** each frontier score is reproduced in under one minute (quick mode).

## 4. Phase 2: dqi (×1.15–1.45 plausible)

What the profile of the frontier circuit already shows:
- It already uses measurement-based uncomputation heavily (79k HMR, 809k conditioned CZ
  fix-ups, only 1.5k CCZ). That easy lever is already used.
- Toffolis by stage: **unextract 94.4k**, unconvert 53.2k, convert 52.6k, extract 53.2k,
  scan 44.3k, unsynthesis 38.2k, synthesis 31.5k.
- Uncomputing "extract" costs **1.8×** computing it. That asymmetry is the first thing to
  explain and the most likely place for a cheap win.

Experiments, in order:
1. **Unextract asymmetry.** Find which sub-block makes unextract more expensive than
   extract, and whether it can become a mirror image of extract or use MBU (a measured
   phase fix-up instead of a full recompute). Upside up to about −40k Toffolis (≈10%).
2. **Field-multiplier cost.** Count GF(2¹²) multiplies per stage and the AND-gate count of
   each. Swap in the cheapest bilinear formula: we already have a 42-term GF(4)-CRT
   multiply and a 54-term tower multiply, versus 63 in the reference. Upside depends on
   how many general multiplies remain; it's measured first.
3. **convert / unconvert basis changes.** Look for linear (CNOT-only) replacements and for
   shared sub-computations between the forward and reverse passes.
4. **Work qubits.** 51 qubits sit on top of the 840-qubit syndrome. Each one removed is
   worth about 0.22%, so 20 fewer gives ×1.05.

**Gate:** a candidate passes the official full verifier with 0 failures and costs less
than the frontier.

## 5. Phase 3: dfthcsos (×1.15–1.4 plausible)

What the verifier's cost model shows for the frontier (R=9, B=143, C=39, brot=8, J=5):
- Toffolis: angle loaders **8,932 (57%)**, inner state prep 4,086, Givens 2,250 (set by
  brot), outer state prep 317.
- Qubits: 125 persistent (40 of them the angle words, J·brot) plus **303 temporary** from the
  inner state-prep lookup. The verifier picks that lookup's size itself, and a shape scan
  over R, B, C and brot keeps it near 270–305. **Q is essentially fixed at ~390–430**, so
  the levers are T and λ_eff.

Experiments:
1. **Cheaper angle loaders** (the largest term). Optimize the ANF / paired-product
   loaders: better pairing certificates, reordering rows so repeated deltas cost nothing,
   and choosing the J/K split. A 30% cut in loader Toffolis gives ×1.21 overall on this task.
2. **brot = 7 with joint refit.** Saves 300 Givens Toffolis, 5 persistent qubits and
   enables more words per loader. It needs the factorization re-optimized at 7 bits
   (marculera's compact pipeline) and must stay under 20 mHa, with E_gap above 3.914898.
   The frontier sits exactly on both limits today, so this is the riskiest step.
3. **Shape search** (R, B, C) with the analytic model first, then refit only the shapes
   that predict a lower λ·T·Q.

**Gate:** both residuals ≤ 20 mHa, gap and Λ constraints pass, and λ·T·Q falls below 313.7M.

## 6. Phase 4: invsqrt (×1.05–1.2 plausible)

- The 62-qubit width is already minimal for this circuit, and the error is at 0.9993 of
  tolerance, so there's no slack to trade.
- **Rotations are 19% of cost.** Fewer phase-kick rotations (say 41 → 30) also loosen the
  per-rotation precision (ε grows as rotations fall), which cuts the T cost of each one.
- Try a cheaper seed/Newton split (lookup size against Newton steps) under the same
  62-qubit cap, scored by Q × C.

**Gate:** full verifier PASS and Q × C below 98,433.

## 7. Phase 5: integrate, verify independently, submit

- Every improved component goes into the notebook behind a switch, with the frontier as
  the automatic fallback (with a loud message, never silent).
- Full official verification of each final circuit, run separately from the build, with
  SHA-256 pinned in `results/`.
- Submit, record the public score, and only then start the next round.

## 8. Order of work and expected outcome

| Step | Effort | Expected overall after step |
|---|---|---|
| Phase 0 fixes | small | 24.49 (on Kaggle, reliably) |
| Phase 1 harness | small–medium | 24.49 |
| Phase 2 dqi | large | 25.7–27.6 |
| Phase 3 chemistry | large | 27–30.5 |
| Phase 4 invsqrt | medium | 27.5–32 |

These ranges are hypotheses to test, not promises. Each phase starts by measuring where
the cost is, and is dropped quickly if the measured opening is small.
