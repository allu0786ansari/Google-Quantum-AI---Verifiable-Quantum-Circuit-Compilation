# Evaluation of the "30+ roadmap" against what has been measured

Baseline: the 24.94 public portfolio (heliosli's release, 18 files hash-pinned). It contains:

| Task | Ratio | Details |
|---|---|---|
| invsqrt | 21.25 | 62 qubits, 1,281.5 Toffolis |
| dqi | 65.22 | 892 qubits, 395,370 Toffolis |
| dfthcsos | 11.20 | (8, 143, 44), λ_eff 46.893, 14,966 Toffolis, 427 qubits |

All numbers below come from `docs/05_results_log.md`, which records the tools and the exact
measurements.

## 1. Two corrections to the roadmap's premises

- **Tasks 12, 18, 19 and 20 target our fallback code, not the submitted circuits.**
  - `DFT_BROT = 12` with refitting disabled and the invsqrt `InvDesign(k=8, frac=(27,19,9), ...)`
    are our own early fallbacks, scoring 1.00 and 7.50.
  - The submitted chemistry already uses 8-bit angles. The submitted invsqrt is a different
    62-qubit architecture (lookup-seeded Newton) at 21.25.
  - Tuning the fallbacks cannot beat the circuits we submit.
- **Task 8 assumes the dqi decoder is extended Euclid. It isn't.**
  - The frontier is a GF(2²⁴) confluent / Hermite rational-interpolation decoder.
  - The published reversible Berlekamp–Massey (Gu et al. 2026) costs 440k Toffolis at 1,837 qubits,
    4.7× worse than the frontier on the dqi cost.

## 2. Task-by-task status

| # | Task | Status | Measured result |
|---|---|---|---|
| 1 | Freeze baseline | done | 24.94 files hash-pinned in the notebook. The simulator reproduces dqi T (395,368.9 vs 395,370) and invsqrt T (1,281.47); exact chemistry metrics reproduced |
| 2 | Score calculator | done | `tools/harness.py` reproduces every official score |
| 3 | Experiment manager | partial | results log plus JSON files; no database (not limiting) |
| 4 | dqi profile | done | mean T per stage. Reverse minus forward: extract +35k, synthesis +10k, convert ≈0; scan 44.3k at its floor. Peak qubits in unextract (companion blocks) |
| 5 | GF(2¹²) multiplication | measured | GF(2²⁴) products 99–105 Toffolis, at best-known bilinear complexity (e.g. 7 GF(64) multiplies × 15) |
| 6 | Inversion / division | not tested in isolation | inverse12 = 66 and inverse24 = 190 Toffolis. Inversions are ~3% of dqi T, so the upside is ~<1% |
| 7 | Liveness / workspace | done | exact chunked simulation: all 891–892 qubits live together at 22 points, ≥ 889 at 186. MBU tails: 44 convertible Toffolis out of 482,802 |
| 8 | Berlekamp–Massey | answered by literature | 440k T at 1,837 Q (Gu et al.), far worse than the frontier |
| 9 | Chien search | done | 4,025 positions × 11 = 44,275 = frontier: at the floor (OR of 12 bits needs 11 ANDs) |
| 10 | Fused phase | already in frontier | phase bit only; no error word is built |
| 11 | Approximate decoding | done | each failure costs ~0.46% in R. Pruning rarely-firing gates: +0.3% on the 893 build, ~0.1% on the 892; larger failure budgets lose |
| 12 | 11-bit chemistry refit | not applicable | already at 8 bits. 7-bit re-rounding gives error_2 455 mHa |
| 13 | Batch size J | analysed | J = 5 is forced by the 128-qubit cap with the current loader (11 selection + 40 words + 77 cache). J = 6 needs a smaller cache, which means a bigger low dimension and more terms |
| 14 | Givens trees | not tested | blocked: any new words need new loaders, and our builder costs ~12% more loader Toffolis than the frontier certificates |
| 15 | Angle-loader compressor | done | dimension-count floor 194 terms vs frontier 205–219 (≤ +2.9% on chemistry). Public builder 249. Our annealing and flip-graph searches couldn't beat the frontier |
| 16 | Joint factorization + loader | done (bounds) | frontier weights KKT-optimal (0.23%); continuous retraining bound 0.2%; grid edits ~0; low-Λ shape (15, 57, 19) fits 2.4× worse |
| 17 | invsqrt local rewrites | not tested | five 7-wire rewrites exist (Kirill Tymchishin's matcher). Larger-window search never tried |
| 18–20 | invsqrt parameters, polynomial degree, multipliers | not applicable | target our 7.50 fallback, not the 62-qubit circuit |
| 21 | Normalization / rotations | analysed | rotations are 19% of C; ×10 margin saves ~19 Toffoli-equivalents (1.2% on invsqrt) |
| 22–24 | Portfolio, acceptance, packaging | done | notebook: hash-pinned sources, fallback chain, official verification, zip checks |

**Summary:**
- Of the roadmap's 24 tasks, 17 are done, measured or answered.
- 4 are not applicable (they target our fallbacks).
- 3 remain untested: 6, 14 and 17.
- None of the measured items has more than ~1% left, except the dqi unextract asymmetry (§3, R4).

## 3. The actual roadmap

Ordered by expected gain per effort. Each item has a kill criterion.

| # | Work | Expected overall | Effort | Kill if |
|---|---|---|---|---|
| R1 | **Exact local-rewrite engine** (Task 17, generalised): enumerate windows of ≤ 8 wires in the 62-qubit invsqrt and in dqi, resynthesize each with fewer Toffolis, prove equivalence exhaustively on the window (values and phase, measurement branches included), then run the official verifier | +0.1–1% | 1–2 sessions | < 10 Toffolis found in invsqrt and < 0.2% in dqi after the first pass |
| R2 | **dqi inversion variants** (Task 6): replace inverse12/inverse24 instances with addition-chain or tower variants at equal width, priced on the full decoder | ≤ +0.5% | 1 session | isolated saving < 5 Toffolis per call |
| R3 | **invsqrt margin retune** (Task 21): spend Toffolis where they enlarge the worst-case phase margin, to cheapen the 41 rotations | ≤ +0.4% | 1 session | net Q·C not lower |
| R4 | **dqi unextract companion redesign**: the only opening above 1%. Make the return chain (companion_return 21.8k, fused_chart 9.2k, normalize 6.3k) cost no more than the forward, at equal width | up to +3% | several sessions | prototype on the 893 source doesn't cut ≥ 10k Toffolis at equal Q. Needs the 892 decoder's source (ask its authors on Kaggle) to reach the submitted circuit |
| R5 | **Frontier-quality loader compiler** (Task 15): prerequisite for any new chemistry tables (Tasks 13, 14, 16) | enables later work; alone ≤ +1% | 2+ sessions | can't reach ≤ 210 terms per slot on the existing tables |

**Expected outcome if R1–R4 all succeed:** about +2–5%, i.e. roughly 25.4–26.2. No measured
opening supports 30. The ×1.75 product the 28–29 leaders imply would need an unpublished
structural idea (most plausibly a different dqi decoder design), which the evidence here doesn't
point to.


## Update: R1-R3 and Tasks 6/14 executed (7 Nov)

All negative; details in `05_results_log.md`:
- R1: 0 Clifford-reducible windows (invsqrt, and 99.6k sampled dqi Toffolis).
- R2: inversions are 3.0% of dqi T.
- Task 14: tree already depth 6.
- R3: rotation angles are already LP-optimal for the verifier's margin.

Open: R4 and R5 only.

- R4 / DQI-A (companion chain) closed after specification and pricing (`09_dqi_companion_spec.md`):
  local alternatives need the chain value again or exceed the qubit budget; the projective-interface
  redesign is net x0.60-0.96 on dqi.
