# Results log

## Public frontier assembly (6 Nov 2026, evening)

After researching the leaderboard (top ≈ 25), we build on the strongest public circuits
(Apache-2.0, credited in the notebook), from ahoirg's *QCC Reproducible Circuit Assets*
release 2. Each is checked with the official verifiers here:

| Task | Component | Q | Toffolis | Score | Ratio | Check |
|---|---|---|---|---|---|---|
| invsqrt | 62-qubit circuit + five proved rewrites (SHA `0ddf5dce…`) | 62 | 1,289.2 | 98,433 | **21.098** | **full: PASS** (200.5M shots, worst 0.99932 of tol) |
| dqi | 891-qubit decoder (SHA `fc3ba816…`) | 891 | 396,921.5 | 3.7538×10¹¹ | **65.106** | **full: PASS, 0 failures** |
| dfthcsos | C39 table + paired-product loaders | 428 | 15,585 | 3.1368×10⁸ | **10.698** | **full: PASS** (error_2 19.9997 mHa) |

Overall from fully verified components: **24.49** (21.098 × 65.106 × 10.698)^(1/3). The source bundle's build reproduced all 17 pinned output hashes.
Our own circuits (below) remain in the notebook as the fallback.

Verified results for each task, newest first. "Quick" means the official verifier on a
subset (invsqrt: x ≤ 65,536, 4 shots each); "full" means exactly as the leaderboard runs it.

## invsqrt

Baseline: Q = 260, 7,735 Toffolis, 37 rotations, C = 7,987.6, **score 2,076,774**
(computed from the reference model's exact error margin).

| Version | Design | Q | Mean Toffolis | Worst error / tol | Est. score | Ratio | Check |
|---|---|---|---|---|---|---|---|
| **v4 (searched)** | k=8, frac (27,19,9), δ bits (13,14), pf 27, qubit cap 128 | **128** | **1,914.0** | 0.957 | **276,794** | **7.50** | **full: PASS** |
| v3 (Booth) | k=8, frac (28,19,9), δ bits (13,14), pf 26 | 136 | 1,882 | 0.717 | 287,065 | 7.23 | quick |
| v2 | radix-2, k=7, frac (29,20,11), pf 26 | 142 | 2,122 | 0.795 | 333,980 | 6.22 | quick |
| v1 | radix-2, full-span subtraction | 156 | 2,478 | 0.887 | 421,437 | 4.93 | quick |

What each step changed:

- **v1:** piecewise quadratic (256-entry lookup) replaces the cubic guess + Newton step;
  5-stage binary normalization (flags are the binary shift amount, cleared by measurement
  with phase fixes that run on half the shots); kickmix barrel shifter for denormalization.
- **v2:** each multiplier row is subtracted only from the low bits it can affect, plus one
  borrow qubit; one controlled decrement fixes the high bits. Saves 14 qubits and
  about 350 Toffolis.
- **v4:** local search over all register widths (single and paired moves), scored by building
  the real circuit; then a 128-qubit cap so kickmix's adders use fewer carry qubits (−0.5%).
- **v3:** radix-4 Booth recoding of the multiplier. Half the rows means half the
  carry-propagating subtractions and half the truncation noise, which also allows narrower
  registers. Negative digits are free (complement, subtract, complement).

v4 full verification (`results/invsqrt_v4_full_verify.json`): PASS on all 3,133,452 inputs ×
64 shots (200,540,928 shots), 34 rotations at ε = 1.42×10⁻¹³, circuit SHA-256
`972d538d…5622f`. The notebook regenerates this exact file.

Design notes:

- An exact numpy model of the fixed-point arithmetic is checked on all 3,133,452 inputs;
  the circuit matches it bit for bit on sampled inputs.
- Coefficients are fitted per segment, then refined by least squares against the exact
  truncated arithmetic, then re-centered.
- Cubic designs (k=5, 6) were worse: more qubits and Horner work than the table saves.
- kickmix adders trade qubits for Toffolis automatically; capping qubits made Q×C worse.

## dqi

Baseline: Q = 1,893, about 5,725,093 mean Toffolis, R = 1.1913 (0 failures in 102,400),
**score 2.4439×10¹³**.

| Version | Change | Q | Toffolis | Est. score | Ratio | Check |
|---|---|---|---|---|---|---|
| **v2** | + GF(4)-CRT multiply (42 terms) | **1,866** | **3,840,178** | **1.5929×10¹³** | **1.534** | **full: PASS, 0 failures in 102,400** |
| v1 | tower-field inverse, divide and multiply | 1,866 | 4,239,733 | 1.7586×10¹³ | 1.39 | quick (1,024 shots, 0 failures) |

v2 full verification (`results/dqi_v2_full_verify.json`): PASS, 0 value and 0 phase errors in
102,400 shots, R = 1.1913, circuit SHA-256 `a8d5ba71…cf45`. The notebook regenerates this file.

v1 rewrites GF(2¹²) as the tower GF(((2³)²)²). Basis changes are free CNOTs.

| Operation | Reference Toffolis | Tower Toffolis | Uses |
|---|---|---|---|
| inverse | 315 | 77 | Forney, 4,095 per shot |
| divide | 378 | 131 (v1), **119** (v2) | Euclid quotients, 444 per pass |
| controlled multiply | 75 | 66 (v1), **54** (v2) | Euclid sweep, 16,204 per pass |

- Multiply formulas (6, 18 and 54 AND terms for GF(8), GF(64), GF(4096)) are found by
  solving a GF(2) linear system for the output wiring of Karatsuba product terms.
- Every intermediate value is bilinear in live registers, so clean-up is X-basis
  measurement plus classically conditioned CZs: zero Toffolis.
- The multiply uses no scratch qubits (forms computed in place, fan-out through the
  output bits), which keeps the qubit peak at 1,866.
- v2 multiplies with GF(2¹²) seen as degree-6 polynomials over GF(4): values at the 4
  points of GF(4), the leading coefficient and remainders mod 3 irreducible quadratics fix
  the product (CRT), so 14 GF(4) multiplies × 3 ANDs = 42 terms.
- Google's reference decoder is imported unchanged; a builder wrapper redirects the four
  field operations, so the notebook carries only the new arithmetic.

## dfthcsos

Baseline: brot 12, J 8, K 10, λ_eff 49.301, T 93,237, Q 730, **score 3,355,599,588**.

Pipeline (notebook section 4): factorization → angle words (two-body rounded, one-body rows
recomputed as eigenvectors) → circuits (largest J under the 128-qubit cap, reference lookup
builder) → official verifier. On the baseline factorization at brot 12 it reproduces the
sample submission byte for byte (all 10 circuits) and its exact score.

Findings so far:

| Experiment | Result |
|---|---|
| Plain rounding at brot 11 / 10 | error_2 25.0 / 40.0 mHa (limit 20) |
| Prune 4 of 60 vectors per rank, no retraining | error_2 261 mHa |
| J = 9 at brot 12 | needs 131 qubits with the standard lookup (cap 128) |
| brot 11: freeze rounded vectors, refit weights, 300 epochs | error_2 21.2 mHa, λ_eff 50.04 |
| same, +1,000 epochs (decaying rate) | plateau at **20.7 mHa** (limit 20), λ_eff 49.18, E_gap 3.9375 |

Weights alone can't absorb the 11-bit rounding error; the next step is giving the refit more
capacity (a few extra copies C, which add weights but no angles) or quantization-aware
training of the vectors themselves.

Training runs at about 1.2 s per epoch on this 2-core CPU; Google's baseline used 50,000
epochs, so structural changes (fewer vectors) need GPU training.


## Phase 0 (pipeline doc 06): Kaggle kickmix fix and data access

- Cause of the Kaggle errors: `import kickmix` picked up the data's `kickmix/` folder as an
  empty namespace package (`kickmix.__file__` was None), so `Circuit`, `CircuitBuilder` and
  `GF2Field` were missing and every later cell failed.
- Fix: `load_kickmix()` checks the imported module for the API and otherwise loads the
  compiled `kickmix*.so` from the data by file path (folder search follows symlinks); a
  wheel is pip-installed as a last resort; failure lists what was found.
- `kagglehub.dataset_download` for `ahoirg/qcc-reproducible-circuit-assets` and
  `marculera/dfthc-compact-pipeline` in the setup cell.
- `REQUIRE_FRONTIER = True`: the notebook stops with a clear error if the frontier circuits
  can't be built, instead of silently submitting a much lower score.
- Tests: (1) simulated shadowing (fake empty `kickmix` module) with the competition's
  `kickmix.so` -> all three frontier rows, estimate 24.494; (2) assets missing -> stops
  with the RuntimeError.


## Phase 1: scoring harness and measured openings

`tools/harness.py` reproduces every frontier score exactly (invsqrt 21.098 in 1 s quick,
dqi 65.107 in 6 s at 256 shots, chemistry 10.698 in 19 s, chemistry analytic model 10.698
in 0.3 s) and prints the overall score with the other two tasks at the frontier.

Openings measured (target: product of ratios ×1.75 to pass 29.5, ×1.84 for 30):

- **Chemistry, lambda_eff.** `lambda_eff = sqrt(E_gap (2 Lambda - E_gap))`. E_gap sits on its
  floor (3.914899), so lambda_eff ~ sqrt(Lambda). Lambda = 284.4 (LCU one-norm), floor 49:
  -20% Lambda -> chemistry x1.12.
- **Chemistry, shape.** Verifier model with the frontier loader cost scaled by output bits
  (same lambda, accuracy NOT checked): (9,143,C=30,brot 7) cost x0.84, (9,143,20,7) x0.78,
  (9,143,39,6) x0.82. At brot 8, C <= 20 doubles Q (the inner lookup jumps to k2 = 5).
- **Chemistry, Q.** Q = 125 persistent + 303 temporary (inner state-prep lookup, size picked by
  the verifier); stays ~390-430 for any frontier-like shape.
- **dqi, failure budget.** Error weights 21-35 carry 99.99% of shots; w <= 20 appears ~5.6 times
  in 102,400. R by failures: 0 -> 1.19126, 1 -> 1.19673, 5 -> 1.21311, 10 -> 1.22901
  (~0.46% per failure), so dropping rare-weight paths pays if it saves more Toffolis/qubits.
- **dqi, structure.** Already MBU-heavy (79k HMR, 809k conditioned CZ); Toffolis by stage:
  unextract 94.4k, unconvert 53.2k, convert 52.6k, extract 53.2k, scan 44.3k, unsynthesis
  38.2k, synthesis 31.5k. Source: 248-file q18 bundle in ahoirg's assets (dqi-source-bundle).


## Experiments toward 30 (round 1): cheap levers on the frontier

| Experiment | Result | Verdict |
|---|---|---|
| Chemistry: refit W (= `H2_I`) only, words/w fixed (error_2 is W-independent), SLSQP on lambda_eff with E_gap floor, analytic gradients (`tools/experiments/chem_w_refit.py`) | lambda_eff 47.025650 -> 47.025592; Lambda split h1 13.13 + h2 271.26 | frontier already W-optimal; no gain |
| Chemistry: loader structure | 3,131 paired-product terms (205-219 per slot) + 78 high-cache + 11 low-cache ANDs per circuit; words carry bit7 = bit6 so one qubit per word is a free ancilla | near the floor for this table; a GF(2) term search could save a few % of loader Toffolis (< 2% on chemistry) |
| invsqrt: rotations | 41 Z_POW angles already individually tuned (not powers of two); worst error 0.99932 of tolerance | margin x10 would save ~19 Toffoli-eq (1.2%) and needs better arithmetic; small |

Conclusion: the frontier is optimal on the inexpensive levers. Reaching 30 needs (a) a new
chemistry factorization (retraining w/u, e.g. brot 7 or smaller C, at the 20 mHa limit; GPU
hours) and/or (b) dqi decoder changes at source level (rare-weight trimming against the
0.46%-per-failure budget, the unextract stage).


## Experiments toward 30 (round 2): chemistry at frozen words, dqi source build

- **Fast exact chemistry metrics** (`tools/experiments/fastfit.py`): with the angle words frozen
  (so the 15 loaders are unchanged), error_2^2 = sum_kl (w_k^T K w_l)^2 - 2 sum_k w_k^T M w_k + ||T'||^2
  with 143x143 blocks K = (U U^T)^2 and M = (bb|b'b') + shift terms. Reproduces the verifier to all
  printed digits (error_2 19.99965 mHa, lambda_eff 47.025650, Lambda 284.3928); milliseconds per
  evaluation, JAX gradients.
- **KKT check** (`kkt.py`): the frontier (w, W) is a KKT point of min lambda_eff s.t. error_2 <= 20 mHa,
  E_gap >= floor (relative residual 0.17%). Multipliers: 0.36 lambda per mHa of error_2 (0.76% per mHa),
  5.7 lambda per Ha of gap.
- **Copy split / drop** (`refit.py --add/--drop`): exact splits (1/sqrt2) to C = 40, 42 and refits give
  lambda_eff 47.0349 / 47.0336 at error_2 19.995: no gain; each copy layer costs 0.56% of Toffolis.
- **7-bit angles** (`chem_reround.py`): plain re-rounding of the frontier words gives error_2 455 mHa
  (limit 20), lambda_eff 74.9: the words are tuned to 8 bits; 7 bits needs a new factorization.
- **Shape model**: C 39 -> 30 at the same lambda lifts the overall score only to 24.93; a new
  factorization would also need new loader certificates (the public generic builder costs ~7% more
  loader Toffolis than the frontier's).
  => Chemistry has ~1-3% overall upside; not worth GPU time now.
- **dqi source build**: ahoirg's `build_dqi.py` rebuilds the 893-qubit hybrid6 decoder from the
  248-file source bundle hash-exactly (fb2da4f2..., 482,802 static Toffolis) in 175 s, 2.5 GB peak.
  Stage plan (playbook/kit/decoder/decode.py): mask projection, initial linear maps, extract (jet
  kernel x289: forward 290, reverse 393 Toffolis at 43 private qubits; the reverse is costlier
  because of the qubit cap), synthesis, convert, scan, then the inverses.


## Experiments toward 30 (round 3): dqi

| Experiment | Result | Verdict |
|---|---|---|
| Generic measurement-based uncomputation: Toffolis in XOR-only qubit tails ending at a reset/end (891 binary, `tails.py`) | 0 CCX (14 tails; 79,399 lifetimes already end in HMR) | closed |
| Same with exact lifetimes: 893 source rebuilt with a builder proxy tagging every `free()` (`build_tagged.py`, `tails_free.py`) | 44 CCX convertible out of 482,802 static (121,376 HMR endings) | closed: MBU is already used everywhere it applies |
| Qubit relabelling: conservative static liveness (`static_live.py`) | peak 891 | inconclusive |
| Exact sampled occupancy: whole 891 circuit run in 20k-line chunks on verifier-style inputs (256 shots, 30 s, `chunk_live.py`) | all 891 qubits nonzero at once at 22 points; >= 889 at 186 points (extract and unextract) | Q = 891 is real; -2 qubits would need ~186 sites changed for +0.15% overall |
| Field arithmetic | GF(2^24) products 99-105 Toffolis (e.g. 7 GF(64) multiplies x 15), GF(2^12) inverse 66, cube 28 | at best-known bilinear complexity |

The 893 source build also differs materially from the 891 binary (CCZ 37,324 vs 1,483 static;
mean T 407,444 vs 396,922), so source-level work starts ~3% behind the submitted circuit.


## Experiments toward 30 (round 4): chemistry loaders and 8-bit words

| Experiment | Result | Verdict |
|---|---|---|
| Loader tensors (`chem_loaders/load.py`) | 15 x 85 x 35 mixed ANF tensors, density 0.45-0.50, all slices full rank; frontier certificates use 410-438 rank-1 pieces (205-219 paired terms); a dimension count gives ~335 pieces (~168 terms) as a floor for generic tensors | room in principle |
| Public builder `tensorqrom.shared` (3 seeds, slots 0 and 5) | 249 terms vs frontier 219 / 205 | public loaders cost ~12% more; a new factorization starts ~7% behind on chemistry Toffolis |
| Annealing on rank-2 groups (`anneal.py`) and Kauers-Moosbauer flip-graph walk on rank-1 pieces (`flipgraph.py`), both from the frontier certificates | no reduction (rank-2 flips raise ranks; only the low factor ever collides, 79 distinct values) | frontier certificates are a strong local optimum for these moves |
| Grid local search on the 8-bit two-body words (`gridmoves.py`): exact change of lambda_eff and error_2 for single +-1 word moves, Lagrangian with the KKT multiplier 0.36 lambda/mHa | 8 of 300 random moves improve, by ~1e-5 lambda | frontier words are grid-optimal; small edits can't buy anything |

A single-row word change alters a slot's tensor by one rank-1 piece (ANF of a point indicator is a
low (x) high product), so small edits would be cheap to patch, but there is nothing to gain from them.
Chemistry now needs a factorization with substantially lower Lambda (>= ~10-15%) to pay for
loaders built without the frontier's certificate optimizer.


## Headroom map (floors vs frontier)

| Component | Frontier | Floor | Max gain on the task | Overall |
|---|---|---|---|---|
| dqi scan (Chien zero tests) | 44,275 | 4,025 positions x 11 (OR of 12 bits needs 11 ANDs) = 44,275 | none | 0 |
| dqi qubits | 891 | ~845 (840 syndrome qubits + minimal scratch) | Q^2 x1.11 | +3.6% |
| dqi jet reverse vs forward | 389 vs 269 Toffolis x 289 calls | reverse = forward | -35k Toffolis (-8.7%) | +2.8% |
| dqi forward locator work | ~289 jets x ~3 GF(2^24) mult-equivalents (~2,600 GF(2^12)) | Berlekamp-Massey count ~2,400 GF(2^12) multiplies | small | ~0 |
| chemistry loaders | 205-219 terms per slot | dimension count: k rank-2 matrices spanning 35 slices of 15x85 needs k >= 194 | x1.029 | +1% |
| chemistry lambda_eff | Lambda 284 | Low et al. reach 180 under a looser accuracy test | up to x1.26 | up to +8% |
| invsqrt rotations | 41 x ~7.2 Toffoli-eq (19% of C) | margin-limited | ~x1.03 | +1% |

Stacking realistic shares of these (dqi x1.15, chemistry x1.2, invsqrt x1.03) gives about 28.
That matches the 28-29 leaders more plausibly than any single breakthrough.


## Public frontier update: heliosli's release (local GM 24.94)

heliosli's `the-dataset-is-public-now-enjoy` (v17), dataset `heliosli/qcc-circuits-live`, portfolio
`v812e4edf89e_y33_ay3`: invsqrt 21.2464 (62 qubits, 1,576.6 Toffoli-eq), dqi 65.2156 (892 qubits,
395,370 mean Toffolis), dfthcsos 11.198 ((8, 143, 44), lambda_eff 46.893, 14,966 Toffolis, 427 qubits);
GM 24.9415. The notebook now prefers these files (each pinned by SHA-256) and keeps ahoirg's build
as the fallback. Tested with a synthetic release in the same manifest format (including an `.xz`
member): files picked by hash, verified, packed. yoyobar's and SomeJoint's public notebooks (24.94 /
24.65) pack prebuilt files only; none contains decoder source or the 28-29 methods.


## Experiments toward 30 (round 5): dqi jet-kernel headroom and distribution-aware pruning

- **Where the qubit peak is** (`dqi/chunk_live_stage.py`, 893 source build, sampled every 2,000 lines):
  unextract/companion_normalize and companion_fused_chart reach 891 (true peak 893 inside chunks);
  jet reverse blocks (reverse multiply v / scale u / scale v) run at ~888; synthesis, convert and
  unconvert at ~881; scan at 867. The jet reverse has no free headroom: extra scratch raises Q
  (each qubit ~0.224% of cost, ~912 Toffolis at T = 407k).
- **Jet kernel anatomy** (playbook/kit/decoder/exp/jet24/pro_jet.py): forward = out-of-place GF(2^24)
  inverse (190) + HMR of the old u (phase debt deferred in a classical token) + v *= u^2 (99) + HMR
  of the old v. The reverse must restore u_old, v_old and erase u' = inv(u_old) by measurement, whose
  repayment is a phase of a linear function of an inverse (not free). That is why reverse (389) >
  forward (269). Potential if the reverse matched the forward: -34.7k Toffolis (8.5%); worth it at up
  to ~38 extra qubits.
- **Bit-sliced simulator** (`dqi_dce/sim.c`, `compile_ops.py`, `run_stats.py`): compiles the 63.2M-op
  circuit to arrays and simulates 4,096 shots in ~10 s with per-gate statistics; reproduces the
  verifier's mean Toffolis exactly (407,444.2) and checks syndromes restored / ancillas zero.
- **Distribution-aware dead-gate removal** (`analyze.py`, `sweep.py`; 102,400 natural shots + 4,096
  per error weight 14-35): 15,838 Toffolis never fire in 4,096 shots, but 13,889 of them fire somewhere
  in 102,400. Truly never firing in 192,512 shots: 1,200 ops = 788 Toffolis/shot (0.19%).
  Trade-off (remove by firing-probability/saving order):

  | expected failures / 102,400 | Toffolis saved | R | net dqi |
  |---|---|---|---|
  | 0 | 1,010 (0.25%) | 1.19126 | x1.0025 |
  | 0.5 | 1,256 (0.31%) | 1.19126 | x1.0031 |
  | 1 | 1,291 | 1.19673 | x0.9986 |
  | 10 | 1,442 | 1.22901 | x0.9727 |

  The lineage has already pruned its rare branches; at most +0.1% overall.


## Correction: where dqi's Toffolis really go (893 build, mean per shot, simulator + stage tags)

| Stage | Forward | Reverse | Reverse - forward |
|---|---|---|---|
| extract / unextract | ~64k (extract 51.8k + inverse24/12, multiply v) | ~99k (unextract 51.6k + companion_return 21.8k + companion_fused_chart 9.2k + companion_normalize 6.3k + reverse jet ~10k) | **+35k** |
| synthesis / unsynthesis | 32.2k | 41.9k | +9.8k |
| convert / unconvert | ~60k (incl. aligned_route 3.0k) | ~61k | ~+1k |
| scan | 44.3k (at floor) | — | — |

The pro_jet kernel (forward 269 / reverse 389) runs only ~30 times in this build (2,970 / 99
Toffolis in 'multiply v'), not 289 as in the older playbook decoder: its asymmetry is worth ~3.6k,
not 35k. The large asymmetry is the companion-based return chain in unextract
(exp/new/a01/companion.py, patched in exp/new/a14/decode_owner.py). Removing all of it would be
-35k Toffolis (8.6% of T): dqi x1.094, about +3% overall, if it can be done without more qubits
(unextract is also where the 893 peak sits).


## Confirmation test for chemistry retraining (negative)

Question: can lambda_eff fall at error_2 = 20 mHa and E_gap = floor if the two-body vectors U move
continuously (no 8-bit rounding, loader cost ignored)? This bounds what any retraining from the
frontier basin could deliver.

- Two-electron tensor low-rank factor: rank 786 at eigenvalue 1e-5 (0.07 mHa Frobenius truncation),
  rank 1,042 at 1e-6 (0.008 mHa); exact metrics reproduced (`chem_pareto/pareto_metrics.py`).
- Adam on an augmented Lagrangian diverges: a 1e-5 step on every parameter raises error_2 from 20 to
  27 mHa (extremely stiff quartic).
- Gradient projection with Newton restoration onto both active constraints (`chem_pareto/gproj.py`):
  monotone descent with constraints held to 1e-8:

  | iteration | lambda_eff | lower than frontier |
  |---|---|---|
  | 10 | 46.9866 | 0.041% |
  | 30 | 46.9733 | 0.111% |
  | 50 | 46.9566 | 0.147% |
  | 80 | 46.9408 | 0.180% (gains per 10 iterations: 0.04% -> 0.007%) |

Verdict: even the continuous upper bound is ~0.2-0.5%, before rounding to 8 bits and before the
~12% loader penalty of re-building circuits. Chemistry retraining in this basin is not confirmed, so
no GPU notebook.


## Checks on the 24.94 baseline (heliosli's release, all 18 files hash-matched)

| Component | Check | Result |
|---|---|---|
| dfthcsos (8, 143, 44) | exact fast metrics: lambda_eff 46.89275, error_2 19.99999998 mHa, error_1 19.9999998 mHa, gap 3.9148982 (floor); KKT test of (w, W) | relative KKT residual 0.23%: already optimal, both error limits active, no slack |
| invsqrt (62 qubits) | exhaustive bit-sliced simulation over all 3,133,452 x x 4 measurement seeds (18 s); sim mean 1,281.474 Toffolis | 29 conditioned gates never fire = 0.056 Toffolis/shot; no control constantly 1: nothing to prune |
| dqi (892 qubits) | 58.27M ops simulated; sim mean 395,368.9 (verifier 395,370); 8,081 Toffolis never fire in 4,096 shots | 7,736 of them fire within 102,400 natural shots; the remainder is ~0.1% of dqi. Low-weight strata (w = 14) fail by design: the decoder spends its failure budget on rare low weights |

Conclusion: the public 24.94 portfolio is at the optimum of every generic lever measured
(weights, words, uncomputation, dead gates, qubit reuse). Further gains need new algorithms in
specific kernels; the largest identified is the dqi unextract companion chain (<= +3% overall),
which can only be built from the 893-qubit source.


## Product reuse and the low-Lambda chemistry shape (both negative)

- **Repeated AND products in the 892 decoder** (`dqi_dce/cse.py`): CCX/CCZ on a control pair whose
  values are unchanged since an earlier product on the same pair: 276 ops, 269 Toffolis/shot (0.07%).
  Product reuse is already applied.
- **Low et al.'s shape (15, 57, 19)** (Lambda 180 in the paper under a CCSD(T) criterion) against the
  current shape (8, 143, 44), same Hamiltonian-derived seed (marculera's `df_init.py`, eigen-leaves,
  OMP bases) and the same 1,500 L-BFGS iterations on error_2:

  | shape | error_2 at 1,500 | slope at the end |
  |---|---|---|
  | (8, 143, 44) | 67.2 mHa | ~4.5 mHa / 100 iterations |
  | (15, 57, 19) | 164.2 mHa | ~2 mHa / 100 iterations |

  The low-Lambda shape fits 2.4x worse and more slowly with fewer parameters; it is not a realistic
  route to error_2 <= 20 mHa (rounding to 8 bits needs ~12-13 mHa continuous). Its verifier cost is
  also worse at small B (Q 662). Chemistry stays with the (8-9, 143, 39-44) family.
