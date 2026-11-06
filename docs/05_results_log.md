# Results log

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
| v1 | tower-field inverse, divide and multiply | 1,866 | 4,239,733 | 1.7586×10¹³ | 1.39 | quick (1,024 shots, 0 failures); full: running |

v1 rewrites GF(2¹²) as the tower GF(((2³)²)²). Basis changes are free CNOTs.

| Operation | Reference Toffolis | Tower Toffolis | Uses |
|---|---|---|---|
| inverse | 315 | 77 | Forney, 4,095 per shot |
| divide | 378 | 131 | Euclid quotients, 444 per pass |
| controlled multiply | 75 | 66 | Euclid sweep, 16,204 per pass |

- Multiply formulas (6, 18 and 54 AND terms for GF(8), GF(64), GF(4096)) are found by
  solving a GF(2) linear system for the output wiring of Karatsuba product terms.
- Every intermediate value is bilinear in live registers, so clean-up is X-basis
  measurement plus classically conditioned CZs: zero Toffolis.
- The multiply uses no scratch qubits (forms computed in place, fan-out through the
  output bits), which keeps the qubit peak at 1,866.
- Google's reference decoder is imported unchanged; a builder wrapper redirects the four
  field operations, so the notebook carries only the new arithmetic.

## dfthcsos

Baseline only so far (score 3,355,599,588).
