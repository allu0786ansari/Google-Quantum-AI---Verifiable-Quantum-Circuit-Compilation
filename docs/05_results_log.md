# Results log

Verified results for each task, newest first. "Quick" means the official verifier on a
subset (invsqrt: x ≤ 65,536, 4 shots each); "full" means exactly as the leaderboard runs it.

## invsqrt

Baseline: Q = 260, 7,735 Toffolis, 37 rotations, C = 7,987.6, **score 2,076,774**
(computed from the reference model's exact error margin).

| Version | Design | Q | Mean Toffolis | Worst error / tol | Est. score | Ratio | Check |
|---|---|---|---|---|---|---|---|
| v4 (searched) | k=8, frac (27,19,9), δ bits (13,14), pf 27, qubit cap 128 | 128 | ≈1,914 | 0.957 | ≈275,100 | ≈7.55 | full: running |
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

Design notes:

- An exact numpy model of the fixed-point arithmetic is checked on all 3,133,452 inputs;
  the circuit matches it bit for bit on sampled inputs.
- Coefficients are fitted per segment, then refined by least squares against the exact
  truncated arithmetic, then re-centered.
- Cubic designs (k=5, 6) were worse: more qubits and Horner work than the table saves.
- kickmix adders trade qubits for Toffolis automatically; capping qubits made Q×C worse.

## dqi

Baseline only so far.

## dfthcsos

Baseline only so far (score 3,355,599,588).
