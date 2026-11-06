# Paper analysis: Khattar et al., arXiv:2610.23832

Notes on the companion paper for the Kaggle competition "Google Quantum AI – Verifiable
Quantum Circuit Compilation". Where the paper and the verifier code disagree, the code wins;
see [02_verifiers_findings.md](02_verifiers_findings.md).

## Scoring

Score = geometric mean over the three tasks of (baseline cost / our cost).

- Each task counts equally in multiplicative terms: halving any one task's cost raises the
  overall score by 2^(1/3) ≈ 1.26×.
- Submitting a task's baseline gives that task a ratio of exactly 1.0.
- All three tasks must be present and must pass.

| Task | Cost |
|---|---|
| invsqrt (fusion) | Q × C, C = Toffolis + T/4 |
| dqi (optimization) | R(ε₊) × Q² × T, R(ε) = ln2 / −ln(½ + 2√ε) |
| dfthcsos (chemistry) | λ_eff × T × Q, λ_eff = √(E_gap(2Λ − E_gap)) |

## kickmix rules (paper §II)

- Instructions: X, CX, CCX, SWAP, CSWAP; Z, CZ, CCZ, Z_POW, NEG; HMR (X-basis measure and
  reset), R (reset); BIT_STORE0/1, BIT_INVERT; PUSH_CONDITION / POP_CONDITION.
- R on a nonzero qubit or an unfixed HMR causes a random sign flip, so sloppy cleanup fails.
- CSWAP costs 1 Toffoli. Measurement-based uncomputation of an AND costs 0 Toffolis
  (HMR plus a classically conditioned CZ).
- Simulator: about 6×10¹⁰ instructions per second. Flame charts show where cost goes.

## Benchmark I: 1/√x phasing (fusion)

- A 22-qubit register `x`, x from 1 to 3,133,452. Inputs above that are never tested.
- Apply phase exp(2πi/√x), tolerance 2⁻²⁴/√x turns (24 bits relative).
- Phase is checked mod 1 turn: bits worth ≥ 1 turn are irrelevant, the ½- and ¼-turn bits
  are free Clifford gates, the ⅛-turn bit costs one T. Finer bits need priced rotations.
- No failures allowed: every x, 64 shots each (about 200.5M shots), all must pass.
- Baseline: Q = 260, 7,735 Toffolis, 1 T, 37 rotations. Worst error 2⁻²⁵·⁰⁷ (about one bit
  of slack).

Baseline Toffolis per forward pass (the reverse pass costs about the same):

| Stage | Toffolis |
|---|---|
| Normalize (leading zero pairs, shift) | 120 |
| Horner cubic, 17 bit | 1,071 |
| y₀² | 442 |
| x_norm · y₀² | 1,043 |
| 3 − z | 28 |
| p = y₀(3 − z) | 867 |
| Shift p ≪ s (CSWAPs) | 300 |

Main levers:
1. Drop the Newton step (about 60% of C with its uncompute). A finer piecewise polynomial
   (for example a cubic on 32 segments, about 2⁻²⁹ approximation error) could reach the
   tolerance on its own. Rough estimate: about 1,000 Toffolis forward instead of 3,870.
2. Fit coefficients to the fixed-point arithmetic, checking exhaustively (fast).
3. Reduce the denormalization shift (300 CSWAPs each way).
4. Cut qubits: Q weighs as much as C. The baseline keeps 30–40-bit intermediates.
5. Possibly fuse the last multiply with the phasing (rows ≥ ¼ turn are free).

## Benchmark II: DQI Reed–Solomon decoder (optimization)

- GF(2¹²) with α¹²+α⁶+α⁴+α+1, little-endian. m = 4095 positions, n = 70 syndromes,
  up to ℓ = 35 errors.
- Input: 840 syndrome qubits; record c of 49,140 classical bits.
- Output: apply (−1)^(c·bin(e(s))), keep syndromes, clean all workspace.
- Error weights follow the principal eigenvector of the DQI tridiagonal matrix: mean about
  31.05, weight 35 occurs about 2.8% of the time, so the maximum weight can't be dropped.
- Baseline: 1,893 qubits, 5,725,461 Toffolis. Synchronized extended Euclid (222 steps),
  then Chien search plus Forney with a field inversion per position.

Failure pricing (ε₊ is Clopper–Pearson at confidence 1 − 2⁻¹²⁸ over 102,400 shots):

| Failures | ε₊ | R |
|---|---|---|
| 0 | 8.7×10⁻⁴ | 1.19 |
| ~5 | 1.0×10⁻³ | 1.21 |
| ~20 | 1.4×10⁻³ | 1.25 |

A perfect decoder still pays about 1.19×, and rare failures cost very little on top.

Main levers:
1. Qubits (squared). The 840 syndrome qubits are a floor. Ideas: leaner Euclid state,
   Berlekamp–Massey, Bernstein–Yang divsteps, in-place syndrome → (σ, Ω).
2. Field representation: a basis change is free (CNOTs). A tower field GF((2⁶)²) makes
   inversion much cheaper.
3. A better GF(2¹²) multiplier than the 63-Toffoli Karatsuba.
4. Use the free failure budget (for example, handle only small quotient degrees).
5. Classical processing of the record c is free; adding a Reed–Solomon dual codeword to the
   mask changes the phase only by a free linear function of the syndromes, so masks can be
   zeroed at chosen positions (small saving, about 1.7%).

## Benchmark III: DF-THC + spectrum amplification (chemistry)

- FeMoCo, N = 76 orbitals, η = 113 electrons.
- We choose R, B, C; the vectors and weights; the Givens tree; the angle precision brot;
  the batch size J; and the angle-loading circuits.
- The symmetry shift (S, β₁) is fixed to the baseline's.
- Gates: one- and two-body residuals ≤ 20 mHa (computed from the loaded, quantized words),
  Λ > 49.0, E_gap > 3.914898 Ha.
- λ_eff ≈ √(2·Λ·E_gap) when Λ ≫ E_gap. At both floors λ_eff ≈ 19.19 Ha.
- E_SOS uses the negative eigenvalues of h1′ (Eq. 57 in the paper has a sign typo; the code
  and Eq. A10 agree).

Main levers:
1. Classical optimization of the factorization, quantization-aware, sweeping R, B, C.
2. Tree shape (shallow trees compound rounding error less).
3. Lookup circuit design.

## Other notes

- Google committed SHA-256 hashes of stronger internal results.
- Related reading: Gu, Kemeny, Yelin & Lukin (arXiv:2610.22367), LLM agents on the same
  problem types; the ecdsa.fail "record-and-replay arithmetic" paper.
- Dates: Milestone prize 20 Nov 2026 (needs a public notebook); Longest Leader cutoff
  4 Jan 2027; entry and merger deadline 13 Jan; final submission 20 Jan 2027.
