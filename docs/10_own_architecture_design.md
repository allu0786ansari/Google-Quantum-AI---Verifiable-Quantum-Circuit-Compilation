# Our own architectures, layer by layer (from the competition paper)

Source: Khattar, Rubin, Gidney, Broughton, Low & Babbush, "Leveraging verifiable structure in quantum
circuits to find more efficient compilations with agentic AI" (the competition paper), §III A–C, plus
our measurements (`05_results_log.md`). No code or circuit from the public notebooks is used. Each
benchmark is split into layers. Each layer has:
- a classical model, checked exhaustively or statistically before any circuit is built;
- a kickmix block with its own unit test;
- a budget of qubits and Toffolis.

Integration happens only after every layer passes its own test. The official verifier is the final
check.

**What to expect.** Our first own circuits will score below the public 24.94; our earlier own
designs reached 7.5, 1.5 and 1.0. The value of this track is a pipeline we fully understand and
control, where every layer can be replaced and measured. The targets below are what each design can
plausibly reach.

---

## Benchmark I: phase by 1/√x (score Q × C)

The paper's pipeline (§III A 2): normalise x by leading zero pairs → initial guess (cubic polynomial
from QROM coefficients) → one Newton step → shift by s → one Z_POW per bit of y1 → uncompute.
Reference: Q 260, C 7,988.

The authors name the open directions:
- a simpler polynomial over a larger range;
- more accurate QROM coefficients;
- "a completely different route": a polynomial over the whole domain.

| Layer | Paper baseline | Our choice | Why | Test |
|---|---|---|---|---|
| I-1 normalisation | s in unary (np = 10 qubits); 2-bit shift per stage | s in **binary** (4 bits), 4 conditional-shift stages (16/8/4/2) with flags computed by one AND tree per stage and cleared by measurement | unary s costs 6 extra qubits; each qubit is worth ~1.6% of score | all 3,133,452 x: x_norm ∈ [2^20, 2^22), s exact |
| I-2 initial guess | cubic per half-interval, coefficients from QROM, three multiplications | **direct seed table**: QROM on the top k bits of x_norm returns y0 with ~(k+1) correct bits; no polynomial | removes all seed multiplications; the Newton step roughly squares accuracy, so ~12–13-bit seeds suffice for 24 bits | exhaustive: worst relative error of y0 per k ∈ {8..12} |
| I-3 Newton step | p = y0(3 − x_norm·y0²)/2 at full width | truncated, error-budgeted multipliers: square y0, multiply by x_norm with only the bits the 2^-24 budget needs, then multiply by y0; Booth or schoolbook chosen per stage | multiplication is the dominant Toffoli cost; truncation widths come from an exhaustive error model | exhaustive model: worst relative error ≤ 0.9·2^-24 |
| I-4 shift and phase | y1 = p << s, then one Z_POW per bit | **no physical shift**: the phase of p·2^s is applied with s-controlled angle sets, or the shift is fused with the rotations | removes the shift workspace; trades against rotation count | exhaustive phase check against isqrt reference |
| I-5 uncompute | full reverse | measurement-based uncomputation (HMR + CZ fix-ups) for every product register; reverse only the in-place parts | halves the cleanup Toffolis where it applies | official verifier (all x, 64 shots) |

Target for our first full version: Q ≈ 70–90, C ≈ 1,800–2,500, ratio ≈ 9–16 (the public frontier is
62 / 1,281, ratio 21.2).

---

## Benchmark II: Reed–Solomon decoding as a phase oracle (score R·Q²·T)

The paper's contract (§III B 1–3): input syndromes s (70 × 12 = 840 qubits) and the classical record
c (49,140 bits). Apply (-1)^(c · bin(e(s))), restore s, clear workspace. Up to 35 errors over
GF(2^12); failures are allowed but priced by R(ε+), and every exact circuit pays R = 1.19.

Reference: EEA (Proos–Zalka synchronised, register-shared) → Chien + Forney with an Itoh–Tsujii
inversion per position → uncompute. Q 1,893, T 5.73M. The paper: "every part of this baseline—explicit
remainder storage, pointer-driven EEA steps, polynomial-basis F_2^12 arithmetic, and per-position
Itoh–Tsujii inversion—can be improved", and implementations "may … compute the required parity
directly".

| Layer | Paper baseline | Our choice | Why | Test |
|---|---|---|---|---|
| II-0 field arithmetic | 63-Toffoli Karatsuba multiply; Itoh–Tsujii inversion 315 Toffolis + 48 scaffold qubits | our **42-term GF(4)-CRT bilinear multiply** (already built and verified on 2,000 random pairs); tower inversion via norm to GF(2^6) (~66 Toffolis, ~12 scratch) | multiply 63 → 42, inversion 315 → ~66, scaffold 48 → ~12 | exhaustive on all field pairs/elements |
| II-1 mask projection | none | add a syndrome-space vector to the mask so ~70 positions get μ_j = 0 and are skipped by the scan; free (linear phase on s) | 70 × 11 Toffolis saved; zero cost | classical: phase identity over random (s, c) |
| II-2 key equation | synchronised EEA, 222 steps, pointer-driven slices | **scale-invariant Berlekamp–Massey**: no inversion per step (both branches share one update; the locator is only needed up to a scalar), fixed 70 iterations, controlled swap for the length change, in-place on a σ register | BM does ~2t² multiplications vs EEA's ~3n²; no pointer arithmetic | classical BM on 102,400 verifier-style syndromes against the true locators |
| II-3 scan + phase | per position: Chien update, Itoh–Tsujii inversion of σ', Forney, Z gates | Chien registers updated by constant multiplies (CNOTs); **zero test 11 Toffolis per position** (the floor); phase Tr(μ_j Ω/σ') through a precomputed reciprocal, so no per-position inversion | scan floor 4,025 × 11 = 44,275; no per-position inversion (removes 4,095 × 351) | phase identity on random shots |
| II-4 cleanup | uncompute EEA | reverse BM with MBU on every product register (phase of a bilinear form is free CZs) | compute + uncompute ~2×; MBU removes the products' share | official verifier, 102,400 shots |

Qubit plan: 840 syndrome + σ (36 × 12) + B (36 × 12) + discrepancy/scratch (~40) ≈ 1,750 before any
in-place trick. In-place reuse of consumed syndromes is the second-generation goal (σ can overwrite
syndrome words once no future discrepancy needs them).

Targets:
- first version: Q ≈ 1,750, T ≈ 450–600k: ratio ≈ 9–12;
- with in-place syndrome reuse (Q ≈ 1,200): ratio ≈ 20–25;
- the public frontier (GF(2^24) confluent interpolation, Q 892, T 395k) is 65.

---

## Benchmark III: DF-THC block encoding (score λ_eff · T_tot · Q_tot)

The paper's pipeline (§III C 2–7), with three user-controlled parts:
- the factorization (R, B, C; vectors u; weights w, W) under ε_fro ≤ 20 mHa (eq. 62), minimising
  λ_eff = √(E_gap(2Λ − E_gap));
- the Givens tree;
- the Boolean angle oracles (QROM), at most 128 qubits each.

State preparation and Givens costs are fixed formulas (`heuristic_costs`).

| Layer | Paper baseline | Our choice | Why | Test |
|---|---|---|---|---|
| III-1 factorization | (65, 60, 60), Λ 308, ε 18.6 mHa | train at (8–9, 143, ~40) (rows R·B ≈ 1,150–1,290, Q stays ~427) with the paper's objective (eq. 62) plus the cost model (eq. 63) in the loop; our exact fast evaluator (milliseconds per evaluation, JAX gradients) and gradient projection with restoration | the shape family that minimises the verifier's cost model (shape scans: small B blows Q up to 600+) | exact metrics vs verifier; continuous ε ≤ ~12–13 mHa before rounding |
| III-2 Givens tree | reference tree | shallow tree (depth ≤ 6, mean path ~3), chosen to minimise rounding error amplification | rounding error compounds along root paths | exact re-encode / decode error |
| III-3 angle words | 12-bit plain rounding | 8 bits, signed gauge (top bits equal, so 7 bits of information and a free ancilla per word), progressive quantization-aware rounding with refits between stages | 8-bit words with J = 5 per oracle fit the cap; plain rounding fails (we measured 455 mHa at 7 bits) | ε_fro and ε_1 after rounding ≤ 20 mHa |
| III-4 angle oracles | plain QROM, L − 1 Toffolis per circuit | ANF loader: split the 11 address bits 4 low / 7 high, cache high monomials (78 ANDs), up to 4 low composites, **paired-product terms** (one Toffoli per rank-2 piece) from our own GF(2) decomposition search | the 15 × 85 × 35 tensor admits ≥ 194 terms (dimension count); good decompositions reach ~205–220 | simulate every address; words exact |
| III-5 weights | — | final (w, W) refit at frozen words (KKT) | spends the remaining error budget on λ | verifier |

Targets: with our own factorization and loaders, ratio ≈ 8–10 (frontier 11.2). The loader compiler
is the critical layer: the public generic builder is ~12% worse than the frontier certificates.

---

## Build order

1. **Benchmark I first.** Smallest circuit, exhaustive checks in seconds, every layer independently
   testable. Order: I-1 → I-2 → I-3 → I-4 → I-5. Decision points: seed width k, multiplier widths.
2. **Benchmark II, layer by layer with classical models first.** II-0 (have), II-2 classical BM on
   verifier syndromes, II-3 classical scan/phase identity, then reversible blocks for each, then the
   in-place generation.
3. **Benchmark III last.** It needs training compute (GPU for the factorization) and the loader
   compiler.

Each layer gets a folder under `ours/<benchmark>/<layer>/` with its classical model, its kickmix
block, a unit test, and a recorded cost. The integrated circuit enters the notebook as the "ours"
tier and is used only if it beats the current component.
