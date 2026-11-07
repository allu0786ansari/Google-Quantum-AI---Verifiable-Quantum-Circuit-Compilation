# Literature review: algorithms for the three tasks, and where 30 could come from

Status: research complete (7 Nov 2026). Numbers are from the papers cited, the official
verifiers, and our own measurements (`docs/05_results_log.md`, rounds 1-3).

## 1. Bottom line

The circuits we submit (the public frontier, 24.49) are **ahead of every published result on
all three tasks**. No paper offers an algorithm that, transplanted as is, beats them. What the
literature does offer is a few directions with measurable upside, all on the order of a few
percent overall, plus one larger but uncertain lever in chemistry (the one-norm Λ).

| Task | Best published | Frontier (ours) | Gap |
|---|---|---|---|
| dqi | Gu et al. 2026 (LLM-discovered): 440k Toffolis, 1,837 qubits; cost ≈ 1.8×10¹² | 396,922 Toffolis, 891 qubits; cost 3.75×10¹¹ | frontier 4.7× better |
| invsqrt | no published circuit for this exact phase function; the closest comparable (Gu et al., fixed-point Gaussian) is 1,583 Toffolis / 139 qubits | 1,289 Toffolis (+41 rotations) / 62 qubits | frontier ahead |
| dfthcsos | Low et al. 2025 (spectrum amplification): DF-THC + BLISS, Λ = 179.7, λ_eff = 43.65 (under a looser accuracy criterion) | Λ = 284.4, λ_eff = 47.03, 15,585 Toffolis, 428 qubits | literature has lower Λ, but under a different accuracy test |

Realistic outlook from everything below: about 25.5–27 with the levers we can identify. A score of
30 needs the product of the three ratios to rise ×1.84, which would take an idea that isn't in
the literature yet.

## 2. dqi: Reed–Solomon decoding as a phase oracle

### What is published

| Work | Decoder | Toffolis | Qubits | Notes |
|---|---|---|---|---|
| Jordan et al., *Optimization by DQI* (Nature 2025) | any efficient decoder | — | — | defines the task |
| Khattar et al., *Verifiable quantum advantage via optimized DQI circuits* (arXiv:2510.10967) | synchronized reversible extended Euclid (Proos–Zalka style), in-place quotients; also an implicit "Dialog" (Bernstein–Yang) variant | 5.72M | 1,885 | the competition baseline; Chien + Forney with an inversion per position |
| Gu et al., *Discovering state-of-the-art quantum circuits with LLMs* (arXiv:2610.22367) | scale-invariant Berlekamp–Massey; one history entry plus a failure budget (≈10⁻⁴); phase bit only; reciprocal polynomial so no division in Chien | 440k (locator solve ~140k) | 1,837 | 12× over baseline in T·Q |
| *Quantum circuit design for DQI* (arXiv:2504.18334) | Gauss–Jordan / lookup | small toy instances | — | not applicable to RS(4095) |
| Ancilla-free GF(2^m) multiply/divide (arXiv:2511.20618) | Karatsuba + Itoh–Tsujii | m = 163: 3,605 Toffolis per multiply | 3m | arithmetic only |
| Sun–Wootters sampler for OPI (arXiv:2607.16541) | list decoding, prime fields | asymptotic only | — | not applicable |

### Why the frontier is hard to beat

- **Qubits:** 840 of 891 hold the syndromes. Exact chunked simulation shows all 891 nonzero at once
  at 22 points (≥ 889 at 186), so −1 qubit needs many sites changed; each qubit is worth 0.22%.
- **Uncompute:** measurement-based uncomputation is already used everywhere it applies: rebuilding
  from source with every `free()` tagged leaves 44 convertible Toffolis out of 482,802. A phase of a
  bilinear form is free (CZs), and the frontier already uses that to MBU every product register.
- **Arithmetic:** GF(2²⁴) products cost 99–105 Toffolis (for example 7 GF(64) multiplies × 15), at
  best-known bilinear complexity.
- **Scan:** Chien search needs a zero test per position (4,095 × ~11 = 44k). Mask projection already
  zeroes the masks at the ~70 positions the syndrome space allows (`SELECTED`). Testing fewer bits
  per position fails: 4,060 positions × 2⁻ᵏ false positives needs k ≈ 28.
- **Alternatives considered and rejected:**
  - Residue / extended-syndrome evaluation of the phase (Σⱼ μⱼ eⱼ as a residue of PΩ/σ): needs the
    LFSR run 4,025 steps × ~35 multiplies, ≈ 5.6M Toffolis.
  - Peterson–Gorenstein–Zierler: O(t³) multiplies.
  - Uncomputing the decoded error by measurement: needs the same kind of phase oracle again.

### Remaining dqi directions (each a few percent at most)

1. **Failure budget on rare weights.** Error weights 21–35 cover 99.99% of shots; w ≤ 20 appears
   ~5.6 times in 102,400 shots. R grows ~0.46% per failure (0 → 1.19126, 5 → 1.21311). Logic that
   exists only for low-weight or degenerate cases ("poles", early groups) can be dropped when it
   saves more than that. wahyurahmaniar's 892-qubit build is an example (−1 qubit, −440 Toffolis).
2. **Reverse jet kernel.** The reverse kernel costs 393 vs 290 forward × 289 calls, because of the
   qubit cap. One extra qubit is worth ~890 Toffolis, so the trade is already near its balance.
3. Our source rebuild is the 893-qubit version, ~3% behind the submitted 891, so source-level
   ideas must be ported or must exceed that gap.

## 3. invsqrt: phase by 1/√x

Relevant methods:
- Häner, Roetteler & Svore, *Optimizing quantum circuits for arithmetic* (arXiv:1805.12445):
  piecewise polynomial evaluation.
- Bhaskar et al. (2016): Newton iterations for inverse square root.
- Gidney: windowed arithmetic and table lookups.

The frontier already combines all of these: a lookup seed, a Newton step, measurement-conditioned
blocks and 62 qubits, minimal under relabelling. Its 41 rotation angles are individually tuned
(worst error 0.99932 of tolerance), and rotations are 19% of the cost.

Directions:
- A wider error margin for cheaper rotations: ×10 margin saves ~19 Toffoli-equivalents (1.2%).
- Fewer phase bits.
- A narrower architecture: 1 qubit ≈ 1.6%.

Ceiling: ~+1–3% overall.

## 4. dfthcsos: DF-THC block encoding with spectrum amplification

Literature:
- Lee et al., *Even more efficient quantum computations of chemistry through tensor
  hypercontraction* (arXiv:2011.03494): FeMoCo THC one-norm ≈ 306.
- Low, King, Berry, Babbush et al., *Fast quantum simulation of electronic structure by spectrum
  amplification* (arXiv:2502.15882; PRX 2025):
  - λ_eff = √(E_gap(2Λ − E_gap)), the formula the verifier uses.
  - FeMoCo-76 with DF-THC + BLISS + SA: ranks (R, B, C) = (15, 57, 19), Λ = 179.73, E_gap = 5.382,
    λ_eff = 43.65; stated lower bound λ_eff ≥ 21.69.
  - Accuracy filter: CCSD(T) correlation-energy error ≤ 0.3 mHa.
- The competition instead requires a Frobenius error ≤ 20 mHa on the decoded two-body tensor, fixes
  the BLISS shift, and floors E_gap at 3.914898.

Our measurements:
- At the frontier's frozen vectors, (w, W) is a KKT point. Error budget is worth 0.36 λ per mHa;
  splitting or dropping copies doesn't pay.
- 7-bit angles: error_2 455 mHa (the words are tuned to 8 bits).
- Cost model (verifier formulas, frontier loader cost per row):

| Scenario | λ_eff | Chemistry ratio | Overall |
|---|---|---|---|
| frontier (9, 143, 39), Λ 284.4 | 47.03 | 10.70 | 24.49 |
| same shape, Λ 240 | 43.17 | 11.65 | 25.20 |
| same shape, Λ 180 (literature level) | 37.3 | ~13.5 | ~26.4 |
| paper shape (15, 57, 19), Λ 180 | 37.31 | 13.33 | 26.36 (Q jumps to 662 at small B) |

- **Λ is the one chemistry lever with real size.** Whether Λ ≪ 284 is reachable under the 20 mHa
  Frobenius test is unknown; it needs new training: a GPU, quantization-aware rounding to 8 bits,
  then new loader certificates.
- The public generic loader builder costs ~7% more angle Toffolis than the frontier's certificates.
  So a new factorization must cut Λ by ~8% just to break even, unless we also build a certificate
  optimizer (GF(2) tensor-decomposition search over the 15 × 85 × 35 mixed tensors, ~205–219 terms per slot).

## 5. Recommended implementation order

| # | Work | Where | Expected overall | Risk |
|---|---|---|---|---|
| 1 | **Chemistry, lower Λ.** Retrain the factorization at the frontier shape with a Λ-weighted objective, starting from the C39 factorization (stays close to the 8-bit words); quantization-aware rounding to 8 bits; then our own loader-certificate search so loaders match the frontier | GPU (Kaggle) for training; certificate search on CPU here | +3–8% if Λ falls 15–35% | high: Λ under the Frobenius test may not move |
| 2 | **dqi, failure-budget trimming** of rare-weight logic, measured against the R penalty; port to the 891 circuit or exceed the 3% gap | CPU here | +1–3% | medium |
| 3 | **invsqrt** margin/rotation trade and narrower design search | CPU here | +1–3% | medium |

Combined, if 1–3 all succeed: about 26.5–28. Beating 29.5 needs, in addition, a new idea, most
plausibly a different dqi locator/scan formulation with far fewer multiplies, or a chemistry
factorization family with a much lower Λ at 20 mHa.

## Sources

- Khattar et al., Verifiable Quantum Advantage via Optimized DQI Circuits — https://arxiv.org/abs/2510.10967
- Gu, Kemeny, Yelin, Lukin, Discovering state-of-the-art quantum circuits with large language models — https://arxiv.org/abs/2610.22367
- Jordan et al., Optimization by decoded quantum interferometry — https://www.nature.com/articles/s41586-025-09527-5
- Quantum Circuit Design for Decoded Quantum Interferometry — https://arxiv.org/abs/2504.18334
- Optimized ancilla-free GF(2^m) multiplication/division — https://arxiv.org/abs/2511.20618
- Efficient exact quantum sampling from the Sun–Wootters distribution for OPI — https://arxiv.org/abs/2607.16541
- Low et al., Fast quantum simulation of electronic structure by spectrum amplification — https://arxiv.org/abs/2502.15882
- Lee et al., Even more efficient quantum computations of chemistry through tensor hypercontraction — https://arxiv.org/abs/2011.03494
- Häner, Roetteler, Svore, Optimizing quantum circuits for arithmetic — https://arxiv.org/abs/1805.12445
