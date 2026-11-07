# Our own architectures: end-to-end results

Every circuit here is built by our code from the competition paper's problem definitions and checked
with the official verifiers. No circuit or factorization from the public notebooks is used. Scores are
improvement ratios over the competition baselines. Higher is better, and 1.0 means the baseline.

## Benchmark I: phase by 1/sqrt(x) (score Q x C)

| Design | Layers | Q | Toffolis | C | Ratio | Check |
|---|---|---|---|---|---|---|
| Layered cubic (this round) | binary s (4 flags), 96-entry cubic, AND-row Horner, CSWAP shift, 36 rotations, MBU uncompute | 145 | 4,271 | 4,500 | 3.18 | official, quick (PASS) |
| **v4 (own, earlier)** | 5-stage normalisation, 512-entry quadratic, radix-4 Booth rows on low-bit spans with a borrow qubit, measurement-cleared row bits | **128** | **1,914** | **2,162** | **7.50** | **official, full (200.5M shots)** |

The quadratic on 512 segments with Booth rows dominates the cubic route. Its segments are non-uniform
(the parity of the shift is part of the address), which places the fine segments where the curvature
is highest, and a lookup costs about one Toffoli per entry whatever the word width.

## Benchmark II: Reed-Solomon decoding phase oracle (score R Q^2 T)

| Design | Q | Toffolis | Failures | R | Ratio | Check |
|---|---|---|---|---|---|---|
| v2 (own, earlier): reference decoder + our tower / GF(4)-CRT arithmetic | 1,866 | 3,840,178 | 0 | 1.19 | 1.53 | full |
| v3: **normal-schedule register-shared Euclid** + Chien/Forney | 1,843 | 919,622 | 0.62% | 1.94 | 4.04 | full |
| **v4**: v3 + Karatsuba-paired MACs, lower scratch peak, MBU zero test | **1,801** | **791,004** | 0.60% | 1.93 | **4.94** | **full (102,400 shots)** |

Key-equation layer (`ours/dqi/eea_model.py`, `decoder_v4.py`):
- Euclid on (x^70, S(x)) almost always has degree-1 quotients until the remainder collapses at round
  k, the number of errors. So 35 fixed rounds of two quotient coefficients each cover 99.4% of shots.
- In that schedule every coefficient position is known in advance. There is no pointer arithmetic, no
  unary iteration and no array shifting: the two zeroed top words of the dividend's remainder are
  relabelled as its cofactor's new top words.
- Quotients are uncomputed from the cofactor's top two words. Rounds after the collapse are frozen by
  a sticky flag that zeroes the quotients. One controlled swap at the end, conditioned on the parity
  of the number of frozen rounds, puts sigma and Omega in fixed words for the scan.
- Remaining failures: a remainder's leading coefficient vanishes before the collapse (about 1/4096 per
  round). Handling them exactly would cut R from 1.93 to 1.19, about 1.6x on this benchmark.

## Benchmark III: DF-THC block encoding (score lambda_eff T Q)

| Design | Shape | brot | lambda_eff | T | Q | Ratio | Check |
|---|---|---|---|---|---|---|---|
| Our pipeline on the provided factorization | (65, 60, 60) | 12 | 49.30 | 93,237 | 730 | 1.00 | full |
| Own factorization, trained from our DF initialisation | (8, 143, 44) | - | 88 (not feasible) | - | - | - | did not reach the 20 mHa gate |

Own compact factorization:
- Initialisation (`ours/chem/df_init.py`): weighted pair-basis eigendecomposition of the
  symmetry-shifted two-body tensor; 352 components split into 8 groups; greedy shared vectors and
  least-squares W per group. Starting error 4.5 Ha.
- Training: the starter kit's trainer on the paper's objective (eq. 48), 2 CPU cores, about 130 ms
  per epoch. After 12,000 epochs: error 64 mHa, Lambda 784 Ha, lambda_eff 88 Ha. Passing needs 20 mHa
  and Lambda near 285 Ha (lambda_eff about 47). The trend points to 10^5 or more epochs, which is not
  feasible on this machine. The run hit the 2-hour background limit before saving.
- The cost model shows what such a factorization is worth with plain lookups, at brot 8, J 13, K 6:
  ratio 6.9 at lambda_eff 46.9 (`tools/harness.chem_model`). The factorization is the binding step,
  and it needs GPU training.

## Combined (geometric mean of the three ratios)

| Portfolio | invsqrt | dqi | dfthcsos | Overall |
|---|---|---|---|---|
| Our own architectures | 7.50 | 4.94 | 1.00 | **3.33** |
| Public 24.94 portfolio (submitted) | 21.25 | 65.22 | 11.20 | **24.94** |
| Best of both, per task | 21.25 | 65.22 | 11.20 | 24.94 |

What the experiment shows:
1. On every benchmark, the public circuits are 3x to 13x better than anything we built from the paper
   in this session. The gaps are structural:
   - invsqrt: the 62-qubit circuit computes in place;
   - dqi: the 892-qubit decoder uses only 52 qubits beyond the 840 syndrome qubits, against our 961;
   - dfthcsos: a GPU-trained compact factorization, with loaders compressed about 5x.
2. Our new key-equation layer for dqi is a real architectural gain over the paper's reference: 7.2x
   fewer Toffolis at slightly fewer qubits, and 3.2x the score of our earlier v2.
3. The ranked next steps for our own pipeline:
   - dqi: exact handling of abnormal Euclid steps (R 1.93 to 1.19, about 1.6x);
   - dfthcsos: a compact factorization trained on a GPU (1.0 to about 5-7);
   - invsqrt: an in-place design that drops the separate coefficient registers (Q 128 to about 80).
