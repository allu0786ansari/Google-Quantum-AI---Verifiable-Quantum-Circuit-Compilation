# Research review, 14 Nov 2026: what it takes to pass 30

Leaderboard top: 29.9. Our submitted score: 24.94 (public portfolio of 7 Nov).

## 1. The public frontier moved to 27.29 in one week

Measured from public repositories and receipts (cloned read-only into `/home/claude/data/pub`):

| Task | 7 Nov public (ours) | 14 Nov public best | Ratio now | Source |
|---|---|---|---|---|
| invsqrt | 62 qubits, cost 97,746 | cost 78,766 (qubit count not confirmed here); a source-built **46-qubit** circuit: 80,163 | 26.37 | leotang09 `vqcc-t7421-submission`; Syed Asad Ali's Apache-2.0 generator (1,511.6 Toffolis, 33 rotations, 1 T) |
| dqi | 892 qubits, 3.748e11 | **889 qubits**, 3.602e11 (384,003 static Toffolis) | 67.85 | Syed Asad Ali's Apache-2.0 generator; an 888-qubit notebook by atakanaldemir exists |
| dfthcsos | 2.997e8 | 2.954e8 (coded-word tables) | 11.36 | maoxin123 `qcc-v25-coded-word-chemistry` |
| **Composite** | **24.94** | **27.29** | | Gloria-Deng reports 27.28997 on the public board, rank 404 of 3,023 on 14 Nov |

So about 400 teams already hold 27.29 or more. The public frontier gained 9.4% in seven days.

## 2. The gap

From 27.29 to 30 the product of the three ratios must grow by 1.33. Three ways to get there:

| Route | invsqrt | dqi | dfthcsos |
|---|---|---|---|
| Even | x1.10 | x1.10 | x1.10 |
| One task | x1.33 | - | - |
| Two tasks | x1.15 | x1.15 | - |

## 3. What the scoring models allow (checked, not assumed)

**dfthcsos**
- lambda_eff = sqrt(E_gap (2 Lambda - E_gap)). The public tables sit on the E_gap floor (3.914898) with
  Lambda about 283. The other branch (E_gap near 2 Lambda) is not reachable: per term, 2 Lambda - E_gap
  grows by 0.5 a^2 + a |d| with a = ||W||_1, and a is fixed by the two-body tensor.
- Cost model (`costs_prepare.py`): the shape (8, 143, 44) at 8 bits is on an edge of the model. The inner
  QROAM block size k = 4 is selected only because k = 5 fails the test `11 + 31*19 < 8*75` by equality.
  B <= 128 flips to k = 5 and 696 qubits (measured: (8,127,44) gives T*Q 1.37x worse). B = 143 is the
  largest B with nine 16-address blocks. The community has already mined this model.
- T = 6,662 fixed + about 8,090 loaders. The loaders are 55% of T: 15 circuits x 2 directions x
  (about 76 address-feature ANDs + about 197 products).
- The two-layer loader form (features, then products of two linear forms) is at its counting floor
  (197 products against 194) and the feature count is pinned by the 128-qubit cap (77 workspace qubits).
- **Open lever, measured here.** The vectors are localised (participation ratio: median 2.8 orbitals),
  so most angles act on little amplitude. Sensitivity of the two-body tensor to an angle at
  (vector a, tree node t): s = r(a,t) * ||sum_c W[a,c] L_c||, with r the amplitude entering the node.
  It is heavy-tailed (median is 0.4% of the maximum). A reverse-water-filling bit allocation frees:

  | Extra quantisation-noise power | Table bits that still matter |
  |---|---|
  | +1% | 81% |
  | +10% | 68% |
  | +50% | 57% |

  The public tables free only bit 0 on some columns with Reed-Muller codes (at most about 9% of the bits);
  products per circuit went from 205-219 in the 7 Nov tables to 197-198 now. If products keep scaling with the bits that matter,
  68% gives about 150 products per circuit: loaders -17%, T -9%, ratio +10%.
  This needs a masked GF(2) completion compiler and a GPU refit of the weights; neither exists yet.

**invsqrt** (46 qubits x 1,742.7)
- Rotations: 33 x about 7.0 Toffoli-equivalents = 231 (13%). No lever: there is no Hadamard gate in
  kickmix, so no phase-gradient state; a larger error margin saves under 0.3%.
- Toffolis (from the public generator's own accounting at 48 qubits): divider 587, product 240,
  loaders 150, coefficient 123, normalisation about 122. Each qubit is 2.2% of the score.
- The public line went 97.7k, 89.8k, 87.8k, 80.2k, 78.8k in about three weeks.

**dqi** (889 qubits, 384k Toffolis)
- Qubits: 49 above the 840-qubit floor; removing all of them is worth 1.12x.
- Scan is at its floor (11 Toffolis per position, 70 positions removed by mask projection).
- Classical computation on the E registers is free, but the only use found is the mask projection
  that is already in place.
- The public line went 897, 893, 892, 891, 890, 889, 888 qubits in three weeks.

## 4. Recommendation

1. **Submit the 27.29 public portfolio now** (three components above). This is +9.4% for no research risk.
2. **Track public releases weekly** with our harness (`tools/harness.py`); the public frontier has been
   moving several percent a week and both leading generators are now Apache-2.0 source.
3. **One differentiating contribution: sensitivity-masked chemistry loaders** (section 3). Expected
   +9-17% on dfthcsos, i.e. +3-5.5% on the composite, independent of what others publish.
   Steps: (a) masked completion compiler on CPU, tested on the public table; (b) GPU refit of (W, H2_I)
   at the masked words on Kaggle; (c) official verifier.
4. **Micro-optimisation on the published invsqrt and dqi generators** only if the builders may be run
   here; expected +1-3% on the composite.

Expected composite: 27.29 now; about 28.5-29 with item 3; passing 30 additionally needs either item 4
or two to three more weeks of public progress. Reaching the top of the board is a different problem:
the leaders will have moved too.

## 5. Limits of this review

- Kaggle is not reachable from this machine (no credentials, proxy refuses), so notebooks and the
  forum were read only through search snippets and through GitHub mirrors.
- Third-party build scripts were not run (blocked by this session's permission policy), so the 46-qubit
  and 889-qubit circuits were not rebuilt or re-verified here; their numbers are the authors' receipts.
- The loader estimate assumes product count stays proportional to the bits that matter. That held for the
  one data point available (bit-0 codes); it is not yet demonstrated for unstructured masks.

## 6. Update, 15 Nov: the public portfolio, fully verified on Kaggle

`notebooks/vqcc_public_best.ipynb` downloads the public sources, verifies every candidate with the official
verifiers and keeps the cheapest verified circuit per task. Kaggle run of 15 Nov (all full verifications passed):

| Task | Score | Ratio | Detail |
|---|---|---|---|
| invsqrt | 78,719.5 | 26.382 | 46 qubits, 1,479.75 Toffolis, 33 rotations |
| dqi | 3.6021e11 | 67.848 | 889 qubits, 382,600 Toffolis, 0 failures in 102,400 shots |
| dfthcsos | 2.95128e8 | 11.370 | lambda_eff 46.897, 14,738 Toffolis, 427 qubits, 8-bit words |
| **Overall (local estimate)** | | **27.30** | |

All three chosen circuits came from the same public archive. The source-built 46-qubit invsqrt (80,171 by the
quick estimate) and the source-built 889-qubit dqi were also found and ranked lower.

## 7. Update, 15 Nov: what the chemistry coding lever is worth (measured)

Experiment `ours/chem/coding_noise.py` on the public (8, 143, 44) table whose bit 0 is already coded
on all 75 columns. Cost of forcing one more bit-plane into RM(1,4) per 16-address block, in units of
the plain 8-bit rounding variance D0:

| Bit-plane | Plain nearest codeword | Sensitivity-weighted codeword | Weighted, vectors sorted by importance |
|---|---|---|---|
| bit 1 | 11.4 | **1.70** | 4.0 |
| bit 2 | 45.1 | 6.9 | 18.1 |

- Weighted decoding cuts the added error variance 6.7x. Sorting vectors by importance makes it worse:
  important vectors should be spread across blocks, a few per block, where the 5 free coefficients of an
  affine code word can match them exactly.
- From the public training logs, lambda_eff moves about 0.34% per mHa of error budget
  (47.255 at 19.49 mHa, 47.214 at 19.75 mHa).
- Coding bit 1 on every column removes about 10% of the products, about 4% of the Toffolis. With the
  weighted decoder the price is roughly 1% of lambda_eff after a refit. Net: about +3% on dfthcsos,
  about +1% on the composite. Bit 2 is too expensive with RM(1,4); it needs a weaker code.

So this lever is real but small: it does not carry the composite past 28 on its own.
The notebook now also discovers the top-scoring public notebooks by itself (`AUTO_DISCOVER`), so a later
run picks up newer public circuits without editing the source list.

## 8. Update, 15 Nov: no downloaded circuits

Decision: we do not submit other teams' finished circuits. `notebooks/vqcc_public_best.ipynb` is removed.
`notebooks/vqcc_from_source.ipynb` builds all three tasks from published generator source at a pinned
commit (46-qubit invsqrt and 889-qubit dqi by Syed Asad Ali, Apache-2.0; chemistry recipe and loader
compiler by Francisco Angulo de Lafuente, MIT or Apache-2.0), with credits in the notebook, then runs the
official verifiers and packages the result. Expected from the authors' receipts: about 27.0.

Tested here: setup, the pinned download, the archive hashes, unpacking, verification and packaging.
Not tested here: the three generator runs themselves (this session's permissions block running them).
The generators are still other people's work; what is ours is listed in the notebook's last section.

## 9. Update, 15 Nov: the bit-1 coding addition was tested and fails without retraining

`ours/chem/code_bit1.py` applies the weighted bit-1 coding to the public table; errors are from our exact
evaluator (`tools/experiments/fastfit.py`, which reproduces the table's own 19.978 mHa and lambda_eff 47.110).

| Columns coded | Two-body error before refit | One-body error | After a weight-only refit |
|---|---|---|---|
| none (public table) | 19.98 mHa | 11.2 mHa | - |
| last 3 batches (15 columns) | 42.6 mHa | 159 mHa | 39.1 mHa, infeasible (limit 20) |
| batches 8-14 (35 columns) | 38.5 mHa (one-body rows less protected) | 156 mHa | not run |
| all 75 columns | 203 mHa | 1,120 mHa | not run |

Why: single words are far more sensitive than section 7 assumed. Changing the most sensitive word by two
steps adds 65 mHa on its own (checked on five entries; the sensitivity model is right to within a factor
1.4). The public table is not a rounded table: it was trained through the rounding, so its words are
near-exact where it matters, and "units of rounding variance" in section 7 meant about 180 mHa, not 6.
The estimate in section 7 (+3% on dfthcsos for about 1% of lambda_eff) is withdrawn.

What would work is what the public pipelines do for bit 0: retrain the vectors in stages with the code as a
constraint. That is a GPU training job with its own trainer, not a post-processing step.
