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
