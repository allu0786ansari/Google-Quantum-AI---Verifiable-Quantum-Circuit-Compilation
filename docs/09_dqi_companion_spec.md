# DQI-A: specification of the unextract companion chain, and whether its excess can go

Source: 893-qubit source build (ahoirg's bundle). Files:
- `exp/new/a01/companion.py` (reverse `external_node_chain`);
- `exp/new/a14/decode_owner.py` (installs it for node groups 8..16);
- `playbook/kit/decoder/final_decoder/decode.py` (forward chain and step kernels).

Costs are mean Toffolis per shot from the bit-sliced simulator, which reproduces the verifier's
407,444.2.

## 1. What the node chain computes

Each extraction node group gi has:
- two GF(2²⁴) words N and L;
- 2·gi "heads" h_i (GF(2²⁴) values that persist through the whole decoder);
- classical constants q_i (nonzero) and c_i.

**Forward** (`external_node_chain`, reverse=False): with D = 1, for each head:

    (N, D) <- (q_i * D, N + h_i * D)      projective Möbius step, matrix M_i = [[0, q_i], [1, h_i]]
    L      <- L + h_i * D^2 (scaled)      the "corrected product"

- Each step costs ~200 Toffolis: two GF(2²⁴) products of 99 each (`mul24_xor99`). Squaring is linear.
- Then: pole = [D = 0], with a cswap.
- Normalise: N' = N / D_f, L' = L / D_f^2 (`external_normalize_paired`).
- D_f is measured away; the debt (-1)^(gauge · D_f) is repaid in the reverse.

With M = Π M_i = [[a0, b0], [c0, d0]] (entries are polynomials in the heads), the forward is the
Möbius map

    N' = (a0 N0 + b0) / (c0 N0 + d0),   and   det M = Π q_i, a classical constant.

**Reverse** (companion version, called 11 times per shot):

| Step | What it does | Toffolis per call |
|---|---|---|
| R <- copy N'; D <- 1 | | 0 |
| reverse chain over heads[1:] | applies adj(M_i) to (N', 1): gives (x, y) = (d0 N' + b0, c0 N' + a0), so N0 = x / y | ~200 per head (in unextract's untagged 51.6k) |
| companion_normalize | divides by y, i.e. normalises N, D, L, R together | 574 |
| companion_fused_chart | settles forward and reverse phase debts, recovers the denominator from its chart | 839 |
| measure old N, restore_direct, pole cleanup | | small |
| companion_return | scalar chain over all 2·gi heads (N, D only): erases the projective scale y | 1,980 (= 99 x ~20 heads) |

Total excess over the forward: 3,393 per call, 37.3k per shot. 21.8k of it is the return pass.

## 2. Why the return pass exists

After the division, the register D still holds the scale y = c0 N' + a0 = det M / D_f(N0). It
depends on the heads through the whole chain, so erasing it needs the chain value again. The
scalar chain costs 99 per head, without the L update.

## 3. Alternatives considered (all priced at 1 qubit ≈ 900 Toffolis in Q²T at Q = 892)

| Idea | Cost | Verdict |
|---|---|---|
| Measure y (MBU) instead of erasing it | debt (-1)^Tr(μ y) needs y again: scalar chain (99/head) plus a phase of an inverse; no saving | no |
| Pay y's debt together with the forward gauge debt | the gauge debt is on D_f = det M / y (payable while y is live, as the frontier does); y's own debt only exists after it is measured, so it cannot be prepaid | no |
| Keep D_f from the forward until the reverse, then start the reverse chain from (N', D_f) so no scale is left | 24 qubits x 11 groups held across synthesis, convert, scan and their reverses: far beyond the qubit budget | no |
| Accumulate M's entries (a0, c0) in registers during the reverse chain, so y = c0 R + a0 can be erased bilinearly (CZ phase, free) | +48 qubits at unextract, which is already the qubit peak: ≈ 43k Toffoli-equivalents against 21.8k saved | no |
| Merge the return pass into the first reverse pass | the return needs N0 = x / y, available only after the first pass and the division: inherently sequential | no |
| Leave N', L' unnormalised (projective) in the forward, so no division and no scale | changes the interface of synthesis/convert (they consume normalised N', L'); a full decoder redesign | open, but not a local change |

## 4. Conclusion

Under the decoder's width constraint (unextract is the qubit peak), the 99-per-head return pass is
the cheapest way to erase the projective scale. Every local alternative either needs the chain value
again or costs more in qubits than it saves in Toffolis. The only route that removes it is a
projective-interface redesign of extract → synthesis → convert (no normalisation, carry
(N, D) pairs). That is a new decoder architecture, not an edit, and it would have to be rebuilt
and verified end to end on the 893 source, and then ported to the 892 lineage, whose source is
unavailable.

Measured ceiling of this track if the projective redesign worked perfectly:
- the companion excess is 37.3k Toffolis (9.2% of T): dqi x1.10, about +3.2% overall;
- minus the cost of carrying projective pairs through synthesis/convert (more multiplies, plus
  24 extra qubits per carried denominator at the stages that are already near the peak).


## 5. Projective-interface redesign priced (7 Nov): closed

Extraction completes all 17 groups before synthesis (`decode`: `extract(b, S)`, then
`synthesize(b, S)` on the transformed bank). A denominator carried projectively must therefore stay
live through synthesis, convert, scan and their reverses. Those stages already run at ~881 of 892
qubits.

| Option | Extra qubits | Q^2 | Toffolis saved at most | Net dqi |
|---|---|---|---|---|
| carry D_f for all 11 companion groups | +264 | x1.68 worse | 37.3k (9.4%) | x0.60 |
| carry D_f for the costliest group only | +24 | x1.055 worse | <= 3.4k (0.9%) | x0.96 |
| fold D_f into scaled node values | 0 | - | - | not equivalent: downstream rational interpolation isn't homogeneous per node |

DQI-A is closed: no local or interface-level change to the companion chain lowers the dqi cost
under the decoder's qubit budget.
