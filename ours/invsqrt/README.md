# Benchmark I (1/sqrt(x) phase): our own architecture, built layer by layer

Design: `docs/10_own_architecture_design.md`. Every layer is checked on all 3,133,452 inputs.

| Layer | File | Status | Cost |
|---|---|---|---|
| I-1 normalisation (binary s, 4 stages of 16/8/4/2-bit conditional shifts) | `layer1.py` | done: exact on all x; round trip clean | 84 Toffolis forward, 155 mean round trip, 37 qubits |
| I-2 seed (one 128-entry lookup of C0, C1; y0 = C0 - C1*t with 8 controlled subtractions; C1 cleared by del_lookup) | `seed_model.py`, `layer2.py` | done: bit-exact vs integer model on all x; round trip clean | +478 Toffolis forward (562 with I-1), 1,129 round trip, 70 qubits |
| I-3 Newton step | `newton_model.py` | modelled, not built (see below) | row-adder form ~2,600-3,200 Toffolis: too expensive |
| I-4 shift and phase | — | — | — |
| I-5 uncompute | — | — | — |

Seed accuracy (exact integer model): worst relative error 2^-13.61, i.e. 2^-26.63 after one exact
Newton step (budget 2^-24). Design correction: a pure lookup seed needs ~6,144 entries (~6,100
Toffolis), so a 96-entry linear seed is used instead (`model_l1_l2.py` compares the options).

Where I-2's cost goes: lookup 126 Toffolis, del_lookup ~17, each controlled 9-bit subtraction
~40 (x 8). The multiply is the next optimisation target.


## Layer I-3 findings (Newton step)

Formulation: u*y0^2 = 1 - eps with |eps| <= 2^-12.6, so y1 = y0 + y0*eps/2. Only a window of
u*y0^2 below its known top bits is needed (taken mod 2^WE, read as the signed eps).

Exact integer model of row-by-row truncated arithmetic over all x (`newton_model.py`):

| FA | WE | F1 | worst error | est. Toffolis |
|---|---|---|---|---|
| 28 | 17 | 28 | 2^-22.8 (fails) | ~2,640 |
| 30 | 19 | 30 | 2^-24.8 (passes) | ~2,980 |
| 32 | 21 | 30 | 2^-25.4 | ~3,240 |

Row-by-row truncation costs ~2 bits against exact Newton (2^-26.6). At ~2.4 Toffolis per target bit
per controlled add, the Newton step alone would exceed the public circuit's whole Toffoli count
(1,281). Cheaper formulation to model next:
- v = u*y0 truncated at 2^-27 (~340 partial products);
- eps window = 1 - v*y0, columns 2^-27..2^-11 only (~310);
- correction y0*eps/2 with ~13-bit operands (~170);
- bias-compensated truncation constants (about 1 bit back);
- carry-save accumulation, where a full adder is 1 Toffoli, instead of a controlled ripple add per row;
- option: apply the phase of y0 and of the correction on separate registers with LP-fitted rotation
  angles, so y1 is never formed.


### Cheaper formulations, measured (exact models over all 3,133,452 x)

Target: worst relative error <= 2^-25.3, which leaves margin for the rotations under the 2^-24 tolerance.

| Design | File | Best passing config | Worst error | Partial products | Est. forward Toffolis | Lookup words |
|---|---|---|---|---|---|---|
| Newton, row adds (above) | `newton_model.py` | FA 30, WE 19, F1 30 | 2^-24.8 | - | ~2,980 | seed 27 bits |
| Newton, column-truncated, band-limited residual, phase split | `newton_band.py` | T1 28, T2 29, T3 30 | 2^-25.40 | 833 (792 ANDs) | ~1,560 | seed 27 bits |
| No Newton: piecewise cubic, K=7 (96 entries), truncated Horner, c0 refit last | `poly_model.py` | T 29, z bits 15/14/11 | 2^-25.29 | 449 | ~850 | 29+23+18+13 = 83 bits |
| Piecewise quadratic, K=8 (192 entries) | `poly_model.py` | T 29-30 | 2^-24.0 (fails) | 317-350 | - | 67 bits |

Error isolation (Newton band): the correction cut T3 dominates (T3 = 28 alone gives 2^-24.2); T1 and T2
cost ~0.3 bit each at 28-29.

**Budget check.** The public circuit is 62 qubits / 1,281 Toffolis (Q*C = 79.4k, ratio 21.25). For ours to
help, Q*C must be below that. Best case for the cubic route:
- C: normalise 155 + lookups ~150 + Horner ~850 + uncompute with MBU ~300 + shift and rotations ~480,
  about 1,900;
- Q: 22 (x) + 4 (s) + 29 (y) + coefficient words + scratch, about 90-140.

That gives Q*C ~ 170-270k, ratio ~8-12. This matches the design doc's 9-16 estimate. It sits below the public
21.25, and replacing it would cost about 13-20% of the overall score.


## Complete circuit, layers I-1 to I-5 (`full_circuit.py`, `layer_poly.py`)

| Layer | Implementation |
|---|---|
| I-1 | binary s (4 flags), CSWAP stages; flags cleared by measurement (`layer1.py`) |
| I-2/I-3 | two lookups (c3,c2 then c1,c0) on 7 address bits; Horner rows `target -= z_j * (A >> (nz-j))`; Y's top bit is a constant 1 |
| I-4 | 39-bit register: Y plus the constant bit plus 10 zeros, CSWAP shift by s (stages 1,2,4,8, data-carrying pairs only, 122 CSWAPs); 36 rotations, 1 T, 2 Clifford |
| I-5 | unshift, reverse rows, measurement-based lookup deletes, MBU denormalise |

Bit-exact with `poly_model` on all 3,133,452 inputs (0 failures, worst error 0.41 of tolerance).

| Row adder | Q | Mean Toffolis | C | Ratio | Check |
|---|---|---|---|---|---|
| controlled adder per row | 126 | 5,349 | 5,594 | 2.95 | all x, 1 shot each |
| AND row into scratch, then plain adder | 145 | 4,271 | 4,500 | **3.18** | **official verifier (quick): PASS** |

For comparison, our earlier own design `notebooks/src/invsqrt.py` v4 scores **7.50**, fully verified on 200.5M
shots. It is a piecewise quadratic on 512 segments (the parity of the shift is part of the address), with
radix-4 Booth rows, rows subtracted only from the low bits they can change plus a borrow qubit, and
row bits cleared by measurement. That design dominates the cubic route on both qubits and Toffolis,
so v4 stays our own Benchmark I circuit.
