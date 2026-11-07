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
