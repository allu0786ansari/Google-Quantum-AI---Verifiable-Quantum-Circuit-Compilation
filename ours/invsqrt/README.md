# Benchmark I (1/sqrt(x) phase): our own architecture, built layer by layer

Design: `docs/10_own_architecture_design.md`. Every layer is checked on all 3,133,452 inputs.

| Layer | File | Status | Cost |
|---|---|---|---|
| I-1 normalisation (binary s, 4 stages of 16/8/4/2-bit conditional shifts) | `layer1.py` | done: exact on all x; round trip clean | 84 Toffolis forward, 155 mean round trip, 37 qubits |
| I-2 seed (one 128-entry lookup of C0, C1; y0 = C0 - C1*t with 8 controlled subtractions; C1 cleared by del_lookup) | `seed_model.py`, `layer2.py` | done: bit-exact vs integer model on all x; round trip clean | +478 Toffolis forward (562 with I-1), 1,129 round trip, 70 qubits |
| I-3 Newton step | — | next | — |
| I-4 shift and phase | — | — | — |
| I-5 uncompute | — | — | — |

Seed accuracy (exact integer model): worst relative error 2^-13.61, i.e. 2^-26.63 after one exact
Newton step (budget 2^-24). Design correction: a pure lookup seed needs ~6,144 entries (~6,100
Toffolis), so a 96-entry linear seed is used instead (`model_l1_l2.py` compares the options).

Where I-2's cost goes: lookup 126 Toffolis, del_lookup ~17, each controlled 9-bit subtraction
~40 (x 8). The multiply is the next optimisation target.
