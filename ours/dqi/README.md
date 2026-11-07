# Benchmark II: our own decoder

| Version | Design | Q | Mean Toffolis | Failures / 102,400 | R | Score | Ratio | Check |
|---|---|---|---|---|---|---|---|---|
| v2 (earlier) | Google reference decoder + our tower / GF(4)-CRT field arithmetic | 1,866 | 3,840,178 | 0 | 1.19 | 1.59e13 | 1.53 | full PASS |
| v3 | **new key-equation layer**: normal-schedule register-shared Euclid (35 fixed rounds, no pointers, freed words relabelled, sticky freeze for k < 35) + Chien/Forney scan | 1,843 | 919,622 | 630 (0.62%) | 1.94 | 6.05e12 | **4.04** | full PASS |
| v4 | v3 + Karatsuba-paired multiply-accumulate (3 products per 2 coefficients), lower scratch peak, measurement-based zero-test cleanup | 1,801 | 791,004 | 619 (0.60%) | 1.93 | 4.94e12 | **4.94** | full PASS |

Files: `eea_model.py` is the word-level model; the circuits mirror it word for word. `decoder_v3.py` and
`decoder_v4.py` are the builders. `tower.py` holds the field arithmetic (multiply 42, divide 119,
inverse 77 Toffolis).

Failures are the abnormal Euclid steps: a remainder's leading coefficient vanishes before the collapse,
about 1/4096 per round. Handling them would bring R from 1.93 to 1.19, a 1.6x gain on the dqi score.
