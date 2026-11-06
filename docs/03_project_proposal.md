# Project proposal: compressing three verifiable quantum circuits

*Kaggle competition: Google Quantum AI – Verifiable Quantum Circuit Compilation.*
*Draft, 6 Nov 2026.*

## 1. Summary

We will build smaller versions of three quantum circuits that Google uses as benchmarks,
and prove they work using Google's own verifiers. Each circuit is scored by how much
cheaper it is than Google's reference version. The overall score is the geometric mean of
the three improvement ratios.

Our first goal is a verified submission that beats all three baselines before the
**20 November 2026 milestone prize**. Our second goal is to keep improving through the
final deadline on 20 January 2027.

Everything needed is already in place: the verifiers, the simulator (kickmix), the sample
submission and the paper. All three baseline scores reproduce exactly on our machine.

## 2. Problem statement

Future error-corrected quantum computers will be limited by two costs:

- **Toffoli gates**, which set the running time, because each one consumes an expensive
  "magic state".
- **Logical qubits**, which set the size of the machine.

Many important quantum algorithms spend most of their cost in subroutines that are really
classical arithmetic (adders, multipliers, table lookups) run on a superposition of inputs.
These subroutines can be checked on a classical computer one input at a time, so any
improvement can be verified automatically. The competition asks for cheaper versions of
three such subroutines.

### 2.1 Task 1: invsqrt (fusion energy)

- **What the circuit does:** given a 22-bit integer x (a squared distance between two
  electrons), apply the phase exp(2πi/√x), then clean up all scratch qubits.
- **Accuracy:** relative error at most 2⁻²⁴ (about 7 significant digits), for every x from
  1 to 3,133,452. Every input is tested 64 times. A single failure rejects the circuit.
- **Cost:** qubits × (Toffolis + T gates / 4). Arbitrary-angle phase rotations are priced
  in T gates.
- **Baseline:** 260 qubits, 7,735 Toffolis, 37 rotations. Cost ≈ 2.07 million.
  It uses a cubic first guess followed by one Newton step, then reverses everything.

### 2.2 Task 2: dqi (combinatorial optimization)

- **What the circuit does:** a reversible Reed–Solomon decoder over the field GF(2¹²).
  Given 70 syndrome values (840 qubits) and a classical record of 4,095 values, it must
  find the error pattern implicitly, apply a ±1 sign that depends on it, and restore the
  syndromes and all workspace.
- **Accuracy:** tested on 102,400 random error patterns (up to 35 errors). Occasional
  failures are allowed but are charged through a repetition factor R(ε).
- **Cost:** R(ε) × qubits² × Toffolis. A perfect decoder still pays R ≈ 1.19.
- **Baseline:** 1,893 qubits, 5.73 million Toffolis. Cost ≈ 2.44×10¹³.
  It uses a synchronized extended Euclidean algorithm, then a per-position Chien search
  with one field inversion at each of the 4,095 positions.

### 2.3 Task 3: dfthcsos (quantum chemistry)

- **What we submit:** a compressed (tensor-factorized) description of the FeMoCo molecule's
  Hamiltonian, plus the lookup circuits that load its rotation angles.
- **Accuracy:** the Hamiltonian rebuilt from the angles our circuits actually load must be
  within 20 mHa of the reference (one-body and two-body residuals). The normalization Λ
  must stay above 49.0 and the energy gap above 3.914898 Ha.
- **Cost:** λ_eff × Toffolis × qubits for the whole block encoding. Only the angle-loading
  circuits are simulated; state preparation and Givens rotations are priced by formula.
- **Baseline:** λ_eff 49.30, 93,237 Toffolis, 730 qubits. Cost ≈ 3.36×10⁹.
  Its two-body residual is already 18.64 mHa, close to the 20 mHa limit.

### 2.4 How the final score works

Overall score = (r₁ · r₂ · r₃)^(1/3), where rᵢ = baseline cost / our cost for task i.

- Each task matters equally in multiplicative terms. Halving any one task's cost raises the
  overall score by about 26%.
- A task left at its baseline scores exactly 1.0, so we can submit early and improve one
  task at a time.
- Separate $10,000 prizes go to the best score on each individual task.

## 3. Why this is hard

**Correctness is all or nothing (invsqrt).** Every one of 3.1 million inputs must pass, 64
times each. The baseline has only about one bit of accuracy to spare, so any
approximation must be checked against every input in the exact fixed-point arithmetic the
circuit uses, not in floating point.

**Everything must be reversible.** A quantum circuit can't overwrite or erase values.
Every intermediate result has to be cleaned up, which usually doubles the cost (compute,
use, then uncompute). Measurement-based tricks can make some clean-up nearly free, but
only for certain kinds of values, and a mistake shows up as a random sign error.

**Qubits and Toffolis trade against each other.** In invsqrt a 10% saving in qubits is
worth as much as a 10% saving in Toffolis. In dqi qubits are squared, so a 10% qubit saving
is worth about 20% in Toffolis. Many classic optimizations save one by spending the other.

**Data-dependent control flow (dqi).** The Euclidean algorithm normally takes a different
path for each input. In a quantum circuit all paths run at once in superposition, so the
circuit must follow a fixed, padded schedule. The baseline does a full multiplication sweep
on every one of 222 iterations, even when no multiplication is needed. Fixing this means
redesigning the algorithm's schedule, not just tuning it.

**Huge circuits (dqi).** The reference decoder has 86 million operations and an 800 MB file.
A full check takes about 25 minutes on our machine, so iteration is slow.

**Coupled classical and quantum design (dfthcsos).** The cheapest circuits need fewer
angles, fewer bits per angle and fewer vectors, but each of those makes the factorization
less accurate, and the baseline is already near the accuracy limit. We tested this: just
rounding the baseline's 12-bit angles to 11 bits raises the one-body error from 3 to
221 mHa. Gains here need re-optimizing the factorization and the circuits together.

**Fixed rules and fresh tests.** Each submission's tests are seeded from a hash of the
circuit plus a secret salt, so we can't tune against known tests. Local checks must be
thorough enough that the leaderboard run passes the first time.

## 4. Solution approach

Our general method is a closed loop:

1. Find where each baseline spends its cost (per-stage Toffoli and qubit profiles).
2. Redesign the most expensive stage.
3. Check it against an exact classical model, then against the official verifier.
4. Keep only verified improvements.

### 4.1 invsqrt: replace the Newton step with a direct polynomial

The baseline computes a rough cubic guess and then refines it with a Newton step. The
Newton step and its uncompute are about 60% of the cost.

Our plan:

- **Piecewise polynomial instead of Newton.** Split the normalized input range into many
  short segments and fit a low-degree polynomial on each. Short segments make a cubic
  (or even a quadratic) accurate enough on its own. Our estimate: 32 segments per octave
  for a cubic, or 128 for a quadratic, reaches the 2⁻²⁴ target with room for rounding.
- **Cheap coefficient lookup.** kickmix's lookup table loads a segment's coefficients for
  about one Toffoli per table entry, and clears them again for very little. We measured
  about 25 Toffolis to clear a 256-entry table, against about 254 to load it.
- **Fit to the real arithmetic.** Because we can check all 3.1 million inputs in seconds
  with a classical model, we will tune coefficients and register widths against the exact
  truncated fixed-point arithmetic the circuit uses.
- **Cheaper denormalization.** The baseline shifts the result back with 300 controlled
  swaps each way. The shift amount can be converted to binary for free (with CNOTs only),
  and a 4-stage binary shifter should need about half as many swaps.
- **Fewer qubits.** Without Newton's wide intermediate registers, we expect to need roughly
  140–170 qubits instead of 260.

Target: at least a 4× improvement on this task. Our rough estimate is 4–7×.

### 4.2 dqi: cheaper field arithmetic first, then a leaner schedule

GF(2¹²) multiplications are about 59% of all Toffolis. That makes field arithmetic the
first target:

- **Cheaper inversion in the Chien/Forney stage.** The baseline does one inversion (5
  multiplications) at each of 4,095 positions, about 1.29 million Toffolis. Rewriting the
  field as GF((2⁶)²) is a free change of basis, and in that form an inversion needs a few
  small GF(2⁶) operations. We estimate about 40% of the current inversion cost, which
  would save roughly 0.7 million Toffolis (about 13% of the task).
- **A cheaper multiplier.** The built-in multiplier uses 63 Toffolis. Known constructions
  need about 54, and better ones may exist. This saves about 5% more.
- **A leaner Euclid schedule (stretch goal).** Cut the multiply sweeps on iterations that
  don't need them, or switch to an algorithm with a fixed schedule. This is the largest
  potential gain and the riskiest.
- **Fewer qubits (stretch goal).** The 840 syndrome qubits are fixed, but the Euclid
  workspace (about 1,050 qubits) has some slack.
- **Use the free failure budget.** About 10⁻⁴ failures barely change R, so rare cases
  (such as unusually long quotients) can be dropped where that saves cost.

Target for the milestone: 1.2–1.3×. Stretch target after the milestone: 2× or more.

### 4.3 dfthcsos: re-optimize the factorization at lower cost

The angle-loading circuits are 85% of the Toffolis, and each one costs about one Toffoli
per address. So the cost is driven by how many vectors there are (R·B), how many bits each
angle has, and how many angles fit per circuit.

Our plan:

- **Quantization-aware training.** Use the starter kit's trainer, but optimize the factorization
  for the precision it will actually be stored at. Then pick the angle words directly,
  rather than rounding.
- **Search over shape and precision.** Try fewer two-body vectors (smaller R or B) and fewer
  bits per angle (brot = 10 or 11), staying under 20 mHa. Our cost model shows halving B
  alone would cut T × Q by about 45%.
- **Lower Λ.** λ_eff grows roughly with √Λ. The gap is already near its floor, so gains in
  λ_eff must come from lowering Λ.
- **Circuit tricks.** Workspace inside the angle circuits is free up to the 128-qubit cap,
  and segments of identical rows cost no Toffolis. Reordering the vectors (which we're free
  to do) can create such segments.

Target for the milestone: 1.2× or more. This task depends on the starter kit, which we
don't have yet.

### 4.4 Expected outcome

| Scenario | invsqrt | dqi | dfthcsos | Overall score |
|---|---|---|---|---|
| Baseline | 1.0 | 1.0 | 1.0 | 1.00 |
| Conservative | 3.0 | 1.15 | 1.0 | 1.51 |
| Milestone target | 4.0 | 1.3 | 1.3 | 1.89 |
| Stretch | 6.0 | 1.6 | 2.0 | 2.68 |

These are projections, not results. Each will be replaced with verified numbers as the
work lands.

## 5. Deliverables

1. Improved circuits for all three tasks, each passing the official verifier locally.
2. A self-contained Kaggle notebook that builds every circuit from source, runs checks and
   writes `submission.zip`. It must be public to qualify for the milestone prize.
3. Documentation in this repository: analysis, verified results and how to reproduce them.

## 6. Risks

| Risk | Effect | Mitigation |
|---|---|---|
| invsqrt fails on a few inputs | Submission rejected | Exact fixed-point model checked on all 3.1M inputs before every build, then the full official verify |
| dqi changes are slow to verify | Fewer iterations | Small-shot checks first; full 102,400-shot runs only for candidates |
| dfthcsos starter kit is delayed | No chemistry gain by the milestone | Submit that task at baseline (ratio 1.0) and keep working on it |
| Leaderboard kickmix differs from ours | Mismatched scores | Already checked: the Kaggle and PyPI builds give identical results |
| Notebook too slow on Kaggle | Can't regenerate submission | Build is fast (seconds to minutes); full verification is optional and switchable |
