# %% [markdown]
# # Verifiable Quantum Circuit Compilation: submission notebook
#
# This notebook builds all three circuits from source, checks them, and writes
# `submission.zip` in the layout the competition expects:
#
# ```
# invsqrt.kmb
# dqi.kmb
# dfthcsos/angle_k_0.kmb ... angle_k_<K-1>.kmb, submission.npz
# ```
#
# | Task | Version in this notebook |
# |---|---|
# | invsqrt | **Improved**: piecewise quadratic replaces the Newton step (section 2) |
# | dqi | Baseline reference decoder (improvement in progress) |
# | dfthcsos | Baseline sample submission (improvement in progress) |
#
# Every task falls back to its baseline circuit if a check fails, so the zip is
# always valid. Settings are in the next cell.

# %%
# ---- Settings ---------------------------------------------------------------
VERIFY = 'quick'      # 'none' | 'quick' (minutes) | 'full' (official verifiers, ~30-60 min)
USE_IMPROVED = {'invsqrt': True, 'dqi': False, 'dfthcsos': False}

# %% [markdown]
# ## 1. Setup
#
# Finds the competition data (on Kaggle it is under `/kaggle/input`; locally set
# the `VQCC_DATA` environment variable), imports kickmix and the official
# verifiers, and prepares the output folder.

# %%
import glob
import io
import json
import math
import os
import shutil
import sys
import time
import zipfile
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path

import numpy as np


def find_data_root() -> Path:
    candidates = []
    if os.environ.get('VQCC_DATA'):
        candidates.append(Path(os.environ['VQCC_DATA']))
    candidates += [Path(p).parents[2] for p in glob.glob(
        '/kaggle/input/**/verifiers/invsqrt/verify_submission_invsqrt.py', recursive=True)]
    for c in candidates:
        if (c / 'verifiers' / 'invsqrt' / 'verify_submission_invsqrt.py').exists():
            return c
    raise FileNotFoundError('Competition data not found. Attach it, or set VQCC_DATA.')


DATA = find_data_root()
VERIFIERS = DATA / 'verifiers'
SAMPLE = DATA / 'sample_submission'
try:
    import kickmix as km
except ImportError:
    sys.path.insert(0, str(DATA / 'kickmix'))
    import kickmix as km
sys.path.insert(0, str(VERIFIERS))

OUT = Path('/kaggle/working') if Path('/kaggle/working').exists() else Path('output').resolve()
SUB = OUT / 'submission'
if SUB.exists():
    shutil.rmtree(SUB)
(SUB / 'dfthcsos').mkdir(parents=True)
REPORT = {}
print('data:', DATA)
print('kickmix:', km.__file__)
print('output:', SUB)

# %% [markdown]
# ## 2. invsqrt: phasing by 1/sqrt(x)
#
# **Contract.** A 22-qubit register `x` holds 1 <= x <= 3,133,452. The circuit
# must apply the phase exp(2 pi i / sqrt(x)) with relative error at most
# 2^-24 on every input, and return all scratch qubits to 0.
# Cost = qubits x (Toffolis + T/4); rotations are priced in T gates.
#
# **Baseline** (260 qubits, 7,735 Toffolis): cubic first guess, one Newton step,
# then a full reverse pass.
#
# **Our design.**
#
# 1. *Normalize* with a 5-stage binary search on leading zeros (shifts of 16, 8,
#    4, 2, 1). Then x_norm = x * 2^t has its top bit set, M = x_norm / 2^21 is
#    in [1, 2), and 1/sqrt(x) = p * 2^(sh - 10) with sh = t >> 1 and
#    p = M^-1/2 (times 1/sqrt(2) when t is even), so p is in (1/2, 1].
#    The stage flags *are* the binary digits of t, so no conversion is needed.
#    Clearing each flag afterwards uses measurement plus a phase fix that only
#    runs on half the shots.
# 2. *Look up* the coefficients of a quadratic for M's segment
#    (address = parity of t and the top k bits of M). kickmix's lookup costs
#    about one Toffoli per entry and clearing it costs about a tenth of that.
# 3. *Evaluate* y = c0 - (c1 - c2*d)*d with truncated fixed-point multipliers,
#    where d is the rest of M's bits.
#    - The multipliers use radix-4 Booth recoding: each row handles two bits of d
#      with a digit in {-2, ..., 2}, so there are half as many carry-propagating
#      subtractions and half as much truncation noise. A negative digit costs
#      nothing extra: complement the target with CNOTs, subtract, complement back.
#    - Each row is subtracted only from the low bits that can change, plus a
#      borrow qubit, followed by one controlled decrement of the high bits.
#    - Row bits are cleared by X-basis measurement with conditioned CZ fixes, so
#      clearing them costs no Toffolis.
# 4. *Shift* the phased part of y left by sh with kickmix's barrel shifter, apply
#    one Z_POW per bit, then undo steps 4 to 1.
#
# Coefficients are fitted per segment, then refined by least squares against the
# exact truncated arithmetic, so the model below reproduces the circuit's
# output bit for bit on all 3.1M inputs.

# %%
INV_WIDTH, INV_MAX_X = 22, 3 * (2**10 - 2)**2


@dataclass(frozen=True)
class InvDesign:
    k: int        # mantissa bits in the lookup address (table has 2^(k+1) entries)
    frac: tuple   # fractional bits of registers c0 (= output y), c1, ..., c_deg
    dbits: tuple  # delta bits used by the Horner step that produces c0, c1, ...
    pf: int       # bits of p below the binary point that get phased

    @property
    def degree(self):
        return len(self.frac) - 1

    @property
    def widths(self):
        return [f + (1 if j == 0 else 0) for j, f in enumerate(self.frac)]


def inv_rows(d, j):
    """(x_norm bit index, weight exponent m) for each multiplier row of step j."""
    top = INV_WIDTH - 2 - d.k
    return [(i, 21 - i) for i in range(top, max(top - d.dbits[j], -1), -1)]


def inv_inputs():
    xs = np.arange(1, INV_MAX_X + 1, dtype=np.int64)
    lz = INV_WIDTH - 1 - np.floor(np.log2(xs)).astype(np.int64)
    xn = xs << lz
    assert ((xn >> (INV_WIDTH - 1)) == 1).all()
    e = lz & 1
    return xs, lz, xn, e


def inv_booth_digits(d, j):
    """Radix-4 Booth digits of step j: (lo, mid, hi, m) with x_norm bit indices
    (None = constant 0); digit = -2*hi + mid + lo, weight 2^-m."""
    used = sorted(i for i, _ in inv_rows(d, j))
    n = len(used)
    bit = lambda t: used[t] if 0 <= t < n else None
    m0 = 21 - used[0]
    out = []
    for jj in range((n + 2) // 2):
        lo, mid, hi = bit(2 * jj - 1), bit(2 * jj), bit(2 * jj + 1)
        if (lo, mid, hi) != (None, None, None):
            out.append((lo, mid, hi, m0 - 2 * jj))
    return out


def _inv_shift(r, s):
    return (r >> s) if s >= 0 else (r << -s)


def inv_steps(d, table, xn, e):
    """Exact fixed-point Horner evaluation. Yields (j, c_j, result_j) per step."""
    addr = (e << d.k) | ((xn >> (INV_WIDTH - 1 - d.k)) & ((1 << d.k) - 1))
    r = table[addr, d.degree].copy()
    for j in range(d.degree - 1, -1, -1):
        acc = table[addr, j].copy()
        fin, fout = d.frac[j + 1], d.frac[j]
        for lo, mid, hi, m in inv_booth_digits(d, j):
            g = lambda i: ((xn >> i) & 1) if i is not None else 0
            dg = -2 * g(hi) + g(mid) + g(lo)
            s = m + fin - fout
            row = np.where(np.abs(dg) == 2, _inv_shift(r, s - 1), np.where(np.abs(dg) == 1, _inv_shift(r, s), 0))
            acc -= np.sign(dg) * row
        yield j, table[addr, j], acc
        r = acc


def inv_evaluate(d, table, xn, e):
    """Returns y (the final Horner result) in units of 2^-frac[0]."""
    for j, c, acc in inv_steps(d, table, xn, e):
        pass
    return acc


def inv_phase_turns(d, y, lz):
    """The phase the circuit applies, in turns (bits below 2^-(pf+10) are not phased)."""
    return (((y << (lz >> 1)) >> (d.frac[0] - d.pf))) / 2.0**(d.pf + 10)


def inv_fit(d, refine_iters=3):
    """Per-segment polynomial fit, then least-squares refinement on the exact arithmetic."""
    xs, lz, xn, e = inv_inputs()
    nseg, h = 1 << d.k, 2.0**-d.k
    table = np.zeros((2 * nseg, d.degree + 1), dtype=np.int64)
    nodes = np.cos(np.pi * (np.arange(4 * d.degree + 8) + 0.5) / (4 * d.degree + 8)) * 0.5 + 0.5
    for ee in (0, 1):
        scale = 1.0 if ee else 1 / math.sqrt(2)
        for s in range(nseg):
            c = np.polynomial.polynomial.polyfit(nodes * h, scale / np.sqrt(1 + s * h + nodes * h), d.degree)
            table[(ee << d.k) | s] = [round(abs(c[j]) * 2**d.frac[j]) for j in range(d.degree + 1)]
    addr = (e << d.k) | ((xn >> (INV_WIDTH - 1 - d.k)) & (nseg - 1))
    dl = (xn & ((1 << (INV_WIDTH - 1 - d.k)) - 1)) / 2.0**21
    target = 1 / np.sqrt(xs.astype(np.float64))
    scale = 2.0**(10 - (lz >> 1))
    order = np.argsort(addr, kind='stable')
    bounds = np.searchsorted(addr[order], np.arange(2 * nseg + 1))
    for it in range(refine_iters + 1):
        err = (inv_phase_turns(d, inv_evaluate(d, table, xn, e), lz) - target) * scale
        for a in range(2 * nseg):
            idx = order[bounds[a]:bounds[a + 1]]
            if len(idx) == 0:
                continue
            if it < refine_iters:
                coef = np.linalg.lstsq(np.vander(dl[idx], d.degree + 1, increasing=True), err[idx], rcond=None)[0]
                for j in range(d.degree + 1):
                    table[a, j] -= round((1 if j % 2 == 0 else -1) * coef[j] * 2**d.frac[j])
            else:
                rel = err[idx] / (target[idx] * scale[idx])
                table[a, 0] -= round((rel.max() + rel.min()) / 2 * (target[idx] * scale[idx]).mean() * 2**d.frac[0])
    y = inv_evaluate(d, table, xn, e)
    phase = inv_phase_turns(d, y, lz)
    stats = dict(
        worst_ratio=float((np.abs(phase - target) / (2.0**-24 * target)).max()),
        min_margin_turns=float((2.0**-24 * target - np.abs(phase - target)).min()),
        low_widths=inv_low_widths(d, table, xn, e),
    )
    return table, y, stats


def inv_low_widths(d, table, xn, e):
    """Bits needed for the total amount P each Horner step subtracts (over all x).
    The borrow trick needs 0 <= P < 2^w and a non-negative result."""
    out = [0] * d.degree
    for j, c, acc in inv_steps(d, table, xn, e):
        P = c - acc
        assert P.min() >= 0 and (acc >= 0).all(), 'Horner range check failed'
        out[j] = int(P.max()).bit_length()
    return out

# %% [markdown]
# ### Circuit builder

# %%
INV_STAGES = (16, 8, 4, 2, 1)   # normalization shifts; flag L means "top L bits of x are 0"


def _inv_flag_zero_top(b, x, L, flag):
    top = x[INV_WIDTH - L:]
    b.x(top)
    if L == 1:
        b.cx(top[0], flag)
    else:
        b.multi_controlled_x(top, flag)
    b.x(top)


def _inv_mbu_clear_flag(b, x, L, flag):
    """Measure flag in the X basis; on outcome 1, apply (-1)^flag computed from x."""
    m = b.alloc_bits(1)[0]
    b.hmr(flag, m)
    top = x[INV_WIDTH - L:]
    b.x(top)
    b.push_condition(m)
    if L == 1:
        b.z(top[0])
    elif L == 2:
        b.cz(top[0], top[1])
    else:
        anc = b.alloc_qubits(L - 2)
        b.init_and(top[0], top[1], anc[0])
        for i in range(1, L - 2):
            b.init_and(anc[i - 1], top[i + 1], anc[i])
        b.cz(anc[L - 3], top[L - 1])
        for i in range(L - 3, 0, -1):
            b.del_and(anc[i - 1], top[i + 1], anc[i])
        b.del_and(top[0], top[1], anc[0])
        b.free(anc)
    b.pop_condition()
    b.x(top)
    b.bit_store(m, False)
    b.free(m)


def inv_normalize(b, x, flags, inverse):
    stages = list(zip(INV_STAGES, flags))
    for L, f in (reversed(stages) if inverse else stages):
        with b.mark(f'{"un" if inverse else ""}norm_{L}'):
            if not inverse:
                _inv_flag_zero_top(b, x, L, f)
                for j in range(INV_WIDTH - 1, L - 1, -1):
                    b.cswap(f, x[j], x[j - L])
            else:
                for j in range(L, INV_WIDTH):
                    b.cswap(f, x[j], x[j - L])
                _inv_mbu_clear_flag(b, x, L, f)


def _inv_booth_row(b, inverse, x, lo, mid, hi, r, s, tgt, scratch, wtop):
    """tgt -= digit * (r >> s) (or += when inverse), digit = -2*hi + mid + lo from x bits.

    Row bit p = one * r[p+s] xor two * r[p+s-1], where one = mid xor lo and
    two = (not one) and (hi xor mid). Bits at or above len(tgt) are dropped:
    tgt is mod-2^len(tgt) arithmetic and the final value is known to fit.
    """
    g = lambda i: x[i] if i is not None else None
    qlo, qmid, qhi = g(lo), g(mid), g(hi)
    R = lambda t: r[t] if 0 <= t < len(r) else None
    o, t = b.alloc_qubits(2)
    for q in (qmid, qlo):                       # o = one
        if q is not None:
            b.cx(q, o)
    both = qhi is not None and qmid is not None
    if both:
        b.cx(qmid, qhi)                         # hi temporarily holds hi xor mid
    hm = qhi if qhi is not None else qmid
    has_two = hm is not None
    if has_two:                                 # t = two
        b.x(o); b.init_and(o, hm, t); b.x(o)
    pos1 = [p for p in range(wtop + 1) if R(p + s) is not None]
    pos2 = [p for p in range(wtop + 1) if R(p + s - 1) is not None] if has_two else []
    positions = sorted(set(pos1) | set(pos2))
    if positions:
        for p in pos1:
            b.init_and(o, R(p + s), scratch[p])
        for p in pos2:
            if p in pos1:
                b.ccx(t, R(p + s - 1), scratch[p])
            else:
                b.init_and(t, R(p + s - 1), scratch[p])
        if both:
            b.cx(qmid, qhi)                     # restore hi: it is the digit's sign
        span = tgt[positions[0]:]
        if qhi is not None:
            b.cx(qhi, span)                     # negative digit: ~(~t - v) = t + v
        if inverse:
            b.iadd(scratch[positions[0]:len(tgt)], target=span)
        else:
            b.isub(scratch[positions[0]:len(tgt)], target=span)
        if qhi is not None:
            b.cx(qhi, span)
        if both:
            b.cx(qmid, qhi)
        mb = b.alloc_bits(len(positions))       # clear row bits by measurement
        for p, m in zip(positions, mb):
            b.hmr(scratch[p], m)
            b.push_condition(m)
            if p in pos1:
                b.cz(o, R(p + s))
            if p in pos2:
                b.cz(t, R(p + s - 1))
            b.pop_condition()
        b.bit_store(mb, False)
        b.free(mb)
    if has_two:
        b.x(o); b.del_and(o, hm, t); b.x(o)
    if both:
        b.cx(qmid, qhi)
    for q in (qlo, qmid):
        if q is not None:
            b.cx(q, o)
    b.free(km.array([o, t]))


def inv_horner(b, d, x, regs, scratch, borrows, low_widths, inverse):
    steps = list(range(d.degree - 1, -1, -1))
    for j in (reversed(steps) if inverse else steps):
        acc, r = regs[j], regs[j + 1]
        fin, fout, wl = d.frac[j + 1], d.frac[j], low_widths[j]
        tgt = km.array(list(acc[:wl]) + [borrows[j]])
        with b.mark(f'horner_{j}'):
            if inverse:
                b.iadd(1, target=acc[wl:], control=borrows[j])
            digits = inv_booth_digits(d, j)
            for lo, mid, hi, m in (reversed(digits) if inverse else digits):
                _inv_booth_row(b, inverse, x, lo, mid, hi, r, m + fin - fout, tgt, scratch, wl)
            if not inverse:
                b.isub(1, target=acc[wl:], control=borrows[j])


def inv_build(d, table, low_widths, compute_only=False, max_qubits=None):
    b = km.CircuitBuilder()
    if max_qubits:
        b.set_max_qubits(max_qubits)   # kickmix adders then use fewer carry qubits
    x = b.create_quantum_register(INV_WIDTH, name='x')
    alloc = (lambda n, name: b.create_quantum_register(n, name=name)) if compute_only else (lambda n, name: b.alloc_qubits(n))
    flags = alloc(len(INV_STAGES), 'flags')          # flags[4] (shift 1) is the parity of t
    widths = d.widths
    coef = alloc(sum(widths), 'coef')
    offs = [sum(widths[:j]) for j in range(len(widths))]
    regs = [coef[offs[j]:offs[j] + widths[j]] for j in range(len(widths))]
    ext = alloc(10, 'ext')
    borrows = alloc(d.degree, 'borrows')
    scratch = alloc(max(low_widths) + 1, 'scratch')
    words = []
    for row in table:
        v = 0
        for j, c in enumerate(row):
            assert 0 <= int(c) < (1 << widths[j]), 'coefficient does not fit its register'
            v |= int(c) << offs[j]
        words.append(v)
    address = km.array([x[INV_WIDTH - 1 - d.k + i] for i in range(d.k)] + [flags[4]])
    sh = km.array([flags[3], flags[2], flags[1], flags[0]])         # sh = t >> 1, little-endian
    yh = km.array(list(regs[0][d.frac[0] - d.pf:]) + list(ext))     # phased part of y, LSB = 2^-pf

    def compute(inverse):
        if not inverse:
            inv_normalize(b, x, flags, False)
            with b.mark('lookup'):
                b.init_lookup(table=words, address=address, target=coef)
            inv_horner(b, d, x, regs, scratch, borrows, low_widths, False)
            with b.mark('denormalize'):
                b.left_rotate(yh, sh)
        else:
            with b.mark('renormalize'):
                b.right_rotate(yh, sh)
            inv_horner(b, d, x, regs, scratch, borrows, low_widths, True)
            with b.mark('unlookup'):
                b.del_lookup(table=words, address=address, target=coef)
            inv_normalize(b, x, flags, True)

    with b.mark('compute'):
        compute(False)
    if compute_only:
        return b.finish_circuit()
    with b.mark('phase'):
        for pos in range(len(yh)):
            w = pos - d.pf - 10                       # this bit is worth 2^w turns
            if w < 0:                                 # whole turns do nothing
                b.z_pow(yh[pos], Fraction(2) ** (w + 1))
    with b.mark('uncompute'):
        compute(True)
    b.free(scratch); b.free(borrows); b.free(ext); b.free(coef); b.free(flags)
    return b.finish_circuit()

# %% [markdown]
# ### Build and check
#
# Three levels of checks, each stricter than the last:
#
# 1. **Model, all inputs.** The exact fixed-point model must stay within 97% of
#    the tolerance on every x (the remaining margin pays for rotation synthesis).
# 2. **Circuit = model.** Simulate the compute half of the circuit on 4,096
#    inputs (all small x, the largest x and random ones) and require the
#    shifted output register to match the model bit for bit.
# 3. **Official verifier.** `quick`: the first 65,536 inputs, 4 shots each.
#    `full`: every input, 64 shots each, exactly as the leaderboard runs it.

# %%
# Found by a local search over widths, scored by the real circuit's qubits x Toffolis.
INV_DESIGN = InvDesign(k=8, frac=(27, 19, 9), dbits=(13, 14), pf=27)
INV_MAX_QUBITS = 128


def inv_check_circuit_matches_model(d, table, y, low_widths):
    xs, lz, xn, e = inv_inputs()
    c = inv_build(d, table, low_widths, compute_only=True)
    rng = np.random.default_rng(0)
    sample = np.unique(np.concatenate([np.arange(1, 2049), rng.integers(1, INV_MAX_X + 1, 2047), [INV_MAX_X]]))
    sim = km.Simulator(batch_size=len(sample))
    sim.use_same_registers_as(c)
    sim.write_across_shots(c.register_data['x'], (sample[None, :] >> np.arange(22)[:, None]) & 1 == 1)
    sim.do(c)
    rd = lambda name: sim.read_across_shots(c.register_data[name]).astype(np.int64)
    coef, ext = rd('coef'), rd('ext')
    g = d.frac[0] - d.pf
    yh_bits = np.concatenate([coef[g:d.widths[0]], ext])
    yh = (yh_bits << np.arange(len(yh_bits))[:, None]).sum(0)
    low = (coef[:g] << np.arange(g)[:, None]).sum(0)
    ys, sh = y[sample - 1], lz[sample - 1] >> 1
    return bool(np.array_equal(yh, (ys >> g) << sh) and np.array_equal(low, ys & ((1 << g) - 1)))


def inv_estimated_score(circuit, stats, mean_toffolis):
    n_rot = INV_DESIGN.pf + 10 - 3
    eps = math.pi * stats['min_margin_turns'] / n_rot
    t = 1 + n_rot * (0.56 * math.log2(1 / eps) + 5.3)
    return circuit.num_qubits * (mean_toffolis + t / 4)


def run_invsqrt():
    from invsqrt import verify_submission_invsqrt as V
    t0 = time.time()
    table, y, stats = inv_fit(INV_DESIGN)
    print(f"model: worst error {stats['worst_ratio']:.3f} x tolerance on all {INV_MAX_X:,} inputs")
    if stats['worst_ratio'] > 0.97:
        raise AssertionError('model error too close to tolerance')
    circuit = inv_build(INV_DESIGN, table, stats['low_widths'], max_qubits=INV_MAX_QUBITS)
    if not inv_check_circuit_matches_model(INV_DESIGN, table, y, stats['low_widths']):
        raise AssertionError('circuit output differs from the model')
    print('circuit matches model bit for bit on 4,096 sample inputs')
    path = SUB / 'invsqrt.kmb'
    circuit.to_file(str(path), format='kmb')
    result = dict(num_qubits=circuit.num_qubits,
                  static_toffolis=circuit.operation_statistics().num_toffoli_like)
    if VERIFY in ('quick', 'full'):
        full = VERIFY == 'full'
        spec = V.ProblemSpecInvSqrt.for_n_p(10) if full else V.ProblemSpecInvSqrt(
            width=22, max_x=65536, beta=Fraction(1), precision_bits=24, shots_per_x=4, system_qubits=0)
        r = V.verify(km.Circuit.from_file(str(path)), spec)
        result.update(verifier='full' if full else 'quick (x <= 65,536)', verifier_status=r['status'],
                      mean_toffolis=r['toffoli_count'])
        result['score'] = r['score'] if full else inv_estimated_score(circuit, stats, r['toffoli_count'])
    else:
        result['score'] = inv_estimated_score(circuit, stats, result['static_toffolis'])
    result['seconds'] = round(time.time() - t0, 1)
    return result


def use_baseline_invsqrt():
    src = SAMPLE / 'invsqrt.kmb'
    if src.exists():
        shutil.copy(src, SUB / 'invsqrt.kmb')
    else:
        from invsqrt.nr_reference_circuit import build_nr_circuit
        build_nr_circuit(10).to_file(str(SUB / 'invsqrt.kmb'), format='kmb')


if USE_IMPROVED['invsqrt']:
    try:
        REPORT['invsqrt'] = dict(version='improved', **run_invsqrt())
    except Exception as ex:
        print('invsqrt: improved circuit rejected, using baseline:', ex)
        use_baseline_invsqrt()
        REPORT['invsqrt'] = dict(version='baseline (fallback)', error=str(ex))
else:
    use_baseline_invsqrt()
    REPORT['invsqrt'] = dict(version='baseline')
print(json.dumps(REPORT['invsqrt'], indent=2))

# %% [markdown]
# ## 3. dqi: Reed-Solomon decoder
#
# Baseline for now: the reference decoder (1,893 qubits, about 5.73M Toffolis).
# It is copied from the sample submission when available, or rebuilt from the
# reference script (about 7 s; the file is about 890 MB before zipping).

# %%
def use_baseline_dqi():
    src = SAMPLE / 'dqi.kmb'
    if src.exists():
        shutil.copy(src, SUB / 'dqi.kmb')
    else:
        from dqi.dqi_reference_circuit import build_reference_circuit
        build_reference_circuit(4095, 70, 12).to_file(str(SUB / 'dqi.kmb'), format='kmb')


t0 = time.time()
use_baseline_dqi()
REPORT['dqi'] = dict(version='baseline', seconds=round(time.time() - t0, 1))
print(REPORT['dqi'])

# %% [markdown]
# ## 4. dfthcsos: chemistry block encoding
#
# Baseline for now: the sample submission's 10 angle circuits and factorization.

# %%
for f in sorted((SAMPLE / 'dfthcsos').iterdir()):
    shutil.copy(f, SUB / 'dfthcsos' / f.name)
REPORT['dfthcsos'] = dict(version='baseline', files=len(list((SUB / 'dfthcsos').iterdir())))
if VERIFY in ('quick', 'full'):
    from dfthcsos.verify import verify as verify_dfthcsos
    r = verify_dfthcsos(str(SUB / 'dfthcsos'))
    REPORT['dfthcsos'].update(verifier_status=r['status'], score=r['score'])
print(REPORT['dfthcsos'])

# %% [markdown]
# ## 5. Package and report
#
# Writes `submission.zip` with the three tasks at the archive root and checks
# its layout. The ratio column is baseline cost / our cost (higher is better);
# the overall score is the geometric mean of the three ratios.

# %%
BASELINE_SCORE = {'invsqrt': 2076774.23, 'dqi': None, 'dfthcsos': 3355599588.3155375}

zip_path = OUT / 'submission.zip'
with zipfile.ZipFile(zip_path, 'w', compression=zipfile.ZIP_DEFLATED, compresslevel=6) as z:
    for f in sorted(SUB.rglob('*')):
        if f.is_file():
            z.write(f, f.relative_to(SUB).as_posix())
with zipfile.ZipFile(zip_path) as z:
    names = set(z.namelist())
need = {'invsqrt.kmb', 'dqi.kmb', 'dfthcsos/submission.npz'}
angles = sorted(n for n in names if n.startswith('dfthcsos/angle_k_'))
assert need <= names, f'missing {need - names}'
assert angles == [f'dfthcsos/angle_k_{i}.kmb' for i in range(len(angles))], 'angle files not contiguous'

ratios = {}
for task, info in REPORT.items():
    if info.get('version', '').startswith('baseline'):
        ratios[task] = 1.0
    elif info.get('score') and BASELINE_SCORE.get(task):
        ratios[task] = BASELINE_SCORE[task] / info['score']
print(f'{"task":10} {"version":22} {"score":>16} {"ratio":>8}')
for task, info in REPORT.items():
    s = info.get('score')
    print(f'{task:10} {info["version"]:22} {(f"{s:,.0f}" if s else "-"):>16} {ratios.get(task, float("nan")):8.2f}')
overall = math.prod(ratios.values()) ** (1 / len(ratios))
print(f'\nestimated overall score: {overall:.3f}')
print(f'{zip_path} ({zip_path.stat().st_size / 1e6:.1f} MB), {len(names)} files')
