"""Our DQI decoder v3: normal-schedule register-shared Euclid + Chien/Forney phase scan.

Layer II-2 (key equation): see eea_model.py, which this builder mirrors word for word.
  35 fixed rounds, each = two quotient coefficients (q1 x + q0). Every position is classically known,
  so there is no pointer logic, no unary iteration and no array shifting: freed words are relabelled.
  Rounds that start after the remainder collapsed (k < 35 errors) are frozen by a sticky flag that
  zeroes the quotients. Abnormal degree drops (about 0.7% of shots) are decoding failures.
Layer II-3 (scan): sigma = A.t, Omega = A.r. Position j: zero test of sigma(a^-j), tower inverse of
  sigma_odd(a^-j), phase (-1)^popcount(M_j & (Omega'(a^-j) / sigma_odd)) via a Toffoli-free bilinear
  phase, gated by the zero flag. Chien updates are constant multiplies (CNOTs).
Layer II-4: the key-equation rounds run in reverse, restoring S and clearing every word.
Field arithmetic: our tower / GF(4)-CRT circuits (multiply 42, divide 119, inverse 77 Toffolis).

    python decoder_v3.py out.kmb
"""
import sys
import kickmix as km
sys.path.insert(0, '/home/claude/google-quantum-ai---verifiable-quantum-circuit-compilation/ours/dqi')
import tower as TW

GF = km.GF2Field(12)
ALPHA = GF.primitive_element
N_SYN, T_ERR, M_POS, DEG = 70, 35, 4095, 12
FLAG_FROM = 19          # freezing is possible from this round on (k >= 18 errors almost surely)


def mul_into(b, lhs, rhs, target):
    TW.tw_ixor_mul(b, lhs, rhs, target)


def del_mul(b, lhs, rhs, target):
    """target holds lhs*rhs (computed into a zero register): clear by measurement, 0 Toffolis."""
    for reg in (lhs, rhs, target):
        TW.tw_apply_linear(b, reg, TW._TW_TO_4)
    TW.tw_del_bilinear(b, TW.TW_MUL_CRT, lhs, rhs, target)
    for reg in (rhs, lhs):
        TW.tw_apply_linear(b, reg, TW._TW_FROM_4)


class Decoder:
    def __init__(self, b, S, E):
        self.b, self.S, self.E = b, S, E

    # ---- quotient handling ----------------------------------------------------------------
    def _make_q(self, num, den, f):
        """q := (num / den) and not f  (f None: ungated). Returns the q register."""
        b = self.b
        ratio, scaf = TW.tw_init_div(b, num, den)
        q = b.alloc_qubits(DEG)
        if f is None:
            b.cx(ratio, q)
        else:
            b.x(f); b.init_and(f, ratio, q); b.x(f)
        TW.tw_del_div(b, num, den, ratio, scaf); b.free(ratio); b.free(scaf)
        return q

    def _clear_q(self, q, num, den, f):
        """q holds (num / den) and not f: clear it."""
        b = self.b
        ratio, scaf = TW.tw_init_div(b, num, den)
        if f is None:
            b.cx(ratio, q)
        else:
            b.x(f); b.del_and(f, ratio, q); b.x(f)
        TW.tw_del_div(b, num, den, ratio, scaf); b.free(ratio); b.free(scaf)
        b.free(q)

    def _clear_q_from_cofactor(self, q1, q0, At, Bt, dt, f):
        """After a round, B.t[dt+1] = q1 A.t[dt] and B.t[dt] = q1 A.t[dt-1] + q0 A.t[dt]."""
        b = self.b
        r1, s1 = TW.tw_init_div(b, Bt[dt + 1], At[dt])
        u = None
        if dt >= 1:
            u = b.alloc_qubits(DEG); mul_into(b, r1, At[dt - 1], u); b.cx(u, Bt[dt])
        self._clear_q(q0, Bt[dt], At[dt], f)
        if u is not None:
            b.cx(u, Bt[dt]); del_mul(b, r1, At[dt - 1], u); b.free(u)
        if f is None:
            b.cx(r1, q1)
        else:
            b.x(f); b.del_and(f, r1, q1); b.x(f)
        TW.tw_del_div(b, Bt[dt + 1], At[dt], r1, s1); b.free(r1); b.free(s1); b.free(q1)

    def _make_q_from_cofactor(self, At, Bt, dt, f):
        """Inverse of _clear_q_from_cofactor: recreate (q1, q0) from the cofactor's top words."""
        b = self.b
        r1, s1 = TW.tw_init_div(b, Bt[dt + 1], At[dt])
        q1 = b.alloc_qubits(DEG)
        if f is None:
            b.cx(r1, q1)
        else:
            b.x(f); b.init_and(f, r1, q1); b.x(f)
        u = None
        if dt >= 1:
            u = b.alloc_qubits(DEG); mul_into(b, r1, At[dt - 1], u); b.cx(u, Bt[dt])
        q0 = self._make_q(Bt[dt], At[dt], f)
        if u is not None:
            b.cx(u, Bt[dt]); del_mul(b, r1, At[dt - 1], u); b.free(u)
        TW.tw_del_div(b, Bt[dt + 1], At[dt], r1, s1); b.free(r1); b.free(s1)
        return q1, q0

    # ---- one round, forward or reverse -------------------------------------------------------
    def _macs_r1(self, q1, Ar, Br, Da):
        for c in range(1, Da + 2):
            mul_into(self.b, q1, Ar[c - 1], Br[c])

    def _macs_r0(self, q0, Ar, Br, Da):
        for c in range(0, Da + 1):
            mul_into(self.b, q0, Ar[c], Br[c])

    def _macs_t(self, q1, q0, At, Bt, dt):
        for c in range(1, dt + 2):
            mul_into(self.b, q1, At[c - 1], Bt[c])
        for c in range(0, dt + 1):
            mul_into(self.b, q0, At[c], Bt[c])

    def _flag(self, i, Ar, Da, prev):
        """f_i = f_{i-1} or [A.r[Da] == 0], in a fresh qubit (reversible: f_{i-1} is kept)."""
        b = self.b
        f = b.alloc_qubits(1)[0]
        if prev is None:
            b.flip_if_equal(Ar[Da], 0, target=f)
        else:
            z = b.alloc_qubits(1)[0]
            b.flip_if_equal(Ar[Da], 0, target=z)
            b.cx(prev, f); b.x(prev); b.ccx(prev, z, f); b.x(prev)
            b.flip_if_equal(Ar[Da], 0, target=z); b.free(z)
        return f

    def _unflag(self, f, Ar, Da, prev):
        b = self.b
        if prev is None:
            b.flip_if_equal(Ar[Da], 0, target=f)
        else:
            z = b.alloc_qubits(1)[0]
            b.flip_if_equal(Ar[Da], 0, target=z)
            b.x(prev); b.ccx(prev, z, f); b.x(prev); b.cx(prev, f)
            b.flip_if_equal(Ar[Da], 0, target=z); b.free(z)
        b.free(f)

    def forward(self):
        b = self.b
        self.one_t = b.alloc_qubits(DEG); b.x(self.one_t[0])
        self.zeros = [b.alloc_qubits(DEG) for _ in range(70)]
        self.one_r = b.alloc_qubits(DEG); b.x(self.one_r[0])
        Ar, At = list(self.S), [self.one_t]
        Br, Bt = list(self.zeros) + [self.one_r], []
        self.flags = []
        prev = None
        for i in range(1, T_ERR + 1):
            Da, dt = 70 - i, i - 1
            with b.mark('round'):
                f = None
                if i >= FLAG_FROM:
                    f = self._flag(i, Ar, Da, prev); prev = f
                self.flags.append(f)
                q1 = self._make_q(Br[Da + 1], Ar[Da], f)
                self._macs_r1(q1, Ar, Br, Da)
                q0 = self._make_q(Br[Da], Ar[Da], f)
                self._macs_r0(q0, Ar, Br, Da)
                top1 = Br.pop(); top0 = Br.pop(); Bt = Bt + [top0, top1]
                self._macs_t(q1, q0, At, Bt, dt)
                self._clear_q_from_cofactor(q1, q0, At, Bt, dt, f)
            Ar, At, Br, Bt = Br, Bt, Ar, At
        # parity of frozen rounds -> one controlled swap brings sigma/Omega to canonical words
        self.parity = b.alloc_qubits(1)[0]
        for f in self.flags:
            if f is not None:
                b.cx(f, self.parity)
        self._parity_swap(Ar, At, Br, Bt)
        self.state = (Ar, At, Br, Bt)

    def _parity_swap(self, Ar, At, Br, Bt):
        pairs = [(Ar[c], Br[c]) for c in range(35)] + [(At[c], Bt[c]) for c in range(35)] + [(At[35], Br[35])]
        lhs = km.array([q for a, _ in pairs for q in a]); rhs = km.array([q for _, c in pairs for q in c])
        self.b.cswap(self.parity, lhs, rhs)

    def reverse(self):
        b = self.b
        Ar, At, Br, Bt = self.state
        self._parity_swap(Ar, At, Br, Bt)
        for f in self.flags:
            if f is not None:
                b.cx(f, self.parity)
        b.free(self.parity)
        for i in range(T_ERR, 0, -1):
            Ar, At, Br, Bt = Br, Bt, Ar, At          # undo the role swap
            Da, dt = 70 - i, i - 1
            f = self.flags[i - 1]
            prev = self.flags[i - 2] if i >= 2 else None
            with b.mark('unround'):
                q1, q0 = self._make_q_from_cofactor(At, Bt, dt, f)
                self._macs_t(q1, q0, At, Bt, dt)
                top0, top1 = Bt[-2], Bt[-1]; Bt = Bt[:-2]; Br = Br + [top0, top1]
                self._macs_r0(q0, Ar, Br, Da)
                self._clear_q(q0, Br[Da], Ar[Da], f)
                self._macs_r1(q1, Ar, Br, Da)
                self._clear_q(q1, Br[Da + 1], Ar[Da], f)
                if f is not None:
                    self._unflag(f, Ar, Da, prev)
        assert Ar == list(self.S) and At == [self.one_t] and Br == self.zeros + [self.one_r] and Bt == []
        b.x(self.one_r[0]); b.x(self.one_t[0])
        for w in self.zeros + [self.one_r, self.one_t]:
            b.free(w)

    # ---- Chien / Forney phase scan ----------------------------------------------------------
    def scan(self):
        b = self.b
        Ar, At, _, _ = self.state
        sig, om = At, Ar                     # sigma_c = At[c] (c <= 35), Omega_c = Ar[c] (c <= 34)
        mult_sig = [GF.pow(ALPHA, (-c) % M_POS) for c in range(36)]
        mult_om = [GF.pow(ALPHA, (-(c + 1)) % M_POS) for c in range(35)]
        # sums: sig[1] <- sigma_odd, sig[0] <- sigma, om[0] <- Omega'(x) = x Omega(x)
        pairs = ([(sig[c], sig[1]) for c in range(3, 36, 2)] + [(sig[c], sig[0]) for c in range(2, 36, 2)]
                 + [(sig[1], sig[0])] + [(om[c], om[0]) for c in range(1, 35)])
        fs = km.array([q for s, _ in pairs for q in s]); fd = km.array([q for _, d in pairs for q in d])
        rs = km.array([q for s, _ in reversed(pairs) for q in s]); rd = km.array([q for _, d in reversed(pairs) for q in d])
        is_zero = b.alloc_qubits(1)[0]; located = b.alloc_qubits(DEG); is_err = b.alloc_bits(1)
        step = b.mark('forney_step')
        for j in range(M_POS):
            with step:
                b.bit_store(is_err[0], True, control=self.E[j])
                b.push_condition(is_err[0])
                b.cx(fs, fd)
                b.flip_if_equal(sig[0], 0, target=is_zero)
                inv, scaf = TW.tw_init_inverse(b, sig[1])
                b.init_and(is_zero, om[0], located)
                b.gf2_phase_by_product(located, inv, mask=self.E[j], field=GF)
                b.del_and(is_zero, om[0], located)
                TW.tw_del_inverse(b, sig[1], inv, scaf); b.free(inv); b.free(scaf)
                b.flip_if_equal(sig[0], 0, target=is_zero)
                b.cx(rs, rd)
                b.pop_condition()
                b.bit_store(is_err, False)
                if j + 1 < M_POS:
                    for c in range(1, 36):
                        b.gf2_imul(mult_sig[c], target=sig[c], field=GF)
                    for c in range(35):
                        b.gf2_imul(mult_om[c], target=om[c], field=GF)
        for c in range(1, 36):
            b.gf2_idiv(GF.pow(mult_sig[c], M_POS - 1), target=sig[c], field=GF)
        for c in range(35):
            b.gf2_idiv(GF.pow(mult_om[c], M_POS - 1), target=om[c], field=GF)
        b.free(is_err); b.free(located); b.free(is_zero)


def build(m=M_POS, n=N_SYN):
    b = km.CircuitBuilder()
    S = [b.create_quantum_register(DEG, name=f'S{i}') for i in range(n)]
    E = [b.create_classical_register(DEG, name=f'E{j}') for j in range(m)]
    d = Decoder(b, S, E)
    with b.mark('key_equation'):
        d.forward()
    with b.mark('scan'):
        d.scan()
    with b.mark('unkey_equation'):
        d.reverse()
    return b.finish_circuit()


if __name__ == '__main__':
    import time
    t0 = time.time()
    c = build()
    print(f'qubits {c.num_qubits}, static Toffolis {c.operation_statistics().num_toffoli_like:,} ({time.time() - t0:.0f}s)')
    c.to_file(sys.argv[1] if len(sys.argv) > 1 else 'decoder_v3.kmb', format='kmb')
