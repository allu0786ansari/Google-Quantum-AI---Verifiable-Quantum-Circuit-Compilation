# %% [markdown]
# ## 3. dqi: Reed-Solomon decoder
#
# **Contract.** Syndromes S0..S69 (12-bit quantum registers over GF(2^12)) and a
# classical record E0..E4094. Decode the error pattern implicitly, apply the sign
# (-1)^(sum popcount(e_j & E_j)), restore the syndromes and clean all workspace.
# Tested on 102,400 random error patterns. Cost = R(eps+) x qubits^2 x Toffolis.
#
# **Baseline** (1,893 qubits, about 5.73M Toffolis): synchronized extended Euclid,
# then a Chien search with one field inversion at each of the 4,095 positions.
# GF(2^12) multiplications (63 Toffolis each) and inversions (315 each) are about
# 59% of all Toffolis.
#
# **Our change: tower-field arithmetic.** We keep Google's reference decoder
# (imported unchanged from the competition data) and redirect its field
# operations to cheaper circuits:
#
# | Operation | Reference | Ours |
# |---|---|---|
# | inverse | 315 | **77** |
# | divide (inverse + multiply) | 378 | **119** |
# | controlled multiply | 75 | **54** |
#
# Changing basis is a fixed linear map, so it costs only CNOTs. Two other views of
# GF(2^12) make the arithmetic cheaper:
#
# - **Inverse:** in the tower GF(((2^3)^2)^2), an inverse is a norm in GF(64), an
#   inverse in GF(8) and two GF(64) multiplies.
# - **Multiply:** as degree-6 polynomials over GF(4), a product is fixed by its
#   values at the 4 points of GF(4), its leading coefficient, and its remainders
#   mod 3 irreducible quadratics (Chinese remainder theorem). That is 14 GF(4)
#   multiplies of 3 ANDs each: 42 AND terms instead of 63.
#
# Each AND term is one Toffoli. Every intermediate value is bilinear in registers
# that are still alive, so it is cleared by X-basis measurement plus CZ fixes, for
# zero Toffolis. The multiply needs no scratch qubits: linear forms are computed in
# place and restored.
#
# Each piece is checked against kickmix's own GF(2^12) arithmetic before the decoder
# is built.

# %%
TW_F = km.GF2Field(12)       # the competition field: alpha^12 = alpha^6 + alpha^4 + alpha + 1


def tw_vec(v, n):
    return np.array([(v >> i) & 1 for i in range(n)], dtype=np.uint8)


def tw_int(v):
    return int(sum(int(b) << i for i, b in enumerate(v)))


def tw_gf2_solve(A, rhs):
    """One solution x of A x = rhs over GF(2), or None."""
    m, n = A.shape
    M = np.concatenate([A % 2, rhs.reshape(-1, 1) % 2], axis=1).astype(np.uint8)
    r, pivots = 0, []
    for c in range(n):
        piv = next((i for i in range(r, m) if M[i, c]), None)
        if piv is None:
            continue
        M[[r, piv]] = M[[piv, r]]
        for i in range(m):
            if i != r and M[i, c]:
                M[i] ^= M[r]
        pivots.append(c)
        r += 1
    if any(M[i, n] for i in range(r, m)):
        return None
    x = np.zeros(n, dtype=np.uint8)
    for i, c in enumerate(pivots):
        x[c] = M[i, n]
    return x


def tw_gf2_inv(M):
    n = M.shape[0]
    cols = [tw_gf2_solve(M, np.eye(n, dtype=np.uint8)[:, j]) for j in range(n)]
    if any(c is None for c in cols):
        raise ValueError('singular matrix')
    return np.array(cols, dtype=np.uint8).T


def _tw_root(coeffs, candidates):
    """A root of sum c_i X^i (field-element coefficients) among candidates."""
    for x in candidates:
        acc, p = 0, 1
        for c in coeffs:
            acc ^= TW_F.mul(c, p)
            p = TW_F.mul(p, x)
        if acc == 0:
            return x
    raise ValueError('no root found')


def _tw_trace(x, degree):
    acc, p = 0, x
    for _ in range(degree):
        acc ^= p
        p = TW_F.mul(p, p)
    return acc


# Tower elements inside GF(4096): GF(8) = {alpha^(585 m)}, GF(64) = {alpha^(65 m)}.
_a = TW_F.primitive_element
_gamma = _tw_root([1, 1, 0, 1], [TW_F.pow(_a, 585 * m) for m in range(7)])   # gamma^3 = gamma + 1
_GAMMA = [1, _gamma, TW_F.mul(_gamma, _gamma)]
_gf8 = lambda bits: (_GAMMA[0] if bits & 1 else 0) ^ (_GAMMA[1] if bits & 2 else 0) ^ (_GAMMA[2] if bits & 4 else 0)
TW_MU_BITS = next(m for m in range(1, 8) if _tw_trace(_gf8(m), 3) == 1)
_z = _tw_root([_gf8(TW_MU_BITS), 1, 1], [TW_F.pow(_a, 65 * m) for m in range(63)])  # z^2 = z + mu
_gf64 = lambda bits: _gf8(bits & 7) ^ TW_F.mul(_gf8(bits >> 3), _z)
TW_NU_BITS = next(m for m in range(1, 64) if _tw_trace(_gf64(m), 6) == 1)
_Z = _tw_root([_gf64(TW_NU_BITS), 1, 1], range(1, 4096))                         # Z^2 = Z + nu

# Tower coordinate bit 6i + 3j + k is the coefficient of gamma^k z^j Z^i.
_BASIS = [TW_F.mul(TW_F.mul(_GAMMA[k], [1, _z][j]), [1, _Z][i]) for i in range(2) for j in range(2) for k in range(3)]
TW_T = np.array([tw_vec(e, 12) for e in _BASIS], dtype=np.uint8).T     # tower -> polynomial basis
TW_TI = tw_gf2_inv(TW_T)                                              # polynomial -> tower

tw_to_tower = lambda x: tw_int(TW_TI @ tw_vec(x, 12) % 2)
tw_from_tower = lambda t: tw_int(TW_T @ tw_vec(t, 12) % 2)
tw_mul = lambda a, b: tw_to_tower(TW_F.mul(tw_from_tower(a), tw_from_tower(b)))   # in tower coordinates


def tw_matrix(fn, n):
    """GF(2) matrix of a linear map on n-bit tower coordinates."""
    return np.array([tw_vec(fn(1 << i), n) for i in range(n)], dtype=np.uint8).T


def tw_solve_formula(products, n, mul_fn):
    """Given product pairs (u_k, v_k), find output masks w_k with
    sum_k (u_k . a)(v_k . b) w_k = mul_fn(a, b). Returns [(u, v, w)]."""
    B = np.zeros((n, n, n), dtype=np.uint8)
    for j, l in itertools.product(range(n), range(n)):
        B[:, j, l] = tw_vec(mul_fn(1 << j, 1 << l), n)
    A = np.array([np.outer(u, v).reshape(-1) % 2 for u, v in products], dtype=np.uint8).T
    W = np.zeros((len(products), n), dtype=np.uint8)
    for i in range(n):
        x = tw_gf2_solve(A, B[i].reshape(-1))
        if x is None:
            raise ValueError('products do not span the multiplication')
        W[:, i] = x
    return [(u, v, W[k]) for k, (u, v) in enumerate(products) if W[k].any()]


def tw_karatsuba(n_sub, pairs):
    """Lift product pairs to the quadratic extension: a0 b0, a1 b1, (a0 + a1)(b0 + b1)."""
    zero, out = np.zeros(n_sub, dtype=np.uint8), []
    for u, v in pairs:
        out += [(np.concatenate([u, zero]), np.concatenate([v, zero])),
                (np.concatenate([zero, u]), np.concatenate([zero, v])),
                (np.concatenate([u, u]), np.concatenate([v, v]))]
    return out


_e3 = lambda *i: np.array([1 if t in i else 0 for t in range(3)], dtype=np.uint8)
_sub = lambda n: (lambda a, b: tw_mul(a, b) & ((1 << n) - 1))
TW_MUL8 = tw_solve_formula([(_e3(0), _e3(0)), (_e3(1), _e3(1)), (_e3(2), _e3(2)),
                            (_e3(0, 1), _e3(0, 1)), (_e3(0, 2), _e3(0, 2)), (_e3(1, 2), _e3(1, 2))], 3, _sub(3))
TW_MUL64 = tw_solve_formula(tw_karatsuba(3, [(u, v) for u, v, _ in TW_MUL8]), 6, _sub(6))
TW_MUL4096 = tw_solve_formula(tw_karatsuba(6, [(u, v) for u, v, _ in TW_MUL64]), 12, tw_mul)
print(f'tower: mu={TW_MU_BITS}, nu={TW_NU_BITS}; multiply terms GF(8) {len(TW_MUL8)}, '
      f'GF(64) {len(TW_MUL64)}, GF(4096) {len(TW_MUL4096)}')


# GF(4096) as polynomials of degree < 6 over GF(4) = {0, 1, w, w^2}, w = alpha^1365.
# Coordinate bit 2i + c is the coefficient of w^c theta^i, with theta = alpha.
_w = TW_F.pow(_a, 1365)
_W4 = [0, 1, _w, _w ^ 1]
_g4 = lambda a, b: _W4.index(TW_F.mul(_W4[a], _W4[b]))          # GF(4) multiply on 2-bit codes
TW_T4 = np.array([tw_vec(TW_F.mul([1, _w][c], TW_F.pow(_a, i)), 12) for i in range(6) for c in range(2)],
                 dtype=np.uint8).T                                # GF(4)-poly coords -> polynomial basis
TW_T4I = tw_gf2_inv(TW_T4)
tw_to4 = lambda x: tw_int(TW_T4I @ tw_vec(x, 12) % 2)
tw_from4 = lambda c: tw_int(TW_T4 @ tw_vec(c, 12) % 2)


def _crt_rows(fn, n_out):
    """GF(2) matrix (n_out x 12) of a GF(4)-linear map given on coefficient i, code c."""
    M = np.zeros((n_out, 12), dtype=np.uint8)
    for bit in range(12):
        i, c = divmod(bit, 2)
        M[:, bit] = fn(i, 1 << c)
    return M


def _crt_eval(t):            # A(t) for t in GF(4): 2 bits
    def f(i, c):
        p = 1
        for _ in range(i):
            p = _g4(p, t)
        v = _g4(c, p)
        return [v & 1, v >> 1]
    return _crt_rows(f, 2)


def _crt_mod(q0, q1):        # A(x) mod (x^2 + q1 x + q0): 4 bits
    def f(i, c):
        r0, r1 = 1, 0
        for _ in range(i):
            r0, r1 = _g4(r1, q0), _g4(r1, q1) ^ r0
        v0, v1 = _g4(c, r0), _g4(c, r1)
        return [v0 & 1, v0 >> 1, v1 & 1, v1 >> 1]
    return _crt_rows(f, 4)


def _gf4_pairs(A, B):        # Karatsuba GF(4) product of two 2-bit linear images: 3 ANDs
    return [(A[0], B[0]), (A[1], B[1]), ((A[0] + A[1]) % 2, (B[0] + B[1]) % 2)]


_lead = np.zeros((2, 12), dtype=np.uint8); _lead[0, 10] = _lead[1, 11] = 1
_irreducible = [(q0, q1) for q0 in range(1, 4) for q1 in range(4)
                if all(_g4(t, t) ^ _g4(q1, t) ^ q0 for t in range(4))]
_mul4 = lambda a, b: tw_to4(TW_F.mul(tw_from4(a), tw_from4(b)))
for _qs in itertools.combinations(_irreducible, 3):
    _pairs = []
    for _t in range(4):
        _pairs += _gf4_pairs(_crt_eval(_t), _crt_eval(_t))
    _pairs += _gf4_pairs(_lead, _lead)
    for _q in _qs:
        P = _crt_mod(*_q)
        _pairs += (_gf4_pairs(P[0:2], P[0:2]) + _gf4_pairs(P[2:4], P[2:4])
                   + _gf4_pairs((P[0:2] + P[2:4]) % 2, (P[0:2] + P[2:4]) % 2))
    try:
        TW_MUL_CRT = tw_solve_formula(_pairs, 12, _mul4)
        break
    except ValueError:
        continue
print(f'GF(4)-CRT multiply: {len(TW_MUL_CRT)} AND terms (quadratics {_qs})')

# %% [markdown]
# ### Circuits
#
# - `tw_apply_linear`: an in-place GF(2) linear map (CNOTs and SWAPs from Gaussian elimination).
# - `tw_xor_bilinear`: out ^= sum of AND terms, one Toffoli per term.
# - `tw_del_bilinear`: clears such a register by measurement plus conditioned CZs.
# - `TowerBuilder`: wraps kickmix's builder and redirects the four GF(2^12) operations
#   the reference decoder uses.

# %%
def tw_cnot_ops(M):
    """CNOT/SWAP ops, in application order, mapping v -> M v in place."""
    n, A, ops = M.shape[0], M.copy() % 2, []
    for c in range(n):
        piv = next((i for i in range(c, n) if A[i, c]), None)
        if piv is None:
            raise ValueError('singular matrix')
        if piv != c:
            A[[c, piv]] = A[[piv, c]]
            ops.append(('swap', c, piv))
        for i in range(n):
            if i != c and A[i, c]:
                A[i] ^= A[c]
                ops.append(('cx', c, i))
    return list(reversed(ops))


_TW_TO_TOWER, _TW_TO_POLY = tw_cnot_ops(TW_TI), tw_cnot_ops(TW_T)
_TW_TO_4, _TW_FROM_4 = tw_cnot_ops(TW_T4I), tw_cnot_ops(TW_T4)


def tw_apply_linear(b, qubits, ops):
    for op, j, i in ops:
        (b.swap if op == 'swap' else b.cx)(qubits[j], qubits[i])


def tw_xor_linear(b, M, src, dst):
    """dst ^= M src."""
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            if M[i, j]:
                b.cx(src[j], dst[i])


def _tw_form(b, mask, qubits):
    idx = [j for j in range(len(mask)) if mask[j]]
    if len(idx) == 1:
        return qubits[idx[0]], None
    p = b.alloc_qubits(1)[0]
    for j in idx:
        b.cx(qubits[j], p)
    return p, idx


def _tw_unform(b, p, idx, qubits):
    if idx is not None:
        for j in idx:
            b.cx(qubits[j], p)
        b.free(p)


def tw_xor_bilinear(b, terms, L1, L2, out):
    """out ^= sum_k (u_k . L1)(v_k . L2) w_k with forms in ancillas (L1 may equal L2)."""
    r = b.alloc_qubits(1)[0]
    for u, v, w in terms:
        p, pi = _tw_form(b, u, L1)
        q, qi = _tw_form(b, v, L2)
        targets = [out[i] for i in range(len(w)) if w[i]]
        if p == q:
            for t in targets:
                b.cx(p, t)
        elif len(targets) == 1:
            b.ccx(p, q, targets[0])
        else:
            b.init_and(p, q, r)
            for t in targets:
                b.cx(r, t)
            b.del_and(p, q, r)
        _tw_unform(b, q, qi, L2)
        _tw_unform(b, p, pi, L1)
    b.free(r)


def tw_xor_bilinear_inplace(b, terms, L1, L2, out):
    """Same, with no ancillas (L1, L2, out disjoint): forms are accumulated into a pivot
    qubit of their own register and restored; a product reaching several outputs is
    fanned out through the first one (t_i ^= t0; t0 ^= P; t_i ^= t0)."""
    for u, v, w in terms:
        ui = [j for j in range(len(u)) if u[j]]
        vi = [j for j in range(len(v)) if v[j]]
        for j in ui[1:]:
            b.cx(L1[j], L1[ui[0]])
        for j in vi[1:]:
            b.cx(L2[j], L2[vi[0]])
        t = [out[i] for i in range(len(w)) if w[i]]
        for ti in t[1:]:
            b.cx(t[0], ti)
        b.ccx(L1[ui[0]], L2[vi[0]], t[0])
        for ti in t[1:]:
            b.cx(t[0], ti)
        for j in reversed(vi[1:]):
            b.cx(L2[j], L2[vi[0]])
        for j in reversed(ui[1:]):
            b.cx(L1[j], L1[ui[0]])


def tw_del_bilinear(b, terms, L1, L2, out):
    """Clear out (holding exactly the bilinear form) by X-basis measurement + CZ fixes."""
    mb = b.alloc_bits(len(out))
    b.hmr(km.array(list(out)), mb)
    for u, v, w in terms:
        p, pi = _tw_form(b, u, L1)
        q, qi = _tw_form(b, v, L2)
        for i in range(len(w)):
            if w[i]:
                b.push_condition(mb[i])
                if p == q:
                    b.z(p)
                else:
                    b.cz(p, q)
                b.pop_condition()
        _tw_unform(b, q, qi, L2)
        _tw_unform(b, p, pi, L1)
    b.bit_store(mb, False)
    b.free(mb)


def _tw_lift(terms, n_u, off_u=0, Mu=None, n_v=None, off_v=0, Mv=None):
    out = []
    for u, v, w in terms:
        uu = (Mu.T @ u % 2) if Mu is not None else u
        vv = (Mv.T @ v % 2) if Mv is not None else v
        U = np.zeros(n_u, dtype=np.uint8); U[off_u:off_u + len(uu)] = uu
        V = np.zeros(n_v or n_u, dtype=np.uint8); V[off_v:off_v + len(vv)] = vv
        out.append((U, V, w))
    return out


def _tw_conj_times(terms, n_half):
    """((x0 + x1) + x1 Y) * y for an extension of degree 2 over terms on n_half bits."""
    zero = np.zeros(n_half, dtype=np.uint8)
    lo = [(np.concatenate([u[:n_half], u[:n_half]]), v, np.concatenate([w, zero])) for u, v, w in terms]
    hi = [(np.concatenate([zero, u[:n_half]]), v, np.concatenate([zero, w])) for u, v, w in terms]
    return lo + hi


TW_SQ64 = tw_matrix(lambda x: tw_mul(x, x) & 63, 6)
TW_SQ8 = tw_matrix(lambda x: tw_mul(x, x) & 7, 3)
_NU = tw_matrix(lambda x: tw_mul(x, TW_NU_BITS) & 63, 6)
_MU = tw_matrix(lambda x: tw_mul(x, TW_MU_BITS) & 7, 3)
# Inverse of a = a0 + a1 Z (a_i in GF(64)), all in tower coordinates:
TW_STEP_N = _tw_lift(TW_MUL64, 12, 0, None, 12, 6)                    # N = a0 a1 + a0^2 + nu a1^2
TW_LIN_N = np.concatenate([TW_SQ64, _NU @ TW_SQ64 % 2], axis=1)
TW_STEP_N6 = _tw_lift(TW_MUL8, 6, 0, None, 6, 3)                      # N6 = norm of N down to GF(8)
TW_LIN_N6 = np.concatenate([TW_SQ8, _MU @ TW_SQ8 % 2], axis=1)
TW_STEP_I3 = _tw_lift(TW_MUL8, 3, 0, TW_SQ8 @ TW_SQ8 % 2, 3, 0, TW_SQ8)  # I3 = N6^6 = N6^-1
TW_STEP_I6 = _tw_conj_times(_tw_lift(TW_MUL8, 6, 0, None, 3, 0), 3)   # I6 = N^-1
TW_STEP_R = _tw_conj_times(_tw_lift(TW_MUL64, 12, 0, None, 6, 0), 6)  # R = a^-1
TW_INV_TOFFOLIS = sum(len(s) for s in (TW_STEP_N, TW_STEP_N6, TW_STEP_I3, TW_STEP_I6, TW_STEP_R))


def tw_init_inverse(b, x, target=None):
    """target := x^-1 (polynomial basis, 0 -> 0). Returns (target, scaffold)."""
    N_, N6, I3, I6 = (b.alloc_qubits(n) for n in (6, 3, 3, 6))
    R = b.alloc_qubits(12) if target is None else target
    tw_apply_linear(b, x, _TW_TO_TOWER)
    tw_xor_bilinear(b, TW_STEP_N, x, x, N_); tw_xor_linear(b, TW_LIN_N, x, N_)
    tw_xor_bilinear(b, TW_STEP_N6, N_, N_, N6); tw_xor_linear(b, TW_LIN_N6, N_, N6)
    tw_xor_bilinear(b, TW_STEP_I3, N6, N6, I3)
    tw_xor_bilinear(b, TW_STEP_I6, N_, I3, I6)
    tw_xor_bilinear(b, TW_STEP_R, x, I6, R)
    tw_apply_linear(b, R, _TW_TO_POLY)
    tw_apply_linear(b, x, _TW_TO_POLY)
    return R, km.array(list(N_) + list(N6) + list(I3) + list(I6))


def tw_del_inverse(b, x, R, scaffold):
    """Clear R and its scaffold with zero Toffolis (frees nothing)."""
    N_, N6, I3, I6 = scaffold[0:6], scaffold[6:9], scaffold[9:12], scaffold[12:18]
    tw_apply_linear(b, x, _TW_TO_TOWER)
    tw_apply_linear(b, R, _TW_TO_TOWER)
    tw_del_bilinear(b, TW_STEP_R, x, I6, R)
    tw_del_bilinear(b, TW_STEP_I6, N_, I3, I6)
    tw_del_bilinear(b, TW_STEP_I3, N6, N6, I3)
    tw_xor_linear(b, TW_LIN_N6, N_, N6); tw_del_bilinear(b, TW_STEP_N6, N_, N_, N6)
    tw_xor_linear(b, TW_LIN_N, x, N_); tw_del_bilinear(b, TW_STEP_N, x, x, N_)
    tw_apply_linear(b, x, _TW_TO_POLY)


def tw_init_div(b, num, den):
    """ratio := num / den (den = 0 gives 0, like kickmix). Returns (ratio, scaffold)."""
    R, s = tw_init_inverse(b, den)
    ratio = b.alloc_qubits(12)
    for reg in (num, R):
        tw_apply_linear(b, reg, _TW_TO_4)
    tw_xor_bilinear_inplace(b, TW_MUL_CRT, num, R, ratio)
    for reg in (ratio, R, num):
        tw_apply_linear(b, reg, _TW_FROM_4)
    return ratio, km.array(list(R) + list(s))


def tw_del_div(b, num, den, ratio, scaffold):
    R, s = scaffold[:12], scaffold[12:]
    for reg in (num, R, ratio):
        tw_apply_linear(b, reg, _TW_TO_4)
    tw_del_bilinear(b, TW_MUL_CRT, num, R, ratio)
    for reg in (R, num):
        tw_apply_linear(b, reg, _TW_FROM_4)
    tw_del_inverse(b, den, R, s)


def tw_ixor_mul(b, lhs, rhs, target, control=True):
    """target ^= control * lhs * rhs: 42 Toffolis, plus 12 when controlled."""
    src, tmp = lhs, None
    if control is not True:
        tmp = b.alloc_qubits(12)
        b.init_and(control, lhs, tmp)
        src = tmp
    for reg in (src, rhs, target):
        tw_apply_linear(b, reg, _TW_TO_4)
    tw_xor_bilinear_inplace(b, TW_MUL_CRT, src, rhs, target)
    for reg in (target, rhs, src):
        tw_apply_linear(b, reg, _TW_FROM_4)
    if tmp is not None:
        b.del_and(control, lhs, tmp)
        b.free(tmp)


class TowerBuilder:
    """A km.CircuitBuilder that sends GF(2^12) inverse, division and multiplication to
    the tower circuits above; every other call goes to the real builder."""

    def __init__(self, inner):
        self._b = inner

    def __getattr__(self, name):
        return getattr(self._b, name)

    def init_gf2_inverse_with_scaffold(self, input, *, target='alloc', field=None):
        assert len(input) == 12
        return tw_init_inverse(self._b, input, None if isinstance(target, str) else target)

    def del_gf2_inverse_with_scaffold(self, input, *, target, scaffold, field=None):
        tw_del_inverse(self._b, input, target, scaffold)

    def init_gf2_div_with_scaffold(self, num, den, *, field=None):
        assert len(num) == 12
        return tw_init_div(self._b, num, den)

    def del_gf2_div_with_scaffold(self, num, den, *, target, scaffold, field=None):
        tw_del_div(self._b, num, den, target, scaffold)

    def ixor_gf2_mul(self, lhs, rhs, *, target, field=None, control=True):
        assert len(lhs) == 12
        tw_ixor_mul(self._b, lhs, rhs, target, control)

# %% [markdown]
# ### Unit checks, build and verify
#
# 1. **Formulas:** GF(8) and GF(64) multiplies on every pair; both GF(4096)
#    multiplies (tower and CRT) on 2,000 random pairs, against kickmix's field.
# 2. **Circuits:** inverse on all 4,096 inputs, and divide and controlled multiply on
#    8,192 random inputs each. Each must give the right value, restore its inputs,
#    leave clean workspace and no phase.
# 3. **Decoder:** `quick` runs the official verifier on 1,024 error patterns and
#    requires zero failures; `full` runs it exactly as the leaderboard does
#    (102,400 patterns, about 30-60 min).

# %%
def tw_check_formulas():
    for n, terms in ((3, TW_MUL8), (6, TW_MUL64)):
        for a in range(1 << n):
            for b_ in range(1 << n):
                av, bv, out = tw_vec(a, n), tw_vec(b_, n), np.zeros(n, dtype=np.uint8)
                for u, v, w in terms:
                    if (u @ av) % 2 and (v @ bv) % 2:
                        out ^= w
                assert tw_int(out) == _sub(n)(a, b_), f'GF(2^{n}) formula'
    rng = np.random.default_rng(1)
    for a, b_ in rng.integers(0, 4096, (2000, 2)):
        av, bv, out = tw_vec(tw_to_tower(int(a)), 12), tw_vec(tw_to_tower(int(b_)), 12), np.zeros(12, dtype=np.uint8)
        for u, v, w in TW_MUL4096:
            if (u @ av) % 2 and (v @ bv) % 2:
                out ^= w
        assert tw_from_tower(tw_int(out)) == TW_F.mul(int(a), int(b_)), 'GF(4096) tower formula'
        av, bv, out = tw_vec(tw_to4(int(a)), 12), tw_vec(tw_to4(int(b_)), 12), np.zeros(12, dtype=np.uint8)
        for u, v, w in TW_MUL_CRT:
            if (u @ av) % 2 and (v @ bv) % 2:
                out ^= w
        assert tw_from4(tw_int(out)) == TW_F.mul(int(a), int(b_)), 'GF(4096) CRT formula'


def tw_check_circuit(kind):
    """Simulate one operation on many inputs; check value, inputs, workspace and phase."""
    b = km.CircuitBuilder()
    regs = {n: b.create_quantum_register(12, n) for n in ('a', 'x', 'out')}
    ctl = b.create_quantum_register(1, 'c')
    if kind == 'inverse':
        R, s = tw_init_inverse(b, regs['a']); b.cx(R, regs['out']); tw_del_inverse(b, regs['a'], R, s)
        b.free(R); b.free(s)
    elif kind == 'divide':
        R, s = tw_init_div(b, regs['a'], regs['x']); b.cx(R, regs['out']); tw_del_div(b, regs['a'], regs['x'], R, s)
        b.free(R); b.free(s)
    else:
        tw_ixor_mul(b, regs['a'], regs['x'], regs['out'], control=ctl[0])
    c = b.finish_circuit()
    rng = np.random.default_rng(2)
    n = 4096 if kind == 'inverse' else 8192
    av = np.arange(4096) if kind == 'inverse' else rng.integers(0, 4096, n)
    xv, ov, cv = rng.integers(0, 4096, n), rng.integers(0, 4096, n), rng.integers(0, 2, n)
    if kind == 'divide':
        xv[:64] = 0
    if kind != 'multiply':
        ov[:] = 0
    sim = km.Simulator(batch_size=n)
    sim.use_same_registers_as(c)
    for name, v, w in (('a', av, 12), ('x', xv, 12), ('out', ov, 12), ('c', cv, 1)):
        sim.write_across_shots(c.register_data[name], (v[None, :] >> np.arange(w)[:, None]) & 1 == 1)
    sim.do(c)
    rd = lambda name: (sim.read_across_shots(c.register_data[name]).astype(np.int64) << np.arange(12)[:, None]).sum(0)
    if kind == 'inverse':
        exp = np.array([TW_F.invert(int(v)) if v else 0 for v in av])
    elif kind == 'divide':
        exp = np.array([TW_F.mul(int(p), TW_F.invert(int(d))) if d else 0 for p, d in zip(av, xv)])
    else:
        exp = ov ^ np.array([TW_F.mul(int(p), int(q)) if cc else 0 for p, q, cc in zip(av, xv, cv)])
    named = set(regs['a']) | set(regs['x']) | set(regs['out']) | set(ctl)
    work = km.array([km.q(i) for i in range(c.num_qubits) if km.q(i) not in named])
    assert np.array_equal(rd('out'), exp), f'{kind}: wrong value'
    assert np.array_equal(rd('a'), av) and np.array_equal(rd('x'), xv), f'{kind}: inputs changed'
    assert not (len(work) and sim.read_across_shots(work).any()), f'{kind}: dirty workspace'
    assert not sim.read_phase_flipped_across_shots().any(), f'{kind}: phase left behind'
    return c.operation_statistics().num_toffoli_like


def build_dqi():
    from dqi.dqi_reference_circuit import ProosAndZalkaDecoder, default_n_iters
    b = km.CircuitBuilder()
    syndromes = [b.create_quantum_register(12, name=f'S{i}') for i in range(70)]
    record = [b.create_classical_register(12, name=f'E{j}') for j in range(4095)]
    ProosAndZalkaDecoder(70, km.GF2Field(12)).build_syndrome_decoding_mbuc(
        TowerBuilder(b), syndromes=syndromes, error_mbuc=record, n_iters=default_n_iters(70))
    return b.finish_circuit()


def load_dqi_verifier():
    """The official dqi verifier, minus its reaction_depth() report (not part of the
    score, and it can run out of memory on small machines)."""
    import types
    path = VERIFIERS / 'dqi' / 'verify_submission_dqi.py'
    mod = types.ModuleType('verify_submission_dqi_local')
    mod.__file__ = str(path)
    sys.modules[mod.__name__] = mod
    exec(compile(path.read_text().replace('circuit.reaction_depth()', '-1'), str(path), 'exec'), mod.__dict__)
    return mod


def dqi_score(V, num_qubits, mean_toffolis, failures=0, shots=102400):
    eps = V.clopper_pearson_upper_bound(failures, shots)
    return V.repetition_overhead(eps) * num_qubits**2 * mean_toffolis


def run_dqi():
    t0 = time.time()
    tw_check_formulas()
    costs = {k: tw_check_circuit(k) for k in ('inverse', 'divide', 'multiply')}
    print('unit checks passed; Toffolis per operation:', costs)
    circuit = build_dqi()
    path = SUB / 'dqi.kmb'
    circuit.to_file(str(path), format='kmb')
    result = dict(num_qubits=circuit.num_qubits, static_toffolis=circuit.operation_statistics().num_toffoli_like,
                  file_mb=round(path.stat().st_size / 1e6))
    print('built:', result)
    if VERIFY in ('quick', 'full'):
        V = load_dqi_verifier()
        full = VERIFY == 'full'
        spec = V.ProblemSpecDQI.benchmark_instance(num_shots=102400 if full else 1024)
        r = V.verify(km.Circuit.from_file(str(path)), spec, max_workers=os.cpu_count() if full else 1,
                     log2_delta=-128.0 if full else -1.0)
        failures = r['value_errors'] + r['phase_errors']
        if failures and not full:
            raise AssertionError(f'{failures} failed shots in the quick check')
        result.update(verifier='full' if full else 'quick (1,024 shots)', verifier_status=r['status'],
                      failures=failures, mean_toffolis=r['toffoli_count'])
        result['score'] = r['score'] if full else dqi_score(V, circuit.num_qubits, r['toffoli_count'])
    else:   # no simulation: estimate from the static count (the per-shot mean is ~100 lower)
        result['score'] = dqi_score(load_dqi_verifier(), circuit.num_qubits, result['static_toffolis'])
    result['seconds'] = round(time.time() - t0, 1)
    return result


def use_baseline_dqi():
    src = SAMPLE / 'dqi.kmb'
    if src.exists():
        shutil.copy(src, SUB / 'dqi.kmb')
    else:
        from dqi.dqi_reference_circuit import build_reference_circuit
        build_reference_circuit(4095, 70, 12).to_file(str(SUB / 'dqi.kmb'), format='kmb')


if USE_IMPROVED['dqi']:
    try:
        REPORT['dqi'] = dict(version='improved', **run_dqi())
    except Exception as ex:
        print('dqi: improved circuit rejected, using baseline:', ex)
        use_baseline_dqi()
        REPORT['dqi'] = dict(version='baseline (fallback)', error=str(ex))
else:
    use_baseline_dqi()
    REPORT['dqi'] = dict(version='baseline')
print(json.dumps(REPORT['dqi'], indent=2))
