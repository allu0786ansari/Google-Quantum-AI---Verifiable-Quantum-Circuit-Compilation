"""Layers I-2/I-3 (polynomial route): y = c0 - z*(c1 - z*(c2 - z*c3)) by truncated Horner rows.

x holds x_norm (after I-1). idx = x[15:22] addresses one lookup that loads all four coefficients
(units 2^-29) into the registers that become the Horner accumulators:
    A2 <- c2 - z*A3,  A1 <- c1 - z*A2,  Y <- c0 - z*A1
Row j of a stage with nz z-bits is `target -= A >> (nz - j)` controlled by z bit j, which is x_norm
bit 15 - nz + j. This is exactly the model in poly_model.fit (band_mul with cut 2^-29).
Y lies in [2^28, 2^29) and every partial result stays there, so its top bit is a constant 1 that is
never stored (Y register = 28 qubits, value Y - 2^28).
"""
import numpy as np
import kickmix as km

T = np.load('/home/claude/google-quantum-ai---verifiable-quantum-circuit-compilation/ours/invsqrt/poly_tables.npz')
C = [[int(v) for v in T[f'C{k}']] for k in range(4)]
NZ = [int(v) for v in T['NZ']]
W = [29, 23, 18, 13]
WY = 28                                     # stored bits of Y (bit 28 is the constant 1)
ADDR = list(range(15, 22))


T32 = [c3 | (c2 << W[3]) for c2, c3 in zip(C[2], C[3])]
T10 = [c1 | ((c0 - (1 << 28) if c0 else 0) << W[1]) for c0, c1 in zip(C[0], C[1])]


def _addr(x):
    return km.array([x[i] for i in ADDR])


MODE = {'rows': 'and'}        # 'ctl': controlled adder per row; 'and': AND the row into scratch, plain adder


def _stage(b, x, target, src, nz, sign):
    """target -= sum_j z_j * (src >> (nz - j))   (sign=-1), or the inverse (sign=+1)."""
    src = list(src)
    rows = range(nz) if sign < 0 else reversed(range(nz))
    op = b.isub if sign < 0 else b.iadd
    for j in rows:
        lo = nz - j
        if lo >= len(src):
            continue
        ctrl = x[15 - nz + j]
        if MODE['rows'] == 'ctl':
            op(km.array(src[lo:]), target=km.array(list(target)), control=ctrl)
        else:
            row = src[lo:]
            tmp = b.alloc_qubits(len(row))
            for a, t in zip(row, tmp):
                b.init_and(ctrl, a, t)
            op(km.array(list(tmp)), target=km.array(list(target)))
            for a, t in zip(row, tmp):
                b.del_and(ctrl, a, t)
            b.free(tmp)


def horner(b, x):
    a3, a2 = b.alloc_qubits(W[3]), b.alloc_qubits(W[2])
    b.init_lookup(table=T32, address=_addr(x), target=km.array(list(a3) + list(a2)))
    _stage(b, x, a2, a3, NZ[3], -1)
    b.del_lookup(table=C[3], address=_addr(x), target=a3); b.free(a3)
    a1, y = b.alloc_qubits(W[1]), b.alloc_qubits(WY)
    b.init_lookup(table=T10, address=_addr(x), target=km.array(list(a1) + list(y)))
    _stage(b, x, a1, a2, NZ[2], -1)
    _stage(b, x, y, a1, NZ[1], -1)
    return y, a1, a2


def unhorner(b, x, y, a1, a2):
    _stage(b, x, y, a1, NZ[1], +1)
    _stage(b, x, a1, a2, NZ[2], +1)
    b.del_lookup(table=T10, address=_addr(x), target=km.array(list(a1) + list(y))); b.free(a1); b.free(y)
    a3 = b.alloc_qubits(W[3])
    b.init_lookup(table=C[3], address=_addr(x), target=a3)
    _stage(b, x, a2, a3, NZ[3], +1)
    b.del_lookup(table=T32, address=_addr(x), target=km.array(list(a3) + list(a2))); b.free(a3); b.free(a2)
