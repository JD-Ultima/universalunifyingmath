"""GU round 3 (G), Task 3: principal symbols of the deformation complex (draft Sec. 10, (10.1)-(10.9)).

Formalization: formalization/geometric-unity/G-full-el-fermions.md, Sec. G.4.
Setting as in D_selfadjoint_shiab.py: d = 14, eta = (+^7, -^7) (and, for the draft's "Euclidean" claim,
eta = (+^14)); gauge algebra g = u(64,64) inside Cl(7,7) (x) C = M_128(C).  Everything below is C-linear
and maps the real form u(64,64) to itself, so ranks over C of the complexified maps on gl(128,C) = Cl (x) C
equal real ranks on u(64,64).

Exact method.  Clifford monomials E_K (K a bitmask of {0..13}) are multiplied symbolically,
E_J E_K = (-1)^{N(J,K)} prod_{i in J&K} eta_i E_{J xor K}.  Each operator below sends a basis vector
e^a (x) E_K to a signed sum of basis vectors, so its matrix is a sparse integer matrix; it splits into small
connected blocks (found by breadth-first search), whose ranks are computed by exact integer elimination
(distinct blocks are cached).  The monomial implementation is cross-checked against the explicit 128 x 128
integer matrices of clifford77.py and the matrix Shiab code of D_selfadjoint_shiab.py.

Operators (all K = Spin(7,7)-equivariant, tilted by eps in the field equations; at a point eps acts by an
invertible conjugation, which does not change ranks):
  W_A  xi = [Phi11 ^ xi]_A,  Phi11 = *Phi3,  Phi3 = sum e^abc (x) i gamma_abc       (Theorem D3; self-adjoint)
  W'_A xi = [Phi11' ^ xi]_A, Phi11' = *(sum e^abc (x) i gamma_abc Gamma)            (self-adjoint)
  RS   xi = [Phi1 ^ *xi]_C - (1/2) *[Phi1 ^ *[Phi2 ^ *xi]_A]_A   (round-1 R-S (C;A,A), c = 1; not self-adjoint;
       implemented as 2*RS so that all coefficients are integers -- ranks are unchanged)
Symbols at a covector xi (Prop. G10):
  sigma(Shiab o d)(xi) : alpha in Lambda^1 (x) g  ->  Shiab(xi ^ alpha) in Lambda^13 (x) g
  sigma(delta_2)(xi)   : (alpha, gamma') -> Shiab(xi ^ alpha) - kappa *(xi (x) gamma')    [T = 0]
  the symbol sequence 0 -> g -> Lambda^1 g (+) g -> Lambda^13 g -> 0 is exact iff rank sigma(delta_2) = 14 dim g.
Tests
  Y0  monomial algebra = explicit matrices (random products); monomial W_A, W'_A, RS = matrix versions on random
      basis inputs (exact).
  Y1  (7,7): ranks at xi = e^0 (eta = +1), e^7 (eta = -1), e^0 + e^1 (non-null, not a basis vector),
      e^0 + e^7 (null), e^0 + e^1 + e^7 + e^8 (null); kappa in {0, 1}.
  Y2  (14,0): ranks at xi = e^0 and e^0 + e^1 (every real covector is non-null).
Exact integer arithmetic throughout (Python ints).
"""
import itertools
import os
import sys
import time
from collections import deque

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from shiab_checks import Ext  # noqa: E402  (exterior algebra with arbitrary coefficients)

N = 14
FULL = (1 << N) - 1
DIMG = 1 << N                      # dim_C gl(128) = 2^14 Clifford monomials


# ------------------------------------------------------------------ Clifford monomials (bitmasks)
class Cliff:
    """Monomial algebra of Cl(eta) over Z: E_J E_K = sign * E_{J ^ K}."""

    def __init__(self, eta):
        self.eta = list(eta)
        self.neg = sum(1 << i for i, e in enumerate(eta) if e == -1)

    def sign(self, J, K):
        s = 0
        k = K
        while k:
            low = k & -k
            i = low.bit_length() - 1
            s += (J >> (i + 1)).bit_count()
            k ^= low
        s += (J & K & self.neg).bit_count()
        return -1 if s & 1 else 1


def mask(idx):
    m = 0
    for i in idx:
        m |= 1 << i
    return m


class CL:
    """Clifford element: dict mask -> coefficient (exact Gaussian integers as Python complex)."""
    __slots__ = ('d',)

    def __init__(self, d=None):
        self.d = d or {}

    def __add__(self, o):
        if isinstance(o, int) and o == 0:
            return self
        out = dict(self.d)
        for k, v in o.d.items():
            out[k] = out.get(k, 0) + v
        return CL({k: v for k, v in out.items() if v != 0})
    __radd__ = __add__

    def __rmul__(self, c):
        return CL({k: c * v for k, v in self.d.items()})


class OP:
    """Linear map X -> sum c E_L X E_R on Cl, stored as dict (L, R) -> c."""
    __slots__ = ('d',)

    def __init__(self, d=None):
        self.d = d or {}

    def __add__(self, o):
        if isinstance(o, int) and o == 0:
            return self
        out = dict(self.d)
        for k, v in o.d.items():
            out[k] = out.get(k, 0) + v
        return OP({k: v for k, v in out.items() if v != 0})
    __radd__ = __add__

    def __rmul__(self, c):
        return OP({k: c * v for k, v in self.d.items()})

    def __neg__(self):
        return (-1) * self

    def __sub__(self, o):
        return self + (-1) * o


IDOP = OP({(0, 0): 1})


def lmul(cl, P, Y):
    out = {}
    for m, cm in P.d.items():
        for (L, R), c in Y.d.items():
            k = (m ^ L, R)
            out[k] = out.get(k, 0) + cm * c * cl.sign(m, L)
    return OP({k: v for k, v in out.items() if v != 0})


def rmul(cl, Y, P):
    out = {}
    for m, cm in P.d.items():
        for (L, R), c in Y.d.items():
            k = (L, R ^ m)
            out[k] = out.get(k, 0) + cm * c * cl.sign(R, m)
    return OP({k: v for k, v in out.items() if v != 0})


def br(cl, s):
    """bracket [P ^ Y]_s with P a Clifford coefficient (form factor 1st), Y an OP coefficient."""
    if s == 'C':
        return lambda P, Y: lmul(cl, P, Y) - rmul(cl, Y, P)
    return lambda P, Y: 1j * (lmul(cl, P, Y) + rmul(cl, Y, P))


def apply_op(cl, op, K):
    """op applied to E_K: dict mask -> coefficient."""
    out = {}
    for (L, R), c in op.d.items():
        LK = L ^ K
        s = cl.sign(L, K) * cl.sign(LK, R)
        k = LK ^ R
        out[k] = out.get(k, 0) + c * s
    return {k: v for k, v in out.items() if v != 0}


# ------------------------------------------------------------------ the Shiab operators on e^b ^ e^a (x) X
def shiab_ops(eta, which):
    """dict (b, a) -> dict m -> OP, where Shiab(e^b ^ e^a (x) X) = sum_m [13-form missing m] (x) OP_m(X)."""
    cl = Cliff(eta)
    E = Ext(eta)
    G = FULL
    if which in ('WA', 'WpA'):
        Phi3 = {abc: CL({mask(abc) ^ (G if which == 'WpA' else 0):
                         1j * (1 if which == 'WA' else cl.sign(mask(abc), G))}) for abc in itertools.combinations(range(N), 3)}
        Phi11 = E.star(Phi3)
        A = br(cl, 'A')

        def S(xi):
            return E.wedge(Phi11, xi, A)
    elif which == 'RS':
        Phi1 = {(m,): CL({1 << m: 1}) for m in range(N)}
        Phi2 = {(m, k): CL({(1 << m) | (1 << k): 1}) for m, k in itertools.combinations(range(N), 2)}
        C_, A_ = br(cl, 'C'), br(cl, 'A')

        def S(xi):
            t1 = E.wedge(Phi1, E.star(xi), C_)
            inner = E.star(E.wedge(Phi2, E.star(xi), A_))
            t2 = E.star(E.wedge(Phi1, inner, A_))
            out = {K: 2 * v for K, v in t1.items()}          # 2 * RS (integer coefficients; rank unchanged)
            for K, v in t2.items():
                out[K] = out.get(K, 0) + (-1) * v
            return out
    else:
        raise ValueError(which)
    ops = {}
    for b in range(N):
        for a in range(N):
            if a == b:
                continue
            s = 1 if b < a else -1
            xi = {tuple(sorted((a, b))): s * IDOP}
            res = S(xi)
            d = {}
            for I, op in res.items():
                assert len(I) == N - 1
                m = next(i for i in range(N) if i not in I)
                op = OP({k: v for k, v in op.d.items() if v != 0})
                if op.d:
                    d[m] = op
            ops[(b, a)] = d
    return cl, E, ops


def combine_ops(o1, o2, x, y):
    """ops of x*S1 + y*S2"""
    out = {}
    for key in set(o1) | set(o2):
        d = {}
        for o, w in ((o1, x), (o2, y)):
            for m, op in o.get(key, {}).items():
                d[m] = d.get(m, 0) + w * op
        out[key] = {m: OP({k: v for k, v in op.d.items() if v != 0}) for m, op in d.items()}
    return out


def exact_coeff(z):
    z = complex(z)
    assert z.imag == 0 and z.real == int(z.real) and abs(z.real) < 2 ** 50, z
    return int(z.real)


# ------------------------------------------------------------------ exact rank of sparse integer matrices
def int_rank(rows):
    """exact rank over Q of an integer matrix (list of lists), integer elimination with gcd reduction."""
    from math import gcd
    M = [list(r) for r in rows if any(r)]
    if not M:
        return 0
    ncol = len(M[0])
    r = 0
    for c in range(ncol):
        piv = None
        for i in range(r, len(M)):
            if M[i][c] != 0:
                piv = i
                break
        if piv is None:
            continue
        M[r], M[piv] = M[piv], M[r]
        p = M[r][c]
        for i in range(r + 1, len(M)):
            q = M[i][c]
            if q:
                row = [p * x - q * y for x, y in zip(M[i], M[r])]
                g = 0
                for x in row:
                    g = gcd(g, x)
                M[i] = [x // g for x in row] if g > 1 else row
        r += 1
        if r == len(M):
            break
    return r


def sparse_rank(rows, cols, vals, nr, nc):
    """exact rank of the integer matrix sum vals[i] e_{rows[i], cols[i]} (duplicates summed): split into
    connected blocks (scipy.sparse.csgraph), exact integer elimination per distinct block."""
    import scipy.sparse as sps
    from scipy.sparse.csgraph import connected_components
    A = sps.coo_matrix((np.asarray(vals, dtype=np.int64), (np.asarray(rows), np.asarray(cols))),
                       shape=(nr, nc)).tocsr()
    A.sum_duplicates()
    A.eliminate_zeros()
    assert np.abs(A.data).max(initial=0) < 2 ** 40
    Ac = A.tocoo()
    G = sps.coo_matrix((np.ones(Ac.nnz), (Ac.row, Ac.col + nr)), shape=(nr + nc, nr + nc))
    ncomp, lab = connected_components(G, directed=False)
    order = np.argsort(lab[Ac.row], kind='stable')
    r_, c_, v_ = Ac.row[order], Ac.col[order], Ac.data[order]
    comp = lab[r_]
    bounds = np.flatnonzero(np.diff(comp)) + 1
    cache, total, nb, mx = {}, 0, 0, 0
    for seg in np.split(np.arange(len(comp)), bounds):
        if len(seg) == 0:
            continue
        rr, cc, vv = r_[seg], c_[seg], v_[seg]
        ur, ri = np.unique(rr, return_inverse=True)
        uc, ci = np.unique(cc, return_inverse=True)
        M = np.zeros((len(ur), len(uc)), dtype=np.int64)
        M[ri, ci] = vv
        key = (M.shape, M.tobytes())
        if key not in cache:
            cache[key] = int_rank(M.tolist())
        total += cache[key]
        nb += 1
        mx = max(mx, M.shape[0], M.shape[1])
    return total, nb, mx, len(cache)


# ------------------------------------------------------------------ symbol matrices (vectorized signs)
_KS = np.arange(DIMG, dtype=np.int64)


def parity_table():
    P = np.zeros(DIMG, dtype=np.int64)
    for k in range(N):
        P |= (np.bitwise_count(_KS >> (k + 1)).astype(np.int64) & 1) << k
    return P


_PJ = parity_table()


def vsign(cl, J, K):
    """vectorised sign(J, K): J scalar or array, K array."""
    s = (np.bitwise_count(_PJ[J] & K).astype(np.int64) + np.bitwise_count(J & K & cl.neg).astype(np.int64)) & 1
    return 1 - 2 * s


def symbol_rank(ops, cl, E, xi, kappa, colmask=None):
    rows, cols, vals, nr, nc = symbol_entries(ops, cl, E, xi, kappa)
    if colmask is not None:
        keep = colmask[cols]
        rows, cols, vals = rows[keep], cols[keep], vals[keep]
    return sparse_rank(rows, cols, vals, nr, nc)


def symbol_entries(ops, cl, E, xi, kappa):
    """rank of (alpha, gamma') -> Shiab(xi ^ alpha) - kappa *(xi (x) gamma') on (Lambda^1 (x) Cl) (+) Cl.
    Columns a*2^14 + K (alpha = e^a (x) E_K) and 14*2^14 + K (gamma' = E_K);
    rows m*2^14 + K ([13-form missing m] (x) E_K)."""
    rows, cols, vals = [], [], []
    for (b, a), d in ops.items():
        if xi[b] == 0:
            continue
        for m, op in d.items():
            for (L, R), c in op.d.items():
                c = exact_coeff(c) * xi[b]
                if c == 0:
                    continue
                LK = L ^ _KS
                sg = vsign(cl, L, _KS) * vsign(cl, LK, R)
                rows.append(m * DIMG + (LK ^ R))
                cols.append(a * DIMG + _KS)
                vals.append(c * sg)
    nc = N * DIMG
    if kappa:
        for b in range(N):
            if xi[b]:
                cb = E.star_coeff((b,))[1]
                rows.append(b * DIMG + _KS)
                cols.append(N * DIMG + _KS)
                vals.append(np.full(DIMG, -kappa * xi[b] * cb, dtype=np.int64))
        nc += DIMG
    return np.concatenate(rows), np.concatenate(cols), np.concatenate(vals), N * DIMG, nc


# ------------------------------------------------------------------ cross-checks against explicit matrices
def crosscheck(rng):
    from clifford77 import clifford_pq
    from shiab_checks import Shiab, comm, acomm
    pos, neg = clifford_pq(7)
    gam = [g.astype(complex) for g in pos + neg]
    eta = [1] * 7 + [-1] * 7
    cl = Cliff(eta)
    Id = np.eye(128, dtype=complex)

    def mono(K):
        M = Id.copy()
        for i in range(N):
            if K >> i & 1:
                M = M @ gam[i]
        return M
    ok = True
    for _ in range(200):
        J, K = int(rng.integers(0, DIMG)), int(rng.integers(0, DIMG))
        ok &= np.array_equal(mono(J) @ mono(K), cl.sign(J, K) * mono(J ^ K))
    print(f"Y0. monomial products E_J E_K = sign * E_(J^K) agree with clifford77 matrices (200 random pairs): {ok}")
    assert ok
    E = Ext(eta)
    Gm = mono(FULL)
    Phi3 = {abc: 1j * mono(mask(abc)) for abc in itertools.combinations(range(N), 3)}
    Phi3p = {abc: 1j * mono(mask(abc)) @ Gm for abc in itertools.combinations(range(N), 3)}
    mats = {'WA': lambda xi: E.wedge(E.star(Phi3), xi, acomm),
            'WpA': lambda xi: E.wedge(E.star(Phi3p), xi, acomm),
            'RS': lambda xi: {K: 2 * v for K, v in Shiab(E, gam).apply(xi, 'C', 'A', 'A', 1.0).items()}}
    for which in ('WA', 'WpA', 'RS'):
        _, _, ops = shiab_ops(eta, which)
        good = True
        for _ in range(4):
            b, a = rng.choice(N, 2, replace=False)
            b, a = int(b), int(a)
            K = int(rng.integers(0, DIMG))
            s = 1 if b < a else -1
            ref = mats[which]({tuple(sorted((a, b))): s * mono(K)})
            got = {}
            for m, op in ops[(b, a)].items():
                I = tuple(i for i in range(N) if i != m)
                got[I] = sum((c * mono(k) for k, c in apply_op(cl, op, K).items()), 0 * Id)
            for I in set(ref) | set(got):
                good &= np.array_equal(np.asarray(ref.get(I, 0 * Id)), np.asarray(got.get(I, 0 * Id)))
        print(f"    monomial {which:3s} = matrix {which:3s} on 4 random inputs e^b ^ e^a (x) E_K (exact): {good}")
        assert good


def validate_rank_pipeline(rng):
    """Y0b: (i) int_rank = sympy rank; (ii) assembled symbol columns = apply_op; (iii) independent rank of
    W_A(xi ^ .) on Lambda^1 (x) span{E_K : |K| <= 1} from the explicit 128 x 128 matrices (exact integer images,
    exact rank of a random integer projection, which can only lower the rank) versus the monomial pipeline."""
    import sympy as sp
    from clifford77 import clifford_pq
    from shiab_checks import acomm
    ok = True
    for _ in range(40):
        r, c, k = (int(x) for x in rng.integers(1, 9, 3))
        Mx = rng.integers(-3, 4, (r, k)) @ rng.integers(-3, 4, (k, c))
        ok &= int_rank(Mx.tolist()) == sp.Matrix(Mx.tolist()).rank()
    print(f"Y0b. (i) int_rank = sympy rank on 40 random low-rank integer matrices: {ok}")
    assert ok
    eta = [1] * 7 + [-1] * 7
    cl, E, ops = shiab_ops(eta, 'WA')
    xi = [0] * N
    xi[0] = xi[7] = 1
    rows, cols, vals, nr, nc = symbol_entries(ops, cl, E, xi, 1)
    import scipy.sparse as sps
    A = sps.coo_matrix((vals, (rows, cols)), shape=(nr, nc)).tocsc()
    good = True
    for _ in range(30):
        a, K = int(rng.integers(0, N)), int(rng.integers(0, DIMG))
        ref = {}
        for b in range(N):
            if xi[b] and (b, a) in ops:
                for m, op in ops[(b, a)].items():
                    for k2, v in apply_op(cl, op, K).items():
                        ref[m * DIMG + k2] = ref.get(m * DIMG + k2, 0) + xi[b] * exact_coeff(v)
        col = A[:, a * DIMG + K].toarray().ravel()
        got = {int(i): int(col[i]) for i in np.flatnonzero(col)}
        good &= got == {k_: v for k_, v in ref.items() if v}
    print(f"     (ii) assembled symbol columns = direct application (30 random columns, xi = e^0+e^7): {good}")
    assert good
    pos, neg = clifford_pq(7)
    gam = [g.astype(np.int64) for g in pos + neg]
    Id = np.eye(128, dtype=np.int64)

    def mono(K):
        M = Id.copy()
        for i in range(N):
            if K >> i & 1:
                M = M @ gam[i]
        return M
    Phi11 = E.star({abc: mono(mask(abc)) for abc in itertools.combinations(range(N), 3)})
    acomm_real = lambda P_, Y: -(P_ @ Y + Y @ P_)       # [i P ^ Y]_A with P real: i * i (PY + YP)
    Ks = [0] + [1 << i for i in range(N)]
    U = rng.integers(-5, 6, (12, 128))
    V = rng.integers(-5, 6, (128, 12))
    for xd in ({0: 1}, {0: 1, 7: 1}):
        xiv = [xd.get(i, 0) for i in range(N)]
        colsimg = []
        for a in range(N):
            for K in Ks:
                al = {(a,): mono(K)}
                xa = E.wedge({(b,): xiv[b] for b in range(N) if xiv[b]}, al, lambda s, v: s * v)
                img = E.wedge(Phi11, xa, acomm_real)
                vec = []
                for m in range(N):
                    I = tuple(i for i in range(N) if i != m)
                    Mm = img.get(I, 0 * Id)
                    vec.extend((U @ Mm @ V).ravel().tolist())
                colsimg.append(vec)
        r_proj = int_rank(colsimg)
        cm = np.zeros(N * DIMG, dtype=bool)
        for a in range(N):
            for K in Ks:
                cm[a * DIMG + K] = True
        r_mono, *_ = symbol_rank(ops, cl, E, xiv, 0, colmask=np.concatenate([cm, np.zeros(DIMG, dtype=bool)]))
        print(f"     (iii) xi = {xd}: rank of W_A(xi ^ .) on Lambda^1 (x) span{{E_K: |K|<=1}} (210 columns): "
              f"matrix route (random projection, lower bound) {r_proj}, monomial route {r_mono}: equal {r_proj == r_mono}")
        assert r_proj == r_mono


def main():
    t0 = time.time()
    rng = np.random.default_rng(20210401)
    crosscheck(rng)
    validate_rank_pipeline(rng)
    full = N * DIMG
    print(f"\n dim g = {DIMG}; Lambda^1 (x) g = Lambda^13 (x) g has dimension {full}; "
          f"exactness of the symbol sequence <=> rank sigma(delta_2) = {full}\n")
    cases77 = [("e^0      (eta = +1)", {0: 1}), ("e^7      (eta = -1)", {7: 1}),
               ("e^0+e^1  (eta = +2)", {0: 1, 1: 1}), ("e^0+e^7  (null)", {0: 1, 7: 1}),
               ("e^0+e^1+e^7+e^8 (null)", {0: 1, 1: 1, 7: 1, 8: 1})]
    print("=== Y1. signature (7,7) ===")
    eta77 = [1] * 7 + [-1] * 7
    base = {w: shiab_ops(eta77, w) for w in ('WA', 'WpA', 'RS')}
    cl, E, _ = base['WA']
    family = [('WA', base['WA'][2]), ('WpA', base['WpA'][2]), ('RS', base['RS'][2]),
              ('WA+WpA', combine_ops(base['WA'][2], base['WpA'][2], 1, 1)),
              ('WA-WpA', combine_ops(base['WA'][2], base['WpA'][2], 1, -1)),
              ('WA+2WpA', combine_ops(base['WA'][2], base['WpA'][2], 1, 2))]
    for which, ops in family:
        for lab, xd in cases77:
            xi = [xd.get(i, 0) for i in range(N)]
            r0, *_ = symbol_rank(ops, cl, E, xi, 0)
            r1, nb, mx, nd = symbol_rank(ops, cl, E, xi, 1)
            print(f"  {which:7s} xi = {lab:24s}: rank Shiab(xi^.) = {r0:6d} (= {full} - {full - r0}); "
                  f"rank sigma(delta_2), kappa=1: {r1:6d} (deficit {full - r1}); exact: {r1 == full}  "
                  f"[{nb} blocks, max size {mx}, {nd} distinct]")
    print("\n=== Y2. signature (14,0) (draft's 'Euclidean' claim, GU-057/GU-072) ===")
    b14 = {w: shiab_ops([1] * 14, w) for w in ('WA', 'WpA', 'RS')}
    cl, E, _ = b14['WA']
    fam14 = [('WA', b14['WA'][2]), ('WpA', b14['WpA'][2]), ('RS', b14['RS'][2]),
             ('WA+WpA', combine_ops(b14['WA'][2], b14['WpA'][2], 1, 1)),
             ('WA+2WpA', combine_ops(b14['WA'][2], b14['WpA'][2], 1, 2))]
    for which, ops in fam14:
        for lab, xd in [("e^0", {0: 1}), ("e^0+e^1", {0: 1, 1: 1})]:
            xi = [xd.get(i, 0) for i in range(N)]
            r0, *_ = symbol_rank(ops, cl, E, xi, 0)
            r1, *_ = symbol_rank(ops, cl, E, xi, 1)
            print(f"  {which:7s} xi = {lab:8s}: rank Shiab(xi^.) = {r0:6d} (= {full} - {full - r0}); "
                  f"rank sigma(delta_2), kappa=1: {r1:6d} (deficit {full - r1}); exact (kappa=1): {r1 == full}")
    print("\n=== Y3. gauge algebra so(64,64) (grades 1,2,5,6,9,10,13,14; preserved by W_A, W'_A) ===")
    so_grades = {1, 2, 5, 6, 9, 10, 13, 14}
    in_so = np.array([bin(K).count('1') in so_grades for K in range(DIMG)])
    dso = int(in_so.sum())
    fso = N * dso
    colmask = np.concatenate([np.tile(in_so, N), in_so])
    for eta, sig, fam_, cl_, E_ in ((eta77, '(7,7)', family, base['WA'][0], base['WA'][1]),
                                    ([1] * 14, '(14,0)', fam14, b14['WA'][0], b14['WA'][1])):
        for which, ops in fam_:
            if which not in ('WA', 'WpA', 'WA+WpA'):
                continue
            for lab, xd in [("e^0", {0: 1}), ("e^0+e^7", {0: 1, 7: 1})]:
                xi = [xd.get(i, 0) for i in range(N)]
                r0, *_ = symbol_rank(ops, cl_, E_, xi, 0, colmask=colmask)
                r1, *_ = symbol_rank(ops, cl_, E_, xi, 1, colmask=colmask)
                print(f"  {sig:6s} {which:7s} xi = {lab:8s}: dim so(64,64) = {dso}; rank Shiab(xi^.) = {r0:6d} "
                      f"(= {fso} - {fso - r0}); rank sigma(delta_2), kappa=1: {r1:6d} (deficit {fso - r1})")
    print(f"\nruntime {time.time() - t0:.0f} s")


if __name__ == '__main__':
    main()
