"""Referee cross-check, round CT-5 (reports/reviews/REVIEW5-CT-2026-10-30.md).

Independent of conjugation.py and records_uniform.py: neither is imported or copied.  The only import from this
directory is config_enumeration.py (unmodified; its printout is captured and compared with config_enumeration.out),
used for Theorem 1's predicates (readable configurations, end rule, two-point definition).

Everything is exact: Gaussian rationals (own class over fractions.Fraction) for Part A, rationals for Part B (all
letters are real rational matrices: Lean's tracelessBasis, i.e. the Gell-Mann / Pauli letters up to factors of i
and of sqrt 3, which change no span, no orthogonality and no (anti)commutation), sympy only for the genuine
Gell-Mann matrices and for the final 1/value numbers.  Real spans of the real Lie algebra u(d1)+u(d2) are computed
on realified vectors.  No randomness.

Part A  conjugation.md (CT-G03/G04)
  A1  relabellings eps: image of g on L, the scalar of Q on L, dependence on e1 e2 (Thm 3); the image of h_W itself
  A2  stabilizer of the pair ket; centres at (2,3) under A and A* in labels and as operator pairs (Thm 4(1), Prop. 6);
      the Supplement's explicit null-space basis (Sec. S-matrices); the 128 of 12288 finite matchings
  A3  matching closure for several [alpha:beta] and in eps-labels (Thm 4(2),(3))
  A4  Lie closure of su+su with m_0 and with the full stabilizer (Thm 5); i 1_L test (Thm 5(4))
  A5  Theorem 1's selection per reading through config_enumeration's predicates (Thm 8, Cor. 9)
  A6  trace calculus spot checks against explicit representation matrices; hypercharges (Thm 8(1), Prop. 1)
  A7  the clause P relative to the Letter's own definition; Prop. 6(4) in labels and in operator form
Part B  records-uniform.md (CT-G09) and the F-CT-7 correction
  B1  model: letters, structure constants, P<=2
  B2  order-sign domain D (Thm 2), one-end anticommutators at d = 2, 3, 4, and a counterexample to Thm 2(3)'s "iff"
  B3  C_L, C, C+ and Lie(G_Omega) (Prop. 1)
  B4  frames: classes and values under the Letter's weights; P2 + C+ under E'; Gell-Mann literal (Prop. 1(4), Sec. 8)
  B5  exchanges: G_Omega covariance on P<=2, involution, algebra type, Schur intertwiners, Fm = Ad_W (Thm 7)
  B6  identification rule: stab_Omega, Comp_f, D cap Comp_f^perp (Thm 6)
  B7  anticommuting letters: Gell-Mann, Lean's basis, Weyl at d = 3 and d = 4 (Prop. 3; F-CT-7)
  B8  collapse kappa(X) = X|Omega> (Prop. 8)
  B9  F-CT-7: class counts by enumeration at several (dQ, dq); the (2,4) values 331.40 / 307.91 / 216.04; (OD)+(ID)
      at (dQ, dq) = (2,3) and the (2,3) value with commutator records at its qubit centre

Run from anywhere:  python3 -I xcheck_round_ct5.py > xcheck_round_ct5.out
"""
import contextlib
import io
import itertools
import os
import sys
import time
from fractions import Fraction as Fr

T0 = time.perf_counter()
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import sympy as sp

NCHK = [0, 0]


def check(cond, msg):
    cond = bool(cond)
    NCHK[0] += 1
    NCHK[1] += cond
    print(f"  [{NCHK[0]:3d}] {msg}: {cond}")
    assert cond, msg


def head(s):
    print()
    print(s)
    print("-" * len(s))


# =====================================================================================================
# exact sparse linear algebra over Q (vectors are dicts key -> Fraction)
# =====================================================================================================
class Span:
    """incremental row echelon basis over Q"""

    def __init__(self, vecs=()):
        self.piv, self.rows = [], []
        for v in vecs:
            self.add(v)

    def reduce(self, v):
        w = {k: Fr(x) for k, x in v.items() if x != 0}
        for p, b in zip(self.piv, self.rows):
            c = w.get(p)
            if c:
                for k, x in b.items():
                    y = w.get(k, 0) - c * x
                    if y:
                        w[k] = y
                    else:
                        w.pop(k, None)
        return w

    def add(self, v):
        w = self.reduce(v)
        if not w:
            return False
        p = min(w)
        c = w[p]
        self.piv.append(p)
        self.rows.append({k: x / c for k, x in w.items()})
        return True

    def contains(self, v):
        return not self.reduce(v)

    @property
    def dim(self):
        return len(self.rows)


def rank(vecs):
    return Span(vecs).dim


def same_space(A, B):
    sa = Span(A)
    return all(sa.contains(b) for b in B) and Span(B).dim == sa.dim


def inter_dim(A, B):
    return rank(A) + rank(B) - rank(list(A) + list(B))


def nullspace(vecs):
    """all coefficient vectors c with sum_i c_i vecs[i] = 0 (as dicts i -> c_i)"""
    piv, rows, out = [], [], []
    for i, v in enumerate(vecs):
        w = {("v",) + (k if isinstance(k, tuple) else (k,)): Fr(x) for k, x in v.items() if x != 0}
        w[("z", i)] = Fr(1)
        for p, b in zip(piv, rows):
            c = w.get(p)
            if c:
                for k, x in b.items():
                    y = w.get(k, 0) - c * x
                    if y:
                        w[k] = y
                    else:
                        w.pop(k, None)
        vk = [k for k in w if k[0] == "v"]
        if not vk:
            out.append({k[1]: x for k, x in w.items()})
            continue
        p = min(vk)
        c = w[p]
        piv.append(p)
        rows.append({k: x / c for k, x in w.items()})
    return out


def combo(coeffs, vecs):
    out = {}
    for i, c in coeffs.items():
        for k, x in vecs[i].items():
            y = out.get(k, 0) + c * x
            if y:
                out[k] = y
            else:
                out.pop(k, None)
    return out


# =====================================================================================================
# Gaussian rationals and small complex matrices (Part A)
# =====================================================================================================
class G:
    __slots__ = ("re", "im")

    def __init__(self, re=0, im=0):
        self.re, self.im = Fr(re), Fr(im)

    def __add__(self, o):
        return G(self.re + o.re, self.im + o.im)

    def __sub__(self, o):
        return G(self.re - o.re, self.im - o.im)

    def __mul__(self, o):
        return G(self.re * o.re - self.im * o.im, self.re * o.im + self.im * o.re)

    def __neg__(self):
        return G(-self.re, -self.im)

    def conj(self):
        return G(self.re, -self.im)

    def iszero(self):
        return self.re == 0 and self.im == 0

    def __eq__(self, o):
        return self.re == o.re and self.im == o.im


ZERO, ONE, IU = G(0), G(1), G(0, 1)


def zm(r, c=None):
    return [[G() for _ in range(c if c is not None else r)] for _ in range(r)]


def eye(n):
    m = zm(n)
    for i in range(n):
        m[i][i] = G(1)
    return m


def Emat(n, a, b, z=ONE):
    m = zm(n)
    m[a][b] = z
    return m


def madd(A, B):
    return [[a + b for a, b in zip(ra, rb)] for ra, rb in zip(A, B)]


def msub(A, B):
    return [[a - b for a, b in zip(ra, rb)] for ra, rb in zip(A, B)]


def msc(z, A):
    return [[z * a for a in r] for r in A]


def mmul(A, B):
    n, k, m = len(A), len(B), len(B[0])
    out = zm(n, m)
    for i in range(n):
        for l in range(k):
            a = A[i][l]
            if a.iszero():
                continue
            for j in range(m):
                b = B[l][j]
                if not b.iszero():
                    out[i][j] = out[i][j] + a * b
    return out


def mconj(A):
    return [[a.conj() for a in r] for r in A]


def mT(A):
    return [list(c) for c in zip(*A)]


def kron(A, B):
    n, m = len(A), len(B)
    out = zm(n * m)
    for i in range(n):
        for j in range(n):
            a = A[i][j]
            if a.iszero():
                continue
            for k in range(m):
                for l in range(m):
                    out[i * m + k][j * m + l] = a * B[k][l]
    return out


def mtr(A):
    t = G()
    for i in range(len(A)):
        t = t + A[i][i]
    return t


def mvec(A, v):
    return [sum((A[i][j] * v[j] for j in range(len(v))), G()) for i in range(len(A))]


def realify_m(A, tag=()):
    out = {}
    for i, r in enumerate(A):
        for j, a in enumerate(r):
            if a.re:
                out[tag + (i, j, 0)] = a.re
            if a.im:
                out[tag + (i, j, 1)] = a.im
    return out


def realify_v(v, tag=()):
    out = {}
    for i, a in enumerate(v):
        if a.re:
            out[tag + (i, 0)] = a.re
        if a.im:
            out[tag + (i, 1)] = a.im
    return out


def diag_of(A):
    return [A[i][i] for i in range(len(A))]


def gstr(z):
    def f(x):
        return str(x)
    if z.im == 0:
        return f(z.re)
    if z.re == 0:
        return f"{f(z.im)}i"
    return f"{f(z.re)}{'+' if z.im > 0 else '-'}{f(abs(z.im))}i"


# ---------------------------------------------------------------- the real Lie algebra g = u(d1) + u(d2)
def u_basis(d):
    B = []
    for j in range(d):
        for k in range(j + 1, d):
            B.append(msub(Emat(d, j, k), Emat(d, k, j)))
            B.append(madd(Emat(d, j, k, IU), Emat(d, k, j, IU)))
    for j in range(d):
        B.append(Emat(d, j, j, IU))
    return B


def su_basis(d):
    B = [b for b in u_basis(d)[: d * (d - 1)]]
    for j in range(d - 1):
        B.append(madd(Emat(d, j, j, IU), Emat(d, j + 1, j + 1, -IU)))
    return B


def g_basis(d1, d2):
    return [(X, zm(d2)) for X in u_basis(d1)] + [(zm(d1), X) for X in u_basis(d2)]


def su_su(d1, d2):
    return [(X, zm(d2)) for X in su_basis(d1)] + [(zm(d1), X) for X in su_basis(d2)]


def gadd(x, y):
    return (madd(x[0], y[0]), madd(x[1], y[1]))


def gsc(c, x):
    return (msc(G(c), x[0]), msc(G(c), x[1]))


def gbr(x, y):
    return (msub(mmul(x[0], y[0]), mmul(y[0], x[0])), msub(mmul(x[1], y[1]), mmul(y[1], x[1])))


def greal(x):
    v = realify_m(x[0], (0,))
    v.update(realify_m(x[1], (1,)))
    return v


def gcombo(coeffs, basis):
    d1, d2 = len(basis[0][0]), len(basis[0][1])
    out = (zm(d1), zm(d2))
    for i, c in coeffs.items():
        out = gadd(out, gsc(c, basis[i]))
    return out


def trs(x):
    """(tr X1 / i, tr X2 / i), real for anti-Hermitian blocks"""
    return mtr(x[0]).im, mtr(x[1]).im


def drho(x, eps=(1, 1)):
    X1 = x[0] if eps[0] == 1 else mconj(x[0])
    X2 = x[1] if eps[1] == 1 else mconj(x[1])
    return madd(kron(X1, eye(len(X2))), kron(eye(len(X1)), X2))


def omega(d1, d2):
    m = min(d1, d2)
    v = [G() for _ in range(d1 * d2)]
    for i in range(m):
        v[i * d2 + i] = G(1)
    return v


def kernel_of_functional(basis, phi):
    vals = [phi(x) for x in basis]
    nz = [j for j, v in enumerate(vals) if v != 0]
    if not nz:
        return list(basis)
    j0 = nz[0]
    return [gadd(basis[j], gsc(-vals[j] / vals[j0], basis[j0])) for j in range(len(basis)) if j != j0]


def subspace_from_null(basis, vec_of):
    ns = nullspace([vec_of(x) for x in basis])
    return [gcombo(c, basis) for c in ns]


def lie_closure(gens):
    S = Span()
    B = []
    for x in gens:
        if S.add(greal(x)):
            B.append(x)
    while True:
        new = []
        for x, y in itertools.combinations(B, 2):
            z = gbr(x, y)
            if S.add(greal(z)):
                new.append(z)
        if not new:
            return B
        B += new


def centre(B):
    """centre of the Lie algebra spanned by the list B (assumed closed)"""
    vecs = []
    for x in B:
        v = {}
        for l, y in enumerate(B):
            for k, val in greal(gbr(x, y)).items():
                v[(l,) + k] = val
        vecs.append(v)
    return [gcombo(c, B) for c in nullspace(vecs)]


def in_real_span(B, x):
    return Span([greal(b) for b in B]).contains(greal(x))


def acting(h, d1, d2, eps=(1, 1)):
    """is i 1_L in d rho_eps(h)?"""
    N = d1 * d2
    S = Span([realify_m(drho(x, eps)) for x in h])
    return S.contains(realify_m(msc(IU, eye(N))))


def h_ab(d1, d2, a, b):
    return kernel_of_functional(g_basis(d1, d2), lambda x: a * trs(x)[0] + b * trs(x)[1])


def h_B(d1, d2, eps=(1, 1)):
    return kernel_of_functional(g_basis(d1, d2), lambda x: mtr(drho(x, eps)).im)


def stab(h, d1, d2, eps=(1, 1)):
    w = omega(d1, d2)
    return subspace_from_null(h, lambda x: realify_v(mvec(drho(x, eps), w)))


def m0(d1, d2):
    """matchings (zeta, conj zeta + 0) (or the mirror image for d1 > d2), zeta in u(m)"""
    m = min(d1, d2)
    out = []
    for z in u_basis(m):
        Z1, Z2 = zm(d1), zm(d2)
        if d1 <= d2:
            src, dst = z, mconj(z)
        else:
            src, dst = mconj(z), z
        for i in range(m):
            for j in range(m):
                Z1[i][j] = src[i][j]
                Z2[i][j] = dst[i][j]
        out.append((Z1, Z2))
    return out


def relabel(x, eps):
    """label of the operator pair x under theta_eps (conjugation is an involution)"""
    return (x[0] if eps[0] == 1 else mconj(x[0]), x[1] if eps[1] == 1 else mconj(x[1]))


# =====================================================================================================
# Part A
# =====================================================================================================
head("A1. Relabellings eps in {+1,-1}^2 (conjugation.md Thm 3)")
TYPES3 = [(2, 2), (2, 3), (3, 3)]
EPS = [(1, 1), (-1, -1), (1, -1), (-1, 1)]
for d1, d2 in TYPES3:
    gb, ss = g_basis(d1, d2), su_su(d1, d2)
    ref = Span([realify_m(drho(x)) for x in gb])
    ref0 = Span([realify_m(drho(x)) for x in ss])
    ok = True
    for e in EPS:
        img = [realify_m(drho(x, e)) for x in gb]
        img0 = [realify_m(drho(x, e)) for x in ss]
        ok &= rank(img) == ref.dim and all(ref.contains(v) for v in img)
        ok &= rank(img0) == ref0.dim and all(ref0.contains(v) for v in img0)
    check(ok and ref.dim == d1 * d1 + d2 * d2 - 1 and ref0.dim == d1 * d1 + d2 * d2 - 2,
          f"({d1},{d2}): d rho_eps(g) (dim {ref.dim}) and d rho_eps(su+su) (dim {ref0.dim}) are the same spaces for all four eps")
print("  Q = i(d2 1, -d1 1) in the Letter's labels h_W; its scalar on L under rho_eps, and acting = (i 1_L in d rho_eps(h_W)):")
ACT_EPS = {}
for e in EPS:
    row = []
    for d1, d2 in TYPES3:
        Q = (msc(G(0, d2), eye(d1)), msc(G(0, -d1), eye(d2)))
        M = drho(Q, e)
        scal = M[0][0]
        isscalar = all(M[i][j] == (scal if i == j else ZERO) for i in range(len(M)) for j in range(len(M)))
        hW = h_ab(d1, d2, 1, 1)
        act = acting(hW, d1, d2, e)
        ACT_EPS[(e, d1, d2)] = act
        check(isscalar and scal == G(0, e[0] * d2 - e[1] * d1) and act == (scal.im != 0),
              f"eps={e}, ({d1},{d2}): Q acts on L as {gstr(scal)} = i(e1 d2 - e2 d1); acting {act}")
for d1, d2 in TYPES3:
    check(ACT_EPS[((1, 1), d1, d2)] == ACT_EPS[((-1, -1), d1, d2)] and ACT_EPS[((1, -1), d1, d2)] == ACT_EPS[((-1, 1), d1, d2)],
          f"({d1},{d2}): the acting predicate depends on eps only through e1 e2")
# the image on L of the identification algebra h_W = s(u(d1)+u(d2)) itself does depend on e1 e2 at d1 = d2
for d1, d2 in TYPES3:
    hW = h_ab(d1, d2, 1, 1)
    im_p = [realify_m(drho(x, (1, 1))) for x in hW]
    im_m = [realify_m(drho(x, (1, -1))) for x in hW]
    same = same_space(im_p, im_m)
    check(same == (d1 != d2), f"({d1},{d2}): d rho_eps(h_W) is the same space for eps = (1,1) and (1,-1): {same} "
          f"(dims {rank(im_p)}, {rank(im_m)}); the image of the identification algebra is not relabelling-invariant at d1 = d2")
# trace functional under theta_eps
for d in (2, 3):
    ok = all(mtr(mconj(X)) == -mtr(X) for X in u_basis(d))
    check(ok, f"tr conj(X) = -tr X on a basis of u({d}), so theta_eps*(tr1 + tr2) = e1 tr1 + e2 tr2")
# transposition at both ends is eps = (-1,-1)
ok = True
for d1, d2 in TYPES3:
    for x in g_basis(d1, d2):
        ok &= mT(drho(x)) == msc(G(-1), drho(x, (-1, -1)))
check(ok, "(X1(x)1 + 1(x)X2)^T = -(conj X1 (x) 1 + 1 (x) conj X2) on bases at (2,2),(2,3),(3,3) (Remark (c))")
# conjugation on U(2): sigma_y conj(U) sigma_y = U on SU(2); inverts the centre
sy = [[G(), G(0, -1)], [G(0, 1), G()]]
a, b = G(Fr(3, 5)), G(0, Fr(4, 5))          # |a|^2 + |b|^2 = 1
U = [[a, -b.conj()], [b, a.conj()]]
check(mmul(mmul(sy, mconj(U)), sy) == U, "sigma_y conj(U) sigma_y = U for U = [[3/5, 4i/5],[4i/5, 3/5]] in SU(2)")
ph = G(Fr(3, 5), Fr(4, 5))
V = msc(ph, U)
check(mmul(mmul(sy, mconj(V)), sy) == msc(ph.conj() * ph.conj(), V),
      "for V = e^{i phi} U: sigma_y conj(V) sigma_y = e^{-2 i phi} V (the central U(1) is inverted)")
D3 = [[IU, ZERO, ZERO], [ZERO, IU, ZERO], [ZERO, ZERO, G(-1)]]
check(mtr(D3) == G(-1, 2) and mtr(mconj(D3)) == G(-1, -2) and (IU * IU * G(-1)) == G(1),
      "U = diag(i,i,-1) in SU(3): tr U = -1+2i != tr conj U = -1-2i, so conj is outer on SU(3) (Remark (b))")

head("A2. Stabilizer of the pair ket; centres at (2,3) under A and A* (Thm 4(1), Prop. 6, Thm S-stabilizer)")
for d1, d2 in [(2, 2), (2, 3), (3, 3), (2, 4), (3, 4)]:
    gb = g_basis(d1, d2)
    S = stab(gb, d1, d2)
    M0 = m0(d1, d2)
    w = omega(d1, d2)
    fix = all(all(z.iszero() for z in mvec(drho(x), w)) for x in M0)
    check(len(S) == d1 * d1 + (d2 - d1) ** 2 and fix and all(in_real_span(S, x) for x in M0),
          f"({d1},{d2}): dim stab_g(omega) = {len(S)} = d1^2 + (d2-d1)^2; m_0 (dim {rank([greal(x) for x in M0])}) fixes omega and lies in it")
d1, d2 = 2, 3
hW = h_ab(d1, d2, 1, 1)                      # the Letter's labels: tr X1 + tr X2 = 0
Uop = (msc(IU, eye(2)), [[-IU, ZERO, ZERO], [ZERO, -IU, ZERO], [ZERO, ZERO, ZERO]])
CENTRES = {}
for e, name in [((1, 1), "A "), ((1, -1), "A*")]:
    S = stab(hW, d1, d2, e)                   # elements of h_W (labels) annihilating omega under rho_eps
    Z = centre(S)
    check(len(S) == 4 and len(Z) == 1, f"{name}: h_W cap stab (labels) has dim {len(S)}; its centre has dim {len(Z)}")
    z = Z[0]
    dg = diag_of(z[0]) + diag_of(z[1])
    off = all(z[k][i][j].iszero() for k in (0, 1) for i in range(len(z[k])) for j in range(len(z[k])) if i != j)
    s = dg[0]
    c = Fr(1) / s.im                          # normalise so that the first entry is i
    lab = [G(x.re * c, x.im * c) for x in dg]
    CENTRES[name] = lab
    opx = relabel(z, e)
    opd = [G(x.re * c, x.im * c) for x in diag_of(opx[0]) + diag_of(opx[1])]
    onL = drho(z, e)
    diagL = [onL[i][i] for i in range(6)]
    prop = all(onL[i][j].iszero() for i in range(6) for j in range(6) if i != j) and \
        all(diagL[i].iszero() for i in range(6) if i % 3 != 2) and not diagL[2].iszero() and diagL[2] == diagL[5]
    print(f"      centre, labels: i diag({', '.join(str(x.im) for x in lab)}); operator pair: i diag({', '.join(str(x.im) for x in opd)})")
    check(off and prop, f"{name}: centre is diagonal and acts on L as a multiple of 1 (x) diag(0,0,1)")
check([x.im for x in CENTRES["A "]] == [1, 1, -1, -1, 0], "A: centre = U = i diag(1,1;-1,-1,0), as Thm S-stabilizer states")
check([x.im for x in CENTRES["A*"]] == [1, 1, 1, 1, -4], "A*: centre = i diag(1,1;1,1,-4) in labels (not U)")
SA = stab(hW, 2, 3, (1, 1))
SAs = stab(hW, 2, 3, (1, -1))
check(same_space([realify_m(drho(x)) for x in SA], [realify_m(drho(x, (1, -1))) for x in SAs]),
      "as operators on L the two stabilizers coincide (Prop. 6(3))")
M0SU = [(z, [[mconj(z)[i][j] if (i < 2 and j < 2) else G() for j in range(3)] for i in range(3)]) for z in su_basis(2)]
UM = (msc(IU, eye(2)), [[-IU if (i == j and i < 2) else G() for j in range(3)] for i in range(3)])
check(all(in_real_span(SAs, relabel(x, (1, -1))) for x in M0SU) and not in_real_span(SAs, relabel(UM, (1, -1))) and
      in_real_span(SA, UM),
      "A*: the label stabilizer contains the su(2) part of the matchings but not the matching phase U_m; A contains U_m (Prop. 6(2))")
# the Supplement's explicit basis (Sec. S-matrices): indices 0,1 qubit block, 2,3,4 qutrit block of C^2 + C^3
def E5(a, b, z=ONE):
    return Emat(5, a, b, z)


def split5(M):
    return ([[M[i][j] for j in range(2)] for i in range(2)], [[M[i][j] for j in range(2, 5)] for i in range(2, 5)])


e = {}
e[1] = msub(E5(0, 1), E5(1, 0)); e[2] = madd(E5(0, 1, IU), E5(1, 0, IU)); e[3] = madd(E5(0, 0, IU), E5(1, 1, -IU))
e[4] = msub(E5(2, 3), E5(3, 2)); e[5] = msub(E5(2, 4), E5(4, 2)); e[6] = msub(E5(3, 4), E5(4, 3))
e[7] = madd(E5(2, 3, IU), E5(3, 2, IU)); e[8] = madd(E5(2, 4, IU), E5(4, 2, IU)); e[9] = madd(E5(3, 4, IU), E5(4, 3, IU))
e[10] = madd(E5(2, 2, IU), E5(3, 3, -IU)); e[11] = madd(E5(3, 3, IU), E5(4, 4, -IU))
e[12] = [[G(0, 3) if (i == j and i < 2) else (G(0, -2) if i == j else G()) for j in range(5)] for i in range(5)]
E12 = [split5(e[k]) for k in range(1, 13)]
check(same_space([greal(x) for x in E12], [greal(x) for x in hW]), "Supplement e_1..e_12 span the Letter's h_W at (2,3)")
lin = lambda terms: split5(sum_m(terms))


def sum_m(terms):
    out = zm(5)
    for c, k in terms:
        out = madd(out, msc(G(c), e[k]))
    return out


NULL = [lin([(1, 1), (1, 4)]), lin([(1, 2), (-1, 7)]), lin([(1, 3), (-1, 10)]), lin([(1, 3), (2, 11), (-1, 12)])]
check(same_space([greal(x) for x in NULL], [greal(x) for x in SA]),
      "the Supplement's null-space basis e1+e4, e2-e7, e3-e10, e3+2e11-e12 spans the A stabilizer (labels = operators)")
w23 = omega(2, 3)
bad = [i for i, x in enumerate(NULL) if not all(z.iszero() for z in mvec(drho(x, (1, -1)), w23))]
check(bad == [1, 2, 3], "under A* (far block acting by its conjugate) e2-e7, e3-e10 and e3+2e11-e12 do not annihilate omega")
zc = lin([(-8, 10), (-16, 11), (8, 12)])
check([x.im for x in diag_of(zc[0]) + diag_of(zc[1])] == [24, 24, -24, -24, 0] and in_real_span(SA, zc),
      "-8e10-16e11+8e12 = 24 U lies in the A stabilizer (the Supplement's centre)")
# finite matching groups: G2 order 32, G3 order 384; pairs fixing the pair ket
def monomials(d):
    out = []
    phases = [G(1), IU, G(-1), -IU]
    for perm in itertools.permutations(range(d)):
        for ph in itertools.product(phases, repeat=d):
            M = zm(d)
            for j in range(d):
                M[perm[j]][j] = ph[j]
            out.append(M)
    return out


G2, G3 = monomials(2), monomials(3)
def det(M):
    if len(M) == 2:
        return M[0][0] * M[1][1] - M[0][1] * M[1][0]
    return (M[0][0] * (M[1][1] * M[2][2] - M[1][2] * M[2][1]) - M[0][1] * (M[1][0] * M[2][2] - M[1][2] * M[2][0])
            + M[0][2] * (M[1][0] * M[2][1] - M[1][1] * M[2][0]))


cntA = cntAs = cntS = 0
for g2 in G2:
    for g3 in G3:
        # (g2 (x) g3) omega = sum_{i<2} g2|i> (x) g3|i>
        okA = okAs = True
        for i in range(2):
            for j in range(3):
                vA = g2[i][0] * g3[j][0] + g2[i][1] * g3[j][1]
                vAs = g2[i][0] * g3[j][0].conj() + g2[i][1] * g3[j][1].conj()
                tgt = G(1) if (i == j) else G(0)
                okA &= vA == tgt
                okAs &= vAs == tgt
        cntA += okA
        cntAs += okAs
        cntS += okA and (det(g2) * det(g3) == G(1))
check(len(G2) == 32 and len(G3) == 384 and cntA == 128 and cntAs == 128,
      f"finite matchings: |G2| = {len(G2)}, |G3| = {len(G3)}; pairs fixing the pair ket: {cntA} (A), {cntAs} (A*), of {len(G2) * len(G3)}")
check(cntS == 32, f"of the {cntA} matching pairs fixing the pair ket (G3 acts on the whole qutrit, the third level's phase is free), "
      f"{cntS} lie in the identification group S(U(2) x U(3)) (det g2 det g3 = 1): the Letter's matchings are not all identifications")

head("A3. Matching closure (Thm 4(2),(3))")
RATIOS = [(1, 1), (2, 2), (-3, -3), (1, -1), (2, 3), (3, 2), (1, 0), (0, 1), (3, -2), (5, 4)]
for d1, d2 in TYPES3 + [(3, 2)]:
    M0 = m0(d1, d2)
    res = []
    for a_, b_ in RATIOS:
        h = h_ab(d1, d2, a_, b_)
        res.append((all(in_real_span(h, x) for x in M0)) == (a_ == b_))
    check(all(res), f"({d1},{d2}): m_0 in h_[a:b] iff a = b, for [a:b] in {RATIOS}")
for d1, d2 in TYPES3:
    hW = h_ab(d1, d2, 1, 1)
    M0 = m0(d1, d2)
    ok = True
    for e in EPS:
        inside = all(in_real_span(hW, relabel(x, e)) for x in M0)
        same_ops = all(drho(relabel(x, e), e) == drho(x) for x in M0)
        ok &= same_ops and (inside == (e[0] == e[1]))
    check(ok, f"({d1},{d2}): in eps-labels the matchings act on L as the same operators and lie in h_W iff e1 = e2")

head("A4. Generation (Thm 5)")
for d1, d2 in [(2, 2), (2, 3), (3, 3)]:
    gens = su_su(d1, d2) + m0(d1, d2)
    B = lie_closure(gens)
    hp = h_ab(d1, d2, 1, 1)
    check(len(B) == d1 * d1 + d2 * d2 - 1 and same_space([greal(x) for x in B], [greal(x) for x in hp]),
          f"({d1},{d2}): <su+su, m_0> has dim {len(B)} and equals s(u({d1})+u({d2})) = h_+")
    Bs = lie_closure(su_su(d1, d2) + stab(g_basis(d1, d2), d1, d2))
    want = d1 * d1 + d2 * d2 - (1 if d1 == d2 else 0)
    check(len(Bs) == want, f"({d1},{d2}): <su+su, full stabilizer> has dim {len(Bs)} ({'h_+' if d1 == d2 else 'all of g'})")
    ph = acting(su_su(d1, d2) + m0(d1, d2), d1, d2)
    check(ph == (d1 != d2), f"({d1},{d2}): i 1_L in d rho(su+su) + d rho(m_0): {ph} (iff d1 != d2)")
# the full stabilizer at (2,3) lies in no invariant hyperplane
S23 = stab(g_basis(2, 3), 2, 3)
tr_img = rank([{0: trs(x)[0], 1: trs(x)[1]} for x in S23])
check(tr_img == 2, "(2,3): the trace map X -> (tr X1, tr X2) has rank 2 on the full stabilizer, so no h_[a:b] contains it")

head("A5. Theorem 1's selection per reading, through config_enumeration.py's predicates (Thm 8, Cor. 9)")
buf = io.StringIO()
with contextlib.redirect_stdout(buf):
    import config_enumeration as CE
with open(os.path.join(HERE, "config_enumeration.out"), encoding="utf-8") as fh:
    check(buf.getvalue() == fh.read(), "config_enumeration.py imported unmodified; its printout equals config_enumeration.out")
TY4 = [(2, 2), (2, 3), (3, 2), (3, 3)]


def selection(act):
    sel = [k for k in CE.readable if CE.end_rule(k) and CE.two_point(k) and all(act[(k[0], e_)] for e_, _ in k[1])]
    return [f"c{k[0]}: " + " ".join(f"{e_}{t}" for e_, t in k[1]) for k in sel]


OBJ, QUB = "c3: 2R 2R 2U", "c2: 2R 2R 2U"
READINGS = {
    "A   h_+ (trace on H1+H2)": (lambda c, e_: h_ab(c, e_, 1, 1), (1, 1)),
    "A*  h_- (trace on H1+conj H2)": (lambda c, e_: h_ab(c, e_, 1, -1), (1, 1)),
    "A*  labels h_W, eps=(1,-1)": (lambda c, e_: h_ab(c, e_, 1, 1), (1, -1)),
    "B   (d rho)^-1 su(L)": (lambda c, e_: h_B(c, e_), (1, 1)),
    "C   all of g": (lambda c, e_: g_basis(c, e_), (1, 1)),
    "MC with the full stabilizer, P at (2,2)": (lambda c, e_: lie_closure(su_su(c, e_) + stab(g_basis(c, e_), c, e_)), (1, 1)),
    "W' = H1+conj H2 (1,-1)": (lambda c, e_: h_ab(c, e_, 1, -1), (1, 1)),
    "Lambda^2 L (d2,d1)(d1d2-1)": (lambda c, e_: h_ab(c, e_, e_ * (c * e_ - 1), c * (c * e_ - 1)), (1, 1)),
    "L' = H1 (x) conj H2 (d2,-d1)": (lambda c, e_: h_ab(c, e_, e_, -c), (1, 1)),
}
EXPECT = {"A   h_+ (trace on H1+H2)": [OBJ], "A*  h_- (trace on H1+conj H2)": [QUB, OBJ], "A*  labels h_W, eps=(1,-1)": [QUB, OBJ],
          "B   (d rho)^-1 su(L)": [], "C   all of g": [QUB, OBJ], "MC with the full stabilizer, P at (2,2)": [OBJ],
          "W' = H1+conj H2 (1,-1)": [QUB, OBJ], "Lambda^2 L (d2,d1)(d1d2-1)": [], "L' = H1 (x) conj H2 (d2,-d1)": [QUB, OBJ]}
for name, (mk, e_) in READINGS.items():
    act = {}
    mc = {}
    for c, f_ in TY4:
        h = mk(c, f_)
        act[(c, f_)] = acting(h, c, f_, e_)
        M0 = [relabel(x, e_) for x in m0(c, f_)]
        mc[(c, f_)] = all(in_real_span(h, x) for x in M0)
    sel = selection(act)
    print(f"      {name:42s} acts (2,2),(3,2): {act[(2, 2)]!s:5s} {act[(3, 2)]!s:5s}  MC (2,2),(3,2): {mc[(2, 2)]!s:5s} {mc[(3, 2)]!s:5s}  selected: {sel}")
    check(sel == EXPECT[name], f"{' '.join(name.split())}: selection {sel}")
# Theorem 8(3): fixed ratios in [-4,4]^2
cnt = {"obj": 0, "both": 0, "qub": 0, "none": 0}
okall = True
for a_ in range(-4, 5):
    for b_ in range(-4, 5):
        if a_ == 0 and b_ == 0:
            continue
        act = {}
        for c, f_ in TY4:
            act[(c, f_)] = acting(h_ab(c, f_, a_, b_), c, f_)
        sel = selection(act)
        key = "obj" if sel == [OBJ] else "both" if sel == [QUB, OBJ] else "qub" if sel == [QUB] else "none"
        cnt[key] += 1
        okall &= (key == "obj") == (a_ == b_) and (key == "qub") == (2 * b_ == 3 * a_)
check(okall and cnt == {"obj": 8, "both": 70, "qub": 2, "none": 0},
      f"fixed [a:b] over the 80 nonzero pairs in [-4,4]^2: {cnt}; object alone iff a = b, qubit configuration alone iff [a:b] = [2:3]")

head("A6. Trace calculus against explicit representation matrices; hypercharges (Thm 8(1), Prop. 1)")


def wedge2(M):
    n = len(M)
    idx = [(a, b) for a in range(n) for b in range(a + 1, n)]
    pos = {p: k for k, p in enumerate(idx)}
    R = zm(len(idx))
    for k, (a, b) in enumerate(idx):
        # X(e_a ^ e_b) = X e_a ^ e_b + e_a ^ X e_b
        for c in range(n):
            for (u, v, z) in ((c, b, M[c][a]), (a, c, M[c][b])):
                if z.iszero() or u == v:
                    continue
                sgn = 1 if u < v else -1
                p = pos[(min(u, v), max(u, v))]
                R[p][k] = R[p][k] + (z if sgn == 1 else -z)
    return R


def blockdiag(A, B):
    n, m = len(A), len(B)
    M = zm(n + m)
    for i in range(n):
        for j in range(n):
            M[i][j] = A[i][j]
    for i in range(m):
        for j in range(m):
            M[n + i][n + j] = B[i][j]
    return M


for d1, d2 in [(2, 2), (2, 3), (3, 3)]:
    n = d1 + d2
    oks = []
    for x in g_basis(d1, d2):
        t1, t2 = trs(x)
        W_ = blockdiag(x[0], x[1])
        Wp = blockdiag(x[0], mconj(x[1]))
        L_ = drho(x)
        Lp = drho(x, (1, -1))
        vals = {"Lambda^2 W": (mtr(wedge2(W_)).im, (n - 1) * (t1 + t2)),
                "Lambda^2 W'": (mtr(wedge2(Wp)).im, (n - 1) * (t1 - t2)),
                "L": (mtr(L_).im, d2 * t1 + d1 * t2),
                "Lambda^2 L": (mtr(wedge2(L_)).im, (d1 * d2 - 1) * (d2 * t1 + d1 * t2)),
                "L'": (mtr(Lp).im, d2 * t1 - d1 * t2)}
        oks.append(all(a_ == b_ for a_, b_ in vals.values()))
    check(all(oks), f"({d1},{d2}): tr of the representation matrices = Thm 8's polynomials for Lambda^2 W, Lambda^2 W', L, Lambda^2 L, L' on a basis of g")
# hypercharges at (2,3): Y = Q/(6i) on W, Lambda^2 W (off-diagonal block) and End W (Hom(H2,H1)) under ad
Q5 = blockdiag(msc(G(0, 3), eye(2)), msc(G(0, -2), eye(3)))
L2 = wedge2(Q5)
idx = [(a_, b_) for a_ in range(5) for b_ in range(a_ + 1, 5)]
offd = [k for k, (a_, b_) in enumerate(idx) if a_ < 2 <= b_]
ys = {L2[k][k].im / 6 for k in offd}
check(ys == {Fr(1, 6)} and all(L2[k][j].iszero() for k in range(10) for j in range(10) if k != j),
      "(2,3): Q acts on the H1(x)H2 block of Lambda^2 W as 6i Y with Y = 1/6 (the 10)")
yad = {(Q5[a_][a_] - Q5[b_][b_]).im / 6 for a_ in range(2) for b_ in range(2, 5)}
check(yad == {Fr(5, 6)}, "(2,3): ad Q on Hom(H2,H1) = E_{a,b}, a in H1, b in H2: Y = 5/6 (the 24)")

head("A7. The clause P relative to the Letter's definition; Prop. 6(4)")
for d1, d2 in TYPES3:
    hW = h_ab(d1, d2, 1, 1)
    proper = len(hW) == d1 * d1 + d2 * d2 - 1
    ok = proper
    for e in EPS:
        mc = all(in_real_span(hW, relabel(x, e)) for x in m0(d1, d2))
        ok &= mc == (e[0] == e[1])
    check(ok, f"({d1},{d2}): h_W is a proper subalgebra (dim {len(hW)}) under every eps, and MC in labels <=> e1 = e2")
sel_ok = True
for e in EPS:
    act = {(c, f_): acting(h_ab(c, f_, 1, 1), c, f_, e) for c, f_ in TY4}
    mc = all(all(in_real_span(h_ab(c, f_, 1, 1), relabel(x, e)) for x in m0(c, f_)) for c, f_ in TY4)
    sel = selection(act)
    sel_ok &= (sel == [OBJ]) == mc
check(sel_ok, "within the Letter's family (h_W, any eps): Theorem 1 selects the object alone <=> MC (P holds automatically)")
# Prop. 6(4) at (2,2): centres of h_+, h_- (operator pairs) and of the label h_W
for name, (a_, b_) in [("h_+", (1, 1)), ("h_-", (1, -1))]:
    Z = centre(h_ab(2, 2, a_, b_))
    z = Z[0]
    x_, y_ = z[0][0][0].im, z[1][0][0].im
    has_equal = (x_ == y_)
    check(len(Z) == 1 and has_equal == (name == "h_-"),
          f"(2,2): centre of {name} (operator pairs) is the line of i({x_}, {y_}); of the form (ix 1, ix 1), x != 0: {has_equal}")
Zl = centre(h_ab(2, 2, 1, 1))[0]
check(Zl[0][0][0].im == -Zl[1][0][0].im and drho(Zl, (1, -1))[0][0].im == 2 * Zl[0][0][0].im,
      "(2,2), A* in labels: the centre of h_W is (i, -i) (no nonzero (ix, ix)), and it acts on L as 2i under rho_(1,-1)")

# =====================================================================================================
# Part B
# =====================================================================================================
head("B1. Model: real letters (Lean's tracelessBasis), structure constants, P<=2")


def basis_d(d):
    """identity first, then Lean's tracelessBasis d: symmetric, antisymmetric off-diagonal, diag(1..1,-k,0..)"""
    def Z():
        return [[Fr(0)] * d for _ in range(d)]
    out = []
    I_ = Z()
    for i in range(d):
        I_[i][i] = Fr(1)
    out.append(I_)
    pairs = [(i, j) for i in range(d) for j in range(i + 1, d)]
    for i, j in pairs:
        m = Z(); m[i][j] = m[j][i] = Fr(1); out.append(m)
    for i, j in pairs:
        m = Z(); m[i][j] = Fr(1); m[j][i] = Fr(-1); out.append(m)
    for k in range(1, d):
        m = Z()
        for a_ in range(k):
            m[a_][a_] = Fr(1)
        m[k][k] = Fr(-k)
        out.append(m)
    return out


def rmul(A, B):
    n = len(A)
    return [[sum(A[i][k] * B[k][j] for k in range(n)) for j in range(n)] for i in range(n)]


def rT(A):
    return [list(c) for c in zip(*A)]


def rpair(A, B):
    return sum(A[i][j] * B[i][j] for i in range(len(A)) for j in range(len(A)))


BAS, NORM, TSG, SC = {}, {}, {}, {}
for d in (2, 3, 4):
    B = basis_d(d)
    BAS[d] = B
    NORM[d] = [rpair(b, b) for b in B]
    TSG[d] = [1 if rT(b) == b else -1 for b in B]
    assert all(rT(b) == [[-x for x in r] for r in b] for b, s in zip(B, TSG[d]) if s == -1)
    orth = all(rpair(B[i], B[j]) == 0 for i in range(len(B)) for j in range(len(B)) if i != j)
    sc = {}
    for i, j in itertools.product(range(len(B)), repeat=2):
        P = rmul(B[i], B[j])
        coef = {k: rpair(B[k], P) / NORM[d][k] for k in range(len(B))}
        coef = {k: c for k, c in coef.items() if c}
        rec = [[sum(c * B[k][r][s] for k, c in coef.items()) for s in range(d)] for r in range(d)]
        assert rec == P
        sc[(i, j)] = coef
    SC[d] = sc
    check(orth and len(B) == d * d, f"d = {d}: identity + {d * d - 1} real letters, orthogonal under tr(A^T B); products expanded exactly")

DIMS = (3, 2, 3, 2)          # Q1, q1, Q2, q2
SLOT_NAMES = ("Q1", "q1", "Q2", "q2")


def nlet(t):
    return sum(1 for x in t if x)


def NQ(t):
    return (t[0] != 0) + (t[2] != 0)


def Nq(t):
    return (t[1] != 0) + (t[3] != 0)


def omul(X, Y):
    out = {}
    for t, x in X.items():
        for u, y in Y.items():
            acc = {(): x * y}
            for s in range(4):
                nxt = {}
                for pre, c in acc.items():
                    for k, cc in SC[DIMS[s]][(t[s], u[s])].items():
                        key = pre + (k,)
                        nxt[key] = nxt.get(key, 0) + c * cc
                acc = nxt
            for k, c in acc.items():
                if c:
                    v = out.get(k, 0) + c
                    if v:
                        out[k] = v
                    else:
                        out.pop(k)
    return out


def oadd(X, Y, c=1):
    out = dict(X)
    for k, y in Y.items():
        v = out.get(k, 0) + c * y
        if v:
            out[k] = v
        else:
            out.pop(k, None)
    return out


def osc(c, X):
    return {k: c * x for k, x in X.items()} if c else {}


def comm(X, Y):
    return oadd(omul(X, Y), omul(Y, X), -1)


def single(s, a):
    t = [0, 0, 0, 0]
    t[s] = a
    return {tuple(t): Fr(1)}


def ident():
    return {(0, 0, 0, 0): Fr(1)}


PLE2 = [t for t in itertools.product(range(9), range(4), range(9), range(4)) if nlet(t) <= 2]
P2T = [t for t in PLE2 if nlet(t) == 2]
VQT = [t for t in PLE2 if nlet(t) == 1 and (t[0] or t[2])]
VqT = [t for t in PLE2 if nlet(t) == 1 and (t[1] or t[3])]
vec = lambda t: {t: Fr(1)}
P2 = [vec(t) for t in P2T]
VQ = [vec(t) for t in VQT]
Vq = [vec(t) for t in VqT]
ID = [ident()]
check(len(PLE2) == 192 and len(P2T) == 169 and len(VQT) == 16 and len(VqT) == 6,
      "P<=2 has dim 192 = 1 + 16 + 6 + 169; slot-pair sizes " +
      str(sorted(sum(1 for t in P2T if tuple(i for i in range(4) if t[i]) == p) for p in itertools.combinations(range(4), 2))))


def weight(t):
    w = Fr(1)
    for s in range(4):
        w *= NORM[DIMS[s]][t[s]]
    return w


def ipair(X, Y):
    return sum(x * Y.get(k, 0) * weight(k) for k, x in X.items())


# exchanges
def Eplain(X):
    return {(t[2], t[3], t[0], t[1]): x for t, x in X.items()}


def tsign(t):
    s = 1
    for i in range(4):
        s *= TSG[DIMS[i]][t[i]]
    return s


def make_f(eQ, eq):
    def f(X):
        return {(t[2], t[3], t[0], t[1]): x * tsign(t) * (eQ ** NQ(t)) * (eq ** Nq(t)) for t, x in X.items()}
    return f


FR = {"E": Eplain, "E'": make_f(1, 1), "Fm": make_f(1, -1), "FQ": make_f(-1, 1), "sE'": make_f(-1, -1)}


def eig_space(f, vecs, sign):
    """sign-eigenspace of the involution f on span(vecs) (f maps the span to itself)"""
    out = [oadd(v, f(v), sign) for v in vecs]
    return [v for v in out if v]


head("B2. The order-sign domain D (records-uniform.md Thm 2)")
LET = [(s, a_) for s in range(4) for a_ in range(1, DIMS[s] ** 2)]
check(len(LET) == 22, "22 letters (8, 3, 8, 3)")
Ssym, Sant = Span(), Span()
MSYM, MANT, MALL = [], [], []
for (s, a_), (u, b_) in itertools.product(LET, repeat=2):
    xy = omul(single(s, a_), single(u, b_))
    yx = omul(single(u, b_), single(s, a_))
    MALL.append(xy)
    MSYM.append(oadd(xy, yx))
    MANT.append(oadd(xy, yx, -1))
dall, dsym, dant = rank(MALL), rank(MSYM), rank(MANT)
SSYM, SANT = Span(MSYM), Span(MANT)
check(dall == 192 and dsym == 186 and dant == 22, f"dim m(V(x)V) = {dall}, dim m(Sym^2 V) = {dsym}, dim m(Lambda^2 V) = {dant}")
check(same_space(MSYM, ID + P2 + VQ) and same_space(MANT, VQ + Vq),
      "m(Sym^2 V) = 1 + P2 + V_Q and m(Lambda^2 V) = V_Q + V_q")
Adim = inter_dim(MSYM, MANT)
check(Adim == 16 and all(SSYM.contains(v) and SANT.contains(v) for v in VQ), "A = m(Sym^2) cap m(Lambda^2) = V_Q (dim 16)")
DT = [t for t in PLE2 if t not in set(VQT)]
check(len(DT) == 176 and all(ipair(vec(t), a_) == 0 for t in DT for a_ in VQ),
      "D = A^perp in P<=2 = 1 + P2 + V_q (dim 176 = 192 - 16, every element orthogonal to A)")
check(all((SSYM.contains(vec(t)) != SANT.contains(vec(t))) and (SSYM.contains(vec(t)) == (nlet(t) % 2 == 0)) for t in DT),
      "each basis operator of D lies in exactly one of m(Sym^2), m(Lambda^2); the order sign equals the own sign (-1)^N")
for d in (2, 3, 4):
    B = BAS[d]
    nonzero = False
    for i, j in itertools.product(range(1, d * d), repeat=2):
        P = rmul(B[i], B[j]); Q_ = rmul(B[j], B[i])
        anti = [[P[r][c] + Q_[r][c] for c in range(d)] for r in range(d)]
        tr_ = sum(anti[r][r] for r in range(d))
        trl = [[anti[r][c] - (tr_ / d if r == c else 0) for c in range(d)] for r in range(d)]
        nonzero |= any(x != 0 for r in trl for x in r)
    check(nonzero == (d != 2), f"d = {d}: the traceless part of {{A,B}} vanishes for all letters: {not nonzero}")
# counterexample to Thm 2(3) sufficiency: B = A + span{1 + sigma_x^(q1)} contains A but the flip does not descend
sx, aY, sz = 1, 2, 3                         # qubit letters: sigma_x, i sigma_y (real), sigma_z
xp = [(Fr(1), (1, sx), (1, sx))]             # symmetric tensor sigma_x (x) sigma_x at q1
xm = [(Fr(-1, 2), (1, aY), (1, sz)), (Fr(1, 2), (1, sz), (1, aY))]   # antisymmetric tensor


def m_of(tensor, flip=False):
    out = {}
    for c, x, y in tensor:
        if flip:
            x, y = y, x
        out = oadd(out, osc(c, omul(single(*x), single(*y))))
    return out


mx = oadd(m_of(xp), m_of(xm))
mtx = oadd(m_of(xp, True), m_of(xm, True))
Bspace = Span(VQ + [oadd(ident(), single(1, sx))])
check(mx == oadd(ident(), single(1, sx)) and mtx == oadd(ident(), single(1, sx), -1),
      "x = sx(x)sx - (iy(x)sz - sz(x)iy)/2 at q1: m(x) = 1 + sigma_x^(q1), m(flip x) = 1 - sigma_x^(q1)")
check(Bspace.contains(mx) and not Bspace.contains(mtx),
      "B = A + span{1 + sigma_x^(q1)} contains A, m(x) in B, m(flip x) not in B: Thm 2(3)'s 'iff B contains A' fails (only B = A, or B/A flip-stable, works)")

head("B3. C_L, C, C+ and Lie(G_Omega) (Prop. 1)")
CL = eig_space(FR["E"], Vq, -1)
CLp = eig_space(FR["E"], Vq, +1)
C = eig_space(FR["E'"], Vq, -1)
Cp = eig_space(FR["E'"], Vq, +1)
check(rank(CL) == rank(C) == rank(Cp) == 3, "C_L = V_q^{E=-1}, C = V_q^{E'=-1}, C+ = V_q^{E'=+1} are 3-dimensional")
check(inter_dim(CL, C) == 2 and inter_dim(CL, Cp) == 1 and inter_dim(C, Cp) == 0,
      f"dim(C_L cap C) = {inter_dim(CL, C)}, dim(C_L cap C+) = {inter_dim(CL, Cp)}, dim(C cap C+) = {inter_dim(C, Cp)}")
Dx = oadd(single(1, sx), single(3, sx), -1); Dz = oadd(single(1, sz), single(3, sz), -1); Dy = oadd(single(1, aY), single(3, aY), -1)
Sx = oadd(single(1, sx), single(3, sx)); Sz = oadd(single(1, sz), single(3, sz)); Sy = oadd(single(1, aY), single(3, aY))
check(same_space(C, [Dx, Sy, Dz]) and same_space(Cp, [Sx, Dy, Sz]) and same_space(CL, [Dx, Dy, Dz]),
      "C = {D_x, S_y, D_z}, C+ = {S_x, D_y, S_z}, C_L = {D_x, D_y, D_z}")
# the Letter's commutator records [s_mu, s_nu]_q1 - [s_mu, s_nu]_q2
LCR = [oadd(comm(single(1, i), single(1, j)), comm(single(3, i), single(3, j)), -1) for i, j in [(1, 2), (1, 3), (2, 3)]]
check(same_space(LCR, CL), "the Letter's three records [s_mu,s_nu]_q1 - [s_mu,s_nu]_q2 span C_L")
LIE = [oadd(single(s, a_), osc(TSG[DIMS[s]][a_], single(s + 2, a_)), -1) for s in (0, 1) for a_ in range(1, DIMS[s] ** 2)]
check(rank(LIE) == 11 and same_space(LIE, eig_space(FR["E'"], VQ + Vq, -1)),
      "Lie(G_Omega) = span{a^(1) - (a^T)^(2)} (dim 11) = the E'-odd single letters")
check(same_space([v for v in LIE if all(k[0] == 0 and k[2] == 0 for k in v)], C), "C is the qubit part of Lie(G_Omega)")
SLIE = Span(LIE)
check(inter_dim(CL, LIE) == 2 and SLIE.contains(Dx) and SLIE.contains(Dz) and not SLIE.contains(Dy) and
      inter_dim(CLp, LIE) == 1 and SLIE.contains(Sy) and not SLIE.contains(Sx) and not SLIE.contains(Sz),
      "plain frame: C_L cap Lie(G_Omega) = {D_x, D_z}, C_L+ cap Lie(G_Omega) = {S_y}")

head("B4. Frames: classes and values under the Letter's weight law (Prop. 1(4), Sec. 8)")


def cls(t):
    n = nlet(t)
    if n == 0:
        return "id"
    if n == 1:
        return "Q" if (t[0] or t[2]) else "q"
    if NQ(t) == 2:
        return "QQ"
    if Nq(t) == 2:
        return "qq"
    return "mix"


def cell_dims(f, tuples):
    """dims of the f-eigenspaces per class, for the span of basis tuples (closed under f up to sign)"""
    out = {}
    seen = set()
    for t in tuples:
        if t in seen:
            continue
        img = f(vec(t))
        (u, s), = img.items()
        c = cls(t)
        if u == t:
            seen.add(t)
            key = (c, 1 if s == 1 else -1)
            out[key] = out.get(key, 0) + 1
        else:
            seen.add(t); seen.add(u)
            out[(c, 1)] = out.get((c, 1), 0) + 1
            out[(c, -1)] = out.get((c, -1), 0) + 1
    return out


R2, R3, R6 = sp.sqrt(2), sp.sqrt(3), sp.sqrt(6)


def inv_value(num):
    a_, b_, c_, d_ = num
    return sp.N(sp.Integer(31104) / (a_ + b_ * R2 + c_ * R3 + d_ * R6), 15)


def classes_of(f, single_space):
    """the Letter's classes for R = P2 + single_space in frame f: (f-odd products, f-even mixed, f-even same-end, f-odd single qubit)"""
    cd = cell_dims(f, P2T)
    odd = sum(v for (c, p), v in cd.items() if p == -1)
    em = cd.get(("mix", 1), 0)
    es = cd.get(("QQ", 1), 0) + cd.get(("qq", 1), 0)
    comm_ = rank(single_space) if single_space and all(f(v) == osc(-1, v) for v in single_space) else 0
    return odd, em, es, comm_


VAL = {}
for fname, sing in [("E", CL), ("E'", C), ("Fm", eig_space(FR["Fm"], Vq, -1)), ("FQ", eig_space(FR["FQ"], Vq, -1)),
                    ("sE'", eig_space(FR["sE'"], Vq, -1))]:
    k = classes_of(FR[fname], sing)
    num = (k[0], k[1], k[2], k[3])          # weights 1, sqrt2, sqrt3, sqrt6
    VAL[fname] = inv_value(num)
    check(k == (79, 48, 42, 3) and abs(VAL[fname] - sp.Float("137.035999189566", 15)) < 1e-9,
          f"R_{fname} = P2 + V_q^({fname}=-1): classes {k}, numerator 79+48r2+42r3+3r6, 1/value = {VAL[fname]}")
check(same_space(eig_space(FR["Fm"], Vq, -1), Cp) and same_space(eig_space(FR["sE'"], Vq, -1), Cp) and
      same_space(eig_space(FR["FQ"], Vq, -1), C), "V_q^{f=-1}: C for E' and (-1)^NQ E'; C+ for Fm and sigma E'")
# P2 + C+ sorted by E': operator form M_L and the word rule
def ML(c, sig, p, nQ, nq):
    S = 1 if sig * p == 1 else 0
    return (3 if S and nQ * nq == 0 else 1) * (2 if S and nq == 1 else 1)


k = classes_of(FR["E'"], [])
# C+ is E'-even single-qubit (sigma = -1, parity +1): S = 0 under M_L (k = 1); word rule: X survives (3), Y does not
nev = rank(Cp) if all(FR["E'"](v) == v for v in Cp) else -1
numML = (k[0] + nev * (ML("q", -1, 1, 0, 1) == 1), k[1], k[2], 0)
numW = (k[0], k[1], k[2] + nev, 0)
check(numML == (82, 48, 42, 0) and str(inv_value(numML))[:7] == "139.712",
      f"P2 + C+ sorted by E', M_L: numerator {numML}, 1/value {inv_value(numML)}")
check(numW == (79, 48, 45, 0) and str(inv_value(numW))[:7] == "138.347",
      f"P2 + C+ sorted by E', word rule: numerator {numW}, 1/value {inv_value(numW)}")
check(ML("q", -1, -1, 0, 1) == 6 and ML("Q", -1, -1, 1, 0) == 3,
      "M_L gives an f-odd single qubit letter k = 6 and an f-odd single qutrit letter k = 3")
for nm, num, want in [("P2 only", (79, 48, 42, 0), "141.621"), ("D_0 under E', M_L", (82, 48, 42, 3), "135.248"),
                      ("D_0, word rule", (79, 48, 45, 3), "133.969"), ("P2 + Lie(G_Omega)", (79, 48, 50, 3), "129.151"),
                      ("P<=2,0, M_L", (90, 48, 50, 3), "123.510"), ("P<=2,0, word", (79, 48, 61, 3), "119.683")]:
    check(str(inv_value(num))[:7] == want, f"variant {nm}: numerator {num}, 1/value {inv_value(num)}")
# Thm 4 cells on D_0 for E: (sigma, f) classes
for fname in ("E", "E'", "Fm"):
    cd = cell_dims(FR[fname], P2T)
    pp = sum(v for (c, p), v in cd.items() if p == 1)
    pm = sum(v for (c, p), v in cd.items() if p == -1)
    mp = rank(eig_space(FR[fname], Vq, 1)); mm = rank(eig_space(FR[fname], Vq, -1))
    check((pp, pm, mp, mm) == (90, 79, 3, 3), f"{fname}: (sigma, f) classes on D_0 = {(pp, pm, mp, mm)}")

head("B5. Exchanges: G_Omega covariance, involution, algebra type, Schur intertwiners, Fm = Ad_W (Thm 7)")
PLE2V = [vec(t) for t in PLE2]
for fname, f in FR.items():
    inv = all(f(f(v)) == v for v in PLE2V)
    cov = True
    for K in LIE:
        for v in PLE2V:
            if f(comm(K, v)) != comm(K, f(v)):
                cov = False
                break
        if not cov:
            break
    singles = [ident()] + [single(s, a_) for s, a_ in LET]
    aut = all(f(omul(x, y)) == omul(f(x), f(y)) for x in singles for y in singles)
    anti = all(f(omul(x, y)) == omul(f(y), f(x)) for x in singles for y in singles)
    typ = "automorphism" if aut else "anti-automorphism" if anti else "neither"
    want = {"E": (False, "automorphism"), "E'": (True, "anti-automorphism")}.get(fname, (True, "neither"))
    check(inv and cov == want[0] and typ == want[1],
          f"{fname}: involution {inv}; commutes with ad_K for all 11 generators K of Lie(G_Omega) on P<=2: {cov}; {typ}")
# Schur: intertwiners V_s -> V_s' with phi(ad_a b) = ad_{-a^T} phi(b)
for s, d in ((1, 2), (0, 3)):
    n = d * d - 1
    eqs = {}
    # unknown phi[i][j]: phi(b_j) = sum_i phi[i][j] b_i
    for a_ in range(1, n + 1):
        for j in range(1, n + 1):
            lhs = SC[d][(a_, j)]; rhs = SC[d][(j, a_)]
            adab = {k: lhs.get(k, 0) - rhs.get(k, 0) for k in set(lhs) | set(rhs)}
            # phi(ad_a b_j) = sum_k adab[k] phi(b_k) = sum_k adab[k] sum_i phi[i][k] b_i
            # ad_{-a^T}(phi(b_j)) = -t_a sum_i phi[i][j] [b_a, b_i]
            for out in range(1, n + 1):
                row = {}
                for kk, c in adab.items():
                    if kk and c:
                        row[(out, kk)] = row.get((out, kk), 0) + c
                for i in range(1, n + 1):
                    l2 = SC[d][(a_, i)]; r2 = SC[d][(i, a_)]
                    c = l2.get(out, 0) - r2.get(out, 0)
                    if c:
                        row[(i, j)] = row.get((i, j), 0) + TSG[d][a_] * c
                eqs[(a_, j, out)] = {k: v for k, v in row.items() if v}
    unknowns = [(i, j) for i in range(1, n + 1) for j in range(1, n + 1)]
    cols = [{key: eq.get(u, 0) for key, eq in eqs.items() if eq.get(u, 0)} for u in unknowns]
    ns = nullspace(cols)
    sol = ns[0] if len(ns) == 1 else {}
    sol = {unknowns[i]: c for i, c in sol.items()}
    c0 = sol.get((1, 1), 0)
    isT = len(ns) == 1 and c0 != 0 and all(sol.get((i, j), 0) == (c0 * TSG[d][i] * TSG[d][1] if i == j else 0) for i, j in unknowns)
    check(isT, f"d = {d}: the G_Omega-intertwiners V_{SLOT_NAMES[s]} -> V_{SLOT_NAMES[s + 2]} form a line, spanned by b -> b^T")
# Fm on qubit-only operators = Ad_W, W = (sigma_y (x) sigma_y) SWAP (real orthogonal 4x4)
B2 = BAS[2]


def k4(A, B):
    return [[A[i // 2][j // 2] * B[i % 2][j % 2] for j in range(4)] for i in range(4)]


SYY = [[0, 0, 0, -1], [0, 0, 1, 0], [0, 1, 0, 0], [-1, 0, 0, 0]]
SWAP = [[1, 0, 0, 0], [0, 0, 1, 0], [0, 1, 0, 0], [0, 0, 0, 1]]
W4 = rmul([[Fr(x) for x in r] for r in SYY], [[Fr(x) for x in r] for r in SWAP])
ok = True
for i, j in itertools.product(range(4), repeat=2):
    lhs = rmul(rmul(W4, k4(B2[i], B2[j])), rT(W4))
    img = FR["Fm"]({(0, i, 0, j): Fr(1)})
    (t, c), = img.items()
    rhs = [[c * x for x in r] for r in k4(B2[t[1]], B2[t[3]])]
    ok &= lhs == rhs
check(ok, "Fm acts on the 16 qubit-only basis operators of P<=2 as Ad_W, W = (sigma_y (x) sigma_y) SWAP")
Vm = [[Fr(0), Fr(1)], [Fr(-1), Fr(0)]]                       # V = i sigma_y in SU(2), real, so conj V = V
VV = k4(Vm, Vm)
check([[Fr(x) for x in r] for r in SYY] == [[-x for x in r] for r in VV],
      "sigma_y (x) sigma_y = -(V (x) conj V), V = i sigma_y: on the read qubits Fm = Ad_(V (x) conj V) o (plain swap), a G_Omega matching after E")

head("B6. Identification rule (Thm 6)")
IB = ID + [single(s, a_) for s, a_ in LET]           # I = P<=1, dim 23
check(rank(IB) == 23, "I = P<=1 has dim 23")
# stab_Omega
def apply_omega(X):
    """X|Omega>, Omega = sum_{i<3, j<2} |i j i j> (slot order Q1 q1 Q2 q2), as dict over (k0,k1,k2,k3)"""
    out = {}
    for t, x in X.items():
        M = [BAS[DIMS[s]][t[s]] for s in range(4)]
        A = rmul(M[0], rT(M[2]))
        Bq = rmul(M[1], rT(M[3]))
        for k0, k2 in itertools.product(range(3), repeat=2):
            if A[k0][k2] == 0:
                continue
            for k1, k3 in itertools.product(range(2), repeat=2):
                c = A[k0][k2] * Bq[k1][k3]
                if c:
                    key = (k0, k1, k2, k3)
                    v = out.get(key, 0) + x * c
                    if v:
                        out[key] = v
                    else:
                        out.pop(key)
    return out


OM = apply_omega(ident())
ns = nullspace([apply_omega(X) for X in IB] + [OM])
STAB = [combo({i: c for i, c in v.items() if i < 23}, IB) for v in ns]
STAB = [v for v in STAB if v]
check(rank(STAB) == 12 and same_space(STAB, ID + LIE), "stab_Omega = {X in I : X Omega in C Omega} = 1 + Lie(G_Omega) (dim 12)")


def comp(f):
    vecs = []
    for X in IB:
        v = {}
        for j, Y in enumerate(PLE2V):
            for k_, c in oadd(f(comm(X, f(Y))), comm(X, Y), -1).items():
                v[(j,) + k_] = c
        vecs.append(v)
    return [combo(c, IB) for c in nullspace(vecs)]


IE = eig_space(FR["E"], IB, +1)
COMP = {}
for fname in FR:
    COMP[fname] = comp(FR[fname])
check(same_space(COMP["E"], IE) and rank(IE) == 12 and not same_space(COMP["E"], STAB),
      "Comp_E = I^E (the E-even identifications, dim 12) != stab_Omega")
check(all(same_space(COMP[f], STAB) for f in ("E'", "Fm", "FQ", "sE'")), "Comp_f = stab_Omega for the four covariant exchanges")
DI = ID + Vq                                         # D cap I; P2 is orthogonal to I automatically


def perp_in(space, against):
    """basis of {v in span(space): <v, a> = 0 for all a in against}"""
    vecs = [{j: ipair(v, a_) for j, a_ in enumerate(against)} for v in space]
    return [combo(c, space) for c in nullspace(vecs)]


check(all(ipair(p, x) == 0 for p in P2 for x in IB), "P2 is orthogonal to I")
for fname, want in [("E", CL), ("E'", Cp), ("Fm", Cp), ("FQ", Cp), ("sE'", Cp)]:
    part = perp_in(DI, COMP[fname])
    check(same_space(part, want) and rank(part) == 3,
          f"{fname}: D cap Comp_f^perp = P2 + {'C_L (the Letter 172)' if fname == 'E' else 'C+'}")
for fname, want, lab in [("E", CL, "C_L"), ("E'", C, "C"), ("Fm", Cp, "C+")]:
    If = eig_space(FR[fname], IB, +1)
    part = perp_in(DI, If)
    check(same_space(part, want), f"{fname}: D cap (I^f)^perp = R_f (single-letter part {lab})")
DIE = eig_space(FR["E'"], DI, +1)                   # 1 + C+
br = [comm(x, y) for x, y in itertools.combinations(Cp, 2)]
check(same_space(br, C), "[C+, C+] = C: D cap I^{E'} = 1 + C+ is not a Lie algebra")
check(Span(C).contains(comm(C[0], C[1])) and all(Span(CLp).contains(comm(x, y)) for x, y in itertools.combinations(CLp, 2)),
      "[C, C] in C (1 + C = D cap I^{Fm} is a Lie algebra) and [C_L+, C_L+] in C_L+ (1 + C_L+ = D cap I^E)")

head("B7. Anticommuting letters (Prop. 3); the qutrit and d = 4")
I_ = sp.I


def gell_mann():
    E_ = lambda a_, b_: sp.Matrix(3, 3, lambda i, j: 1 if (i, j) == (a_, b_) else 0)
    return [E_(0, 1) + E_(1, 0), -I_ * E_(0, 1) + I_ * E_(1, 0), E_(0, 0) - E_(1, 1), E_(0, 2) + E_(2, 0),
            -I_ * E_(0, 2) + I_ * E_(2, 0), E_(1, 2) + E_(2, 1), -I_ * E_(1, 2) + I_ * E_(2, 1),
            (E_(0, 0) + E_(1, 1) - 2 * E_(2, 2)) / sp.sqrt(3)]


GM = gell_mann()
pairs = [(i + 1, j + 1) for i, j in itertools.combinations(range(8), 2) if sp.simplify(GM[i] * GM[j] + GM[j] * GM[i]) == sp.zeros(3)]
check(pairs == [(1, 2), (1, 3), (2, 3), (4, 5), (6, 7)], f"Gell-Mann: anticommuting pairs {pairs}")
cm = [sp.simplify(GM[i - 1] * GM[j - 1] - GM[j - 1] * GM[i - 1]) for i, j in pairs]
flat = lambda Ms: sp.Matrix([[x for x in M] for M in Ms])
check(flat(cm).rank() == 4 and flat(cm + [GM[0], GM[1], GM[2], GM[7]]).rank() == 4,
      "their commutators span {lambda_1, lambda_2, lambda_3, lambda_8} (dim 4)")
check([sp.simplify(GM[i].det()) == 0 for i in range(8)] == [True] * 7 + [False], "lambda_1..lambda_7 singular, lambda_8 invertible")
nqut = flat(cm).rank()
nGM = (79, 48, 42 + nqut, 3)               # the qutrit differences: own sign -1, E-odd, k = 3 under M_L(E) and the word rule
check(nGM == (79, 48, 46, 3) and str(inv_value(nGM))[:7] == "132.977",
      f"Gell-Mann literal: R_E + {nqut} E-odd qutrit differences at k = 3: numerator {nGM}, 1/value {inv_value(nGM)}")
check(ML("Q", -1, -1, 1, 0) == 3, "an E-odd single qutrit letter has k = 3 under M_L (Pi_Y = 0: the centre is not moved)")
lean_pairs = [(i, j) for i, j in itertools.combinations(range(1, 9), 2)
              if all(x == 0 for r in [[a_ + b_ for a_, b_ in zip(r1, r2)] for r1, r2 in zip(rmul(BAS[3][i], BAS[3][j]), rmul(BAS[3][j], BAS[3][i]))] for x in r)]
check(len(lean_pairs) == 5, f"Lean's tracelessBasis 3: {len(lean_pairs)} anticommuting pairs {lean_pairs}")


def weyl(d):
    w = sp.exp(2 * sp.pi * I_ / d)
    X = sp.Matrix(d, d, lambda i, j: 1 if i == (j + 1) % d else 0)
    Z = sp.diag(*[sp.nsimplify(sp.expand_complex(w ** j)) for j in range(d)])
    return [(a_, b_, sp.simplify(X ** a_ * Z ** b_)) for a_ in range(d) for b_ in range(d) if (a_, b_) != (0, 0)]


for d in (3, 4):
    Wl = weyl(d)
    ap = [((a1, b1), (a2, b2)) for (a1, b1, M1), (a2, b2, M2) in itertools.combinations(Wl, 2)
          if sp.simplify(M1 * M2 + M2 * M1) == sp.zeros(d)]
    cms = [sp.simplify(M1 * M2 - M2 * M1) for (a1, b1, M1), (a2, b2, M2) in itertools.combinations(Wl, 2)
           if sp.simplify(M1 * M2 + M2 * M1) == sp.zeros(d)]
    rk = flat(cms).rank() if cms else 0
    WEYL_RANK = rk
    if d == 3:
        check(len(ap) == 0, "Weyl letters X^a Z^b at d = 3: no anticommuting pair")
    else:
        check(len(ap) > 0 and (((2, 0), (0, 1)) in ap or ((0, 1), (2, 0)) in ap) and rk == 15,
              f"Weyl letters at d = 4: {len(ap)} anticommuting pairs (e.g. X^2, Z); their commutators span {rk} = all traceless")
# d = 4, Lean's basis: commutator span and anticommuting pairs
B4 = BAS[4]
cm4 = []
anti4 = 0
for i, j in itertools.combinations(range(1, 16), 2):
    P = rmul(B4[i], B4[j]); Q_ = rmul(B4[j], B4[i])
    cm4.append({(r, c): P[r][c] - Q_[r][c] for r in range(4) for c in range(4) if P[r][c] != Q_[r][c]})
    anti4 += all(P[r][c] + Q_[r][c] == 0 for r in range(4) for c in range(4))
check(len(cm4) == 105 and rank(cm4) == 15,
      f"d = 4: the 105 commutators of the 15 letters span {rank(cm4)} dims, so 105 independent commutator records do not exist")
print(f"      (Lean's tracelessBasis 4 has {anti4} anticommuting pairs; weylIsSign 4 = false in Lean gives no commutator records)")

head("B8. Collapse kappa(X) = X|Omega> (Prop. 8)")
check(all(not apply_omega(v) for v in C), "kappa(C) = 0")
check(rank([apply_omega(v) for v in Cp]) == 3 and rank([apply_omega(v) for v in CL]) == 1 and rank([apply_omega(v) for v in CLp]) == 2,
      "dim kappa(C+) = 3; plain frame: dim kappa(C_L) = 1, dim kappa(C_L+) = 2")
KP2 = Span([apply_omega(v) for v in P2])
check(all(KP2.contains(apply_omega(v)) for v in Cp) and rank([apply_omega(v) for v in PLE2V]) == 36,
      "kappa(C+) lies in kappa(P2); rank kappa on P<=2 = 36")

head("B9. F-CT-7: the construction at other (dQ, dq), by enumeration")


def fct7(dQ, dq, ncomm):
    a_, b_ = dQ * dQ - 1, dq * dq - 1
    dims = (dQ, dq, dQ, dq)
    n = {"odd": 0, "em": 0, "es": 0}
    seen = set()
    for s1, s2 in itertools.combinations(range(4), 2):
        for i in range(1, dims[s1] ** 2):
            for j in range(1, dims[s2] ** 2):
                t = [0, 0, 0, 0]; t[s1] = i; t[s2] = j
                t = tuple(t)
                if t in seen:
                    continue
                u = (t[2], t[3], t[0], t[1])
                same = (s1 % 2) == (s2 % 2)
                seen.add(t); seen.add(u)
                key = "es" if same else "em"
                if u == t:
                    n[key] += 1
                else:
                    n[key] += 1
                    n["odd"] += 1
    num = n["odd"] + n["em"] * sp.sqrt(dq) + n["es"] * sp.sqrt(dQ) + ncomm * sp.sqrt(dQ * dq)
    N = a_ * b_ * (dQ * dq) ** 4
    return n, num, N, sp.N(N / num, 12)


def approx(v, target):
    nd = len(target.split(".")[1])
    return round(float(v), nd) == float(target)


n32, num32, N32, v32 = fct7(3, 2, 3)
check(n32 == {"odd": 79, "em": 48, "es": 42} and N32 == 31104 and str(v32).startswith("137.035999"),
      f"(3,2): classes {n32}, N = {N32}, 1/value = {v32}")
for (dQ, dq, nc, want) in [(2, 2, 3, "31.81"), (2, 3, 0, "140.40"), (3, 3, 0, "791.66"), (4, 2, 3, "314.67"),
                           (2, 4, 0, "331.40"), (2, 4, WEYL_RANK, "307.91"), (2, 4, 105, "216.04"), (3, 4, 0, "2215.38"),
                           (4, 4, 0, "7228.24"), (5, 4, 0, "16551.02")]:
    n_, num_, N_, v_ = fct7(dQ, dq, nc)
    check(approx(v_, want), f"({dQ},{dq}) with {nc:3d} commutator records: numerator {sp.nsimplify(num_)}, N = {N_}, 1/value = {v_}")

# (OD) + (ID) of records-uniform.md applied at (dQ, dq) = (2,3): slot dimensions (2,3,2,3)
DIMS = (2, 3, 2, 3)
LET23 = [(s_, a_) for s_ in range(4) for a_ in range(1, DIMS[s_] ** 2)]
MS23, MA23 = [], []
for (s_, a_), (u_, b_) in itertools.product(LET23, repeat=2):
    xy = omul(single(s_, a_), single(u_, b_)); yx = omul(single(u_, b_), single(s_, a_))
    MS23.append(oadd(xy, yx)); MA23.append(oadd(xy, yx, -1))
Vc23 = [single(s_, a_) for s_, a_ in LET23 if s_ in (0, 2)]     # letters at the qubit centre copies
Vf23 = [single(s_, a_) for s_, a_ in LET23 if s_ in (1, 3)]     # letters at the qutrit far ends
SS23, SA23 = Span(MS23), Span(MA23)
check(inter_dim(MS23, MA23) == 16 and all(SS23.contains(v) and SA23.contains(v) for v in Vf23)
      and not any(SS23.contains(v) for v in Vc23),
      "(2,3): A = m(Sym^2) cap m(Lambda^2) = the far qutrit letters (16); the 6 qubit-centre letters are order-definite")
IE23 = eig_space(Eplain, [ident()] + Vc23 + Vf23, +1)
part23 = perp_in([ident()] + Vc23, IE23)                        # D cap I is 1 + centre letters
diffs23 = [oadd(single(0, a_), single(2, a_), -1) for a_ in range(1, 4)]
check(same_space(part23, diffs23) and rank(part23) == 3,
      "(2,3): (OD)+(ID) with the plain exchange keeps 3 single-letter records, the centre differences")
DIMS = (3, 2, 3, 2)

# (2,3) is (3,2) with sqrt2 <-> sqrt3 on the even classes; with three commutator records at its qubit centre
# (supplement.tex:310 'formed at the qubit end'; (OD)+(ID) at a qubit slot) it lands within 1% of 137.036
ALPHA = 137.035999177
n23, num23, N23, v23 = fct7(2, 3, 0)
check(n23 == n32 and N23 == N32, "(2,3) has the same classes (79, 48, 42) and N = 31104 as (3,2)")
for wname, w in [("sqrt2 (M_L at a qubit centre: Pi_X = 1, Pi_Y = 0)", sp.sqrt(2)), ("sqrt3", sp.sqrt(3)),
                 ("sqrt6 (the (3,2) commutator weight)", sp.sqrt(6))]:
    v = sp.N(N23 / (num23 + 3 * w), 12)
    rel = abs(float(v) - ALPHA) / ALPHA
    check(rel < 0.01, f"(2,3) with 3 commutator records at the qubit centre, weight {wname}: 1/value = {v}, {100 * rel:.2f}% from 1/alpha (< 1%)")
rel0 = abs(float(v23) - ALPHA) / ALPHA
check(abs(rel0 - 0.0246) < 0.0005, f"(2,3) with none (far-end placement, F-CT-7): 1/value = {v23}, {100 * rel0:.2f}% from 1/alpha")

print()
print(f"SUMMARY: {NCHK[0]} checks, {NCHK[1]} True, {NCHK[0] - NCHK[1]} False")
print(f"runtime: {time.perf_counter() - T0:.1f} s")
