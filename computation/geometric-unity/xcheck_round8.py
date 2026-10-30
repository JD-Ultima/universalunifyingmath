"""Referee cross-check of GU round 8 (formalization/geometric-unity/L-nonradial-z2.md), REVIEW6-GU.

Imports only xcheck_round6.py (and through it xcheck_round5.py), unmodified: the round-5/round-6 referee construction
(own Clifford algebra, the four Shiab templates from their definitions, the bilinear 2B, the vector-sector patterns
P_k and the grade extraction), used here ONLY to obtain the J1/J2 tables S, M, l, b and c1, c13 independently of
J_slaved_sector.py.  Nothing from L_nonradial.py, J_slaved_sector.py, I_elliptic_reduction.py, G_symbols.py or
F_jets.py is imported.  Written here from the definitions:
  - an own bookkeeping of the odd equation on the vector sector (Sigma = -M^-1 S; Hessian, gradient and Lambda^2
    coefficients), and a second, direct one along curves: for A_k = a_k(rho) d rho the traceless part of E^(i) is
    (sum_k l_k a_k) (Hess rho)_0 + (sum_k l_k a_k' + sum_{k<=l} w_kl (b1 + b2) a_k a_l) (d rho (x) d rho)_0;
  - own truncated-polynomial jets and an own textbook Ricci routine (Christoffel symbols and their first derivatives
    from a metric 2-jet), applied to Z2 at points OFF the section over a random polynomial base metric of degree 3
    in general coordinates;
  - the fibre connection Gamma(h, k) = -(1/2)(h y^-1 k + k y^-1 h), checked to be the Levi-Civita connection of
    c_V G_lambda (metric, symmetric) for lambda = 0, 1/2, 1; its curvature and Ricci tensor recomputed from Gamma;
  - the referee's own derivation of the order-1 integrability condition (C1) of Lemma L7 (REVIEW6-GU, section 2), its
    A-term from the definition of (AU, AW) by exact differentiation, and the pairings over own dual bases of p;
  - an own brute-force Taylor solution of the fibre equation Hess_F w = mu(w)(Ric_F + s (A., A.)) at y = eta to order
    eps^2 (1000 unknowns), with the control scale s* derived by hand from (C1).
Exact arithmetic only (Fractions, Python/numpy integers, sympy rationals and Gaussian rationals); fixed seeds; one core.
Every check is asserted and printed.  Sections X0-X6 (X6a-c) are cited as xcheck section Xk in the report.
Run (from computation/geometric-unity): python3 xcheck_round8.py > xcheck_round8.out 2>&1   (~2 min, 70 checks)
A section can be run alone for testing: python3 xcheck_round8.py x2 x6c (x1, x4, x5 also run X0).
"""
import contextlib
import io
import itertools
import os
import random
import sys
import time
from fractions import Fraction as Fr

for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import numpy as np  # noqa: E402
import sympy as sp  # noqa: E402
from sympy.polys.matrices import DomainMatrix  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

T0 = time.time()
RES = []
I = sp.I
x, y, z, w, kap, t, u = sp.symbols('x y z w kappa t u')
sig, cc = sp.symbols('sigma c')
KS = (0, 2, 12, 14)                      # A, C, D, B
KN = {0: 'A', 2: 'C', 12: 'D', 14: 'B'}
LAB = {1: '(14,0)', -1: '(7,7)'}


def say(*a):
    print(*a)
    sys.stdout.flush()


def check(label, ok):
    ok = bool(ok)
    RES.append((label, ok))
    say(f"  [{len(RES):3d}] {label}: {ok}")
    assert ok, label


def section(s):
    say(f"\n=== {s}")


def tick(t1):
    return f"[{time.time() - t1:.0f} s]"


def zero(e):
    return sp.expand(e) == 0


def gpoly(e, gens):
    return sp.Poly(sp.expand(e), *gens, domain='QQ_I')


# =============================================================================================== X0 inputs
def x0():
    section("X0  Inputs: S, M, l, b, c1, c13 from the round-6 referee construction (xcheck_round6 X0-X2, own Clifford "
            "algebra), not from J_slaved_sector.py or L_nonradial.py")
    t1 = time.time()
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        import xcheck_round6 as R6
        tabs = {}
        for eta in (R6.R5.EUC, R6.R5.E77):
            C = R6.Ctx(eta)
            c1, c13 = R6.x0(C)
            d1 = R6.x1(C)
            L, B = R6.x2(C)
            tabs[C.eps] = dict(c1=c1, c13=c13, S=d1['S'], M=d1['M'], L=L, B=B)
    labs = [lab for lab, _ in R6.RES]
    check(f"xcheck_round6 X0-X2 rerun as a library, both signatures: {len(R6.RES)} internal checks, all True (J1 fits on "
          f"all 14 covectors, J2 on all 196 tensors and sampled pairs; tables equal J's printed tables) {tick(t1)}",
          len(R6.RES) == 38 and all(ok for _, ok in R6.RES) and
          all(sum((lb + ' ' + key) in s_ for s_ in labs) == 1 for lb in ('(14,0)', '(7,7)')
              for key in ('J1(a) S equals', 'J1(b) M equals', 'J2 derivative table equals', 'J2 quadratic table equals')))
    for ep in (1, -1):
        T = tabs[ep]
        say(f"     {LAB[ep]}: S = {T['S'].tolist()}")
        say(f"     {LAB[ep]}: c1 = {T['c1']};  c13 = {T['c13']}")
    return tabs


# =============================================================================================== own bookkeeping
def lrow(T, nm):
    return [T['L'][(nm, k)] for k in KS]


def assemble(T, Sig, dSt, dSu):
    """own bookkeeping: Hessian H[i][phi], symmetric gradient coefficients Qs[i][(phi,psi)] (of dphi (x) dpsi, symmetrised),
    trace coefficients Qp[i][(phi,psi)] (of <dphi, dpsi> g), Lambda^2 coefficients W[nm] (of dt ^ du = dt(x)du - du(x)dt
    in the antisymmetric part of E1, E13, and in w3, w11).  Sigma[k][phi], dSt = d_t Sigma, dSu = d_u Sigma."""
    dS = (dSt, dSu)
    H, Qs, Qp, W = {}, {}, {}, {}
    for nm in ('E1', 'E13'):
        l = lrow(T, nm)
        H[nm] = [sum(l[k] * Sig[k, a] for k in range(4)) for a in range(2)]
        Qf = {}
        for a in range(2):
            for b in range(2):
                v = sum(l[k] * dS[b][k, a] for k in range(4))       # l_k d_psi Sigma_k,phi  (N^T)
                vp = -v
                for (ik, k), (il, ll) in itertools.combinations_with_replacement(list(enumerate(KS)), 2):
                    b1, b2, b3 = T['B'][(nm, k, ll)]
                    wkl = sp.Rational(1, 2) if k == ll else 1
                    v += wkl * (b1 * Sig[ik, a] * Sig[il, b] + b2 * Sig[il, a] * Sig[ik, b])
                    vp += wkl * b3 * Sig[ik, a] * Sig[il, b]
                Qf[(a, b)] = v
                Qp.setdefault(nm, {})[(a, b)] = vp
        Qs[nm] = {(a, b): (Qf[(a, b)] + Qf[(b, a)]) / 2 for a in range(2) for b in range(2)}
        W[nm] = (Qf[(0, 1)] - Qf[(1, 0)]) / 2
    for nm in ('w3', 'w11'):
        l = lrow(T, nm)
        v = sum(l[k] * (dSt[k, 1] - dSu[k, 0]) for k in range(4))
        for (ik, k), (il, ll) in itertools.combinations(list(enumerate(KS)), 2):
            v += T['B'][(nm, k, ll)] * (Sig[ik, 0] * Sig[il, 1] - Sig[ik, 1] * Sig[il, 0])
        W[nm] = v
    return H, Qs, Qp, W


def slaving(T, sub):
    S = T['S'].subs(sub)
    M = T['M'].subs(sub)
    D = sp.expand(M.det())
    NS = (-M.adjugate() * S).applyfunc(sp.expand)
    return S, M, D, NS


# =============================================================================================== X1 Lemma L1
def x1(tabs):
    section("X1  Lemma L1: the Lambda^2 coefficients varpi (own assembly from the X0 tables)")
    t1 = time.time()
    for ep in (1, -1):
        T = tabs[ep]
        lab = LAB[ep]
        # x W_A, exact, x symbolic
        sub = {y: 0, z: 0, w: 0}
        S, M, D, NS = slaving(T, sub)
        Sig = (NS / D).applyfunc(sp.cancel)
        H, Qs, Qp, W = assemble({**T, 'L': {k: v.subs(sub) for k, v in T['L'].items()},
                                 'B': {k: (tuple(a.subs(sub) for a in v) if isinstance(v, tuple) else v.subs(sub))
                                       for k, v in T['B'].items()}},
                                Sig, Sig.diff(t), Sig.diff(u))
        dA = kap ** 2 + 224 * kap * t * x - 3328 * ep * u ** 2 * x ** 2
        okW = sp.cancel(W['E1'] - ep * 1916928 * kap * u * x ** 4 / dA ** 2) == 0 and \
            all(sp.cancel(W[nm]) == 0 for nm in ('E13', 'w3', 'w11'))
        check(f"{lab} x W_A: varpi(E1) = {'+' if ep == 1 else '-'}1916928 kappa u x^4/(kappa^2 + 224 kappa t x "
              f"{'-' if ep == 1 else '+'} 3328 u^2 x^2)^2, the other three vanish (nonzero rational function for x != 0)", okW)
        # x W_A + z W_C AdG: varpi(E13) = x z^2 P / (det M)^2
        sub = {y: 0, w: 0}
        Tz = {**T, 'L': {k: v.subs(sub) for k, v in T['L'].items()},
              'B': {k: (tuple(a.subs(sub) for a in v) if isinstance(v, tuple) else v.subs(sub)) for k, v in T['B'].items()}}
        S, M, D, NS = slaving(T, sub)
        Sig = NS / D
        H, Qs, Qp, W = assemble(Tz, Sig, Sig.diff(t), Sig.diff(u))
        P13 = sp.cancel(W['E13'] * D ** 2 / (x * z ** 2))
        n_, d_ = sp.fraction(sp.together(P13))
        check(f"{lab} x W_A + z W_C AdG: varpi(E13) (det M)^2 = x z^2 P(x, z, kappa, t, u) with P a nonzero polynomial, "
              f"P(kappa, 0, 0) = 5 kappa^4 * 9216 kappa^2", sp.Poly(d_, x, z, kap, t, u).is_ground and n_ != 0 and
              zero(P13.subs({t: 0, u: 0}) - 9216 * 5 * kap ** 6))
        # every member, leading order kappa^-2 of varpi(w3), varpi(w11)  (Neumann series of M^-1)
        S = T['S']
        M1 = T['M'] - kap * sp.eye(4)
        e_ = sp.Symbol('e_')
        Sig2 = -S * e_ + M1 * S * e_ ** 2                       # e_ = 1/kappa
        H, Qs, Qp, W = assemble(T, Sig2, Sig2.diff(t), Sig2.diff(u))
        ell = x * z + ep * y * w
        lead = {}
        for nm in ('w3', 'w11', 'E1', 'E13'):
            Pw = sp.Poly(sp.expand(W[nm]), e_)
            lead[nm] = sp.factor(Pw.coeff_monomial(e_ ** 2))
        okL = zero(lead['w3'] - (-ep) * 88320 * I * y * ell) and zero(lead['w11'] - 88320 * I * x * ell) and \
            zero(lead['E1'] - (-46080 * w * ell if ep == 1 else 46080 * w * ell)) and zero(lead['E13'] - 46080 * z * ell)
        check(f"{lab} every member, coefficient of kappa^-2: varpi(w3) = {'-' if ep == 1 else '+'}88320 i y l, varpi(w11) "
              f"= 88320 i x l, varpi(E1) = {'-' if ep == 1 else '+'}46080 w l, varpi(E13) = 46080 z l (l = xz "
              f"{'+' if ep == 1 else '-'} yw): some varpi is nonzero whenever l != 0", okL)
        # observation: members with l = 0 outside L1's list (True = the referee's statement is confirmed)
        obs = []
        if ep == 1:
            mems = [((0, 1, 0, 0), "W'_A", True), ((1, 0, 0, 1), "W_A + W'_C AdG", True),
                    ((0, 1, 1, 0), "W'_A + W_C AdG", True), ((1, 1, 0, 0), "W_A + W'_A", True),
                    ((0, 0, 1, 0), "W_C AdG", False), ((0, 0, 0, 1), "W'_C AdG", False)]
        else:
            mems = [((0, 1, 0, 0), "W'_A", True), ((1, 0, 0, 1), "W_A + W'_C AdG", True),
                    ((0, 1, 1, 0), "W'_A + W_C AdG", True), ((2, 1, 1, 2), "(2, 1, 1, 2)", True),
                    ((1, 1, 0, 0), "W_A + W'_A", False), ((1, 1, 1, 1), "(1, 1, 1, 1)", False),
                    ((0, 0, 1, 0), "W_C AdG", False), ((0, 0, 0, 1), "W'_C AdG", False)]
        okO = True
        for mem, nmem, expect in mems:
            sub = dict(zip((x, y, z, w), mem))
            Tm = {**T, 'L': {k: v.subs(sub) for k, v in T['L'].items()},
                  'B': {k: (tuple(a.subs(sub) for a in v) if isinstance(v, tuple) else v.subs(sub))
                        for k, v in T['B'].items()}}
            S_, M_, D_, NS_ = slaving(T, sub)
            Sg = NS_ / D_
            _, _, _, Wm = assemble(Tm, Sg, Sg.diff(t), Sg.diff(u))
            nz = [nm for nm in ('E1', 'E13', 'w3', 'w11') if sp.cancel(Wm[nm]) != 0]
            okO &= mem[0] * mem[2] + ep * mem[1] * mem[3] == 0 and (len(nz) > 0) == expect
            obs.append(f"{nmem}: {nz if nz else 'all four vanish identically'}")
        check(f"{lab} observation (not in the note), members with l = 0 outside L1's list, exact varpi: " + "; ".join(obs) +
              ". So L1's Lambda^2 step extends to the first ones, but NOT to the last ones, for which dt ^ du = 0 on open "
              "solution sets is not established", okO)
    # Shiab#: exact varpi(E13) nonzero, numerator 9216 kappa^2 (5 kappa^4 - 1600 kappa^3 t - ...)
    for ep in (1, -1):
        T = tabs[ep]
        sub = {x: 1, y: 0, z: 1, w: 0}
        Tm = {**T, 'L': {k: v.subs(sub) for k, v in T['L'].items()},
              'B': {k: (tuple(a.subs(sub) for a in v) if isinstance(v, tuple) else v.subs(sub)) for k, v in T['B'].items()}}
        S_, M_, D_, NS_ = slaving(T, sub)
        Sg = NS_ / D_
        _, _, _, Wm = assemble(Tm, Sg, Sg.diff(t), Sg.diff(u))
        num = sp.expand(sp.cancel(Wm['E13'] * D_ ** 2))
        Pn = sp.Poly(num, kap, t, u)
        check(f"{LAB[ep]} Shiab#: varpi(E13) (det M)^2 = 9216 kappa^2 (5 kappa^4 - 1600 kappa^3 t + ...), a nonzero "
              f"polynomial (J Cor. J2.1(e))", Pn.coeff_monomial(kap ** 6) == 9216 * 5 and
              Pn.coeff_monomial(kap ** 5 * t) == -9216 * 1600)
    say(f"     {tick(t1)}")


# =============================================================================================== jets
class TP:
    """truncated polynomials in nv variables of total degree <= deg: dict monomial-index -> Fraction"""

    def __init__(self, nv, deg):
        self.nv, self.deg = nv, deg
        self.mons = []
        for d in range(deg + 1):
            for comb in itertools.combinations_with_replacement(range(nv), d):
                e = [0] * nv
                for c in comb:
                    e[c] += 1
                self.mons.append(tuple(e))
        self.ix = {m: i for i, m in enumerate(self.mons)}
        self.dg = [sum(m) for m in self.mons]
        self.mt = {}
        for i, a in enumerate(self.mons):
            for j, b in enumerate(self.mons):
                if self.dg[i] + self.dg[j] <= deg:
                    self.mt[(i, j)] = self.ix[tuple(p + q for p, q in zip(a, b))]

    def c(self, v):
        return {0: Fr(v)} if v != 0 else {}

    def var(self, k):
        e = [0] * self.nv
        e[k] = 1
        return {self.ix[tuple(e)]: Fr(1)}

    def add(self, a, b, s=1):
        o = dict(a)
        for k, v in b.items():
            nv_ = o.get(k, 0) + s * v
            if nv_:
                o[k] = nv_
            else:
                o.pop(k, None)
        return o

    def scl(self, a, s):
        return {k: v * s for k, v in a.items()} if s else {}

    def mul(self, a, b):
        o = {}
        mt = self.mt
        for i, va in a.items():
            for j, vb in b.items():
                k = mt.get((i, j))
                if k is not None:
                    o[k] = o.get(k, 0) + va * vb
        return {k: v for k, v in o.items() if v}

    def ssum(self, lst):
        o = {}
        for p in lst:
            o = self.add(o, p)
        return o

    def deriv(self, a, i):
        o = {}
        for k, v in a.items():
            m = self.mons[k]
            if m[i]:
                mm = list(m)
                mm[i] -= 1
                o[self.ix[tuple(mm)]] = o.get(self.ix[tuple(mm)], 0) + v * m[i]
        return o

    def embed(self, a, src, pos):
        """polynomial a of TP src -> self, variable j of src -> variable pos[j] of self (degree truncated)"""
        o = {}
        for k, v in a.items():
            m = src.mons[k]
            if sum(m) > self.deg:
                continue
            e = [0] * self.nv
            for j, p in enumerate(m):
                e[pos[j]] += p
            o[self.ix[tuple(e)]] = o.get(self.ix[tuple(e)], 0) + v
        return o


def pm_mul(T, A, B):
    n, m, r = len(A), len(B), len(B[0])
    return [[T.ssum([T.mul(A[i][k], B[k][j]) for k in range(m) if A[i][k] and B[k][j]]) for j in range(r)]
            for i in range(n)]


def pm_add(T, A, B, s=1):
    return [[T.add(A[i][j], B[i][j], s) for j in range(len(A[0]))] for i in range(len(A))]


def pm_scl(T, A, s):
    return [[T.scl(A[i][j], s) for j in range(len(A[0]))] for i in range(len(A))]


def pm_tr(T, A):
    return T.ssum([A[i][i] for i in range(len(A))])


def pm_trprod(T, A, B):
    n = len(A)
    return T.ssum([T.mul(A[i][j], B[j][i]) for i in range(n) for j in range(n) if A[i][j] and B[j][i]])


def pm_const(T, A):
    return [[T.c(A[i][j]) for j in range(len(A[0]))] for i in range(len(A))]


def finv(Mx):
    """exact inverse of a Fraction matrix (Gauss-Jordan)"""
    n = len(Mx)
    A = [[Fr(Mx[i][j]) for j in range(n)] + [Fr(int(i == j)) for j in range(n)] for i in range(n)]
    for c in range(n):
        p = next(r for r in range(c, n) if A[r][c] != 0)
        A[c], A[p] = A[p], A[c]
        pv = A[c][c]
        A[c] = [v / pv for v in A[c]]
        for r in range(n):
            if r != c and A[r][c] != 0:
                f = A[r][c]
                A[r] = [a - f * b for a, b in zip(A[r], A[c])]
    return [row[n:] for row in A]


def fmul(A, B):
    return [[sum(A[i][k] * B[k][j] for k in range(len(B))) for j in range(len(B[0]))] for i in range(len(A))]


def ftr(A):
    return sum(A[i][i] for i in range(len(A)))


def fT(A):
    return [list(r) for r in zip(*A)]


def fadd(A, B, s=1):
    return [[A[i][j] + s * B[i][j] for j in range(len(A[0]))] for i in range(len(A))]


def fscl(A, s):
    return [[s * v for v in r] for r in A]


def fbr(A, B):
    return fadd(fmul(A, B), fmul(B, A), -1)


def ricci_from_jet(T, g):
    """Ricci tensor at the origin of the metric with 2-jet g[a][b] (TP polynomials, deg 2): own textbook routine.
    R^k_lij = d_i Gam^k_jl - d_j Gam^k_il + Gam^k_im Gam^m_jl - Gam^k_jm Gam^m_il,  Ric_lj = R^k_lkj."""
    n = len(g)
    g0 = [[g[a][b].get(0, Fr(0)) for b in range(n)] for a in range(n)]
    gi = finv(g0)
    lin = [T.ix[tuple(int(i == m) for i in range(T.nv))] for m in range(n)]
    dg = [[[g[a][b].get(lin[m], Fr(0)) for b in range(n)] for a in range(n)] for m in range(n)]

    def ddg(m, p, a, b):
        e = [0] * T.nv
        e[m] += 1
        e[p] += 1
        v = g[a][b].get(T.ix[tuple(e)], Fr(0))
        return 2 * v if m == p else v
    # Gamma0^k_ij
    Tl = [[[dg[i][j][l] + dg[j][i][l] - dg[l][i][j] for j in range(n)] for i in range(n)] for l in range(n)]
    G0 = [[[sum(gi[k][l] * Tl[l][i][j] for l in range(n)) / 2 for j in range(n)] for i in range(n)] for k in range(n)]
    dgi = [[[-sum(gi[k][a] * dg[m][a][b] * gi[b][l] for a in range(n) for b in range(n) if dg[m][a][b])
             for l in range(n)] for k in range(n)] for m in range(n)]

    def dGam(m, k, i, j):
        s1 = sum(dgi[m][k][l] * Tl[l][i][j] for l in range(n))
        s2 = sum(gi[k][l] * (ddg(m, i, j, l) + ddg(m, j, i, l) - ddg(m, l, i, j)) for l in range(n) if gi[k][l])
        return (s1 + s2) / 2
    Ric = [[Fr(0)] * n for _ in range(n)]
    for l in range(n):
        for j in range(l, n):
            v = sum(dGam(k, k, j, l) for k in range(n)) - sum(dGam(j, k, k, l) for k in range(n))
            v += sum(G0[k][k][m] * G0[m][j][l] - G0[k][j][m] * G0[m][k][l] for k in range(n) for m in range(n))
            Ric[l][j] = Ric[j][l] = v
    return Ric, G0, gi, g0


def riemann_from_jet(T, g):
    """R^k_lij at the origin (base metric, small n)"""
    n = len(g)
    g0 = [[g[a][b].get(0, Fr(0)) for b in range(n)] for a in range(n)]
    gi = finv(g0)
    lin = [T.ix[tuple(int(i == m) for i in range(T.nv))] for m in range(n)]
    dg = [[[g[a][b].get(lin[m], Fr(0)) for b in range(n)] for a in range(n)] for m in range(n)]

    def ddg(m, p, a, b):
        e = [0] * T.nv
        e[m] += 1
        e[p] += 1
        v = g[a][b].get(T.ix[tuple(e)], Fr(0))
        return 2 * v if m == p else v
    Tl = [[[dg[i][j][l] + dg[j][i][l] - dg[l][i][j] for j in range(n)] for i in range(n)] for l in range(n)]
    G0 = [[[sum(gi[k][l] * Tl[l][i][j] for l in range(n)) / 2 for j in range(n)] for i in range(n)] for k in range(n)]
    dgi = [[[-sum(gi[k][a] * dg[m][a][b] * gi[b][l] for a in range(n) for b in range(n)) for l in range(n)]
            for k in range(n)] for m in range(n)]
    dG = [[[[(sum(dgi[m][k][l] * Tl[l][i][j] for l in range(n)) +
              sum(gi[k][l] * (ddg(m, i, j, l) + ddg(m, j, i, l) - ddg(m, l, i, j)) for l in range(n))) / 2
             for j in range(n)] for i in range(n)] for k in range(n)] for m in range(n)]
    R = [[[[dG[i][k][j][l] - dG[j][k][i][l] + sum(G0[k][i][m] * G0[m][j][l] - G0[k][j][m] * G0[m][i][l]
                                                     for m in range(n))
            for j in range(n)] for i in range(n)] for l in range(n)] for k in range(n)]
    return R


# =============================================================================================== X2 O'Neill off s(X)
SYMIDX = [(i, j) for i in range(4) for j in range(i, 4)]


def Emat(A):
    i, j = SYMIDX[A]
    Mx = [[Fr(0)] * 4 for _ in range(4)]
    Mx[i][j] = Fr(1)
    Mx[j][i] = Fr(1)
    return Mx


def coords(Mx):
    return [Mx[i][j] for (i, j) in SYMIDX]


def x2():
    section("X2  Lemma L2 inputs OFF the section: brute-force Ricci of Z2 (own jets, own Ricci routine) at (0, y0), "
            "y0 != g(0), random polynomial base metric of degree 3 in general coordinates")
    rng = random.Random(20261108)
    # control: round 3-sphere, g = 4/(1+|x|^2)^2 delta, 2-jet 4 delta - 8 |x|^2 delta at 0, Ric = 2 g
    T3 = TP(3, 2)
    g = [[T3.add(T3.c(4), T3.ssum([T3.scl(T3.mul(T3.var(m), T3.var(m)), -8) for m in range(3)])) if a == b else {}
          for b in range(3)] for a in range(3)]
    Ric, _, _, g0 = ricci_from_jet(T3, g)
    check("control: own Ricci routine on the round 3-sphere (2-jet at 0) gives Ric = 2 g", all(
        Ric[a][b] == 2 * g0[a][b] for a in range(3) for b in range(3)))
    cases = [((1, 1, 1, 1), Fr(0), Fr(3), Fr(2)), ((1, -1, -1, -1), Fr(1, 2), Fr(3), Fr(2)),
             ((-1, 1, 1, 1), Fr(1), Fr(-3), Fr(2)), ((1, 1, -1, -1), Fr(1, 2), Fr(3), Fr(-2))]
    for eta, lam, cH, cV in cases:
        t1 = time.time()
        sigl = f"({eta.count(1)},{eta.count(-1)})"
        Tb = TP(4, 3)
        # random symmetric polynomial base metric g(x) = eta + deg 1..3 terms (general coordinates)
        gb = [[None] * 4 for _ in range(4)]
        for a in range(4):
            for b in range(a, 4):
                p = Tb.c(eta[a]) if a == b else {}
                for k in range(1, len(Tb.mons)):
                    if rng.random() < 0.35:
                        p = Tb.add(p, {k: Fr(rng.randint(-3, 3), rng.choice((5, 7, 10)))})
                gb[a][b] = gb[b][a] = p
        # inverse to degree 2: eta - eta P eta + eta P eta P eta
        etaM = [[Fr(eta[i]) if i == j else Fr(0) for j in range(4)] for i in range(4)]
        Pb = [[Tb.add(gb[a][b], Tb.c(eta[a]) if a == b else {}, -1) for b in range(4)] for a in range(4)]
        Ec = pm_const(Tb, etaM)
        EP = pm_mul(Tb, Ec, Pb)
        EPE = pm_mul(Tb, EP, Ec)
        gbi = pm_add(Tb, pm_add(Tb, Ec, EPE, -1), pm_mul(Tb, EP, EPE))
        gbi = [[{k: v for k, v in e.items() if Tb.dg[k] <= 2} for e in r] for r in gbi]
        # Gamma^k_{mu i}(x) to degree 2
        dgb = [[[Tb.deriv(gb[a][b], m) for b in range(4)] for a in range(4)] for m in range(4)]
        Gam = [[[{kk: v for kk, v in Tb.scl(Tb.ssum([Tb.mul(gbi[k][l], Tb.add(Tb.add(dgb[m][i][l], dgb[i][m][l]),
                                                                                 dgb[l][m][i], -1)) for l in range(4)]),
                                              Fr(1, 2)).items() if Tb.dg[kk] <= 2}
                 for i in range(4)] for m in range(4)] for k in range(4)]          # Gam[k][m][i]
        # base curvature at 0
        T4 = TP(4, 2)
        gb2 = [[T4.embed(gb[a][b], Tb, list(range(4))) for b in range(4)] for a in range(4)]
        Rb = riemann_from_jet(T4, gb2)                 # Rb[k][l][i][j] = R^k_{l i j}
        Ricb = [[sum(Rb[k][l][k][j] for k in range(4)) for j in range(4)] for l in range(4)]
        scal = sum(Fr(eta[l]) * Ricb[l][l] for l in range(4))
        Bend = {(i, j): [[Rb[k][l][i][j] for l in range(4)] for k in range(4)] for i in range(4) for j in range(4)}
        okso = all(all(fadd(fmul(fT(Bm), etaM), fmul(etaM, Bm))[a][b] == 0 for a in range(4) for b in range(4))
                   for Bm in Bend.values())
        # y0 = P^T eta P, P unimodular integer (same signature, off the section)
        while True:
            Pm = [[Fr(int(i == j)) for j in range(4)] for i in range(4)]
            for _ in range(5):
                i, j = rng.sample(range(4), 2)
                cf = rng.choice((-1, 1, 2))
                Pm = [[Pm[r][c] + (cf * Pm[j][c] if r == i else 0) for c in range(4)] for r in range(4)]
            y0 = fmul(fmul(fT(Pm), etaM), Pm)
            if y0 != etaM:
                break
        # Z2 metric 2-jet in 14 variables (x^0..x^3, Y_0..Y_9) at (0, y0)
        T14 = TP(14, 2)
        posx = [0, 1, 2, 3]
        GamE = [[[T14.embed(Gam[k][m][i], Tb, posx) for i in range(4)] for k in range(4)] for m in range(4)]  # [m][k][i]
        Yv = [[T14.add(T14.c(y0[i][j]), T14.var(4 + SYMIDX.index((min(i, j), max(i, j))))) for j in range(4)]
              for i in range(4)]
        y0i = finv(y0)
        Yc = pm_const(T14, y0i)
        Dl = [[T14.var(4 + SYMIDX.index((min(i, j), max(i, j)))) for j in range(4)] for i in range(4)]
        A1 = pm_mul(T14, Yc, pm_mul(T14, Dl, Yc))
        A2 = pm_mul(T14, A1, pm_mul(T14, Dl, Yc))
        yi = pm_add(T14, pm_add(T14, Yc, A1, -1), A2)        # y^-1 to degree 2
        K = []
        for m in range(4):
            Gm_ = GamE[m]                                      # (Gamma_mu)[k][i]
            GmT = [[Gm_[i][k] for i in range(4)] for k in range(4)]
            K.append(pm_add(T14, pm_mul(T14, GmT, Yv), pm_mul(T14, Yv, Gm_)))
        Pmu = [pm_mul(T14, yi, K[m]) for m in range(4)]
        QA = [pm_mul(T14, yi, pm_const(T14, Emat(A))) for A in range(10)]

        def Gl(P1, P2):
            return T14.add(pm_trprod(T14, P1, P2), T14.mul(pm_tr(T14, P1), pm_tr(T14, P2)), -lam)
        gz = [[None] * 14 for _ in range(14)]
        for m in range(4):
            for n_ in range(m, 4):
                gz[m][n_] = gz[n_][m] = T14.add(T14.scl(T14.embed(gb[m][n_], Tb, posx), cH), T14.scl(Gl(Pmu[m], Pmu[n_]), cV))
            for A in range(10):
                gz[m][4 + A] = gz[4 + A][m] = T14.scl(Gl(QA[A], Pmu[m]), -cV)
        for A in range(10):
            for B in range(A, 10):
                gz[4 + A][4 + B] = gz[4 + B][4 + A] = T14.scl(Gl(QA[A], QA[B]), cV)
        Ric, _, _, g0 = ricci_from_jet(T14, gz)

        # ---- O'Neill predictions
        def G0f(h, k, l_=lam):
            return ftr(fmul(fmul(y0i, h), fmul(y0i, k))) - l_ * ftr(fmul(y0i, h)) * ftr(fmul(y0i, k))

        def Rdot(Bm, yy):
            return fscl(fadd(fmul(fT(Bm), yy), fmul(yy, Bm)), -1)
        RY = {ij: Rdot(Bm, y0) for ij, Bm in Bend.items()}
        Am = [[cV ** 2 / (4 * cH ** 2) * sum(eta[i] * eta[j] * G0f(Emat(A), RY[(i, j)]) * G0f(Emat(B), RY[(i, j)])
                                             for i in range(4) for j in range(4)) for B in range(10)] for A in range(10)]
        RicF = [[-G0f(Emat(A), Emat(B), Fr(1, 4)) for B in range(10)] for A in range(10)]
        okVV = all(Ric[4 + A][4 + B] == RicF[A][B] + Am[A][B] for A in range(10) for B in range(10))
        Kc = [coords([[K[m][i][j].get(0, Fr(0)) for j in range(4)] for i in range(4)]) for m in range(4)]
        Xh = [[Fr(int(a == m)) for a in range(4)] + Kc[m] for m in range(4)]

        def bil(Q, X1, X2):
            return sum(X1[a] * Q[a][b] * X2[b] for a in range(14) for b in range(14) if X1[a] and X2[b])
        Phi = [[sum(eta[i] * ftr(fmul(fmul(y0i, Rdot(Bend[(m, i)], y0)), fmul(y0i, Rdot(Bend[(n_, i)], y0))))
                    for i in range(4)) for n_ in range(4)] for m in range(4)]
        okHH = all(bil(Ric, Xh[m], Xh[n_]) == Ricb[m][n_] - cV / (2 * cH) * Phi[m][n_] for m in range(4) for n_ in range(4))
        okHm = all(bil(g0, Xh[m], Xh[n_]) == cH * (eta[m] if m == n_ else 0) for m in range(4) for n_ in range(4))
        check(f"base {sigl}, lambda = {lam}, c_H = {cH}, c_V = {cV}, y0 = P^T eta P != eta: B_ij = R_g(e_i, e_j) in "
              f"so(eta) ({okso}); VV block of Ric(Z2) = Ric_F + (AU, AW), (AU, AW) = c_V^2/(4 c_H^2) sum eps_i eps_j "
              f"G(U, R_ij.y) G(W, R_ij.y), Ric_F = -G_1/4 (100 components)", okso and okVV)
        check(f"base {sigl}: horizontal lifts X_mu = d_mu + K_mu(0, y0) are g_aleph-orthogonal to V with norm c_H g, and "
              f"HH block = Ric_g - (c_V/(2 c_H)) Phi_y0 (16 components)", okHm and okHH and
              all(bil(g0, Xh[m], [Fr(0)] * 4 + [Fr(int(b == A)) for b in range(10)]) == 0 for m in range(4) for A in range(10)))
        Ev = [Fr(0)] * 4 + coords(y0)
        okE = all(sum(Ric[a][b] * Ev[a] for a in range(14)) == 0 for b in range(14)) and \
            bil(g0, Ev, Ev) == 4 * cV * (1 - 4 * lam)
        okAE = all(sum(Am[A][B] * Ev[4 + A] for A in range(10)) == 0 for B in range(10))
        check(f"base {sigl}: Euler field E = y d_y: Ric(E, .) = 0 (all 14 components), g(E, E) = 4 c_V (1 - 4 lambda) = "
              f"{4 * cV * (1 - 4 * lam)}, and (AE, A.) = 0", okE and okAE)
        gV = [[g0[4 + A][4 + B] for B in range(10)] for A in range(10)]
        gVi = finv(gV)
        trV = sum(gVi[A][B] * Ric[4 + B][4 + A] for A in range(10) for B in range(10))
        A2n = cV / (4 * cH ** 2) * sum(eta[i] * eta[j] * G0f(RY[(i, j)], RY[(i, j)]) for i in range(4) for j in range(4))
        trVA = sum(gVi[A][B] * Am[B][A] for A in range(10) for B in range(10))
        trH = sum(Fr(eta[m]) / cH * bil(Ric, Xh[m], Xh[m]) for m in range(4))
        check(f"base {sigl}: L2 steps 1-2 off the section: tr_V Ric(Z2) = -9/c_V + |A|^2 and tr_H Ric(Z2) = scal_g/c_H - "
              f"2 |A|^2, with |A|^2 = c_V/(4 c_H^2) sum eps_i eps_j G(R_ij.y, R_ij.y) = tr_V (A., A.) = {A2n} != 0 "
              f"{tick(t1)}", trV == -9 / cV + A2n and trH == scal / cH - 2 * A2n and trVA == A2n and A2n != 0)


# =============================================================================================== X3 L2, L3 algebra
def x3():
    section("X3  Lemma L2 bookkeeping; Lemma L3 (commutant, [p, p]); L7 order 0 (commuting pair)")
    al, A2, sc, ch, cv = sp.symbols('alpha A2 scal c_H c_V')
    sol = sp.solve([9 * al - (-9 / cv + A2), 4 * al - (sc / ch - 2 * A2)], [al, A2], dict=True)[0]
    okal = zero(sol[al] - sp.Rational(1, 11) * (sc / (2 * ch) - 9 / cv))
    s0 = sp.solve(sp.Eq(sol[A2], 0), sc)[0]
    check("L2: 9 alpha = -9/c_V + |A|^2, 4 alpha = scal/c_H - 2|A|^2 give alpha = (scal/(2 c_H) - 9/c_V)/11; |A|^2 = 0 "
          "gives scal = -4 c_H/c_V and alpha = -1/c_V", okal and zero(s0 + 4 * ch / cv) and
          zero(sol[al].subs(sc, s0) + 1 / cv))
    for eta in ((1, 1, 1, 1), (1, -1, -1, -1), (-1, 1, 1, 1), (1, 1, -1, -1)):
        sigl = f"({eta.count(1)},{eta.count(-1)})"
        E = sp.diag(*eta)
        pb = [sp.diag(*[1 if k == i else (-1 if k == i + 1 else 0) for k in range(4)]) for i in range(3)]
        for i in range(4):
            for j in range(i + 1, 4):
                Mx = sp.zeros(4)
                Mx[i, j], Mx[j, i] = eta[i], eta[j]
                pb.append(Mx)
        kb = []
        for i in range(4):
            for j in range(i + 1, 4):
                Mx = sp.zeros(4)
                Mx[i, j], Mx[j, i] = eta[i], -eta[j]
                kb.append(Mx)
        okp = all((E * X).T == E * X and X.trace() == 0 for X in pb) and all((E * K).T == -(E * K) for K in kb)
        brs = [X * Y - Y * X for X in pb for Y in pb]
        rk = sp.Matrix([list(B) for B in brs]).rank()
        rkk = sp.Matrix([list(B) for B in brs + kb]).rank()
        # commutant of so(eta) in gl(4)
        cs = sp.symbols('m0:16')
        Mg = sp.Matrix(4, 4, cs)
        eqs = []
        for K in kb:
            eqs += list(Mg * K - K * Mg)
        solc = sp.linsolve(eqs, cs)
        (vec,) = solc
        free = set().union(*[sp.sympify(v).free_symbols for v in vec])
        okc = len(free) == 1 and sp.Matrix(4, 4, list(vec)) == list(free)[0] * sp.eye(4)
        # {v in p : [[X, Y], v] = 0 for all X, Y in p} = 0
        vs = sp.symbols('v0:9')
        V = sum((vs[i] * pb[i] for i in range(9)), sp.zeros(4))
        eq2 = []
        for B in brs:
            eq2 += list(B * V - V * B)
        (vv,) = sp.linsolve(eq2, vs)
        # order 0 of L7: Y generic in p, Z = Y^2 - tr(Y^2)/4 commutes with Y and is independent of it
        rng = random.Random(7 + eta.count(1))
        Yg = sum((rng.randint(-3, 3) * pb[i] for i in range(9)), sp.zeros(4))
        Zg = Yg * Yg - (Yg * Yg).trace() / 4 * sp.eye(4)
        okY = Yg * Zg == Zg * Yg and sp.Matrix([list(Yg), list(Zg)]).rank() == 2 and (E * Zg).T == E * Zg
        check(f"base {sigl}: own bases of p (9, eta-self-adjoint traceless) and so(eta) (6); [p, p] has rank {rk} and "
              f"equals so(eta) (rank with so(eta) {rkk}); the commutant of so(eta) in gl(4) is R 1 (absolutely "
              f"irreducible, not only irreducible); [[X, Y], v] = 0 for all X, Y forces v = 0; a random Y and "
              f"Z = Y^2 - tr(Y^2)/4 are commuting, independent elements of p", okp and rk == 6 and rkk == 6 and okc and
              all(c == 0 for c in vv) and okY)


# =============================================================================================== X4 L5
def x4(tabs):
    section("X4  Proposition L5: the curvature-free channel E' = y E1 - x E13 (own assembly; residual computed twice)")
    gens = (x, y, z, w, t, u)
    for ep in (1, -1):
        t1 = time.time()
        T = tabs[ep]
        lab = LAB[ep]
        ell = x * z + ep * y * w
        N = x ** 2 + ep * y ** 2
        de = t * y - u * x
        S, M, D, NS = slaving(T, {kap: 1})
        lE1, lE13 = lrow(T, 'E1'), lrow(T, 'E13')
        lam = [y * lE1[k] - x * lE13[k] for k in range(4)]                      # coefficients of E'
        nh = [sp.expand(sum(lam[k] * NS[k, a] for k in range(4))) for a in range(2)]   # h = nh / D
        vt, rt = sp.div(gpoly(nh[1], gens), gpoly(576 * ell, gens))
        vu, ru = sp.div(gpoly(-nh[0], gens), gpoly(576 * ell, gens))
        okpoly = rt.is_zero and ru.is_zero
        v = (vt.as_expr(), vu.as_expr())
        # direct route: integral curve gamma' = v, a_k = Sigma_k . v = n_k / D, a_k' = v . grad a_k
        nk = [sp.expand(NS[k, 0] * v[0] + NS[k, 1] * v[1]) for k in range(4)]
        Dp = gpoly(D, gens)

        def vgrad(e):
            return sp.expand(v[0] * sp.diff(e, t) + v[1] * sp.diff(e, u))
        vD = vgrad(D)
        hv = sp.expand(sum(lam[k] * nk[k] for k in range(4)))
        num = 0
        for k in range(4):
            if lam[k] != 0:
                num += lam[k] * (vgrad(nk[k]) * D - nk[k] * vD)
        for (ik, k), (il, ll) in itertools.combinations_with_replacement(list(enumerate(KS)), 2):
            b1, b2, _ = T['B'][('E1', k, ll)]
            c1_, c2_, _ = T['B'][('E13', k, ll)]
            beta = y * (b1 + b2) - x * (c1_ + c2_)
            if beta != 0:
                wkl = sp.Rational(1, 2) if k == ll else 1
                num += wkl * beta * nk[ik] * nk[il]
        target = -1916928 * N * ell ** 2 * de
        okR = gpoly(num, gens) == gpoly(target, gens) * Dp ** 2
        check(f"{lab} kappa = 1, symbolic in x, y, z, w, t, u: v = det M (h_u, -h_t)/(576 l) is a polynomial ({okpoly}); "
              f"along an integral curve of v the (Hess rho)_0 coefficient h.v of E' vanishes, and its (d rho (x) d rho)_0 "
              f"coefficient (direct route, sum_k l'_k a_k' + sum w_kl (b1+b2)' a_k a_l) is "
              f"R' = -1916928 (x^2 {'+' if ep == 1 else '-'} y^2)(xz {'+' if ep == 1 else '-'} yw)^2 (ty - ux) "
              f"{tick(t1)}", okpoly and hv == 0 and okR)
        # the note's route: R' = h.(Dv)v + v^T Q''_sym v, with own H, Q bookkeeping on numerators:
        # Sigma = NS/D, so H D = assemble(NS) and Q D^2 = assemble(NS, d NS D - NS d D)
        t1 = time.time()
        Tk = {**T, 'L': {k: vv.subs(kap, 1) for k, vv in T['L'].items()},
              'B': {k: (tuple(a.subs(kap, 1) for a in vv) if isinstance(vv, tuple) else vv.subs(kap, 1))
                    for k, vv in T['B'].items()}}
        dNt = (NS.diff(t) * D - NS * sp.diff(D, t)).applyfunc(sp.expand)
        dNu = (NS.diff(u) * D - NS * sp.diff(D, u)).applyfunc(sp.expand)
        Hn, Qn, _, _ = assemble(Tk, NS, dNt, dNu)
        hrow = [sp.expand(y * Hn['E1'][a] - x * Hn['E13'][a]) for a in range(2)]          # = h D
        okh = all(zero(hrow[a] - nh[a]) for a in range(2))
        Dv = [[sp.diff(v[a], (t, u)[b]) for b in range(2)] for a in range(2)]
        hDvv = sum(hrow[a] * Dv[a][b] * v[b] for a in range(2) for b in range(2))         # (h.(Dv)v) D
        Qpp = {ab: y * Qn['E1'][ab] - x * Qn['E13'][ab] for ab in Qn['E1']}             # Q''_sym D^2
        vQv = sum(v[a] * Qpp[(a, b)] * v[b] for a in range(2) for b in range(2))
        okN = gpoly(hDvv * D + vQv, gens) == gpoly(target, gens) * Dp ** 2
        check(f"{lab} the note's route R' := h.(Dv)v + v^T Q''_sym v with own H, Q: h = (y, -x) H equals the direct row, "
              f"and R' is the same polynomial (L5(a), L5(d): q^ = mu^2 R') {tick(t1)}", okh and okN)
        if ep == 1:
            # L.0 / L8 step 3: in the u(128) real form ((x, y, z, w) = -i(x', ..), t, u, gamma', gamma'' in iR) the
            # reduced coefficients h^ = h.gamma' and q^ = h.gamma'' + gamma'^T Q'' gamma' are real, so F''/F' is real
            rng = random.Random(41)
            okre = True
            for _ in range(3):
                pt = {x: -I * sp.Rational(rng.randint(-5, 5), 3), y: -I * sp.Rational(rng.randint(-5, 5), 3),
                      z: -I * sp.Rational(rng.randint(1, 5), 3), w: -I * sp.Rational(rng.randint(1, 5), 3),
                      t: I * sp.Rational(rng.randint(-3, 3), 50), u: I * sp.Rational(rng.randint(-3, 3), 50)}
                g1 = [I * sp.Rational(rng.randint(-4, 4), 5) for _ in range(2)]
                g2 = [I * sp.Rational(rng.randint(-4, 4), 5) for _ in range(2)]
                Dn = D.subs(pt)
                hh = [nh[a].subs(pt) / Dn for a in range(2)]
                QQn = {ab: Qpp[ab].subs(pt) / Dn ** 2 for ab in Qpp}
                hhat = sp.expand_complex(sum(hh[a] * g1[a] for a in range(2)))
                qhat = sp.expand_complex(sum(hh[a] * g2[a] for a in range(2)) +
                                         sum(g1[a] * QQn[(a, b)] * g1[b] for a in range(2) for b in range(2)))
                okre &= sp.im(hhat) == 0 and sp.im(qhat) == 0 and hhat != 0 and Dn != 0
            check("(14,0) real form: at three random exact points with (x, y, z, w) = -i(x', y', z', w'), t, u in iR and "
                  "curve jets gamma', gamma'' in iR^2, the reduced coefficients h^ = h.gamma' and q^ = h.gamma'' + "
                  "gamma'^T Q'' gamma' of E' are real (so F''/F' = q^/h^ is real, L.0 and L8 step 3)", okre)
        # invariance of the line
        vl = sp.expand((y * v[0] - x * v[1]).subs({t: sig * x, u: sig * y}))
        check(f"{lab} on the line (t, u) = sigma (x, y): v.grad(ty - ux) = -(xw - yz)(160 N sigma - 1) (kappa = 1), so the "
              f"line is v-invariant iff xw = yz", zero(vl + (x * w - y * z) * (160 * N * sig - 1)))
        # general kappa at random exact members: R' = -1916928 kappa^3 N l^2 delta, v.grad delta = -kappa^2 (..)(160 N s - kappa)
        rng = random.Random(31 + ep)
        okk = True
        for _ in range(2):
            vals = {x: sp.Rational(rng.randint(-5, 5), rng.randint(1, 4)), y: sp.Rational(rng.randint(-5, 5), rng.randint(1, 4)),
                    z: sp.Rational(rng.randint(-5, 5), rng.randint(1, 4)), w: sp.Rational(rng.randint(-5, 5), rng.randint(1, 4)),
                    kap: sp.Rational(rng.randint(2, 9), rng.randint(1, 5))}
            Tn = {**T, 'S': T['S'].subs(vals), 'M': T['M'].subs(vals)}
            S_, M_, D_, NS_ = slaving(Tn, {})
            l1 = [e.subs(vals) for e in lE1]
            l13 = [e.subs(vals) for e in lE13]
            lm = [vals[y] * l1[k] - vals[x] * l13[k] for k in range(4)]
            ellv, Nv = ell.subs(vals), N.subs(vals)
            nh_ = [sp.expand(sum(lm[k] * NS_[k, a] for k in range(4))) for a in range(2)]
            v_ = (sp.cancel(nh_[1] / (576 * ellv)), sp.cancel(-nh_[0] / (576 * ellv)))
            nk_ = [sp.expand(NS_[k, 0] * v_[0] + NS_[k, 1] * v_[1]) for k in range(4)]

            def vg(e):
                return sp.expand(v_[0] * sp.diff(e, t) + v_[1] * sp.diff(e, u))
            vD_ = vg(D_)
            nm_ = sum(lm[k] * (vg(nk_[k]) * D_ - nk_[k] * vD_) for k in range(4))
            for (ik, k), (il, ll) in itertools.combinations_with_replacement(list(enumerate(KS)), 2):
                b1, b2, _ = T['B'][('E1', k, ll)]
                c1_, c2_, _ = T['B'][('E13', k, ll)]
                beta = (vals[y] * (b1 + b2) - vals[x] * (c1_ + c2_)).subs(vals)
                if beta != 0:
                    nm_ += (sp.Rational(1, 2) if k == ll else 1) * beta * nk_[ik] * nk_[il]
            tg = -1916928 * vals[kap] ** 3 * Nv * ellv ** 2 * (t * vals[y] - u * vals[x])
            okk &= sp.Poly(sp.expand(nm_), t, u, domain='QQ_I') == sp.Poly(sp.expand(tg * D_ ** 2), t, u, domain='QQ_I')
            vl_ = sp.expand((vals[y] * v_[0] - vals[x] * v_[1]).subs({t: sig * vals[x], u: sig * vals[y]}))
            okk &= zero(vl_ + vals[kap] ** 2 * (vals[x] * vals[w] - vals[y] * vals[z]) * (160 * Nv * sig - vals[kap]))
        check(f"{lab} general kappa, two random exact members (t, u symbolic): R' = -1916928 kappa^3 N l^2 delta and, on "
              f"the line, v.grad delta = -kappa^2 (xw - yz)(160 N sigma - kappa), zero only at sigma* = kappa/(160 N)", okk)
        # L5(c): the family x W_A + z W_C AdG (kappa symbolic); kernel field v_c of the E13 channel and its residual
        sub = {y: 0, w: 0}
        S_, M_, D_, NS_ = slaving(T, sub)
        vc = (kap ** 2 + 224 * kap * t * x - 3328 * ep * u ** 2 * x ** 2, 384 * kap * u * x)
        l13 = [e.subs(sub) for e in lE13]
        nkc = [sp.expand(NS_[k, 0] * vc[0] + NS_[k, 1] * vc[1]) for k in range(4)]

        def vgc(e):
            return sp.expand(vc[0] * sp.diff(e, t) + vc[1] * sp.diff(e, u))
        hvc = sp.expand(sum(l13[k] * nkc[k] for k in range(4)))
        numc = sum(l13[k] * (vgc(nkc[k]) * D_ - nkc[k] * vgc(D_)) for k in range(4))
        for (ik, k), (il, ll) in itertools.combinations_with_replacement(list(enumerate(KS)), 2):
            c1_, c2_, _ = T['B'][('E13', k, ll)]
            be = (c1_ + c2_).subs(sub)
            if be != 0:
                numc += (sp.Rational(1, 2) if k == ll else 1) * be * nkc[ik] * nkc[il]
        okc = hvc == 0 and gpoly(numc, (x, z, kap, t, u)) == gpoly(-1916928 * kap * u * x ** 4 * D_ ** 2, (x, z, kap, t, u))
        # L5(a)'s kernel field for this family (kappa symbolic): v_a = det M (h_u, -h_t)/(576 l), h = (0, -x) H
        nha = [sp.expand(sum(-x * l13[k] * NS_[k, a] for k in range(4))) for a in range(2)]
        va = (sp.cancel(nha[1] / (576 * x * z)), sp.cancel(-nha[0] / (576 * x * z)))
        okc &= all(sp.cancel(va[i] + kap * z * vc[i]) == 0 for i in range(2))
        check(f"{lab} L5(c), family x W_A + z W_C AdG, kappa symbolic: v_c = (kappa^2 + 224 kappa t x "
              f"{'-' if ep == 1 else '+'} 3328 u^2 x^2, 384 kappa u x) is a kernel field of the E13 row (H13.v_c = 0) and "
              f"the E13 residual along it is -1916928 kappa u x^4 (independent of z); L5(a)'s kernel field is v_a = -kappa z v_c, "
              f"so L5(a)'s R' = -x (kappa z)^2 R13 (consistent with -1916928 kappa^3 N l^2 delta at y = w = 0)", okc)
        # x W_A: Sigma, H13 = 0, E13 = c13 g - 16 x Sigma_C Sigma_B (dt (x) dt - |dt|^2 g)
        sub = {y: 0, z: 0, w: 0}
        S_, M_, D_, NS_ = slaving(T, sub)
        Sg = (NS_ / D_).applyfunc(sp.cancel)
        Ts = {**T, 'L': {k: vv.subs(sub) for k, vv in T['L'].items()},
              'B': {k: (tuple(a.subs(sub) for a in vv) if isinstance(vv, tuple) else vv.subs(sub)) for k, vv in T['B'].items()}}
        H, Qs, Qp, W = assemble(Ts, Sg, Sg.diff(t), Sg.diff(u))
        dA = kap ** 2 + 224 * kap * t * x - 3328 * ep * u ** 2 * x ** 2
        SC, SB = 24 * kap * x / dA, 4992 * u * x ** 2 / dA
        okWA = sp.cancel(Sg[1, 0] - SC) == 0 and sp.cancel(Sg[3, 0] - SB) == 0 and \
            all(sp.cancel(Sg[k, a]) == 0 for k in range(4) for a in range(2) if (k, a) not in ((1, 0), (3, 0)))
        okWA &= all(sp.cancel(H['E13'][a]) == 0 for a in range(2))
        okWA &= sp.cancel(Qs['E13'][(0, 0)] + 16 * x * SC * SB) == 0 and sp.cancel(Qp['E13'][(0, 0)] - 16 * x * SC * SB) == 0
        okWA &= all(sp.cancel(Qs['E13'][ab]) == 0 and sp.cancel(Qp['E13'][ab]) == 0 for ab in ((0, 1), (1, 0), (1, 1)))
        check(f"{lab} members x W_A: Sigma = (0, Sigma_C, 0, Sigma_B) dt with Sigma_C = 24 kappa x/d, Sigma_B = 4992 u x^2/d; "
              f"H13 = 0; E13 = c13 g - 16 x Sigma_C Sigma_B (dt (x) dt - |dt|^2 g) (L5(c))", okWA)


# =============================================================================================== X5 L8 step 5
def x5(tabs):
    section("X5  Theorem L8 step 5: the line (t, u) = sigma (x, y), (z, w) = c (x, y) (own direct assembly, symbolic in "
            "x, y, c, kappa, sigma)")
    for ep in (1, -1):
        t1 = time.time()
        T = tabs[ep]
        lab = LAB[ep]
        N = x ** 2 + ep * y ** 2
        sub = {z: cc * x, w: cc * y}
        Tm = {**T, 'S': T['S'].subs(sub), 'M': T['M'].subs(sub),
              'L': {k: vv.subs(sub) for k, vv in T['L'].items()},
              'B': {k: (tuple(a.subs(sub) for a in vv) if isinstance(vv, tuple) else vv.subs(sub)) for k, vv in T['B'].items()}}
        line = {t: sig * x, u: sig * y}
        M_ = Tm['M'].subs(line)
        D_ = sp.factor(M_.det())
        Sg = (-M_.adjugate() * Tm['S'] / D_).applyfunc(sp.cancel)
        a = [sp.cancel(Sg[k, 0] * x + Sg[k, 1] * y) for k in range(4)]          # A_k = a_k(sigma) d sigma
        ap = [sp.cancel(sp.diff(a[k], sig)) for k in range(4)]
        d1 = kap + 224 * N * sig
        ok = zero(sp.expand(D_) - sp.expand(kap ** 2 * d1 * (kap + 544 * N * sig)))
        okc = True
        for nm, xi in (('E1', x), ('E13', y)):
            l = lrow(Tm, nm)
            Hc = sp.cancel(sum(l[k] * a[k] for k in range(4)))
            Qc = sum(l[k] * ap[k] for k in range(4))
            Qpc = -sum(l[k] * ap[k] for k in range(4))
            for (ik, k), (il, ll) in itertools.combinations_with_replacement(list(enumerate(KS)), 2):
                b1, b2, b3 = Tm['B'][(nm, k, ll)]
                wkl = sp.Rational(1, 2) if k == ll else 1
                Qc += wkl * (b1 + b2) * a[ik] * a[il]
                Qpc += wkl * b3 * a[ik] * a[il]
            Qc, Qpc = sp.cancel(Qc), sp.cancel(Qpc)
            ok &= sp.cancel(Hc - xi * 576 * N / d1) == 0
            ok &= sp.cancel(Qc + xi * 258048 * N ** 2 / d1 ** 2) == 0
            ok &= sp.cancel(Qpc - xi * 78336 * N ** 2 / d1 ** 2) == 0
            okc &= cc not in (Hc.free_symbols | Qc.free_symbols | Qpc.free_symbols)
            ci = T['c1'] if nm == 'E1' else T['c13']
            ok &= zero(ci.subs(line) - xi * (kap * sig + 312 * N * sig ** 2))
        # antisymmetric and Lambda^2 parts vanish identically on the line (A_k proportional to d sigma)
        F = -1 / (224 * N * d1)
        Fp, Fpp = sp.diff(F, sig), sp.diff(F, sig, 2)
        Hn, Qn = 576 * N / d1, -258048 * N ** 2 / d1 ** 2
        okF = sp.cancel(Fpp / Fp - Qn / Hn) == 0
        mu = sp.cancel(Fp / Hn)
        okF &= sp.cancel(mu + sp.Rational(7, 18) * F) == 0
        check(f"{lab} on the line: det M = kappa^2 (kappa + 224 N sigma)(kappa + 544 N sigma); E_i = x_i E~ with Hessian "
              f"coefficient 576 N/(kappa + 224 N sigma), gradient coefficients -258048 N^2/(..)^2, 78336 N^2/(..)^2 and "
              f"zeroth order (kappa sigma + 312 N sigma^2) g, N = x^2 {'+' if ep == 1 else '-'} y^2, all independent of c "
              f"{tick(t1)}", ok and okc)
        check(f"{lab} r = F(sigma) = -1/(224 N (kappa + 224 N sigma)) solves F''/F' = Q~/H~, and (Hess r)_0 = mu (Ric)_0 with "
              f"mu = F'/H~ = -(7/18) r exactly (sign: E~ = -G + H~ (Hess sigma - Box sigma g) + ...)", okF)


# =============================================================================================== X6 L7
def pq_bases(eta):
    n = 4
    pb = []
    for i in range(3):
        Mx = np.zeros((n, n), dtype=np.int64)
        Mx[i, i], Mx[i + 1, i + 1] = 1, -1
        pb.append(Mx)
    for i in range(n):
        for j in range(i + 1, n):
            Mx = np.zeros((n, n), dtype=np.int64)
            Mx[i, j], Mx[j, i] = eta[i], eta[j]
            pb.append(Mx)
    kb = []
    for i in range(n):
        for j in range(i + 1, n):
            Mx = np.zeros((n, n), dtype=np.int64)
            Mx[i, j], Mx[j, i] = eta[i], -eta[j]
            kb.append(Mx)
    return np.array(pb), np.array(kb)


def x6a():
    section("X6a Lemma L7: the referee's (C1) and its pairings (own dual bases of p; exact integers)")
    for eta in ((1, 1, 1, 1), (1, -1, -1, -1), (-1, 1, 1, 1), (1, 1, -1, -1)):
        sigl = f"({eta.count(1)},{eta.count(-1)})"
        Xb, Kb = pq_bases(eta)
        G = [[Fr(int(np.trace(Xb[a] @ Xb[b]))) for b in range(9)] for a in range(9)]
        Gi = finv(G)
        den = 1
        for r in Gi:
            for v in r:
                den = den * v.denominator // np.gcd(den, v.denominator)
        Xd = np.array([sum(int(Gi[a][b] * den) * Xb[b] for b in range(9)) for a in range(9)])   # den * dual basis
        okd = all(int(np.trace(Xb[a] @ Xd[b])) == den * (a == b) for a in range(9) for b in range(9))

        def comm(A, B):
            return np.einsum('...ij,...jk->...ik', A, B) - np.einsum('...ij,...jk->...ik', B, A)

        def c0t(Xs):            # 4 c0(X_a, X_b, X_c) = -[[X_b, X_c], X_a]  -> (9, 9, 9, 4, 4)
            BC = comm(Xs[None, :, None], Xs[None, None, :])
            return -comm(BC, Xs[:, None, None])

        def c1t(Xs):            # c1 = tr(X_a X_b) X_c - tr(X_a X_c) X_b
            g = np.einsum('aij,bji->ab', Xs, Xs)
            return np.einsum('ab,cij->abcij', g, Xs) - np.einsum('ac,bij->abcij', g, Xs)

        def qt(Xs, K, Kp):      # 4 q(K (x) K') with 4 d(X, Y) Z = tr(Z [K, X]) [K', Y] + tr(Z [K, Y]) [K', X]
            P = comm(K[None], Xs)            # [K, X_a]
            Q = comm(Kp[None], Xs)           # [K', X_b]
            tau = np.einsum('cij,aji->ca', Xs, P)          # tr(X_c [K, X_a])
            dd = np.einsum('ca,bij->abcij', tau, Q) + np.einsum('cb,aij->abcij', tau, Q)   # 4 d(X_a, X_b) X_c
            return dd - np.transpose(dd, (0, 2, 1, 3, 4))  # 4 q = 4 d(a,b)c - 4 d(a,c)b
        C0, C1 = c0t(Xb), c1t(Xb)
        C0d, C1d = c0t(Xd), c1t(Xd)

        def pair(E, Fd, scale):
            v = int(np.einsum('abcij,abcji->', E, Fd))
            return Fr(v, scale)
        D3 = den ** 3
        p00 = pair(C0, C0d, 16 * D3)
        p01 = pair(C0, C1d, 4 * D3)
        p11 = pair(C1, C1d, D3)
        okq = True
        for K in Kb:
            for Kp in Kb:
                Qt = qt(Xb, K, Kp)
                assert np.abs(Qt).max() < 2 ** 40
                tkk = int(np.trace(K @ Kp))
                okq &= pair(Qt, C1d, 4 * D3) == -9 * tkk and pair(Qt, C0d, 16 * D3) == Fr(-27, 4) * tkk
                Tq = 4 * pair(Qt, C0d, 16 * D3) - 3 * pair(Qt, C1d, 4 * D3)
                okq &= Tq == 0
        check(f"base {sigl}: dual basis of p (x {den}); <c0, c0> = {p00}, <c0, c1> = {p01}, <c1, c1> = {p11}; on all 36 "
              f"basis elements K (x) K' of so (x) so: <q, c1> = -9 tr(KK'), <q, c0> = -(27/4) tr(KK'); so with "
              f"T = 4 c0 - 3 c1: <q(B), T> = 0 for every B, <c0, T> = {4 * p00 - 3 * p01}, <c1, T> = {4 * p01 - 3 * p11}",
              okd and p00 == Fr(27, 2) and p01 == 18 and p11 == 144 and okq and 4 * p00 - 3 * p01 == 0 and
              4 * p01 - 3 * p11 == -360)
        # Casimir structure (justifies the control of X6c): q(Cas) = -3 c0, sum eps eps R1_ij (x) R1_ij = -4 Cas
        kg = [[Fr(int(np.trace(Kb[a] @ Kb[b]))) for b in range(6)] for a in range(6)]
        kgi = finv(kg)
        q_cas = sum((qt(Xb, Kb[a], Kb[b]).astype(object) * kgi[a][b] for a in range(6) for b in range(6) if kgi[a][b]),
                    np.zeros((9, 9, 9, 4, 4), dtype=object))
        okcas = all(v == 0 for v in (q_cas * Fr(1, 4) + C0.astype(object) * Fr(3, 4)).flatten())
        R1 = {(i, j): np.array([[(int(k == i) * eta[j] * int(l == j) - int(k == j) * eta[i] * int(l == i))
                                 for l in range(4)] for k in range(4)], dtype=np.int64) for i in range(4) for j in range(4)}
        lhs = sum(eta[i] * eta[j] * np.einsum('ij,kl->ijkl', R1[(i, j)], R1[(i, j)]) for i in range(4) for j in range(4))
        cas = sum(np.einsum('ij,kl->ijkl', Kb[a], Kb[b]).astype(object) * kgi[a][b] for a in range(6) for b in range(6))
        okR1 = all(v == 0 for v in (lhs.astype(object) + 4 * cas).flatten())
        check(f"base {sigl}: for B = Cas (the so(eta) Casimir w.r.t. tr), q(Cas) = -3 c0; and for the unit constant-curvature "
              f"tensor R1 (R(X,Y)Z = g(Y,Z)X - g(X,Z)Y), sum eps_i eps_j R1_ij (x) R1_ij = -4 Cas, so q(B_1) = 12 (c_V/c_H)^2 "
              f"c0 and (C1) at mu' = 0 is solvable iff the A-term is scaled by s* = -(c_H/c_V)^2/12", okcas and okR1)


def rand_curv(rng, eta):
    """random algebraic curvature tensor R_abcd (sum of Kulkarni-Nomizu products of random symmetric forms)"""
    n = 4
    R = [[[[Fr(0)] * n for _ in range(n)] for _ in range(n)] for _ in range(n)]
    for _ in range(3):
        h = [[0] * n for _ in range(n)]
        k = [[0] * n for _ in range(n)]
        for i in range(n):
            for j in range(i, n):
                h[i][j] = h[j][i] = rng.randint(-2, 2)
                k[i][j] = k[j][i] = rng.randint(-2, 2)
        for a, b, c, d in itertools.product(range(n), repeat=4):
            R[a][b][c][d] += Fr(h[a][c] * k[b][d] + h[b][d] * k[a][c] - h[a][d] * k[b][c] - h[b][c] * k[a][d])
    return R


def endos(R, eta):
    """B_cd = R(e_c, e_d) as matrices: (B_cd)^a_b = eta_a R_abcd"""
    return {(c, d): [[Fr(eta[a]) * R[a][b][c][d] for b in range(4)] for a in range(4)] for c in range(4) for d in range(4)}


def x6b():
    section("X6b Lemma L7 inputs from definitions: the fibre connection, its curvature (c0), Ric_F, and the A-term of (C1)")
    rng = random.Random(20261109)
    for eta in ((1, -1, -1, -1), (1, 1, -1, -1)):
        sigl = f"({eta.count(1)},{eta.count(-1)})"
        etaM = [[Fr(eta[i]) if i == j else Fr(0) for j in range(4)] for i in range(4)]
        lam, cH, cV = Fr(1, 2), Fr(3), Fr(2)
        # random point of the fibre: y = P^T eta P
        while True:
            Pm = [[Fr(rng.randint(-2, 2) + (3 if i == j else 0)) for j in range(4)] for i in range(4)]
            if sp.Matrix(Pm).det() != 0:
                break
        yy = fmul(fmul(fT(Pm), etaM), Pm)
        yi = finv(yy)

        def Gam(h, k, yinv=yi):
            return fscl(fadd(fmul(fmul(h, yinv), k), fmul(fmul(k, yinv), h)), Fr(-1, 2))

        def gF(h, k, yinv=yi, lmb=lam):
            return cV * (ftr(fmul(fmul(yinv, h), fmul(yinv, k))) - lmb * ftr(fmul(yinv, h)) * ftr(fmul(yinv, k)))

        def rsym():
            Mx = [[Fr(0)] * 4 for _ in range(4)]
            for i in range(4):
                for j in range(i, 4):
                    Mx[i][j] = Mx[j][i] = Fr(rng.randint(-3, 3))
            return Mx
        U, h, k = rsym(), rsym(), rsym()
        # D_U g(h, k) at y
        yUy = fmul(fmul(yi, U), yi)
        okLC = Gam(U, h) == Gam(h, U)
        for lmb in (Fr(0), lam, Fr(1)):
            dg = cV * (-ftr(fmul(fmul(yUy, h), fmul(yi, k))) - ftr(fmul(fmul(yi, h), fmul(yUy, k))) +
                       lmb * (ftr(fmul(yUy, h)) * ftr(fmul(yi, k)) + ftr(fmul(yi, h)) * ftr(fmul(yUy, k))))
            okLC &= dg == gF(Gam(U, h), k, lmb=lmb) + gF(h, Gam(U, k), lmb=lmb)

        def DGam(Uv, h_, k_):
            return fscl(fadd(fmul(fmul(fmul(h_, yi), fmul(Uv, yi)), k_), fmul(fmul(fmul(k_, yi), fmul(Uv, yi)), h_)), Fr(1, 2))

        def Rc(Uv, Vv, Wv):
            return fadd(fadd(DGam(Uv, Vv, Wv), DGam(Vv, Uv, Wv), -1), fadd(Gam(Uv, Gam(Vv, Wv)), Gam(Vv, Gam(Uv, Wv)), -1))
        # Ric_F(V, W) = tr(U -> R(U, V) W) in the coordinate basis E_A (coordinates = upper triangle)
        V, Wm = rsym(), rsym()
        ricF = sum(coords(Rc(Emat(A), V, Wm))[A] for A in range(10))
        G14 = ftr(fmul(fmul(yi, V), fmul(yi, Wm))) - Fr(1, 4) * ftr(fmul(yi, V)) * ftr(fmul(yi, Wm))
        okRic = ricF == -G14
        # curvature at eta in the endomorphism picture: R(eta Y, eta Z)(eta X) = eta (-(1/4)[[Y, Z], X])
        yi_e = etaM

        def DGe(Uv, h_, k_):
            return fscl(fadd(fmul(fmul(fmul(h_, yi_e), fmul(Uv, yi_e)), k_), fmul(fmul(fmul(k_, yi_e), fmul(Uv, yi_e)), h_)),
                        Fr(1, 2))

        def Ge(h_, k_):
            return Gam(h_, k_, yi_e)

        def Re(Uv, Vv, Wv):
            return fadd(fadd(DGe(Uv, Vv, Wv), DGe(Vv, Uv, Wv), -1), fadd(Ge(Uv, Ge(Vv, Wv)), Ge(Vv, Ge(Uv, Wv)), -1))
        Xb, Kb = pq_bases(eta)
        Xf = [[[Fr(int(v)) for v in r] for r in Mx] for Mx in Xb]
        Kf = [[[Fr(int(v)) for v in r] for r in Mx] for Mx in Kb]

        def rp():
            cf_ = [Fr(rng.randint(-2, 2)) for _ in range(9)]
            return [[sum(cf_[a] * Xf[a][i][j] for a in range(9)) for j in range(4)] for i in range(4)]
        Xr, Yr, Zr = rp(), rp(), rp()
        lhs = Re(fmul(etaM, Yr), fmul(etaM, Zr), fmul(etaM, Xr))
        okc0 = fmul(etaM, lhs) == fscl(fbr(fbr(Yr, Zr), Xr), Fr(-1, 4))
        check(f"base {sigl}: Gamma(h, k) = -(1/2)(h y^-1 k + k y^-1 h) is metric for c_V G_lambda (lambda = 0, 1/2, 1) "
              f"and symmetric at a random y = P^T eta P, hence the Levi-Civita connection for every lambda tested; its curvature gives Ric_F = -G_1/4 "
              f"there, and at eta, in H = y^-1 h, R(Y, Z) X = -(1/4)[[Y, Z], X] (the c0 term of (C1))", okLC and okRic and okc0)
        # A-term from the definition: f_ij(U)(y) = G_y(U, R_ij.y); d/ds at eta along y = eta + s eta X, U = eta Z
        R = rand_curv(rng, eta)
        Bc = endos(R, eta)

        def Bdot(Bm, m_):
            return fscl(fadd(fmul(fT(Bm), m_), fmul(m_, Bm)), -1)

        def fder(U_, Xd_, Bm):
            # exact derivative at s = 0 (product rule) of f(s) = tr(y^-1 U y^-1 (B.y)) - lam tr(y^-1 U) tr(y^-1 B.y),
            # y = eta + s H, H = eta X:  d(y^-1) = -eta H eta
            Hm = fmul(etaM, Xd_)
            Y1 = fscl(fmul(fmul(etaM, Hm), etaM), -1)
            B0, B1 = Bdot(Bm, etaM), Bdot(Bm, Hm)
            return (ftr(fmul(fmul(Y1, U_), fmul(etaM, B0))) + ftr(fmul(fmul(etaM, U_), fmul(Y1, B0))) +
                    ftr(fmul(fmul(etaM, U_), fmul(etaM, B1))) -
                    lam * (ftr(fmul(Y1, U_)) * ftr(fmul(etaM, B0)) +
                           ftr(fmul(etaM, U_)) * (ftr(fmul(Y1, B0)) + ftr(fmul(etaM, B1)))))

        def fder_sym(U_, Xd_, Bm):
            # the same derivative by sympy differentiation of the exact rational function (spot check of fder)
            s_ = sp.Symbol('s')
            ys = sp.Matrix(4, 4, lambda i, j: sp.Rational(etaM[i][j]) + s_ * sp.Rational(fmul(etaM, Xd_)[i][j]))
            yinv = ys.inv()
            Bs = sp.Matrix(4, 4, lambda i, j: sp.Rational(Bm[i][j]))
            By = -(Bs.T * ys + ys * Bs)
            Us = sp.Matrix(4, 4, lambda i, j: sp.Rational(U_[i][j]))
            f = (yinv * Us * yinv * By).trace() - sp.Rational(lam) * (yinv * Us).trace() * (yinv * By).trace()
            return sp.Rational(sp.diff(f, s_).subs(s_, 0))
        # Hess (AU, AW) at eta in directions X, Y: c sum eps eps [df(U)(X) df(W)(Y) + df(U)(Y) df(W)(X)]
        cA = cV ** 2 / (4 * cH ** 2)
        Wb = []                       # basis of eta-self-adjoint matrices (tangent vectors U = eta W~)
        for i in range(4):
            for j in range(i, 4):
                Mx = [[Fr(0)] * 4 for _ in range(4)]
                Mx[i][j] = Fr(eta[i])
                Mx[j][i] = Fr(eta[j])
                Wb.append(Mx)
        dfX = {}
        for ij, Bm in Bc.items():
            for key, D_ in (('X', Xr), ('Y', Yr)):
                for wi, Wt in enumerate(Wb + [Zr]):
                    dfX[(ij, key, wi)] = fder(fmul(etaM, Wt), D_, Bm)
        rhs = []
        for wi in range(10):
            rhs.append(cA * sum(eta[i] * eta[j] * (dfX[((i, j), 'X', 10)] * dfX[((i, j), 'Y', wi)] +
                                                   dfX[((i, j), 'Y', 10)] * dfX[((i, j), 'X', wi)])
                                for i in range(4) for j in range(4)))
        # solve g_F(D, W~) = c_V [tr(D W~) - lam tr D tr W~] = rhs for D = sum d_m Wb_m
        Gm = [[cV * (ftr(fmul(Wb[m], Wb[wi])) - lam * ftr(Wb[m]) * ftr(Wb[wi])) for m in range(10)] for wi in range(10)]
        Gmi = finv(Gm)
        dco = [sum(Gmi[m][wi] * rhs[wi] for wi in range(10)) for m in range(10)]
        Dend = [[sum(dco[m] * Wb[m][i][j] for m in range(10)) for j in range(4)] for i in range(4)]
        lhsq = fscl(Dend, cV)                      # c_V A''(X, Y) Z
        # the note's closed form with B = (c_V/c_H)^2 sum eps eps R_ij (x) R_ij
        f2 = (cV / cH) ** 2
        cf = [[Fr(0)] * 4 for _ in range(4)]
        for i in range(4):
            for j in range(4):
                Bm = Bc[(i, j)]
                tX = ftr(fmul(Zr, fbr(Bm, Xr)))
                tY = ftr(fmul(Zr, fbr(Bm, Yr)))
                cf = fadd(cf, fscl(fadd(fscl(fbr(Bm, Yr), tX), fscl(fbr(Bm, Xr), tY)), f2 * eta[i] * eta[j] / 4))
        check(f"base {sigl}: the A-term of (C1) from the definition of (AU, AW) (exact derivative of G_y(U, R_ij.y) at "
              f"eta, random curvature tensor, random X, Y, Z in p): c_V (g_F^-1 Hess (A., A.))(X, Y) Z = (1/4) sum_s "
              f"[tr(Z[K_s, X])[K'_s, Y] + tr(Z[K_s, Y])[K'_s, X]] with B = (c_V/c_H)^2 sum eps_i eps_j R_ij (x) R_ij",
              lhsq == cf and any(v != 0 for r in cf for v in r) and
              fder(fmul(etaM, Zr), Xr, Bc[(0, 1)]) == fder_sym(fmul(etaM, Zr), Xr, Bc[(0, 1)]))


def fibre_system(eta, R, cH, cV, s_fac, mu1, lam=Fr(1, 2)):
    """Taylor system at y = eta (coordinates Y_A, y = eta + sum Y_A E_A) for
       Hess_F w = (1 + mu1 (w - 1)) (Ric_F + s_fac (A., A.)),  w(eta) = 1,  to order eps^2 (w of degree <= 4).
       Returns (rank of the coefficient matrix, rank of the augmented matrix, number of unknowns)."""
    T = TP(10, 2)
    etaM = [[Fr(eta[i]) if i == j else Fr(0) for j in range(4)] for i in range(4)]
    Ec = pm_const(T, etaM)
    Dl = [[T.var(SYMIDX.index((min(i, j), max(i, j)))) for j in range(4)] for i in range(4)]
    A1 = pm_mul(T, Ec, pm_mul(T, Dl, Ec))
    A2 = pm_mul(T, A1, pm_mul(T, Dl, Ec))
    yi = pm_add(T, pm_add(T, Ec, A1, -1), A2)
    Yv = pm_add(T, Ec, Dl)
    EA = [pm_const(T, Emat(A)) for A in range(10)]
    QA = [pm_mul(T, yi, EA[A]) for A in range(10)]
    Gam = {}
    for A in range(10):
        for B in range(A, 10):
            Mx = pm_scl(T, pm_add(T, pm_mul(T, EA[A], QA[B]), pm_mul(T, EA[B], QA[A])), Fr(-1, 2))
            Gam[(A, B)] = [Mx[i][j] for (i, j) in SYMIDX]          # coordinates Gamma^C_AB
    trQ = [pm_tr(T, QA[A]) for A in range(10)]
    Bc = endos(R, eta)
    fA = {}
    for ij, Bm in Bc.items():
        BmP = pm_const(T, Bm)
        BmT = pm_const(T, fT(Bm))
        By = pm_scl(T, pm_add(T, pm_mul(T, BmT, Yv), pm_mul(T, Yv, BmP)), -1)
        yBy = pm_mul(T, yi, By)
        for A in range(10):
            fA[(ij, A)] = T.add(pm_trprod(T, QA[A], yBy), T.mul(trQ[A], pm_tr(T, yBy)), -lam)
    cA = cV ** 2 / (4 * cH ** 2) * s_fac
    Bt = {}
    for A in range(10):
        for B in range(A, 10):
            ric = T.scl(T.add(pm_trprod(T, QA[A], QA[B]), T.mul(trQ[A], trQ[B]), Fr(-1, 4)), -1)
            at = T.ssum([T.scl(T.mul(fA[((i, j), A)], fA[((i, j), B)]), cA * eta[i] * eta[j])
                         for i in range(4) for j in range(4) if fA[((i, j), A)] and fA[((i, j), B)]])
            Bt[(A, B)] = T.add(ric, at)
    # unknowns: coefficients of w of degree 1..4
    W4 = TP(10, 4)
    unk = [k for k in range(len(W4.mons)) if 1 <= W4.dg[k] <= 4]
    uix = {k: i for i, k in enumerate(unk)}
    rows, rhs = [], []
    low = {T.mons[k]: k for k in range(len(T.mons))}
    for A in range(10):
        for B in range(A, 10):
            for n_ in T.mons:
                row = {}
                # d_A d_B w at monomial n_: coefficient of Y^(n_ + e_A + e_B)
                m = list(n_)
                m[A] += 1
                m[B] += 1
                fac = (m[A] * (m[A] - 1) if A == B else m[A] * m[B])
                row[uix[W4.ix[tuple(m)]]] = Fr(fac)
                # - Gamma^C_AB d_C w
                for C in range(10):
                    gc = Gam[(A, B)][C]
                    for kk, gv in gc.items():
                        n1 = T.mons[kk]
                        n2 = [p - q for p, q in zip(n_, n1)]
                        if min(n2) < 0:
                            continue
                        m2 = list(n2)
                        m2[C] += 1
                        j = uix[W4.ix[tuple(m2)]]
                        row[j] = row.get(j, 0) - gv * m2[C]
                # - mu1 (w - 1) Bt
                if mu1:
                    for kk, bv in Bt[(A, B)].items():
                        n1 = T.mons[kk]
                        n2 = tuple(p - q for p, q in zip(n_, n1))
                        if min(n2) < 0 or sum(n2) == 0:
                            continue
                        j = uix[W4.ix[n2]]
                        row[j] = row.get(j, 0) - mu1 * bv
                rows.append(row)
                rhs.append(Bt[(A, B)].get(low[n_], Fr(0)))
    nunk = len(unk)
    dA_, dAb = {}, {}
    for r, row in enumerate(rows):
        ent = {j: sp.QQ(v.numerator, v.denominator) for j, v in row.items() if v}
        if ent:
            dA_[r] = dict(ent)
        if rhs[r]:
            ent[nunk] = sp.QQ(rhs[r].numerator, rhs[r].denominator)
        if ent:
            dAb[r] = ent
    rA = DomainMatrix(dA_, (len(rows), nunk), sp.QQ).rank()
    rAb = DomainMatrix(dAb, (len(rows), nunk + 1), sp.QQ).rank()
    return rA, rAb, nunk, len(rows)


def x6c():
    section("X6c Lemma L7 by brute force: own Taylor solution of Hess_F w = mu(w)(Ric_F + s (A., A.)) at y = eta to order "
            "eps^2, mu(w) = 1 + mu' (w - 1), w(eta) = 1 (own Gamma, Ric_F = -G_1/4, (A., A.) from its definition)")
    rng = random.Random(20261110)
    cH, cV = Fr(3), Fr(2)
    sstar = -(cH / cV) ** 2 / 12
    R1 = {e: [[[[Fr(int(a == c) * e[a] * int(b == d) * e[b] - int(a == d) * e[a] * int(b == c) * e[b])
                 for d in range(4)] for c in range(4)] for b in range(4)] for a in range(4)]
          for e in ((1, -1, -1, -1), (1, 1, -1, -1))}
    zeroR = [[[[Fr(0)] * 4 for _ in range(4)] for _ in range(4)] for _ in range(4)]
    runs = [((1, -1, -1, -1), 'random R', None, Fr(1), Fr(-7, 18), False),
            ((1, -1, -1, -1), 'R = 0 (flat base)', zeroR, Fr(1), Fr(-7, 18), False),
            ((1, -1, -1, -1), 'R = R1, s = s*', 'R1', sstar, Fr(0), True),
            ((1, -1, -1, -1), 'R = R1, s = s*', 'R1', sstar, Fr(-7, 18), False),
            ((1, -1, -1, -1), 'R = R1, s = 1', 'R1', Fr(1), Fr(0), False),
            ((1, 1, -1, -1), 'random R', None, Fr(1), Fr(-7, 18), False),
            ((1, 1, -1, -1), 'R = R1, s = s*', 'R1', sstar, Fr(0), True),
            ((1, 1, 1, 1), 'random R', None, Fr(1), Fr(-7, 18), False)]
    for eta, nmR, Rsel, sfac, mu1, expect in runs:
        t1 = time.time()
        sigl = f"({eta.count(1)},{eta.count(-1)})"
        if Rsel is None:
            R = rand_curv(rng, eta)
        elif Rsel == 'R1':
            R = R1[eta]
        else:
            R = Rsel
        rA, rAb, nu, nr = fibre_system(eta, R, cH, cV, sfac, mu1)
        cons = rA == rAb
        check(f"base {sigl}, {nmR}, c_V/c_H = 2/3, mu' = {mu1}: Taylor system ({nr} equations, {nu} unknowns) is "
              f"{'consistent' if cons else 'inconsistent'} (ranks {rA}, {rAb}), as (C1) predicts "
              f"({'solvable' if expect else 'not solvable'}; s* = {sstar}) {tick(t1)}", cons == expect)


# =============================================================================================== main
def main():
    say("Referee cross-check of GU round 8 (L-nonradial-z2.md), REVIEW6-GU.  Exact arithmetic; fixed seeds; one core.")
    tabs = x0()
    x1(tabs)
    x2()
    x3()
    x4(tabs)
    x5(tabs)
    x6a()
    x6b()
    x6c()
    nf = [lab for lab, ok in RES if not ok]
    say(f"\nSUMMARY: {len(RES)} checks, {len(RES) - len(nf)} True, {len(nf)} False")
    say(f"runtime {time.time() - T0:.0f} s")


if __name__ == '__main__':
    if len(sys.argv) > 1:
        sel = sys.argv[1:]
        tabs = x0() if any(s in sel for s in ('x1', 'x4', 'x5')) else None
        for s in sel:
            f = globals()[s]
            f(tabs) if s in ('x1', 'x4', 'x5') else f()
        say(f"\n{len(RES)} checks; runtime {time.time() - T0:.0f} s")
    else:
        main()
