"""GU round 3 (G), Task 1: the complete Euler-Lagrange system of the bosonic action (9.4) and the Noether
identity of the tau'(H) gauge symmetry (draft (9.5)-(9.6): Xi = D Upsilon).

Formalization: formalization/geometric-unity/G-full-el-fermions.md, Sec. G.2 (Theorems G1-G3).
Full PDE treatment in d = 14 is infeasible; the statements G1-G3 are proved for every dimension, Shiab operator
and Ad-invariant pairing.  This script checks them EXACTLY in a non-abelian model where the Euler operator can be
evaluated without approximation:

  model M3: d = 3, eta = diag(1, 1, -1) (Hodge star), structure group H = SL(2,R), ad = sl(2,R), pairing Tr;
     background A0 = constant sl(2)-valued 1-form with F0 = A0 ^ A0 != 0 (plays nabla^0);
     eps(x) = U(f) L(g) U(k)  (U, L unipotent; eps^-1 = U(-k) L(-g) U(-f) is polynomial);
     Shiab (d-1 = 2, so Omega^2 -> Omega^2):  Sh_eps(w (x) X) = w (x) Ad_eps^-1 M Ad_eps X  with
        M_sa X = [Q,[Q,X]]  (Tr-symmetric: the self-adjoint case of Theorem D2; Psi in Lambda^0 (x) S^2 g*)
        M_ns X = [Q,[R,X]]  (not Tr-symmetric).
  Action density, draft (9.4) under H5 of D.1:  L = Tr(T ^ Sh_eps(F_B + 1/2 d_B T + 1/3 T^T)) + kappa/2 Tr(T ^ *T),
     B = A0 + eps^-1 d_A0 eps,  T = varpi - eps^-1 d_A0 eps,  (T^T)_ab = [T_a, T_b].
  Fields are random integer polynomials (degree 5 for f, g, k; 4 for varpi); all quantities are carried as exact
  truncated Taylor polynomials (degree <= 3 at x = 0) with 66 first-order dual directions = the jet coordinates
  of the variations gamma (delta eps = eps gamma, |N| <= 2) and alpha (delta varpi, |N| <= 1).  The Euler operator
  E_u = sum_N (-D)^N dL/du_N is then evaluated exactly (Python integers; the density is scaled by 6).
Tests
  N1  pointwise factorisation (Theorem G1): L(eps, varpi) = L_J(S) as Taylor polynomials, S = pi'(eps, varpi)
      = eps varpi eps^-1 - (d_A0 eps) eps^-1, L_J(S) = Tr(S ^ M(F0 + 1/2 d_0 S + 1/3 S^S)) + kappa/2 Tr(S ^ *S).
  N2  Noether identity (Theorem G2): Xi := E_eps (left-trivialised) equals d_A Upsilon, Upsilon := E_varpi,
      A = A0 + varpi, at x = 0; both for M_sa and M_ns.
  N3  (M_sa, M_id) Upsilon = Sh F_A + kappa *T - 1/3 h(T) + 1/2 D(T), h(T) = 2 M_eps(T^T) - M_eps(T)^T - T^M_eps(T)
      (Corollary D5) and D(T) = d_B(Sh_eps T) - Sh_eps(d_B T), which vanishes when Sh_eps is d_B-parallel (true in GU,
      where the Phi's are nabla^0-parallel; false in M3 for M_sa because Q is not A0-parallel); Upsilon != printed
      Upsilon^pr = Sh F_A + kappa *T (9.7) and Xi != d_A Upsilon^pr; control M_id = Id (Chern-Simons, parallel):
      h = D = 0 and Upsilon = Upsilon^pr.
  N4  covariance: the tau'-transformed configuration (eps h, h^-1 varpi h + h^-1 d_A0 h) has the same L_J density.
"""
import itertools
import os
import sys
import time
from math import comb, factorial

import numpy as np
import sympy as sp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from shiab_checks import Ext, perm_sign  # noqa: E402

X = sp.symbols('x0:3')
D = 3                                   # Taylor truncation degree
MONS = [m for dg in range(D + 1) for m in itertools.product(range(D + 1), repeat=3) if sum(m) == dg]
MIDX = {m: i for i, m in enumerate(MONS)}
NC = len(MONS)
TRIPLES = [(MIDX[p], MIDX[q], MIDX[tuple(a + b for a, b in zip(p, q))]) for p in MONS for q in MONS
           if sum(p) + sum(q) <= D]
ETA = [1, 1, -1]
E3 = Ext(ETA)
TAU = [sp.Matrix([[0, 1], [0, 0]]), sp.Matrix([[0, 0], [1, 0]]), sp.Matrix([[1, 0], [0, -1]])]
TAU_DUAL = [TAU[1], TAU[0], TAU[2] / 2]          # Tr(tau_i tau^j) = delta_ij
UNIT = [tuple(int(i == j) for i in range(3)) for j in range(3)]

# dual directions: gamma (i, N) |N|<=2  and  alpha (a, i, N) |N|<=1
MI2 = [m for m in MONS if sum(m) <= 2]
MI1 = [m for m in MONS if sum(m) <= 1]
DUALS = [('g', i, Nn) for i in range(3) for Nn in MI2] + [('w', a, i, Nn) for a in range(3) for i in range(3) for Nn in MI1]
DIDX = {d: k for k, d in enumerate(DUALS)}
ND = len(DUALS)


# ------------------------------------------------------------------ truncated Taylor polynomials with duals
class TP:
    """value = sum_k rows[k] t_k (t_0 = 1, t_k^2 = t_k t_l = 0), each row a Taylor polynomial (NC coefficients)."""
    __slots__ = ('a',)

    def __init__(self, a):
        self.a = a

    @staticmethod
    def zero(nd):
        return TP(np.zeros((1 + nd, NC), dtype=object) * 0)

    @staticmethod
    def const(c, nd):
        t = TP.zero(nd)
        t.a[0, 0] = c
        return t

    def __add__(self, o):
        return TP(self.a + o.a)

    def __sub__(self, o):
        return TP(self.a - o.a)

    def __neg__(self):
        return TP(-self.a)

    def scale(self, c):
        return TP(self.a * c)

    def __mul__(self, o):
        a, b = self.a, o.a
        out = np.zeros_like(a)
        dual = a.shape[0] > 1
        for p, q, r in TRIPLES:
            ap, bq = a[0, p], b[0, q]
            if ap:
                out[0, r] += ap * bq
                if dual:
                    out[1:, r] += ap * b[1:, q]
            if dual and bq:
                out[1:, r] += a[1:, p] * bq
        return TP(out)


def mzero(nd):
    return [[TP.zero(nd), TP.zero(nd)], [TP.zero(nd), TP.zero(nd)]]


def madd(A, B_):
    return [[A[i][j] + B_[i][j] for j in range(2)] for i in range(2)]


def msub(A, B_):
    return [[A[i][j] - B_[i][j] for j in range(2)] for i in range(2)]


def mscale(A, c):
    return [[A[i][j].scale(c) for j in range(2)] for i in range(2)]


def mmul(A, B_):
    return [[A[i][0] * B_[0][j] + A[i][1] * B_[1][j] for j in range(2)] for i in range(2)]


def mconst(Mc, nd):
    """constant integer matrix as TP matrix"""
    return [[TP.const(int(Mc[i, j]), nd) for j in range(2)] for i in range(2)]


def cmul(Mc, A, left=True):
    """constant integer matrix times TP matrix (cheap)"""
    out = []
    for i in range(2):
        row = []
        for j in range(2):
            if left:
                acc = A[0][j].scale(int(Mc[i, 0])) + A[1][j].scale(int(Mc[i, 1]))
            else:
                acc = A[i][0].scale(int(Mc[0, j])) + A[i][1].scale(int(Mc[1, j]))
            row.append(acc)
        out.append(row)
    return out


def mtr(A):
    return A[0][0] + A[1][1]


def mcomm(A, B_):
    return msub(mmul(A, B_), mmul(B_, A))


# ------------------------------------------------------------------ configuration
def rand_poly(rng, deg, r=1):
    return sum(int(rng.integers(-r, r + 1)) * X[0] ** e0 * X[1] ** e1 * X[2] ** e2
               for e0 in range(deg + 1) for e1 in range(deg + 1) for e2 in range(deg + 1) if e0 + e1 + e2 <= deg)


def trunc(expr, deg=D + 2):
    """drop monomials of total degree > deg (only Taylor terms of degree <= D of derivatives of order <= 2
    at x = 0 are ever used, so this is exact for every quantity computed)"""
    P = sp.Poly(sp.expand(expr), *X)
    return sum((c * X[0] ** m[0] * X[1] ** m[1] * X[2] ** m[2] for m, c in P.terms() if sum(m) <= deg), sp.Integer(0))


def tmat(Mx):
    return Mx.applyfunc(trunc)


def taylor_row(expr):
    """Taylor coefficients (degree <= D) at 0 of a sympy polynomial expression -> integer row"""
    P = sp.Poly(sp.expand(expr), *X)
    row = np.zeros(NC, dtype=object) * 0
    for mon, c in P.terms():
        if sum(mon) <= D:
            row[MIDX[tuple(mon)]] = int(c)
    return row


def dpoly(expr, Nn):
    for i, k in enumerate(Nn):
        if k:
            expr = sp.diff(expr, X[i], k)
    return expr


def binom_mi(M, Nn):
    out = 1
    for m, k in zip(M, Nn):
        out *= comb(m, k)
    return out


def leq(Nn, M):
    return all(a <= b for a, b in zip(Nn, M))


def msub_mi(M, Nn):
    return tuple(a - b for a, b in zip(M, Nn))


def tpm_from_sym(Msym, nd):
    return [[TP(np.vstack([taylor_row(Msym[i, j])] + [np.zeros(NC, dtype=object) * 0] * nd)) for j in range(2)]
            for i in range(2)]


def add_dual(Mtp, k, Mrow):
    """add the (constant-in-t) TP matrix Mrow (base rows only) into dual row k of Mtp"""
    for i in range(2):
        for j in range(2):
            Mtp[i][j].a[1 + k] += Mrow[i][j].a[0]


def build_inputs(cfg, nd_on):
    """jets of eps (|M|<=2), eps^-1 (|M|<=1), varpi_a (|M|<=1) as TP matrices; with duals if nd_on."""
    nd = ND if nd_on else 0
    eps, epsi, w = cfg['eps'], cfg['epsi'], cfg['w']
    epsJ = {M: tpm_from_sym(eps.applyfunc(lambda e: dpoly(e, M)), nd) for M in MI2}
    epsiJ = {M: tpm_from_sym(epsi.applyfunc(lambda e: dpoly(e, M)), nd) for M in MI1}
    wJ = {(a, M): tpm_from_sym(w[a].applyfunc(lambda e: dpoly(e, M)), nd) for a in range(3) for M in MI1}
    if nd_on:
        base_eps = {M: tpm_from_sym(eps.applyfunc(lambda e: dpoly(e, M)), 0) for M in MI2}
        base_epsi = {M: tpm_from_sym(epsi.applyfunc(lambda e: dpoly(e, M)), 0) for M in MI1}
        for M in MI2:                    # d_M (eps gamma) = sum_{N<=M} binom(M,N) d_{M-N} eps . d_N gamma
            for Nn in MI2:
                if leq(Nn, M):
                    for i in range(3):
                        term = cmul(TAU[i], base_eps[msub_mi(M, Nn)], left=False)
                        add_dual(epsJ[M], DIDX[('g', i, Nn)], mscale(term, binom_mi(M, Nn)))
        for M in MI1:                    # d_M (-gamma eps^-1)
            for Nn in MI1:
                if leq(Nn, M):
                    for i in range(3):
                        term = cmul(TAU[i], base_epsi[msub_mi(M, Nn)], left=True)
                        add_dual(epsiJ[M], DIDX[('g', i, Nn)], mscale(term, -binom_mi(M, Nn)))
        for a in range(3):
            for M in MI1:
                for i in range(3):
                    add_dual(wJ[(a, M)], DIDX[('w', a, i, M)], mconst(TAU[i], 0))
    return epsJ, epsiJ, wJ, nd


def lagrangian6(epsJ, epsiJ, wJ, nd, A0, Mop, kappa, want=False):
    """6 x density of (9.4); also returns pointwise objects if want."""
    O = (0, 0, 0)
    e0, ei = epsJ[O], epsiJ[O]
    A0t = [mconst(A0[a], nd) for a in range(3)]
    de = [epsJ[UNIT[a]] for a in range(3)]
    dei = [epsiJ[UNIT[b]] for b in range(3)]
    dde = {(a, b): epsJ[tuple(UNIT[a][i] + UNIT[b][i] for i in range(3))] for a in range(3) for b in range(3)}
    Bf = [madd(mmul(ei, de[a]), mmul(ei, mmul(A0t[a], e0))) for a in range(3)]
    dB = {(b, a): madd(madd(mmul(dei[b], de[a]), mmul(ei, dde[(a, b)])),
                       madd(mmul(dei[b], mmul(A0t[a], e0)), mmul(ei, mmul(A0t[a], de[b]))))
          for a in range(3) for b in range(3)}
    T = [madd(msub(wJ[(a, O)], Bf[a]), A0t[a]) for a in range(3)]
    dT = {(b, a): msub(wJ[(a, UNIT[b])], dB[(b, a)]) for a in range(3) for b in range(3)}
    pairs = [(0, 1), (0, 2), (1, 2)]
    FB = {(a, b): madd(msub(dB[(a, b)], dB[(b, a)]), mcomm(Bf[a], Bf[b])) for a, b in pairs}
    dBT = {(a, b): madd(msub(dT[(a, b)], dT[(b, a)]), madd(mcomm(Bf[a], T[b]), mcomm(T[a], Bf[b]))) for a, b in pairs}
    TT = {(a, b): mcomm(T[a], T[b]) for a, b in pairs}
    Xf = {p: madd(madd(mscale(FB[p], 6), mscale(dBT[p], 3)), mscale(TT[p], 2)) for p in pairs}

    def shiab(Y):
        return mmul(ei, mmul(Mop(mmul(e0, mmul(Y, ei))), e0))
    L = TP.zero(nd)
    for a in range(3):
        for (b, c) in pairs:
            if a in (b, c):
                continue
            L = L + mtr(mmul(T[a], shiab(Xf[(b, c)]))).scale(perm_sign((a, b, c)))
    for a in range(3):
        L = L + mtr(mmul(T[a], T[a])).scale(3 * kappa * ETA[a])
    if want:
        return L, dict(T=T, B=Bf, FB=FB, dBT=dBT, TT=TT, shiab=shiab)
    return L


def lagrangianJ6(cfg, A0, Mop, kappa):
    """6 x L_J(S), S = pi'(eps, varpi) = eps varpi eps^-1 - (d_A0 eps) eps^-1, assembled from the jets of the
    configuration (no duals) -- an independent formula: no B, no T, no tilt."""
    epsJ, epsiJ, wJ, _ = build_inputs(cfg, False)
    O = (0, 0, 0)
    e0, ei = epsJ[O], epsiJ[O]
    de = [epsJ[UNIT[a]] for a in range(3)]
    dei = [epsiJ[UNIT[b]] for b in range(3)]
    dde = {(a, b): epsJ[tuple(UNIT[a][i] + UNIT[b][i] for i in range(3))] for a in range(3) for b in range(3)}
    A0t = [mconst(A0[a], 0) for a in range(3)]
    dA0e = [madd(de[a], mcomm(A0t[a], e0)) for a in range(3)]                     # d_A0 eps
    St = [msub(mmul(e0, mmul(wJ[(a, O)], ei)), mmul(dA0e[a], ei)) for a in range(3)]
    dS = {(b, a): msub(madd(madd(mmul(de[b], mmul(wJ[(a, O)], ei)), mmul(e0, mmul(wJ[(a, UNIT[b])], ei))),
                            mmul(e0, mmul(wJ[(a, O)], dei[b]))),
                       madd(mmul(madd(dde[(a, b)], mcomm(A0t[a], de[b])), ei), mmul(dA0e[a], dei[b])))
          for a in range(3) for b in range(3)}
    pairs = [(0, 1), (0, 2), (1, 2)]
    F0 = {(a, b): mcomm(A0t[a], A0t[b]) for a, b in pairs}
    d0S = {(a, b): madd(msub(dS[(a, b)], dS[(b, a)]), madd(mcomm(A0t[a], St[b]), mcomm(St[a], A0t[b])))
           for a, b in pairs}
    SS = {(a, b): mcomm(St[a], St[b]) for a, b in pairs}
    L = TP.zero(0)
    for a in range(3):
        for (b, c) in pairs:
            if a in (b, c):
                continue
            Y = madd(madd(mscale(F0[(b, c)], 6), mscale(d0S[(b, c)], 3)), mscale(SS[(b, c)], 2))
            L = L + mtr(mmul(St[a], Mop(Y))).scale(perm_sign((a, b, c)))
    for a in range(3):
        L = L + mtr(mmul(St[a], St[a])).scale(3 * kappa * ETA[a])
    return L


def dcoef(row, M):
    """partial_M of a Taylor row at 0"""
    if sum(M) > D:
        raise ValueError
    f = 1
    for k in M:
        f *= factorial(k)
    return row[MIDX[M]] * f


def to_sym(Mtp_base):
    """value at x = 0 of a TP matrix (base row) -> sympy Matrix"""
    return sp.Matrix(2, 2, lambda i, j: sp.Integer(Mtp_base[i][j].a[0, 0]))


def main():
    t0 = time.time()
    rng = np.random.default_rng(20210401)
    Q = sp.Matrix([[1, 2], [1, -1]])
    R = sp.Matrix([[0, 1], [3, 0]])
    A0 = [sp.Matrix([[1, 1], [0, -1]]), sp.Matrix([[0, -1], [2, 0]]), sp.Matrix([[2, 0], [1, -2]])]
    kappa = 2
    Qn, Rn = np.array(Q.tolist(), dtype=object), np.array(R.tolist(), dtype=object)

    def comm_const(Cn, Y):
        return msub(cmul(Cn, Y, left=True), cmul(Cn, Y, left=False))
    Mops = {'M_sa = ad_Q ad_Q': lambda Y: comm_const(Qn, comm_const(Qn, Y)),
            'M_ns = ad_Q ad_R': lambda Y: comm_const(Qn, comm_const(Rn, Y)),
            'M_id = Id (CS)': lambda Y: Y}
    Msym = {'M_sa = ad_Q ad_Q': lambda Y: Q * (Q * Y - Y * Q) - (Q * Y - Y * Q) * Q,
            'M_ns = ad_Q ad_R': lambda Y: Q * (R * Y - Y * R) - (R * Y - Y * R) * Q}
    trsym = all(sp.expand((Y1 * Msym['M_sa = ad_Q ad_Q'](Y2)).trace() - (Y2 * Msym['M_sa = ad_Q ad_Q'](Y1)).trace()) == 0
                for Y1, Y2 in [(TAU[0], TAU[1]), (TAU[0], TAU[2]), (TAU[1], TAU[2]), (TAU[2], TAU[2])])
    trns = all(sp.expand((Y1 * Msym['M_ns = ad_Q ad_R'](Y2)).trace() - (Y2 * Msym['M_ns = ad_Q ad_R'](Y1)).trace()) == 0
               for Y1, Y2 in [(TAU[0], TAU[1]), (TAU[0], TAU[2]), (TAU[1], TAU[2])])
    F0 = [A0[a] * A0[b] - A0[b] * A0[a] for a, b in [(0, 1), (0, 2), (1, 2)]]
    print(f"setup: d = 3, eta = {ETA}, sl(2) basis tau_i; M_sa Tr-symmetric: {trsym}; M_ns Tr-symmetric: {trns}; "
          f"F0 = A0^A0 != 0: {any(f != sp.zeros(2) for f in F0)}; kappa = {kappa}; {ND} dual directions, "
          f"Taylor degree <= {D}")
    for trial in range(2):
        f, g, k = (rand_poly(rng, 5) for _ in range(3))
        U = lambda s: sp.Matrix([[1, s], [0, 1]])
        Lw = lambda s: sp.Matrix([[1, 0], [s, 1]])
        eps = tmat(U(f) * Lw(g) * U(k))
        epsi = tmat(U(-k) * Lw(-g) * U(-f))
        assert tmat(eps * epsi) == sp.eye(2)          # exact modulo degree > 5 (all that is used)
        w = [sp.Matrix(sum((rand_poly(rng, 4) * TAU[i] for i in range(3)), sp.zeros(2))) for a in range(3)]
        cfg = dict(eps=eps, epsi=epsi, w=w)
        print(f"\n--- configuration {trial}: eps(0) = {eps.subs({x: 0 for x in X}).tolist()}, "
              f"varpi(0) = {[wa.subs({x: 0 for x in X}).tolist() for wa in w]}")
        for mname, Mop in Mops.items():
            epsJ, epsiJ, wJ, nd = build_inputs(cfg, True)
            L6, P = lagrangian6(epsJ, epsiJ, wJ, nd, A0, Mop, kappa, want=True)
            LJ6 = lagrangianJ6(cfg, A0, Mop, kappa)
            fact = all(L6.a[0, i] == LJ6.a[0, i] for i in range(NC))
            print(f"  [{mname}] N1 factorisation L(eps,varpi) = L_J(pi'(eps,varpi)) for all {NC} Taylor coefficients "
                  f"(degree <= {D}) at 0: {fact}")
            assert fact
            # Euler operator
            rows = L6.a

            def P_(d):
                return rows[1 + DIDX[d]]
            Xi = sp.zeros(2)
            for i in range(3):
                Ei = sum((-1) ** sum(Nn) * dcoef(P_(('g', i, Nn)), Nn) for Nn in MI2)
                Xi += sp.Integer(Ei) * TAU_DUAL[i]
            Ut, dUt = [sp.zeros(2) for _ in range(3)], [[sp.zeros(2) for _ in range(3)] for _ in range(3)]
            for a in range(3):
                for i in range(3):
                    E0 = sum((-1) ** sum(Nn) * dcoef(P_(('w', a, i, Nn)), Nn) for Nn in MI1)
                    Ut[a] += sp.Integer(E0) * TAU_DUAL[i]
                    for j in range(3):
                        Ej = sum((-1) ** sum(Nn) * dcoef(P_(('w', a, i, Nn)), tuple(Nn[q] + UNIT[j][q] for q in range(3)))
                                 for Nn in MI1)
                        dUt[j][a] += sp.Integer(Ej) * TAU_DUAL[i]
            # Upsilon_{12} = Ut_0, Upsilon_{02} = -Ut_1, Upsilon_{01} = Ut_2  (and derivatives)
            Ups = {(1, 2): Ut[0], (0, 2): -Ut[1], (0, 1): Ut[2]}
            dUps = {j: {(1, 2): dUt[j][0], (0, 2): -dUt[j][1], (0, 1): dUt[j][2]} for j in range(3)}
            A = [A0[a] + w[a].subs({x: 0 for x in X}) for a in range(3)]

            def dA_top(U2, dU2):
                d = dU2[0][(1, 2)] - dU2[1][(0, 2)] + dU2[2][(0, 1)]
                wedge = A[0] * U2[(1, 2)] - A[1] * U2[(0, 2)] + A[2] * U2[(0, 1)]
                wedge2 = U2[(1, 2)] * A[0] - U2[(0, 2)] * A[1] + U2[(0, 1)] * A[2]
                return d + wedge - wedge2
            dAU = dA_top(Ups, dUps)
            noe = sp.simplify(Xi - dAU) == sp.zeros(2)
            print(f"  [{mname}] N2 Noether identity  Xi = d_A Upsilon  at x = 0 (exact; values x6): {noe};  "
                  f"Xi = {Xi.tolist()}, Upsilon_01 = {Ups[(0, 1)].tolist()}")
            assert noe
            nonzero = Xi != sp.zeros(2)
            print(f"      Xi != 0 off-shell (the identity is not vacuous): {nonzero}")
            if not mname.startswith('M_ns'):
                # pointwise objects at x = 0 and their 1-jets (no duals) for printed Upsilon
                epsJ0, epsiJ0, wJ0, _ = build_inputs(cfg, False)
                _, Q0 = lagrangian6(epsJ0, epsiJ0, wJ0, 0, A0, Mop, kappa, want=True)
                Tm = Q0['T']
                FA = {p: madd(madd(Q0['FB'][p], Q0['dBT'][p]), Q0['TT'][p]) for p in Q0['FB']}
                shiab = Q0['shiab']

                def star1(Tl):
                    out = {}
                    for a in range(3):
                        Ic, c = E3.star_coeff((a,))
                        out[Ic] = madd(out[Ic], mscale(Tl[a], c)) if Ic in out else mscale(Tl[a], c)
                    return out
                sT = star1(Tm)
                Ms = lambda Y: shiab(Y)        # tilted M on a TP matrix
                MT = [Ms(Tm[a]) for a in range(3)]
                h = {(a, b): msub(msub(mscale(Ms(mcomm(Tm[a], Tm[b])), 2), msub(mmul(MT[a], Tm[b]), mmul(MT[b], Tm[a]))),
                                  msub(mmul(Tm[a], MT[b]), mmul(Tm[b], MT[a]))) for a, b in FA}
                Upr = {p: madd(shiab(mscale(FA[p], 6)), mscale(sT[p], 6 * kappa)) for p in FA}
                Ud5 = {p: msub(Upr[p], mscale(h[p], 2)) for p in FA}
                # non-parallelism term D(T) = d_B(Sh T) - Sh(d_B T) (zero when Sh_eps is d_B-parallel, as in GU)
                Bm = Q0['B']
                Y = [Ms(Tm[a]) for a in range(3)]
                dY = lambda a, b: sp.Matrix(2, 2, lambda r, s: sp.Integer(dcoef(Y[b][r][s].a[0], UNIT[a])))
                DT = {(a, b): dY(a, b) - dY(b, a) + to_sym(mcomm(Bm[a], Y[b])) + to_sym(mcomm(Y[a], Bm[b]))
                      - to_sym(Ms(Q0['dBT'][(a, b)])) for a, b in FA}
                ok5 = all(to_sym(Ud5[p]) == Ups[p] for p in FA)
                okgen = all(to_sym(Ud5[p]) + 3 * DT[p] == Ups[p] for p in FA)
                okpr = all(to_sym(Upr[p]) == Ups[p] for p in FA)
                diff = max(abs(int(x)) for p in FA for x in (to_sym(Upr[p]) - Ups[p]))
                hz = all(to_sym(h[p]) == sp.zeros(2) for p in FA)
                Dz = all(DT[p] == sp.zeros(2) for p in FA)
                print(f"  [{mname}] N3 Upsilon (Euler operator) = Sh F_A + kappa *T - (1/3) h(T) + (1/2) D(T) at 0: {okgen}"
                      f"  [D(T) = d_B(Sh T) - Sh(d_B T) = 0: {Dz}; h(T) = 0: {hz}]")
                print(f"      without the D-term (Corollary D5, which assumes d_B Sh_eps = 0): {ok5}; "
                      f"= printed Sh F_A + kappa *T (9.7): {okpr} (max |difference| x6 = {diff})")
                # d_A of printed Upsilon needs its 1-jet
                dUpr = {j: {p: sp.Matrix(2, 2, lambda r, s: sp.Integer(dcoef(Upr[p][r][s].a[0], UNIT[j])))
                            for p in FA} for j in range(3)}
                Upr0 = {p: to_sym(Upr[p]) for p in FA}
                mism = sp.simplify(Xi - dA_top(Upr0, dUpr)) != sp.zeros(2)
                print(f"      Xi != d_A Upsilon^pr (the redundancy (9.6) fails for the printed Upsilon): {mism}")
            # N4 covariance under tau'
        hpoly = [rand_poly(rng, 4) for _ in range(3)]
        U = lambda s: sp.Matrix([[1, s], [0, 1]])
        Lw = lambda s: sp.Matrix([[1, 0], [s, 1]])
        h = tmat(U(hpoly[0]) * Lw(hpoly[1]) * U(hpoly[2]))
        hi = tmat(U(-hpoly[2]) * Lw(-hpoly[1]) * U(-hpoly[0]))
        A0s = [sp.Matrix(A0[a]) for a in range(3)]
        w2 = [tmat(tmat(hi * w[a]) * h + hi * (h.diff(X[a]) + A0s[a] * h - h * A0s[a])) for a in range(3)]
        cfg2 = dict(eps=tmat(eps * h), epsi=tmat(hi * epsi), w=w2)
        Mop = Mops['M_sa = ad_Q ad_Q']
        epsJ2, epsiJ2, wJ2, _ = build_inputs(cfg2, False)
        L2 = lagrangian6(epsJ2, epsiJ2, wJ2, 0, A0, Mop, kappa)
        epsJ1, epsiJ1, wJ1, _ = build_inputs(cfg, False)
        L1 = lagrangian6(epsJ1, epsiJ1, wJ1, 0, A0, Mop, kappa)
        inv = all(L1.a[0, i] == L2.a[0, i] for i in range(NC))
        print(f"  N4 tau'-invariance of the density: L(eps h, h^-1 varpi h + h^-1 d_A0 h) = L(eps, varpi) "
              f"(all Taylor coefficients, random polynomial h): {inv}")
        assert inv
    print(f"\nruntime {time.time() - t0:.0f} s")


if __name__ == '__main__':
    main()
