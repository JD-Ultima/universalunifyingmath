"""Referee cross-check of GU round 11 (formalization/geometric-unity/O-gimel-variation.md), REVIEW9-GU.

Nothing from O_gimel.py, F_jets.py or I_elliptic_reduction.py is imported.  The only repository import is
xcheck_round5.py (the round-5 referee's own Clifford algebra and Shiab templates, unmodified), used in X7 only.
Everything else is built here from the definitions:
  - own truncated Taylor jets (order 2, exact Fractions) and a textbook scalar-curvature routine at a point;
  - own Zorro metrics Z1, Z2 on Y = Met(R^n) from R-Z1/R-Z2 of F (Dy = dy - (Gamma^T y + y Gamma) dx);
  - own algebraic curvature tensors, Schur traces, boost matrices on Sym^2, the |A|^2 formula from A = -(1/2) R.y;
  - own sympy reductions of the ray-weight and smooth-weight functionals, the BG variant, the constant modes and
    the round-10 on-shell elimination.
Sections:
  X0  inputs: F5 (Lambda_A, Lambda_P), I6 (c1, c13); the potential P~ from Euler homogeneity (O4).
  X1  O2: Sym^2 Jacobian, Diff-naturality of mu_can, det G_lambda, det g_aleph with Gamma != 0, rho.
  X2  own brute force (n = 3, also n = 2): controls; F3/F4(b) on the section; Lemma O8 (degree <= 2 in eps,
      constant sigma_1/c_V, linear tr(y^-1 Ric)/c_H, quadratic != 0 off the ray, = 0 on the ray, 3-jet irrelevant,
      c_H-scaling, invariance under linear and cubic coordinate changes); the negative ray for p = q; Z2 off the
      section (O'Neill with the |A|^2 formula); a null curvature tensor with |A|^2 = 0 on the whole fibre.
  X3  Theorem O6 inputs: boost spectra on Sym^2, diagonalisability, joint kernel R g, [p, p] = k; (4,0) rotations.
  X4  Lemma O9 (Schur traces, n = 4), the eigen-formula of |A|^2, Gauss-Bonnet split; null tensors in (1,3), (2,2);
      the kernel of B -> f_B (|A|^2 as a function of sum eps eps R_ij (x) R_ij) is R * (Hodge star) for n = 4, 0 for n = 3.
  X5  Theorem O10 / Cor. O11 / Prop. O12 / Remark O14 (sympy, own derivations).
  X6  Corollary O15: the elimination, det M, and the double-root form.
  X7  Clifford (xcheck_round5 templates): O4's <T,T> and cubic term; O13(b) grade-2 identity; Psi(T ^ xi ^ T) = 0.
Exact arithmetic (Fractions, sympy rationals).  Seeds: random.Random(20261111 + k), k listed at each use; the
round-5 module is used with random.seed(20261130).  Every check is asserted.  "observation" marks checks of the
referee's own claims (REVIEW9-GU items), not of the note.
Run (from computation/geometric-unity): python3 xcheck_round11.py > xcheck_round11.out 2>&1   (~30 s on one core;
52 checks)
"""
import itertools
import os
import random
import sys
import time
from fractions import Fraction as Fr

import sympy as sp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

T0 = time.time()
NCHK = [0]
NTRUE = [0]


def P(*a):
    print(*a)
    sys.stdout.flush()


def check(label, ok):
    NCHK[0] += 1
    ok = bool(ok)
    NTRUE[0] += ok
    P(f"  [{NCHK[0]:3d}] {label}: {ok}")
    assert ok, label


def section(s):
    P(f"\n=== {s}   [{time.time() - T0:.0f} s]")


def Q(v):
    return sp.Rational(v.numerator, v.denominator) if isinstance(v, Fr) else sp.Rational(v)


# =============================================================================================== jets and curvature
class JetSpace:
    """truncated Taylor polynomials of total order <= 2 in m variables, exact Fractions.
    layout: [c] + [x_i] (m) + [x_i x_j, i <= j]"""

    def __init__(self, m):
        self.m = m
        self.quad = [(i, j) for i in range(m) for j in range(i, m)]
        self.qidx = {p: m + 1 + k for k, p in enumerate(self.quad)}
        self.L = 1 + m + len(self.quad)

    def const(self, c):
        v = [Fr(0)] * self.L
        v[0] = Fr(c)
        return v

    def var(self, i, c0=0):
        v = self.const(c0)
        v[1 + i] = Fr(1)
        return v

    def zero(self):
        return [Fr(0)] * self.L

    def add(self, a, b):
        return [p + q for p, q in zip(a, b)]

    def sub(self, a, b):
        return [p - q for p, q in zip(a, b)]

    def scale(self, a, s):
        return [p * s for p in a]

    def mul(self, a, b):
        m = self.m
        a0, b0 = a[0], b[0]
        out = [a0 * bi + b0 * ai for ai, bi in zip(a, b)]
        out[0] = a0 * b0
        la, lb = a[1:1 + m], b[1:1 + m]
        if any(la) and any(lb):
            for (i, j), k in self.qidx.items():
                out[k] += la[i] * lb[i] if i == j else la[i] * lb[j] + la[j] * lb[i]
        return out

    def d1(self, a, i):
        return a[1 + i]

    def d2(self, a, i, j):
        k = self.qidx[(min(i, j), max(i, j))]
        return 2 * a[k] if i == j else a[k]


def jsum(J, lst):
    v = J.zero()
    for w in lst:
        v = J.add(v, w)
    return v


def finv(M):
    n = len(M)
    A = [list(map(Fr, row)) + [Fr(int(i == j)) for j in range(n)] for i, row in enumerate(M)]
    for c in range(n):
        p = next(r for r in range(c, n) if A[r][c] != 0)
        A[c], A[p] = A[p], A[c]
        pv = A[c][c]
        A[c] = [v / pv for v in A[c]]
        for r in range(n):
            if r != c and A[r][c] != 0:
                f = A[r][c]
                A[r] = [v - f * w for v, w in zip(A[r], A[c])]
    return [row[n:] for row in A]


def fnull(rows):
    """exact basis of the null space of a list of Fraction rows (reduced row echelon form)"""
    A = [list(map(Fr, r)) for r in rows]
    m = len(A[0])
    piv, rk = [], 0
    for c in range(m):
        p = next((r for r in range(rk, len(A)) if A[r][c] != 0), None)
        if p is None:
            continue
        A[rk], A[p] = A[p], A[rk]
        pv = A[rk][c]
        A[rk] = [v / pv for v in A[rk]]
        for r in range(len(A)):
            if r != rk and A[r][c] != 0:
                f = A[r][c]
                A[r] = [v - f * w for v, w in zip(A[r], A[rk])]
        piv.append(c)
        rk += 1
    out = []
    for fc in (c for c in range(m) if c not in piv):
        v = [Fr(0)] * m
        v[fc] = Fr(1)
        for i, pc in enumerate(piv):
            v[pc] = -A[i][fc]
        out.append(v)
    return out


def fmat(A, B):
    return [[sum((A[i][k] * B[k][j] for k in range(len(B))), Fr(0)) for j in range(len(B[0]))] for i in range(len(A))]


def ftr(A):
    return sum((A[i][i] for i in range(len(A))), Fr(0))


def scalar_at_origin(J, G):
    """scalar curvature at 0 of the metric with order-2 jet components G[a][b]; Ric_bc = d_a Gam^a_bc - d_c Gam^a_ab
    + Gam^a_ad Gam^d_bc - Gam^a_cd Gam^d_ab (unit sphere: Ric = (n-1) g)."""
    m = J.m
    g0 = [[G[a][b][0] for b in range(m)] for a in range(m)]
    gi = finv(g0)
    dg = [[[J.d1(G[a][b], c) for b in range(m)] for a in range(m)] for c in range(m)]
    Gl = [[[(dg[b][e][c] + dg[c][e][b] - dg[e][b][c]) / 2 for c in range(m)] for b in range(m)] for e in range(m)]
    ddg = {}
    for c in range(m):
        for d in range(c, m):
            ddg[(c, d)] = [[J.d2(G[a][b], c, d) for b in range(m)] for a in range(m)]

    def dd(c, d, a, b):
        return ddg[(min(c, d), max(c, d))][a][b]
    Gam = [[[sum(gi[a][e] * Gl[e][b][c] for e in range(m)) for c in range(m)] for b in range(m)] for a in range(m)]
    dgi = []
    for d in range(m):
        tmp = [[sum(dg[d][f][h] * gi[h][e] for h in range(m)) for e in range(m)] for f in range(m)]
        dgi.append([[-sum(gi[a][f] * tmp[f][e] for f in range(m)) for e in range(m)] for a in range(m)])

    def dGam(d, a, b, c):
        s = Fr(0)
        for e in range(m):
            s += gi[a][e] * (dd(d, b, e, c) + dd(d, c, e, b) - dd(d, e, b, c)) / 2 + dgi[d][a][e] * Gl[e][b][c]
        return s
    R = Fr(0)
    trG = [sum(Gam[a][a][d] for a in range(m)) for d in range(m)]
    for b in range(m):
        for c in range(m):
            if gi[b][c] == 0:
                continue
            ric = Fr(0)
            for a in range(m):
                ric += dGam(a, a, b, c) - dGam(c, a, b, a)
            for d in range(m):
                ric += trG[d] * Gam[d][b][c]
                for a in range(m):
                    ric -= Gam[a][c][d] * Gam[d][b][a]
            R += gi[b][c] * ric
    return R


# ---- polynomials in x (dict exponent tuple -> Fraction)
def pderiv(p, i):
    out = {}
    for e, c in p.items():
        if e[i]:
            f = list(e)
            f[i] -= 1
            out[tuple(f)] = out.get(tuple(f), Fr(0)) + c * e[i]
    return {k: v for k, v in out.items() if v}


def padd(*ps):
    out = {}
    for p in ps:
        for e, c in p.items():
            out[e] = out.get(e, Fr(0)) + c
    return {k: v for k, v in out.items() if v}


def pscale(p, s):
    return {e: c * s for e, c in p.items() if c * s}


def pmul(p, q):
    out = {}
    for e, c in p.items():
        for f, d in q.items():
            k = tuple(x + y for x, y in zip(e, f))
            out[k] = out.get(k, Fr(0)) + c * d
    return {k: v for k, v in out.items() if v}


def plinsub(p, Pm, n):
    """p(P x'): substitute x_a = sum_b P_ab x'_b"""
    lin = []
    for a in range(n):
        lin.append({tuple(int(k == b) for k in range(n)): Fr(Pm[a][b]) for b in range(n) if Pm[a][b]})
    out = {}
    for e, c in p.items():
        term = {tuple([0] * n): c}
        for a in range(n):
            for _ in range(e[a]):
                term = pmul(term, lin[a])
        out = padd(out, term)
    return out


def pjet(J, p, n):
    v = J.zero()
    for e, c in p.items():
        d = sum(e)
        if d == 0:
            v[0] += c
        elif d == 1:
            v[1 + e.index(1)] += c
        elif d == 2:
            idx = [i for i in range(n) for _ in range(e[i])]
            v[J.qidx[(idx[0], idx[1])]] += c
    return v


def sym_index(n):
    return [(i, j) for i in range(n) for j in range(i, n)]


def Emat(n, A):
    i, j = A
    M = [[0] * n for _ in range(n)]
    M[i][j] = 1
    M[j][i] = 1
    return M


def mono(n, idx):
    e = [0] * n
    for z in idx:
        e[z] += 1
    return tuple(e)


def rand_quad(n, rng, r=2):
    return {(i, j): {mono(n, (a, b)): Fr(c) for (a, b) in sym_index(n) for c in [rng.randint(-r, r)] if c}
            for (i, j) in sym_index(n)}


def rand_cub(n, rng, r=2):
    return {(i, j): {mono(n, abc): Fr(c) for abc in itertools.combinations_with_replacement(range(n), 3)
                     for c in [rng.randint(-r, r)] if c} for (i, j) in sym_index(n)}


def base_metric(n, g0, quad, cub, eps):
    """g_ij(x) = g0_ij + eps (quad_ij + cub_ij)(x); dg(0) = 0"""
    g = {}
    for i in range(n):
        for j in range(n):
            key = (min(i, j), max(i, j))
            p = {tuple([0] * n): Fr(g0[i][j])} if g0[i][j] else {}
            g[(i, j)] = padd(p, pscale(quad.get(key, {}), eps), pscale(cub.get(key, {}), eps))
    return g


def base_christoffel(n, g0, g):
    """Gamma^k_{mu i}(x) up to order 2 (valid as g - g0 = O(x^2)): (g0^-1)^{kl} Gamma_{l, mu i}"""
    g0i = finv(g0)
    dg = {(c, i, j): pderiv(g[(i, j)], c) for c in range(n) for i in range(n) for j in range(n)}
    low = {(l, mu, i): pscale(padd(dg[(mu, l, i)], dg[(i, l, mu)], pscale(dg[(l, mu, i)], -1)), Fr(1, 2))
           for l in range(n) for mu in range(n) for i in range(n)}
    return {(k, mu, i): padd(*[pscale(low[(l, mu, i)], g0i[k][l]) for l in range(n) if g0i[k][l]])
            for k in range(n) for mu in range(n) for i in range(n)}


def base_curvature(n, g0, Gam):
    """R^a_{bcd}(0) = d_c Gamma^a_{db} - d_d Gamma^a_{cb}; Ric_bd = R^a_{bad}; scal = g0^{bd} Ric_bd"""
    def d0(p, c):
        return p.get(mono(n, (c,)), Fr(0))
    Rm = {}
    for a, b, c, d in itertools.product(range(n), repeat=4):
        v = d0(Gam[(a, d, b)], c) - d0(Gam[(a, c, b)], d)
        if v:
            Rm[(a, b, c, d)] = v
    Ric = [[sum((Rm.get((a, b, a, d), Fr(0)) for a in range(n)), Fr(0)) for d in range(n)] for b in range(n)]
    g0i = finv(g0)
    scal = sum(g0i[b][d] * Ric[b][d] for b in range(n) for d in range(n))
    return Rm, Ric, scal


def zorro_scalar(n, g, Gam, y0, lam, cH, cV, variant):
    """scalar curvature of Z1 / Z2 at (x = 0, y0), from R-Z1/R-Z2 of F, own construction"""
    N = n * (n + 1) // 2
    m = n + N
    J = JetSpace(m)
    SI = sym_index(n)
    Aidx = {A: k for k, A in enumerate(SI)}
    y = [[J.var(n + Aidx[(min(i, j), max(i, j))], y0[i][j]) for j in range(n)] for i in range(n)]
    D = [[J.var(n + Aidx[(min(i, j), max(i, j))]) for j in range(n)] for i in range(n)]
    Y0i = finv(y0)

    def mm(A, B):
        return [[jsum(J, [J.mul(A[i][k], B[k][j]) for k in range(n)]) for j in range(n)] for i in range(n)]
    Y0iJ = [[J.const(Y0i[i][j]) for j in range(n)] for i in range(n)]
    T1 = mm(mm(Y0iJ, D), Y0iJ)
    T2 = mm(mm(T1, D), Y0iJ)
    yinv = [[J.add(J.sub(Y0iJ[i][j], T1[i][j]), T2[i][j]) for j in range(n)] for i in range(n)]
    MA = []
    for A in SI:
        E = Emat(n, A)
        MA.append([[jsum(J, [J.scale(yinv[i][k], E[k][j]) for k in range(n) if E[k][j]]) for j in range(n)]
                   for i in range(n)])
    trA = [jsum(J, [MA[a][i][i] for i in range(n)]) for a in range(N)]
    Gm = [[None] * N for _ in range(N)]
    for a in range(N):
        for b in range(a, N):
            v = jsum(J, [J.mul(MA[a][i][k], MA[b][k][i]) for i in range(n) for k in range(n)])
            Gm[a][b] = Gm[b][a] = J.sub(v, J.scale(J.mul(trA[a], trA[b]), lam))
    GamJ = {key: pjet(J, p, n) for key, p in Gam.items()}
    K = [[jsum(J, [J.mul(GamJ[(k, mu, i)], y[k][j]) for k in range(n)]
               + [J.mul(y[i][k], GamJ[(k, mu, j)]) for k in range(n)]) for mu in range(n)] for (i, j) in SI]
    GK = [[jsum(J, [J.mul(Gm[a][b], K[b][mu]) for b in range(N)]) for mu in range(n)] for a in range(N)]
    G = [[None] * m for _ in range(m)]
    for mu in range(n):
        for nu in range(mu, n):
            hh = J.scale(y[mu][nu], cH) if variant == 'Z1' else J.scale(pjet(J, g[(mu, nu)], n), cH)
            G[mu][nu] = G[nu][mu] = J.add(hh, J.scale(jsum(J, [J.mul(K[a][mu], GK[a][nu]) for a in range(N)]), cV))
        for a in range(N):
            G[mu][n + a] = G[n + a][mu] = J.scale(GK[a][mu], -cV)
    for a in range(N):
        for b in range(N):
            G[n + a][n + b] = J.scale(Gm[a][b], cV)
    return scalar_at_origin(J, G)


def rand_y(n, g0, rng, r=2, d=3):
    """random y = P^T g0 P of the signature of g0 (P invertible, exact)"""
    while True:
        Pm = [[Fr(rng.randint(-r, r)) + (Fr(d) if i == j else 0) for j in range(n)] for i in range(n)]
        if sp.Matrix(n, n, lambda i, j: Q(Pm[i][j])).det() != 0:
            return fmat(fmat([[Pm[j][i] for j in range(n)] for i in range(n)], g0), Pm)


def sigma1(n, lam):
    return Fr(-n * (n - 1) * (n + 2), 8) - Fr(n, 4) * ((1 + lam) / (1 - n * lam) + Fr(n + 1, 2))


def sF(n):
    return Fr(-n * (n - 1) * (n + 2), 8)


def Gy(y, h, k, lam):
    yi = finv(y)
    A, B = fmat(yi, h), fmat(yi, k)
    return ftr(fmat(A, B)) - lam * ftr(A) * ftr(B)


def Rij_endo(Rm, n, i, j):
    """(R_ij)^k_l = R^k_{l i j}"""
    return [[Rm.get((k, l, i, j), Fr(0)) for l in range(n)] for k in range(n)]


def act(Rend, y):
    """R . y = -(R^T y + y R)"""
    n = len(y)
    RT = [[Rend[j][i] for j in range(n)] for i in range(n)]
    a, b = fmat(RT, y), fmat(y, Rend)
    return [[-(a[i][j] + b[i][j]) for j in range(n)] for i in range(n)]


def A2(Rm, eta, y, lam, cH, cV):
    """|A|^2 = (c_V/(4 c_H^2)) sum_{i,j} eps_i eps_j G_y(R_ij . y, R_ij . y)  (A_X Y = -(1/2)(R(X,Y).y)^vert)"""
    n = len(y)
    s = Fr(0)
    for i in range(n):
        for j in range(n):
            if i != j:
                h = act(Rij_endo(Rm, n, i, j), y)
                s += eta[i] * eta[j] * Gy(y, h, h, lam)
    return cV / (4 * cH * cH) * s


def lagrange_coeffs(xs, ys):
    """exact coefficients of the interpolating polynomial through (xs, ys)"""
    E = sp.Symbol('E')
    pol = sp.expand(sp.interpolate(list(zip([Q(x) for x in xs], [Q(v) for v in ys])), E))
    return [pol.coeff(E, k) for k in range(len(xs))]


# =============================================================================================== X0
c_t, c_u, kap, xx, yy = sp.symbols('t u kappa x y')
cH_, cV_, lam_ = sp.symbols('c_H c_V lambda')
SG = {1: '(14,0)', -1: '(7,7)'}


def x0():
    section("X0. Inputs: F5 (Lambda_A, Lambda_P), I6 (c1, c13); the potential P~ from Euler homogeneity (O4)")
    n = 4
    s1 = -sp.Rational(n * (n - 1) * (n + 2), 8) - sp.Rational(n, 4) * ((1 + lam_) / (1 - n * lam_) + sp.Rational(n + 1, 2))
    LA1 = -cH_ * s1 / (2 * cV_)
    LA1_F5 = cH_ * n / (16 * cV_) * (n ** 2 + 2 * n - 1 + 2 * (1 + lam_) / (1 - n * lam_))
    LP1 = LA1 + cH_ / (4 * cV_ * (n * lam_ - 1))
    LA2 = cH_ * n * (n - 1) * (n + 2) / (16 * cV_)
    vals = [(LA1.subs({lam_: v, cH_: 1, cV_: 1}), LP1.subs({lam_: v, cH_: 1, cV_: 1})) for v in (0, sp.Rational(1, 2), 1)]
    check(f"F5: Lambda_A(Z1) = -c_H sigma_1/(2 c_V) = F5's general-n form = (5 c_H/4 c_V)(5 - 18 lam)/(1 - 4 lam);"
          f" Lambda_P(Z1) = (3 c_H/2 c_V)(4 - 15 lam)/(1 - 4 lam); Lambda_A(Z2) = 9 c_H/(2 c_V); (A, P) at lam = 0, 1/2, 1:"
          f" {vals}; zeros 5/18 (A), 4/15 (P)",
          sp.simplify(LA1 - LA1_F5) == 0
          and sp.simplify(LA1 - sp.Rational(5, 4) * cH_ / cV_ * (5 - 18 * lam_) / (1 - 4 * lam_)) == 0
          and sp.simplify(LP1 - sp.Rational(3, 2) * cH_ / cV_ * (4 - 15 * lam_) / (1 - 4 * lam_)) == 0
          and LA2 == sp.Rational(9, 2) * cH_ / cV_
          and vals[1] == (5, sp.Rational(21, 4))
          and sp.solve(sp.numer(sp.cancel(sp.together(LA1))), lam_) == [sp.Rational(5, 18)]
          and sp.solve(sp.numer(sp.cancel(sp.together(LP1))), lam_) == [sp.Rational(4, 15)])
    inv = {}
    for chi in (1, -1):
        c1 = kap * c_t + xx * (312 * c_t ** 2 + chi * 104 * c_u ** 2) + chi * 208 * yy * c_t * c_u   # I6, transcribed
        c13 = kap * c_u + 208 * xx * c_t * c_u + yy * (104 * c_t ** 2 + chi * 312 * c_u ** 2)
        D = sp.expand(c_t * c1 + chi * c_u * c13)
        s_ = sp.Symbol('s_')
        Dd = sp.Poly(sp.expand(D.subs({c_t: s_ * c_t, c_u: s_ * c_u}, simultaneous=True)), s_)
        D2, D3 = Dd.coeff_monomial(s_ ** 2), Dd.coeff_monomial(s_ ** 3)
        Pt = sp.expand(D2 / 2 + D3 / 3)             # Euler: the degree-k part of the action is (1/k) p(T, Upsilon_k)
        a = xx * c_t + chi * yy * c_u
        ok = sp.expand(Pt - (c_t ** 2 + chi * c_u ** 2) * (kap / 2 + 104 * a)) == 0
        ok &= sp.expand(sp.diff(Pt, c_t) - c1) == 0 and sp.expand(sp.diff(Pt, c_u) - chi * c13) == 0
        check(f"{SG[chi]}: from I6's c1, c13 and Euler homogeneity, the T-part of the density is 14 P~ with"
              f" P~ = D_2/2 + D_3/3 = (t^2 + chi u^2)(kappa/2 + 104 a), a = {a}; grad P~ = (c1, chi c13); with"
              f" <Phi^1, G e(x)gamma> = tr G = -6 s_Y the curvature part is 6 a s_Y (O4)", ok)
        inv[chi] = dict(c1=c1, c13=c13, D=D, Pt=Pt, a=a)
    return inv, dict(s1=s1, LA1=LA1, LP1=LP1, LA2=LA2)


# =============================================================================================== X1
def x1():
    section("X1. Proposition O2: Sym^2 Jacobian, naturality of mu_can, det G_lambda, det g_aleph with Gamma != 0, rho")
    rng = random.Random(20261111 + 1)
    ok = True
    # (a) Jacobian of y -> P^T y P on Sym^2 (coordinates y_A, A = (i <= j))
    for n in (2, 3, 4):
        SI = sym_index(n)
        for _ in range(2):
            Pm = [[Fr(rng.randint(-3, 3), rng.randint(1, 3)) for _ in range(n)] for _ in range(n)]
            if sp.Matrix(n, n, lambda i, j: Q(Pm[i][j])).det() == 0:
                continue
            PT = [[Pm[j][i] for j in range(n)] for i in range(n)]
            Jm = []
            for A in SI:
                img = fmat(fmat(PT, [[Fr(v) for v in row] for row in Emat(n, A)]), Pm)
                Jm.append([img[i][j] for (i, j) in SI])
            detJ = sp.Matrix(len(SI), len(SI), lambda r, c: Q(Jm[c][r])).det()
            ok &= detJ == sp.Matrix(n, n, lambda i, j: Q(Pm[i][j])).det() ** (n + 1)
    # symbolic n = 2
    p = sp.symbols('p0:4')
    Pm = sp.Matrix(2, 2, p)
    ys = sp.symbols('y0:3')
    Y = sp.Matrix([[ys[0], ys[1]], [ys[1], ys[2]]])
    Yp = Pm.T * Y * Pm
    Jac = sp.Matrix([[sp.diff(Yp[i, j], v) for v in ys] for (i, j) in sym_index(2)])
    ok &= sp.expand(Jac.det() - Pm.det() ** 3) == 0
    check("O2(a): Jacobian of y -> P^T y P on Sym^2 is det(P)^(n+1) (n = 2 symbolic; n = 3, 4 random exact P)", ok)
    # naturality: (x, y) -> (A x, A^-T y A^-1) preserves |det y|^(-n/2) d^n x d^N y
    ok = True
    for n in (3, 4):
        A = [[Fr(rng.randint(-2, 2) + (3 if i == j else 0)) for j in range(n)] for i in range(n)]
        Ai = finv(A)
        AiT = [[Ai[j][i] for j in range(n)] for i in range(n)]
        SI = sym_index(n)
        Jm = []
        for E in SI:
            img = fmat(fmat(AiT, [[Fr(v) for v in row] for row in Emat(n, E)]), Ai)
            Jm.append([img[i][j] for (i, j) in SI])
        dA = sp.Matrix(n, n, lambda i, j: Q(A[i][j])).det()
        detJy = sp.Matrix(len(SI), len(SI), lambda r, c: Q(Jm[c][r])).det()
        y0 = [[Fr(rng.randint(-2, 2)) + (Fr(4) if i == j else 0) for j in range(n)] for i in range(n)]
        for i in range(n):
            for j in range(i):
                y0[i][j] = y0[j][i]
        y1 = fmat(fmat(AiT, y0), Ai)
        d0, d1 = sp.Matrix(n, n, lambda i, j: Q(y0[i][j])).det(), sp.Matrix(n, n, lambda i, j: Q(y1[i][j])).det()
        # density transforms by |dA| |det Jy| (|d1|/|d0|)^(-n/2): check its square = 1
        ok &= (dA * detJy) ** 2 * (d1 / d0) ** (-n) == 1
    check("O2(a): mu_can = |det y|^(-n/2) d^n x d^N y is invariant under (x, y) -> (A x, A^-T y A^-1) (n = 3, 4,"
          " random exact A, y)", ok)
    # (b) det G_lambda in the basis E_ii, E_ij + E_ji, lambda symbolic, signatures (n,0), (1,n-1), (2,n-2)
    ok = True
    for n in (2, 3, 4):
        SI = sym_index(n)
        N = len(SI)
        for sig in sorted({n, 1, 2 if n > 2 else 1}):
            eta = [1] * sig + [-1] * (n - sig)
            Pm = [[Fr(rng.randint(-2, 2)) + (Fr(3) if i == j else 0) for j in range(n)] for i in range(n)]
            y = fmat(fmat([[Pm[j][i] for j in range(n)] for i in range(n)],
                          [[Fr(eta[i]) if i == j else Fr(0) for j in range(n)] for i in range(n)]), Pm)
            yi = sp.Matrix(n, n, lambda i, j: Q(y[i][j])).inv()
            Es = [sp.Matrix(Emat(n, A)) for A in SI]
            Gm = sp.Matrix(N, N, lambda a, b: (yi * Es[a] * yi * Es[b]).trace()
                           - lam_ * (yi * Es[a]).trace() * (yi * Es[b]).trace())
            dety = sp.Matrix(n, n, lambda i, j: Q(y[i][j])).det()
            ok &= sp.expand(Gm.det() - (1 - n * lam_) * 2 ** (N - n) * dety ** (-(n + 1))) == 0
    check("O2(b): det G_lambda|_y = (1 - n lambda) 2^(N-n) det(y)^-(n+1) (n = 2, 3, 4; lambda symbolic; y of"
          " signature (n,0), (1,n-1), (2,n-2))", ok)
    # (c) det g_aleph with arbitrary K = Gamma^T y + y Gamma (Gamma != 0), Z1 and Z2
    ok = True
    for n, eta, lam in ((2, (1, -1), Fr(1, 3)), (3, (1, 1, 1), Fr(0)), (3, (1, -1, -1), Fr(1, 2)), (4, (1, -1, -1, -1), Fr(1, 2))):
        SI = sym_index(n)
        N = len(SI)
        cH, cV = Fr(2), Fr(3)
        consts = set()
        for rep in range(2):
            Pm = [[Fr(rng.randint(-2, 2)) + (Fr(3) if i == j else 0) for j in range(n)] for i in range(n)]
            y = fmat(fmat([[Pm[j][i] for j in range(n)] for i in range(n)],
                          [[Fr(eta[i]) if i == j else Fr(0) for j in range(n)] for i in range(n)]), Pm)
            g = [[Fr(eta[i]) if i == j else Fr(0) for j in range(n)] for i in range(n)]
            Gam = [[[Fr(rng.randint(-2, 2)) for _ in range(n)] for _ in range(n)] for _ in range(n)]   # Gam[mu][k][i]
            Kc = []
            for (i, j) in SI:
                Kc.append([sum(Gam[mu][k][i] * y[k][j] + y[i][k] * Gam[mu][k][j] for k in range(n)) for mu in range(n)])
            Gm = [[Gy(y, [[Fr(v) for v in row] for row in Emat(n, A)], [[Fr(v) for v in row] for row in Emat(n, B)], lam)
                   for B in SI] for A in SI]
            for var in ('Z1', 'Z2'):
                h = y if var == 'Z1' else g
                M = [[Fr(0)] * (n + N) for _ in range(n + N)]
                for mu in range(n):
                    for nu in range(n):
                        M[mu][nu] = cH * h[mu][nu] + cV * sum(Kc[a][mu] * Gm[a][b] * Kc[b][nu]
                                                              for a in range(N) for b in range(N))
                    for a in range(N):
                        M[mu][n + a] = M[n + a][mu] = -cV * sum(Gm[a][b] * Kc[b][mu] for b in range(N))
                for a in range(N):
                    for b in range(N):
                        M[n + a][n + b] = cV * Gm[a][b]
                dM = sp.Matrix(n + N, n + N, lambda r, c: Q(M[r][c])).det()
                dh = sp.Matrix(n, n, lambda i, j: Q(h[i][j])).det()
                dG = sp.Matrix(N, N, lambda r, c: Q(cV * Gm[r][c])).det()
                dy = sp.Matrix(n, n, lambda i, j: Q(y[i][j])).det()
                ok &= dM == Q(cH) ** n * dh * dG
                if var == 'Z1':
                    consts.add(dM * dy ** n)        # = C^2 sign: independent of y, Gamma and the base
                    ok &= dM / (dy * dG) == Q(cH) ** n * 1           # rho_Z1^2 = c_H^n det(u), det g = 1 here
                else:
                    ok &= dM / (sp.Matrix(n, n, lambda i, j: Q(g[i][j])).det() * dG) == Q(cH) ** n
        C2 = Q(cH) ** n * Q(cV) ** N * (1 - n * Q(lam)) * 2 ** (N - n)
        ok &= consts == {C2}
    check("O2(c): det g_aleph = c_H^n det(h_H) det(c_V G_lambda) at random (y, Gamma != 0), Z1 (h_H = y) and Z2"
          " (h_H = g); Z1: det g_aleph (det y)^n = c_H^n c_V^N (1 - n lam) 2^(N-n), the same for every y and Gamma"
          " (so dvol = C |det y|^(-n/2) d^n x d^N y, gimel- and Gamma-free); rho^2 = c_H^n det u (Z1), c_H^n (Z2)"
          " (n = 2, 3, 4; four signatures)", ok)


# =============================================================================================== X2
def x2():
    section("X2. Own brute force on Y = Met(R^n): controls, F3/F4(b), Lemma O8, the negative ray, Z2 off the section")
    # controls
    J = JetSpace(3)
    G = [[J.zero() for _ in range(3)] for _ in range(3)]
    for a in range(3):
        v = J.const(4)
        for i in range(3):
            v[J.qidx[(i, i)]] = Fr(-8)
        G[a][a] = v
    J2 = JetSpace(2)
    H = [[J2.zero() for _ in range(2)] for _ in range(2)]
    for a in range(2):
        v = J2.const(4)
        for i in range(2):
            v[J2.qidx[(i, i)]] = Fr(8)
        H[a][a] = v
    check("control: own scalar-curvature routine gives 6 for S^3 (4|dx|^2/(1+|x|^2)^2) and -2 for H^2",
          scalar_at_origin(J, G) == 6 and scalar_at_origin(J2, H) == -2)
    n = 3
    cH, cV = Fr(2), Fr(3)
    runs = [((1, 1, 1), Fr(0), 3), ((1, 1, -1), Fr(1, 2), 4), ((1, -1, -1), Fr(1, 2), 5)]
    for eta, lam, k in runs:
        t1 = time.time()
        rng = random.Random(20261111 + k)
        g0 = [[Fr(eta[i]) if i == j else Fr(0) for j in range(n)] for i in range(n)]
        quad, cub = rand_quad(n, rng), rand_cub(n, rng)
        # off-section y0 of the right signature
        Pm = [[Fr(rng.randint(-1, 1)) + (Fr(2) if i == j else 0) for j in range(n)] for i in range(n)]
        y0 = fmat(fmat([[Pm[j][i] for j in range(n)] for i in range(n)], g0), Pm)
        sig = (eta.count(1), eta.count(-1))
        # section (F3, F4(b)) and Lemma O8 off the section, eps in {0, ..., 4} (+ one control value)
        g1 = base_metric(n, g0, quad, cub, Fr(1))
        Gm1 = base_christoffel(n, g0, g1)
        Rm1, Ric1, scal1 = base_curvature(n, g0, Gm1)
        s_sec1 = zorro_scalar(n, g1, Gm1, g0, lam, cH, cV, 'Z1')
        s_sec2 = zorro_scalar(n, g1, Gm1, g0, lam, cH, cV, 'Z2')
        check(f"n = 3, base {sig}, lam = {lam}, (c_H, c_V) = (2, 3), random quadratic + cubic base metric: on the"
              f" section s_Y^Z1 = R/c_H + sigma_1/c_V (F4(b)) and s_Y^Z2 = R/c_H + s_F (F3), R = {scal1}",
              s_sec1 == scal1 / cH + sigma1(n, lam) / cV and s_sec2 == scal1 / cH + sF(n) / cV)
        epss = [Fr(e) for e in (0, 1, 2, 3, 4)]
        vals, vals_nc, vals_ray, vals_cH = [], [], [], []
        yr = [[2 * v for v in row] for row in g0]
        for e in epss + [Fr(-3)]:
            ge = base_metric(n, g0, quad, cub, e)
            Ge = base_christoffel(n, g0, ge)
            vals.append(zorro_scalar(n, ge, Ge, y0, lam, cH, cV, 'Z1'))
            vals_ray.append(zorro_scalar(n, ge, Ge, yr, lam, cH, cV, 'Z1'))
            vals_cH.append(zorro_scalar(n, ge, Ge, y0, lam, 2 * cH, cV, 'Z1'))
            if e in (1, 2):
                gn = base_metric(n, g0, quad, {}, e)
                vals_nc.append(zorro_scalar(n, gn, base_christoffel(n, g0, gn), y0, lam, cH, cV, 'Z1'))
        co = lagrange_coeffs(epss, vals[:5])
        cr = lagrange_coeffs(epss, vals_ray[:5])
        cc = lagrange_coeffs(epss, vals_cH[:5])
        E = sp.Symbol('E')
        polv = sum(co[k] * E ** k for k in range(3))
        y0i = finv(y0)
        trR = sum(y0i[a][b] * Ric1[b][a] for a in range(n) for b in range(n))
        ok = co[3] == 0 and co[4] == 0 and polv.subs(E, -3) == Q(vals[5])
        ok &= co[0] == Q(sigma1(n, lam) / cV) and co[1] == Q(trR / cH) and co[2] != 0
        okr = cr[3] == 0 and cr[4] == 0 and cr[2] == 0 and cr[0] == Q(sigma1(n, lam) / cV) and cr[1] == Q(scal1 / (2 * cH))
        okc = cc[3] == 0 and cc[4] == 0 and cc[0] == co[0] and cc[1] == co[1] / 2 and cc[2] == co[2] / 4
        oknc = vals_nc == vals[1:3]
        check(f"Lemma O8, n = 3, base {sig}, lam = {lam}: s_Y^Z1(eps) at y0 off the ray is the quadratic"
              f" {co[0]} + ({co[1]}) eps + ({co[2]}) eps^2 (degree-3, -4 coefficients 0; a 6th value fits):"
              f" constant sigma_1/c_V, linear tr(y0^-1 Ric_1)/c_H, quadratic != 0", ok)
        check(f"Lemma O8, base {sig}: on the ray y0 = 2 g(0) the quadratic term vanishes (q_cg = 0): s_Y = sigma_1/c_V"
              f" + R eps/(2 c_H); doubling c_H halves the linear and quarters the quadratic coefficient; dropping the cubic"
              f" terms of g changes nothing (no d^2 Gamma)   [{time.time() - t1:.0f} s]", okr and okc and oknc)
    # coordinate invariance: linear (GL) and cubic changes; s depends on the base only through Rm(0)
    rng = random.Random(20261111 + 6)
    eta, lam = (1, 1, -1), Fr(1, 2)
    g0 = [[Fr(eta[i]) if i == j else Fr(0) for j in range(n)] for i in range(n)]
    quad, cub = rand_quad(n, rng), rand_cub(n, rng)
    Pm = [[Fr(rng.randint(-1, 1)) + (Fr(2) if i == j else 0) for j in range(n)] for i in range(n)]
    y0 = fmat(fmat([[Pm[j][i] for j in range(n)] for i in range(n)], g0), Pm)
    g = base_metric(n, g0, quad, cub, Fr(1))
    s_ref = zorro_scalar(n, g, base_christoffel(n, g0, g), y0, lam, cH, cV, 'Z1')
    # linear change x = L x': g'(x') = L^T g(L x') L, y0' = L^T y0 L
    L = [[Fr(1), Fr(1), Fr(0)], [Fr(0), Fr(2), Fr(-1)], [Fr(1), Fr(0), Fr(1)]]
    LT = [[L[j][i] for j in range(n)] for i in range(n)]
    gp = {}
    for i in range(n):
        for j in range(n):
            acc = {}
            for a in range(n):
                for b in range(n):
                    if L[a][i] and L[b][j]:
                        acc = padd(acc, pscale(plinsub(g[(a, b)], L, n), L[a][i] * L[b][j]))
            gp[(i, j)] = acc
    g0p = fmat(fmat(LT, g0), L)
    y0p = fmat(fmat(LT, y0), L)
    s_lin = zorro_scalar(n, gp, base_christoffel(n, g0p, gp), y0p, lam, cH, cV, 'Z1')
    # cubic change x = x' + phi(x'), phi^a = cubic: g'_ij = g_ij + g0_aj d_i phi^a + g0_ia d_j phi^a + O(x^3)
    phi = [{mono(n, abc): Fr(rng.randint(-2, 2)) for abc in itertools.combinations_with_replacement(range(n), 3)}
           for _ in range(n)]
    gc = {}
    for i in range(n):
        for j in range(n):
            extra = {}
            for a in range(n):
                extra = padd(extra, pscale(pderiv(phi[a], i), g0[a][j]), pscale(pderiv(phi[a], j), g0[i][a]))
            gc[(i, j)] = padd({e: c for e, c in g[(i, j)].items() if sum(e) <= 2}, extra)
    Gc = base_christoffel(n, g0, gc)
    s_cub = zorro_scalar(n, gc, Gc, y0, lam, cH, cV, 'Z1')
    Rm_a, _, _ = base_curvature(n, g0, base_christoffel(n, g0, g))
    Rm_b, _, _ = base_curvature(n, g0, Gc)
    dGam_changed = any(Gc[k] != base_christoffel(n, g0, g)[k] for k in Gc)
    check("O8 premises: s_Y^Z1(0, y0) is unchanged under a linear change of coordinates (y0 -> L^T y0 L) and under a"
          " cubic change x = x' + phi(x') that alters dGamma(0) but not Rm(0) (base (2,1)): s depends on the base only"
          " through Rm(0)", s_lin == s_ref and s_cub == s_ref and Rm_a == Rm_b and dGam_changed)
    # the negative ray for p = q (n = 2, (1,1)): Z1 at y = -2 g: s_Y = R(-2g)/c_H + sigma_1/c_V = -R/(2 c_H) + ...
    n2 = 2
    rng = random.Random(20261111 + 7)
    eta2 = (1, -1)
    g02 = [[Fr(eta2[i]) if i == j else Fr(0) for j in range(n2)] for i in range(n2)]
    ok = True
    for lam2 in (Fr(1, 3), Fr(1)):
        sc2 = 0
        while sc2 == 0:                                  # a base with R(0) != 0
            q2, c2 = rand_quad(n2, rng), rand_cub(n2, rng)
            gq = base_metric(n2, g02, q2, c2, Fr(1))
            Gq = base_christoffel(n2, g02, gq)
            _, _, sc2 = base_curvature(n2, g02, Gq)
        for cval in (Fr(2), Fr(-2), Fr(-1, 3)):
            yv = [[cval * v for v in row] for row in g02]
            sv = zorro_scalar(n2, gq, Gq, yv, lam2, cH, cV, 'Z1')
            ok &= sv == sc2 / (cval * cH) + sigma1(n2, lam2) / cV
            sz2 = zorro_scalar(n2, gq, Gq, yv, lam2, cH, cV, 'Z2')
            ok &= sz2 == sc2 / cH + sF(n2) / cV
        ok &= sc2 != 0
    check("observation (R11-2): base (1,1), Z1 on the ray y = c g for c = 2, -2, -1/3: s_Y = R/(c c_H) + sigma_1/c_V,"
          " so on the negative ray (p = q) the R-coefficient changes sign; Z2: s_Y = R/c_H + s_F for every c (lam = 1/3, 1)",
          ok)
    # Z2 off the section: O'Neill s_Y = R/c_H + s_F - |A|^2 with |A|^2 = (c_V/4c_H^2) sum eps eps G(R_ij.y, R_ij.y)
    ok = True
    for eta, lam, k in (((1, 1, 1), Fr(0), 8), ((1, 1, -1), Fr(1, 2), 9)):
        rng = random.Random(20261111 + k)
        g0 = [[Fr(eta[i]) if i == j else Fr(0) for j in range(n)] for i in range(n)]
        quad, cub = rand_quad(n, rng), rand_cub(n, rng)
        g = base_metric(n, g0, quad, cub, Fr(1))
        Gm_ = base_christoffel(n, g0, g)
        Rm, _, scal = base_curvature(n, g0, Gm_)
        for _ in range(2):
            Pm = [[Fr(rng.randint(-1, 1)) + (Fr(2) if i == j else 0) for j in range(n)] for i in range(n)]
            y0 = fmat(fmat([[Pm[j][i] for j in range(n)] for i in range(n)], g0), Pm)
            s2 = zorro_scalar(n, g, Gm_, y0, lam, cH, cV, 'Z2')
            a2 = A2(Rm, eta, y0, lam, cH, cV)
            ok &= s2 == scal / cH + sF(n) / cV - a2 and a2 != 0
        # Lemma O9's eigen-formula (definite): at u = diag, |A|^2 = (c_V/4c_H^2) sum R_klij^2 (u_k/u_l + u_l/u_k - 2)
        if eta == (1, 1, 1):
            ud = [Fr(1, 2), Fr(2), Fr(3)]
            yd = [[ud[i] if i == j else Fr(0) for j in range(n)] for i in range(n)]
            Rlow = {(k_, l, i, j): Rm.get((k_, l, i, j), Fr(0)) for k_, l, i, j in itertools.product(range(n), repeat=4)}
            eig = cV / (4 * cH * cH) * sum(Rlow[(k_, l, i, j)] ** 2 * (ud[k_] / ud[l] + ud[l] / ud[k_] - 2)
                                           for k_, l, i, j in itertools.product(range(n), repeat=4))
            s2d = zorro_scalar(n, g, Gm_, yd, lam, cH, cV, 'Z2')
            ok &= s2d == scal / cH + sF(n) / cV - eig
    check("Z2 off the section (n = 3, bases (3,0) lam = 0 and (2,1) lam = 1/2, two y0 each): brute-force s_Y = R/c_H + s_F"
          " - |A|^2 with |A|^2 = (c_V/4c_H^2) sum eps_i eps_j G_y(R_ij.y, R_ij.y) != 0; Euclidean, y = diag(1/2, 2, 3):"
          " the eigen-formula of Lemma O9", ok)
    # a null curvature tensor in (2,1): R = omega (x) omega, omega = l ^ m, l null, m spacelike; |A|^2 == 0 on F_x
    eta = (1, 1, -1)
    g0 = [[Fr(eta[i]) if i == j else Fr(0) for j in range(n)] for i in range(n)]
    lup, mup = [Fr(1), Fr(0), Fr(1)], [Fr(0), Fr(1), Fr(0)]
    ll = [eta[i] * lup[i] for i in range(n)]
    ml = [eta[i] * mup[i] for i in range(n)]
    om = [[ll[a] * ml[b] - ml[a] * ll[b] for b in range(n)] for a in range(n)]
    quad = {}
    for (i, j) in sym_index(n):          # g_ij = eta_ij - (1/3) R_iajb x^a x^b with R_iajb = om_ia om_jb
        p = {}
        for a in range(n):
            for b in range(n):
                c = -Fr(1, 3) * om[i][a] * om[j][b]
                if c:
                    p[mono(n, (a, b))] = p.get(mono(n, (a, b)), Fr(0)) + c
        quad[(i, j)] = {e: c for e, c in p.items() if c}
    g = base_metric(n, g0, quad, {}, Fr(1))
    Gm_ = base_christoffel(n, g0, g)
    Rm, Ric, scal = base_curvature(n, g0, Gm_)
    Rlow = {(a, b, c, d): eta[a] * Rm.get((a, b, c, d), Fr(0)) for a, b, c, d in itertools.product(range(n), repeat=4)}
    prop = all(Rlow[(a, b, c, d)] == om[a][b] * om[c][d] for a, b, c, d in itertools.product(range(n), repeat=4))
    ok = prop and scal == 0 and any(Ric[i][j] for i in range(n) for j in range(n))
    rng = random.Random(20261111 + 10)
    for _ in range(3):
        Pm = [[Fr(rng.randint(-2, 2)) + (Fr(3) if i == j else 0) for j in range(n)] for i in range(n)]
        y0 = fmat(fmat([[Pm[j][i] for j in range(n)] for i in range(n)], g0), Pm)
        s2 = zorro_scalar(n, g, Gm_, y0, Fr(1, 2), cH, cV, 'Z2')
        ok &= A2(Rm, eta, y0, Fr(1, 2), cH, cV) == 0 and s2 == sF(n) / cV
    check("observation (R11-1): base (2,1), Rm(0) = omega (x) omega with omega = l ^ m, l null: Rm != 0, Ric = <m,m> l (x) l"
          " != 0, R = 0, and |A|^2 = 0 at three random off-section y0; brute-force Z2 s_Y = s_F there. So for Rm != 0,"
          " |A|^2 need not be unbounded (or nonzero) on the fibre in indefinite signature", ok)


# =============================================================================================== X3
def x3():
    section("X3. Theorem O6 inputs: boosts on Sym^2 (spectra, diagonalisable, joint kernel R g), [p, p] = k; (4,0)")
    lamb = sp.Symbol('L')
    for p_, q_ in ((1, 1), (2, 1), (1, 3), (3, 1), (2, 2), (2, 3)):
        n = p_ + q_
        eta = [1] * p_ + [-1] * q_
        g = sp.diag(*eta)
        SI = sym_index(n)
        N = len(SI)
        Es = [sp.Matrix(Emat(n, A)) for A in SI]

        def coords(Mm):
            return [Mm[i, j] for (i, j) in SI]

        def actmat(X):
            return sp.Matrix(N, N, lambda r, c: coords(X.T * Es[c] + Es[c] * X)[r])
        boosts = []
        for i in range(p_):
            for a in range(p_, n):
                X = sp.zeros(n)
                X[i, a] = X[a, i] = 1
                boosts.append(X)
        rots = []
        for i, j in itertools.combinations(range(n), 2):
            if (i < p_) == (j < p_):
                X = sp.zeros(n)
                X[i, j], X[j, i] = 1, -1
                rots.append(X)
        ok = all((X.T * g + g * X) == sp.zeros(n) for X in boosts + rots)
        Ms = [actmat(X) for X in boosts]
        I_ = sp.eye(N)
        for M in Ms:
            roots = set(sp.roots(M.charpoly(lamb).as_expr(), lamb).keys())
            ok &= roots <= {0, 1, -1, 2, -2}
            ok &= (M * (M - I_) * (M + I_) * (M - 2 * I_) * (M + 2 * I_)) == sp.zeros(N)
        ns = sp.Matrix.vstack(*Ms).nullspace()
        ok &= len(ns) == 1 and coords(g) == list(ns[0] / ns[0][[k for k in range(N) if ns[0][k] != 0][0]]
                                                 * coords(g)[[k for k in range(N) if ns[0][k] != 0][0]])
        # [p, p] = k
        br = [X * Y - Y * X for X, Y in itertools.combinations(boosts, 2)]
        kdim = len(rots)
        rk = sp.Matrix([list(B) for B in br]).rank() if br else 0
        ok &= rk == kdim and all(((B.T * g + g * B) == sp.zeros(n)) and B == -B.T for B in br)
        check(f"base ({p_},{q_}): {len(boosts)} elementary boosts act on Sym^2 diagonalisably with eigenvalues in"
              f" {{0, +-1, +-2}}; joint kernel = R g; [p, p] spans k (dim {kdim})", ok)
    # (4,0): the action of so(4) preserves the Frobenius form tr(y1 y2) -> skew, purely imaginary spectrum
    n = 4
    SI = sym_index(n)
    N = len(SI)
    Es = [sp.Matrix(Emat(n, A)) for A in SI]
    Gram = sp.Matrix(N, N, lambda a, b: (Es[a] * Es[b]).trace())
    ok = True
    for i, j in itertools.combinations(range(n), 2):
        X = sp.zeros(n)
        X[i, j], X[j, i] = 1, -1
        M = sp.Matrix(N, N, lambda r, c: [(X.T * Es[c] + Es[c] * X)[a, b] for (a, b) in SI][r])
        ok &= (Gram * M + (Gram * M).T) == sp.zeros(N) and Gram.det() != 0
    check("base (4,0): every rotation acts on Sym^2 skew-symmetrically for the positive form tr(y1 y2): no real nonzero"
          " eigenvalue, no boost; K compact", ok)


# =============================================================================================== X4
def curv_space(n):
    """basis of algebraic curvature tensors R_{abcd} (all-lower), as dicts, from symmetric maps on Lambda^2 + Bianchi"""
    pairs = list(itertools.combinations(range(n), 2))
    params = [(A, B) for A in range(len(pairs)) for B in range(A, len(pairs))]

    def tensor(vec):
        Rt = {}
        for (A, B), v in zip(params, vec):
            if v == 0:
                continue
            for (I, Jp) in (((A, B), (B, A)) if A != B else ((A, B),)):
                (k, l), (i, j) = pairs[I], pairs[Jp]
                for (kk, ll, sk) in ((k, l, 1), (l, k, -1)):
                    for (ii, jj, si) in ((i, j, 1), (j, i, -1)):
                        Rt[(kk, ll, ii, jj)] = Rt.get((kk, ll, ii, jj), 0) + sk * si * v
        return Rt
    units = [tensor([1 if q == p else 0 for q in range(len(params))]) for p in range(len(params))]
    rows = []
    for a, b, c, d in itertools.product(range(n), repeat=4):
        row = [Rt.get((a, b, c, d), 0) + Rt.get((a, c, d, b), 0) + Rt.get((a, d, b, c), 0) for Rt in units]
        if any(row):
            rows.append(row)
    ns = sp.Matrix(rows).nullspace()
    return [tensor(list(v)) for v in ns]


def kn(h, k, n):
    Rt = {}
    for a, b, c, d in itertools.product(range(n), repeat=4):
        v = h[a][c] * k[b][d] + h[b][d] * k[a][c] - h[a][d] * k[b][c] - h[b][c] * k[a][d]
        if v != 0:
            Rt[(a, b, c, d)] = v
    return Rt


def x4():
    section("X4. Lemma O9 (Schur traces, n = 4), the |A|^2 formula, Gauss-Bonnet; null tensors; the kernel of B -> f_B")
    n = 4
    basis = curv_space(n)
    ip = lambda R1, R2: sum(v * R2.get(k, 0) for k, v in R1.items())       # noqa: E731
    I4 = [[int(i == j) for j in range(n)] for i in range(n)]
    scal_t = [kn(I4, I4, n)]
    ric0 = []
    for i, j in sym_index(n):
        h = [[0] * n for _ in range(n)]
        if i == j:
            if i == n - 1:
                continue
            h[i][i], h[n - 1][n - 1] = 1, -1
        else:
            h[i][j] = h[j][i] = 1
        ric0.append(kn(h, I4, n))
    # Weyl: totally trace-free elements of the space
    rows = []
    for b, d in itertools.product(range(n), repeat=2):
        rows.append([sum(Rt.get((a, b, a, d), 0) for a in range(n)) for Rt in basis])
    wco = sp.Matrix(rows).nullspace()
    weyl = []
    for v in wco:
        Rt = {}
        for c, Rb in zip(v, basis):
            for k, val in Rb.items():
                Rt[k] = Rt.get(k, 0) + c * val
        weyl.append({k: val for k, val in Rt.items() if val != 0})
    dims = (len(basis), len(weyl), len(ric0), len(scal_t))
    rng = random.Random(20261111 + 11)
    ok = dims == (20, 10, 9, 1)
    ratios = []
    for _ in range(3):
        u = [Fr(rng.randint(1, 6), rng.randint(1, 6)) for _ in range(n)]
        f = {(k, l): u[k] / u[l] + u[l] / u[k] - 2 for k in range(n) for l in range(n)}
        Sf = sum(f[(k, l)] for k in range(n) for l in range(n) if k != l)
        trs = []
        for V in (weyl, ric0, scal_t):
            Gm = sp.Matrix(len(V), len(V), lambda a, b: ip(V[a], V[b]))
            Qm = sp.Matrix(len(V), len(V), lambda a, b: sum(Q(f[(k, l)]) * V[a].get((k, l, i, j), 0) * V[b].get((k, l, i, j), 0)
                                                             for k, l, i, j in itertools.product(range(n), repeat=4)))
            trs.append((Gm.inv() * Qm).trace() / len(V))
        ok &= all(t_ == Q(Sf) / 12 for t_ in trs)
        ok &= Sf == 2 * (sum(u) * sum(1 / v for v in u) - 16)
        ratios.append(trs[0])
    check(f"Lemma O9: algebraic curvature tensors in dim 4: dims (all, Weyl, Ric_0, scalar) = {dims}; for Q_u(R) ="
          f" sum R_klij^2 f_kl at three random diagonal u, tr(Q|_V)/dim V = sum_(k!=l) f_kl/12 on all three parts"
          f" ({ratios}), and sum_(k!=l) f_kl = 2 (tr u tr u^-1 - 16)", ok)
    # eigen-formula of |A|^2 against the definition G_y(R_ij.y, R_ij.y), random R and diagonal u (n = 4, Euclidean)
    ok = True
    for _ in range(2):
        coeffs = [rng.randint(-2, 2) for _ in basis]
        R = {}
        for c, Rb in zip(coeffs, basis):
            for k, v in Rb.items():
                R[k] = R.get(k, 0) + c * v
        Rup = {k: Fr(v) for k, v in R.items() if v}               # Euclidean: index position irrelevant
        u = [Fr(rng.randint(1, 5), rng.randint(1, 5)) for _ in range(n)]
        yd = [[u[i] if i == j else Fr(0) for j in range(n)] for i in range(n)]
        a2 = A2(Rup, (1,) * n, yd, Fr(1, 3), Fr(2), Fr(3))
        eig = Fr(3) / 16 * sum(Rup.get((k, l, i, j), 0) ** 2 * (u[k] / u[l] + u[l] / u[k] - 2)
                                for k, l, i, j in itertools.product(range(n), repeat=4))
        Rsq = sum(v * v for v in R.values())
        Ric = [[sum(R.get((a, b, a, d), 0) for a in range(n)) for d in range(n)] for b in range(n)]
        sc = sum(Ric[i][i] for i in range(n))
        ric0sq = sum((Ric[i][j] - (Fr(sc, 4) if i == j else 0)) ** 2 for i in range(n) for j in range(n))
        ric2 = sum(Ric[i][j] ** 2 for i in range(n) for j in range(n))
        EGB = Rsq - 4 * ric2 + sc * sc
        ok &= a2 == eig and Rsq == EGB + 4 * ric0sq
    check("O9 / O10(b): at u = diag (Euclidean, n = 4) |A|^2 = (c_V/4c_H^2) sum R_klij^2 4 sinh^2((a_k - a_l)/2) from the"
          " definition; |Rm|^2 = E_GB + 4 |Ric_0|^2 (random algebraic R)", ok)
    # null tensors in (1,3) and (2,2): |A|^2 == 0 on F_x; the kernel of B -> f_B mod constants
    for eta, lup, mup in (((1, -1, -1, -1), [1, 1, 0, 0], [0, 0, 1, 0]), ((1, 1, -1, -1), [1, 0, 1, 0], [0, 1, 0, 1])):
        ll = [eta[i] * lup[i] for i in range(n)]
        ml = [eta[i] * mup[i] for i in range(n)]
        om = [[ll[a] * ml[b] - ml[a] * ll[b] for b in range(n)] for a in range(n)]
        norm_l = sum(eta[i] * lup[i] ** 2 for i in range(n))
        lm = sum(eta[i] * lup[i] * mup[i] for i in range(n))
        Rup = {(a, b, c, d): Fr(eta[a] * om[a][b] * om[c][d]) for a, b, c, d in itertools.product(range(n), repeat=4)
               if om[a][b] * om[c][d]}
        ok = norm_l == 0 and lm == 0 and bool(Rup)
        g0 = [[Fr(eta[i]) if i == j else Fr(0) for j in range(n)] for i in range(n)]
        for _ in range(3):
            y = rand_y(n, g0, rng)
            ok &= A2(Rup, eta, y, Fr(1, 2), Fr(2), Fr(3)) == 0
        # f_B(y) = sum_{a <= b} B_ab G_y(w_a.y, w_b.y) over a basis w of so(eta); kernel of B -> f_B mod constants
        ok &= kernel_is_star(n, eta, rng)
        check(f"observation (R11-1): base {(eta.count(1), eta.count(-1))}: R = omega (x) omega (omega = l ^ m, l null) has"
              f" |A|^2 = 0 at three random y in F_x (Rm != 0); the linear map B -> f_B on Sym^2(so) (dim 21), f_B(y) ="
              f" sum B_ab G_y(w_a.y, w_b.y), has a 1-dimensional kernel modulo constants, spanned by the Hodge star"
              f" (eps_ijkl pattern). So |A|^2 == 0 on F_x iff sum eps_i eps_j R_ij (x) R_ij is a multiple of the Hodge"
              f" star; with Bianchi this forces (R^)^2 = 0 and scal = 0 (REVIEW9-GU R11-1)", ok)
    rng = random.Random(20261111 + 12)
    ok = all(kernel_is_star(3, eta, rng, expect=0) for eta in ((1, 1, 1), (1, 1, -1)))
    check("observation (R11-1): n = 3, bases (3,0), (2,1): the map B -> f_B on Sym^2(so) (dim 6) is injective modulo"
          " constants", ok)


def kernel_is_star(n, eta, rng, expect=1):
    g0 = [[Fr(eta[i]) if i == j else Fr(0) for j in range(n)] for i in range(n)]
    pairs = list(itertools.combinations(range(n), 2))
    so = []
    for i, j in pairs:
        Xm = [[Fr(0)] * n for _ in range(n)]
        Xm[i][j] = Fr(1)
        Xm[j][i] = Fr(-eta[i] * eta[j])                # X^T g + g X = 0; X_ij = eta_jj e_i ^ e_j
        so.append(Xm)
    pairsB = [(a, b) for a in range(len(so)) for b in range(a, len(so))]

    def fvec(y):
        acts = [act(Xm, y) for Xm in so]
        return [Gy(y, acts[a], acts[b], Fr(0)) * (1 if a == b else 2) for (a, b) in pairsB]
    ys_ = [rand_y(n, g0, rng) for _ in range(len(pairsB) + 6)]
    F0 = fvec(ys_[0])
    ns = fnull([[v - w for v, w in zip(fvec(y), F0)] for y in ys_[1:]])
    if len(ns) != expect:
        return False
    if expect == 0:
        return True
    v = ns[0]
    star = []
    for (a, b) in pairsB:
        (i, j), (k, l) = pairs[a], pairs[b]
        idx = (i, j, k, l)
        if len(set(idx)) < 4 or a == b:
            star.append(Fr(0))
            continue
        perm_sign = 1
        lst = list(idx)
        for p_ in range(4):
            for q_ in range(p_ + 1, 4):
                if lst[p_] > lst[q_]:
                    perm_sign = -perm_sign
        star.append(Fr(perm_sign * eta[j] * eta[l]))
    k0 = next(k for k in range(len(v)) if star[k] != 0)
    return all(v[k] * star[k0] == star[k] * v[k0] for k in range(len(v)))


# =============================================================================================== X5
def x5(inv, F5):
    section("X5. Theorem O10, Cor. O11, Prop. O12, Remark O14: own reductions (sympy)")
    s1, LA1, LP1, LA2 = F5['s1'], F5['LA1'], F5['LP1'], F5['LA2']
    Rs, tau = sp.symbols('R tau')
    w = sp.symbols('w0:3', positive=True)
    tk = sp.symbols('tau0:3')
    sFv = sp.Rational(-9) / cV_
    out = {}
    for chi in (1, -1):
        d = inv[chi]
        a, Pt = d['a'], d['Pt']
        LamO1 = LA1 - sp.Rational(7, 6) * cH_ * Pt / a
        LamO2 = LA2 - sp.Rational(7, 6) * cH_ * Pt / a
        # Z1 ray: y = e^tau g; s = tau; rho = c_H^2 e^{2 tau}; s_Y = e^-tau R/c_H + sigma_1/c_V; mu = sum w_k delta_{tau_k}
        dens1 = cH_ ** 2 * sp.exp(2 * tau) * (6 * a * (sp.exp(-tau) * Rs / cH_ + s1 / cV_) + 14 * Pt)
        I1 = sum(w[k] * dens1.subs(tau, tk[k]) for k in range(3))
        al1 = sp.diff(I1, Rs)
        be1 = sp.expand(I1 - Rs * al1)
        Lam1 = -be1 / (2 * al1)
        m1 = sum(w[k] * sp.exp(tk[k]) for k in range(3))
        m2 = sum(w[k] * sp.exp(2 * tk[k]) for k in range(3))
        ok1 = sp.simplify(Lam1 - m2 / m1 * LamO1) == 0
        # Z2 ray: s_Y = R/c_H + s_F, rho = c_H^2
        dens2 = cH_ ** 2 * (6 * a * (Rs / cH_ + sFv) + 14 * Pt)
        I2 = sum(w[k] * dens2 for k in range(3))
        al2 = sp.diff(I2, Rs)
        Lam2 = -sp.expand(I2 - Rs * al2) / (2 * al2)
        ok2 = sp.simplify(Lam2 - LamO2) == 0
        # scaling consistency (L-F6 trace part): for F(c) = int (alpha R(c g0) + beta) dvol_{c g0} at Einstein g0
        # (R(c g0) = R0/c, dvol scales by c^2), F'(1) = 0 iff R0 = 4 Lambda with Lambda = -beta/(2 alpha)
        cc, R0, al, be = sp.symbols('c R_0 alpha beta')
        Fc = (al * R0 / cc + be) * cc ** 2
        oksc = sp.solve(sp.diff(Fc, cc).subs(cc, 1), R0) == [-2 * be / al]
        # closed form and T -> 0
        e_ = sp.Symbol('e_')
        closed = LA2 - cH_ * (c_t ** 2 + chi * c_u ** 2) * (7 * kap / (12 * a) + sp.Rational(364, 3))
        okc = sp.simplify(LamO2 - closed) == 0
        lead = sp.limit(LamO2.subs({c_t: e_ * c_t, c_u: e_ * c_u}, simultaneous=True), e_, 0)
        ser = sp.series(LamO2.subs({c_t: e_ * c_t, c_u: e_ * c_u}, simultaneous=True), e_, 0, 2).removeO()
        first = sp.simplify(sp.expand(ser - LA2).coeff(e_) + sp.Rational(7, 12) * cH_ * kap * (c_t ** 2 + chi * c_u ** 2) / a)
        check(f"{SG[chi]}: ray weight mu = sum w_k delta(tau_k) (three atoms, symbolic): Z1 reduced functional"
              f" int (alpha R + beta), Lambda = -beta/(2 alpha) = (m_2/m_1)(Lambda_A - (7 c_H/6) P~/a), m_k = int e^(k tau)"
              f" d mu; Z2: Lambda = Lambda_A - (7 c_H/6) P~/a for every mu; scaling consistency F'(1) = 0 <=> R_0 = 4 Lambda;"
              f" closed form -c_H (t^2 + chi u^2)(7 kappa/(12 a) + 364/3); T -> 0: Lambda_A, first order"
              f" -(7 c_H kappa/12)(t^2 + chi u^2)/a",
              ok1 and ok2 and oksc and okc and sp.simplify(lead - LA2) == 0 and first == 0)
        # BG identification (Z2), own derivation: E^psi = E + p(delta_psi T, Upsilon); delta_psi T = -(1/2) dg . (t e(x)gamma
        # + u e(x)gamma Gamma); Upsilon = -*((x G^Y - c1 g) e(x)gamma + (y G^Y - c13 g) e(x)gamma Gamma);
        # on the ray dg_HH = c_H d(gimel), G^Y_HH = G - (1/2) c_H s_F gimel, g_HH = c_H gimel.  Coefficients of
        # <d gimel, G> and <d gimel, gimel> (per unit m and dvol): main -alpha (G + Lambda_O gimel), alpha = 6 a c_H.
        Dd = d['D']
        # pairing (per unit m, with rho = c_H^2): c_H^2 (1/(2 c_H)) <d gimel, a G - c_H (a s_F/2 + D) gimel>
        cG = -6 * a * cH_ + cH_ ** 2 / (2 * cH_) * a
        cg = -6 * a * cH_ * LamO2 + cH_ ** 2 / (2 * cH_) * (-cH_ * (a * sFv / 2 + Dd))
        LamBG = sp.simplify(cg / cG)
        okbg = sp.simplify(LamBG - (LA2 - cH_ * (14 * Pt - Dd) / (11 * a))) == 0
        okbg &= sp.simplify(LamBG - (LA2 - cH_ * (c_t ** 2 + chi * c_u ** 2) * (6 * kap / (11 * a) + 104))) == 0
        serB = sp.series(LamBG.subs({c_t: e_ * c_t, c_u: e_ * c_u}, simultaneous=True), e_, 0, 2).removeO()
        okbg &= sp.simplify(serB.coeff(e_, 0) - LA2) == 0
        okbg &= sp.simplify(serB.coeff(e_, 1) + sp.Rational(6, 11) * cH_ * kap * (c_t ** 2 + chi * c_u ** 2) / a) == 0
        check(f"{SG[chi]}: Prop. O12(b) re-derived by hand (Z2, Bourguignon-Gauduchon): G + Lambda' g = 0 with Lambda' ="
              f" Lambda_A - c_H (14 P~ - D)/(11 a) = Lambda_A - c_H (t^2 + chi u^2)(6 kappa/(11 a) + 104); equal to"
              f" Lambda_O at leading order, first order -(6/11) c_H kappa (t^2 + chi u^2)/a (not -7/12)", okbg)
        out[chi] = dict(LamO1=LamO1, LamO2=LamO2)
    # Theorem O10(b), Z2, smooth natural weight on a (4,0) base: <|A|^2>_w = (c_V/(24 c_H^2)) F' |Rm|^2 (O9),
    # |Rm|^2 = E_GB + 4 |Ric_0|^2; EL of int (alpha R + beta + gamma |Ric_0|^2): -alpha (G + Lambda g) + gamma E = 0
    m0, Fp, Ric0sq, EGB = sp.symbols("m_0 F' Ric0sq E_GB", positive=True)
    a, Pt = inv[1]['a'], inv[1]['Pt']
    Lw = cH_ ** 2 * (6 * a * (m0 * Rs / cH_ + m0 * sFv - cV_ / (24 * cH_ ** 2) * Fp * (EGB + 4 * Ric0sq)) + 14 * Pt * m0)
    alw, gaw = sp.diff(Lw, Rs), sp.diff(Lw, Ric0sq)
    bew = sp.expand(Lw - Rs * alw - Ric0sq * gaw - EGB * sp.diff(Lw, EGB))
    okw = sp.simplify(-bew / (2 * alw) - out[1]['LamO2']) == 0
    okw &= sp.simplify(-gaw / alw - cV_ * Fp / (6 * cH_ * m0)) == 0
    rng = random.Random(20261111 + 13)
    for _ in range(5):
        uu_ = [Fr(rng.randint(1, 9), rng.randint(1, 9)) for _ in range(4)]
        okw &= sum(uu_) * sum(1 / v for v in uu_) >= 16
    okw &= Fr(4) * Fr(4) == 16 and sum([Fr(3)] * 4) * sum([Fr(1, 3)] * 4) == 16
    # O11(c): Lambda_O is homogeneous of degree 1 in c_H: r Lambda_O(c_H) = Lambda_O(r c_H)
    rr = sp.Symbol('r', positive=True)
    okw &= all(sp.simplify(out[chi][k].subs(cH_, rr * cH_) - rr * out[chi][k]) == 0 for chi in (1, -1)
               for k in ('LamO1', 'LamO2'))
    check("Theorem O10(b), Z2, smooth natural weight (base (4,0)): reduced Lagrangian alpha R + beta + gamma |Ric_0|^2"
          " + (Gauss-Bonnet) with O9's average; EL G + Lambda g + kappa_w E_|Ric_0|^2 = 0, Lambda = Lambda_O(Z2) and"
          " kappa_w = c_V F'/(6 c_H m_0), F' = int (tr u tr u^-1 - 16) w >= 0 (Cauchy-Schwarz; = 0 on the ray);"
          " Cor. O11(c): r Lambda_O(c_H) = Lambda_O(r c_H) in Z1 and Z2", okw)
    # p = q: the negative ray enters m_1 with a minus sign (X2), so r_mu can be negative; the sign of Lambda flips
    a, Pt = inv[-1]['a'], inv[-1]['Pt']
    wp, wm = sp.symbols('w_p w_m', nonnegative=True)
    dens_p = cH_ ** 2 * (6 * a * (Rs / cH_ + s1 / cV_) + 14 * Pt)            # mu_+ = w_p delta_0
    dens_m = cH_ ** 2 * (6 * a * (-Rs / cH_ + s1 / cV_) + 14 * Pt)           # mu_- = w_m delta_0 on y = -g
    Itot = wp * dens_p + wm * dens_m
    alt = sp.diff(Itot, Rs)
    Lt = sp.simplify(-sp.expand(Itot - Rs * alt) / (2 * alt))
    okpq = sp.simplify(Lt - (wp + wm) / (wp - wm) * out[-1]['LamO1']) == 0
    okpq &= sp.simplify(Lt.subs({wp: 0, wm: 1}) + out[-1]['LamO1']) == 0
    okpq &= sp.simplify(alt.subs({wp: 1, wm: 1})) == 0
    check("observation (R11-2): base (2,2), Z1, mu_+ = w_p delta_0, mu_- = w_m delta_0 (both natural): Lambda ="
          " ((w_p + w_m)/(w_p - w_m)) Lambda_O; mu = delta on the negative ray gives -Lambda_O (sign reversed);"
          " w_p = w_m gives alpha = 0 (no Einstein term)", okpq)
    # Remark O14 (Z2, ray weight, y = 0): constant-mode equations d_t I = d_u I = 0 at R = 4 Lambda_O
    for chi in (1, -1):
        d = inv[chi]
        a, Pt, c1, c13 = d['a'], d['Pt'], d['c1'], d['c13']
        Lam = out[chi]['LamO2']
        dens = 6 * a * (Rs / cH_ + sFv) + 14 * Pt
        Et = sp.together((sp.diff(dens, c_t)).subs(Rs, 4 * Lam).subs(yy, 0))
        Eu = sp.factor((sp.diff(dens, c_u)).subs(Rs, 4 * Lam).subs(yy, 0))
        okm = sp.simplify(Eu - chi * 14 * c_u * (kap + 208 * c_t * xx)) == 0
        brA = sp.solve(sp.numer(Et.subs(c_t, -kap / (208 * xx))), c_u)
        u2 = sp.simplify(brA[0] ** 2)
        okm &= sp.simplify(u2 - (-chi) * (kap ** 2 / (43264 * xx ** 2) + sp.Rational(27, 728) / cV_)) == 0
        okm &= sp.simplify(Lam.subs({yy: 0, c_t: -kap / (208 * xx)}) - LA2) == 0
        brB = [v for v in sp.solve(sp.numer(Et.subs(c_u, 0)), c_t) if v != 0]          # a = x t != 0
        okm &= {sp.simplify(v ** 2) for v in brB} == {-sp.Rational(27, 728) / cV_}
        LamB = sp.simplify(sp.expand(Lam.subs({yy: 0, c_u: 0})).subs(c_t ** 2, -sp.Rational(27, 728) / cV_))
        okm &= sp.simplify(LamB - (9 * cH_ / cV_ - 7 * cH_ * kap * c_t / (12 * xx))) == 0
        check(f"{SG[chi]}: Remark O14 re-derived (Z2, ray weight, y = 0, R = 4 Lambda_O): u-equation 14 chi u (kappa +"
              f" 208 t x) = 0; branch t = -kappa/(208 x): u^2 = {-chi}(kappa^2/(43264 x^2) + 27/(728 c_V)), Lambda ="
              f" Lambda_A; branch u = 0: t^2 = -27/(728 c_V), Lambda = 9 c_H/c_V - 7 c_H kappa t/(12 x)", okm)
    # Cor. O11(d), (e): comparison with (P) of I7(a) at y = 0; P~ = 0 on a = -kappa/208 and t = +-u (7,7)
    ok = True
    for chi in (1, -1):
        d = inv[chi]
        e_ = sp.Symbol('e_')
        LP4 = sp.Symbol('L_P') - cH_ * d['c1'] / xx
        fP = sp.series(LP4.subs({c_t: e_ * c_t, c_u: e_ * c_u, yy: 0}, simultaneous=True), e_, 0, 2).removeO().coeff(e_)
        fO = sp.series(out[chi]['LamO2'].subs({c_t: e_ * c_t, c_u: e_ * c_u, yy: 0}, simultaneous=True), e_, 0, 2) \
            .removeO().coeff(e_)
        ok &= sp.simplify(fP + cH_ * kap * c_t / xx) == 0 and sp.simplify(fP - fO) != 0
        ok &= sp.expand(d['Pt'].subs(yy, 0).subs(c_t, -kap / (208 * xx))) == 0
        if chi == -1:
            ok &= sp.expand(d['Pt'].subs(c_u, c_t)) == 0 and sp.expand(d['Pt'].subs(c_u, -c_t)) == 0
    check("Cor. O11(d), (e): at y = 0 the (P) value Lambda_P - c_H c1/x has first-order term -c_H kappa t/x, different"
          " from the gimel value; P~ = 0 on t = -kappa/(208 x) (y = 0), and in (7,7) on t = +-u", ok)
    return out


# =============================================================================================== X6
def x6(inv, F5):
    section("X6. Corollary O15: Upsilon = 0 at the round-10 flat-base solutions together with the gimel-equation")
    s1, LA1 = F5['s1'], F5['LA1']
    l5 = sp.Rational(1, 5)
    s15 = s1.subs(lam_, l5)
    LA15 = sp.simplify(LA1.subs(lam_, l5))
    # flat Z1 at lambda = 1/(n+1): Ric = -(n+1)/(4 c_V) g_aleph (F6/F7); G = Ric - s/2 g, d = 14
    ricc = -sp.Rational(5, 4) / cV_
    sY = 14 * ricc
    GY = ricc - sY / 2
    ok = sp.simplify(GY - sp.Rational(15, 2) / cV_) == 0 and sp.simplify(sY - s15 / cV_) == 0 and s15 == -sp.Rational(35, 2)
    ok &= sp.simplify(LA15 - sp.Rational(35, 4) * cH_ / cV_) == 0
    check("inputs: sigma_1(1/5) = -35/2; flat Z1 at lambda = 1/5 is Einstein with Ric = -(5/(4 c_V)) g_aleph, so s_Y ="
          " -35/(2 c_V) = sigma_1/c_V and G = (15/(2 c_V)) g_aleph; Lambda_A(Z1, 1/5) = 35 c_H/(4 c_V)", ok)
    # (14,0) real form: x = -i x', t = i s, u = 0, y = 0
    xp, s = sp.symbols("x' s", real=True)
    d = inv[1]
    sub = {xx: -sp.I * xp, yy: 0, c_t: sp.I * s, c_u: 0}
    ups = sp.expand((d['c1'] - xx * GY).subs(sub) / sp.I)          # x G = c1 g  (I6, constant t, u, T_e = 0)
    ell = sp.expand((6 * d['a'] * s15 / cV_ + 14 * d['Pt']).subs(sub))   # beta / c_H^2 m_2: the flat-base gimel-equation
    ok = sp.expand(ups - (312 * xp * s ** 2 + kap * s + sp.Rational(15, 2) * xp / cV_)) == 0
    disc = sp.discriminant(312 * xp * s ** 2 + kap * s + sp.Rational(15, 2) * xp / cV_, s)
    ok &= sp.simplify(disc - (kap ** 2 - 9360 * xp ** 2 / cV_)) == 0
    sol = sp.solve([ups, ell], [s, cV_], dict=True)
    sol = [so for so in sol if so[s] != 0]
    ok &= len(sol) == 1
    so = sol[0]
    ok &= sp.simplify(so[cV_] * kap ** 2 - 12480 * xp ** 2) == 0 and sp.simplify(so[s] + kap / (416 * xp)) == 0
    ok &= sp.simplify(12480 - 9360) > 0
    detM = kap ** 2 * (kap + 224 * xp * so[s]) * (kap + 544 * xp * so[s])
    ok &= sp.simplify(detM + sp.Rational(24, 169) * kap ** 4) == 0
    ok &= sp.Rational(12480) != sp.Rational(277440, 29)
    # double-root form: beta = 0 is ell = 0 at the solution; Upsilon = 0 traced is d_t ell = 0 (6 x s_Y + 14 c1 = 0)
    lt = sp.expand((6 * xx * c_t * s15 / cV_ + 14 * d['Pt']).subs({yy: 0, c_u: 0}))
    dl = sp.expand(sp.diff(lt, c_t))
    tr_ups = sp.expand((-6 * xx * sY - 14 * d['c1']).subs({yy: 0, c_u: 0}))
    ok &= sp.expand(dl + tr_ups) == 0
    check("Cor. O15: Upsilon = 0 (x G = c1 g) is 312 x' s^2 + kappa s + 15 x'/(2 c_V) = 0, real iff c_V kappa^2 >="
          " 9360 x'^2; with the flat-base gimel-equation beta = 0 (6 a sigma_1/c_V + 14 P~ = 0) the only solution with"
          " s != 0 is c_V kappa^2 = 12480 x'^2, s = -kappa/(416 x'); there det M = kappa^2 (kappa + 224 x' s)(kappa +"
          " 544 x' s) = -(24/169) kappa^4 != 0 (and 12480 != 277440/29); equivalently the cubic l(t) has a double"
          " root (l = 0 and, by the trace of Upsilon, d_t l = 0)", ok)


# =============================================================================================== X7
def x7():
    section("X7. Clifford checks with the round-5 referee templates: O4 terms, O13(b) grade-2 identity, Psi(T^xi^T) = 0")
    import xcheck_round5 as R5
    random.seed(20261130)
    for eta in (R5.EUC, R5.E77):
        S = R5.Shiab(eta)
        cl = S.cl
        N = R5.N
        chi = 1 if eta == R5.EUC else -1
        Phi1 = [{1 << a: 1} for a in range(N)]
        Phi1G = [cl.mul({1 << a: 1}, S.G) for a in range(N)]
        G2 = cl.mul(S.G, S.G)
        okc = G2 == {0: -chi}
        # <T, T> = sum_a eta_a Tr(T_a T_a) = 14 (t^2 + chi u^2); (1/3) Psi(T^T^T) per operator
        ok = okc
        for (tv, uv) in ((1, 0), (0, 1), (2, 1), (1, -3)):
            T = [R5.add(R5.scl(Phi1[a], tv), R5.scl(Phi1G[a], uv)) for a in range(N)]
            nrm = sum(eta[a] * cl.tr(T[a], T[a]) for a in range(N))
            ok &= nrm == 14 * (tv ** 2 + chi * uv ** 2)
            WW = S.wedge(T, T)
            for op in R5.OPS:
                cub = R5.PHASE[op] * S.psi(op, T, WW)
                if op == 'W_A':
                    pred = 3 * 14 * 104 * tv * (tv ** 2 + chi * uv ** 2)
                elif op == "W'_A":
                    pred = 3 * 14 * 104 * chi * uv * (tv ** 2 + chi * uv ** 2)
                else:
                    pred = 0
                ok &= sp.simplify(cub - pred) == 0
        check(f"{R5.SIG[eta]}: Gamma^2 = -chi; for T = t Phi^1 + u Phi^1 Gamma at four (t, u): <T, T> = 14 (t^2 + chi u^2)"
              f" and (1/3) Psi(T^T^T) = 14 * 104 a (t^2 + chi u^2) per operator (W_A: a = t, W'_A: a = chi u, twisted: 0),"
              f" in the round-5 referee's normalisation (O4)", ok)
        # template identity: {g Gamma, Y Gamma} = chi {g, Y} for odd g (grade 3) and odd Y
        ok = True
        for _ in range(20):
            abc = tuple(random.sample(range(N), 3))
            Y = R5.relt(R5.ODD, 4)
            lhs = S.M("W'_A", abc, cl.mul(Y, S.G))
            rhs = R5.scl(S.M('W_A', abc, Y), chi)
            ok &= R5.add(lhs, rhs, -1) == {}
        check(f"{R5.SIG[eta]}: observation: the W'_A template on Y Gamma equals chi times the W_A template on Y for every odd"
              f" Y (20 random), so W'_A(xi ^ Phi^1 Gamma) = chi W_A(xi ^ Phi^1) exactly, in all grades", ok)
        # O13(b): grade-2 part of Shiab(xi ^ T) = a W_A(xi ^ Phi^1)_2, twisted operators no grade 2
        ok = True
        for k in (0, 5, 9):
            xi = [({0: 1} if a == k else {}) for a in range(N)]
            Dp = {op: S.dual(op, S.wedge(xi, Phi1)) for op in R5.OPS}
            DpG = {op: S.dual(op, S.wedge(xi, Phi1G)) for op in R5.OPS}
            g2 = lambda D: [R5.grade_part(X, lambda q: q == 2) for X in D]        # noqa: E731
            WA2 = g2(Dp['W_A'])
            ok &= any(WA2)
            ok &= all(not any(g2(Dp[op])) and not any(g2(DpG[op])) for op in ('W_C.AdG', "W'_C.AdG"))
            ok &= not any(g2(Dp["W'_A"])) and not any(g2(DpG['W_A']))
            ok &= all(R5.add(Xa, R5.scl(Xb, chi), -1) == {} for Xa, Xb in zip(g2(DpG["W'_A"]), WA2))
        check(f"{R5.SIG[eta]}: Prop. O13(b): for xi = e^0, e^5, e^9 the grade-2 part of Shiab(xi ^ (t Phi^1 + u Phi^1 Gamma))"
              f" is (x t + chi y u) W_A(xi ^ Phi^1)_2 (z, w and the cross terms give no grade 2), and W_A(xi ^ Phi^1)_2 != 0",
              ok)
        # O13(a): Psi(T ^ xi ^ T) = 0 for the four templates, random sparse T and xi; control Psi(alpha ^ xi ^ T) != 0
        ok = True
        nz = 0
        for _ in range(2):
            T = [R5.relt(R5.ALL, 2) for _ in range(N)]
            al = [R5.relt(R5.ALL, 2) for _ in range(N)]
            xi = [({0: random.randint(-2, 2)} if random.random() < 0.5 else {}) for _ in range(N)]
            xi = [{k: v for k, v in X.items() if v} for X in xi]
            W = S.wedge(xi, T)
            for op in R5.OPS:
                D = S.dual(op, W)
                ok &= S.pair(T, D) == 0
                nz += S.pair(al, D) != 0
        ok &= nz > 0
        check(f"{R5.SIG[eta]}: Prop. O13(a) input: Psi(T ^ xi ^ T) = 0 for the four templates on random T, xi (two samples);"
              f" Psi(alpha ^ xi ^ T) != 0 in {nz} of 8 controls", ok)


def main():
    inv, F5 = x0()
    x1()
    x2()
    x3()
    x4()
    x5(inv, F5)
    x6(inv, F5)
    x7()
    P(f"\nSUMMARY: {NCHK[0]} checks, {NTRUE[0]} True, {NCHK[0] - NTRUE[0]} False")
    P(f"runtime {time.time() - T0:.0f} s")


if __name__ == '__main__':
    main()
