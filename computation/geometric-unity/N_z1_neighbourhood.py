"""GU round 10 (N): Upsilon = 0 on a neighbourhood of a point of s(X) in the Zorro metric Z1.

Formalization: formalization/geometric-unity/N-z1-neighbourhood.md.  Builds on F (F1, F4, F6, F7), K (K1, K2, Theorem K),
J (J1, J2, J4), J5d, L (L1, L3, L5, L8).  Imports F_jets.py and L_nonradial.py (and through it J_slaved_sector.py, ...)
unmodified.  Exact arithmetic throughout: Python Fractions (also as coefficients of truncated power series, class T below),
sympy rationals.  Fixed seeds.  One core.  Every check is asserted.

The Ricci tensor of Z1 at points near a point p of s(X) is computed by brute force with F_jets.total_metric /
F_jets.ricci_at_point, using as coefficient ring the truncated polynomials T in some of the coordinates of Y (fibre
coordinates y - g(x0), and for the full 1-jet also a shift xi of the base point).  Since every quantity is a polynomial
or a convergent series in those coordinates, the truncation is exact to the stated order.

Sections (labels = results of the note)
  N0  the line equation in the variable r: coefficients from L's own `reduced` (Shiab#, W_A, both real forms) equal the
      line formulas; reparametrisation r = -1/(224 N (kappa + 224 N sigma)); mu_1 = -7/18; the trace equation and nu;
      Phi(r) has a 1/r term with nonzero residue A.
  N1  Lemma N1 (Euler field of Z1), off the section, curved bases in general coordinates: T = 0, E^flat = c ds,
      nabla E = (1/2) P_H, |E|^2 = c, Ric(E, .) = -(1/c) E^flat.
  N2  Lemma N2 (flat model): Ric = a g + b (g_H + c ds (x) ds) at points off the section; b = 0 iff lambda = 1/5; Einstein
      flat model at lambda = 1/5 (F7) with G = (15/(2 cV)) g.
  N3  Lemma N5 (fibre jets of Ric(Z1)): along the fibre the curvature corrections of Z1 equal those of Z2 (O'Neill
      A-terms) to order 2 in the HH and VV blocks, to order 1 in HV; they differ at order 3.
  N4  Lemma N6 (data at p): R(.,.)E, Ric, nabla Ric at p in closed form (brute-force 1-jet of Ric in all 14 directions).
  N5  Theorem N7 (degree-1 integrability of the line equation at p): rank tests of the linear system (*) in (v, n);
      predicted solutions; Einstein bases with sigma^2 = 10(1-5l)(13-56l)/(7 cV^2 (1-4l)^2) pass at one point and
      fail at the point (x0, 2 g(x0)) of the radial set Sigma.
  N6  Proposition N8 (lambda = 13/56 excluded by the trace equation along Sigma).
  N7  Proposition N9 / N10 (lambda = 1/5: the order-eps^2 fibre pairing <q, T'> = alpha |W|^2, alpha != 0).
  N8  Proposition N11 (flat model: the radial ODE; the trace equation kills lambda = 13/56).
  N9  Theorem N12(c) (existence: flat base, lambda = 1/5, (14,0), T_e = 0; I7(b) algebra; det M).
Run (from computation/geometric-unity): python3 N_z1_neighbourhood.py > N_z1_neighbourhood.out 2>&1
"""
import itertools
import math
import os
import random
import sys
import time
from fractions import Fraction as Fr

for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import sympy as sp  # noqa: E402
from sympy.polys.matrices import DomainMatrix  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import F_jets as FJ  # noqa: E402
import L_nonradial as L  # noqa: E402

T0 = time.time()
NCHK = [0]
n = 4
IDX = FJ.sym_index(n)
MU1 = Fr(-7, 18)


def P(*a):
    print(*a)
    sys.stdout.flush()


def check(label, ok):
    NCHK[0] += 1
    P(f"  [{NCHK[0]:3d}] {label}: {bool(ok)}")
    assert ok, label


def section(s):
    P(f"\n=== {s}   [{time.time() - T0:.0f} s]")


def sig_of(eta):
    return (eta.count(1), eta.count(-1))


# ====================================================================== truncated power series with Fraction coefficients
class Alg:
    """truncated polynomial algebra in nv variables; total degree <= K"""

    def __init__(self, nv, K):
        self.nv, self.K = nv, K
        mons = []
        for d in range(K + 1):
            for comb in itertools.combinations_with_replacement(range(nv), d):
                e = [0] * nv
                for c in comb:
                    e[c] += 1
                mons.append(tuple(e))
        self.mons = mons
        self.idx = {m: i for i, m in enumerate(mons)}
        self.deg = [sum(m) for m in mons]


class T:
    """element of Alg: sparse dict monomial-index -> Fraction"""
    __slots__ = ('c', 'A')

    def __init__(self, A, c):
        self.A = A
        self.c = c

    @staticmethod
    def const(A, v):
        v = Fr(v)
        return T(A, {0: v} if v else {})

    def _co(self, o):
        return o if isinstance(o, T) else T.const(self.A, o)

    def __add__(self, o):
        o = self._co(o)
        d = dict(self.c)
        for k, v in o.c.items():
            w = d.get(k, 0) + v
            if w:
                d[k] = w
            else:
                d.pop(k, None)
        return T(self.A, d)
    __radd__ = __add__

    def __neg__(self):
        return T(self.A, {k: -v for k, v in self.c.items()})

    def __sub__(self, o):
        return self + (-self._co(o))

    def __rsub__(self, o):
        return self._co(o) - self

    def __mul__(self, o):
        if not isinstance(o, T):
            o = Fr(o)
            return T(self.A, {k: v * o for k, v in self.c.items()}) if o else T(self.A, {})
        if len(self.c) == 1 and 0 in self.c:
            return o * self.c[0]
        if len(o.c) == 1 and 0 in o.c:
            return self * o.c[0]
        A = self.A
        d = {}
        K, deg, idx, mons = A.K, A.deg, A.idx, A.mons
        for i, a in self.c.items():
            di = deg[i]
            for j, b in o.c.items():
                if di + deg[j] <= K:
                    k = idx[tuple(p + q for p, q in zip(mons[i], mons[j]))] if (i and j) else (i or j)
                    d[k] = d.get(k, 0) + a * b
        return T(A, {k: v for k, v in d.items() if v})
    __rmul__ = __mul__

    def inv(self):
        a0 = self.c.get(0, 0)
        if a0 == 0:
            raise ZeroDivisionError('non-unit in truncated series')
        e = self * (1 / a0) - 1
        acc = T.const(self.A, 1)
        p = T.const(self.A, 1)
        for _ in range(self.A.K):
            p = p * (-e)
            acc = acc + p
        return acc * (1 / a0)

    def __truediv__(self, o):
        return self * (1 / Fr(o)) if not isinstance(o, T) else self * o.inv()

    def __rtruediv__(self, o):
        return self._co(o) * self.inv()

    def __eq__(self, o):
        if isinstance(o, T):
            return self.c == o.c
        o = Fr(o)
        return self.c == ({0: o} if o else {})

    def __ne__(self, o):
        return not self.__eq__(o)

    def __hash__(self):
        return hash(tuple(sorted(self.c.items())))

    def coef(self, m):
        return self.c.get(self.A.idx[m], Fr(0))

    def part(self, d):
        return {k: v for k, v in self.c.items() if self.A.deg[k] == d}


class FCtx:
    """drop-in for F_jets.Ctx: Fraction coefficients (A = None) or truncated series T"""

    def __init__(self, A=None):
        self.A = A
        self.zero = Fr(0) if A is None else T.const(A, 0)
        self.one = Fr(1) if A is None else T.const(A, 1)
        self.gens = ()

    def c(self, x):
        if isinstance(x, T):
            return x
        return Fr(x) if self.A is None else T.const(self.A, x)


FC = FCtx()


def mon(a, nv):
    return tuple(1 if k == a else 0 for k in range(nv))


def QQ(v):
    return sp.QQ(v.numerator, v.denominator)


# ====================================================================== curvature tensors of the base
def kn(h, k):
    return {(a, b, c, d): h[a][c] * k[b][d] + h[b][d] * k[a][c] - h[a][d] * k[b][c] - h[b][c] * k[a][d]
            for a, b, c, d in itertools.product(range(n), repeat=4)}


def radd(*Rs):
    out = {k: Fr(0) for k in itertools.product(range(n), repeat=4)}
    for coef, R in Rs:
        for k, v in R.items():
            out[k] += coef * v
    return out


def gmat(eta, s=1):
    return [[Fr(s * eta[i]) if i == j else Fr(0) for j in range(n)] for i in range(n)]


def ric_of(R, eta):
    return [[sum(Fr(eta[a]) * R[(a, b, a, d)] for a in range(n)) for d in range(n)] for b in range(n)]


def scal_of(R, eta):
    Ric = ric_of(R, eta)
    return sum(Fr(eta[a]) * Ric[a][a] for a in range(n))


def weyl_of(R, eta):
    g = gmat(eta)
    Ric = ric_of(R, eta)
    s = scal_of(R, eta)
    Pm = [[(Ric[i][j] - s / 6 * g[i][j]) / 2 for j in range(n)] for i in range(n)]
    return radd((1, R), (-1, kn(Pm, g)))


def R1(eta):
    g = gmat(eta)
    return radd((Fr(1, 2), kn(g, g)))


def rand_R(rng, eta, k=3):
    out = radd()
    for _ in range(k):
        h = [[None] * n for _ in range(n)]
        for i in range(n):
            for j in range(i, n):
                v = Fr(rng.randint(-3, 3), rng.randint(1, 2))
                h[i][j] = v
                h[j][i] = v
        out = radd((1, out), (rng.choice([1, -1]), kn(h, h)))
    return out


def norm2(R, eta):
    return sum(Fr(eta[a] * eta[b] * eta[c] * eta[d]) * R[(a, b, c, d)] ** 2
               for a, b, c, d in itertools.product(range(n), repeat=4))


def weyl_pm(W, eta):
    """W^+- = (W +- *W)/2 on the first pair (Hodge star of eta); only for eta with *^2 = 1 on 2-forms"""
    perm_sign = {}
    for p in itertools.permutations(range(4)):
        s = 1
        for i in range(4):
            for j in range(i + 1, 4):
                if p[i] > p[j]:
                    s = -s
        perm_sign[p] = s
    vol = 1
    for e in eta:
        vol *= e
    star = {}
    for a, b, c, d in itertools.product(range(n), repeat=4):
        v = Fr(0)
        for e, f in itertools.product(range(n), repeat=2):
            eps = perm_sign.get((a, b, e, f), 0)
            if eps:
                v += Fr(eps, 2) * Fr(eta[e] * eta[f]) * W[(e, f, c, d)]
        star[(a, b, c, d)] = v
    return ({k: (W[k] + star[k]) / 2 for k in W}, {k: (W[k] - star[k]) / 2 for k in W})


def normal_base(eta, R, S=None, scale=1):
    """polynomial base metric scale*(eta - 1/3 R_{mu a nu b} x^a x^b [- 1/6 S_{c mu a nu b} x^a x^b x^c]); Gamma(0) = 0"""
    Q = {}
    for a in range(n):
        for b in range(a, n):
            Q[(a, b)] = [[-Fr(1, 3) * scale * (R[(mu, a, nu, b)] + (R[(mu, b, nu, a)] if a != b else 0))
                          for nu in range(n)] for mu in range(n)]
    C = {}
    if S is not None:
        for a, b, c in itertools.combinations_with_replacement(range(n), 3):
            C[(a, b, c)] = [[-Fr(1, 6) * scale * sum(S[(r, mu, p, nu, q)] for (p, q, r) in set(itertools.permutations((a, b, c))))
                             for nu in range(n)] for mu in range(n)]
    return (tuple(scale * e for e in eta), Q, C)


def rand_cubic(rng):
    def rsym():
        Mx = [[None] * n for _ in range(n)]
        for i in range(n):
            for j in range(i, n):
                v = Fr(rng.randint(-2, 2), rng.randint(1, 3))
                Mx[i][j] = v
                Mx[j][i] = v
        return Mx
    return {(a, b, c): rsym() for a in range(n) for b in range(a, n) for c in range(b, n)}


def base_jets_at(pb, ctx, xi):
    """3-jet at x = xi of the polynomial base pb = (diag g0, Q, C) (exact: degree <= 3); xi may carry T variables"""
    g0d, Q, C = pb
    D = 3
    X = [FJ.jadd(FJ.jvar(mu, D, ctx), FJ.jconst(ctx.c(xi[mu]), D)) for mu in range(n)]
    G = [[FJ.jconst(ctx.c(Fr(g0d[i]) if i == j else 0), D) for j in range(n)] for i in range(n)]
    for (a, b), Mx in Q.items():
        xx = FJ.jmul(X[a], X[b], D)
        for i in range(n):
            for j in range(n):
                if Mx[i][j] != 0:
                    G[i][j] = FJ.jadd(G[i][j], FJ.jscale(xx, Mx[i][j]))
    for (a, b, c), Mx in C.items():
        xx = FJ.jmul(FJ.jmul(X[a], X[b], D), X[c], D)
        for i in range(n):
            for j in range(n):
                if Mx[i][j] != 0:
                    G[i][j] = FJ.jadd(G[i][j], FJ.jscale(xx, Mx[i][j]))
    return G


# ====================================================================== Z1 / Z2 near a point of s(X)
def fibre_point(ctx, y0d, A=None, offset=0):
    """y = y0 (+ fibre variables of A, numbered from offset)"""
    y0 = [[ctx.c(Fr(y0d[i]) if i == j else 0) for j in range(n)] for i in range(n)]
    if A is not None:
        for a, (i, j) in enumerate(IDX):
            e = T(A, {A.idx[mon(offset + a, A.nv)]: Fr(1)})
            y0[i][j] = y0[i][j] + e
            if i != j:
                y0[j][i] = y0[j][i] + e
    return y0


def fibre_series(pb, lam, cH, cV, K, variant):
    """Ric of Z1/Z2 at (0, g(0) + h), h the 10 fibre variables, to order K (adapted frame = coordinate frame there)"""
    gj = base_jets_at(pb, FC, [0] * n)
    base = FJ.base_geometry(gj, n, FC)
    A = Alg(10, K)
    ctx = FCtx(A)
    y0 = fibre_point(ctx, pb[0], A)
    G, m, _, Psi0 = FJ.total_metric(gj, base, n, y0, ctx.c(lam), ctx.c(cH), ctx.c(cV), variant, ctx)
    assert all(v == 0 for row in Psi0 for v in row)
    Ric, _, _ = FJ.ricci_at_point(G, m, ctx)
    return dict(A=A, ctx=ctx, base=base, Ric=Ric, y0=y0)


def dewitt_at(ctx, y0, lam):
    Yj = [[FJ.jconst(y0[i][j], 0) for j in range(n)] for i in range(n)]
    W = FJ.jinv(Yj, 0, ctx)
    Gm = FJ.dewitt_matrix(W, n, ctx.c(lam), 0, ctx)
    return [[Gm[a][b][0][()] for b in range(10)] for a in range(10)]


def point_data(pb, lam, cH, cV):
    """at p = (0, g(0)): metric, inverse, Christoffels, Riemann (Fractions) and nabla Ric (from the 1-jet of Ric in all 14
    directions: shifted base point xi and fibre variables; exact because the base is polynomial of degree <= 3)"""
    m = 14
    gj = base_jets_at(pb, FC, [0] * n)
    base = FJ.base_geometry(gj, n, FC)
    y0 = fibre_point(FC, pb[0])
    G, _, _, Psi0 = FJ.total_metric(gj, base, n, y0, lam, cH, cV, 'Z1', FC)
    assert all(v == 0 for row in Psi0 for v in row)
    Rm, g0, g0i = FJ.riemann_at_point(G, m, FC)
    _, dG, _ = FJ.metric_point_data(G, m, FC)
    _, Gam = FJ.christoffel(g0i, dG, m, FC)
    A = Alg(14, 1)
    ctx = FCtx(A)
    xi = [T(A, {A.idx[mon(k, 14)]: Fr(1)}) for k in range(n)]
    gjs = base_jets_at(pb, ctx, xi)
    bs = FJ.base_geometry(gjs, n, ctx)
    y0s = fibre_point(ctx, pb[0], A, offset=n)
    Gs, _, _, _ = FJ.total_metric(gjs, bs, n, y0s, ctx.c(lam), ctx.c(cH), ctx.c(cV), 'Z1', ctx)
    Rics, _, _ = FJ.ricci_at_point(Gs, m, ctx)
    Ric = [[Rics[a][b].coef((0,) * 14) for b in range(m)] for a in range(m)]
    dRic = [[[Rics[a][b].coef(mon(c, 14)) for b in range(m)] for a in range(m)] for c in range(m)]
    nRic = [[[dRic[c][a][b] - sum(Gam[d][c][a] * Ric[d][b] + Gam[d][c][b] * Ric[a][d] for d in range(m))
              for b in range(m)] for a in range(m)] for c in range(m)]
    return dict(g0=g0, g0i=g0i, R=Rm, Ric=Ric, nRic=nRic, base=base, Gam=Gam)


def E_vec(eta, scale=1):
    E = [Fr(0)] * 14
    for A_, (i, j) in enumerate(IDX):
        if i == j:
            E[n + A_] = Fr(scale * eta[i])
    return E


def star_system(D, mu1=MU1):
    """(*)  R(X,Y)v - mu1[<v,X> Ric^Y - <v,Y> Ric^X] - [<n,X>Y - <n,Y>X] = (nabla_X Ric^)Y - (nabla_Y Ric^)X,
    unknowns (v, n) = nabla r(p)/(mu1 r0), nabla nu(p)/(mu1 r0); rows (a<b, output d)"""
    m = 14
    g0, g0i, R, Ric, nRic = D['g0'], D['g0i'], D['R'], D['Ric'], D['nRic']
    Rh = [[sum(g0i[d][e] * Ric[e][b] for e in range(m)) for b in range(m)] for d in range(m)]
    nRh = [[[sum(g0i[d][e] * nRic[c][e][b] for e in range(m)) for b in range(m)] for d in range(m)] for c in range(m)]
    rows, rhs = [], []
    for a in range(m):
        for b in range(a + 1, m):
            for d in range(m):
                row = [Fr(0)] * (2 * m)
                for c in range(m):
                    row[c] += sum(g0i[d][e] * R[(e, c, a, b)] for e in range(m) if g0i[d][e] != 0)
                    row[c] -= mu1 * (g0[a][c] * Rh[d][b] - g0[b][c] * Rh[d][a])
                    row[m + c] -= (g0[a][c] * (d == b) - g0[b][c] * (d == a))
                rows.append(row)
                rhs.append(nRh[a][d][b] - nRh[b][d][a])
    return rows, rhs


def rank_test(rows, rhs):
    Mx = DomainMatrix([[QQ(v) for v in r] for r in rows], (len(rows), len(rows[0])), sp.QQ)
    Ma = DomainMatrix([[QQ(v) for v in r] + [QQ(b)] for r, b in zip(rows, rhs)], (len(rows), len(rows[0]) + 1), sp.QQ)
    return Mx.rank(), Ma.rank()


def consts(lam, cV):
    c = 4 * cV * (1 - 4 * lam)
    b = (1 - 5 * lam) / (cV * (1 - 4 * lam))
    bS = Fr(-5, 4) / cV
    return c, b, bS


# ====================================================================== N0
def n0():
    section("N0. The line equation in r (L5): coefficients, reparametrisation, mu_1, trace equation, Phi")
    k_, N_, s_, r_ = sp.symbols('kappa N sigma r')
    Ht = 576 * N_ / (k_ + 224 * N_ * s_)
    Qt = -258048 * N_ ** 2 / (k_ + 224 * N_ * s_) ** 2
    Qpt = 78336 * N_ ** 2 / (k_ + 224 * N_ * s_) ** 2
    ct = k_ * s_ + 312 * N_ * s_ ** 2
    ok = True
    for eps in (1, -1):
        for nm, sub in (('Shiab#', {L.x: 1, L.y: 0, L.z: 1, L.w: 0}), ('W_A', {L.x: 1, L.y: 0, L.z: 0, L.w: 0})):
            rr = L.reduced(L.subs_tables(L.tables(eps), sub))
            vals = [sp.sympify(e).subs(L.u, 0) for e in (rr['H'][0, 0], rr['Q'][('E1', 't', 't')], rr['Qp'][('E1', 't', 't')])]
            ref = [e.subs({N_: 1, s_: L.t, k_: L.kap}) for e in (Ht, Qt, Qpt)]
            ok &= all(sp.cancel(a - b) == 0 for a, b in zip(vals, ref))
    c1 = (L.kap * L.t + L.x * (312 * L.t ** 2 + 104 * L.u ** 2) + 208 * L.y * L.t * L.u).subs({L.x: 1, L.y: 0, L.u: 0})
    ok &= sp.expand(c1 - ct.subs({N_: 1, s_: L.t, k_: L.kap})) == 0
    check("L's `reduced` (Shiab# and W_A, (14,0) and (7,7)) at u = 0: Hessian 576/(kappa + 224 t), gradient coefficients"
          " -258048/(..)^2, 78336/(..)^2, c1 = kappa t + 312 t^2: the line formulas of L5 at N = 1", ok)
    Z = -1 / (224 * N_ * r_)
    sig_r = (Z - k_) / (224 * N_)
    sp1, sp2 = sp.diff(sig_r, r_), sp.diff(sig_r, r_, 2)
    sub = {s_: sig_r}
    hess = sp.simplify((Ht * sp1).subs(sub))
    grad = sp.simplify((Ht * sp2 + Qt * sp1 ** 2).subs(sub))
    trgr = sp.simplify((-Ht * sp2 + Qpt * sp1 ** 2).subs(sub))
    check(f"with r = -1/(224 N (kappa + 224 N sigma)): H~ dsigma/dr = {hess}, the dr (x) dr coefficient H~ sigma'' + Q~ sigma'^2 ="
          f" {grad}, and the |dr|^2 g coefficient -H~ sigma'' + Q~' sigma'^2 = {trgr}",
          sp.simplify(hess + sp.Rational(18, 7) / r_) == 0 and grad == 0 and sp.simplify(trgr + sp.Rational(351, 98) / r_ ** 2) == 0)
    # E~ = -G + c~ g - 18/(7r)(Hess r - Box r g) - 351/(98 r^2)|dr|^2 g; with Hess r = mu1 r Ric + nu g, d = 14
    mu1, s_Y, nu, box, dr2 = sp.symbols('mu1 s_Y nu Box dr2')
    Phi = sp.Rational(7, 18) * r_ * ct.subs(sub)
    # traceless part: -Ric_0 - 18/(7r) (mu1 r Ric_0) = 0  -> mu1 = -7/18
    check("traceless part: (Hess r)_0 = mu_1 r (Ric)_0 with mu_1 = -7/18 (L8 step 5)",
          sp.solve(-1 - sp.Rational(18, 7) * mu1, mu1) == [sp.Rational(-7, 18)])
    # trace: 6 s_Y + 14 c~ - 18/(7r)(box - 14 box) - 14*351/(98 r^2) dr2 = 0 ; box = mu1 r s_Y + 14 nu
    tr = 6 * s_Y + 14 * ct.subs(sub) + sp.Rational(18, 7) / r_ * 13 * box - 14 * sp.Rational(351, 98) / r_ ** 2 * dr2
    tr = tr.subs(box, sp.Rational(-7, 18) * r_ * s_Y + 14 * nu)
    nus = sp.solve(tr, nu)[0]
    ref = sp.Rational(7, 468) * r_ * s_Y - Phi / 13 + sp.Rational(3, 28) * dr2 / r_
    check("trace part: nu = (7/468) r s_Y - Phi(r)/13 + (3/28)|dr|^2/r with Phi(r) = (7/18) r c~(sigma(r))",
          sp.simplify(nus - ref) == 0)
    rPhi = sp.expand(sp.cancel(r_ * Phi))
    A0 = rPhi.subs(r_, 0)
    check(f"r Phi(r) is a polynomial in r with constant term A = {A0} != 0 (N != 0)",
          sp.denom(sp.together(rPhi)).free_symbols <= {N_} and sp.simplify(A0 - sp.Rational(7 * 312, 18) / (224 ** 4 * N_ ** 3)) == 0)
    return dict(Phi=Phi, r=r_, N=N_, k=k_)


# ====================================================================== N1
def rand_general_base(rng, eta):
    def rsym(s=2):
        Mx = [[None] * n for _ in range(n)]
        for i in range(n):
            for j in range(i, n):
                v = Fr(rng.randint(-s, s), rng.randint(1, 3))
                Mx[i][j] = v
                Mx[j][i] = v
        return Mx
    g0 = gmat(eta)
    lin = [rsym() for _ in range(n)]
    quad = {(a, b): rsym() for a in range(n) for b in range(a, n)}
    cub = {(a, b, c): rsym() for a in range(n) for b in range(a, n) for c in range(b, n)}
    return FJ.base_from_coeffs(g0, lin, quad, cub, n, FC)


def n1():
    section("N1. Lemma N1 (Euler field of Z1) at points OFF the section, curved bases in general coordinates")
    m = 14
    for eta, lam, seed in (((1, 1, 1, 1), Fr(0), 11), ((1, -1, -1, -1), Fr(1, 2), 12), ((1, 1, -1, -1), Fr(2, 3), 13)):
        t1 = time.time()
        rng = random.Random(seed)
        gj = rand_general_base(rng, eta)
        base = FJ.base_geometry(gj, n, FC)
        cH, cV = Fr(2), Fr(3)
        y0 = [[Fr(eta[i]) * (i == j) + (Fr(rng.randint(-2, 2), 5) if i <= j else 0) for j in range(n)] for i in range(n)]
        for i in range(n):
            for j in range(i):
                y0[i][j] = y0[j][i]
        G, _, _, Psi0 = FJ.total_metric(gj, base, n, y0, lam, cH, cV, 'Z1', FC)
        g0, dG, _ = FJ.metric_point_data(G, m, FC)
        g0i = FJ.finv(g0, FC)
        _, Gam = FJ.christoffel(g0i, dG, m, FC)
        Ric, _, _ = FJ.ricci_at_point(G, m, FC)
        c, _, _ = consts(lam, cV)
        E = [Fr(0)] * n + [y0[i][j] for (i, j) in IDX]
        okT = all(Gam[mu][n + A_][n + B] == 0 for mu in range(n) for A_ in range(10) for B in range(10))
        PH = [[Fr(0)] * m for _ in range(m)]          # PH[a][b] = (P_H)^a_b
        for mu in range(n):
            PH[mu][mu] = Fr(1)
            for A_ in range(10):
                PH[n + A_][mu] = Psi0[A_][mu]
        nE = [[(Fr(1) if (a == b and a >= n) else Fr(0)) + sum(Gam[a][b][cc] * E[cc] for cc in range(m)) for b in range(m)]
              for a in range(m)]
        okN = all(nE[a][b] == PH[a][b] / 2 for a in range(m) for b in range(m))
        # ds, s = 1/4 log|det(g(x)^-1 y)|: d_A s = 1/4 tr(y^-1 E_A), d_mu s = -1/4 tr(g^-1 d_mu g)
        Ym = sp.Matrix(n, n, lambda i, j: sp.Rational(y0[i][j]))
        Yi = Ym.inv()
        gm0 = sp.Matrix(n, n, lambda i, j: sp.Rational(FJ.jval(gj[i][j], FC)))
        gmi = gm0.inv()
        ds = []
        for mu in range(n):
            dg = sp.Matrix(n, n, lambda i, j: sp.Rational(FJ.jgrad(gj[i][j], mu, FC)))
            ds.append(-sp.Rational(1, 4) * (gmi * dg).trace())
        for (i, j) in IDX:
            Em = sp.zeros(n, n)
            Em[i, j] = 1
            Em[j, i] = 1
            ds.append(sp.Rational(1, 4) * (Yi * Em).trace())
        Eflat = [sum(g0[a][b] * E[b] for b in range(m)) for a in range(m)]
        okF = all(sp.Rational(Eflat[a]) == sp.Rational(c) * ds[a] for a in range(m))
        EE = sum(Eflat[a] * E[a] for a in range(m))
        RicE = [sum(Ric[a][b] * E[b] for b in range(m)) for a in range(m)]
        okR = all(RicE[a] == -Eflat[a] / c for a in range(m))
        check(f"base {sig_of(eta)}, lambda = {lam}, cH = 2, cV = 3, y0 != g(x0), Gamma(x0) != 0: fibres totally geodesic"
              f" ({okT}); E^flat = c ds ({okF}); nabla E = (1/2) P_H ({okN}); |E|^2 = c = 4 cV (1 - 4 lambda) = {c} ({EE == c});"
              f" Ric(E, .) = -(1/c) E^flat ({okR})  [{time.time() - t1:.0f} s]", okT and okF and okN and EE == c and okR)


# ====================================================================== N2
def n2():
    section("N2. Lemma N2 (flat model): Ric = a g + b (g_H + c ds (x) ds) off the section; b = 0 iff lambda = 1/5")
    m = 14
    eta = (1, 1, 1, 1)
    rng = random.Random(21)
    zero = {(a, b): [[Fr(0)] * n for _ in range(n)] for a in range(n) for b in range(a, n)}
    gj = FJ.base_from_coeffs(gmat(eta), [[[Fr(0)] * n for _ in range(n)] for _ in range(n)], zero, {}, n, FC)
    base = FJ.base_geometry(gj, n, FC)
    y0 = [[Fr(eta[i]) * (i == j) + (Fr(rng.randint(-2, 2), 5) if i <= j else 0) for j in range(n)] for i in range(n)]
    for i in range(n):
        for j in range(i):
            y0[i][j] = y0[j][i]
    Ym = sp.Matrix(n, n, lambda i, j: sp.Rational(y0[i][j]))
    Yi = Ym.inv()
    ds = [Fr(0)] * n
    for (i, j) in IDX:
        Em = sp.zeros(n, n)
        Em[i, j] = 1
        Em[j, i] = 1
        v = sp.Rational(1, 4) * (Yi * Em).trace()
        ds.append(Fr(int(v.p), int(v.q)))
    cH, cV = Fr(2), Fr(3)
    for lam in (Fr(0), Fr(1, 5), Fr(13, 56), Fr(1, 2), Fr(1)):
        G, _, _, _ = FJ.total_metric(gj, base, n, y0, lam, cH, cV, 'Z1', FC)
        Ric, g0, _ = FJ.ricci_at_point(G, m, FC)
        c, b, _ = consts(lam, cV)
        a = Fr(-5, 4) / cV
        ok = all(Ric[p][q] == a * g0[p][q] + b * ((cH * y0[p][q] if p < n and q < n else 0) + c * ds[p] * ds[q])
                 for p in range(m) for q in range(m))
        sY = sum(FJ.finv(g0, FC)[p][q] * Ric[p][q] for p in range(m) for q in range(m))
        extra = ''
        if lam == Fr(1, 5):
            Gt = [[Ric[p][q] - sY / 2 * g0[p][q] for q in range(m)] for p in range(m)]
            okG = all(Gt[p][q] == Fr(15, 2) / cV * g0[p][q] for p in range(m) for q in range(m))
            extra = f"; Einstein, Ric = -(5/(4 cV)) g and G = (15/(2 cV)) g ({okG})"
            ok &= okG and b == 0
        check(f"flat base (4,0), y0 != eta, lambda = {lam}: Ric(Z1) = a g + b (g_H + c ds (x) ds), a = -5/(4 cV), b = (1 - 5 lambda)/"
              f"(cV (1 - 4 lambda)) = {b}{extra}", ok)


# ====================================================================== N3
def n3():
    section("N3. Lemma N5 (fibre jets): Ric(Z1) - Ric(Z2) - [flat difference] along the fibre over x0, in powers of h = y - g(x0)")
    rng = random.Random(31)
    m = 14
    for eta, lam, K in (((1, 1, 1, 1), Fr(1, 7), 3), ((1, -1, -1, -1), Fr(1, 2), 2), ((1, 1, -1, -1), Fr(2, 3), 2)):
        t1 = time.time()
        R = rand_R(rng, eta)
        pb = normal_base(eta, R, None)
        pb = (pb[0], pb[1], rand_cubic(rng))
        cH, cV = Fr(2), Fr(3)
        S1 = fibre_series(pb, lam, cH, cV, K, 'Z1')
        S2 = fibre_series(pb, lam, cH, cV, K, 'Z2')
        A, ctx = S1['A'], S1['ctx']
        c, _, _ = consts(lam, cV)
        G15, G14 = dewitt_at(ctx, S1['y0'], Fr(1, 5)), dewitt_at(ctx, S1['y0'], Fr(1, 4))
        res = {}
        nzA = 0
        for blk, cond in (('HH', lambda a, b: a < n and b < n), ('HV', lambda a, b: (a < n) != (b < n)),
                          ('VV', lambda a, b: a >= n and b >= n)):
            oks = []
            for o in range(K + 1):
                okb = True
                for a in range(m):
                    for b in range(m):
                        if not cond(a, b):
                            continue
                        if a < n and b < n:
                            fl = -(1 / c) * cH * S1['y0'][a][b]
                        elif a >= n and b >= n:
                            fl = Fr(-5, 4) * G15[a - n][b - n] + G14[a - n][b - n]
                        else:
                            fl = ctx.zero
                        if (S1['Ric'][a][b] - S2['Ric'][a][b] - fl).part(o):
                            okb = False
                        if o == 2 and blk != 'HV':
                            nzA += len(S2['Ric'][a][b].part(2))
                oks.append(okb)
            res[blk] = oks
        want = {'HH': [True, True, True] + ([False] if K == 3 else []), 'HV': [True, True] + [False] * (K - 1),
                'VV': [True, True, True] + ([False] if K == 3 else [])}
        check(f"base {sig_of(eta)}, lambda = {lam}, random R and cubic terms (nabla R != 0): agreement order by order"
              f" {res} (orders 0..{K}); the order-2 A-terms are nonzero ({nzA} coefficients)  [{time.time() - t1:.0f} s]",
              res == want and nzA > 0)


# ====================================================================== N4
def traceless_V0(eta):
    """basis of V_0 at eta: off-diagonal E_ij and eta_i E_ii - eta_j E_jj"""
    out = []
    for A_, (i, j) in enumerate(IDX):
        if i != j:
            v = [Fr(0)] * 10
            v[A_] = Fr(1)
            out.append(v)
    diag = [A_ for A_, (i, j) in enumerate(IDX) if i == j]
    for k in range(1, n):
        v = [Fr(0)] * 10
        v[diag[0]] = Fr(eta[0])
        v[diag[k]] = -Fr(eta[k])
        out.append(v)
    return out


def symmat(v):
    M = [[Fr(0)] * n for _ in range(n)]
    for A_, (i, j) in enumerate(IDX):
        M[i][j] += v[A_]
        if i != j:
            M[j][i] += v[A_]
    return M


def n4():
    section("N4. Lemma N6 (data at p in closed form): R(.,.)E, Ric and nabla Ric of Z1 at a point p of s(X)")
    rng = random.Random(41)
    m = 14
    out = []
    for eta, lam in (((1, 1, 1, 1), Fr(1, 7)), ((1, -1, -1, -1), Fr(1, 2)), ((1, 1, -1, -1), Fr(2, 3))):
        cH, cV = Fr(2), Fr(3)
        R = rand_R(rng, eta)
        W = weyl_of(R, eta)
        Ric_g = ric_of(R, eta)
        s = scal_of(R, eta)
        Ric0 = any(Ric_g[i][j] - s / 4 * Fr(eta[i]) * (i == j) != 0 for i in range(n) for j in range(n))
        if eta != (1, -1, -1, -1):
            Wp, Wm = weyl_pm(W, eta)
            gen = any(v != 0 for v in Wp.values()) and any(v != 0 for v in Wm.values())
        else:
            gen = any(v != 0 for v in W.values())
        check(f"base {sig_of(eta)}: the random curvature tensor is generic for Schur (W^+, W^- (or W in (1,3)), traceless Ricci and"
              f" scalar parts all nonzero)", gen and Ric0 and s != 0)
        for cub in (False, True):
            t1 = time.time()
            pb = normal_base(eta, R, None)
            if cub:
                pb = (pb[0], pb[1], rand_cubic(rng))
            D = point_data(pb, lam, cH, cV)
            g0, R4, nRic, base = D['g0'], D['R'], D['nRic'], D['base']
            okRb = all(base['R'][k] == R[k] for k in R)
            c, b, bS = consts(lam, cV)
            mu = lam / (1 - 4 * lam)
            g = gmat(eta)
            E = E_vec(eta)
            # R(e_c, e_d) E = R^a_{b c d} E^b
            g0i = D['g0i']

            def RE(cc, dd):
                return [sum(g0i[a][e] * R4[(e, bb, cc, dd)] * E[bb] for e in range(m) for bb in range(m) if E[bb] != 0)
                        for a in range(m)]
            okE = True
            for X in range(n):
                for Y in range(n):
                    okE &= all(v == 0 for v in RE(X, Y))
                for A_ in range(10):
                    Um = symmat([Fr(1) if k == A_ else 0 for k in range(10)])
                    pred = [-Fr(1, 4) * sum(Fr(eta[i]) * Um[i][X] * (i == a) for i in range(n)) if a < n else Fr(0) for a in range(m)]
                    okE &= RE(X, n + A_) == pred
            for A_ in range(10):
                for B in range(10):
                    okE &= all(v == 0 for v in RE(n + A_, n + B))
            # R(E, X) E = X/4
            for X in range(n):
                v = [sum(E[cc] * x for cc, x in enumerate(col)) for col in zip(*[RE(cc, X) for cc in range(m)])]
                okE &= v == [Fr(1, 4) * (a == X) for a in range(m)]
            # nabla Ric predictions
            okH = True   # (nabla_X Ric)(Y,Z) = nabla Ric_g
            okHV = True  # (nabla_X Ric)(Y,U) = -Ric_VV(II(X,Y),U) - Ric_HH(Y, U~X/2)
            okXVV = True
            okU = True   # (nabla_U Ric)_VV = 0 and (nabla_U Ric)_HH = -sym(Ric_g o U~), U in V_0
            okEd = True  # (nabla_E Ric): HH = -Ric_g, HV = 0, VV = 0
            G15 = [[Fr(0)] * 10 for _ in range(10)]
            Et = sp.diag(*eta)
            Ems = []
            for (i, j) in IDX:
                Em = sp.zeros(n, n)
                Em[i, j] = 1
                Em[j, i] = 1
                Ems.append(Em)
            for A_ in range(10):
                for B in range(10):
                    a_, b_ = Et.inv() * Ems[A_], Et.inv() * Ems[B]
                    v = (a_ * b_).trace() - sp.Rational(1, 5) * a_.trace() * b_.trace()
                    G15[A_][B] = Fr(int(v.p), int(v.q))

            def tovec(M):
                return [M[i][j] for (i, j) in IDX]
            RicHH = [[Ric_g[i][j] - cH / c * g[i][j] for j in range(n)] for i in range(n)]
            for X in range(n):
                for Y in range(n):
                    for Z in range(n):
                        okH &= nRic[X][Y][Z] == base['nablaRic'][X][Y][Z]
                    IIm = [[-(cH / (2 * cV)) * ((g[i][X] * g[j][Y] + g[i][Y] * g[j][X]) / 2 + mu * g[X][Y] * g[i][j])
                            for j in range(n)] for i in range(n)]
                    IIv = tovec(IIm)
                    for A_ in range(10):
                        Um = symmat([Fr(1) if k == A_ else 0 for k in range(10)])
                        UX = [sum(Fr(eta[i]) * Um[i][X] * (i == a) for i in range(n)) for a in range(n)]   # U~X
                        pred = -Fr(-5, 4) * sum(IIv[P_] * G15[P_][A_] for P_ in range(10)) - sum(RicHH[Y][a] * UX[a] / 2 for a in range(n))
                        okHV &= nRic[X][Y][n + A_] == pred
                for A_ in range(10):
                    for B in range(10):
                        okXVV &= nRic[X][n + A_][n + B] == 0
            for u in traceless_V0(eta):
                Um = symmat(u)
                Ut = [[Fr(eta[i]) * Um[i][j] for j in range(n)] for i in range(n)]
                for X in range(n):
                    for Z in range(n):
                        val = sum(u[A_] * nRic[n + A_][X][Z] for A_ in range(10))
                        pred = -(sum(Ric_g[X][a] * Ut[a][Z] for a in range(n)) + sum(Ric_g[Z][a] * Ut[a][X] for a in range(n))) / 2
                        okU &= val == pred
                for B in range(10):
                    for C in range(10):
                        okU &= sum(u[A_] * nRic[n + A_][n + B][n + C] for A_ in range(10)) == 0
            for X in range(n):
                for Z in range(n):
                    okEd &= sum(E[cc] * nRic[cc][X][Z] for cc in range(m)) == -Ric_g[X][Z]
                for B in range(10):
                    okEd &= sum(E[cc] * nRic[cc][X][n + B] for cc in range(m)) == 0
            for B in range(10):
                for C in range(10):
                    okEd &= sum(E[cc] * nRic[cc][n + B][n + C] for cc in range(m)) == 0
            ok = okRb and okE and okH and okHV and okXVV and okU and okEd
            out.append(ok)
            check(f"base {sig_of(eta)}, lambda = {lam}, {'with' if cub else 'without'} cubic terms: base curvature = R ({okRb});"
                  f" R(X,Y)E = 0, R(X,U)E = -(1/4) U~X, R(U,W)E = 0, R(E,X)E = X/4 ({okE}); (nabla_X Ric)_HH = nabla Ric_g ({okH});"
                  f" (nabla_X Ric)(Y,U) = -Ric(II(X,Y),U) - Ric(Y, U~X/2) ({okHV}); (nabla_X Ric)_VV = 0 ({okXVV});"
                  f" U in V_0: (nabla_U Ric)_VV = 0, (nabla_U Ric)(X,Z) = -sym Ric_g(X, U~Z) ({okU}); nabla_E Ric = -Ric_g on HH,"
                  f" 0 on HV, VV ({okEd})  [{time.time() - t1:.0f} s]", ok)


# ====================================================================== N5
def n5():
    section("N5. Theorem N7: the degree-1 integrability condition (*) of the line equation at a point p of s(X)")
    cH, cV = Fr(2), Fr(3)
    eta = (1, 1, 1, 1)
    flat = (eta, {}, {})
    res = {}
    for lam in (Fr(13, 56), Fr(1, 5), Fr(0), Fr(1, 7), Fr(1, 2), Fr(2, 3), Fr(1)):
        rows, rhs = star_system(point_data(flat, lam, cH, cV))
        res[str(lam)] = rank_test(rows, rhs)
    check(f"flat base: rank (L, [L | Omega]) = {res}; solvable exactly at lambda = 13/56 and 1/5",
          all((v[0] == v[1]) == (k in ('13/56', '1/5')) for k, v in res.items()) and all(v[0] == 28 for v in res.values()))
    rng = random.Random(51)
    for eta in ((1, 1, 1, 1), (1, -1, -1, -1), (1, 1, -1, -1)):
        t1 = time.time()
        W = weyl_of(rand_R(rng, eta), eta)
        pbW = normal_base(eta, W)
        out = {}
        oksol = True
        for lam in (Fr(13, 56), Fr(1, 5), Fr(1, 7), Fr(1, 2)):
            D = point_data(pbW, lam, cH, cV)
            rows, rhs = star_system(D)
            out[str(lam)] = rank_test(rows, rhs)
            if lam in (Fr(13, 56), Fr(1, 5)):
                x_, ok_, rk_ = FJ.solve_fit([[r[k] for r in rows] for k in range(28)], rhs, FC)
                c, b, bS = consts(lam, cV)
                E = E_vec(eta)
                oksol &= ok_ and rk_ == 28 and x_[:14] == [2 * b * e for e in E] and x_[14:] == [-MU1 * bS * 2 * b * e for e in E]
        check(f"Ricci-flat base {sig_of(eta)} (random Weyl tensor, nabla R(x0) = 0): rank (L, [L|Omega]) = {out}; at 13/56 and"
              f" 1/5 the unique solution is v' = 2b E, n' = -mu_1 beta_S v' ({oksol})  [{time.time() - t1:.0f} s]",
              oksol and out['13/56'][0] == out['13/56'][1] and out['1/5'][0] == out['1/5'][1]
              and out['1/7'][1] > out['1/7'][0] and out['1/2'][1] > out['1/2'][0])
    # special Einstein bases: sigma^2 = 10(1-5l)(13-56l)/(7 cV^2 (1-4l)^2)
    l_, cv_ = sp.symbols('lambda c_V')
    c_s = 4 * cv_ * (1 - 4 * l_)
    b_s = (1 - 5 * l_) / (cv_ * (1 - 4 * l_))
    m1 = sp.Rational(-7, 18)
    sig2 = sp.factor(4 * b_s * (4 * m1 * b_s * c_s - 1) / (m1 * c_s))
    check(f"hand solution of (*) for an Einstein base (v = v_E E, v_0 = 0, n_E = -mu_1 beta_S v_E, v'_E = 2b - sigma/2):"
          f" sigma^2 = 4b(4 mu_1 b c - 1)/(mu_1 c) = {sig2}",
          sp.simplify(sig2 - 10 * (1 - 5 * l_) * (13 - 56 * l_) / (7 * cv_ ** 2 * (1 - 4 * l_) ** 2)) == 0)
    eta = (1, 1, 1, 1)
    for lam in (Fr(1, 7), Fr(13, 30)):
        v = 70 * (1 - 5 * lam) * (13 - 56 * lam)
        sq = Fr(math.isqrt(v.numerator), math.isqrt(v.denominator))
        sigma = sq / (7 * cV * (1 - 4 * lam))
        assert sigma ** 2 == 10 * (1 - 5 * lam) * (13 - 56 * lam) / (7 * cV ** 2 * (1 - 4 * lam) ** 2)
        rk = {}
        oksol = True
        for tag, sg, scale in (('sigma', sigma, 1), ('-sigma', -sigma, 1), ('2 sigma', 2 * sigma, 1), ('sigma at (x0, 2g)', sigma, 2)):
            K = sg * cH / 12
            pb = normal_base(eta, radd((K, R1(eta))), None, scale)
            D = point_data(pb, lam, cH, cV)
            rows, rhs = star_system(D)
            rk[tag] = rank_test(rows, rhs)
            if tag in ('sigma', '-sigma'):
                x_, ok_, _ = FJ.solve_fit([[r[k] for r in rows] for k in range(28)], rhs, FC)
                c, b, bS = consts(lam, cV)
                oksol &= ok_ and x_[:14] == [(2 * b - sg / 2) * e for e in E_vec(eta)]
        check(f"constant curvature bases, lambda = {lam} ({'(14,0) side' if lam < Fr(1, 4) else '(7,7) side'}), cH = 2, cV = 3,"
              f" sigma = scal/cH = +-{abs(sigma)}: ranks {rk}; solvable at +-sigma with v'_E = 2b - sigma/2 ({oksol}), not at 2 sigma, and"
              f" not at the point (x0, 2 g(x0)) of Sigma, where sigma halves",
              oksol and rk['sigma'][0] == rk['sigma'][1] and rk['-sigma'][0] == rk['-sigma'][1]
              and rk['2 sigma'][1] > rk['2 sigma'][0] and rk['sigma at (x0, 2g)'][1] > rk['sigma at (x0, 2g)'][0])
    # other bases fail
    rng = random.Random(52)
    eta = (1, 1, 1, 1)
    Rr = rand_R(rng, eta)
    W = weyl_of(Rr, eta)
    g = gmat(eta)
    P0 = [[Fr(rng.randint(-2, 2)) if i <= j else 0 for j in range(n)] for i in range(n)]
    for i in range(n):
        for j in range(i):
            P0[i][j] = P0[j][i]
    trP = sum(Fr(eta[i]) * P0[i][i] for i in range(n))
    P0 = [[P0[i][j] - trP / 4 * g[i][j] for j in range(n)] for i in range(n)]
    cases = [('random R', normal_base(eta, Rr)), ('Einstein W + R_1', normal_base(eta, radd((1, W), (1, R1(eta))))),
             ('traceless Ricci only', normal_base(eta, kn(P0, g))),
             ('Weyl + random cubic (nabla Ric not Codazzi)', (eta, normal_base(eta, W)[1], rand_cubic(rng)))]
    out = {}
    for nm, pb in cases:
        rows, rhs = star_system(point_data(pb, Fr(13, 56), cH, cV))
        a1 = rank_test(rows, rhs)
        rows, rhs = star_system(point_data(pb, Fr(1, 5), cH, cV))
        a2 = rank_test(rows, rhs)
        out[nm] = (a1, a2)
    check(f"base (4,0), lambda = 13/56 and 1/5: not solvable for {list(out)}: ranks {list(out.values())}",
          all(a[1] > a[0] and b_[1] > b_[0] for a, b_ in out.values()))


# ====================================================================== N6
def n6(n0res):
    section("N6. Proposition N8: lambda = 13/56 is excluded by the trace equation along Sigma = {(x, e^s g(x))}")
    r_, N_, k_ = n0res['r'], n0res['N'], n0res['k']
    Phi = n0res['Phi']
    cv = sp.Symbol('c_V', positive=True)
    lam = sp.Rational(13, 56)
    c = 4 * cv * (1 - 4 * lam)
    m1 = sp.Rational(-7, 18)
    b = (1 - 5 * lam) / (cv * (1 - 4 * lam))
    check(f"lambda = 13/56: c = {sp.simplify(c)}, 4 mu_1 b c = {sp.simplify(4 * m1 * b * c)}; along Sigma (Ricci-flat base) E r = r/2,"
          f" |dr|^2 = r^2/(4c), s_Y = sigma_1/cV", sp.simplify(4 * m1 * b * c - 1) == 0)
    sigma1 = -sp.Rational(23, 2) - (1 + lam) / (1 - 4 * lam)
    sY = sigma1 / cv
    nu_E = r_ * (sp.Rational(1, 4) + m1) / c                     # from E(Er) = c nu - mu_1 r with Er = r/2
    nu_Tr = sp.Rational(7, 468) * r_ * sY - Phi / 13 + sp.Rational(3, 28) * (r_ ** 2 / (4 * c)) / r_
    Dm = sp.expand(sp.cancel(r_ * (nu_Tr - nu_E)))
    cst = sp.simplify(Dm.subs(r_, 0))
    check(f"r (nu_trace - nu_E) is a polynomial in r = r(x0, e^s g(x0)) (which varies, = r0 e^(s/2)) with constant term {cst} != 0:"
          f" no solution on any interval of s", sp.denom(sp.together(Dm)).free_symbols <= {N_, cv} and cst != 0
          and sp.simplify(cst + sp.Rational(7 * 312, 18 * 13) / (224 ** 4 * N_ ** 3)) == 0)


# ====================================================================== N7
def n7():
    section("N7. Propositions N9, N10: lambda = 1/5, order eps^2 on the slice S through p: <q, T'> = alpha |W|^2")
    lam = Fr(1, 5)
    eta = (1, 1, 1, 1)
    alphas = {}
    for (cH, cV), seeds in (((Fr(2), Fr(3)), (71, 72)), ((Fr(1), Fr(1)), (73,)), ((Fr(3), Fr(1, 2)), (74,))):
        for seed in seeds:
            t1 = time.time()
            rng = random.Random(seed)
            W = weyl_of(rand_R(rng, eta), eta)
            S = fibre_series(normal_base(eta, W), lam, cH, cV, 2, 'Z1')
            A, ctx = S['A'], S['ctx']
            GF = dewitt_at(ctx, S['y0'], lam)
            c, _, _ = consts(lam, cV)
            kap = -1 / c
            Bm = [[S['Ric'][n + a][n + b] - kap * cV * GF[a][b] for b in range(10)] for a in range(10)]
            GJ = [[FJ.jconst(cV * GF[a][b], 0) for b in range(10)] for a in range(10)]
            GFi = FJ.jinv(GJ, 0, ctx)
            GFi = [[GFi[a][b][0][()] for b in range(10)] for a in range(10)]
            hatB = [[sum((GFi[a][k] * Bm[k][b] for k in range(10)), ctx.zero) for b in range(10)] for a in range(10)]
            o0 = all(not hatB[a][b].part(0) for a in range(10) for b in range(10))
            o1 = all(not hatB[a][b].part(1) for a in range(10) for b in range(10))
            D2 = {}
            for k in range(10):
                for l_ in range(10):
                    mm = [0] * 10
                    mm[k] += 1
                    mm[l_] += 1
                    fac = 2 if k == l_ else 1
                    D2[(k, l_)] = sp.Matrix(10, 10, lambda a, b: sp.Rational(hatB[a][b].coef(tuple(mm)) * fac))
            Et = sp.diag(*eta)
            Ems = []
            for (i, j) in IDX:
                Em = sp.zeros(n, n)
                Em[i, j] = 1
                Em[j, i] = 1
                Ems.append(Em)
            Gm = sp.Matrix(10, 10, lambda a, b: sp.Rational(cV) * ((Et.inv() * Ems[a] * Et.inv() * Ems[b]).trace()
                                                                    - sp.Rational(lam) * (Et.inv() * Ems[a]).trace() * (Et.inv() * Ems[b]).trace()))
            tr = [(Et.inv() * Ems[a]).trace() for a in range(10)]
            Bs = sp.Matrix.vstack(*[k_.T for k_ in sp.Matrix([tr]).nullspace()])
            gS = Bs * Gm * Bs.T
            gSi = gS.inv()
            Pr = Bs.T * gSi * Bs * Gm

            def tomat(v):
                return sum((Ems[a] * v[a] for a in range(10)), sp.zeros(n, n))

            def tovec(M):
                return sp.Matrix([M[i, j] for (i, j) in IDX])

            def q(X, Y, Z):
                dXY = sum((D2[(k, l_)] * X[k] * Y[l_] for k in range(10) for l_ in range(10)), sp.zeros(10, 10))
                dXZ = sum((D2[(k, l_)] * X[k] * Z[l_] for k in range(10) for l_ in range(10)), sp.zeros(10, 10))
                return Pr * (dXY * Z - dXZ * Y)

            def c0(X, Y, Z):
                Xt, Yt, Zt = [Et.inv() * tomat(v) for v in (X, Y, Z)]
                br = lambda a, b: a * b - b * a  # noqa: E731
                return tovec(Et * (-sp.Rational(1, 4) * br(br(Yt, Zt), Xt)))

            def c1(X, Y, Z):
                return (X.T * Gm * Y)[0] * Z - (X.T * Gm * Z)[0] * Y
            Nr = sp.Matrix(9, 9, lambda i, j: 0)
            rr = random.Random(seed + 100)
            for i in range(9):
                for j in range(i, 9):
                    v = sp.Rational(rr.randint(-3, 3), rr.randint(1, 3))
                    Nr[i, j] = v
                    Nr[j, i] = v
            # make N traceless w.r.t. g_S
            Nr = Nr - (gSi * Nr).trace() / 9 * gS
            basis = [Bs[i, :].T for i in range(9)]
            dual = [sum((gSi[i, j] * basis[j] for j in range(9)), sp.zeros(10, 1)) for i in range(9)]
            coords = lambda X: gSi * (Bs * Gm * X)  # noqa: E731

            def Nw(X, Y, Z):
                return (coords(X).T * Nr * coords(Y))[0] * Z - (coords(X).T * Nr * coords(Z))[0] * Y
            trip = list(itertools.product(range(9), repeat=3))

            def tab(F, vecs):
                return {t_: F(vecs[t_[0]], vecs[t_[1]], vecs[t_[2]]) for t_ in trip}
            Fb = {nm_: tab(F, basis) for nm_, F in (('c0', c0), ('c1', c1), ('q', q), ('N', Nw))}
            Fd = {nm_: tab(F, dual) for nm_, F in (('c0', c0), ('c1', c1))}

            def pair(Fa, Fdual):
                return sp.nsimplify(sum((Fa[t_].T * Gm * Fdual[t_])[0] for t_ in trip))
            p00, p01, p11 = pair(Fb['c0'], Fd['c0']), pair(Fb['c0'], Fd['c1']), pair(Fb['c1'], Fd['c1'])
            beta = p01 / p11
            Tpd = {t_: Fd['c0'][t_] - beta * Fd['c1'][t_] for t_ in trip}
            pq = pair(Fb['q'], Tpd)
            pN = pair(Fb['N'], Tpd)
            p1 = pair(Fb['c1'], Tpd)
            W2 = sp.Rational(norm2(W, eta))
            alpha = sp.nsimplify(pq / W2)
            alphas.setdefault((cH, cV), []).append(alpha)
            check(f"(cH, cV) = ({cH}, {cV}), random Weyl W (seed {seed}): hatA = P(G_F^-1 Ric_VV - kappa)P vanishes to orders 0, 1"
                  f" ({o0}, {o1}); <c0,c0>, <c0,c1>, <c1,c1> = {p00}, {p01}, {p11}; T' = c0 - (<c0,c1>/<c1,c1>) c1 kills c1 ({p1})"
                  f" and N^ for a random g_S-traceless N ({pN}); <q, T'> = {pq} = alpha |W|^2 with alpha = {alpha}"
                  f"  [{time.time() - t1:.0f} s]", o0 and o1 and p1 == 0 and pN == 0 and alpha != 0)
    a23 = alphas[(Fr(2), Fr(3))]
    allv = [a for v in alphas.values() for a in v]
    fit = all(a == sp.Rational(45, 8) / (sp.Rational(ch) ** 2 * sp.Rational(cv)) for (ch, cv), v in alphas.items() for a in v)
    check(f"alpha is the same for the two Weyl tensors at (cH, cV) = (2, 3) (Schur: an O(4)-invariant quadratic form on Weyl tensors"
          f" is a multiple of |W|^2): {[str(x) for x in a23]}; alpha = 45/(8 cH^2 cV) for all three parameter sets ({fit}), never 0",
          a23[0] == a23[1] and all(a != 0 for a in allv) and fit)


# ====================================================================== N8
def n8(n0res):
    section("N8. Proposition N11 (flat model): radial reduction r = r(s) and the trace equation")
    l_, cv, s, C = sp.symbols('lambda c_V s C')
    c = 4 * cv * (1 - 4 * l_)
    b = (1 - 5 * l_) / (cv * (1 - 4 * l_))
    m1 = sp.Rational(-7, 18)
    f = sp.Function('f')
    # Hess f(s) = f'' ds^2 + f' g_H/(2c);  mu1 f Ric + nu g, Ric = a g + b (g_H + c ds^2): g_H coefficient f'/(2c) = mu1 b f,
    # ds^2 coefficient f'' = mu1 b c f  (the a g + nu g part is absorbed by nu = -mu1 a f)
    sol = sp.solve(sp.numer(sp.together(2 * c * m1 * b * 2 * c * m1 * b - m1 * b * c)), l_)
    check(f"radial r = f(s) on the flat model: f' = 2 c mu_1 b f and f'' = mu_1 b c f are compatible iff (2 c mu_1 b)^2 = mu_1 b c,"
          f" i.e. lambda in {sorted(sol)}; at lambda = 1/5 (b = 0) f' = 0, so f is constant; at 13/56, f = C e^(s/2)",
          set(sol) == {sp.Rational(13, 56), sp.Rational(1, 5)} and sp.solve(b, l_) == [sp.Rational(1, 5)])
    r_, N_ = n0res['r'], n0res['N']
    Phi = n0res['Phi']
    lam = sp.Rational(13, 56)
    cc = c.subs(l_, lam)
    sigma1 = -sp.Rational(23, 2) - (1 + lam) / (1 - 4 * lam)
    sY = sigma1 / cv
    # flat model, r = C e^{s/2}: Box r = 5 r/(4c), |dr|^2 = r^2/(4c); trace equation 13 Box r + (7/3) r s_Y + 14 Phi - (39/2)|dr|^2/r
    tr = 13 * 5 * r_ / (4 * cc) + sp.Rational(7, 3) * r_ * sY + 14 * Phi - sp.Rational(39, 2) * (r_ ** 2 / (4 * cc)) / r_
    cst = sp.simplify(sp.expand(sp.cancel(r_ * tr)).subs(r_, 0))
    check(f"lambda = 13/56, r = C e^(s/2): r (trace equation) is a polynomial in r with constant term 14 A = {cst} != 0, so the"
          f" trace equation fails on every interval of s", cst != 0)


# ====================================================================== N9
def n9():
    section("N9. Theorem N12(c): existence for a flat base, lambda = 1/5, (14,0), T_e = 0 (F7 + I7(b))")
    s, k, cv = sp.symbols('s kappa c_V', positive=True)
    # I7(b): 312 s^2 + kappa s + 104 r^2 + 15/(2 cV) = 0 with r = 0
    disc = sp.discriminant(312 * s ** 2 + k * s + sp.Rational(15, 2) / cv, s)
    check(f"I7(b) with r = 0: 312 s^2 + kappa s + 15/(2 cV) = 0 has a real root iff its discriminant {sp.factor(disc)} >= 0,"
          f" i.e. cV kappa^2 >= 9360", sp.simplify(disc - (k ** 2 - 9360 / cv)) == 0)
    # det M on r = 0 (J4(d)): kappa^2 (kappa + 224 s)(kappa + 544 s); roots s = -kappa/224, -kappa/544
    v224 = sp.simplify((312 * s ** 2 + k * s + sp.Rational(15, 2) / cv).subs(s, -k / 224))
    v544 = sp.solve(sp.Eq((312 * s ** 2 + k * s + sp.Rational(15, 2) / cv).subs(s, -k / 544), 0), cv)
    check(f"det M = kappa^2 (kappa + 224 s)(kappa + 544 s) on r = 0: s = -kappa/224 is never a root (value {v224} > 0);"
          f" s = -kappa/544 is a root only for cV = {v544}, i.e. cV kappa^2 = {sp.nsimplify(v544[0] * k ** 2)}, where the other root has"
          f" det M != 0", sp.simplify(v224 - (sp.Rational(88, 50176) * k ** 2 + sp.Rational(15, 2) / cv)) == 0 and len(v544) == 1)


def main():
    n0res = n0()
    n1()
    n2()
    n3()
    n4()
    n5()
    n6(n0res)
    n7()
    n8(n0res)
    n9()
    P(f"\n{NCHK[0]} checks, all True (each asserted).")
    P(f"runtime {time.time() - T0:.0f} s")


if __name__ == '__main__':
    main()
