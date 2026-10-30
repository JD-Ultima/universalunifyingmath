"""Referee cross-check of GU round 10 (formalization/geometric-unity/N-z1-neighbourhood.md), REVIEW8-GU.

Independent construction.  Nothing is imported from N_z1_neighbourhood.py, L_nonradial.py, J_slaved_sector.py,
I_elliptic_reduction.py, G_symbols.py or F_jets.py.  The only repository import is xcheck_round6.py (and through it
xcheck_round5.py), unmodified, the round-5/6 referee construction (own Clifford algebra and Shiab templates), used in X5
only, to obtain c1, c13 and the slaving block M of the (14,0) real form.

Written from the definitions here:
  - truncated multivariate jets with exact Fraction coefficients, and multi-dual numbers a + sum_k a_k eps_k
    (eps_i eps_j = 0) as a coefficient ring, so that a Ricci tensor computed at the shifted point p + eps_k d_k gives the
    exact first derivative d_k Ric at p in every direction k at once;
  - a textbook Ricci routine and a Riemann routine from a metric 2-jet (control: the round 3-sphere and the hyperbolic
    3-space);
  - the Zorro metric Z1 = c_H y(dx,dx) + c_V G_lambda|_y(Dy, Dy), Dy = dy - (Gamma^T y + y Gamma) dx, as 2-jets on
    Y = Met(R^4) at an arbitrary point (x_P, y_P), over polynomial base metrics of degree <= 3 in general coordinates;
  - an own assembly of the degree-1 system (*) of Theorem N7;
  - the line-equation algebra of N3(d), N7(c), N8, N11 (sympy) from the line coefficients of L8 step 5.
Sections
  X0  controls (Ricci and Riemann routines; multi-dual derivative)
  X1  Lemma N1 OFF the section (curved cubic bases in general coordinates, Gamma(x0) != 0; three base signatures)
  X2  Lemma N2 (flat model) at y0 != eta, six lambda, three signatures, c_H of both signs; Hess s = g_H/(2c)
  X3  Lemma N6 and N1(e) at p in s(X) (curved cubic bases in general coordinates, generic R, nabla R != 0):
      (a)-(e) in closed form; N6(f) at eight values of lambda in (4,0) (enough for the degree bound of REVIEW8 R10-1),
      three in (1,3), two further (c_H, c_V), one with c_H < 0, and a second random base
  X4  Theorem N7: the 1274 x 28 system (*), own assembly; flat base at seven lambda; constant-curvature bases at
      lambda = 1/7 ((4,0), sigma = +-10/9) and lambda = 13/30 ((1,3), sigma = +-65/33); 2 sigma; the point (x0, 2 g(x0))
  X5  Theorem N12(c): c1, c13, M from the round-6 referee construction; I7(b) roots, c_V kappa^2 >= 9360, det M at the
      roots; the multiples x' W_A and the line (t, u) || (x, y) (observation for REVIEW8 R10-3)
  X6  line-equation algebra (sympy): (star), (Tr), A; N7(c) from the component equations; N8 and N11 constant terms
  X7  Propositions N9, N10: <q, T'> = alpha |W|^2 at lambda = 1/5 from exact second fibre derivatives (nested multi-dual
      shift along the slice), alpha = 45/(8 c_H^2 c_V) at (2, 3) and (1, 1) and for a Lorentzian base, and its independence
      of nabla R(0) (observation for REVIEW8 R10-2)
Exact arithmetic throughout (Fractions; sympy rationals).  Fixed seeds random.Random(20261111..20261116).  One core.
Every check is asserted.  Run from computation/geometric-unity:
    python3 xcheck_round10.py > xcheck_round10.out 2>&1
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
import sympy as sp  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

T0 = time.time()
RES = []


def say(*a):
    print(*a)
    sys.stdout.flush()


def check(label, ok):
    ok = bool(ok)
    RES.append((label, ok))
    say(f"  [{len(RES):3d}] {label}: {ok}")
    assert ok, label


def section(s):
    say(f"\n=== {s}   [{time.time() - T0:.0f} s]")


def tick(t1):
    return f"[{time.time() - t1:.0f} s]"


# ---------------------------------------------------------------- multi-dual numbers  a0 + sum_k a_k eps_k, eps_i eps_j = 0
class MD:
    __slots__ = ('a', 'd')

    def __init__(self, a, d=None):
        self.a = a
        self.d = d if d is not None else {}

    def __add__(self, o):
        if isinstance(o, MD):
            d = dict(self.d)
            for k, v in o.d.items():
                w = d.get(k, 0) + v
                if w:
                    d[k] = w
                else:
                    d.pop(k, None)
            return MD(self.a + o.a, d)
        return MD(self.a + o, self.d)
    __radd__ = __add__

    def __neg__(self):
        return MD(-self.a, {k: -v for k, v in self.d.items()})

    def __sub__(self, o):
        return self + (-o)

    def __rsub__(self, o):
        return (-self) + o

    def __mul__(self, o):
        if isinstance(o, MD):
            a, b = self.a, o.a
            d = {}
            if b:
                for k, v in self.d.items():
                    d[k] = v * b
            if a:
                for k, v in o.d.items():
                    w = d.get(k, 0) + a * v
                    if w:
                        d[k] = w
                    else:
                        d.pop(k, None)
            return MD(a * b, d)
        if not o:
            return MD(Fr(0))
        return MD(self.a * o, {k: v * o for k, v in self.d.items()})
    __rmul__ = __mul__

    def inv(self):
        a = self.a
        if isinstance(a, MD):
            ia = a.inv()
        else:
            if a == 0:
                raise ZeroDivisionError
            ia = 1 / Fr(a)
        ia2 = ia * ia
        return MD(ia, {k: -(v * ia2) for k, v in self.d.items()})

    def __truediv__(self, o):
        if isinstance(o, MD):
            return self * o.inv()
        return self * (1 / Fr(o))

    def __rtruediv__(self, o):
        return self.inv() * o

    def __bool__(self):
        return bool(self.a) or bool(self.d)

    def __eq__(self, o):
        if isinstance(o, MD):
            return self.a == o.a and self.d == o.d
        return self.a == o and not self.d

    def __hash__(self):
        return hash(self.a)

    def __repr__(self):
        return f"MD({self.a}, {self.d})"


def val(x):
    return x.a if isinstance(x, MD) else x


def base_val(x):
    while isinstance(x, MD):
        x = x.a
    return x


def dpart(x, k):
    return x.d.get(k, Fr(0)) if isinstance(x, MD) else Fr(0)


# ---------------------------------------------------------------- truncated polynomial jets in nv variables, degree <= D
class Space:
    def __init__(self, nv, D):
        self.nv, self.D = nv, D
        mons = []
        for d in range(D + 1):
            for comb in itertools.combinations_with_replacement(range(nv), d):
                e = [0] * nv
                for c in comb:
                    e[c] += 1
                mons.append(tuple(e))
        self.mons = mons
        self.idx = {m: i for i, m in enumerate(mons)}
        self.deg = [sum(m) for m in mons]
        self.var = [self.idx[tuple(1 if k == v else 0 for k in range(nv))] for v in range(nv)]
        # multiplication table
        self.mt = []
        for i, mi in enumerate(mons):
            row = {}
            for j, mj in enumerate(mons):
                if self.deg[i] + self.deg[j] <= D:
                    row[j] = self.idx[tuple(p + q for p, q in zip(mi, mj))]
            self.mt.append(row)


def jconst(S, c):
    return {0: c} if c else {}


def jvar(S, v, c0=0):
    out = {S.var[v]: Fr(1)}
    if c0:
        out[0] = c0
    return out


def jadd(a, b, s=1):
    out = dict(a)
    for k, v in b.items():
        w = out.get(k, 0) + (v if s == 1 else s * v)
        if w:
            out[k] = w
        else:
            out.pop(k, None)
    return out


def jscale(a, c):
    if not c:
        return {}
    return {k: v * c for k, v in a.items()}


def jmul(S, a, b):
    if not a or not b:
        return {}
    out = {}
    mt = S.mt
    for i, x in a.items():
        row = mt[i]
        for j, y in b.items():
            k = row.get(j)
            if k is None:
                continue
            w = out.get(k)
            p = x * y
            out[k] = p if w is None else w + p
    return {k: v for k, v in out.items() if v}


def jdiff(S, a, v):
    """partial derivative wrt variable v (degree drops by one; result is exact to degree D-1)"""
    out = {}
    for i, x in a.items():
        m = S.mons[i]
        if m[v]:
            mm = list(m)
            mm[v] -= 1
            out[S.idx[tuple(mm)]] = x * m[v]
    return out


def jembed(S_from, S_to, a, varmap):
    out = {}
    for i, x in a.items():
        m = S_from.mons[i]
        mm = [0] * S_to.nv
        for v, e in enumerate(m):
            if e:
                mm[varmap[v]] += e
        if sum(mm) <= S_to.D:
            out[S_to.idx[tuple(mm)]] = x
    return out


# ---------------------------------------------------------------- dense linear algebra over a field-like ring (Fr or MD)
def minv(M):
    n = len(M)
    A = [list(row) + [Fr(1) if i == j else Fr(0) for j in range(n)] for i, row in enumerate(M)]
    for col in range(n):
        piv = next(r for r in range(col, n) if base_val(A[r][col]) != 0)
        A[col], A[piv] = A[piv], A[col]
        p = A[col][col]
        ip = (1 / p) if not isinstance(p, MD) else p.inv()
        A[col] = [x * ip if x else x for x in A[col]]
        for r in range(n):
            if r != col and A[r][col]:
                f = A[r][col]
                A[r] = [x - f * y if y else x for x, y in zip(A[r], A[col])]
    return [row[n:] for row in A]


def mmul(A, B):
    m = len(B[0])
    out = []
    for row in A:
        r = [Fr(0)] * m
        for k, a in enumerate(row):
            if a:
                Bk = B[k]
                for j in range(m):
                    if Bk[j]:
                        r[j] = r[j] + a * Bk[j]
        out.append(r)
    return out


def jmat_inv(S, M):
    """inverse of a matrix of jets with invertible constant part"""
    n = len(M)
    M0 = [[M[i][j].get(0, Fr(0)) for j in range(n)] for i in range(n)]
    M0i = minv(M0)
    # N = M - M0 (no constant terms); M^-1 = sum_k (-M0i N)^k M0i
    Nn = [[{k: v for k, v in M[i][j].items() if k != 0} for j in range(n)] for i in range(n)]
    A = [[jconst(S, -M0i[i][j]) for j in range(n)] for i in range(n)]
    Q = [[jmatsum(S, [jmul(S, A[i][k], Nn[k][j]) for k in range(n)]) for j in range(n)] for i in range(n)]  # -M0i N
    out = [[jconst(S, M0i[i][j]) for j in range(n)] for i in range(n)]
    term = [[jconst(S, M0i[i][j]) for j in range(n)] for i in range(n)]
    for _ in range(S.D):
        term = [[jmatsum(S, [jmul(S, Q[i][k], term[k][j]) for k in range(n)]) for j in range(n)] for i in range(n)]
        out = [[jadd(out[i][j], term[i][j]) for j in range(n)] for i in range(n)]
    return out


def jmatsum(S, lst):
    out = {}
    for a in lst:
        if a:
            out = jadd(out, a)
    return out


# ---------------------------------------------------------------- curvature from a metric 2-jet (generic ring)
def two_jet(S, G):
    """values, first and second derivatives at 0 of a matrix of jets"""
    m = len(G)
    nv = S.nv
    g = [[G[a][b].get(0, Fr(0)) for b in range(m)] for a in range(m)]
    dg = [[[G[a][b].get(S.var[c], Fr(0)) for b in range(m)] for a in range(m)] for c in range(nv)]
    ddg = [[[[None] * m for _ in range(m)] for _ in range(nv)] for _ in range(nv)]
    for c in range(nv):
        for d in range(nv):
            e = [0] * nv
            e[c] += 1
            e[d] += 1
            k = S.idx[tuple(e)]
            f = 2 if c == d else 1
            for a in range(m):
                for b in range(m):
                    v = G[a][b].get(k, Fr(0))
                    ddg[c][d][a][b] = v * f if v else Fr(0)
    return g, dg, ddg


def christoffel(g, dg, gi=None):
    m = len(g)
    if gi is None:
        gi = minv(g)
    Gl = [[[(dg[b][a][c] + dg[c][a][b] - dg[a][b][c]) * Fr(1, 2) for c in range(m)] for b in range(m)] for a in range(m)]
    Gu = [[[sum((gi[a][d] * Gl[d][b][c] for d in range(m) if gi[a][d] and Gl[d][b][c]), Fr(0)) for c in range(m)]
           for b in range(m)] for a in range(m)]
    return gi, Gl, Gu


def ricci(g, dg, ddg, want_gamma=False):
    """Ric_bd = d_a Gamma^a_bd - d_d Gamma^a_ab + Gamma^a_ae Gamma^e_bd - Gamma^a_de Gamma^e_ab  (unit sphere: Ric = (n-1) g)"""
    m = len(g)
    gi, Gl, Gu = christoffel(g, dg)
    # dgi[e][a][b] = - gi dg[e] gi
    dgi = []
    for e in range(m):
        Mx = mmul(gi, dg[e])
        dgi.append([[-x for x in row] for row in mmul(Mx, gi)])
    w = [sum((dgi[a][a][e] for a in range(m) if dgi[a][a][e]), Fr(0)) for e in range(m)]
    tG = [sum((Gu[a][a][e] for a in range(m) if Gu[a][a][e]), Fr(0)) for e in range(m)]
    Ric = [[None] * m for _ in range(m)]
    for b in range(m):
        for d in range(b, m):
            t1 = Fr(0)
            for e in range(m):
                if w[e] and Gl[e][b][d]:
                    t1 = t1 + w[e] * Gl[e][b][d]
            for a in range(m):
                for e in range(m):
                    x = gi[a][e]
                    if x:
                        y = (ddg[a][b][e][d] + ddg[a][d][e][b] - ddg[a][e][b][d]) * Fr(1, 2)
                        if y:
                            t1 = t1 + x * y
            t2 = Fr(0)
            for a in range(m):
                for e in range(m):
                    if dgi[d][a][e] and dg[b][a][e]:
                        t2 = t2 + dgi[d][a][e] * dg[b][a][e]
                    if gi[a][e] and ddg[d][b][a][e]:
                        t2 = t2 + gi[a][e] * ddg[d][b][a][e]
            t2 = t2 * Fr(1, 2)
            t3 = Fr(0)
            for e in range(m):
                if tG[e] and Gu[e][b][d]:
                    t3 = t3 + tG[e] * Gu[e][b][d]
            t4 = Fr(0)
            for a in range(m):
                for e in range(m):
                    if Gu[a][d][e] and Gu[e][b][a]:
                        t4 = t4 + Gu[a][d][e] * Gu[e][b][a]
            Ric[b][d] = t1 - t2 + t3 - t4
            Ric[d][b] = Ric[b][d]
    if want_gamma:
        return Ric, gi, Gu
    return Ric


def riemann(g, dg, ddg):
    """R^a_{bcd} with R(d_c, d_d) d_b = R^a_{bcd} d_a  (Fractions)"""
    m = len(g)
    gi, Gl, Gu = christoffel(g, dg)
    dgi = []
    for e in range(m):
        Mx = mmul(gi, dg[e])
        dgi.append([[-x for x in row] for row in mmul(Mx, gi)])
    # dGu[e][a][b][c] = d_e Gamma^a_bc
    dGl = [[[[(ddg[e][b][a][c] + ddg[e][c][a][b] - ddg[e][a][b][c]) * Fr(1, 2) for c in range(m)] for b in range(m)]
             for a in range(m)] for e in range(m)]
    dGu = [[[[sum((dgi[e][a][d] * Gl[d][b][c] for d in range(m) if dgi[e][a][d] and Gl[d][b][c]), Fr(0)) +
              sum((gi[a][d] * dGl[e][d][b][c] for d in range(m) if gi[a][d] and dGl[e][d][b][c]), Fr(0))
              for c in range(m)] for b in range(m)] for a in range(m)] for e in range(m)]
    R = {}
    for a in range(m):
        for b in range(m):
            for c in range(m):
                for d in range(c + 1, m):
                    v = dGu[c][a][d][b] - dGu[d][a][c][b]
                    for e in range(m):
                        if Gu[a][c][e] and Gu[e][d][b]:
                            v += Gu[a][c][e] * Gu[e][d][b]
                        if Gu[a][d][e] and Gu[e][c][b]:
                            v -= Gu[a][d][e] * Gu[e][c][b]
                    if v:
                        R[(a, b, c, d)] = v
                        R[(a, b, d, c)] = -v
    return R, gi, Gu


# ---------------------------------------------------------------- the Zorro metric Z1 (and Z2) on Y = Met(R^4)
n = 4
IDX = [(i, j) for i in range(n) for j in range(i, n)]          # fibre coordinates y_A = y_ij, i <= j
AOF = {}
for _A, (_i, _j) in enumerate(IDX):
    AOF[(_i, _j)] = _A
    AOF[(_j, _i)] = _A
S4 = Space(4, 3)          # base jets
S14 = Space(14, 2)        # jets on Y


def Emat(A):
    i, j = IDX[A]
    M = [[Fr(0)] * n for _ in range(n)]
    M[i][j] = Fr(1)
    M[j][i] = Fr(1)
    return M


def vec_of(M):
    """symmetric matrix -> vertical vector (components in d/dy_A):  M = sum_A U^A E_A  =>  U^A = M_ij"""
    return [M[i][j] for (i, j) in IDX]


def mat_of(U):
    M = [[Fr(0)] * n for _ in range(n)]
    for A, (i, j) in enumerate(IDX):
        M[i][j] = M[i][j] + U[A]
        if i != j:
            M[j][i] = M[j][i] + U[A]
    return M


def base_jets(base, xP):
    """base = {exponent tuple (4): 4x4 coefficient matrix}; jets of g_ij at x = xP + xi (S4, degree 3)"""
    X = [jvar(S4, mu, xP[mu]) for mu in range(n)]
    cache = {(0, 0, 0, 0): jconst(S4, Fr(1))}

    def mono(e):
        if e in cache:
            return cache[e]
        k = next(i for i in range(n) if e[i])
        e2 = list(e)
        e2[k] -= 1
        r = jmul(S4, mono(tuple(e2)), X[k])
        cache[e] = r
        return r
    G = [[{} for _ in range(n)] for _ in range(n)]
    for e, C in base.items():
        me = mono(e)
        for i in range(n):
            for j in range(n):
                if C[i][j]:
                    G[i][j] = jadd(G[i][j], jscale(me, C[i][j]))
    return G


def base_gamma(G):
    """Gamma^k_{mu i} as S4 jets (exact to degree 2)"""
    Gi = jmat_inv(S4, G)
    dG = [[[jdiff(S4, G[i][j], mu) for j in range(n)] for i in range(n)] for mu in range(n)]
    Gam = [[[None] * n for _ in range(n)] for _ in range(n)]
    for k in range(n):
        for mu in range(n):
            for i in range(mu, n):
                acc = {}
                for l in range(n):
                    t = jadd(jadd(dG[mu][l][i], dG[i][l][mu]), dG[l][mu][i], -1)
                    acc = jadd(acc, jmul(S4, Gi[k][l], t))
                acc = jscale(acc, Fr(1, 2))
                Gam[k][mu][i] = acc
                Gam[k][i][mu] = acc
    return Gam


def zorro_jets(base, xP, yP, lam, cH, cV, variant='Z1'):
    """the metric on Y as S14 jets at the point (xP, yP); coordinates (x^0..3, y_A)"""
    Gb = base_jets(base, xP)
    Gam4 = base_gamma(Gb)
    vm = list(range(4))
    Gam = [[[jembed(S4, S14, Gam4[k][mu][i], vm) for i in range(n)] for mu in range(n)] for k in range(n)]
    gx = [[jembed(S4, S14, Gb[i][j], vm) for j in range(n)] for i in range(n)]
    Y = [[jadd(jconst(S14, yP[i][j]), {S14.var[4 + AOF[(i, j)]]: Fr(1)}) for j in range(n)] for i in range(n)]
    Yi = jmat_inv(S14, Y)
    # traces tr(Yi E_A) and tr(Yi E_A Yi E_B)
    trA = []
    for A, (i, j) in enumerate(IDX):
        trA.append(jscale(Yi[i][j], Fr(2)) if i != j else Yi[i][i])

    def trAB(A, B):
        i, j = IDX[A]
        k, l = IDX[B]
        pa = [(i, j), (j, i)] if i != j else [(i, i)]
        pb = [(k, l), (l, k)] if k != l else [(k, k)]
        acc = {}
        for (p, q) in pa:          # E_A contains e_p e_q^T
            for (r, s) in pb:      # E_B contains e_r e_s^T ;  tr(Yi e_p e_q^T Yi e_r e_s^T) = Yi[s][p] Yi[q][r]
                acc = jadd(acc, jmul(S14, Yi[s][p], Yi[q][r]))
        return acc
    Gd = [[None] * 10 for _ in range(10)]
    for A in range(10):
        for B in range(A, 10):
            v = jadd(trAB(A, B), jmul(S14, trA[A], trA[B]), -lam) if lam else trAB(A, B)
            Gd[A][B] = v
            Gd[B][A] = v
    # Psi_{A mu} = (Gamma_mu^T y + y Gamma_mu)_ij
    Psi = [[None] * n for _ in range(10)]
    for A, (i, j) in enumerate(IDX):
        for mu in range(n):
            acc = {}
            for k in range(n):
                acc = jadd(acc, jmul(S14, Gam[k][mu][i], Y[k][j]))
                acc = jadd(acc, jmul(S14, Y[i][k], Gam[k][mu][j]))
            Psi[A][mu] = acc
    m = 14
    M = [[None] * m for _ in range(m)]
    for A in range(10):
        for B in range(10):
            M[4 + A][4 + B] = jscale(Gd[A][B], cV)
    GP = [[jmatsum(S14, [jmul(S14, Gd[A][B], Psi[B][mu]) for B in range(10)]) for mu in range(n)] for A in range(10)]
    for A in range(10):
        for mu in range(n):
            v = jscale(GP[A][mu], -cV)
            M[mu][4 + A] = v
            M[4 + A][mu] = v
    for mu in range(n):
        for nu in range(mu, n):
            hor = jscale(Y[mu][nu] if variant == 'Z1' else gx[mu][nu], cH)
            acc = jadd(hor, jscale(jmatsum(S14, [jmul(S14, Psi[A][mu], GP[A][nu]) for A in range(10)]), cV))
            M[mu][nu] = acc
            M[nu][mu] = acc
    Psi0 = [[Psi[A][mu].get(0, Fr(0)) for mu in range(n)] for A in range(10)]
    return M, Psi0, Gb, Psi


def frame(Psi0):
    """columns: e_mu = d_mu + Psi_{A mu} d_A, then d_A"""
    m = 14
    F = [[Fr(0)] * m for _ in range(m)]
    for mu in range(n):
        F[mu][mu] = Fr(1)
        for A in range(10):
            F[4 + A][mu] = Psi0[A][mu]
    for A in range(10):
        F[4 + A][4 + A] = Fr(1)
    return F


# ====================================================================================== helpers
MU1 = Fr(-7, 18)


def gmat(eta, s=1):
    return [[Fr(s * eta[i]) if i == j else Fr(0) for j in range(n)] for i in range(n)]


def rsym(rng, s=2, d=3):
    M = [[None] * n for _ in range(n)]
    for i in range(n):
        for j in range(i, n):
            v = Fr(rng.randint(-s, s), rng.randint(1, d))
            M[i][j] = v
            M[j][i] = v
    return M


def rand_base(rng, eta, degs=(1, 2, 3)):
    """g(x) = eta + random linear, quadratic, cubic terms (general coordinates: Gamma(0) != 0 if 1 in degs)"""
    base = {(0, 0, 0, 0): gmat(eta)}
    for d in degs:
        for comb in itertools.combinations_with_replacement(range(n), d):
            e = [0] * n
            for c in comb:
                e[c] += 1
            base[tuple(e)] = rsym(rng)
    return base


def cc_base(eta, K):
    """(1 - (K/2) eta(x, x)) eta: conformally flat, Gamma(0) = 0, constant curvature K at 0, nabla R(0) = 0"""
    base = {(0, 0, 0, 0): gmat(eta)}
    for mu in range(n):
        e = [0] * n
        e[mu] = 2
        base[tuple(e)] = gmat(eta, -K * eta[mu] / 2)
    return base


def rand_point(rng, eta):
    """y0 = P^T eta P with P = 1 + small random (same signature as eta, y0 != eta)"""
    P = [[Fr(int(i == j)) + Fr(rng.randint(-1, 1), rng.randint(2, 4)) for j in range(n)] for i in range(n)]
    E = gmat(eta)
    return mmul([list(r) for r in zip(*P)], mmul(E, P))


def mT(A):
    return [list(r) for r in zip(*A)]


def mtr(A):
    return sum((A[i][i] for i in range(len(A))), Fr(0))


def msub(A, B, s=1):
    return [[a - s * b for a, b in zip(ra, rb)] for ra, rb in zip(A, B)]


def mscal(A, c):
    return [[a * c for a in r] for r in A]


def br(A, B):
    return msub(mmul(A, B), mmul(B, A))


def consts(lam, cV):
    c = 4 * cV * (1 - 4 * lam)
    b = (1 - 5 * lam) / (cV * (1 - 4 * lam))
    bS = Fr(-5, 4) / cV
    return c, b, bS


def dewitt(y, lam):
    """G_lambda|_y(E_A, E_B) (10 x 10)"""
    yi = minv(y)
    out = [[None] * 10 for _ in range(10)]
    for A in range(10):
        a = mmul(yi, Emat(A))
        for B in range(10):
            bb = mmul(yi, Emat(B))
            out[A][B] = mtr(mmul(a, bb)) - lam * mtr(a) * mtr(bb)
    return out


def base_geometry(base):
    """R^a_bcd, Ric, scal, Christoffels and nabla Ric of the base metric at x = 0 (own 4D routines)"""
    Gb = base_jets(base, [Fr(0)] * n)
    g, dg, ddg = two_jet(S4, Gb)
    R, gi, Gu = riemann(g, dg, ddg)
    Ric = ricci(g, dg, ddg)
    ok_contr = all(sum((R.get((a, b, a, d), Fr(0)) for a in range(n)), Fr(0)) == Ric[b][d] for b in range(n) for d in range(n))
    scal = sum(gi[a][b] * Ric[a][b] for a in range(n) for b in range(n))
    GbM = base_jets(base, [MD(Fr(0), {mu: Fr(1)}) for mu in range(n)])
    gM, dgM, ddgM = two_jet(S4, GbM)
    RicM = ricci(gM, dgM, ddgM)
    nRic = [[[dpart(RicM[a][b], k) - sum(Gu[c][k][a] * Ric[c][b] + Gu[c][k][b] * Ric[a][c] for c in range(n))
              for b in range(n)] for a in range(n)] for k in range(n)]
    Rl = {(a, b, c, d): sum(g[a][e] * R.get((e, b, c, d), Fr(0)) for e in range(n))
          for a, b, c, d in itertools.product(range(n), repeat=4)}
    return dict(g=g, gi=gi, Gu=Gu, R=R, Rl=Rl, Ric=Ric, scal=scal, nRic=nRic, ok_contr=ok_contr)


def kn(h, k):
    return {(a, b, c, d): h[a][c] * k[b][d] + h[b][d] * k[a][c] - h[a][d] * k[b][c] - h[b][c] * k[a][d]
            for a, b, c, d in itertools.product(range(n), repeat=4)}


def generic_for_schur(B, eta):
    """W^+, W^- (or W in (1,3)), traceless Ricci and scalar parts of R all nonzero"""
    g = gmat(eta)
    Ric, s = B['Ric'], B['scal']
    Pm = [[(Ric[i][j] - s / 6 * g[i][j]) / 2 for j in range(n)] for i in range(n)]
    KN = kn(Pm, g)
    W = {k: B['Rl'][k] - KN[k] for k in KN}
    okW = all(sum(Fr(eta[a]) * W[(a, b, a, d)] for a in range(n)) == 0 for b in range(n) for d in range(n))
    sgn = {}
    for p in itertools.permutations(range(4)):
        sg = 1
        for i in range(4):
            for j in range(i + 1, 4):
                if p[i] > p[j]:
                    sg = -sg
        sgn[p] = sg
    star = {}
    for a, b, c, d in itertools.product(range(n), repeat=4):
        v = Fr(0)
        for e, f in itertools.product(range(n), repeat=2):
            ep = sgn.get((a, b, e, f), 0)
            if ep:
                v += Fr(ep, 2) * Fr(eta[e] * eta[f]) * W[(e, f, c, d)]
        star[(a, b, c, d)] = v
    Wnz = any(W.values())
    if eta.count(-1) % 2 == 1:      # *^2 = -1: W^+- = (W -+ i*W)/2 are both nonzero iff W != 0
        gen_w = Wnz
    else:
        gen_w = any(W[k] + star[k] for k in W) and any(W[k] - star[k] for k in W)
    tl = any(Ric[i][j] - s / 4 * g[i][j] for i in range(n) for j in range(n))
    return okW and gen_w and tl and s != 0


def point_data(base, yP0, lam, cH, cV, keys=tuple(range(14)), riem=False, variant='Z1'):
    """Z1 at (0, yP0): metric, Christoffels, Ric and (from the multi-dual shift in the directions `keys`) nabla Ric"""
    xP = [MD(Fr(0), {mu: Fr(1)}) if mu in keys else Fr(0) for mu in range(n)]
    yP = [[None] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            A = AOF[(i, j)]
            yP[i][j] = MD(yP0[i][j], {4 + A: Fr(1)}) if (4 + A) in keys else yP0[i][j]
    M, Psi0, Gb, Psi = zorro_jets(base, xP, yP, lam, cH, cV, variant)
    gM, dgM, ddgM = two_jet(S14, M)
    RicM = ricci(gM, dgM, ddgM)
    m = 14
    g = [[val(x) for x in r] for r in gM]
    dg = [[[val(x) for x in r] for r in M2] for M2 in dgM]
    ddg = [[[[val(x) for x in r] for r in M2] for M2 in M3] for M3 in ddgM]
    gi, Gl, Gu = christoffel(g, dg)
    Ric = [[val(x) for x in r] for r in RicM]
    nRic = {}
    for k in keys:
        nRic[k] = [[dpart(RicM[a][b], k) - sum((Gu[c][k][a] * Ric[c][b] + Gu[c][k][b] * Ric[a][c] for c in range(m)
                                                 if Gu[c][k][a] or Gu[c][k][b]), Fr(0))
                    for b in range(m)] for a in range(m)]
    Psi0v = [[val(x) for x in r] for r in Psi0]
    Psi1 = [[[val(Psi[A][mu].get(S14.var[a], Fr(0))) for a in range(m)] for mu in range(n)] for A in range(10)]
    out = dict(g=g, dg=dg, ddg=ddg, gi=gi, Gl=Gl, Gu=Gu, Ric=Ric, nRic=nRic, Psi0=Psi0v, Psi1=Psi1, Gb=Gb)
    if riem:
        R, _, _ = riemann(g, dg, ddg)
        RA = [[[[Fr(0)] * m for _ in range(m)] for _ in range(m)] for _ in range(m)]   # RA[c][d][a][b] = R^a_bcd
        for (a, b, c, d), v in R.items():
            RA[c][d][a][b] = v
        out['R'] = R
        out['RA'] = RA
    return out


def Rvec(RA, Zc, Zd, Zb):
    m = 14
    out = [Fr(0)] * m
    for c in range(m):
        if not Zc[c]:
            continue
        for d in range(m):
            if not Zd[d] or c == d:
                continue
            f = Zc[c] * Zd[d]
            Mx = RA[c][d]
            for a in range(m):
                row = Mx[a]
                acc = Fr(0)
                for b in range(m):
                    if row[b] and Zb[b]:
                        acc += row[b] * Zb[b]
                if acc:
                    out[a] += f * acc
    return out


def to_frame_vec(Fi, Z):
    return [sum((Fi[i][a] * Z[a] for a in range(14) if Fi[i][a] and Z[a]), Fr(0)) for i in range(14)]


def frame_inv(Psi0):
    m = 14
    Fi = [[Fr(int(i == j)) for j in range(m)] for i in range(m)]
    for A in range(10):
        for mu in range(n):
            Fi[4 + A][mu] = -Psi0[A][mu]
    return Fi


def tensor2_frame(F, T):
    return mmul(mT(F), mmul(T, F))


def nabla_frame(F, nRic, keys):
    """(nabla_{F_i} Ric)(F_j, F_k) for the frame directions i whose coordinate support lies in `keys`"""
    m = 14
    lowered = {k: tensor2_frame(F, nRic[k]) for k in keys}
    out = {}
    for i in range(m):
        sup = [a for a in range(m) if F[a][i]]
        if all(a in keys for a in sup):
            acc = [[Fr(0)] * m for _ in range(m)]
            for a in sup:
                acc = [[x + F[a][i] * y for x, y in zip(r1, r2)] for r1, r2 in zip(acc, lowered[a])]
            out[i] = acc
    return out


def V0_basis(g0):
    """traceless symmetric U (tr(g0^-1 U) = 0), as vertical vectors"""
    gi = minv(g0)
    out = []
    for A in range(10):
        i, j = IDX[A]
        if i != j:
            v = [Fr(0)] * 10
            v[A] = Fr(1)
            out.append(v)
    for k in range(1, n):
        Mx = [[Fr(0)] * n for _ in range(n)]
        Mx[0][0] = 1 / gi[0][0]
        Mx[k][k] = -1 / gi[k][k]
        out.append(vec_of(Mx))
    assert all(mtr(mmul(gi, mat_of(u))) == 0 for u in out)
    return out


# ====================================================================================== X0 controls
def x0():
    section("X0  Controls: own Ricci and Riemann routines; multi-dual derivatives")
    S3 = Space(3, 2)
    ok = True
    for K in (1, -1):
        G = [[{} for _ in range(3)] for _ in range(3)]
        for i in range(3):
            j = {0: Fr(4)}
            for v in range(3):
                e = [0] * 3
                e[v] = 2
                j[S3.idx[tuple(e)]] = Fr(-8 * K)
            G[i][i] = j
        g, dg, ddg = two_jet(S3, G)
        Ric = ricci(g, dg, ddg)
        R, _, _ = riemann(g, dg, ddg)
        ok &= all(Ric[a][b] == 2 * K * g[a][b] for a in range(3) for b in range(3))
        ok &= all(sum((R.get((a, b, a, d), Fr(0)) for a in range(3)), Fr(0)) == Ric[b][d] for b in range(3) for d in range(3))
        # sectional curvature K: R(d0, d1) d1 = K (g11 d0 - g01 d1)
        ok &= R.get((0, 1, 0, 1), Fr(0)) == K * g[1][1]
    check("4|dx|^2/(1 + K|x|^2)^2, K = +-1, dim 3: Ric = 2K g, R^a_{bad} = Ric_bd, R(d0,d1)d1 = K g11 d0 (unit sphere has"
          " Ric = (n-1) g, the convention of F)", ok)
    x = MD(Fr(2), {0: Fr(1)})
    f = (x * x * x + 3) / (x * x - 1)          # d/dx at 2: ((3x^2)(x^2-1) - (x^3+3) 2x)/(x^2-1)^2 = (36 - 44)/9
    check("multi-dual arithmetic: d/dx (x^3 + 3)/(x^2 - 1) at x = 2 is -8/9 (exact)", f.a == Fr(11, 3) and f.d == {0: Fr(-8, 9)})


# ====================================================================================== X1 Lemma N1 off the section
def x1():
    section("X1  Lemma N1 at points OFF the section: curved cubic bases in general coordinates (Gamma(x0) != 0), y0 != g(x0)")
    rng = random.Random(20261111)
    m = 14
    cases = (((1, 1, 1, 1), Fr(0), Fr(2), Fr(3)), ((1, -1, -1, -1), Fr(1, 2), Fr(3), Fr(2)),
             ((1, 1, -1, -1), Fr(2, 3), Fr(-3, 2), Fr(5, 2)), ((1, 1, 1, 1), Fr(1, 5), Fr(1), Fr(-2)))
    for eta, lam, cH, cV in cases:
        t1 = time.time()
        base = rand_base(rng, eta)
        y0 = rand_point(rng, eta)
        M, Psi0, Gb, Psi = zorro_jets(base, [Fr(0)] * n, y0, lam, cH, cV)
        g, dg, ddg = two_jet(S14, M)
        Ric, gi, Gu = ricci(g, dg, ddg, want_gamma=True)
        _, Gl, _ = christoffel(g, dg, gi)
        c, b, bS = consts(lam, cV)
        mu_ = lam / (1 - 4 * lam)
        F = frame(Psi0)
        Psi1 = [[[Psi[A][mu].get(S14.var[a], Fr(0)) for a in range(m)] for mu in range(n)] for A in range(10)]
        gamma0 = any(Psi0[A][mu] for A in range(10) for mu in range(n))
        # (a) T = 0: g(nabla_{d_A} d_B, e_mu) = 0
        okT = all(sum(F[cc][mu] * Gl[cc][4 + A][4 + B] for cc in range(m)) == 0
                  for mu in range(n) for A in range(10) for B in range(10))
        # (a) symmetric part of V nabla_{e_mu} e_nu and its horizontal trace
        def nab(mu, nu):
            Z = [Fr(0)] * m
            for cc in range(m):
                acc = Fr(0)
                for a in range(m):
                    if F[a][mu]:
                        de = Psi1[cc - 4][nu][a] if cc >= 4 else Fr(0)
                        acc += F[a][mu] * de
                        for bb in range(m):
                            if F[bb][nu] and Gu[cc][a][bb]:
                                acc += F[a][mu] * F[bb][nu] * Gu[cc][a][bb]
                Z[cc] = acc
            return [Z[4 + A] - sum(Psi0[A][r] * Z[r] for r in range(n)) for A in range(10)]   # vertical frame components
        okS = True
        trH = [Fr(0)] * 10
        ghi = minv(mscal(y0, cH))
        for mu in range(n):
            for nu in range(mu, n):
                sym = [(p + q) / 2 for p, q in zip(nab(mu, nu), nab(nu, mu))]
                Xf = [y0[i][mu] for i in range(n)]
                Yf = [y0[i][nu] for i in range(n)]
                pred = vec_of([[-(cH / (2 * cV)) * ((Xf[i] * Yf[j] + Yf[i] * Xf[j]) / 2 + mu_ * y0[mu][nu] * y0[i][j])
                                for j in range(n)] for i in range(n)])
                okS &= sym == pred
                for A in range(10):
                    trH[A] += ghi[mu][nu] * sym[A] * (1 if mu == nu else 2)
        okTr = trH == [-e / (2 * cV * (1 - 4 * lam)) for e in vec_of(y0)]
        # (b) E^flat = c ds, |E|^2 = c, nabla E = P_H / 2, L_E g = g_H
        E = [Fr(0)] * n + vec_of(y0)
        Ef = [sum(g[a][bb] * E[bb] for bb in range(m)) for a in range(m)]
        gb0 = [[Gb[i][j].get(0, Fr(0)) for j in range(n)] for i in range(n)]
        gb0i = minv(gb0)
        y0i = minv(y0)
        ds = [-Fr(1, 4) * sum(gb0i[i][j] * Gb[j][i].get(S4.var[mu], Fr(0)) for i in range(n) for j in range(n)) for mu in range(n)]
        ds += [Fr(1, 4) * mtr(mmul(y0i, Emat(A))) for A in range(10)]
        okF = all(Ef[a] == c * ds[a] for a in range(m))
        EE = sum(Ef[a] * E[a] for a in range(m))
        PH = [[Fr(0)] * m for _ in range(m)]
        for mu in range(n):
            for a in range(m):
                PH[a][mu] = F[a][mu]
        nE = [[Fr(int(a == bb and a >= 4)) + sum(Gu[a][bb][cc] * E[cc] for cc in range(m)) for bb in range(m)] for a in range(m)]
        okN = all(nE[a][bb] == PH[a][bb] / 2 for a in range(m) for bb in range(m))
        LE = [[sum(E[cc] * dg[cc][a][bb] for cc in range(m)) + (g[a][bb] if a >= 4 else 0) + (g[a][bb] if bb >= 4 else 0)
               for bb in range(m)] for a in range(m)]
        gH = mmul(mT(PH), mmul(g, PH))
        okL = LE == gH and any(LE[a][bb] for a in range(m) for bb in range(m))
        # (d) Ric(E, .) = -(1/c) E^flat
        RicE = [sum(Ric[a][bb] * E[bb] for bb in range(m)) for a in range(m)]
        okR = all(RicE[a] == -Ef[a] / c for a in range(m))
        sig = (eta.count(1), eta.count(-1))
        check(f"base {sig}, lambda = {lam}, (cH, cV) = ({cH}, {cV}), y0 != g(x0), Gamma(x0) != 0 ({gamma0}): T = 0 ({okT});"
              f" sym V nabla_X Y = -(cH/2cV)[X^b.Y^b + mu y(X,Y) y] ({okS}), horizontal trace -E/(2cV(1-4lambda)) ({okTr});"
              f" E^flat = c ds ({okF}); |E|^2 = c = {c} ({EE == c}); nabla E = P_H/2 ({okN}); L_E g = g_H != 0 ({okL});"
              f" Ric(E,.) = -E^flat/c ({okR}) {tick(t1)}",
              gamma0 and okT and okS and okTr and okF and EE == c and okN and okL and okR)


# ====================================================================================== X2 Lemma N2
def x2():
    section("X2  Lemma N2 (flat model) at y0 != eta: Ric = beta_S g + b (g_H + c ds ds), Hess s = g_H/(2c); Einstein iff lambda = 1/5")
    rng = random.Random(20261112)
    m = 14
    for eta, cH, cV in (((1, 1, 1, 1), Fr(2), Fr(3)), ((1, -1, -1, -1), Fr(-2), Fr(3)), ((1, 1, -1, -1), Fr(3), Fr(-1, 2))):
        t1 = time.time()
        base = {(0, 0, 0, 0): gmat(eta)}
        y0 = rand_point(rng, eta)
        y0i = minv(y0)
        ds = [Fr(0)] * n + [Fr(1, 4) * mtr(mmul(y0i, Emat(A))) for A in range(10)]
        dds = [[Fr(0)] * m for _ in range(m)]
        for A in range(10):
            for B in range(10):
                dds[4 + A][4 + B] = -Fr(1, 4) * mtr(mmul(mmul(y0i, Emat(A)), mmul(y0i, Emat(B))))
        res = {}
        okall = True
        for lam in (Fr(0), Fr(1, 7), Fr(1, 5), Fr(13, 56), Fr(1, 2), Fr(1)):
            M, Psi0, Gb, Psi = zorro_jets(base, [Fr(0)] * n, y0, lam, cH, cV)
            g, dg, ddg = two_jet(S14, M)
            Ric, gi, Gu = ricci(g, dg, ddg, want_gamma=True)
            c, b, bS = consts(lam, cV)
            gH = [[cH * y0[p][q] if p < n and q < n else Fr(0) for q in range(m)] for p in range(m)]
            ok = all(Ric[p][q] == bS * g[p][q] + b * (gH[p][q] + c * ds[p] * ds[q]) for p in range(m) for q in range(m))
            ok &= bS + b == -1 / c and b * c == 4 * (1 - 5 * lam)
            hess = [[dds[p][q] - sum(Gu[r][p][q] * ds[r] for r in range(m)) for q in range(m)] for p in range(m)]
            ok &= all(hess[p][q] == gH[p][q] / (2 * c) for p in range(m) for q in range(m))
            sY = sum(gi[p][q] * Ric[p][q] for p in range(m) for q in range(m))
            ein = all(Ric[p][q] == sY / 14 * g[p][q] for p in range(m) for q in range(m))
            ok &= ein == (lam == Fr(1, 5))
            if lam == Fr(1, 5):
                ok &= all(Ric[p][q] == Fr(-5, 4) / cV * g[p][q] and Ric[p][q] - sY / 2 * g[p][q] == Fr(15, 2) / cV * g[p][q]
                          for p in range(m) for q in range(m))
            res[str(lam)] = ok
            okall &= ok
        sig = (eta.count(1), eta.count(-1))
        check(f"flat base {sig}, (cH, cV) = ({cH}, {cV}), y0 != eta: the N2 formula, beta_S + b = -1/c, bc = 4(1-5lambda) and"
              f" Hess s = g_H/(2c) hold, and Ric is proportional to g exactly at lambda = 1/5 (where Ric = -(5/(4cV)) g,"
              f" G = (15/(2cV)) g, the F7 constant): {res} {tick(t1)}", okall)


# ====================================================================================== X3 Lemma N6, N1(e)
def n6_checks(base, eta, lam, cH, cV, B, full=True):
    m = 14
    g0 = gmat(eta)
    keys = tuple(range(14)) if full else tuple(range(4, 14))
    D = point_data(base, g0, lam, cH, cV, keys=keys, riem=full)
    F = frame(D['Psi0'])
    Fi = frame_inv(D['Psi0'])
    gF = tensor2_frame(F, D['g'])
    RicF = tensor2_frame(F, D['Ric'])
    nF = nabla_frame(F, D['nRic'], keys)
    c, b, bS = consts(lam, cV)
    mu_ = lam / (1 - 4 * lam)
    g0i = minv(g0)
    Ricg = B['Ric']
    E = [Fr(0)] * n + vec_of(g0)
    out = {}
    V0 = V0_basis(g0)
    # N6(f)
    okfV, okfH, hv_nz = True, True, False
    for u in V0:
        Ut = mmul(g0i, mat_of(u))
        T = [[sum(u[A] * nF[4 + A][j][k] for A in range(10)) for k in range(m)] for j in range(m)]
        okfV &= all(T[j][k] == 0 for j in range(4, m) for k in range(4, m))
        for X in range(n):
            for Z in range(n):
                pred = -(sum(Ricg[X][a] * Ut[a][Z] for a in range(n)) + sum(Ricg[Z][a] * Ut[a][X] for a in range(n))) / 2
                okfH &= T[X][Z] == pred
        hv_nz |= any(T[j][k] for j in range(n) for k in range(4, m))
    out['f'] = okfV and okfH
    out['hv_nz'] = hv_nz
    if not full:
        return out
    # gF block form
    G = dewitt(g0, lam)
    okg = all(gF[i][j] == (cH * g0[i][j] if i < n and j < n else (cV * G[i - 4][j - 4] if i >= 4 and j >= 4 else 0))
              for i in range(m) for j in range(m))
    # N6(a)
    G15 = dewitt(g0, Fr(1, 5))
    oka = all(RicF[i][j] == Ricg[i][j] - cH / c * g0[i][j] for i in range(n) for j in range(n))
    oka &= all(RicF[i][j] == 0 for i in range(n) for j in range(4, m))
    oka &= all(RicF[4 + A][4 + Bb] == Fr(-5, 4) * G15[A][Bb] for A in range(10) for Bb in range(10))
    RE = [sum(RicF[i][j] * E[j] for j in range(m)) for i in range(m)]
    gE = [sum(gF[i][j] * E[j] for j in range(m)) for i in range(m)]
    oka &= RE == [-x / c for x in gE]
    for u in V0:
        U = [Fr(0)] * n + u
        oka &= all(sum(RicF[i][j] * U[j] for j in range(m)) == bS * sum(gF[i][j] * U[j] for j in range(m)) for i in range(4, m))
    # N6(b)
    okb = True
    for X in range(n):
        for Y in range(n):
            for Z in range(n):
                okb &= nF[X][Y][Z] == B['nRic'][X][Y][Z]
            II = vec_of([[-(cH / (2 * cV)) * ((g0[i][X] * g0[j][Y] + g0[i][Y] * g0[j][X]) / 2 + mu_ * g0[X][Y] * g0[i][j])
                          for j in range(n)] for i in range(n)])
            for A in range(10):
                UX = [mmul(g0i, Emat(A))[r][X] for r in range(n)]
                pred = -sum(II[P_] * RicF[4 + P_][4 + A] for P_ in range(10)) - sum(RicF[Y][r] * UX[r] for r in range(n)) / 2
                okb &= nF[X][Y][4 + A] == pred
        okb &= all(nF[X][i][j] == 0 for i in range(4, m) for j in range(4, m))
    # N6(c)
    nE = [[sum(E[4 + A] * nF[4 + A][i][j] for A in range(10)) for j in range(m)] for i in range(m)]
    okc = all(nE[i][j] == (-Ricg[i][j] if i < n and j < n else 0) for i in range(m) for j in range(m))
    # N6(d): (nabla_zeta Ric)(E, .) = -(1/c) g(nabla_zeta E, .) - Ric(nabla_zeta E, .), nabla E = P_H/2
    okd = True
    for i in range(m):
        lhs = [sum(nF[i][j][k] * E[j] for j in range(m)) for k in range(m)]
        pred = [(-gF[i][k] / (2 * c) - RicF[i][k] / 2) if i < n else Fr(0) for k in range(m)]
        okd &= lhs == pred
    # Riemann: N1(e) and N6(e)
    RA = D['RA']

    def Rf(i, j, Zb):          # R(F_i, F_j) Zb (Zb in frame components) -> frame components
        Zbc = [sum(F[a][k] * Zb[k] for k in range(m)) for a in range(m)]
        return to_frame_vec(Fi, Rvec(RA, [F[a][i] for a in range(m)], [F[a][j] for a in range(m)], Zbc))
    oke = True
    for X in range(n):
        for Y in range(X + 1, n):
            oke &= all(v == 0 for v in Rf(X, Y, E))
        for A in range(10):
            UX = [mmul(g0i, Emat(A))[r][X] for r in range(n)]
            oke &= Rf(X, 4 + A, E) == [-v / 4 for v in UX] + [Fr(0)] * 10
    for A in range(10):
        for Bb in range(A + 1, 10):
            oke &= all(v == 0 for v in Rf(4 + A, 4 + Bb, E))
    for X in range(n):
        ex = [Fr(int(k == X)) for k in range(m)]
        Ec = [sum(F[a][k] * E[k] for k in range(m)) for a in range(m)]
        v = to_frame_vec(Fi, Rvec(RA, Ec, [F[a][X] for a in range(m)], Ec))
        oke &= v == [Fr(1, 4) * e for e in ex]
    okS = True
    for i1, u in enumerate(V0):
        for i2, w in enumerate(V0):
            if i2 <= i1:
                continue
            Ut, Wt = mmul(g0i, mat_of(u)), mmul(g0i, mat_of(w))
            for up in V0[:5]:
                Upt = mmul(g0i, mat_of(up))
                pred = vec_of(mmul(g0, mscal(br(br(Ut, Wt), Upt), Fr(-1, 4))))
                Zc = [Fr(0)] * n + u
                Zd = [Fr(0)] * n + w
                Zb = [Fr(0)] * n + up
                okS &= to_frame_vec(Fi, Rvec(RA, Zc, Zd, Zb)) == [Fr(0)] * n + pred
    out.update(g=okg, a=oka, b=okb, c=okc, d=okd, e1=oke, e6=okS)
    return out


def x3():
    section("X3  Lemma N6 and N1(e) at p = (0, g(0)) in s(X): curved cubic bases in general coordinates (Gamma(0) != 0), generic R,"
            " nabla R != 0")
    rng = random.Random(20261113)
    for eta, lam_main, lam_more, extra in (((1, 1, 1, 1), Fr(1, 7), (Fr(0), Fr(1, 5), Fr(13, 56), Fr(1, 3), Fr(1, 2), Fr(1), Fr(-1)),
                                            ((Fr(-3), Fr(2), Fr(0)),)),
                                           ((1, -1, -1, -1), Fr(1, 2), (Fr(1, 7), Fr(2)), ((Fr(1, 2), Fr(5), Fr(1, 2)),))):
        t1 = time.time()
        sig = (eta.count(1), eta.count(-1))
        base = rand_base(rng, eta)
        B = base_geometry(base)
        gen = generic_for_schur(B, eta)
        nR = any(B['nRic'][k][a][b2] for k in range(n) for a in range(n) for b2 in range(n))
        check(f"base {sig}: own 4D routines (R^a_bad = Ric_bd: {B['ok_contr']}); R(0) is generic for Schur (W^+, W^- (or W),"
              f" traceless Ricci, scalar all nonzero): {gen}; nabla Ric_g(0) != 0: {nR}", B['ok_contr'] and gen and nR)
        cH, cV = Fr(2), Fr(3)
        o = n6_checks(base, eta, lam_main, cH, cV, B, full=True)
        check(f"base {sig}, lambda = {lam_main}, (cH, cV) = (2, 3), Z1 at p (full 1-jet of Ric in all 14 directions):"
              f" metric in the adapted frame = cH g + cV G_lambda ({o['g']}); N6(a) blocks, E-eigenvalue -1/c, V_0-eigenvalue"
              f" beta_S ({o['a']}); N6(b) ({o['b']}); N6(c) ({o['c']}); N6(d) ({o['d']}); N1(e) R(X,Y)E = 0, R(X,U)E = -U~X/4,"
              f" R(U,W)E = 0, R(E,X)E = X/4 ({o['e1']}); N6(e) R(U,W)U' = -[[U~,W~],U~']/4 on V_0 ({o['e6']}); N6(f)"
              f" (nabla_U Ric)_VV = 0 and (nabla_U Ric)(X,Z) = -(Ric_g(X,U~Z) + Ric_g(Z,U~X))/2 ({o['f']}) {tick(t1)}",
              all(o[k] for k in ('g', 'a', 'b', 'c', 'd', 'e1', 'e6', 'f')))
        res = {str(lam_main): o['f']}
        hv = o['hv_nz']
        t1 = time.time()
        for lam in lam_more:
            o2 = n6_checks(base, eta, lam, cH, cV, B, full=False)
            res[str(lam)] = o2['f']
            hv |= o2['hv_nz']
        check(f"base {sig}: N6(f) at (cH, cV) = (2, 3) for lambda in {list(res)}: {res} (observation: the HV block of nabla_U Ric,"
              f" not claimed by N6(f), is nonzero for this base with nabla R != 0: {hv}) {tick(t1)}", all(res.values()))
        t1 = time.time()
        for cH2, cV2, lam2 in extra:
            o3 = n6_checks(base, eta, lam2, cH2, cV2, B, full=False)
            check(f"base {sig}: N6(f) at (cH, cV) = ({cH2}, {cV2}), lambda = {lam2}: {o3['f']} {tick(t1)}", o3['f'])
        t1 = time.time()
        base2 = rand_base(rng, eta)
        B2 = base_geometry(base2)
        o4 = n6_checks(base2, eta, Fr(1, 5), Fr(2), Fr(3), B2, full=False)
        check(f"base {sig}, a second random cubic base (generic: {generic_for_schur(B2, eta)}): N6(f) at lambda = 1/5,"
              f" (cH, cV) = (2, 3): {o4['f']} {tick(t1)}", o4['f'] and generic_for_schur(B2, eta))


# ====================================================================================== X4 Theorem N7: the system (*)
def star_system(D, keys=tuple(range(14))):
    m = 14
    g, gi, Ric, nRic, RA = D['g'], D['gi'], D['Ric'], D['nRic'], D['RA']
    Rh = mmul(gi, Ric)                                         # Rh[d][b] = (Ric^)^d_b
    nRh = {k: mmul(gi, nRic[k]) for k in keys}
    rows, rhs = [], []
    for a in range(m):
        for b in range(a + 1, m):
            Mx = RA[a][b]                                       # Mx[d][c] = R^d_{c a b}
            for d in range(m):
                row = [Fr(0)] * (2 * m)
                for cc in range(m):
                    row[cc] = Mx[d][cc] - MU1 * (g[a][cc] * Rh[d][b] - g[b][cc] * Rh[d][a])
                    row[m + cc] = -(g[a][cc] * (d == b) - g[b][cc] * (d == a))
                rows.append(row)
                rhs.append(nRh[a][d][b] - nRh[b][d][a])
    return rows, rhs


def solve_exact(rows, rhs):
    """rank of A and of [A | b]; the solution if consistent and unique"""
    A = [list(r) + [bb] for r, bb in zip(rows, rhs)]
    ncol = len(rows[0])
    piv = []
    r = 0
    for col in range(ncol + 1):
        p = next((i for i in range(r, len(A)) if A[i][col] != 0), None)
        if p is None:
            continue
        A[r], A[p] = A[p], A[r]
        inv = 1 / A[r][col]
        A[r] = [x * inv for x in A[r]]
        for i in range(len(A)):
            if i != r and A[i][col] != 0:
                f = A[i][col]
                A[i] = [x - f * y for x, y in zip(A[i], A[r])]
        piv.append(col)
        r += 1
    rkA = sum(1 for cpos in piv if cpos < ncol)
    rkAb = len(piv)
    sol = None
    if rkA == rkAb == ncol:
        sol = [A[i][ncol] for i in range(ncol)]
    return rkA, rkAb, sol


def x4():
    section("X4  Theorem N7: own assembly of the system (*) (1274 equations, unknowns (v, n)/(mu_1 r0) in R^28)")
    cH, cV = Fr(2), Fr(3)
    t1 = time.time()
    eta = (1, 1, 1, 1)
    flat = {(0, 0, 0, 0): gmat(eta)}
    res = {}
    oksol = True
    for lam in (Fr(13, 56), Fr(1, 5), Fr(0), Fr(1, 7), Fr(1, 2), Fr(2, 3), Fr(1)):
        D = point_data(flat, gmat(eta), lam, cH, cV, riem=True)
        rows, rhs = star_system(D)
        rkA, rkAb, sol = solve_exact(rows, rhs)
        res[str(lam)] = (len(rows), rkA, rkAb)
        if lam in (Fr(13, 56), Fr(1, 5)):
            c, b, bS = consts(lam, cV)
            E = [Fr(0)] * n + vec_of(gmat(eta))
            oksol &= sol is not None and sol[:14] == [2 * b * e for e in E] and sol[14:] == [-MU1 * bS * 2 * b * e for e in E]
    check(f"flat base (4,0): (rows, rank A, rank [A|b]) = {res}; consistent exactly at lambda = 13/56 and 1/5, with the unique"
          f" solution v' = 2b E, n' = -mu_1 beta_S v' ({oksol}) {tick(t1)}",
          oksol and all(v[1] == 28 for v in res.values()) and all((v[1] == v[2]) == (k in ('13/56', '1/5')) for k, v in res.items()))
    for eta, lam in (((1, 1, 1, 1), Fr(1, 7)), ((1, -1, -1, -1), Fr(13, 30))):
        t1 = time.time()
        c, b, bS = consts(lam, cV)
        sig2 = 10 * (1 - 5 * lam) * (13 - 56 * lam) / (7 * cV ** 2 * (1 - 4 * lam) ** 2)
        num, den = sig2.numerator, sig2.denominator
        sq = Fr(int(round(num ** 0.5)), int(round(den ** 0.5)))
        assert sq * sq == sig2
        res = {}
        oks = True
        for tag, K, scale in (('+sigma', sq * cH / 12, 1), ('-sigma', -sq * cH / 12, 1), ('2 sigma', 2 * sq * cH / 12, 1),
                              ('+sigma at (x0, 2g)', sq * cH / 12, 2)):
            base = cc_base(eta, K)
            B = base_geometry(base)
            sigma = B['scal'] / cH
            okK = all(B['Rl'][(a, b2, c2, d2)] == K * (gmat(eta)[a][c2] * gmat(eta)[b2][d2] - gmat(eta)[a][d2] * gmat(eta)[b2][c2])
                      for a, b2, c2, d2 in itertools.product(range(n), repeat=4))
            yP = gmat(eta, scale)
            D = point_data(base, yP, lam, cH, cV, riem=True)
            rows, rhs = star_system(D)
            rkA, rkAb, sol = solve_exact(rows, rhs)
            res[tag] = (rkA, rkAb)
            if tag in ('+sigma', '-sigma'):
                E = [Fr(0)] * n + vec_of(yP)
                vE = 2 * b - sigma / 2
                oks &= okK and sol is not None and sol[:14] == [vE * e for e in E] and sol[14:] == [-MU1 * bS * vE * e for e in E]
            else:
                oks &= okK and sol is None
        sig = (eta.count(1), eta.count(-1))
        check(f"constant-curvature base {sig}, lambda = {lam}, (cH, cV) = (2, 3), N7(c) value sigma_g = +-{sq}: (rank A, rank"
              f" [A|b]) = {res}; solvable at +-sigma with v' = (2b - sigma/2) E, n' = -mu_1 beta_S v'; not at 2 sigma; not at the"
              f" point (x0, 2 g(x0)) of Sigma: {oks} {tick(t1)}", oks and res['+sigma'][0] == 28)


# ====================================================================================== X5 Theorem N12(c)
def x5():
    section("X5  Theorem N12(c): c1, c13 and M of the (14,0) real form from the round-6 referee construction (xcheck_round6"
            " X0-X1, own Clifford algebra), not from J_slaved_sector.py")
    t1 = time.time()
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        import xcheck_round6 as R6
        C = R6.Ctx(R6.R5.EUC)
        c1, c13 = R6.x0(C)
        d1 = R6.x1(C)
    M = d1['M']
    check(f"xcheck_round6 X0-X1 rerun as a library in (14,0): {len(R6.RES)} internal checks, all True {tick(t1)}",
          len(R6.RES) > 0 and all(ok for _, ok in R6.RES))
    x, y, z, w, kap, t, u = R6.x, R6.y, R6.z, R6.w, R6.kap, R6.t, R6.u
    s, r = sp.symbols('s r', real=True)
    cV = sp.symbols('c_V', positive=True)
    xp, yp, zp, wp = sp.symbols("x' y' z' w'", real=True)
    sub = {x: -sp.I * xp, y: -sp.I * yp, z: -sp.I * zp, w: -sp.I * wp, t: sp.I * s, u: sp.I * r}
    e1 = sp.expand(-sp.I * c1.subs(sub))
    e13 = sp.expand(-sp.I * c13.subs(sub))
    # Einstein flat model at lambda = 1/5: G = (15/(2 cV)) g, odd equation x' G + e1 g = 0, y' G + e13 g = 0
    sh = {xp: 1, yp: 0}
    eq1 = sp.expand(e1.subs(sh) + sp.Rational(15, 2) / cV)
    eq13 = sp.factor(e13.subs(sh))
    ok1 = sp.expand(eq1.subs(r, 0) - (312 * s ** 2 + kap * s + sp.Rational(15, 2) / cV)) == 0 and sp.expand(eq13 - r * (kap + 208 * s)) == 0
    disc = sp.factor(sp.discriminant(eq1.subs(r, 0), s))
    ok2 = sp.simplify(disc - (kap ** 2 - 9360 / cV)) == 0
    br2 = sp.expand(eq1.subs(s, -kap / 208))
    ok3 = sp.expand(br2 - (kap ** 2 / 416 + 104 * r ** 2 + sp.Rational(15, 2) / cV)) == 0
    check(f"(14,0) real form (x = -i x', t = i s, u = i r): e1 = {e1}, e13 = {e13}; Shiab# and W_A (x' = 1, y' = 0) on the"
          f" Einstein flat model: {sp.factor(eq1)} = 0 and {eq13} = 0; on r = 0 this is 312 s^2 + kappa s + 15/(2 cV) = 0"
          f" ({ok1}), discriminant {disc} >= 0 iff cV kappa^2 >= 9360 ({ok2}); branch s = -kappa/208 gives {br2} > 0, no real"
          f" solution ({ok3})", ok1 and ok2 and ok3)
    detM = sp.factor(sp.expand(M.det()))
    ok = True
    vals = {}
    for nm, mem in (('Shiab#', {xp: 1, yp: 0, zp: 1, wp: 0}), ('W_A', {xp: 1, yp: 0, zp: 0, wp: 0})):
        dm = sp.factor(sp.expand(detM.subs(sub).subs(mem).subs(r, 0)))
        ok &= sp.expand(dm - kap ** 2 * (kap + 224 * s) * (kap + 544 * s)) == 0
        vals[nm] = dm
    roots = sp.solve(312 * s ** 2 + kap * s + sp.Rational(15, 2) / cV, s)
    q224 = sp.simplify((312 * s ** 2 + kap * s + sp.Rational(15, 2) / cV).subs(s, -kap / 224))
    cv544 = sp.solve(sp.Eq((312 * s ** 2 + kap * s + sp.Rational(15, 2) / cV).subs(s, -kap / 544), 0), cV)
    ok &= sp.simplify(q224 - (sp.Rational(88, 50176) * kap ** 2 + sp.Rational(15, 2) / cV)) == 0
    ok &= len(cv544) == 1 and sp.simplify(cv544[0] * kap ** 2 - sp.Rational(277440, 29)) == 0
    # at cV kappa^2 = 277440/29 the other root
    kk = sp.symbols('k', positive=True)
    cvv = sp.Rational(277440, 29) / kk ** 2
    rts = sp.solve(312 * s ** 2 + kk * s + sp.Rational(15, 2) / cvv, s)
    other = [rt for rt in rts if sp.simplify(rt + kk / 544) != 0]
    ok &= len(other) == 1 and sp.simplify((kk + 224 * other[0]) * (kk + 544 * other[0])) != 0
    check(f"det M (from the round-6 M, the (14,0) real form) on r = 0: Shiab# {vals['Shiab#']}, W_A {vals['W_A']}; s = -kappa/224"
          f" is never a root ({q224} > 0); s = -kappa/544 is a root only for cV kappa^2 = 277440/29 ({cv544}), and then the other"
          f" root {other} has det M != 0: {ok}", ok)
    # observation (REVIEW8 R10-1): multiples x' W_A, and the line (s, r) = rho (x', y')
    X_ = sp.symbols('X', positive=True)
    eqX = sp.expand(e1.subs({xp: X_, yp: 0}).subs(r, 0) + 15 * X_ / (2 * cV))
    dX = sp.factor(sp.discriminant(eqX, s))
    rho, Nn = sp.symbols('rho N', positive=True)
    lin = {s: rho * xp, r: rho * yp}
    l1 = sp.factor(sp.expand(e1.subs(lin)))
    l13 = sp.factor(sp.expand(e13.subs(lin)))
    okl = sp.expand(l1 - xp * (kap * rho + 312 * (xp ** 2 + yp ** 2) * rho ** 2)) == 0 and \
        sp.expand(l13 - yp * (kap * rho + 312 * (xp ** 2 + yp ** 2) * rho ** 2)) == 0
    dl = sp.factor(sp.discriminant(kap * rho + 312 * Nn * rho ** 2 + sp.Rational(15, 2) / cV, rho))
    check(f"observation (REVIEW8 R10-3): for x' W_A (x' = X > 0) the r = 0 equation is {eqX} = 0 with discriminant {dX}, i.e. a real root iff"
          f" cV kappa^2 >= 9360 X^2; on the line (s, r) = rho (x', y') e1 = {l1}, e13 = {l13} ({okl}), so there a root needs"
          f" discriminant {dl} >= 0, i.e. cV kappa^2 >= 9360 (x'^2 + y'^2)",
          okl and sp.simplify(dX - (kap ** 2 - 9360 * X_ ** 2 / cV)) == 0 and sp.simplify(dl - (kap ** 2 - 9360 * Nn / cV)) == 0)


# ====================================================================================== X6 line-equation algebra
def x6():
    section("X6  Line-equation algebra (sympy) from the line coefficients of L8 step 5; N7(c) from the component equations; N8, N11")
    k, N_, sg, r = sp.symbols('kappa N sigma r')
    Dd = k + 224 * N_ * sg
    Ht, Qt, Qpt = 576 * N_ / Dd, -258048 * N_ ** 2 / Dd ** 2, 78336 * N_ ** 2 / Dd ** 2
    ct = k * sg + 312 * N_ * sg ** 2
    sig_r = sp.solve(sp.Eq(r, -1 / (224 * N_ * Dd)), sg)[0]
    s1, s2 = sp.diff(sig_r, r), sp.diff(sig_r, r, 2)
    sub = {sg: sig_r}
    hcoef = sp.simplify((Ht * s1).subs(sub))
    ddcoef = sp.simplify((Ht * s2 + Qt * s1 ** 2).subs(sub))
    gcoef = sp.simplify((-Ht * s2 + Qpt * s1 ** 2).subs(sub))
    ok1 = sp.simplify(hcoef + sp.Rational(18, 7) / r) == 0 and ddcoef == 0 and sp.simplify(gcoef + sp.Rational(351, 98) / r ** 2) == 0
    # E~ = -G + c~ g + hcoef (Hess r - Box r g) + gcoef |dr|^2 g ; traceless: -Ric_0 + hcoef (Hess r)_0 = 0
    mu1 = sp.simplify(1 / (hcoef * r))
    # trace (dim 14): 6 sY + 14 c~ + hcoef (Box r - 14 Box r) + 14 gcoef |dr|^2 = 0 with Box r = mu1 r sY + 14 nu
    sY, nu, dr2 = sp.symbols('s_Y nu dr2')
    tr = 6 * sY + 14 * ct.subs(sub) - 13 * hcoef * (mu1 * r * sY + 14 * nu) + 14 * gcoef * dr2
    nus = sp.solve(tr, nu)[0]
    Phi = sp.Rational(7, 18) * r * ct.subs(sub)
    ok2 = mu1 == sp.Rational(-7, 18) and sp.simplify(nus - (sp.Rational(7, 468) * r * sY - Phi / 13 + sp.Rational(3, 28) * dr2 / r)) == 0
    rPhi = sp.expand(sp.cancel(r * Phi))
    A = rPhi.subs(r, 0)
    ok3 = sp.simplify(A - sp.Rational(13, 269746176) / N_ ** 3) == 0 and sp.Poly(rPhi, r).degree() == 2
    check(f"with r = -1/(224 N (kappa + 224 N sigma)): Hessian coefficient {hcoef}, dr dr coefficient {ddcoef}, |dr|^2 g coefficient"
          f" {gcoef} ({ok1}); mu_1 = {mu1}, nu = (7/468) r s_Y - Phi/13 + (3/28)|dr|^2/r ({ok2}); r Phi(r) is a quadratic polynomial"
          f" in r with constant term A = {A} ({ok3})", ok1 and ok2 and ok3)
    # N7 component equations (steps 3-5) for an Einstein base Ric_g = (scal/4) g, sigma = scal/cH
    lam, cV, cH, r0, sgm, vE, nE = sp.symbols('lambda c_V c_H r0 sigma v_E n_E')
    c = 4 * cV * (1 - 4 * lam)
    b = (1 - 5 * lam) / (cV * (1 - 4 * lam))
    bS = sp.Rational(-5, 4) / cV
    m1 = sp.Rational(-7, 18)
    ok4 = sp.simplify(bS + b + 1 / c) == 0 and sp.simplify(b * c - 4 * (1 - 5 * lam)) == 0
    # step 3: -(cH/4) vE U(X,Z) = mu1 r0 [ -(cH b/2) U(X,Z) + (1/2) Ric_g(X, U~Z) ],  Ric_g(X,U~Z) = (sigma cH/4) U(X,Z)
    st3 = sp.solve(sp.Eq(-cH / 4 * vE, m1 * r0 * (-cH * b / 2 + sgm * cH / 8)), vE)[0]
    # step 4: -mu1 vE c beta_S - nE c = 0
    st4 = sp.solve(sp.Eq(-m1 * vE * c * bS - nE * c, 0), nE)[0]
    # step 5: R(X,E)v + mu1 vE c Ric^X + nE c X = mu1 r0 [(nabla_X Ric^)E - (nabla_E Ric^)X]
    lhs5 = -vE / 4 + m1 * vE * c * (sgm / 4 - 1 / c) + nE * c
    rhs5 = m1 * r0 * (-(1 / (2 * c)) - (sgm / 4 - 1 / c) / 2 + sgm / 4)
    eq5 = sp.simplify((lhs5 - rhs5).subs(nE, st4).subs(vE, st3))
    sols = sp.solve(sp.numer(sp.together(eq5)), sgm)
    pred = 10 * (1 - 5 * lam) * (13 - 56 * lam) / (7 * cV ** 2 * (1 - 4 * lam) ** 2)
    ok5 = sp.simplify(st3 - m1 * r0 * (2 * b - sgm / 2)) == 0 and len(sols) == 2 and all(sp.simplify(sv ** 2 - pred) == 0 for sv in sols)
    roots = sp.solve(sp.numer(sp.together(b * (4 * m1 * b * c - 1))), lam)
    ok6 = set(roots) == {sp.Rational(1, 5), sp.Rational(13, 56)} and sp.simplify(4 * m1 * b * c + sp.Rational(56, 9) * (1 - 5 * lam)) == 0
    check(f"N7: beta_S + b = -1/c, bc = 4(1 - 5 lambda) ({ok4}); my reading of steps 3-5 gives v_E = {sp.factor(st3)}, n_E ="
          f" -mu_1 beta_S v_E, and sigma_g^2 = {sp.factor(pred)} ({ok5}); 4 mu_1 b c = -(56/9)(1 - 5 lambda), and b(4 mu_1 b c - 1)"
          f" = 0 iff lambda in {sorted(roots)} ({ok6})", ok4 and ok5 and ok6)
    # N8: lambda = 13/56 along Sigma: E r = r/2, nu = r (1/4 + mu1)/c, |dr|^2 = r^2/(4c), s_Y = sigma_1/cV
    l0 = sp.Rational(13, 56)
    cc = c.subs(lam, l0)
    sigma1 = -sp.Rational(23, 2) - (1 + l0) / (1 - 4 * l0)
    nuE = r * (sp.Rational(1, 4) + m1) / cc
    nuT = (sp.Rational(7, 468) * r * sY - Phi / 13 + sp.Rational(3, 28) * dr2 / r).subs({sY: sigma1 / cV, dr2: r ** 2 / (4 * cc)})
    Dm = sp.expand(sp.cancel(r * (nuT - nuE)))
    cst8 = Dm.subs(r, 0)
    ok7 = sp.simplify(sp.simplify(4 * m1 * b * c).subs(lam, l0) - 1) == 0 and sp.simplify(cst8 + A / 13) == 0 and \
        sp.Poly(Dm, r).degree() == 2
    # N11: flat model, lambda = 13/56, r = C e^{s/2}: nu = -mu1 beta_S r, s_Y = 14 beta_S + 5 b, |dr|^2 = r^2/(4c)
    bb = b.subs(lam, l0)
    nuF = -m1 * bS * r
    nuT2 = (sp.Rational(7, 468) * r * sY - Phi / 13 + sp.Rational(3, 28) * dr2 / r).subs({sY: 14 * bS + 5 * bb, dr2: r ** 2 / (4 * cc)})
    Dm2 = sp.expand(sp.cancel(r * (nuT2 - nuF)))
    cst11 = Dm2.subs(r, 0)
    ode = sp.solve(sp.numer(sp.together((2 * c * m1 * b) ** 2 - m1 * b * c)), lam)
    ok8 = sp.simplify(cst11 + A / 13) == 0 and set(ode) == {sp.Rational(1, 5), sp.Rational(13, 56)} and \
        sp.simplify(-182 * cst11 - 14 * A) == 0
    check(f"N8 (Sigma, Ricci-flat base) and N11 (flat model), lambda = 13/56: r (nu_Tr - nu) is a quadratic polynomial in r with"
          f" constant term -A/13 = {cst8} in both cases ({ok7}); N11's printed '14 A' is the same constant in the normalisation"
          f" (7r/18) x (trace equation) = -182 r (nu_Tr - nu); the radial ODE f' = 2c mu_1 b f, f'' = mu_1 b c f is compatible iff"
          f" lambda in {sorted(ode)} ({ok8})", ok7 and ok8)


# ====================================================================================== X7 Propositions N9, N10: alpha
def weyl_from(R, eta):
    g = gmat(eta)
    Ric = [[sum(Fr(eta[a]) * R[(a, b, a, d)] for a in range(n)) for d in range(n)] for b in range(n)]
    s = sum(Fr(eta[a]) * Ric[a][a] for a in range(n))
    Pm = [[(Ric[i][j] - s / 6 * g[i][j]) / 2 for j in range(n)] for i in range(n)]
    KN = kn(Pm, g)
    return {k: R[k] - KN[k] for k in R}


def rand_curv(rng, eta, k=3):
    out = {key: Fr(0) for key in itertools.product(range(n), repeat=4)}
    for _ in range(k):
        h = rsym(rng, 3, 2)
        K = kn(h, h)
        sg = rng.choice([1, -1])
        for key in out:
            out[key] += sg * K[key]
    return out


def weyl_base(W, eta, rng=None):
    """eta - (1/3) W_{mu a nu b} x^a x^b [+ random cubic terms]: Gamma(0) = 0, R(0) = W, nabla R(0) != 0 with cubic terms"""
    base = {(0, 0, 0, 0): gmat(eta)}
    for a in range(n):
        for b in range(a, n):
            e = [0] * n
            e[a] += 1
            e[b] += 1
            base[tuple(e)] = [[-Fr(1, 3) * (W[(mu, a, nu, b)] + (W[(mu, b, nu, a)] if a != b else 0)) for nu in range(n)]
                              for mu in range(n)]
    if rng is not None:
        for comb in itertools.combinations_with_replacement(range(n), 3):
            e = [0] * n
            for cc in comb:
                e[cc] += 1
            base[tuple(e)] = rsym(rng)
    return base


def alpha_pairing(base, eta, cH, cV):
    """<q, T'> at p for lambda = 1/5, from the exact second fibre derivatives of hatA (nested multi-dual shift along V_0)"""
    lam = Fr(1, 5)
    g0 = gmat(eta)
    V0 = V0_basis(g0)
    k9 = range(9)
    yP = [[None] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            A = AOF[(i, j)]
            inner = MD(g0[i][j], {jj: V0[jj][A] for jj in k9 if V0[jj][A]})
            yP[i][j] = MD(inner, {ii: V0[ii][A] for ii in k9 if V0[ii][A]})
    M, Psi0, Gb, Psi = zorro_jets(base, [Fr(0)] * n, yP, lam, cH, cV)
    gM, dgM, ddgM = two_jet(S14, M)
    RicM = ricci(gM, dgM, ddgM)
    c, b, bS = consts(lam, cV)
    kE = -1 / c
    gVV = [[gM[4 + A][4 + B] for B in range(10)] for A in range(10)]
    Bm = [[RicM[4 + A][4 + B] - kE * gVV[A][B] for B in range(10)] for A in range(10)]
    hatB = mmul(minv(gVV), Bm)
    o0 = all(base_val(x) == 0 for r in hatB for x in r)
    o1 = all((not isinstance(x, MD)) or (all(base_val(v) == 0 for v in x.d.values()) and
                                          (not isinstance(x.a, MD) or not any(x.a.d.values()))) for r in hatB for x in r)
    D2 = [[[[(x.d[i].d.get(j, Fr(0)) if isinstance(x, MD) and i in x.d and isinstance(x.d[i], MD) else Fr(0))
             for x in r] for r in hatB] for j in k9] for i in k9]
    G = dewitt(g0, lam)
    gV = [[cV * G[A][B] for B in range(10)] for A in range(10)]
    E = vec_of(g0)
    EE = sum(gV[A][B] * E[A] * E[B] for A in range(10) for B in range(10))
    Ef = [sum(gV[A][B] * E[B] for B in range(10)) for A in range(10)]
    Pr = [[Fr(int(A == B)) - E[A] * Ef[B] / EE for B in range(10)] for A in range(10)]
    GS = [[sum(V0[a][A] * gV[A][B] * V0[bb][B] for A in range(10) for B in range(10)) for bb in k9] for a in k9]
    GSi = minv(GS)
    g0i = minv(g0)
    Xt = [mmul(g0i, mat_of(v)) for v in V0]

    def mv(Mx, v):
        return [sum(Mx[A][B] * v[B] for B in range(10)) for A in range(10)]

    def q(a, bb, cc):
        u1 = mv(D2[a][bb], V0[cc])
        u2 = mv(D2[a][cc], V0[bb])
        return mv(Pr, [x - y for x, y in zip(u1, u2)])

    def c0(a, bb, cc):
        return vec_of(mmul(g0, mscal(br(br(Xt[bb], Xt[cc]), Xt[a]), Fr(-1, 4))))

    def c1(a, bb, cc):
        return [GS[a][bb] * x - GS[a][cc] * y for x, y in zip(V0[cc], V0[bb])]
    rr = random.Random(20261116)
    Nr = rsym9(rr)
    trN = sum(GSi[i][j] * Nr[j][i] for i in k9 for j in k9)
    Nr = [[Nr[i][j] - trN / 9 * GS[i][j] for j in k9] for i in k9]

    def Nw(a, bb, cc):
        return [Nr[a][bb] * x - Nr[a][cc] * y for x, y in zip(V0[cc], V0[bb])]
    trip = list(itertools.product(k9, repeat=3))
    tabs = {nm: {t_: F(*t_) for t_ in trip} for nm, F in (('q', q), ('c0', c0), ('c1', c1), ('N', Nw))}

    def dual(tab):
        out = tab
        for slot in range(3):
            new = {}
            for t_ in trip:
                acc = [Fr(0)] * 10
                for k in k9:
                    if GSi[t_[slot]][k]:
                        tt = list(t_)
                        tt[slot] = k
                        v = out[tuple(tt)]
                        f = GSi[t_[slot]][k]
                        acc = [x + f * y for x, y in zip(acc, v)]
                new[t_] = acc
            out = new
        return out
    dl = {nm: dual(tabs[nm]) for nm in ('c0', 'c1')}

    def pair(e, d):
        return sum(sum(e[t_][A] * Ef_ for A, Ef_ in enumerate(mv(gV, d[t_]))) for t_ in trip)
    p00, p01, p11 = pair(tabs['c0'], dl['c0']), pair(tabs['c0'], dl['c1']), pair(tabs['c1'], dl['c1'])
    beta = p01 / p11
    Td = {t_: [x - beta * y for x, y in zip(dl['c0'][t_], dl['c1'][t_])] for t_ in trip}
    return dict(o0=o0, o1=o1, pq=pair(tabs['q'], Td), pc1=pair(tabs['c1'], Td), pN=pair(tabs['N'], Td), p=(p00, p01, p11))


def rsym9(rng):
    M = [[None] * 9 for _ in range(9)]
    for i in range(9):
        for j in range(i, 9):
            v = Fr(rng.randint(-3, 3), rng.randint(1, 3))
            M[i][j] = v
            M[j][i] = v
    return M


def x7():
    section("X7  Propositions N9, N10: lambda = 1/5, Weyl base, <q, T'> from the exact second fibre derivatives of hatA along the"
            " slice S (nested multi-dual shift), with and without nabla R(0); bases (4,0) and (1,3)")
    rng = random.Random(20261115)
    E4, E13 = (1, 1, 1, 1), (1, -1, -1, -1)
    vals = {}
    weyls = {'W1': weyl_from(rand_curv(rng, E4), E4), 'W2': weyl_from(rand_curv(rng, E4), E4),
             'W3': weyl_from(rand_curv(rng, E13), E13)}
    for tag, eta, (cH, cV), cub in (('W1 + cubic', E4, (Fr(2), Fr(3)), True), ('W1', E4, (Fr(2), Fr(3)), False),
                                    ('W2 + cubic', E4, (Fr(1), Fr(1)), True), ('W3 + cubic', E13, (Fr(2), Fr(3)), True)):
        t1 = time.time()
        W = weyls[tag.split()[0]]
        base = weyl_base(W, eta, random.Random(20261114) if cub else None)
        B = base_geometry(base)
        okB = all(B['Rl'][k] == W[k] for k in W) and all(B['Ric'][i][j] == 0 for i in range(n) for j in range(n))
        nR = any(B['nRic'][k][a][b2] for k in range(n) for a in range(n) for b2 in range(n))
        W2 = sum(Fr(eta[k[0]] * eta[k[1]] * eta[k[2]] * eta[k[3]]) * W[k] ** 2 for k in W)
        o = alpha_pairing(base, eta, cH, cV)
        alpha = o['pq'] / W2
        pred = Fr(45, 8) / (cH ** 2 * cV)
        vals[tag] = o['pq']
        sig = (eta.count(1), eta.count(-1))
        check(f"{tag}, base {sig}, (cH, cV) = ({cH}, {cV}): base R(0) = W, Ric_g(0) = 0 ({okB}); nabla Ric_g(0) != 0: {nR};"
              f" |W|^2 = {W2}; hatA(p) = 0 and"
              f" d hatA(p) = 0 ({o['o0']}, {o['o1']}); <c0,c0>, <c0,c1>, <c1,c1> = {', '.join(str(v) for v in o['p'])}; <c1, T'> = {o['pc1']}, <N^, T'> ="
              f" {o['pN']} for a random traceless N; <q, T'> = {o['pq']} = alpha |W|^2 with alpha = {alpha} (45/(8 cH^2 cV) ="
              f" {pred}) {tick(t1)}", okB and nR == cub and W2 != 0 and o['o0'] and o['o1'] and o['pc1'] == 0 and o['pN'] == 0
              and alpha == pred)
    check(f"<q, T'> does not depend on nabla R(0) (W1 with and without random cubic terms: {vals['W1 + cubic']} = {vals['W1']}),"
          f" as the structure argument for the VV block requires (REVIEW8 R10-2); the identity <q, T'> = alpha |W|^2 also holds"
          f" for a Lorentzian base (W3), as N.3 assumes for indefinite bases", vals['W1 + cubic'] == vals['W1'])


def main():
    x0()
    x1()
    x2()
    x3()
    x4()
    x5()
    x6()
    x7()
    nt = sum(1 for _, ok in RES if ok)
    say(f"\nSUMMARY: {len(RES)} checks, {nt} True, {len(RES) - nt} False")
    say(f"runtime {time.time() - T0:.0f} s")


if __name__ == '__main__':
    main()
