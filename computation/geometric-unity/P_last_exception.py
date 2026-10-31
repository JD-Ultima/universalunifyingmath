"""GU round 12 (P): the last Z2 exception of M10(b), i.e. the exceptional set Z of M9.

Formalization: formalization/geometric-unity/P-z2-last-exception.md.  Builds on M (M3, M4, M6 steps 1-2, M9, M10), L (L2,
L6, L.0) and J5d (L0(a)).  Imports, unmodified: L_nonradial.py (rand_curv, fibre_data, hessian_system, dm_from_cols: M9's
own pipeline) and the referee scripts xcheck_round9.py / xcheck_round8.py (relaxed_system, ranks, c0t, qt, curv_op,
BR_coeff, pq_bases, rand_curv, endos).  Exact arithmetic throughout (Python Fractions, sympy rationals, sympy QQ[k]);
fixed seeds; one core.  Every check is asserted and printed.

Notation.  Fibre F = Met_(p,q)(R^4) at a base point with g = eta, metric G = G_lambda (c_V absorbed), Levi-Civita
Gamma(h, k) = -(1/2)(h y^-1 k + k y^-1 h), Ric_F = -G_1/4.  For K in gl(4) the Killing field xi_K(y) = -(K^T y + y K)
(= K.y of L-F5).  For b in Sym^2(so(eta)) the A-term A_b(U, W) = sum b^(ab) G(U, xi_a) G(W, xi_b); the base curvature
gives b_R = c_A sum_(c,d) eps_c eps_d R_cd (x) R_cd, c_A = c_V^2/(4 c_H^2), and A_x = A_(b_R) (J5d L0(a)).  The relaxed
system Sigma_y(b) at a fibre point y: unknowns V (10) and xi (10),
     G(R(X,Y) V, Z) - xi(X) B(Y,Z) + xi(Y) B(X,Z) = (nabla_X A)(Y,Z) - (nabla_Y A)(X,Z),     B = -G_1/4 + A_b,
the Ricci identity of Hess w = mu(w) B at y with V = grad w/mu, xi = mu' dw/mu (450 rows).  M9's order-eps Taylor
system (295 unknowns) is equivalent: its extra 275 unknowns (2- and 3-jets) are always solvable for, so
ranks (294, 295) <-> (19, 20).  Cas = kg^-1 in the K-basis K_(ab) (K[a][b] = eta_a, K[b][a] = -eta_b), kg = trace form.

Sections (labels = results of the note)
  P0  the construction: ingredients at random fibre points (Gamma metric, Ric_F = -G_1/4 parallel, xi_K Killing, the
      Killing-field formula for nabla A against an exact derivative), a manufactured solution of the Ricci identity,
      lambda-independence; M9's own samples reproduced (19, 20)/(19, 19); the referee's direct construction at the same
      points (294, 295).
  P1  Proposition P1: A_x depends on R only through b_R = -4 c_A (curv. op.)^2 Cas; null curvatures (b_R = 0, R != 0)
      in bases (1,3), (3,1), (2,2); the trace form is definite in base (4,0).
  P2  Lemma P2 (rank structure) is checked inside P0 and P3 (s-kernel, rank 19 at the section, independent of b).
  P3  Lemma P3: first order in the distance from the section: b in Z forces b in R Cas (own construction; and the
      referee's (C1) tensors: the relaxed (C1) holds only on R Cas).
  P4  Lemma P4: the line R Cas: third order at the section forces k = 0; globally, 20-minors at two off-section points
      have gcd k^2; at the other roots three constructions agree (inconsistent).
Run (from computation/geometric-unity): python3 P_last_exception.py > P_last_exception.out 2>&1   (~12 min on one core, 39 checks)
"""
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
import L_nonradial as L  # noqa: E402
import xcheck_round8 as X8  # noqa: E402
import xcheck_round9 as X9  # noqa: E402

T0 = time.time()
NCHK = [0]
SIGS = ((1, 1, 1, 1), (1, -1, -1, -1), (-1, 1, 1, 1), (1, 1, -1, -1))
SYM = [(i, j) for i in range(4) for j in range(i, 4)]
PAIRS = [(a, b) for a in range(4) for b in range(a + 1, 4)]
kk = sp.Symbol('k')


def P(*a):
    print(*a)
    sys.stdout.flush()


def check(label, ok):
    NCHK[0] += 1
    P(f"  [{NCHK[0]:3d}] {label}: {ok}")
    assert ok, label


def section(s):
    P(f"\n=== {s}   [{time.time() - T0:.0f} s]")


def tk(t1):
    return f"[{time.time() - t1:.0f} s]"


def sl(eta):
    return f"({eta.count(1)},{eta.count(-1)})"


def lam_of(eta):
    return Fr(0) if -1 not in eta else Fr(1, 2)        # the side each real form requires (A4, A5), as in L, M


# ---------------------------------------------------------------- 4x4 Fraction matrices
def mz():
    return [[Fr(0)] * 4 for _ in range(4)]


def eye():
    return [[Fr(1) if i == j else Fr(0) for j in range(4)] for i in range(4)]


def diag(eta):
    return [[Fr(eta[i]) if i == j else Fr(0) for j in range(4)] for i in range(4)]


def mmul(A, B):
    return [[sum(A[i][k] * B[k][j] for k in range(4)) for j in range(4)] for i in range(4)]


def madd(A, B, s=1):
    return [[A[i][j] + s * B[i][j] for j in range(4)] for i in range(4)]


def msc(A, s):
    return [[A[i][j] * s for j in range(4)] for i in range(4)]


def mT(A):
    return [[A[j][i] for j in range(4)] for i in range(4)]


def mtr(A):
    return sum(A[i][i] for i in range(4))


def minv(A):
    M = sp.Matrix(4, 4, lambda i, j: sp.Rational(A[i][j].numerator, A[i][j].denominator)).inv()
    return [[Fr(int(sp.numer(M[i, j])), int(sp.denom(M[i, j]))) for j in range(4)] for i in range(4)]


def coords(M):
    return [M[i][j] for (i, j) in SYM]


def Emat(A):
    i, j = SYM[A]
    M = mz()
    M[i][j] = Fr(1)
    M[j][i] = Fr(1)
    return M


EB = [Emat(A) for A in range(10)]


def killing(K, y):
    """xi_K(y) = -(K^T y + y K)"""
    return msc(madd(mmul(mT(K), y), mmul(y, K)), -1)


def comm(A, B):
    return madd(mmul(A, B), mmul(B, A), -1)


def kbasis(eta):
    out = []
    for a, b in PAIRS:
        K = mz()
        K[a][b] = Fr(eta[a])
        K[b][a] = Fr(-eta[b])
        out.append(K)
    return out


def trace_form(Kb):
    return sp.Matrix(6, 6, lambda a, b: sp.Rational(mtr(mmul(Kb[a], Kb[b]))))


def cas_bc(eta, scale=Fr(1)):
    """the Casimir Cas = kg^-1 (K basis) as a coefficient list [(c, a, b)]"""
    Kb = kbasis(eta)
    kgi = trace_form(Kb).inv()
    return [(scale * Fr(int(sp.numer(kgi[a, b])), int(sp.denom(kgi[a, b]))), a, b)
            for a in range(6) for b in range(6) if kgi[a, b] != 0]


def curv_Ks(Rl, eta, cA):
    """K_s = R_cd (c != d), coefficients c_A eps_c eps_d: b_R"""
    Ks, bc = [], []
    for c in range(4):
        for d in range(4):
            if c != d:
                Ks.append([[Fr(eta[a]) * Fr(Rl(a, b, c, d)) for b in range(4)] for a in range(4)])
                bc.append((cA * eta[c] * eta[d], len(Ks) - 1, len(Ks) - 1))
    return Ks, bc


def R_of(Rnest):
    return lambda a, b, c, d: Rnest[a][b][c][d]  # noqa: E731


def R1_fun(eta):
    return lambda a, b, c, d: Fr((eta[a] if a == c else 0) * (eta[b] if b == d else 0)  # noqa: E731
                                 - (eta[a] if a == d else 0) * (eta[b] if b == c else 0))


def nest(Rl):
    return [[[[Fr(Rl(a, b, c, d)) for d in range(4)] for c in range(4)] for b in range(4)] for a in range(4)]


# ---------------------------------------------------------------- the fibre at a point
class FibrePoint:
    def __init__(self, y0, lam):
        self.y0, self.lam = y0, lam
        self.yi = minv(y0)
        self.Eh = [mmul(self.yi, EB[A]) for A in range(10)]       # hat E_A = y^-1 E_A (endomorphism picture)
        self.Gm = [[self.G(EB[a], EB[b]) for b in range(10)] for a in range(10)]
        Gs = sp.Matrix(10, 10, lambda a, b: sp.Rational(self.Gm[a][b].numerator, self.Gm[a][b].denominator)).inv()
        self.Gi = [[Fr(int(sp.numer(Gs[a, b])), int(sp.denom(Gs[a, b]))) for b in range(10)] for a in range(10)]
        self.G14 = [[self.G(EB[a], EB[b], Fr(1, 4)) for b in range(10)] for a in range(10)]
        self.GamB = [[self.Gam(EB[a], EB[b]) for b in range(10)] for a in range(10)]
        # G(R(E_X, E_Y) E_A, E_Z) with y^-1 R(X,Y)Z = -(1/4)[[X^, Y^], Z^] (F1; checked against the Gamma formula in P0);
        # the lambda-term drops since the commutator is traceless
        Eh = self.Eh
        self.GR = {}
        for X, Y in itertools.combinations(range(10), 2):
            C = comm(Eh[X], Eh[Y])
            rows = []
            for A in range(10):
                Rh = msc(comm(C, Eh[A]), Fr(-1, 4))
                rows.append([mtr(mmul(Rh, Eh[Z])) for Z in range(10)])
            self.GR[(X, Y)] = rows

    def G(self, h, k, lam=None):
        lam = self.lam if lam is None else lam
        a, b = mmul(self.yi, h), mmul(self.yi, k)
        return mtr(mmul(a, b)) - lam * mtr(a) * mtr(b)

    def Gam(self, h, k):
        yi = self.yi
        return msc(madd(mmul(mmul(h, yi), k), mmul(mmul(k, yi), h)), Fr(-1, 2))

    def DGam(self, X, h, k):
        yi = self.yi
        return msc(madd(mmul(mmul(mmul(mmul(h, yi), X), yi), k), mmul(mmul(mmul(mmul(k, yi), X), yi), h)), Fr(1, 2))

    def Rm(self, X, Y, Z):
        return madd(madd(self.DGam(X, Y, Z), self.DGam(Y, X, Z), -1),
                    madd(self.Gam(X, self.Gam(Y, Z)), self.Gam(Y, self.Gam(X, Z)), -1))

    def kdata(self, Ks):
        xs = [killing(K, self.y0) for K in Ks]
        Gx = [[self.G(EB[W], xs[s]) for W in range(10)] for s in range(len(Ks))]
        Gn = [[[self.G(madd(killing(K, EB[X]), self.Gam(EB[X], xs[s])), EB[Y]) for Y in range(10)] for X in range(10)]
              for s, K in enumerate(Ks)]
        return Gx, Gn

    def aterm(self, Ks, bc):
        """B = -G_1/4 + A_b and nB[X][Y][Z] = (nabla_X A_b)(Y, Z) (Killing-field formula)"""
        Gx, Gn = self.kdata(Ks)
        A = [[sum(c * Gx[s][a] * Gx[s2][b] for c, s, s2 in bc) for b in range(10)] for a in range(10)]
        B = [[-self.G14[a][b] + A[a][b] for b in range(10)] for a in range(10)]
        nB = [[[sum(c * (Gn[s][X][Y] * Gx[s2][Z] + Gx[s][Y] * Gn[s2][X][Z]) for c, s, s2 in bc) for Z in range(10)]
                for Y in range(10)] for X in range(10)]
        return A, B, nB


def ricci_rows(pt, B, nB):
    """relaxed Ricci-identity rows: unknowns V (0..9), xi (10..19)"""
    rows, rhs = [], []
    for X, Y in itertools.combinations(range(10), 2):
        GR = pt.GR[(X, Y)]
        for Z in range(10):
            row = {A: GR[A][Z] for A in range(10) if GR[A][Z]}
            if B[Y][Z]:
                row[10 + X] = row.get(10 + X, 0) - B[Y][Z]
            if B[X][Z]:
                row[10 + Y] = row.get(10 + Y, 0) + B[X][Z]
            rows.append({k_: v for k_, v in row.items() if v})
            rhs.append(nB[X][Y][Z] - nB[Y][X][Z])
    return rows, rhs


def relaxed(pt, Ks, bc):
    A, B, nB = pt.aterm(Ks, bc)
    rows, rhs = ricci_rows(pt, B, nB)
    return rows, rhs, A


def qq(v):
    return sp.QQ(v.numerator, v.denominator)


def ranks(rows, rhs, n=20):
    dA, dAb = {}, {}
    for r, row in enumerate(rows):
        ent = {j: qq(v) for j, v in row.items()}
        if ent:
            dA[r] = dict(ent)
        if rhs is not None and rhs[r]:
            ent[n] = qq(rhs[r])
        if ent:
            dAb[r] = ent
    ra = DomainMatrix(dA, (len(rows), n), sp.QQ).rank()
    rb = DomainMatrix(dAb, (len(rows), n + 1), sp.QQ).rank()
    return ra, rb


def apply_rows(rows, vec):
    return [sum(v * vec[j] for j, v in row.items()) for row in rows]


def randP(rng, den=3):
    while True:
        Pm = [[Fr(rng.randint(-2, 2), den) + (1 if i == j else 0) for j in range(4)] for i in range(4)]
        if sp.Matrix(Pm).det() != 0:
            return Pm


def rand_sym(rng, lo=-3, hi=3):
    S = mz()
    for i in range(4):
        for j in range(i, 4):
            S[i][j] = S[j][i] = Fr(rng.randint(lo, hi), rng.randint(1, 3))
    return S


# ---------------------------------------------------------------- dual numbers for exact first derivatives
def dmul(A, B):
    return (mmul(A[0], B[0]), madd(mmul(A[0], B[1]), mmul(A[1], B[0])))


def dtr(A):
    return (mtr(A[0]), mtr(A[1]))


def dG(yi, h, k, lam):
    a, b = dmul(yi, h), dmul(yi, k)
    t1 = dtr(dmul(a, b))
    ta, tb = dtr(a), dtr(b)
    return (t1[0] - lam * ta[0] * tb[0], t1[1] - lam * (ta[0] * tb[1] + ta[1] * tb[0]))


# ====================================================================== P0
def p0():
    section("P0. The construction: ingredients, a manufactured solution, lambda-independence, M9's samples reproduced")
    rng = random.Random(20261201)
    for eta in SIGS:
        t1 = time.time()
        lam = lam_of(eta)
        P0m = randP(rng)
        y0 = mmul(mmul(P0m, diag(eta)), mT(P0m))
        pt = FibrePoint(y0, lam)
        X = rand_sym(rng)
        Z0 = [[Fr(0)] * 4 for _ in range(4)]
        yid = (pt.yi, msc(mmul(mmul(pt.yi, X), pt.yi), -1))            # y^-1 along y0 + eps X
        cst = lambda M: (M, Z0)  # noqa: E731
        ok_metric = all(dG(yid, cst(EB[a]), cst(EB[b]), lam)[1] == pt.G(pt.Gam(X, EB[a]), EB[b]) + pt.G(EB[a], pt.Gam(X, EB[b]))
                        for a in range(10) for b in range(10))
        ok_par = all(dG(yid, cst(EB[a]), cst(EB[b]), Fr(1, 4))[1] == pt.G(pt.Gam(X, EB[a]), EB[b], Fr(1, 4))
                     + pt.G(EB[a], pt.Gam(X, EB[b]), Fr(1, 4)) for a in range(10) for b in range(10))
        ok_ric = all(sum(coords(pt.Rm(EB[A], EB[a], EB[b]))[A] for A in range(10)) == -pt.G14[a][b]
                     for a in range(10) for b in range(10))
        ok_ric &= all(pt.GR[(X_, Y_)][A_][Z_] == pt.G(pt.Rm(EB[X_], EB[Y_], EB[A_]), EB[Z_])
                      for X_, Y_ in itertools.combinations(range(10), 2) for A_ in range(10) for Z_ in range(10))
        Kr = [[Fr(rng.randint(-3, 3)) for _ in range(4)] for _ in range(4)]
        xr = killing(Kr, y0)
        nx = [madd(killing(Kr, EB[a]), pt.Gam(EB[a], xr)) for a in range(10)]
        ok_kill = all(pt.G(nx[a], EB[b]) + pt.G(EB[a], nx[b]) == 0 for a in range(10) for b in range(10))
        # nabla A by the Killing formula against the exact derivative along X
        Rr = X8.rand_curv(rng, eta)
        Ks, bc = curv_Ks(R_of(Rr), eta, Fr(1, 4))
        A, B, nB = pt.aterm(Ks, bc)
        xd = [(killing(K, y0), killing(K, X)) for K in Ks]
        DA = [[sum(c * (lambda p, q: p[0] * q[1] + p[1] * q[0])(dG(yid, cst(EB[a]), xd[s], lam), dG(yid, cst(EB[b]), xd[s2], lam))
                   for c, s, s2 in bc) for b in range(10)] for a in range(10)]
        cX = coords(X)

        def Abil(h, k):
            ch, ck = coords(h), coords(k)
            return sum(ch[a] * ck[b] * A[a][b] for a in range(10) for b in range(10) if ch[a] and ck[b])
        nAX = [[sum(cX[C] * nB[C][a][b] for C in range(10)) for b in range(10)] for a in range(10)]
        ok_nA = all(DA[a][b] - Abil(pt.Gam(X, EB[a]), EB[b]) - Abil(EB[a], pt.Gam(X, EB[b])) == nAX[a][b]
                    for a in range(10) for b in range(10))
        # manufactured solution: w = tr(Ay) tr(By) tr(Cy) + tr(M y N y) + tr(D y); B := Hess w, nB := nabla Hess w
        s = sp.symbols('s0:10')
        ys = sp.Matrix(4, 4, lambda i, j: sp.Rational(y0[i][j]) + sum(s[A_] * EB[A_][i][j] for A_ in range(10)))
        Am = [sp.Matrix(4, 4, lambda i, j: rng.randint(-2, 2)) for _ in range(5)]
        wexp = sp.expand((Am[0] * ys).trace() * (Am[1] * ys).trace() * (Am[2] * ys).trace() + (Am[3] * ys * Am[4] * ys).trace()
                         + (Am[0] * Am[1] * ys).trace())
        pw = sp.Poly(wexp, *s).as_dict()                  # derivatives at s = 0 from the monomial coefficients

        def der(*ix):
            e = [0] * 10
            for i_ in ix:
                e[i_] += 1
            c_ = pw.get(tuple(e), 0)
            for m_ in e:
                c_ *= sp.factorial(m_)
            return L.fr(c_)
        d1 = [der(a) for a in range(10)]
        d2 = [[der(a, b) for b in range(10)] for a in range(10)]
        d3 = [[[der(a, b, c) for c in range(10)] for b in range(10)] for a in range(10)]
        dw = lambda M: sum(d1[C] * v for C, v in enumerate(coords(M)))  # noqa: E731
        d2v = lambda a, M: sum(d2[a][C] * v for C, v in enumerate(coords(M)))  # noqa: E731
        H = [[d2[a][b] - dw(pt.GamB[a][b]) for b in range(10)] for a in range(10)]
        Hb = lambda M, b: sum(v * H[C][b] for C, v in enumerate(coords(M)))  # noqa: E731
        DH = [[[d3[c][a][b] - d2v(c, pt.GamB[a][b]) - dw(pt.DGam(EB[c], EB[a], EB[b])) for b in range(10)] for a in range(10)]
              for c in range(10)]
        nH = [[[DH[c][a][b] - Hb(pt.GamB[c][a], b) - Hb(pt.GamB[c][b], a) for b in range(10)] for a in range(10)]
              for c in range(10)]
        rows, rhs = ricci_rows(pt, H, nH)
        Vw = [sum(pt.Gi[a][b] * d1[b] for b in range(10)) for a in range(10)] + [Fr(0)] * 10
        ok_man = all(v == r for v, r in zip(apply_rows(rows, Vw), rhs)) and any(rhs)
        # lambda-independence of the relaxed system
        rowsL, rhsL, _ = relaxed(pt, Ks, bc)
        same = True
        for lam2 in (Fr(-1), Fr(2)):
            pt2 = FibrePoint(y0, lam2)
            r2, h2, _ = relaxed(pt2, Ks, bc)
            same &= r2 == rowsL and h2 == rhsL
        check(f"base {sl(eta)}, random fibre point y0 = P eta P^T, lambda = {lam}: Gamma is G-metric "
              f"(exact derivative along a random X): {ok_metric}; G_1/4 is parallel: {ok_par}; Ric of R(X,Y)Z = "
              f"D_X Gamma(Y,Z) - D_Y Gamma(X,Z) + [Gamma_X, Gamma_Y] Z equals -G_1/4 (F1), and G(R(X,Y)A, Z) equals the closed form -(1/4) tr([[X^,Y^],A^] Z^) used below (hat = y^-1 .): {ok_ric}; xi_K is Killing for a "
              f"random K in gl(4): {ok_kill}; the Killing-field formula for nabla A_(b_R) equals D_X A - A(Gamma .,.) - "
              f"A(., Gamma .) for a random R: {ok_nA}; manufactured solution w (cubic in y) with B := Hess w: (V, xi) = "
              f"(grad w, 0) solves all 450 Ricci-identity rows (nonzero rhs): {ok_man}; the relaxed system is identical for "
              f"lambda = {lam}, -1, 2: {same}  {tk(t1)}", ok_metric and ok_par and ok_ric and ok_kill and ok_nA and ok_man and same)
    # M9's own samples (M_remaining.m8: random.Random(5) for P, L.rand_curv(Random(90/91 + sum eta)), R_1, control)
    rng = random.Random(5)
    for eta in ((1, -1, -1, -1), (-1, 1, 1, 1), (1, 1, 1, 1), (1, 1, -1, -1)):
        t1 = time.time()
        R1 = R1_fun(eta)
        samples = [('random R #1', L.rand_curv(random.Random(90 + sum(eta)))), ('random R #2', L.rand_curv(random.Random(91 + sum(eta)))),
                   ('R = R_1', R1), ('control P = 1 (section)', L.rand_curv(random.Random(74)))]
        res, okk = [], True
        pt = FibrePoint(diag(eta), lam_of(eta))
        for lab, Rl in samples:
            Pm = sp.eye(4) if lab.startswith('control') else sp.Matrix(4, 4, lambda i, j: (1 if i == j else 0) + sp.Rational(rng.randint(-2, 2), 5))
            Pi = Pm.T.inv()
            Ks, bc = [], []
            for c in range(4):
                for d in range(4):
                    if c != d:
                        Bm = sp.Matrix(4, 4, lambda a, b: eta[a] * sp.Rational(Rl(a, b, c, d)))
                        Bp = Pm.T * Bm * Pi
                        Ks.append([[L.fr(Bp[a, b]) for b in range(4)] for a in range(4)])
                        bc.append((Fr(1, 4) * eta[c] * eta[d], len(Ks) - 1, len(Ks) - 1))
            rows, rhs, _ = relaxed(pt, Ks, bc)
            r = ranks(rows, rhs)
            kv = coords(pt.y0) + [Fr(0)] * 10                 # the dilation (grad s direction) with xi = 0
            okk &= all(v == 0 for v in apply_rows(rows, kv))
            okk &= r == ((19, 19) if lab.startswith('control') else (19, 20))
            res.append(r)
        check(f"base {sl(eta)}: M9's own samples (same seeds and pull-back y0 = P eta P^T -> eta as M_remaining section M8): "
              f"pointwise relaxed system (20 unknowns, 450 rows) has (rank, augmented rank) = {res} for "
              f"{[lab for lab, _ in samples]}; M_remaining.out [32]-[35] has (294, 295) x 3 and (294, 294) for the same samples, "
              f"i.e. the same ranks + 275 (the 2- and 3-jets); the dilation field with xi = 0 is in the kernel  {tk(t1)}", okk)
    # the referee's direct construction (xcheck_round9.relaxed_system) at the same (R, y0) without pull-back
    rng = random.Random(20261202)
    for eta in ((1, 1, 1, 1), (-1, 1, 1, 1)):         # xcheck_round9 section Z8 has (1,3) and (2,2)
        t1 = time.time()
        R = X8.rand_curv(rng, eta)
        P1 = randP(rng)
        y0 = mmul(mmul(P1, diag(eta)), mT(P1))
        pt = FibrePoint(y0, lam_of(eta))
        Ks, bc = curv_Ks(R_of(R), eta, Fr(1, 4))
        rows, rhs, _ = relaxed(pt, Ks, bc)
        mine = ranks(rows, rhs)
        rr, rh, okr, nu = X9.relaxed_system(eta, R, y0, Fr(1, 4), lam_of(eta) if -1 in eta else Fr(0))
        ref = X9.ranks(rr, rh, nu)
        check(f"base {sl(eta)}, random R, random y0 = P eta P^T (direct, no pull-back): own pointwise system {mine}, the "
              f"referee's order-eps Taylor system (xcheck_round9.relaxed_system, 295 unknowns) {ref}: both inconsistent, ranks "
              f"differ by 275  {tk(t1)}", mine == (19, 20) and ref == (294, 295) and okr)


# ====================================================================== P1
def ptwave(eta):
    """null curvatures R = sum H_ij alpha_i (x) alpha_j, alpha_i = l ^ e_i, l null, e_i orthogonal to l"""
    if eta == (1, -1, -1, -1):
        l, es = (1, 1, 0, 0), [(0, 0, 1, 0), (0, 0, 0, 1)]
    elif eta == (-1, 1, 1, 1):
        l, es = (1, 1, 0, 0), [(0, 0, 1, 0), (0, 0, 0, 1)]
    else:                                             # (2,2): l = e0 + e2; e1 spacelike, m = e1 + e3 null
        l, es = (1, 0, 1, 0), [(0, 1, 0, 0), (0, 1, 0, 1)]
    low = lambda v: [eta[a] * v[a] for a in range(4)]  # noqa: E731
    al = [[[Fr(low(l)[a] * low(e)[b] - low(l)[b] * low(e)[a]) for b in range(4)] for a in range(4)] for e in es]
    Hm = [[Fr(1), Fr(2)], [Fr(2), Fr(-3)]]
    R = [[[[sum(Hm[i][j] * al[i][a][b] * al[j][c][d] for i in range(2) for j in range(2)) for d in range(4)] for c in range(4)]
          for b in range(4)] for a in range(4)]
    return R


def p1():
    section("P1. Proposition P1: A_x depends on R only through b_R = -4 c_A R^2 Cas; null curvatures with b_R = 0")
    rng = random.Random(20261203)
    for eta in SIGS:
        t1 = time.time()
        R = X8.rand_curv(rng, eta)
        Rm = X9.curv_op(R, eta)
        b = X9.BR_coeff(R, eta)                        # coefficients of sum eps eps R_cd (x) R_cd in the K basis
        Kb = kbasis(eta)
        kg = trace_form(Kb)
        okid = b == -4 * Rm * Rm * kg.inv()
        P0m = randP(rng)
        y0 = mmul(mmul(P0m, diag(eta)), mT(P0m))
        pt = FibrePoint(y0, lam_of(eta))
        Ks, bc = curv_Ks(R_of(R), eta, Fr(1, 4))
        A1, B1, n1 = pt.aterm(Ks, bc)
        bcK = [(Fr(1, 4) * L.fr(b[a, c]), a, c) for a in range(6) for c in range(6) if b[a, c] != 0]
        A2, B2, n2 = pt.aterm(Kb, bcK)
        check(f"base {sl(eta)}, random R: sum_cd eps_c eps_d R_cd (x) R_cd = -4 R^2 kg^-1 (referee's curv_op, BR_coeff; "
              f"R = curvature operator): {okid}; at a random fibre point the A-term and its covariant derivative built from the "
              f"twelve R_cd equal those built from b_R in the K basis: {A1 == A2 and n1 == n2}  {tk(t1)}", okid and A1 == A2 and n1 == n2)
    for eta in ((1, -1, -1, -1), (-1, 1, 1, 1), (1, 1, -1, -1)):
        t1 = time.time()
        R = ptwave(eta)
        okb = X9.bianchi_ok(R) and any(R[a][b][c][d] for a, b, c, d in itertools.product(range(4), repeat=4))
        ric = X9.ricci(R, eta)
        Rm = X9.curv_op(R, eta)
        b = X9.BR_coeff(R, eta)
        okn = Rm * Rm == sp.zeros(6, 6) and b == sp.zeros(6, 6) and Rm != sp.zeros(6, 6)
        Ks, bc = curv_Ks(R_of(R), eta, Fr(1, 4))
        okA, rk = True, []
        for _ in range(2):
            P0m = randP(rng)
            y0 = mmul(mmul(P0m, diag(eta)), mT(P0m))
            pt = FibrePoint(y0, lam_of(eta))
            rows, rhs, A = relaxed(pt, Ks, bc)
            okA &= all(v == 0 for r_ in A for v in r_) and not any(rhs)
            rk.append(ranks(rows, rhs))
        check(f"base {sl(eta)}: R = sum H_ij (l ^ e_i) (x) (l ^ e_j), l null, H = [[1,2],[2,-3]]: Bianchi and R != 0: {okb}; "
              f"Ric != 0: {any(v for r_ in ric for v in r_)}; R^2 = 0 and b_R = 0 (R != 0): {okn}; the A-term vanishes identically "
              f"at two random fibre points, the relaxed system is consistent there {rk} (so R lies in M9's set Z)  {tk(t1)}",
              okb and okn and okA and all(r == (19, 19) for r in rk) and any(v for r_ in ric for v in r_))
    kg = trace_form(kbasis((1, 1, 1, 1)))
    check(f"base (4,0): the trace form on so(4) is {kg.diagonal().tolist()} (negative definite), so a self-adjoint curvature "
          f"operator with R^2 = 0 vanishes: in base (4,0) b_R = 0 iff R = 0", kg.is_diagonal() and all(v < 0 for v in kg.diagonal()))


# ====================================================================== P3
def eps_jets(pt, Kb, X, order):
    """at eta with K'(eps) = (1 + eps X^T) K (1 + eps X^T)^-1 = sum eps^n K_n: G(E_W, xi^n), G(nabla_Y xi^n, E_Z)"""
    E = pt.y0
    XT = mT(X)
    out_x, out_n = [], []
    for K in Kb:
        C = comm(XT, K)
        Kn = [K]
        for n in range(1, order + 1):
            Kn.append(msc(C, Fr((-1) ** (n + 1))))
            C = mmul(C, XT)
        xi = [killing(Kk, E) for Kk in Kn]
        out_x.append([[pt.G(EB[W], xi[n]) for W in range(10)] for n in range(order + 1)])
        out_n.append([[[pt.G(madd(killing(Kn[n], EB[Y]), pt.Gam(EB[Y], xi[n])), EB[Z]) for Z in range(10)] for Y in range(10)]
                      for n in range(order + 1)])
    return out_x, out_n


def rhs_jet(Gx, Gn, bc, m):
    vec = []
    for Y, Z in itertools.combinations(range(10), 2):
        for W in range(10):
            v = Fr(0)
            for c, a, b in bc:
                for n in range(m + 1):
                    v += c * (Gn[a][n][Y][Z] * Gx[b][m - n][W] + Gx[a][n][Z] * Gn[b][m - n][Y][W]
                              - Gn[a][n][Z][Y] * Gx[b][m - n][W] - Gx[a][n][Y] * Gn[b][m - n][Z][W])
            vec.append(v)
    return vec


def A0_rows(pt):
    zero = [[Fr(0)] * 10 for _ in range(10)]
    B = [[-pt.G14[a][b] for b in range(10)] for a in range(10)]
    rows, _ = ricci_rows(pt, B, [zero] * 10)
    return rows


def p3():
    section("P3. Lemma P3: first order in the distance from the section; b in Z forces b in R Cas")
    for eta in SIGS:
        t1 = time.time()
        pt = FibrePoint(diag(eta), lam_of(eta))
        Kb = kbasis(eta)
        A0 = A0_rows(pt)
        rA0 = ranks(A0, None)[0]
        okker = all(v == 0 for v in apply_rows(A0, coords(pt.y0) + [Fr(0)] * 10))
        sbasis = [(a, b) for a in range(6) for b in range(a, 6)]
        Xs = []
        for i in range(4):
            for j in range(4):
                X = mz()
                X[i][j] = Fr(1)
                Xs.append(X)
        nunk = 21 + 20 * len(Xs)
        rows = []
        for cX, X in enumerate(Xs):
            Gx, Gn = eps_jets(pt, Kb, X, 1)
            cols = []
            for a, b in sbasis:
                bc = [(Fr(1), a, b)] + ([(Fr(1), b, a)] if a != b else [])
                cols.append(rhs_jet(Gx, Gn, bc, 1))
            for r_, row in enumerate(A0):
                ent = {21 + 20 * cX + j: v for j, v in row.items()}
                for ip in range(21):
                    if cols[ip][r_]:
                        ent[ip] = -cols[ip][r_]
                rows.append(ent)
        M = DomainMatrix({i: {j: qq(v) for j, v in r_.items()} for i, r_ in enumerate(rows) if r_}, (len(rows), nunk), sp.QQ)
        ns = M.nullspace().to_Matrix()
        bp = ns[:, :21]
        rk = bp.rank()
        kgi = trace_form(Kb).inv()
        casv = sp.Matrix([[kgi[a, b] for (a, b) in sbasis]])
        okc = rk == 1 and sp.Matrix.vstack(bp.rref()[0][:1, :], casv).rank() == 1
        check(f"base {sl(eta)}: at the section A_0 has rank {rA0} (independent of b; kernel = (dilation, 0): {okker}); the "
              f"first-order system b_1(b, X) = A_0 (V_X, xi_X) for all 16 X in gl(4) ({len(rows)} rows, {nunk} unknowns: b in "
              f"Sym^2 so(eta), (V_X, xi_X)) has nullity {ns.shape[0]} = 16 + 1 and b-projection of rank {rk}, spanned by "
              f"Cas = kg^-1: {okc}  {tk(t1)}", rA0 == 19 and okker and ns.shape[0] == 17 and okc)
    # the referee's (C1) tensors: relaxed (C1)  q(b) in span{c0(V X_a, X_b, X_c)} + span{xi_a(X_b) X_c - xi_a(X_c) X_b}
    for eta in SIGS:
        t1 = time.time()
        Xb, Kb = X8.pq_bases(eta)
        C0 = X9.c0t(Xb)
        pairs = [(a, b) for a in range(6) for b in range(a, 6)]
        gens = []
        for a, b in pairs:
            q = X9.qt(Xb, Kb[a], Kb[b]) + (X9.qt(Xb, Kb[b], Kb[a]) if a != b else 0)
            gens.append(q.flatten())
        for e in range(9):
            for a in range(9):
                T = np.zeros((9, 9, 9, 4, 4), dtype=np.int64)
                T[a] = C0[e]
                gens.append(-T.flatten())
        for a in range(9):
            for f in range(9):
                T = np.zeros((9, 9, 9, 4, 4), dtype=np.int64)
                for c in range(9):
                    T[a, f, c] += Xb[c]
                    T[a, c, f] -= Xb[c]
                gens.append(-T.flatten())
        Mt = np.array(gens).T
        Mt = Mt[np.nonzero(np.any(Mt != 0, axis=1))[0]]
        dm = DomainMatrix([[sp.QQ(int(v)) for v in r_] for r_ in Mt.tolist()], Mt.shape, sp.QQ)
        ns = dm.nullspace().to_Matrix()
        bp = ns[:, :21]
        kg = sp.Matrix(6, 6, lambda a, b: int(np.trace(Kb[a] @ Kb[b])))
        cas = kg.inv()
        casv = sp.Matrix([[cas[a, b] for (a, b) in pairs]])
        okc = bp.rank() == 1 and sp.Matrix.vstack(bp.rref()[0][:1, :], casv).rank() == 1
        check(f"base {sl(eta)}, referee's (C1) tensors (xcheck_round9 c0t, qt; xcheck_round8 pq_bases): the RELAXED (C1), "
              f"q(b)(X,Y,Z) = c0(V X, Y, Z) + xi_X(Y) Z - xi_X(Z) Y with V in End(p), xi in Hom(p, p*) free ({Mt.shape[0]} "
              f"equations, 21 + 81 + 81 unknowns), has nullity {ns.shape[0]}, with b in R Cas: {okc}  {tk(t1)}",
              ns.shape[0] == 1 and okc)


# ====================================================================== P4
def solve20(rows, vec):
    M = sp.Matrix(len(rows), 20, lambda i, j: sp.Rational(rows[i].get(j, 0)))
    sol, params = M.gauss_jordan_solve(sp.Matrix([sp.Rational(v) for v in vec]))
    return list(sol.subs({p_: 0 for p_ in params}))


def rank_with(rows, vecs):
    M = sp.Matrix(len(rows), 20, lambda i, j: sp.Rational(rows[i].get(j, 0)))
    for v in vecs:
        M = M.row_join(sp.Matrix([sp.Rational(x_) for x_ in v]))
    return M.rank()


def p4():
    section("P4. Lemma P4: the line b = k Cas; third order at the section forces k = 0; global and cross-construction checks")
    rng = random.Random(20261204)
    for eta in SIGS:
        t1 = time.time()
        pt = FibrePoint(diag(eta), lam_of(eta))
        Kb = kbasis(eta)
        bc = cas_bc(eta)
        A0 = A0_rows(pt)
        X = [[Fr(rng.randint(-3, 3)) for _ in range(4)] for _ in range(4)]
        Gx, Gn = eps_jets(pt, Kb, X, 3)
        b1, b2, b3 = (rhs_jet(Gx, Gn, bc, m) for m in (1, 2, 3))
        z1 = solve20(A0, b1)
        ok1 = [sp.Rational(v) for v in apply_rows(A0, [L.fr(v) for v in z1])] == [sp.Rational(v) for v in b1]
        # A_2: the A-term at order eps^2 enters only the xi-columns
        A2m = [[sum(c * (Gx[a][1][i] * Gx[b][1][j]) for c, a, b in bc) for j in range(10)] for i in range(10)]
        zero = [[Fr(0)] * 10 for _ in range(10)]
        A2rows, _ = ricci_rows(pt, A2m, [zero] * 10)
        A2rows = [{j: v for j, v in r_.items() if j >= 10} for r_ in A2rows]   # drop the curvature part (order 0 only)
        A2z1 = apply_rows(A2rows, [L.fr(v) for v in z1])
        r0, r_b2, r_a2, r_b3 = 19, rank_with(A0, [b2]), rank_with(A0, [A2z1]), rank_with(A0, [b3])
        r_both = rank_with(A0, [b3, A2z1])
        okk = all(v == 0 for v in apply_rows(A2rows, coords(pt.y0) + [Fr(0)] * 10))
        check(f"base {sl(eta)}, b = Cas, a random X (seeded): b_1 in Im A_0 (Cas passes P3): {ok1}; rank [A_0 | b_2] = {r_b2} "
              f"(order 2 passes); rank [A_0 | A_2 zeta_1] = {r_a2}; rank [A_0 | b_3] = {r_b3}; rank [A_0 | b_3 | A_2 zeta_1] = "
              f"{r_both}; A_2 kills ker A_0: {okk}.  For b = k Cas the order-3 equation A_0 zeta_3 = k b_3 - k^2 A_2 zeta_1 then "
              f"forces k b_3 in Im A_0, i.e. k = 0  {tk(t1)}", ok1 and r_b2 == r0 and r_a2 == r0 and r_b3 == 20 and okk)
    # global: two off-section points, 20-minors of [A(k) | b(k)] with b = k Cas
    rng = random.Random(20261205)
    for eta in SIGS:
        t1 = time.time()
        Kb = kbasis(eta)
        polys, fls, roots, pts = [], [], [], []
        for trial in range(2):
            Pm = randP(rng)
            y0 = mmul(mmul(Pm, diag(eta)), mT(Pm))
            pt = FibrePoint(y0, lam_of(eta))
            rows0, _, _ = relaxed(pt, Kb, [(Fr(0), 0, 0)])
            rows1, rhs1, _ = relaxed(pt, Kb, cas_bc(eta))
            Mk = []
            for r0_, r1_, bv in zip(rows0, rows1, rhs1):
                row = [sp.Integer(0)] * 21
                for j in set(r0_) | set(r1_):
                    a0, a1 = sp.Rational(r0_.get(j, 0)), sp.Rational(r1_.get(j, 0))
                    row[j] = a0 + kk * (a1 - a0)
                row[20] = kk * sp.Rational(bv)
                Mk.append(row)
            Mn = sp.Matrix(Mk).subs(kk, sp.Rational(3, 7))
            _, piv = Mn.T.rref()
            sel = list(piv[:20])
            _, pc = Mn.extract(sel, list(range(21))).rref()
            Dm = sp.Matrix([[Mk[i][j] for j in pc] for i in sel])
            det = DomainMatrix.from_Matrix(Dm).convert_to(sp.QQ[kk]).det()
            p_ = sp.Poly(sp.QQ[kk].to_sympy(det), kk)
            polys.append(p_)
            cst, fl = sp.factor_list(p_.as_expr())
            fls.append(' * '.join(f"({sp.expand(f_)})^{m_}" if m_ > 1 else f"({sp.expand(f_)})" for f_, m_ in fl))
            if trial == 0:
                roots = sorted({r_ for r_ in sp.Poly(p_, kk).ground_roots() if r_ != 0})
            pts.append((Pm, y0, len(piv), 20 in pc))
        g = sp.gcd(polys[0], polys[1])
        check(f"base {sl(eta)}, b = k Cas at two seeded off-section points y = P eta P^T: rank of [A | b] at k = 3/7: "
              f"{[p_[2] for p_ in pts]}; a 20-minor of [A(k) | b(k)] (containing the rhs column) at each point is, up to a "
              f"constant, {fls[0]}  and  {fls[1]}; gcd = {sp.factor(g.as_expr())}: for every k != 0 the system is inconsistent "
              f"at one of the two points  {tk(t1)}", sp.expand(g.monic().as_expr() - kk ** 2) == 0
              and all(p_[2] == 20 and p_[3] for p_ in pts))
        # at the nonzero roots of the first minor: own, referee's direct, and L's pull-back constructions at the first point
        t1 = time.time()
        Pm, y0 = pts[0][0], pts[0][1]
        pt = FibrePoint(y0, lam_of(eta))
        R1 = R1_fun(eta)
        out, ok = [], True
        for r_ in roots:
            kv = Fr(int(sp.numer(r_)), int(sp.denom(r_)))
            rows, rhs, _ = relaxed(pt, Kb, cas_bc(eta, kv))
            Ks, bcR = curv_Ks(R1, eta, -kv / 4)            # b_(R_1) = -4 c_A Cas, so c_A = -k/4 gives k Cas
            rowsR, rhsR, _ = relaxed(pt, Ks, bcR)
            mine = ranks(rows, rhs)
            rr, rh, okr, nu = X9.relaxed_system(eta, nest(R1), y0, -kv / 4, lam_of(eta) if -1 in eta else Fr(0))
            ref = X9.ranks(rr, rh, nu)
            PmS = sp.Matrix(4, 4, lambda i, j: sp.Rational(Pm[i][j]))
            Pi = PmS.T.inv()
            Bp = {(c, d): PmS.T * sp.Matrix(4, 4, lambda a, b: eta[a] * sp.Rational(R1(a, b, c, d))) * Pi for c in range(4) for d in range(4)}
            Rp = lambda a, b, c, d, Bp=Bp: L.fr(Bp[(c, d)][a, b] * eta[a])  # noqa: E731
            fd = L.fibre_data(eta, Rp, lam=lam_of(eta), maxorder=1, afac=-kv)
            unk, rowkeys, c0, _, rhsL = L.hessian_system(fd, 1, 3, mu0=Fr(1), mu1=Fr(0))
            _, _, c1, _, _ = L.hessian_system(fd, 1, 3, mu0=Fr(1), mu1=Fr(1))
            deg1 = [j for j, mn in enumerate(unk) if sum(mn) == 1]
            rest = [j for j, mn in enumerate(unk) if sum(mn) > 1]
            dcols = [{i: c1[j].get(i, 0) - c0[j].get(i, 0) for i in set(c0[j]) | set(c1[j]) if c1[j].get(i, 0) != c0[j].get(i, 0)}
                     for j in deg1]
            cols = [c0[j] for j in rest] + [c0[j] for j in deg1] + dcols
            lr = (L.dm_from_cols(len(rowkeys), cols).rank(), L.dm_from_cols(len(rowkeys), cols, rhsL).rank())
            ok &= rows == rowsR and rhs == rhsR and mine == (19, 20) and ref == (294, 295) and lr == (294, 295) and okr
            out.append((str(r_), mine, ref, lr))
        check(f"base {sl(eta)}: at the nonzero roots k of the first minor, at the first point: own system for k Cas (equal row "
              f"by row to the one for R = R_1, c_A = -k/4), the referee's direct Taylor system and M9's pull-back Taylor system "
              f"(L.fibre_data, afac = -k): (k, own, referee, M9) = {out}: inconsistent in all three  {tk(t1)}", ok and len(roots) > 0)
    # the value of M5(d) / R_L: k = 1/12
    t1 = time.time()
    eta = (1, -1, -1, -1)
    RL = X9.R_L()
    b = X9.BR_coeff(RL, eta)
    kgi = trace_form(kbasis(eta)).inv()
    okv = b * sp.Rational(1, 48) == kgi / 12
    check(f"base (1,3): for the referee's R_L (curvature operator^2 = -Id) at (c_V/c_H)^2 = 1/12 (c_A = 1/48), b_R = Cas/12, "
          f"i.e. the point k = 1/12 of the line, which is M5(d)'s (C1)_0 value; P4 excludes it on every open subset of the "
          f"fibre (xcheck_round9 section Z8 found it inconsistent at one point)  {tk(t1)}", okv)


def main():
    p0()
    p1()
    p3()
    p4()
    P(f"\n{NCHK[0]} checks, all True (each asserted).")
    P(f"runtime {time.time() - T0:.0f} s")


if __name__ == '__main__':
    main()
