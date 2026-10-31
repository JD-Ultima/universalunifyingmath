"""REVIEW11-GU: independent cross-check of GU round 12 (formalization/geometric-unity/P-z2-last-exception.md).

Referee construction, written for this review.  It does NOT import P_last_exception.py, L_nonradial.py, M_remaining.py,
J_slaved_sector.py or any xcheck_round* script.  The only repository import is F_jets.py, used once (section Y1) for a
brute-force Ricci tensor of the Z2 metric over a null-curvature base.

How it differs from the round-12 code:
  * the fibre geometry is built from the metric alone: G_lam(h,k) = tr(y^-1 h y^-1 k) - lam tr(y^-1 h) tr(y^-1 k) in the
    coordinates y_A of Sym^2(R^4); Christoffel symbols by Koszul from exact first derivatives of the metric, the Riemann
    tensor from the coordinate formula with exact second derivatives (not the closed forms Gamma = -(1/2)(h y^-1 k + ...)
    and R = -(1/4)[[X,Y],Z]);
  * Ric_F is the Ricci tensor of that Riemann tensor (F1 is checked, not assumed);
  * the A-term and its covariant derivative are computed by the coordinate product rule and the Koszul Christoffels (not
    the Killing-field formula);
  * expansions about the section are made directly at y(eps), as exact truncated power series in eps, along straight lines
    eta + eps h and along orbit curves (1 + eps X) eta (1 + eps X)^T (not by pulling the system back to eta);
  * so(eta) basis E_cd = R_1(e_c, e_d); lambda = -1 (base (4,0)) and lambda = 1 (indefinite bases); own seeds and points.

Notation as in the note: Sigma_y(b) has unknowns (V, xi) (20) and the 450 rows
     G(R(X,Y)V, Z) - xi(X) B(Y,Z) + xi(Y) B(X,Z) = (nabla_X A_b)(Y,Z) - (nabla_Y A_b)(X,Z),    B = Ric_F + A_b,
     A_b(U,W) = sum b^(st) G(U, xi_s) G(W, xi_t),   xi_K(y) = -(K^T y + y K),   b_R = c_A sum_(c != d) eps_c eps_d R_cd (x) R_cd.

Sections
  Y0  the construction: F1 (Ric_F = -G_1/4), G_1/4 parallel, derivative formulas against exact eps-series derivatives,
      a manufactured solution with mu' != 0 (tests every column), lambda-independence, the dilation kernel.
  Y1  Proposition P1: b_R = -4 c_A R^2 Cas; null curvatures (own examples): R != 0, R^2 = 0, b_R = 0, A-term zero, A != 0,
      relaxed system consistent with kernel = dilation only (the flat case); an element of Z outside the pp-wave family of
      P1(b) (observation for R12-3); brute-force O'Neill (F_jets) over a pp-wave
      base off the section: Ric^Y on VV equals Ric_F although A != 0; trace form definite in (4,0).
  Y2  Lemma P2: at random off-section points rank A = 19 and, for random R, (19, 20); at the section rank 19, A-term and
      its derivative zero.
  Y3  Lemma P3: first order along the ten coordinate lines through eta, b in Sym^2 so(eta) free: b-projection = R Cas.
  Y4  Lemma P4: on b = k Cas, third order along a line and an orbit curve (direct): orders 1, 2 pass, order 3 forces
      k = 0; the full truncated Taylor systems at five values of k; 20-minors in Q[k] at own off-section points, gcd k^2;
      M5(d)'s value k = 1/12 (R_L) and R_L off the section.
Exact arithmetic (Python Fractions, sympy rationals); fixed seeds random.Random(20261241..20261246); one core.
Run (from computation/geometric-unity): python3 xcheck_round12.py > xcheck_round12.out 2>&1   (~6 min)
"""
import itertools
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

T0 = time.time()
NCHK = [0, 0]
SIGS = ((1, 1, 1, 1), (1, -1, -1, -1), (-1, 1, 1, 1), (1, 1, -1, -1))
SYM = [(i, j) for i in range(4) for j in range(i, 4)]
NA = 10


def P(*a):
    print(*a)
    sys.stdout.flush()


def check(label, ok):
    NCHK[0] += 1
    NCHK[1] += bool(ok)
    P(f"  [{NCHK[0]:2d}] {label}: {bool(ok)}")
    assert ok, label


def section(s):
    P(f"\n=== {s}   [{time.time() - T0:.0f} s]")


def tk(t1):
    return f"[{time.time() - t1:.0f} s]"


def sl(eta):
    return f"({eta.count(1)},{eta.count(-1)})"


def lam_of(eta):
    return Fr(-1) if -1 not in eta else Fr(1)     # lambda < 1/4 for (14,0), > 1/4 for (7,7) (A4, A5); values not those of P


# ---------------------------------------------------------------- truncated power series in eps (lists of Fractions)
def s_add(a, b):
    return [x + y for x, y in zip(a, b)]


def s_sub(a, b):
    return [x - y for x, y in zip(a, b)]


def s_mul(a, b):
    N = len(a)
    return [sum(a[i] * b[n - i] for i in range(n + 1) if a[i] and b[n - i]) for n in range(N)]


def s_scal(a, c):
    return [x * c for x in a]


def s_zero(N):
    return [Fr(0)] * N


def s_const(c, N):
    return [Fr(c)] + [Fr(0)] * (N - 1)


# ---------------------------------------------------------------- 4x4 Fraction matrices and series of them
def mz():
    return [[Fr(0)] * 4 for _ in range(4)]


def mmul(A, B):
    return [[sum(A[i][k] * B[k][j] for k in range(4) if A[i][k] and B[k][j]) for j in range(4)] for i in range(4)]


def madd(A, B, s=1):
    return [[A[i][j] + s * B[i][j] for j in range(4)] for i in range(4)]


def mT(A):
    return [[A[j][i] for j in range(4)] for i in range(4)]


def minv(A):
    n = 4
    M = [list(map(Fr, A[i])) + [Fr(int(i == j)) for j in range(n)] for i in range(n)]
    for c in range(n):
        p = next(r for r in range(c, n) if M[r][c] != 0)
        M[c], M[p] = M[p], M[c]
        pv = M[c][c]
        M[c] = [x / pv for x in M[c]]
        for r in range(n):
            if r != c and M[r][c] != 0:
                f = M[r][c]
                M[r] = [x - f * y for x, y in zip(M[r], M[c])]
    return [row[n:] for row in M]


def EB(A):
    i, j = SYM[A]
    M = mz()
    M[i][j] = Fr(1)
    M[j][i] = Fr(1)
    return M


EBS = [EB(A) for A in range(NA)]


def coords(M):
    return [M[i][j] for (i, j) in SYM]


# series matrix: list of N Fraction matrices (coefficients of eps^0..eps^(N-1))
def pm_mul(A, B):
    N = len(A)
    out = []
    for n in range(N):
        acc = mz()
        for i in range(n + 1):
            acc = madd(acc, mmul(A[i], B[n - i]))
        out.append(acc)
    return out


def pm_tr(A):
    return [sum(C[i][i] for i in range(4)) for C in A]


def pm_trmul(A, B):
    """tr(A B) as a series, without forming A B"""
    N = len(A)
    out = []
    for n in range(N):
        s = Fr(0)
        for i in range(n + 1):
            Ai, Bj = A[i], B[n - i]
            for p in range(4):
                for q in range(4):
                    if Ai[p][q] and Bj[q][p]:
                        s += Ai[p][q] * Bj[q][p]
        out.append(s)
    return out


def pm_inv(Y):
    """inverse of a series matrix with invertible constant term (Neumann series)"""
    N = len(Y)
    Y0i = minv(Y[0])
    D = [mz()] + [mmul(Y0i, Y[n]) for n in range(1, N)]        # Y0^-1 (Y - Y0)
    # X = - D ; inverse = sum_k X^k Y0^-1
    term = [[[Fr(int(i == j)) for j in range(4)] for i in range(4)]] + [mz() for _ in range(N - 1)]
    acc = [r for r in term]
    negD = [[[-x for x in row] for row in M] for M in D]
    for _ in range(1, N):
        term = pm_mul(negD, term)
        acc = [madd(a, t) for a, t in zip(acc, term)]
    return [mmul(a, Y0i) for a in acc]


# ---------------------------------------------------------------- the fibre geometry at a point y (series in eps)
class Geo:
    """Fibre Met(R^4) with G_lam(h,k) = tr(y^-1 h y^-1 k) - lam tr(y^-1 h) tr(y^-1 k), at y (a series matrix),
       coordinates y_A on the basis EBS.  Christoffel symbols by Koszul from the metric derivatives,
       Riemann tensor from the coordinate formula R^a_bcd = d_c G^a_db - d_d G^a_cb + G^a_ce G^e_db - G^a_de G^e_cb."""

    def __init__(self, Y, lam, curv=True):
        self.N = N = len(Y)
        self.Y = Y
        self.lam = lam = Fr(lam)
        self.Yi = Yi = pm_inv(Y)
        self.M = M = [pm_mul(Yi, [EBS[A]] + [mz()] * (N - 1)) for A in range(NA)]       # M_A = y^-1 E_A
        self.t1 = t1 = [pm_tr(M[A]) for A in range(NA)]
        P2 = {(a, b): pm_mul(M[a], M[b]) for a in range(NA) for b in range(NA)}
        self.P2 = P2
        T2 = {(a, b): pm_tr(P2[(a, b)]) for a in range(NA) for b in range(NA)}
        T3 = {(a, b, c): pm_trmul(P2[(a, b)], M[c]) for a in range(NA) for b in range(NA) for c in range(NA)}
        self.T2, self.T3 = T2, T3
        g = [[s_sub(T2[(a, b)], s_scal(s_mul(t1[a], t1[b]), lam)) for b in range(NA)] for a in range(NA)]
        self.g = g
        # d_C g_AB
        dg = [[[None] * NA for _ in range(NA)] for _ in range(NA)]
        for C in range(NA):
            for A in range(NA):
                for B in range(A, NA):
                    v = s_scal(s_add(T3[(C, A, B)], T3[(A, C, B)]), -1)
                    v = s_add(v, s_scal(s_add(s_mul(T2[(C, A)], t1[B]), s_mul(t1[A], T2[(C, B)])), lam))
                    dg[C][A][B] = dg[C][B][A] = v
        self.dg = dg
        # inverse metric (series): g^-1 = sum (-g0^-1 D)^k g0^-1
        g0 = [[g[a][b][0] for b in range(NA)] for a in range(NA)]
        g0i = finv(g0)
        self.ginv = ser_inv10(g, g0i, N)
        # Koszul: Gamma_{E,CD} and Gamma^A_{CD}
        G1 = [[[s_scal(s_sub(s_add(dg[C][E][D], dg[D][E][C]), dg[E][C][D]), Fr(1, 2)) for D in range(NA)]
               for C in range(NA)] for E in range(NA)]
        self.G1 = G1
        gi = self.ginv
        Gam = [[[None] * NA for _ in range(NA)] for _ in range(NA)]
        for A in range(NA):
            for C in range(NA):
                for D in range(C, NA):
                    acc = s_zero(N)
                    for E in range(NA):
                        acc = s_add(acc, s_mul(gi[A][E], G1[E][C][D]))
                    Gam[A][C][D] = Gam[A][D][C] = acc
        self.Gam = Gam
        if curv:
            self._curv()

    def _curv(self):
        N, M, t1, T2, T3, lam, P2 = self.N, self.M, self.t1, self.T2, self.T3, self.lam, self.P2
        T4 = {}

        def t4(a, b, c, d):
            key = (a, b, c, d)
            if key not in T4:
                # cyclic invariance
                for r in ((b, c, d, a), (c, d, a, b), (d, a, b, c)):
                    if r in T4:
                        T4[key] = T4[r]
                        return T4[key]
                T4[key] = pm_trmul(P2[(a, b)], P2[(c, d)])
            return T4[key]
        ddg = {}
        for D in range(NA):
            for C in range(D, NA):
                for A in range(NA):
                    for B in range(A, NA):
                        v = s_zero(N)
                        for tup in ((D, C, A, B), (C, D, A, B), (C, A, D, B), (D, A, C, B), (A, D, C, B), (A, C, D, B)):
                            v = s_add(v, t4(*tup))
                        lt = s_mul(s_add(T3[(D, C, A)], T3[(C, D, A)]), t1[B])
                        lt = s_add(lt, s_mul(T2[(C, A)], T2[(D, B)]))
                        lt = s_add(lt, s_mul(T2[(D, A)], T2[(C, B)]))
                        lt = s_add(lt, s_mul(t1[A], s_add(T3[(D, C, B)], T3[(C, D, B)])))
                        v = s_sub(v, s_scal(lt, lam))
                        ddg[(D, C, A, B)] = ddg[(C, D, A, B)] = ddg[(D, C, B, A)] = ddg[(C, D, B, A)] = v
        self.ddg = ddg
        gi, dg, G1, Gam, g = self.ginv, self.dg, self.G1, self.Gam, self.g
        # d_B g^{AE} = -(g^-1 d_B g g^-1)^{AE}
        dgi = []
        for B in range(NA):
            X = [[None] * NA for _ in range(NA)]
            tmp = [[None] * NA for _ in range(NA)]
            for A in range(NA):
                for G in range(NA):
                    acc = s_zero(N)
                    for F in range(NA):
                        acc = s_add(acc, s_mul(gi[A][F], dg[B][F][G]))
                    tmp[A][G] = acc
            for A in range(NA):
                for E in range(NA):
                    acc = s_zero(N)
                    for G in range(NA):
                        acc = s_add(acc, s_mul(tmp[A][G], gi[G][E]))
                    X[A][E] = s_scal(acc, -1)
            dgi.append(X)
        # d_B Gamma^A_{CD}
        dGam = {}
        for B in range(NA):
            for A in range(NA):
                for C in range(NA):
                    for D in range(C, NA):
                        acc = s_zero(N)
                        for E in range(NA):
                            dG1 = s_scal(s_sub(s_add(ddg[(B, C, E, D)], ddg[(B, D, E, C)]), ddg[(B, E, C, D)]), Fr(1, 2))
                            acc = s_add(acc, s_add(s_mul(dgi[B][A][E], G1[E][C][D]), s_mul(gi[A][E], dG1)))
                        dGam[(B, A, C, D)] = dGam[(B, A, D, C)] = acc
        self.dGam = dGam
        # R^A_{BCD} (C < D) and R_{ZBCD} = g(R(d_C, d_D) d_B, d_Z)
        Rup = {}
        for A in range(NA):
            for B in range(NA):
                for C, D in itertools.combinations(range(NA), 2):
                    v = s_sub(dGam[(C, A, D, B)], dGam[(D, A, C, B)])
                    for E in range(NA):
                        v = s_add(v, s_sub(s_mul(Gam[A][C][E], Gam[E][D][B]), s_mul(Gam[A][D][E], Gam[E][C][B])))
                    Rup[(A, B, C, D)] = v
        Rlo = {}
        for Z in range(NA):
            for B in range(NA):
                for C, D in itertools.combinations(range(NA), 2):
                    acc = s_zero(N)
                    for A in range(NA):
                        acc = s_add(acc, s_mul(g[Z][A], Rup[(A, B, C, D)]))
                    Rlo[(Z, B, C, D)] = acc
        self.Rup, self.Rlo = Rup, Rlo
        # Ricci: Ric_{BD} = R^A_{BAD}
        Ric = [[None] * NA for _ in range(NA)]
        for B in range(NA):
            for D in range(NA):
                acc = s_zero(N)
                for A in range(NA):
                    if A < D:
                        acc = s_add(acc, Rup[(A, B, A, D)])
                    elif A > D:
                        acc = s_sub(acc, Rup[(A, B, D, A)])
                Ric[B][D] = acc
        self.Ric = Ric

    def G(self, h, k, lam=None):
        """G_lam(h, k) at y for series matrices h, k"""
        lam = self.lam if lam is None else Fr(lam)
        a, b = pm_mul(self.Yi, h), pm_mul(self.Yi, k)
        return s_sub(pm_trmul(a, b), s_scal(s_mul(pm_tr(a), pm_tr(b)), lam))

    def dG(self, C, h, k, lam=None):
        """(d_C G)(h, k) at y (the metric's coordinate derivative, h, k fixed)"""
        lam = self.lam if lam is None else Fr(lam)
        a, b = pm_mul(self.Yi, h), pm_mul(self.Yi, k)
        MC = self.M[C]
        v = s_scal(s_add(pm_trmul(pm_mul(MC, a), b), pm_trmul(a, pm_mul(MC, b))), -1)
        v = s_add(v, s_scal(s_add(s_mul(pm_trmul(MC, a), pm_tr(b)), s_mul(pm_tr(a), pm_trmul(MC, b))), lam))
        return v


def finv(A):
    n = len(A)
    M = [list(map(Fr, A[i])) + [Fr(int(i == j)) for j in range(n)] for i in range(n)]
    for c in range(n):
        p = next(r for r in range(c, n) if M[r][c] != 0)
        M[c], M[p] = M[p], M[c]
        pv = M[c][c]
        M[c] = [x / pv for x in M[c]]
        for r in range(n):
            if r != c and M[r][c] != 0:
                f = M[r][c]
                M[r] = [x - f * y for x, y in zip(M[r], M[c])]
    return [row[n:] for row in M]


def ser_inv10(g, g0i, N):
    n = len(g)
    if N == 1:
        return [[[g0i[a][b]] for b in range(n)] for a in range(n)]
    # D_k = g0^-1 g_k  (k >= 1); inverse = sum_m (-D)^m g0^-1 as series
    Dk = [None] + [[[sum(g0i[a][c] * g[c][b][k] for c in range(n)) for b in range(n)] for a in range(n)] for k in range(1, N)]
    # series matrix X (coefficient list): X_0 = I, X_n = -sum_{k=1..n} D_k X_{n-k}
    X = [[[Fr(int(a == b)) for b in range(n)] for a in range(n)]]
    for m in range(1, N):
        acc = [[Fr(0)] * n for _ in range(n)]
        for k in range(1, m + 1):
            Dm, Xp = Dk[k], X[m - k]
            for a in range(n):
                for b in range(n):
                    acc[a][b] -= sum(Dm[a][c] * Xp[c][b] for c in range(n))
        X.append(acc)
    out = [[[sum(X[m][a][c] * g0i[c][b] for c in range(n)) for m in range(N)] for b in range(n)] for a in range(n)]
    return out


# ---------------------------------------------------------------- so(eta), Killing fields, the A-term, the relaxed rows
PAIRS = [(c, d) for c in range(4) for d in range(c + 1, 4)]


def so_basis(eta):
    """E_cd = R_1(e_c, e_d): (E_cd)^a_b = delta^a_c eta_db - delta^a_d eta_cb"""
    out = []
    for c, d in PAIRS:
        K = mz()
        K[c][d] = Fr(eta[d])
        K[d][c] = Fr(-eta[c])
        out.append(K)
    return out


def killing(K, y):
    return [[-(sum(K[k][i] * y[k][j] for k in range(4)) + sum(y[i][k] * K[k][j] for k in range(4))) for j in range(4)]
            for i in range(4)]


def ser(Mx, N):
    return [Mx] + [mz() for _ in range(N - 1)]


def kill_data(geo, Ks):
    N = geo.N
    xs = [[killing(K, Yn) for Yn in geo.Y] for K in Ks]
    Gx = [[geo.G(ser(EBS[A], N), xs[s]) for A in range(NA)] for s in range(len(Ks))]
    dGx = [[[s_add(geo.dG(C, ser(EBS[A], N), xs[s]), geo.G(ser(EBS[A], N), ser(killing(Ks[s], EBS[C]), N)))
             for A in range(NA)] for C in range(NA)] for s in range(len(Ks))]
    return Gx, dGx


def aterm(geo, Ks, beta, kd=None):
    """A_AB = sum beta[s,t] G(E_A, xi_s) G(E_B, xi_t) at y, d_C A_AB by the product rule (coordinate derivative of the
       metric), and (nabla_C A)_AB with the Koszul Christoffels.  beta: dict (s,t) -> c.  Returns A, nA, dA."""
    N = geo.N
    Gx, dGx = kill_data(geo, Ks) if kd is None else kd
    A = [[s_zero(N) for _ in range(NA)] for _ in range(NA)]
    dA = [[[s_zero(N) for _ in range(NA)] for _ in range(NA)] for _ in range(NA)]
    for (s, t), c in beta.items():
        if not c:
            continue
        for a in range(NA):
            for b in range(NA):
                A[a][b] = s_add(A[a][b], s_scal(s_mul(Gx[s][a], Gx[t][b]), c))
                for C in range(NA):
                    dA[C][a][b] = s_add(dA[C][a][b], s_scal(s_add(s_mul(dGx[s][C][a], Gx[t][b]), s_mul(Gx[s][a], dGx[t][C][b])), c))
    Gam = geo.Gam
    nA = [[[None] * NA for _ in range(NA)] for _ in range(NA)]
    for C in range(NA):
        for a in range(NA):
            for b in range(NA):
                v = dA[C][a][b]
                for D in range(NA):
                    v = s_sub(v, s_add(s_mul(Gam[D][C][a], A[D][b]), s_mul(Gam[D][C][b], A[a][D])))
                nA[C][a][b] = v
    return A, nA, dA


def rhs_of(nA):
    return [s_sub(nA[X][Y][Z], nA[Y][X][Z]) for X, Y in itertools.combinations(range(NA), 2) for Z in range(NA)]


def xi_rows(A, N):
    """the A-term's contribution to the xi-columns: row (X<Y, Z): col 10+X: -A_YZ, col 10+Y: +A_XZ"""
    rows = []
    for X, Y in itertools.combinations(range(NA), 2):
        for Z in range(NA):
            row = [s_zero(N) for _ in range(2 * NA)]
            row[10 + X] = s_scal(A[Y][Z], -1)
            row[10 + Y] = A[X][Z]
            rows.append(row)
    return rows


def relaxed_rows(geo, A, nA):
    """rows (X<Y, Z): sum_A V^A R_{Z A X Y} - xi_X B_YZ + xi_Y B_XZ = (nabla_X A)_YZ - (nabla_Y A)_XZ, B = Ric_F + A.
       Each row: list of 20 series; rhs: series.  A, nA may be None (A-term off)."""
    N = geo.N
    rows, rhs = [], []
    for X, Y in itertools.combinations(range(NA), 2):
        for Z in range(NA):
            row = [geo.Rlo[(Z, Aa, X, Y)] for Aa in range(NA)] + [s_zero(N) for _ in range(NA)]
            BYZ = geo.Ric[Y][Z] if A is None else s_add(geo.Ric[Y][Z], A[Y][Z])
            BXZ = geo.Ric[X][Z] if A is None else s_add(geo.Ric[X][Z], A[X][Z])
            row[10 + X] = s_scal(BYZ, -1)
            row[10 + Y] = BXZ
            rows.append(row)
            rhs.append(s_zero(N) if nA is None else s_sub(nA[X][Y][Z], nA[Y][X][Z]))
    return rows, rhs


def order(rows, n):
    return [[e[n] for e in r] for r in rows]


def orderv(v, n):
    return [e[n] for e in v]


# ---------------------------------------------------------------- exact linear algebra over Q
def qq(v):
    return sp.QQ(v.numerator, v.denominator)


def dm(rows, ncols):
    d = {}
    for i, r in enumerate(rows):
        e = {j: qq(v) for j, v in enumerate(r) if v}
        if e:
            d[i] = e
    return DomainMatrix(d, (len(rows), ncols), sp.QQ)


def rank(rows):
    return dm(rows, len(rows[0])).rank()


def augment(rows, *vecs):
    return [list(r) + [v[i] for v in vecs] for i, r in enumerate(rows)]


def solve(rows, vec):
    """one exact solution of rows x = vec (free variables 0), or None if inconsistent"""
    n = len(rows[0])
    R, piv = dm(augment(rows, vec), n + 1).to_Matrix().rref()
    if n in piv:
        return None
    x = [Fr(0)] * n
    for i, p in enumerate(piv):
        x[p] = Fr(int(R[i, n].p), int(R[i, n].q))
    return x


def matvec(rows, x):
    return [sum(a * b for a, b in zip(r, x) if a and b) for r in rows]


def fdet(M):
    M = [list(r) for r in M]
    n, d = len(M), Fr(1)
    for c in range(n):
        p = next((r for r in range(c, n) if M[r][c] != 0), None)
        if p is None:
            return Fr(0)
        if p != c:
            M[c], M[p] = M[p], M[c]
            d = -d
        d *= M[c][c]
        for r in range(c + 1, n):
            if M[r][c] != 0:
                f = M[r][c] / M[c][c]
                M[r] = [x - f * y for x, y in zip(M[r], M[c])]
    return d


def pick_rows(M, k):
    """indices of the first k linearly independent rows of M (greedy elimination)"""
    basis, piv, sel = [], [], []
    for i, r in enumerate(M):
        v = list(r)
        for pc, br in zip(piv, basis):
            if v[pc] != 0:
                f = v[pc] / br[pc]
                v = [x - f * y for x, y in zip(v, br)]
        nz = next((j for j, x in enumerate(v) if x != 0), None)
        if nz is not None:
            basis.append(v)
            piv.append(nz)
            sel.append(i)
            if len(sel) == k:
                break
    return sel


# ---------------------------------------------------------------- base curvature tensors
def diag(eta):
    return [[Fr(eta[i]) if i == j else Fr(0) for j in range(4)] for i in range(4)]


def mtr(A):
    return sum(A[i][i] for i in range(4))


def from_pairs(S):
    R = [[[[Fr(0)] * 4 for _ in range(4)] for _ in range(4)] for _ in range(4)]
    for p, (a, b) in enumerate(PAIRS):
        for q, (c, d) in enumerate(PAIRS):
            for aa, bb, s1 in ((a, b, 1), (b, a, -1)):
                for cc, dd, s2 in ((c, d, 1), (d, c, -1)):
                    R[aa][bb][cc][dd] = s1 * s2 * S[p][q]
    return R


def rand_curv(rng):
    """random algebraic curvature tensor R_abcd (all indices down): symmetric pair matrix, first Bianchi imposed on the
       single totally antisymmetric component (pairs ordered 01, 02, 03, 12, 13, 23)"""
    S = [[None] * 6 for _ in range(6)]
    for i in range(6):
        for j in range(i, 6):
            S[i][j] = S[j][i] = Fr(rng.randint(-3, 3), rng.randint(1, 3))
    S[2][3] = S[3][2] = S[1][4] - S[0][5]
    return from_pairs(S)


def null_curv(eta, l, es, H):
    """R = sum H_ij alpha_i (x) alpha_j, alpha_i = l ^ e_i (indices lowered with eta)"""
    low = lambda v: [eta[a] * Fr(v[a]) for a in range(4)]  # noqa: E731
    al = [[[low(l)[a] * low(e)[b] - low(l)[b] * low(e)[a] for b in range(4)] for a in range(4)] for e in es]
    n = len(es)
    return [[[[sum(Fr(H[i][j]) * al[i][a][b] * al[j][c][d] for i in range(n) for j in range(n)) for d in range(4)]
              for c in range(4)] for b in range(4)] for a in range(4)]


def two_form(entries):
    F = [[Fr(0)] * 4 for _ in range(4)]
    for a, b, v in entries:
        F[a][b] += v
        F[b][a] -= v
    return F


def ip2(F, G, eta):
    """(1/2) F_ab G^ab"""
    return sum(F[a][b] * G[a][b] * eta[a] * eta[b] for a in range(4) for b in range(4)) / 2


def wedge2(F, G):
    """coefficient of e^0123 in F ^ G = (1/4) sum over permutations sgn(p) F_p0p1 G_p2p3"""
    s = Fr(0)
    for p in itertools.permutations(range(4)):
        inv = sum(1 for i in range(4) for j in range(i + 1, 4) if p[i] > p[j])
        s += (-1) ** inv * F[p[0]][p[1]] * G[p[2]][p[3]]
    return s / 4


def R_L():
    """REVIEW7-GU's Lorentzian tensor: R_(01,12) = R_(02,23) = R_(03,31) = 1, other pair components 0"""
    S = [[Fr(0)] * 6 for _ in range(6)]
    pi = {p: i for i, p in enumerate(PAIRS)}

    def setp(a, b, c, d, v):
        sa, sb = (1 if a < b else -1), (1 if c < d else -1)
        i, j = pi[tuple(sorted((a, b)))], pi[tuple(sorted((c, d)))]
        S[i][j] = S[j][i] = Fr(v * sa * sb)
    setp(0, 1, 1, 2, 1)
    setp(0, 2, 2, 3, 1)
    setp(0, 3, 3, 1, 1)
    return from_pairs(S)


def bianchi(R):
    return all(R[a][b][c][d] + R[a][c][d][b] + R[a][d][b][c] == 0 for a, b, c, d in itertools.product(range(4), repeat=4))


def nonzero(R):
    return any(R[a][b][c][d] for a, b, c, d in itertools.product(range(4), repeat=4))


def ricci(R, eta):
    return [[sum(Fr(eta[a]) * R[a][b][a][d] for a in range(4)) for d in range(4)] for b in range(4)]


def Rcd_list(R, eta):
    """(c, d) and the endomorphism R_cd = R(e_c, e_d), (R_cd)^a_b = eta^aa R_abcd, for the 12 ordered pairs c != d"""
    return [((c, d), [[Fr(eta[a]) * R[a][b][c][d] for b in range(4)] for a in range(4)])
            for c in range(4) for d in range(4) if c != d]


def curv_op(R, eta):
    """matrix (columns) of E_cd -> R_cd in the basis so_basis(eta); E_ab[a][b] = eta_b"""
    out = sp.zeros(6, 6)
    for q, (c, d) in enumerate(PAIRS):
        for p, (a, b) in enumerate(PAIRS):
            out[p, q] = sp.Rational(Fr(eta[a]) * R[a][b][c][d] / eta[b])
    return out


def vec16(K):
    return sp.Matrix([sp.Rational(K[i][j]) for i in range(4) for j in range(4)])


def BR16(R, eta):
    """sum_(c != d) eps_c eps_d R_cd (x) R_cd as a 16 x 16 matrix"""
    M = sp.zeros(16, 16)
    for (c, d), K in Rcd_list(R, eta):
        v = vec16(K)
        M += eta[c] * eta[d] * v * v.T
    return M


def kg_of(eta):
    Ks = so_basis(eta)
    return sp.Matrix(6, 6, lambda p, q: sp.Rational(mtr(mmul(Ks[p], Ks[q]))))


def opCas16(Op, eta):
    """(Op (x) 1) Cas as a 16 x 16 matrix, Cas = sum (kg^-1)^pq E_p (x) E_q, Op a 6 x 6 matrix on so(eta)"""
    Ks = so_basis(eta)
    kgi = kg_of(eta).inv()
    Kv = [vec16(K) for K in Ks]
    M = sp.zeros(16, 16)
    for p in range(6):
        img = sum((Op[r, p] * Kv[r] for r in range(6)), sp.zeros(16, 1))
        for q in range(6):
            if kgi[p, q] != 0:
                M += kgi[p, q] * img * Kv[q].T
    return M


def cas_beta(eta, scale=Fr(1)):
    kgi = kg_of(eta).inv()
    return {(p, q): scale * Fr(int(kgi[p, q].p), int(kgi[p, q].q)) for p in range(6) for q in range(6) if kgi[p, q] != 0}


def R_beta(R, eta, cA):
    """(Ks, beta) with A_x = A_beta: the twelve R_cd with coefficients c_A eps_c eps_d"""
    lst = Rcd_list(R, eta)
    return [K for _, K in lst], {(s, s): cA * eta[c] * eta[d] for s, ((c, d), _) in enumerate(lst)}


def randP(rng, num=2, dens=(2, 3)):
    while True:
        Pm = [[Fr(rng.randint(-num, num), rng.choice(dens)) + (1 if i == j else 0) for j in range(4)] for i in range(4)]
        if fdet(Pm) != 0:
            return Pm


def offpoint(rng, eta):
    Pm = randP(rng)
    return mmul(mmul(Pm, diag(eta)), mT(Pm))


def system_at(y0, lam, Ks, beta):
    """the relaxed system at the exact point y0: (A rows, rhs, A-term, nabla A-term, geometry)"""
    geo = Geo([y0], lam)
    A, nA, _ = aterm(geo, Ks, beta)
    rows, rhs = relaxed_rows(geo, A, nA)
    return order(rows, 0), orderv(rhs, 0), A, nA, geo


# ====================================================================== Y0
def manufactured(geo, rng):
    """w = l1 l2 l3 + tr(Q y Q' y) (l_i linear), mu(t) = 1 + t; B := Hess w / mu(w).  Then (V, xi) = (grad w/mu, dw/mu)
       must solve  G(R(X,Y)V,Z) - xi_X B_YZ + xi_Y B_XZ = (nabla_X B)_YZ - (nabla_Y B)_XZ  at y (all 450 rows)."""
    y0 = geo.Y[0]
    L = [[[Fr(rng.randint(-2, 2)) for _ in range(4)] for _ in range(4)] for _ in range(3)]
    Q1 = [[Fr(rng.randint(-2, 2)) for _ in range(4)] for _ in range(4)]
    Q2 = [[Fr(rng.randint(-2, 2)) for _ in range(4)] for _ in range(4)]
    lv = [mtr(mmul(Li, y0)) for Li in L]
    ld = [[mtr(mmul(Li, EBS[A])) for A in range(NA)] for Li in L]
    w = lv[0] * lv[1] * lv[2] + mtr(mmul(mmul(Q1, y0), mmul(Q2, y0)))
    perm = list(itertools.permutations(range(3)))
    d1 = [ld[0][A] * lv[1] * lv[2] + lv[0] * ld[1][A] * lv[2] + lv[0] * lv[1] * ld[2][A]
          + mtr(mmul(mmul(Q1, EBS[A]), mmul(Q2, y0))) + mtr(mmul(mmul(Q1, y0), mmul(Q2, EBS[A]))) for A in range(NA)]
    d2 = [[sum(ld[p[0]][A] * ld[p[1]][B] * lv[p[2]] for p in perm)
           + mtr(mmul(mmul(Q1, EBS[A]), mmul(Q2, EBS[B]))) + mtr(mmul(mmul(Q1, EBS[B]), mmul(Q2, EBS[A])))
           for B in range(NA)] for A in range(NA)]
    d3 = [[[sum(ld[p[0]][A] * ld[p[1]][B] * ld[p[2]][C] for p in perm) for C in range(NA)] for B in range(NA)]
          for A in range(NA)]
    mu, mup = 1 + w, Fr(1)
    if mu == 0:
        return False
    Gm = lambda a, b, c: geo.Gam[a][b][c][0]  # noqa: E731
    dGm = lambda D, a, b, c: geo.dGam[(D, a, b, c)][0]  # noqa: E731
    H = [[d2[A][B] - sum(Gm(C, A, B) * d1[C] for C in range(NA)) for B in range(NA)] for A in range(NA)]
    dH = [[[d3[D][A][B] - sum(dGm(D, C, A, B) * d1[C] + Gm(C, A, B) * d2[D][C] for C in range(NA)) for B in range(NA)]
           for A in range(NA)] for D in range(NA)]
    Bm = [[H[A][B] / mu for B in range(NA)] for A in range(NA)]
    dB = [[[dH[D][A][B] / mu - H[A][B] * mup * d1[D] / mu ** 2 for B in range(NA)] for A in range(NA)] for D in range(NA)]
    nB = [[[dB[D][A][B] - sum(Gm(E, D, A) * Bm[E][B] + Gm(E, D, B) * Bm[A][E] for E in range(NA)) for B in range(NA)]
           for A in range(NA)] for D in range(NA)]
    gi = [[geo.ginv[a][b][0] for b in range(NA)] for a in range(NA)]
    V = [sum(gi[A][C] * d1[C] for C in range(NA)) / mu for A in range(NA)]
    xi = [mup * d1[A] / mu for A in range(NA)]
    ok, nz = True, False
    for X, Y in itertools.combinations(range(NA), 2):
        for Z in range(NA):
            lhs = sum(geo.Rlo[(Z, A, X, Y)][0] * V[A] for A in range(NA)) - xi[X] * Bm[Y][Z] + xi[Y] * Bm[X][Z]
            r = nB[X][Y][Z] - nB[Y][X][Z]
            ok &= lhs == r
            nz |= r != 0
    return ok and nz and any(xi)


def y0sec():
    section("Y0. The construction from the metric alone: F1, derivatives, a manufactured solution, lambda-independence")
    rng = random.Random(20261241)
    for eta in SIGS:
        t1 = time.time()
        lam = lam_of(eta)
        y0 = offpoint(rng, eta)
        geo = Geo([y0], lam)
        s1 = lambda M: ser(M, 1)  # noqa: E731
        okF1 = all(geo.Ric[a][b][0] == -geo.G(s1(EBS[a]), s1(EBS[b]), Fr(1, 4))[0] for a in range(NA) for b in range(NA))
        okpar = True
        for C in range(NA):
            for a in range(NA):
                for b in range(a, NA):
                    v = geo.dG(C, s1(EBS[a]), s1(EBS[b]), Fr(1, 4))[0]
                    for D in range(NA):
                        v -= (geo.Gam[D][C][a][0] * geo.G(s1(EBS[D]), s1(EBS[b]), Fr(1, 4))[0]
                              + geo.Gam[D][C][b][0] * geo.G(s1(EBS[a]), s1(EBS[D]), Fr(1, 4))[0])
                    okpar &= v == 0
        # derivative formulas against exact eps-series derivatives along E_C and along E_C + E_D
        R = rand_curv(rng)
        Ks, beta = R_beta(R, eta, Fr(1, 3))
        A, nA, dA = aterm(geo, Ks, beta)
        C, D = rng.sample(range(NA), 2)
        gC = Geo([y0, EBS[C]], lam, curv=False)
        AC, _, _ = aterm(gC, Ks, beta)
        okd = all(gC.g[a][b][1] == geo.dg[C][a][b][0] for a in range(NA) for b in range(NA))
        okd &= all(AC[a][b][1] == dA[C][a][b][0] for a in range(NA) for b in range(NA))
        gCD = Geo([y0, madd(EBS[C], EBS[D]), mz()], lam, curv=False)
        dd = lambda p, q, a, b: geo.ddg[(p, q, a, b)][0]  # noqa: E731
        okd &= all(2 * gCD.g[a][b][2] == dd(C, C, a, b) + 2 * dd(C, D, a, b) + dd(D, D, a, b) for a in range(NA) for b in range(NA))
        okm = manufactured(geo, rng)
        # lambda-independence of the relaxed system (A-term of the random R) and the dilation kernel
        rows, rhs, _, _, _ = system_at(y0, lam, Ks, beta)
        rows2, rhs2, _, _, _ = system_at(y0, Fr(1, 3), Ks, beta)
        okl = rows == rows2 and rhs == rhs2
        okk = all(v == 0 for v in matvec(rows, coords(y0) + [Fr(0)] * 10))
        check(f"base {sl(eta)}, random fibre point y0 = P eta P^T, lambda = {lam}: Ric of the Koszul/coordinate Riemann "
              f"tensor = -G_1/4 (F1): {okF1}; G_1/4 parallel: {okpar}; d_C g, d_C d_D g and d_C A (product rule) equal the "
              f"exact eps-series derivatives along E_C, E_C + E_D: {okd}; manufactured solution (w cubic + quadratic in y, "
              f"mu = 1 + w, B := Hess w / mu): (grad w/mu, dw/mu) solves all 450 rows, xi != 0, rhs != 0: {okm}; relaxed "
              f"system identical at lambda = {lam} and 1/3: {okl}; dilation (y, 0) in the kernel: {okk}  {tk(t1)}",
              okF1 and okpar and okd and okm and okl and okk)
    # F1 as a series along a line (to order eps^3)
    t1 = time.time()
    eta = (1, 1, -1, -1)
    h = [[Fr(rng.randint(-2, 2)) for _ in range(4)] for _ in range(4)]
    h = [[h[i][j] + h[j][i] for j in range(4)] for i in range(4)]
    geo = Geo([diag(eta), h, mz(), mz()], lam_of(eta))
    sN = lambda M: ser(M, 4)  # noqa: E731
    ok = all(geo.Ric[a][b] == s_scal(geo.G(sN(EBS[a]), sN(EBS[b]), Fr(1, 4)), -1) for a in range(NA) for b in range(NA))
    check(f"base (2,2): along the line eta + eps h (h random) the Koszul/coordinate Ricci tensor equals -G_1/4 as an exact "
          f"power series to order eps^3  {tk(t1)}", ok)


# ====================================================================== Y1
def y1sec():
    section("Y1. Proposition P1: b_R = -4 c_A R^2 Cas; null curvatures; O'Neill by brute force over a pp-wave base")
    rng = random.Random(20261242)
    for eta in SIGS:
        t1 = time.time()
        R = rand_curv(rng)
        Rm = curv_op(R, eta)
        kg = kg_of(eta)
        sa = (kg * Rm).T == kg * Rm
        idn = BR16(R, eta) == -4 * opCas16(Rm * Rm, eta)
        # the A-term of the twelve R_cd equals the A-term of b_R written in the so(eta) basis
        y0 = offpoint(rng, eta)
        geo = Geo([y0], lam_of(eta), curv=False)
        Ks, beta = R_beta(R, eta, Fr(1, 4))
        A1, n1, _ = aterm(geo, Ks, beta)
        bR = -Rm * Rm * kg.inv()                   # b_R / c_A = -4 R^2 Cas, c_A = 1/4
        beta2 = {(p, q): Fr(int(bR[p, q].p), int(bR[p, q].q)) for p in range(6) for q in range(6) if bR[p, q] != 0}
        A2, n2, _ = aterm(geo, so_basis(eta), beta2)
        okA = A1 == A2 and n1 == n2 and any(v != [0] for r in A1 for v in r)
        check(f"base {sl(eta)}, random R: Bianchi {bianchi(R)}; curvature operator self-adjoint for the trace form: {sa}; "
              f"sum_(c != d) eps_c eps_d R_cd (x) R_cd = -4 (R^2 (x) 1) Cas as 16 x 16 matrices: {idn}; at a random fibre "
              f"point the A-term and its covariant derivative from the twelve R_cd equal those from b_R = -4 c_A R^2 Cas: "
              f"{okA}  {tk(t1)}", bianchi(R) and sa and idn and okA)
    cases = [((1, -1, -1, -1), (1, Fr(3, 5), Fr(4, 5), 0), [(0, Fr(4, 5), Fr(-3, 5), 0), (0, 0, 0, 1)], [[2, -1], [-1, 5]]),
             ((-1, 1, 1, 1), (1, Fr(3, 5), Fr(4, 5), 0), [(0, Fr(4, 5), Fr(-3, 5), 0), (0, 0, 0, 1)], [[2, -1], [-1, 5]]),
             ((1, 1, -1, -1), (Fr(3, 5), Fr(4, 5), 1, 0), [(Fr(4, 5), Fr(-3, 5), 0, 0), (0, 0, 0, 1)], [[2, -1], [-1, 5]]),
             ((1, 1, -1, -1), (1, 0, 1, 0), [(0, 1, 0, 1)], [[1]])]
    for eta, l, es, H in cases:
        t1 = time.time()
        lnull = sum(eta[a] * Fr(l[a]) ** 2 for a in range(4)) == 0
        perp = all(sum(eta[a] * Fr(l[a]) * Fr(e[a]) for a in range(4)) == 0 for e in es)
        R = null_curv(eta, l, es, H)
        Rm = curv_op(R, eta)
        ric = ricci(R, eta)
        ricnz = any(v for r in ric for v in r)
        ok0 = bianchi(R) and nonzero(R) and Rm != sp.zeros(6, 6) and Rm * Rm == sp.zeros(6, 6) and BR16(R, eta) == sp.zeros(16, 16)
        Ks, beta = R_beta(R, eta, Fr(1, 4))
        res = []
        okf = True
        for _ in range(2):
            y0 = offpoint(rng, eta)
            rows, rhs, A, nA, _ = system_at(y0, lam_of(eta), Ks, beta)
            Anz = any(any(x for r in killing(K, y0) for x in r) for K in Ks)        # A_X Y ~ (R(X,Y).y)^vert != 0
            okf &= all(v == [0] for r in A for v in r) and not any(rhs) and Anz
            ra = rank(rows)
            okf &= ra == 19 and all(v == 0 for v in matvec(rows, coords(y0) + [Fr(0)] * 10))
            res.append((ra, rank(augment(rows, rhs))))
        okf &= all(r == (19, 19) for r in res)
        lab = 'l ^ m, l, m null and orthogonal (totally null plane)' if len(es) == 1 else 'sum H_ij (l ^ e_i) (x) (l ^ e_j), e_i orthogonal to l'
        check(f"base {sl(eta)}: R = {lab}, l = {tuple(str(x) for x in l)} null: {lnull and perp}; Bianchi, R != 0, "
              f"(curvature operator)^2 = 0, b_R = 0 as a tensor (so A_x = 0 on the whole fibre): {ok0}; Ric != 0: {ricnz}; "
              f"at two random fibre points xi_(R_cd)(y) != 0 (A != 0) but the A-term and the rhs vanish, and the relaxed "
              f"system has ranks {res} with kernel exactly the dilation (so grad w in R d_s and mu' dw = 0: M9's flat case "
              f"applies verbatim)  {tk(t1)}", lnull and perp and ok0 and okf and (ricnz or len(es) == 1))
    # observation (REVIEW11-GU R12-3): Z \ {0} is larger than the pp-wave family of P1(b) in Lorentzian bases
    t1 = time.time()
    res, ok = [], True
    for eta in ((1, -1, -1, -1), (-1, 1, 1, 1)):
        F1 = two_form([(0, 1, 1), (2, 3, 1)])           # e^01 + e^23  (E = B = e_1)
        F2 = two_form([(0, 2, 1), (3, 1, 1)])           # e^02 + e^31  (E = B = e_2)
        R = [[[[F1[a][b] * F1[c][d] - F2[a][b] * F2[c][d] for d in range(4)] for c in range(4)] for b in range(4)]
             for a in range(4)]
        Rm = curv_op(R, eta)
        ric = ricci(R, eta)
        gram_tr = [[ip2(F, G, eta) for G in (F1, F2)] for F in (F1, F2)]
        gram_w = [[wedge2(F, G) for G in (F1, F2)] for F in (F1, F2)]
        ok &= bianchi(R) and nonzero(R) and Rm * Rm == sp.zeros(6, 6) and BR16(R, eta) == sp.zeros(16, 16)
        ok &= gram_tr == [[0, 0], [0, 0]] and gram_w == [[2, 0], [0, 2]]          # Im R totally null, no simple element
        rk_ric = sp.Matrix(4, 4, lambda i, j: sp.Rational(ric[i][j])).rank()
        ok &= rk_ric == 2                                                          # Ric not of the form f l (x) l
        Ks, beta = R_beta(R, eta, Fr(1, 4))
        y0 = offpoint(rng, eta)
        rows, rhs, _, _, _ = system_at(y0, lam_of(eta), Ks, beta)
        r = (rank(rows), rank(augment(rows, rhs)))
        ok &= r == (19, 19)
        res.append((sl(eta), [str(ric[i][i]) for i in range(4)], r))
    check(f"observation (R12-3), bases (1,3), (3,1): R = F1 (x) F1 - F2 (x) F2, F1 = e^01 + e^23, F2 = e^02 + e^31: Bianchi, "
          f"R != 0, (curvature operator)^2 = 0, b_R = 0, so R is in Z; but Im R = span(F1, F2) is totally null with "
          f"F ^ F = 2(a^2 + b^2) vol for F = a F1 + b F2 (no simple element, unlike l ^ e_i), and Ric (diagonal "
          f"{[x[1] for x in res]}) has rank 2, not f l (x) l: Z \\ {{0}} is not only the pp-wave type of P1(b); relaxed system "
          f"at a random fibre point {[x[2] for x in res]}  {tk(t1)}", ok)
    # brute-force O'Neill over a pp-wave base (F_jets): Ric^Y on VV = Ric_F at a point off the section although A != 0
    t1 = time.time()
    eta = (1, -1, -1, -1)
    R = null_curv(eta, *cases[0][1:])
    ctx = FJ.Ctx('q')
    g0 = [[ctx.c(Fr(eta[i]) if i == j else Fr(0)) for j in range(4)] for i in range(4)]
    lin = [[[ctx.zero] * 4 for _ in range(4)] for _ in range(4)]
    quad = {}
    for k in range(4):
        for l_ in range(k, 4):
            M = [[ctx.c(Fr(-(R[i][k][j][l_] + (R[i][l_][j][k] if k != l_ else 0)), 3)) for j in range(4)] for i in range(4)]
            if any(x != 0 for r in M for x in r):
                quad[(k, l_)] = M
    gj = FJ.base_from_coeffs(g0, lin, quad, {}, 4, ctx)
    base = FJ.base_geometry(gj, 4, ctx)
    okR = all(base['R'][(a, b, c, d)] == ctx.c(R[a][b][c][d]) for a, b, c, d in itertools.product(range(4), repeat=4))
    Pm = [[Fr(1), Fr(1, 3), Fr(0), Fr(-1, 2)], [Fr(0), Fr(1), Fr(2, 3), Fr(0)], [Fr(1, 4), Fr(0), Fr(1), Fr(1, 3)],
          [Fr(0), Fr(-1, 3), Fr(0), Fr(1)]]
    y0 = mmul(mmul(Pm, diag(eta)), mT(Pm))
    lam, cH, cV = Fr(1, 2), Fr(-2), Fr(3)
    G, m, idx, Psi0 = FJ.total_metric(gj, base, 4, [[ctx.c(y0[i][j]) for j in range(4)] for i in range(4)], ctx.c(lam),
                                      ctx.c(cH), ctx.c(cV), 'Z2', ctx)
    Ric, _, _ = FJ.ricci_at_point(G, m, ctx)
    Tb = FJ.frame_blocks(Ric, 4, idx, Psi0, ctx)
    geo = Geo([y0], lam, curv=False)
    okVV = all(Tb[4 + a][4 + b] == ctx.c(-geo.G(ser(EBS[a], 1), ser(EBS[b], 1), Fr(1, 4))[0]) for a in range(NA) for b in range(NA))
    ric = ricci(R, eta)
    HHdiff = any(Tb[i][j] != ctx.c(ric[i][j]) for i in range(4) for j in range(4))
    check(f"base (1,3), F_jets brute force: base metric eta - (1/3) R_ikjl x^k x^l with R the pp-wave tensor above (base "
          f"curvature at 0 equals R: {okR}); Z2 metric (c_H, c_V, lambda) = (-2, 3, 1/2) at (0, y0), y0 = P eta P^T != eta: "
          f"Ric^Y on VV equals Ric_F = -G_1/4 exactly (A-term zero): {okVV}; Ric^Y on HH differs from Ric_g, i.e. "
          f"(A_X, A_X') != 0 and A != 0 there: {HHdiff}  {tk(t1)}", okR and okVV and HHdiff)
    kg = kg_of((1, 1, 1, 1))
    check(f"base (4,0): the trace form on so(4) in the basis E_cd is {kg.diagonal().tolist()[0]} (negative definite), so a "
          f"self-adjoint curvature operator with square 0 vanishes: b_R = 0 iff R = 0 there",
          kg.is_diagonal() and all(v < 0 for v in kg.diagonal()))


# ====================================================================== Y2
def y2sec():
    section("Y2. Lemma P2: rank structure off and at the section; generic inconsistency")
    rng = random.Random(20261243)
    for eta in SIGS:
        t1 = time.time()
        lam = lam_of(eta)
        res, ok = [], True
        for _ in range(2):
            R = rand_curv(rng)
            Ks, beta = R_beta(R, eta, Fr(rng.randint(1, 5), rng.randint(1, 5)))
            y0 = offpoint(rng, eta)
            rows, rhs, _, _, _ = system_at(y0, lam, Ks, beta)
            r = (rank(rows), rank(augment(rows, rhs)))
            res.append(r)
            ok &= r == (19, 20) and all(v == 0 for v in matvec(rows, coords(y0) + [Fr(0)] * 10))
        # at the section: random b in Sym^2 so(eta)
        geo = Geo([diag(eta)], lam)
        b = {}
        for p in range(6):
            for q in range(p, 6):
                v = Fr(rng.randint(-3, 3), rng.randint(1, 3))
                b[(p, q)] = v
                b[(q, p)] = v
        A, nA, dA = aterm(geo, so_basis(eta), b)
        okz = all(v == [0] for r in A for v in r) and all(v == [0] for C in nA for r in C for v in r)
        rows, rhs = relaxed_rows(geo, A, nA)
        r0 = (rank(order(rows, 0)), rank(augment(order(rows, 0), orderv(rhs, 0))))
        check(f"base {sl(eta)}: two random R (random c_A), two random off-section points: (rank A, rank [A|b]) = {res} "
              f"(dilation in the kernel; inconsistent); at the section y = eta with a random b in Sym^2 so(eta): A_b(eta) = 0, "
              f"nabla A_b(eta) = 0, ranks {r0}  {tk(t1)}", ok and okz and r0 == (19, 19))


# ====================================================================== Y3
def y3sec():
    section("Y3. Lemma P3: first order at the section along the ten coordinate lines; b in Sym^2 so(eta) free")
    for eta in SIGS:
        t1 = time.time()
        lam = lam_of(eta)
        Ks = so_basis(eta)
        sb = [(p, q) for p in range(6) for q in range(p, 6)]
        betas = [({(p, q): Fr(1), (q, p): Fr(1)} if p != q else {(p, p): Fr(1)}) for p, q in sb]
        kgi = kg_of(eta).inv()
        cas = sp.Matrix([[kgi[p, q] for p, q in sb]])
        geo0 = Geo([diag(eta)], lam)
        rows0, _ = relaxed_rows(geo0, None, None)
        A0 = order(rows0, 0)
        rA0 = rank(A0)
        b1, okO = [], True
        for j in range(NA):
            gj = Geo([diag(eta), EBS[j]], lam, curv=False)
            kd = kill_data(gj, Ks)
            bj = []
            for beta in betas:
                A, nA, _ = aterm(gj, Ks, beta, kd)
                okO &= all(A[a][b][0] == 0 and A[a][b][1] == 0 for a in range(NA) for b in range(NA))
                r = rhs_of(nA)
                okO &= all(v[0] == 0 for v in r)
                bj.append(orderv(r, 1))
            b1.append(bj)
        nunk = 21 + 20 * NA
        big = []
        for j in range(NA):
            for r in range(450):
                row = [Fr(0)] * nunk
                for m in range(21):
                    row[m] = -b1[j][m][r]
                for c in range(20):
                    row[21 + 20 * j + c] = A0[r][c]
                big.append(row)
        ns = dm(big, nunk).nullspace().to_Matrix()
        bp = ns[:, :21]
        rk = bp.rank()
        okc = rk == 1 and sp.Matrix.vstack(bp, cas).rank() == 1
        check(f"base {sl(eta)}: A_0 at eta has rank {rA0}; along eta + eps E_j (j = 1..10) the A-term is O(eps^2) and b(0) = 0 "
              f"for every b in Sym^2 so(eta): {okO}; the first-order conditions b_1(b, E_j) in Im A_0 for all j "
              f"({len(big)} rows, {nunk} unknowns: b and zeta_1(E_j)) have nullity {ns.shape[0]} = 10 + 1 and b-projection "
              f"of rank {rk}, spanned by Cas = kg^-1: {okc}  {tk(t1)}", rA0 == 19 and okO and ns.shape[0] == 11 and okc)


# ====================================================================== Y4
def curve_test(eta, Yser, lam):
    N = len(Yser)
    Ks = so_basis(eta)
    geo = Geo(Yser, lam)
    rows0, _ = relaxed_rows(geo, None, None)                 # A^0(eps): B = Ric_F
    A, nA, _ = aterm(geo, Ks, cas_beta(eta))
    okA = all(A[a][b][0] == 0 and A[a][b][1] == 0 for a in range(NA) for b in range(NA))
    rowsC = xi_rows(A, N)                                     # k A^C(eps): the Casimir A-term in the xi-columns
    bC = rhs_of(nA)                                           # k b^C(eps)
    okA &= all(v[0] == 0 for v in bC)
    A0, A1, A2, A2C = order(rows0, 0), order(rows0, 1), order(rows0, 2), order(rowsC, 2)
    b1, b2, b3 = (orderv(bC, n) for n in (1, 2, 3))
    dser = [coords(Yn) + [Fr(0)] * 10 for Yn in Yser]
    okd = True
    for n in range(N):
        v = [Fr(0)] * 450
        for i in range(n + 1):
            v = [a + b + c for a, b, c in zip(v, matvec(order(rows0, i), dser[n - i]), matvec(order(rowsC, i), dser[n - i]))]
        okd &= all(x == 0 for x in v)
    z1 = solve(A0, b1)
    z2 = solve(A0, [a - b for a, b in zip(b2, matvec(A1, z1))]) if z1 is not None else None
    u = [a - b - c for a, b, c in zip(b3, matvec(A1, z2), matvec(A2, z1))] if z2 is not None else None
    v = [-x for x in matvec(A2C, z1)] if z1 is not None else None
    ru = (rank(augment(A0, u)), rank(augment(A0, v)), rank(augment(A0, u, v))) if u is not None else None
    direct = {}
    for k in (Fr(0), Fr(1, 12), Fr(1), Fr(-3), Fr(7, 5)):
        big, rhs = [], []
        for n in range(N):
            for r in range(450):
                row = [Fr(0)] * (20 * N)
                for i in range(n + 1):
                    for c in range(20):
                        val = rows0[r][c][i] + k * rowsC[r][c][i]
                        if val:
                            row[20 * (n - i) + c] += val
                big.append(row)
                rhs.append(k * bC[r][n])
        direct[str(k)] = (rank(big), rank(augment(big, rhs)))
    return okA, okd, z1 is not None, z2 is not None, ru, direct


def y4sec():
    section("Y4. Lemma P4: the line b = k Cas; order 3 at the section (direct, two curves); 20-minors in Q[k]; R_L")
    rng = random.Random(20261244)
    for eta in SIGS:
        lam = lam_of(eta)
        h = [[Fr(rng.randint(-3, 3)) for _ in range(4)] for _ in range(4)]
        h = [[h[i][j] + h[j][i] for j in range(4)] for i in range(4)]
        X = [[Fr(rng.randint(-3, 3)) for _ in range(4)] for _ in range(4)]
        e = diag(eta)
        for lab, Yser in (('line eta + eps h', [e, h, mz(), mz()]),
                          ('orbit (1 + eps X) eta (1 + eps X)^T', [e, madd(mmul(X, e), mmul(e, mT(X))), mmul(mmul(X, e), mT(X)), mz()])):
            t1 = time.time()
            okA, okd, s1, s2, ru, direct = curve_test(eta, Yser, lam)
            okdir = direct['0'][0] == direct['0'][1] and all(direct[k][1] == direct[k][0] + 1 for k in direct if k != '0')
            check(f"base {sl(eta)}, {lab} (h, X random), b = k Cas, computed at y(eps) as power series to eps^3: A-term "
                  f"O(eps^2), b(0) = 0: {okA}; dilation y(eps) in the kernel to order 3: {okd}; order 1 solvable: {s1}; "
                  f"order 2 solvable: {s2}; with zeta_1 = k z1, zeta_2 = k z2 the order-3 equation is A_0 zeta_3 = k u + k^2 v, "
                  f"ranks [A_0|u], [A_0|v], [A_0|u|v] = {ru}, so k = 0; the full truncated Taylor system (orders 0..3, "
                  f"zeta_0..zeta_3 free) has (rank, augmented rank) {direct}: consistent only at k = 0  {tk(t1)}",
                  okA and okd and s1 and s2 and ru == (20, 19, 20) and okdir)
    kk = sp.Symbol('k')
    rng = random.Random(20261245)
    for eta in SIGS:
        t1 = time.time()
        lam = lam_of(eta)
        Ks = so_basis(eta)
        polys, fls, okall = [], [], True
        for _ in range(2):
            y0 = offpoint(rng, eta)
            geo = Geo([y0], lam)
            rows0, _ = relaxed_rows(geo, None, None)
            A, nA, _ = aterm(geo, Ks, cas_beta(eta))
            M0 = [[e_[0] for e_ in r] + [Fr(0)] for r in rows0]
            bC = rhs_of(nA)
            M1 = [[e_[0] for e_ in r] + [bC[i][0]] for i, r in enumerate(xi_rows(A, 1))]
            d = coords(y0)
            j0 = next(j for j in range(NA) if d[j] != 0)
            cols = [j for j in range(21) if j != j0]                 # 19 A-columns (dilation removed) + the rhs column
            k0 = Fr(5, 11)
            sel = pick_rows([[M0[i][j] + k0 * M1[i][j] for j in cols] for i in range(450)], 20)
            okall &= len(sel) == 20
            pts = [Fr(i + 1, 3) for i in range(13)]
            vals = [fdet([[M0[i][j] + kv * M1[i][j] for j in cols] for i in sel]) for kv in pts]
            poly = sp.Poly(sp.expand(sp.interpolate([(sp.Rational(p_.numerator, p_.denominator), sp.Rational(v_.numerator, v_.denominator))
                                                      for p_, v_ in zip(pts, vals)], kk)), kk)
            okall &= poly.degree() <= 11
            for kv in (Fr(-7, 2), Fr(13, 9)):                        # the interpolant is the minor (degree <= 11)
                okall &= fdet([[M0[i][j] + kv * M1[i][j] for j in cols] for i in sel]) == Fr(str(poly.eval(sp.Rational(kv.numerator, kv.denominator))))
            polys.append(poly)
            _, fl = sp.factor_list(poly.as_expr())
            fls.append(' * '.join(f"({sp.expand(f_)})^{m_}" if m_ > 1 else f"({sp.expand(f_)})" for f_, m_ in fl))
        g = sp.gcd(polys[0], polys[1])
        okg = sp.expand(sp.Poly(g, kk).monic().as_expr() - kk ** 2) == 0
        check(f"base {sl(eta)}, b = k Cas at two own off-section points: a 20-minor of [A(k) | b(k)] (19 A-columns and the rhs "
              f"column; 20 rows independent at k = 5/11) is, up to a constant, {fls[0]}  and  {fls[1]}; gcd = "
              f"{sp.factor(g.as_expr())}: for every k != 0 one of the two minors is nonzero, so k Cas is not in Z_k "
              f"(P2(d)(i))  {tk(t1)}", okall and okg)
    # M5(d)'s value: R_L at (c_V/c_H)^2 = 1/12
    t1 = time.time()
    rng = random.Random(20261246)
    RL = R_L()
    res, ok = [], True
    for eta in ((1, -1, -1, -1), (-1, 1, 1, 1)):
        Rm = curv_op(RL, eta)
        ric = ricci(RL, eta)
        ok &= bianchi(RL) and Rm * Rm == -sp.eye(6) and any(v for r in ric for v in r)
        ok &= Fr(1, 48) * BR16(RL, eta) == opCas16(sp.eye(6), eta) / 12           # b_(R_L) = Cas/12 at c_A = 1/48
        Ks, beta = R_beta(RL, eta, Fr(1, 48))
        y0 = offpoint(rng, eta)
        rows, rhs, _, _, _ = system_at(y0, lam_of(eta), Ks, beta)
        r = (rank(rows), rank(augment(rows, rhs)))
        res.append((sl(eta), r))
        ok &= r == (19, 20)
    check(f"R_L (R_(01,12) = R_(02,23) = R_(03,31) = 1) in bases (1,3), (3,1): Bianchi, (curvature operator)^2 = -Id, Ric != 0, "
          f"and at c_A = 1/48, i.e. (c_V/c_H)^2 = 1/12, b_(R_L) = Cas/12 (M5(d)'s point k = 1/12 of the line); relaxed system "
          f"at a random off-section point: {res} (inconsistent)  {tk(t1)}", ok)


def main():
    y0sec()
    y1sec()
    y2sec()
    y3sec()
    y4sec()
    P(f"\nSUMMARY: {NCHK[0]} checks, {NCHK[1]} True, {NCHK[0] - NCHK[1]} False")
    P(f"runtime {time.time() - T0:.0f} s")


if __name__ == '__main__':
    main()
