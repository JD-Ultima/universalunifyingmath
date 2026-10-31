"""GU round 15 (S): Upsilon = 0 in the Zorro metric Z1 on open sets of Y that need not meet s(X).

Formalization: formalization/geometric-unity/S-z1-off-section.md.  Builds on N (N1, N2, N4, N6, N7, N9-N11), Q (Q1-Q8, Q12,
Q15), P (P1, P2(d)) and F (F2.4, F4, F6).  Imports F_jets.py and N_z1_neighbourhood.py (and through it L_nonradial.py,
J_slaved_sector.py, ...) unmodified.  Exact arithmetic throughout: Python Fractions (also as coefficients of N's truncated
power series), sympy rationals.  Fixed seeds.  One core.  Every check is asserted.

Sections (labels = results of the note)
  S0  Lemma S1 (the Euler field off the section): nabla_E X^H = X^H/2, nabla_E d_A = -d_A at points off s(X) over a curved base
      in general coordinates; the fibre equation alone has solutions over the flat model for every lambda.
  S1  Lemma S2 / Corollary S3 (the parallel set Pi): generic curvature has only R g as invariant symmetric form; product and
      pp-wave bases have more, and at those points Ric(Z1) is F4(c) for the parallel metric through the point.
  S2  Lemma S4 (fibre polynomiality, dichotomy): det(y)^3 Ric(Z1) is polynomial along a fibre line; the degree-1 system (*)
      has rank 28 at section and radial points, is inconsistent at a generic point off the section over a curved base, and is
      consistent exactly as predicted over the flat base.
  S3  Theorem S5, Lemma S6 (Einstein pieces, F9): the order-2 VV fibre coefficient of Ric(Z1) is Z2's A-term for Weyl tensors
      in bases (1,3) and (2,2) (Schur data); the order-2 A-term determines b; null (type N and its (2,2) analogues) curvatures
      fail the Einstein condition off the section.
  S4  Propositions S7-S9 (linear mu): rank 28 of (*) at section and radial points (several bases, mu_1 = -7/18, -5/12);
      the algebra of the trace-equation continuation (lambda = 13/56) and of lambda = 1/5.
  S5  Lemma S13, Theorem S14 (nonlinear mu, the twisted cone): the relaxed system (42 unknowns): rank 42 at section and radial
      points, its two degeneracies, the Einstein solution xi = k v (k of Q6), inconsistency off the section; the primitive
      G(t) of H/x and the derivation of Q7's (D) for a complex branch.
Run (from computation/geometric-unity): python3 S_z1_off_section.py > S_z1_off_section.out 2>&1
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
import N_z1_neighbourhood as N  # noqa: E402

T0 = time.time()
NCHK = [0]
n = 4
IDX = FJ.sym_index(n)
FC = N.FC
T, Alg, FCtx, mon = N.T, N.Alg, N.FCtx, N.mon
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


def tk(t1):
    return f"  [{time.time() - t1:.0f} s]"


def sig_of(eta):
    return (eta.count(1), eta.count(-1))


def QQ(v):
    return sp.QQ(v.numerator, v.denominator)


def rank_of(rows):
    return DomainMatrix([[QQ(v) for v in r] for r in rows], (len(rows), len(rows[0])), sp.QQ).rank()


def rand_sym(rng, lo=-2, hi=2, den=5):
    M = [[Fr(0)] * n for _ in range(n)]
    for i in range(n):
        for j in range(i, n):
            v = Fr(rng.randint(lo, hi), den)
            M[i][j] = v
            M[j][i] = v
    return M


def pert(eta, rng):
    """diag(eta) + a symmetric perturbation with entries in {-1, 0, 1}/10 (Gershgorin: same signature, det != 0)"""
    return madd(diag(eta), rand_sym(rng, -1, 1, 10))


def madd(A, B, s=1):
    return [[A[i][j] + s * B[i][j] for j in range(n)] for i in range(n)]


def diag(v):
    return [[Fr(v[i]) if i == j else Fr(0) for j in range(n)] for i in range(n)]


def det4(M):
    return Fr(sp.Matrix(n, n, lambda i, j: sp.Rational(M[i][j])).det())


# ====================================================================== geometry of Z1 at a point (x0 = 0, y0 arbitrary)
def fibre_point_gen(ctx, y0, A=None, offset=0):
    Y = [[ctx.c(y0[i][j]) for j in range(n)] for i in range(n)]
    if A is not None:
        for a, (i, j) in enumerate(IDX):
            e = T(A, {A.idx[mon(offset + a, A.nv)]: Fr(1)})
            Y[i][j] = Y[i][j] + e
            if i != j:
                Y[j][i] = Y[j][i] + e
    return Y


def base_of(pb):
    gj = N.base_jets_at(pb, FC, [0] * n)
    return gj, FJ.base_geometry(gj, n, FC)


def ric_point(gj, base, y0, lam, cH, cV, variant='Z1'):
    G, m, idx, Psi0 = FJ.total_metric(gj, base, n, y0, lam, cH, cV, variant, FC)
    Ric, g0, g0i = FJ.ricci_at_point(G, m, FC)
    return Ric, g0, Psi0, G


def point_data_y(pb, y0, lam, cH, cV):
    """N.point_data at the point (0, y0) with y0 an arbitrary symmetric matrix (normal-type base: Gamma(0) = 0)"""
    m = 14
    gj, base = base_of(pb)
    G, _, _, Psi0 = FJ.total_metric(gj, base, n, fibre_point_gen(FC, y0), lam, cH, cV, 'Z1', FC)
    assert all(v == 0 for row in Psi0 for v in row)
    Rm, g0, g0i = FJ.riemann_at_point(G, m, FC)
    _, dG, _ = FJ.metric_point_data(G, m, FC)
    _, Gam = FJ.christoffel(g0i, dG, m, FC)
    A = Alg(14, 1)
    ctx = FCtx(A)
    xi = [T(A, {A.idx[mon(k, 14)]: Fr(1)}) for k in range(n)]
    gjs = N.base_jets_at(pb, ctx, xi)
    bs = FJ.base_geometry(gjs, n, ctx)
    Gs, _, _, _ = FJ.total_metric(gjs, bs, n, fibre_point_gen(ctx, y0, A, offset=n), ctx.c(lam), ctx.c(cH), ctx.c(cV),
                                  'Z1', ctx)
    Rics, _, _ = FJ.ricci_at_point(Gs, m, ctx)
    Ric = [[Rics[a][b].coef((0,) * 14) for b in range(m)] for a in range(m)]
    dRic = [[[Rics[a][b].coef(mon(c, 14)) for b in range(m)] for a in range(m)] for c in range(m)]
    nRic = [[[dRic[c][a][b] - sum(Gam[d][c][a] * Ric[d][b] + Gam[d][c][b] * Ric[a][d] for d in range(m))
              for b in range(m)] for a in range(m)] for c in range(m)]
    return dict(g0=g0, g0i=g0i, R=Rm, Ric=Ric, nRic=nRic, base=base, Gam=Gam)


def relaxed_system(D):
    """(*) with mu'(r) <v, X> replaced by a free covector xi(X): unknowns (v^c, xi_c, n^c), 42 in all"""
    m = 14
    g0, g0i, R, Ric, nRic = D['g0'], D['g0i'], D['R'], D['Ric'], D['nRic']
    Rh = [[sum(g0i[d][e] * Ric[e][b] for e in range(m)) for b in range(m)] for d in range(m)]
    nRh = [[[sum(g0i[d][e] * nRic[c][e][b] for e in range(m)) for b in range(m)] for d in range(m)] for c in range(m)]
    rows, rhs = [], []
    for a in range(m):
        for b in range(a + 1, m):
            for d in range(m):
                row = [Fr(0)] * (3 * m)
                for c in range(m):
                    row[c] += sum(g0i[d][e] * R[(e, c, a, b)] for e in range(m) if g0i[d][e] != 0)
                    row[m + c] -= ((a == c) * Rh[d][b] - (b == c) * Rh[d][a])
                    row[2 * m + c] -= (g0[a][c] * (d == b) - g0[b][c] * (d == a))
                rows.append(row)
                rhs.append(nRh[a][d][b] - nRh[b][d][a])
    return rows, rhs


def solve(rows, rhs, K):
    return FJ.solve_fit([[r[k] for r in rows] for k in range(K)], rhs, FC)


def E_of(y0):
    return [Fr(0)] * n + [y0[i][j] for (i, j) in IDX]


def lower(v, g0):
    return [sum(g0[a][b] * v[b] for b in range(14)) for a in range(14)]


# ====================================================================== curvature tensors
def two_form(l, mm, eta):
    lo = [Fr(eta[a]) * l[a] for a in range(n)]
    mo = [Fr(eta[a]) * mm[a] for a in range(n)]
    return [[lo[a] * mo[b] - lo[b] * mo[a] for b in range(n)] for a in range(n)]


def R_forms(terms):
    out = {k: Fr(0) for k in itertools.product(range(n), repeat=4)}
    for s, al in terms:
        for a, b, c, d in itertools.product(range(n), repeat=4):
            out[(a, b, c, d)] += s * al[a][b] * al[c][d]
    return out


def bianchi_ok(R):
    return all(R[(a, b, c, d)] + R[(a, c, d, b)] + R[(a, d, b, c)] == 0 for a, b, c, d in itertools.product(range(n), repeat=4))


PERM = {}
for _p in itertools.permutations(range(4)):
    _s = 1
    for _i in range(4):
        for _j in range(_i + 1, 4):
            if _p[_i] > _p[_j]:
                _s = -_s
    PERM[_p] = _s


def star_form(al, eta):
    return [[sum(Fr(PERM.get((a, b, c, d), 0), 2) * eta[c] * eta[d] * al[c][d] for c in range(n) for d in range(n))
             for b in range(n)] for a in range(n)]


def star_R(W, eta):
    return {(a, b, c, d): sum(Fr(PERM.get((a, b, e, f), 0), 2) * eta[e] * eta[f] * W[(e, f, c, d)]
                              for e in range(n) for f in range(n)) for a, b, c, d in itertools.product(range(n), repeat=4)}


def pair4(A, B, eta):
    return sum(Fr(eta[a] * eta[b] * eta[c] * eta[d]) * A[(a, b, c, d)] * B[(a, b, c, d)]
               for a, b, c, d in itertools.product(range(n), repeat=4))


def op_sq_zero(R, eta):
    """the curvature operator on Lambda^2, M[(ab),(ef)] = eta^e eta^f R_{abef}; returns M^2 == 0"""
    pairs = [(a, b) for a in range(n) for b in range(a + 1, n)]
    M = [[Fr(eta[e] * eta[f]) * R[(a, b, e, f)] for (e, f) in pairs] for (a, b) in pairs]
    M2 = [[sum(M[i][k] * M[k][j] for k in range(6)) for j in range(6)] for i in range(6)]
    return all(v == 0 for r in M2 for v in r), any(v != 0 for r in M for v in r)


def ric_zero(R, eta):
    return all(v == 0 for r in N.ric_of(R, eta) for v in r)


# ====================================================================== S0
def s0():
    section("S0. Lemma S1: the Euler field off the section (what replaces L6); the fibre equation alone over the flat model")
    m = 14
    for eta, lam, seed in (((1, 1, 1, 1), Fr(1, 7), 101), ((1, -1, -1, -1), Fr(1, 2), 102)):
        t1 = time.time()
        rng = random.Random(seed)
        gj = N.rand_general_base(rng, eta)
        base = FJ.base_geometry(gj, n, FC)
        y0 = pert(eta, rng)
        cH, cV = Fr(2), Fr(3)
        G, _, _, Psi0 = FJ.total_metric(gj, base, n, y0, lam, cH, cV, 'Z1', FC)
        g0, dG, _ = FJ.metric_point_data(G, m, FC)
        g0i = FJ.finv(g0, FC)
        _, Gam = FJ.christoffel(g0i, dG, m, FC)
        E = E_of(y0)
        ok_h = True
        for mu in range(n):
            X = [Fr(1) if a == mu else (Psi0[a - n][mu] if a >= n else Fr(0)) for a in range(m)]
            dEX = [Psi0[a - n][mu] if a >= n else Fr(0) for a in range(m)]       # E(Psi) = Psi (Psi linear in y)
            nab = [dEX[a] + sum(Gam[a][b][c] * E[b] * X[c] for b in range(m) for c in range(m)) for a in range(m)]
            ok_h &= nab == [x / 2 for x in X]
        ok_v = all([sum(Gam[a][b][n + A_] * E[b] for b in range(m)) for a in range(m)]
                   == [Fr(-1) if a == n + A_ else Fr(0) for a in range(m)] for A_ in range(10))
        nzPsi = any(v != 0 for r in Psi0 for v in r)
        check(f"base {sig_of(eta)} in general coordinates (Gamma(x0) != 0: {nzPsi}), y0 != g(x0), lambda = {lam}: nabla_E X^H = X^H/2"
              f" for the four horizontal lifts ({ok_h}) and nabla_E d_A = -d_A for the ten vertical coordinate fields ({ok_v}); hence"
              f" Hess f(E, X^H) = E(X^H f) - (X^H f)/2 and Hess f(E, d_A) = d_A(E f) for every f{tk(t1)}", nzPsi and ok_h and ok_v)
    l_, cv = sp.symbols('lambda c_V')
    c = 4 * cv * (1 - 4 * l_)
    b = (1 - 5 * l_) / (cv * (1 - 4 * l_))
    bS = sp.Rational(-5, 4) / cv
    s, nu = sp.symbols('s nu')
    f = sp.exp(sp.sqrt(MU1_s * b * c) * s)
    # g_F = g_{V0} + c ds^2 ; Ric_VV (flat model, N2) = beta_S g_F + b c ds^2 ; Hess^F f(s) = f'' ds^2
    eqV0 = sp.simplify(MU1_s * f * bS + nu)
    nu_s = sp.solve(eqV0, nu)[0]
    eqs = sp.simplify(sp.diff(f, s, 2) - (MU1_s * f * (bS * c + b * c) + nu_s * c))
    check(f"flat model, every lambda: w = exp(sqrt(mu_1 b c) s) solves the FIBRE equation Hess^F w = mu_1 w Ric_VV + nu g_F with"
          f" nu = -mu_1 beta_S w (V_0 block {sp.simplify(eqV0.subs(nu, nu_s))}, ds^2 block {eqs}); so no fibre-only obstruction"
          f" exists in Z1 (on Y, N11 allows only lambda = 13/56, 1/5)", sp.simplify(eqV0.subs(nu, nu_s)) == 0 and eqs == 0)


MU1_s = sp.Rational(-7, 18)


# ====================================================================== S1
def invariants_dim(R, eta, rng=None):
    """dimension of {y in Sym^2 : R(e_c, e_d) . y = 0 for all c < d} (R(X,Y) as endomorphisms, R^e_{b cd} = eta^e R_{ebcd})"""
    rows = []
    for c in range(n):
        for d in range(c + 1, n):
            Rm = [[Fr(eta[e]) * R[(e, b, c, d)] for b in range(n)] for e in range(n)]    # Rm[e][b] = R^e_b
            for i in range(n):
                for j in range(i, n):
                    row = []
                    for (p, q) in IDX:
                        Y = [[Fr(0)] * n for _ in range(n)]
                        Y[p][q] = Fr(1)
                        Y[q][p] = Fr(1)
                        # (R . y)_{ij} = -(R^T y + y R)_{ij} = -(sum_e R^e_i y_ej + y_ie R^e_j)
                        row.append(-sum(Rm[e][i] * Y[e][j] + Y[i][e] * Rm[e][j] for e in range(n)))
                    rows.append(row)
    return 10 - rank_of(rows)


def f4c_check(gj, base, y0, lam, cH, cV):
    """Ric(Z1) at (0, y0) against F4(c) for a metric g' parallel for the base connection with g'(0) = y0 (Gamma(0) = 0)"""
    Ric, g0, Psi0, _ = ric_point(gj, base, y0, lam, cH, cV)
    assert all(v == 0 for r in Psi0 for v in r)
    c = 4 * cV * (1 - 4 * lam)
    G15 = N.dewitt_at(FC, y0, Fr(1, 5))
    okH = all(Ric[i][j] == base['Ric'][i][j] - cH / c * y0[i][j] for i in range(n) for j in range(n))
    okHV = all(Ric[i][n + A_] == 0 for i in range(n) for A_ in range(10))
    okV = all(Ric[n + A_][n + B] == Fr(-5, 4) * G15[A_][B] for A_ in range(10) for B in range(10))
    return okH, okHV, okV, Ric, g0


def s1():
    section("S1. Lemma S2, Corollary S3: the parallel set Pi (points on sections of metrics with the same Levi-Civita connection)")
    rng = random.Random(111)
    dims = {}
    for eta in ((1, 1, 1, 1), (1, -1, -1, -1), (1, 1, -1, -1)):
        dims[sig_of(eta)] = invariants_dim(N.rand_R(rng, eta), eta)
    check(f"generic curvature: the symmetric forms annihilated by every R(X,Y) have dimension {dims} (= R g(x0)), so Pi_x is"
          f" the radial set Sigma_x (and -Sigma_x in base (2,2)) for every base whose curvature at x spans so(g)",
          all(v == 1 for v in dims.values()))
    # product base: g = g1(x0,x1) + g2(x2,x3) (exactly, from normal_base), y0 = g1(0) + a g2(0)
    eta = (1, 1, 1, 1)
    Rp = {k: Fr(0) for k in itertools.product(range(n), repeat=4)}
    for (a, b, K) in ((0, 1, Fr(1)), (2, 3, Fr(-2))):
        for (p, q, r, s, sg) in ((a, b, a, b, 1), (b, a, b, a, 1), (a, b, b, a, -1), (b, a, a, b, -1)):
            Rp[(p, q, r, s)] = sg * K
    pb = N.normal_base(eta, Rp)
    gj, base = base_of(pb)
    prod = all(M[i][j] == 0 or ((i < 2) == (j < 2) == (key[0] < 2) == (key[1] < 2))
               for key, M in pb[1].items() for i in range(n) for j in range(n)) and not pb[2]
    lam, cH, cV = Fr(1, 7), Fr(2), Fr(3)
    y0 = diag((1, 1, 3, 3))
    okH, okHV, okV, _, _ = f4c_check(gj, base, y0, lam, cH, cV)
    y1 = madd(y0, [[Fr(1, 5) if {i, j} == {0, 2} else Fr(0) for j in range(n)] for i in range(n)])
    cH1, cHV1, cV1, _, _ = f4c_check(gj, base, y1, lam, cH, cV)
    check(f"product base (curvatures 1 and -2 on the factors; exact product: {prod}): invariant forms of dimension"
          f" {invariants_dim(Rp, eta)}; at y0 = g1(0) + 3 g2(0) (on the section of the parallel metric g1 + 3 g2), Ric(Z1) is F4(c)"
          f" for that metric: HH ({okH}), HV = 0 ({okHV}), VV = -(5/4) G_1/5 ({okV}); control y0 + (1/5)(E_02 + E_20) (not"
          f" invariant): ({cH1}, {cHV1}, {cV1})",
          prod and invariants_dim(Rp, eta) == 2 and okH and okHV and okV and not (cH1 and cHV1 and cV1))
    # pp-wave (plane wave) base (1,3): g = -(2 du dv + (x^2 - y^2 + x y) du^2 + dx^2 + dy^2); l = du is parallel
    z = [[Fr(0)] * n for _ in range(n)]
    g0 = [row[:] for row in z]
    g0[0][1] = g0[1][0] = Fr(-1)
    g0[2][2] = g0[3][3] = Fr(-1)
    quad = {}
    for key, v in (((2, 2), Fr(-1)), ((3, 3), Fr(1)), ((2, 3), Fr(-1))):
        M = [row[:] for row in z]
        M[0][0] = v
        quad[key] = M
    gjw = FJ.base_from_coeffs(g0, [[row[:] for row in z] for _ in range(n)], quad, {}, n, FC)
    basew = FJ.base_geometry(gjw, n, FC)
    Rw = basew['R']
    # invariants of the curvature endomorphisms R^e_{bcd} = g^{ea} R_{abcd} (g(0) not diagonal: use g0i)
    g0i = basew['g0i']
    rows = []
    for c in range(n):
        for d in range(c + 1, n):
            Rm = [[sum(g0i[e][a] * Rw[(a, b, c, d)] for a in range(n)) for b in range(n)] for e in range(n)]
            for i in range(n):
                for j in range(i, n):
                    row = []
                    for (p, q) in IDX:
                        Y = [[Fr(0)] * n for _ in range(n)]
                        Y[p][q] = Fr(1)
                        Y[q][p] = Fr(1)
                        row.append(-sum(Rm[e][i] * Y[e][j] + Y[i][e] * Rm[e][j] for e in range(n)))
                    rows.append(row)
    dimw = 10 - rank_of(rows)
    ricflat = all(v == 0 for r in basew['Ric'] for v in r)
    out = {}
    for lam in (Fr(1, 2), Fr(1, 5)):
        for a in (Fr(2), Fr(-3, 2)):
            y0 = [row[:] for row in g0]
            y0[0][0] += a
            okH, okHV, okV, Ric, gg = f4c_check(gjw, basew, y0, lam, cH, cV)
            c = 4 * cV * (1 - 4 * lam)
            ein = all(Ric[p][q] == -gg[p][q] / c for p in range(14) for q in range(14))
            out[(str(lam), str(a))] = (okH and okHV and okV, ein)
    ok = all(v[0] for v in out.values()) and all(v[1] == (k[0] == '1/5') for k, v in out.items())
    check(f"plane-wave base (1,3), Ricci-flat ({ricflat}), l = du parallel: invariant forms of dimension {dimw} (g and l (x) l); at"
          f" y0 = g(0) + a l (x) l (a = 2, -3/2) Ric(Z1) is F4(c) for the parallel metric g + a l (x) l, and Z1 is Einstein there"
          f" exactly at lambda = 1/5 (F6 transferred): {out}", ricflat and dimw == 2 and ok)
    # Z1 depends on g only through aleph: base g and 3 g give the same Ricci at the same point of Y
    rng = random.Random(112)
    eta = (1, -1, -1, -1)
    R = N.rand_R(rng, eta)
    gjA, bA = base_of(N.normal_base(eta, R))
    gjB, bB = base_of(N.normal_base(eta, R, None, 3))
    y0 = pert(eta, rng)
    RA, _, _, GA = ric_point(gjA, bA, y0, Fr(1, 2), cH, cV)
    RB, _, _, GB = ric_point(gjB, bB, y0, Fr(1, 2), cH, cV)
    check("F2.4: the bases g and 3g (same Levi-Civita connection) give the same Z1 2-jet and Ricci tensor at the same point of Y",
          RA == RB and all(FJ.jval(GA[a][b], FC) == FJ.jval(GB[a][b], FC) for a in range(14) for b in range(14)))


# ====================================================================== S2
def s2():
    section("S2. Lemma S4: Ric(Z1) along a fibre is polynomial in (y, 1/det y); the dichotomy for the degree-1 system (*)")
    t1 = time.time()
    rng = random.Random(121)
    eta = (1, 1, 1, 1)
    R = N.rand_R(rng, eta)
    nb = N.normal_base(eta, R)
    pb = (nb[0], nb[1], N.rand_cubic(rng))
    gj, base = base_of(pb)
    lam, cH, cV = Fr(1, 7), Fr(2), Fr(3)
    y0 = diag(eta)
    h = rand_sym(rng, -1, 1, 10)
    ts = [Fr(k, 7) for k in range(-6, 14)]
    vals = []
    for t in ts:
        y = madd(y0, h, t)
        Ric, _, _, _ = ric_point(gj, base, y, lam, cH, cV)
        d3 = det4(y) ** 3
        vals.append([[d3 * Ric[a][b] for b in range(14)] for a in range(14)])

    def mindeg(seq):
        dd = list(seq)
        for k in range(len(ts)):
            if all(v == 0 for v in dd):
                return k - 1
            dd = [(dd[i + 1] - dd[i]) / (ts[i + k + 1] - ts[i]) for i in range(len(dd) - 1)]
        return None
    degs = [mindeg([vals[k][a][b] for k in range(len(ts))]) for a in range(14) for b in range(a, 14)]
    mx = max(d for d in degs if d is not None) if all(d is not None for d in degs) else None
    check(f"base (4,0) with random R and cubic terms, lambda = 1/7, the line y = g(x0) + t h of the fibre (20 exact nodes): every"
          f" component of det(y)^3 Ric(Z1) is a polynomial in t, of degree <= {mx} (<= 13 confirmed at >= 6 extra nodes){tk(t1)}",
          mx is not None and mx <= 13)
    # the degree-1 system (*) at section, radial and off-section points of a curved base
    rng = random.Random(122)
    eta = (1, -1, -1, -1)
    R = N.rand_R(rng, eta)
    lam = Fr(1, 2)
    res = {}
    for tag, pbx, y in (('section', N.normal_base(eta, R), diag(eta)), ('(x0, 2 g(x0)) in Sigma', N.normal_base(eta, R, None, 2),
                                                                       diag([2 * e for e in eta])),
                        ('off Pi, y0 = diag(2,-1,-3,-5)', N.normal_base(eta, R), diag((2, -1, -3, -5)))):
        t1 = time.time()
        D = point_data_y(pbx, y, lam, cH, cV)
        rows, rhs = N.star_system(D, MU1)
        res[tag] = N.rank_test(rows, rhs)
    check(f"curved non-Einstein base (1,3), lambda = 1/2, mu_1 = -7/18: rank (L, [L | Omega]) of (*) = {res}: full rank 28 at the"
          f" section and radial points (N7's homogeneous argument), inconsistent everywhere sampled (case (i) of the dichotomy)",
          all(v[0] == 28 and v[1] == 29 for v in res.values()))
    # flat base: the radial solutions r = e^{s/2} (13/56) and r = const (1/5) give consistency off the section; 1/7 does not
    eta = (1, 1, 1, 1)
    flat = (eta, {}, {})
    y = diag((2, 1, 3, 5))
    out = {}
    oks = True
    for lam in (Fr(13, 56), Fr(1, 5), Fr(1, 7)):
        D = point_data_y(flat, y, lam, cH, cV)
        rows, rhs = N.star_system(D, MU1)
        out[str(lam)] = N.rank_test(rows, rhs)
        if lam != Fr(1, 7):
            x, ok, rk = solve(rows, rhs, 28)
            c, b, bS = N.consts(lam, cV)
            E = E_of(y)
            vp = [e / (2 * MU1 * c) for e in E] if lam == Fr(13, 56) else [Fr(0)] * 14
            oks &= ok and rk == 28 and x[:14] == vp and x[14:] == [-MU1 * bS * v for v in vp]
    check(f"flat base (4,0), point y0 = diag(2,1,3,5): (*) has ranks {out}; at 13/56 the unique solution is v' = E/(2 mu_1 c),"
          f" n' = -mu_1 beta_S v' (the radial solution r = e^(s/2) of N11) and at 1/5 it is 0 (r constant) ({oks}); 1/7 inconsistent",
          oks and out['13/56'] == (28, 28) and out['1/5'] == (28, 28) and out['1/7'] == (28, 29))


# ====================================================================== S3
def weyl_pm_13(W, eta):
    """real and imaginary parts of the complex invariant I+ ~ |W|^2 - i <W, *W> in a Lorentzian base"""
    return pair4(W, W, eta), pair4(W, star_R(W, eta), eta)


def s3():
    section("S3. Theorem S5, Lemma S6: Einstein pieces and F9 for Z1 (n = 4) on every open set")
    lam, cH, cV = Fr(1, 5), Fr(2), Fr(3)
    for eta, seeds in (((1, -1, -1, -1), (131, 132)), ((1, 1, -1, -1), (133, 134))):
        t1 = time.time()
        Ws = [N.weyl_of(N.rand_R(random.Random(sd), eta), eta) for sd in seeds]
        if eta == (1, -1, -1, -1):
            inv = [weyl_pm_13(W, eta) for W in Ws]
            gen = all(any(v != 0 for v in W.values()) for W in Ws)
        else:
            pm = [N.weyl_pm(W, eta) for W in Ws]
            inv = [(N.norm2(p[0], eta), N.norm2(p[1], eta)) for p in pm]
            gen = all(any(v != 0 for v in p[0].values()) and any(v != 0 for v in p[1].values()) for p in pm)
        indep = inv[0][0] * inv[1][1] - inv[0][1] * inv[1][0] != 0
        nonnull = all(not op_sq_zero(W, eta)[0] for W in Ws)
        res = []
        for W in (Ws[0], {k: 2 * v for k, v in Ws[0].items()}, Ws[1]):
            pb = N.normal_base(eta, W)
            S1 = N.fibre_series(pb, lam, cH, cV, 2, 'Z1')
            S2 = N.fibre_series(pb, lam, cH, cV, 2, 'Z2')
            ctx = S1['ctx']
            G15, G14 = N.dewitt_at(ctx, S1['y0'], Fr(1, 5)), N.dewitt_at(ctx, S1['y0'], Fr(1, 4))
            agree = all(not (S1['Ric'][n + a][n + b] - S2['Ric'][n + a][n + b] - (Fr(-5, 4) * G15[a][b] + G14[a][b])).part(o)
                        for a in range(10) for b in range(10) for o in range(3))
            nzA = sum(len((S2['Ric'][n + a][n + b] + G14[a][b]).part(2)) for a in range(10) for b in range(10))
            res.append((agree, nzA > 0))
        check(f"Lemma S6, base {sig_of(eta)}, lambda = 1/5: for W1, 2 W1, W2 (random Weyl, nabla R = 0) the VV block of Ric(Z1) along"
              f" the fibre equals -(5/4) G_1/5 + A_x (Z2's O'Neill term) to order 2 ({[r[0] for r in res]}); the order-2 A-terms are"
              f" nonzero ({[r[1] for r in res]}); Schur data: W1, W2 have nonzero self-dual and anti-self-dual parts ({gen}), are not"
              f" null ({nonnull}), and their invariant pairs {[tuple(str(x) for x in i) for i in inv]} are independent ({indep})"
              f"{tk(t1)}", all(r[0] and r[1] for r in res) and gen and indep and nonnull)
    # the order-2 A-term determines b: K -> G_eta(., xi_K(h)) is injective on so(eta) for generic h
    rng = random.Random(135)
    ranks = {}
    for eta in ((1, 1, 1, 1), (1, -1, -1, -1), (1, 1, -1, -1), (1, 1, 1, -1)):
        h = rand_sym(rng, -3, 3, 1)
        Et = sp.diag(*eta)
        Ks = []
        for a in range(n):
            for b in range(a + 1, n):
                K = sp.zeros(n, n)
                K[a, b] = eta[a]
                K[b, a] = -eta[b]
                Ks.append(K)
        H = sp.Matrix(n, n, lambda i, j: sp.Rational(h[i][j]))
        rows = []
        for K in Ks:
            xi = -(K.T * H + H * K)
            rows.append([((Et * Em).T * Et * xi).trace() for Em in
                         [sp.Matrix(n, n, lambda i, j: 1 if (i, j) in ((p, q), (q, p)) else 0) for (p, q) in IDX]])
        ranks[sig_of(eta)] = sp.Matrix(rows).rank()
    check(f"for a random h the six covectors G_eta(., xi_K(h)), K in so(eta), are independent (ranks {ranks}); so the order-2"
          f" coefficient of A_b vanishes only for b = 0, i.e. (P1) only for R^2 = 0", all(v == 6 for v in ranks.values()))
    # null curvatures: type N in (1,3) and the eight (2,2) orbits fail the Einstein condition at a point of the fibre
    rng = random.Random(136)
    eta = (1, -1, -1, -1)
    al = two_form([1, 1, 0, 0], [0, 0, 1, 0], eta)
    be = two_form([1, 1, 0, 0], [0, 0, 0, 1], eta)
    RN = R_forms([(Fr(1), al), (Fr(-1), be)])
    sq, nz = op_sq_zero(RN, eta)
    prop = bianchi_ok(RN) and ric_zero(RN, eta) and sq and nz
    pb = N.normal_base(eta, RN)
    gj, base = base_of(pb)
    y0 = pert(eta, rng)
    defs = {}
    for (ch, cv) in ((Fr(2), Fr(3)), (Fr(2), Fr(-3))):
        c = 4 * cv * (1 - 4 * lam)
        Ric, gg, _, _ = ric_point(gj, base, y0, lam, ch, cv)
        defs[(str(ch), str(cv))] = sum(1 for a in range(n, 14) for b in range(n, 14) if Ric[a][b] + gg[a][b] / c != 0)
        Ric0, gg0, _, _ = ric_point(gj, base, diag(eta), lam, ch, cv)
        assert all(Ric0[a][b] + gg0[a][b] / c == 0 for a in range(14) for b in range(14))
    check(f"type N in (1,3): R = (l^e2)(x)(l^e2) - (l^e3)(x)(l^e3) is a Ricci-flat curvature tensor with R^2 = 0, R != 0 ({prop});"
          f" Z1 at lambda = 1/5 is Einstein at the section point but not at a point of the fibre: nonzero VV components of"
          f" Ric + g/c for (cH, cV) = (2, 3), (2, -3): {defs}", prop and all(v > 0 for v in defs.values()))
    eta = (1, 1, -1, -1)
    ap = two_form([1, 0, 1, 0], [0, 1, 0, 1], eta)
    am = two_form([1, 0, 1, 0], [0, 1, 0, -1], eta)
    sp_ = star_form(ap, eta)
    sm_ = star_form(am, eta)
    duality = (sp_ == ap and sm_ == [[-v for v in r] for r in am]) or (sp_ == [[-v for v in r] for r in ap] and sm_ == am)
    y0 = pert(eta, rng)
    defs = {}
    props = True
    for sP, sM in itertools.product((-1, 0, 1), repeat=2):
        if sP == 0 and sM == 0:
            continue
        Rn = R_forms([(Fr(sP), ap), (Fr(sM), am)])
        sq, nz = op_sq_zero(Rn, eta)
        props &= bianchi_ok(Rn) and ric_zero(Rn, eta) and sq and nz
        gj, base = base_of(N.normal_base(eta, Rn))
        for (ch, cv) in ((Fr(2), Fr(3)), (Fr(2), Fr(-3))):
            c = 4 * cv * (1 - 4 * lam)
            Ric, gg, _, _ = ric_point(gj, base, y0, lam, ch, cv)
            defs[(sP, sM, str(cv))] = sum(1 for a in range(n, 14) for b in range(n, 14) if Ric[a][b] + gg[a][b] / c != 0)
    check(f"base (2,2): l^m+ and l^m- (l, m+- null, totally null planes of the two families) are self-dual and anti-self-dual"
          f" ({duality}); the eight sign classes s+ a+ (x) a+ + s- a- (x) a- are Ricci-flat curvature tensors with R^2 = 0 ({props});"
          f" none is Einstein at the fibre point (nonzero VV components for cV = +-3: {sorted(set(defs.values()))})",
          duality and props and all(v > 0 for v in defs.values()))


# ====================================================================== S4
def s4():
    section("S4. Propositions S7-S9 (linear mu): rank of (*) at section and radial points; the continued trace equation")
    cH, cV = Fr(2), Fr(3)
    out = {}
    rng = random.Random(141)
    cases = []
    R40 = N.rand_R(rng, (1, 1, 1, 1))
    W40 = N.weyl_of(N.rand_R(rng, (1, 1, 1, 1)), (1, 1, 1, 1))
    R22 = N.rand_R(rng, (1, 1, -1, -1))
    cases.append(('(4,0) random, 1/7, scale 1/2', (1, 1, 1, 1), R40, Fr(1, 7), Fr(1, 2), MU1))
    cases.append(('(4,0) Weyl, 13/56, scale 2', (1, 1, 1, 1), W40, Fr(13, 56), Fr(2), MU1))
    cases.append(('(4,0) Weyl, 1/5, scale 3', (1, 1, 1, 1), W40, Fr(1, 5), Fr(3), MU1))
    cases.append(('(2,2) random, 2/3, scale 2, mu_1 = -5/12', (1, 1, -1, -1), R22, Fr(2, 3), Fr(2), Fr(-5, 12)))
    t1 = time.time()
    for tag, eta, R, lam, sc, mu1 in cases:
        D = point_data_y(N.normal_base(eta, R, None, sc), diag([sc * e for e in eta]), lam, cH, cV)
        rows, _ = N.star_system(D, mu1)
        out[tag] = rank_of(rows)
    check(f"rank of the coefficient matrix of (*) at radial points (x0, s g(x0)): {out} (all 28: the homogeneous system has only the"
          f" trivial solution, so consistency there is equivalent to the vanishing of the 29-minors){tk(t1)}",
          all(v == 28 for v in out.values()))
    # the trace-equation continuation at lambda = 13/56
    k_, N_, r = sp.symbols('kappa N r')
    Z = -1 / (224 * N_ * r)
    sig = (Z - k_) / (224 * N_)
    Phi = sp.Rational(7, 18) * r * (k_ * sig + 312 * N_ * sig ** 2)
    P0 = sp.expand(Phi / r)
    A = sp.expand(P0 * r ** 2).subs(r, 0)
    cv = sp.Symbol('c_V', nonzero=True)
    lam = sp.Rational(13, 56)
    c = 4 * cv * (1 - 4 * lam)
    b = (1 - 5 * lam) / (cv * (1 - 4 * lam))
    m1 = sp.Rational(-7, 18)
    dlog = sp.simplify(2 * m1 * b * c)
    nup = sp.simplify((1 + 4 * m1 * b ** 2 * c ** 2) / c)
    nuN8 = sp.simplify((sp.Rational(1, 4) + m1) / (m1 * c))
    sigma1 = -sp.Rational(23, 2) - (1 + lam) / (1 - 4 * lam)
    rho = 13 * (sp.Rational(7, 468) * sigma1 / cv + sp.Rational(3, 28) * m1 ** 2 * 4 * b ** 2 * c - m1 * nup)
    cst = sp.simplify(sp.expand(r ** 2 * (P0 - rho)).subs(r, 0))
    check(f"Prop. S8 (lambda = 13/56, Ricci-flat base): along Sigma the continued data are v' = 2bE, d log r/d theta = mu_1 v'(E)"
          f" = {dlog}, nu' = (1 + 4 mu_1 b^2 c^2)/c (= N8's value (1/4 + mu_1)/(mu_1 c): {sp.simplify(nup - nuN8) == 0}) and"
          f" rho = 13((7/468) s_Y + (3/28) mu_1^2 |v'|^2 - mu_1 nu') constant; P0(r) = Phi(r)/r has r^-2 coefficient A = {A}, so"
          f" r^2 (P0(r) - rho) has constant term {cst} != 0 and P0(r) = rho cannot hold on a curve r = r* e^(theta/2)",
          dlog == sp.Rational(1, 2) and sp.simplify(nup - nuN8) == 0 and cst != 0 and sp.simplify(cst - A) == 0
          and sp.simplify(A - sp.Rational(7 * 312, 18) / (224 ** 4 * N_ ** 3)) == 0)
    # lambda = 1/5: at p, v' = 0, nu'(p) = 1/c, and the fibre Hessian of a local solution vanishes at p (Einstein point)
    rng = random.Random(142)
    W = N.weyl_of(N.rand_R(rng, (1, 1, 1, 1)), (1, 1, 1, 1))
    D = point_data_y(N.normal_base((1, 1, 1, 1), W), diag((1, 1, 1, 1)), Fr(1, 5), cH, cV)
    rows, rhs = N.star_system(D, MU1)
    x, ok, rk = solve(rows, rhs, 28)
    c5 = 4 * cV * (1 - 4 * Fr(1, 5))
    ein = all(D['Ric'][a][b] == -D['g0'][a][b] / c5 for a in range(14) for b in range(14))
    zsol = ok and rk == 28 and all(v == 0 for v in x)
    check(f"Prop. S9 (lambda = 1/5, Weyl base (4,0)): (*) at p is consistent with the unique solution v' = n' = 0 ({zsol}),"
          f" and Ric(p) = -g/c ({ein}), so with nu'(p) = 1/c the fibre Hessian mu_1 r (Ric + nu' g)"
          f" of a local solution of d log r = mu_1 v' vanishes at p (N9 step 1)", ok and rk == 28 and all(v == 0 for v in x) and ein)


# ====================================================================== S5
def s5():
    section("S5. Lemma S13, Theorem S14 (nonlinear mu, the twisted cone): the relaxed system with a free covector xi")
    cH, cV = Fr(2), Fr(3)
    rng = random.Random(151)
    out = {}
    R13 = N.rand_R(rng, (1, -1, -1, -1))
    R40 = N.rand_R(rng, (1, 1, 1, 1))
    W40 = N.weyl_of(N.rand_R(rng, (1, 1, 1, 1)), (1, 1, 1, 1))
    eta13 = (1, -1, -1, -1)
    c12, b12, _ = N.consts(Fr(1, 2), cV)
    Kspec = -cH * b12 / 3                               # R = K R1 has Ric_g = 3K g = -cH b g
    g4 = N.gmat((1, 1, 1, 1))
    om = [Fr(1), Fr(2), Fr(0), Fr(-1)]
    w2 = sum(o * o for o in om)
    Rk1 = N.kn([[(om[i] * om[j] - w2 * g4[i][j] / 6) / 2 for j in range(n)] for i in range(n)], g4)   # Ric_g = om (x) om
    rk1 = sp.Matrix(N.ric_of(Rk1, (1, 1, 1, 1))).rank()
    t1 = time.time()
    for tag, eta, R, lam, sc in (('(1,3) random, 1/2', eta13, R13, Fr(1, 2), 1), ('(1,3) random, 1/2, scale 2', eta13, R13, Fr(1, 2), 2),
                                 ('(4,0) random, 1/7', (1, 1, 1, 1), R40, Fr(1, 7), 1),
                                 ('(4,0) random, 1/5', (1, 1, 1, 1), R40, Fr(1, 5), 1),
                                 ('(1,3) Ric = -cH b g, 1/2', eta13, N.radd((Kspec, N.R1(eta13))), Fr(1, 2), 1),
                                 ('(4,0) rank-1 Ricci, 1/5', (1, 1, 1, 1), Rk1, Fr(1, 5), 1),
                                 ('(4,0) Weyl, 1/5', (1, 1, 1, 1), W40, Fr(1, 5), 1)):
        D = point_data_y(N.normal_base(eta, R, None, sc), diag([sc * e for e in eta]), lam, cH, cV)
        rows, rhs = relaxed_system(D)
        out[tag] = N.rank_test(rows, rhs)
    gen = ['(1,3) random, 1/2', '(1,3) random, 1/2, scale 2', '(4,0) random, 1/7', '(4,0) random, 1/5']
    check(f"Lemma S13: ranks (L, [L|Omega]) of the relaxed system (42 unknowns) at section and radial points: {out}. Full rank 42"
          f" for the random (non-Einstein) bases, where it is inconsistent (relaxed consistency at Sigma forces an Einstein base);"
          f" the predicted degeneracies: rank 41 at Ric_g = -cH b g (b != 0) and at lambda = 1/5 with Ric_g of rank {rk1}, and a"
          f" larger kernel at lambda = 1/5 with Ric_g = 0{tk(t1)}",
          all(out[k] == (42, 43) for k in gen) and out['(1,3) Ric = -cH b g, 1/2'][0] == 41 and rk1 == 1
          and out['(4,0) rank-1 Ricci, 1/5'][0] == 41 and out['(4,0) Weyl, 1/5'][0] < 41)
    res = {}
    for tag, eta, lam, K in (('(1,3), 1/2, K = 1/2', eta13, Fr(1, 2), Fr(1, 2)), ('(4,0), 1/7, K = 1/3', (1, 1, 1, 1), Fr(1, 7), Fr(1, 3))):
        t1 = time.time()
        c, b, bS = N.consts(lam, cV)
        D = point_data_y(N.normal_base(eta, N.radd((K, N.R1(eta)))), diag(eta), lam, cH, cV)
        rows, rhs = relaxed_system(D)
        x, ok, rk = solve(rows, rhs, 42)
        sigma = 12 * K / cH
        k = -4 * b / (c * (sigma ** 2 - 16 * b ** 2))
        E = E_of(diag(eta))
        vE = 2 * b - sigma / 2
        Ef = lower(E, D['g0'])
        good = ok and rk == 42 and x[:14] == [vE * e for e in E] and x[14:28] == [k * vE * e for e in Ef] \
            and x[28:] == [-k * bS * vE * e for e in E]
        res[tag] = (ok, str(k), good)
    check(f"constant-curvature (Einstein) bases at the section: the relaxed system is consistent, and its unique solution is Q6's,"
          f" v' = (2b - sigma/2) E, xi = k v'^flat, n' = -k beta_S v' with k c sigma^2 = 4b(4kbc - 1): {res}",
          all(v[2] for v in res.values()))
    t1 = time.time()
    D = point_data_y(N.normal_base(eta13, R13), diag((2, -1, -3, -5)), Fr(1, 2), cH, cV)
    rows, rhs = relaxed_system(D)
    rr = N.rank_test(rows, rhs)
    check(f"relaxed system at an off-section point over the curved (1,3) base, lambda = 1/2: rank (L, [L|Omega]) = {rr}"
          f" (generic rank 42; inconsistent: case (i) of the dichotomy){tk(t1)}", rr == (42, 43))
    # the primitive of H/x and the quadratic kappa(t) = k
    tau, A, eps, kk = sp.symbols('tau A epsilon k')
    D1 = (1 + eps * A) + 544 * (1 - eps * A) * tau
    Q2 = 1 + 768 * tau + (121856 - 147456 * eps * A) * tau ** 2
    N2 = (7 - 3 * eps * A) + (7616 - 3008 * eps * A - 768 * A ** 2) * tau + (2071552 - 3226624 * eps * A + 1155072 * A ** 2) * tau ** 2
    res_sum = sp.limit(576 * D1 / Q2 * tau, tau, sp.oo)
    quad = sp.Poly(sp.expand(18 * kk * D1 ** 2 + N2), tau)
    ok_res = sp.simplify(res_sum - 576 * 544 * (1 - eps * A) / (121856 - 147456 * eps * A)) == 0
    check(f"dG/dtau = 576 D1/Q2 (G = int H dt / x): the sum of its residues is {sp.factor(res_sum)}, nonzero unless eps A = 1, so G"
          f" has logarithmic singularities at the zeros of Q2 (det M = 0) and at infinity (eps A != 1); kappa(tau) = k is the"
          f" quadratic 18 k D1^2 + N2 = 0 in tau (degree {quad.degree()})", ok_res and quad.degree() == 2)
    B, x_, kv, kp, Hv, rv = sp.symbols('B x kappa kappa_p H rho', nonzero=True)
    tp = x_ * (B / 2) * (1 - rv) / Hv                          # (B): (H/x) dt/dtheta = c v_E = (B/2)(1 - rho)
    lhs = B * kp * tp * (1 - rv ** 2) + B * kv * (-2 * rv * (-rv))     # d/dtheta [B kappa(t)(1 - rho^2)], rho' = -rho
    Phi = (x_ * kp / Hv) / kv ** 2                              # mu mu''/mu'^2 with mu mu'' = x kappa'(t)/H(t)
    target = (B * kv / 2) * (Phi * (1 - rv) + 4 * rv ** 2)
    Bsol = 1 / (kv * (1 - rv ** 2))                             # (C'): B kappa (1 - rho^2) = 1
    diffe = sp.simplify((lhs - target).subs(B, Bsol))
    check(f"Q7's (D) for any branch t(theta) (complex allowed): d/dtheta [B kappa(t)(1 - rho^2)] with (B) and rho' = -rho equals"
          f" (B kappa/2)(Phi (1 - rho) + 4 rho^2) on (C'), Phi = (x kappa'/H)/kappa^2 = mu mu''/mu'^2 (difference {diffe})",
          diffe == 0)


def main():
    s0()
    s1()
    s2()
    s3()
    s4()
    s5()
    P(f"\n{NCHK[0]} checks, all True (each asserted).")
    P(f"runtime {time.time() - T0:.0f} s")


if __name__ == '__main__':
    main()
