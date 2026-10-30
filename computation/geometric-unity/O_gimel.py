"""GU round 11 (O): the variation of the action (9.4) in the section gimel (the metric on X).

Formalization: formalization/geometric-unity/O-gimel-variation.md.  Builds on F (R-Z2: Z1/Z2; F1-F5), G (G1: unitary
gauge; G.2.5), I (I5/I6: the invariant sector T = t Phi^1 + u Phi^1 Gamma), K (Theorem K), L (O'Neill blocks off the
section).  Imports F_jets.py and I_elliptic_reduction.py (and through it G_symbols.py, H_ellipticity.py, shiab_checks.py)
unmodified.  Exact arithmetic throughout (sympy rationals / rational-function fields, Fractions, Gaussian rationals);
fixed seeds; one core.  Every check is asserted.

Sections (labels = results of the note)
  O0  inputs: c1, c13 of I6 recomputed by I's own i5_invariant; integrability with chi = -Gamma^2; the potential P~ and
      the reduced density l = 6 a s_Y + 14 P~ with a = x t + chi y u (O4).
  O1  the action density on the invariant sector evaluated directly in the Clifford model (O4): (1/3) Psi(T^T^T),
      (kappa/2) <T,T>, Psi(T ^ F_R) for a Levi-Civita spin curvature, and the derivative term (zero pointwise).
  O2  ingredients and measures (O1, O2): Sym^2 Jacobian det(P)^(n+1); det G_lambda; det of the Zorro metric at a point
      (Gamma-free; Z1 gimel-free); Z1 scaling = rescaling of c_H; Z2 scale invariance (brute-force scalar curvature).
  O3  Lemma O8: the Z1 scalar curvature off the section is sigma_1/c_V + tr(y^-1 Ric)/c_H + q_y(Rm/c_H), q quadratic,
      zero on the scaling ray, independent of the 3-jet (brute force, n = 3 in three signatures and n = 4, (1,3)).
  O4  the scaling ray y = c g (brute force, n = 4, general coordinates): s_Y and Ricci blocks in Z1 and Z2, Ric_HV = 0.
  O5  Lemma O9 (Euclidean Z2): |A|^2 in eigen-coordinates against brute force (n = 3); the O(4)-average of |A|^2 is
      (c_V/(24 c_H^2)) (tr u tr u^-1 - 16) |Rm|^2 (Schur traces on Weyl, Ric_0, scalar parts, n = 4).
  O6  Theorem O6 (indefinite no-go) inputs: boosts act on Sym^2 with real eigenvalues; their joint kernel is R g;
      regular orbits have finite stabiliser; rotations of a definite metric have no real nonzero eigenvalue on Sym^2;
      Lorentzian traces unbounded below; Weyl Jacobian for n = 2.
  O7  Theorem O10 / Cor. O11 / Prop. O12 / Remark O14: the gimel-Euler-Lagrange equations (sympy): Lambda for ray
      weights in Z1 and Z2, BG identification, leading order, comparison with Lambda_P, Lambda_A and I7(a); the branch
      t = -kappa/(208 x); the constant-mode system.
  O8  Proposition O13: the weight in the varpi-equation: toy Euler-Lagrange check; Psi(T ^ xi ^ T) = 0 for all 48
      types of S_sa; the grade-2 part of Shiab(xi ^ T) is a(t, u) W_A(xi ^ Phi^1)_2.
  O9  Proposition O3(b): the gimel-derivative of pi' at fixed (epsilon, varpi) is -(1 - Ad_eps) dA_0.
Run (from computation/geometric-unity): python3 O_gimel.py > O_gimel.out 2>&1   (~3 min on one core; 59 checks)
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

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import F_jets as FJ  # noqa: E402
import I_elliptic_reduction as IR  # noqa: E402

T0 = time.time()
NCHK = [0]
NTRUE = [0]
t, u, kap = sp.symbols('t u kappa')
x, y = IR.xs, IR.ys
cH, cV, lam = sp.symbols('c_H c_V lambda')
SG = {1: '(14,0)', -1: '(7,7)'}          # chi = -Gamma^2 = +1 in (14,0), -1 in (7,7)


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


def QR(v):
    """exact conversion of int / Fraction / sympy number to a sympy Rational (never nsimplify)"""
    if isinstance(v, Fr):
        return sp.Rational(v.numerator, v.denominator)
    if isinstance(v, int):
        return sp.Integer(v)
    return sp.Rational(sp.sympify(v))


def tofr(v):
    return QR(v.as_expr()) if hasattr(v, 'as_expr') and not isinstance(v, (int, Fr)) else QR(v)


def sigma1(n, l_):
    return (-sp.Rational(n * (n - 1) * (n + 2), 8)
            - sp.Rational(n, 4) * ((1 + l_) / (1 - n * l_) + sp.Rational(n + 1, 2)))


def sF(n):
    return -sp.Rational(n * (n - 1) * (n + 2), 8)          # times 1/c_V


# ====================================================================== O0
def potential(c1, c13, chi):
    s_ = sp.Symbol('s_')
    integrand = t * c1.subs({t: s_ * t, u: s_ * u}, simultaneous=True) + \
        chi * u * c13.subs({t: s_ * t, u: s_ * u}, simultaneous=True)
    return sp.expand(sp.integrate(sp.expand(integrand), (s_, 0, 1)))


def o0():
    section("O0. Inputs: c1, c13 of I6 (recomputed by I's i5_invariant); potential P~ and the reduced density l")
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        res = IR.i5_invariant()
    trans = {1: (kap * t + x * (312 * t ** 2 + 104 * u ** 2) + 208 * y * t * u,
                 kap * u + 208 * x * t * u + y * (104 * t ** 2 + 312 * u ** 2)),
             -1: (kap * t + x * (312 * t ** 2 - 104 * u ** 2) - 208 * y * t * u,
                  kap * u + 208 * x * t * u + y * (104 * t ** 2 - 312 * u ** 2))}
    tsub = {sp.Symbol('t'): t, sp.Symbol('u'): u, sp.Symbol('kappa'): kap}
    out = {}
    for chi in (1, -1):
        C1, C13 = [sp.expand(c.subs(tsub)) for c in res[chi]]
        c1, c13 = trans[chi]
        check(f"{SG[chi]}: I's i5_invariant reproduces the transcription c1 = {sp.expand(c1)}, c13 = {sp.expand(c13)}",
              sp.expand(C1 - c1) == 0 and sp.expand(C13 - c13) == 0)
        integ = sp.expand(sp.diff(c1, u) - chi * sp.diff(c13, t))
        Pt = potential(c1, c13, chi)
        okP = sp.expand(sp.diff(Pt, t) - c1) == 0 and sp.expand(sp.diff(Pt, u) - chi * c13) == 0
        P2 = sp.expand(Pt.subs({x: 0, y: 0}))
        P3 = sp.expand(Pt - P2)
        Dl = sp.expand(t * c1 + chi * u * c13)
        a = x * t + chi * y * u
        okf = sp.expand(Pt - (t ** 2 + chi * u ** 2) * (kap / 2 + 104 * a)) == 0 and \
            sp.expand(Dl - (t ** 2 + chi * u ** 2) * (kap + 312 * a)) == 0
        check(f"{SG[chi]}: d_u c1 = chi d_t c13 ({integ} = 0), so (c1, chi c13) is a gradient; P~ = {P2} + {sp.factor(P3)}"
              f" = (t^2 + chi u^2)(kappa/2 + 104 a), a = {a}; Euler: D := t c1 + chi u c13 = 2 P~_2 + 3 P~_3"
              f" = (t^2 + chi u^2)(kappa + 312 a)", integ == 0 and okP and okf and
              sp.expand(Dl - 2 * P2 - 3 * P3) == 0)
        out[chi] = dict(c1=c1, c13=c13, Pt=Pt, P2=P2, P3=P3, a=a, D=Dl)
    return out


# ====================================================================== O1
def o1(inv):
    section("O1. The action density on the invariant sector, evaluated directly in the Clifford model (I's evaluator)")
    N = IR.N
    for eta in (IR.EUC, IR.E77):
        chi = eta[7]
        Al = IR.Alg(eta)
        Phi1 = [Al.gamma(a_) for a_ in range(N)]
        Phi1G = [Al.gammaG(a_) for a_ in range(N)]
        rng = np.random.default_rng(1)
        R = IR.random_curvature(N, rng, terms=1).astype(np.int64)
        etam = np.diag(eta).astype(np.int64)
        Ric = np.einsum('ac,abcd->bd', etam, R)
        scal = int(np.einsum('bd,bd->', etam, Ric))
        ptR = Al.psi_form(IR.z4_form(Al, R, half=True))
        ok_cub = ok_kap = ok_cur = True
        pts = ((1, 0), (0, 1), (1, 1), (1, -1), (2, 1), (-3, 2))
        for tv, uv in pts:
            T = [IR.add(IR.scale(Phi1[a_], tv), IR.scale(Phi1G[a_], uv)) for a_ in range(N)]
            pf = Al.psi_form(Al.wedge11(T, T))
            P3v = inv[chi]['P3'].subs({t: tv, u: uv})
            for op, cx in (('W_A', x), ("W'_A", y), ('W_C.AdG', None), ("W'_C.AdG", None)):
                unit, sv = IR.OPS[op]
                val = IR.unit_sp(unit) * QR(Al.pair_form(IR.combine(pf, sv), T)) / 3
                exp = 14 * sp.expand(P3v).coeff(cx) if cx is not None else 0
                ok_cub &= sp.expand(val - exp) == 0
                cur = IR.unit_sp(unit) * QR(Al.pair_form(IR.combine(ptR, sv), T))
                expc = 6 * scal * (tv if op == 'W_A' else chi * uv if op == "W'_A" else 0)
                ok_cur &= sp.expand(cur - expc) == 0
            TT = sum(eta[a_] * Al.tr(T[a_], T[a_]) for a_ in range(N))
            ok_kap &= sp.Rational(TT) == 14 * (tv ** 2 + chi * uv ** 2)
        check(f"{SG[chi]}: (1/3) Psi(T^T^T) = 14 P~_3 for T = t Phi^1 + u Phi^1 Gamma at {len(pts)} points, per operator"
              f" (W_A: x-part, W'_A: y-part, twisted: 0)", ok_cub)
        check(f"{SG[chi]}: <T, T> = 14 (t^2 + chi u^2), so (kappa/2) <T,T> = 14 P~_2", ok_kap)
        check(f"{SG[chi]}: Psi(T ^ F_R) = 6 (x t + chi y u) s for the Levi-Civita spin curvature of a random algebraic"
              f" R (s = {scal}); twisted operators 0", ok_cur)
        random.seed(31 + chi)
        okd = True
        for _ in range(3):
            xi1 = [random.randint(-2, 2) for _ in range(N)]
            xi2 = [random.randint(-2, 2) for _ in range(N)]
            zt = {}
            for a_, b_ in itertools.combinations(range(N), 2):
                v = IR.add(IR.add(IR.scale(Phi1[b_], xi1[a_]), IR.scale(Phi1[a_], -xi1[b_])),
                           IR.add(IR.scale(Phi1G[b_], xi2[a_]), IR.scale(Phi1G[a_], -xi2[b_])))
                if v:
                    zt[(a_, b_)] = v
            pz = Al.psi_form(zt)
            for tv, uv in ((1, 2), (-1, 3)):
                T = [IR.add(IR.scale(Phi1[a_], tv), IR.scale(Phi1G[a_], uv)) for a_ in range(N)]
                okd &= all(Al.pair_form(pz[ty], T) == 0 for ty in pz)
        check(f"{SG[chi]}: the derivative term Psi(T ^ d_0 T), d_0 T = dt ^ Phi^1 + du ^ Phi^1 Gamma, vanishes pointwise"
              f" for all 48 types (random dt, du)", okd)


# ====================================================================== O2
def sym_basis(n):
    idx = FJ.sym_index(n)
    E = []
    for (i, j) in idx:
        M = sp.zeros(n, n)
        M[i, j] = 1
        M[j, i] = 1
        E.append(M)
    return idx, E


def sym_action_matrix(n, Pm):
    """matrix of y -> P^T y P in the coordinates y_A (i <= j)"""
    idx, E = sym_basis(n)
    cols = []
    for Eb in E:
        Yb = Pm.T * Eb * Pm          # E^B = d/dy_B (y_ij = y_ji = 1 for i < j); coordinates of Yb are Yb[i, j], i <= j
        cols.append([Yb[i, j] for (i, j) in idx])
    return sp.Matrix(cols).T


def make_base(n, eta, rng, ctx, eps=None, with_lin=True, with_cub=True, scale=2):
    def rsym(s=scale):
        M = [[None] * n for _ in range(n)]
        for i in range(n):
            for j in range(i, n):
                v = ctx.c(Fr(rng.randint(-s, s), rng.randint(1, 3)))
                M[i][j] = v
                M[j][i] = v
        return M
    e = ctx.one if eps is None else eps
    g0 = [[ctx.c(eta[i]) if i == j else ctx.zero for j in range(n)] for i in range(n)]
    lin = [rsym() if with_lin else [[ctx.zero] * n for _ in range(n)] for _ in range(n)]
    quad = {(a_, b_): [[v * e for v in row] for row in rsym()] for a_ in range(n) for b_ in range(a_, n)}
    cub = {(a_, b_, c_): [[v * e for v in row] for row in rsym()] for a_ in range(n) for b_ in range(a_, n)
           for c_ in range(b_, n)} if with_cub else {}
    return g0, lin, quad, cub


def rand_y0(n, eta, rng, ctx, g0, den=3, r=2):
    """a random point y0 != g0 of the fibre over x0: nondegenerate with the signature of eta (numeric test only for the
    choice; everything computed at y0 is exact)"""
    while True:
        y0 = [[g0[i][j] + (ctx.c(Fr(rng.randint(-r, r), den)) if i <= j else ctx.zero) for j in range(n)]
              for i in range(n)]
        for i in range(n):
            for j in range(i):
                y0[i][j] = y0[j][i]
        Ym = np.array([[float(tofr(y0[i][j])) for j in range(n)] for i in range(n)])
        ev = np.linalg.eigvalsh(Ym)
        if (min(abs(ev)) > 1e-3 and sum(ev > 0) == eta.count(1)
                and any(tofr(y0[i][j]) != tofr(g0[i][j]) for i in range(n) for j in range(n))):
            return y0


def scalar_Y(gj, base, n, y0, lamv, cHv, cVv, variant, ctx, blocks=False):
    G, m, idx, Psi0 = FJ.total_metric(gj, base, n, y0, ctx.c(lamv), ctx.c(cHv), ctx.c(cVv), variant, ctx)
    Ric, G0, G0i = FJ.ricci_at_point(G, m, ctx)
    s = FJ.trace(Ric, G0i, ctx)
    if blocks:
        return s, FJ.frame_blocks(Ric, n, idx, Psi0, ctx), G0
    return s


def o2():
    section("O2. Ingredients and measures: Sym^2 Jacobian, det G_lambda, det of the Zorro metric (Gamma-free; Z1"
            " gimel-free), the scaling y -> c y")
    # (a) Jacobian of y -> P^T y P
    ps = sp.symbols('p0:4')
    Pm = sp.Matrix(2, 2, ps)
    ok2 = sp.expand(sym_action_matrix(2, Pm).det() - Pm.det() ** 3) == 0
    rng = random.Random(7)
    okr = True
    for n in (3, 4):
        for _ in range(2):
            Pm = sp.Matrix(n, n, lambda i, j: sp.Rational(rng.randint(-3, 3), rng.randint(1, 2)))
            if Pm.det() == 0:
                Pm += sp.eye(n)
            okr &= sym_action_matrix(n, Pm).det() == Pm.det() ** (n + 1)
    check("Jacobian of y -> P^T y P on Sym^2(R^n) is det(P)^(n+1) (n = 2 symbolic; n = 3, 4 random rational P); hence"
          " |det y|^(-n/2) d^n x d^N y is invariant under linear changes of x (exponent -(n+1) + 1 + n = 0)", ok2 and okr)
    # (b) det of the DeWitt Gram matrix in the F_jets basis
    ctxl = FJ.Ctx('lam')
    L_ = ctxl.gens[0]
    okG = True
    rng = random.Random(8)
    for n, eta in ((2, (1, -1)), (3, (1, 1, -1)), (4, (1, -1, -1, -1))):
        N = n * (n + 1) // 2
        Y0 = sp.Matrix(n, n, lambda i, j: 0)
        while Y0.det() == 0:
            Y0 = sp.diag(*eta) + sp.Matrix(n, n, lambda i, j: 0)
            for i in range(n):
                for j in range(i, n):
                    v = sp.Rational(rng.randint(-2, 2), 3)
                    Y0[i, j] += v
                    if i != j:
                        Y0[j, i] += v
        Wi = Y0.inv()
        W = [[FJ.jconst(ctxl.c(Fr(int(Wi[i, j].p), int(Wi[i, j].q))), 0) for j in range(n)] for i in range(n)]
        Gm = FJ.dewitt_matrix(W, n, L_, 0, ctxl)
        M = sp.Matrix(N, N, lambda A, B: Gm[A][B][0].get((), ctxl.zero).as_expr())
        lsym = sp.Symbol('lam')
        okG &= sp.simplify(M.det() - (1 - n * lsym) * 2 ** (N - n) * Y0.det() ** (-(n + 1))) == 0
    check("det G_lambda|_y (basis E_ii, E_ij + E_ji) = (1 - n lambda) 2^(N - n) det(y)^-(n+1) (n = 2, 3, 4, lambda"
          " symbolic): the fibre-invariant measure is |det y|^-(n+1)/2 d^N y", okG)
    # (c) det of the Zorro metric at a point with Gamma(0) != 0
    ctx = FJ.Ctx('q')
    okd = True
    rng = random.Random(9)
    for n, eta, lamv in ((2, (1, -1), Fr(1, 2)), (3, (1, 1, 1), Fr(0)), (3, (1, -1, -1), Fr(1, 2))):
        N = n * (n + 1) // 2
        dets = {}
        for rep in range(2):
            g0, lin, quad, cub = make_base(n, eta, rng, ctx)
            gj = FJ.base_from_coeffs(g0, lin, quad, cub, n, ctx)
            base = FJ.base_geometry(gj, n, ctx)
            gam0 = any(base['Gam0'][a_][b_][c_] != 0 for a_ in range(n) for b_ in range(n) for c_ in range(n))
            y0 = [[g0[i][j] + (ctx.c(Fr(1, 3)) if i == j == 0 else (ctx.c(Fr(1, 5)) if (i, j) in ((0, 1), (1, 0))
                                                                     else ctx.zero)) for j in range(n)]
                  for i in range(n)]
            Y0 = sp.Matrix(n, n, lambda i, j: tofr(y0[i][j]))
            for var in ('Z1', 'Z2'):
                G, m, idx, Psi0 = FJ.total_metric(gj, base, n, y0, ctx.c(lamv), ctx.c(2), ctx.c(3), var, ctx)
                G0 = sp.Matrix(m, m, lambda a_, b_: tofr(FJ.jval(G[a_][b_], ctx)))
                h = Y0 if var == 'Z1' else sp.diag(*eta)
                pred = 2 ** n * h.det() * 3 ** N * (1 - n * sp.Rational(lamv)) * 2 ** (N - n) * Y0.det() ** (-(n + 1))
                okd &= gam0 and G0.det() == pred
                dets[(rep, var)] = G0.det()
        okd &= dets[(0, 'Z1')] == dets[(1, 'Z1')]
    check("det g_aleph at (x0, y0) with Gamma(x0) != 0 equals c_H^n det(h_H) det(c_V G_lambda): Z1 (h_H = y0) and Z2"
          " (h_H = g(x0)); for Z1 it is the same for two different random base metrics (n = 2, 3; three signatures)", okd)
    # (d) scaling: Z1(c_H) at c y0 = Z1(c c_H) at y0; Z2 invariant
    ok1 = ok2 = True
    rng = random.Random(10)
    for n, eta, lamv in ((3, (1, 1, 1), Fr(0)), (3, (1, 1, -1), Fr(1, 2))):
        g0, lin, quad, cub = make_base(n, eta, rng, ctx)
        gj = FJ.base_from_coeffs(g0, lin, quad, cub, n, ctx)
        base = FJ.base_geometry(gj, n, ctx)
        y0 = [[g0[i][j] + (ctx.c(Fr(1, 2)) if (i, j) in ((0, 1), (1, 0)) else ctx.zero) for j in range(n)]
              for i in range(n)]
        y2 = [[v * 2 for v in row] for row in y0]
        ok1 &= scalar_Y(gj, base, n, y2, lamv, 2, 3, 'Z1', ctx) == scalar_Y(gj, base, n, y0, lamv, 4, 3, 'Z1', ctx)
        ok2 &= scalar_Y(gj, base, n, y2, lamv, 2, 3, 'Z2', ctx) == scalar_Y(gj, base, n, y0, lamv, 2, 3, 'Z2', ctx)
    check("brute-force scalar curvature (n = 3, (3,0) and (2,1), y0 off the section): s_Y^Z1(c_H = 2) at 2 y0 equals"
          " s_Y^Z1(c_H = 4) at y0; s_Y^Z2 at 2 y0 equals s_Y^Z2 at y0", ok1 and ok2)


# ====================================================================== O3
def o3():
    section("O3. Lemma O8: Z1 scalar curvature off the section = sigma_1/c_V + tr(y^-1 Ric)/c_H + q_y(Rm/c_H)")
    ctx = FJ.Ctx('eps')
    eps = ctx.gens[0]
    E = sp.Symbol('eps')
    runs = [(3, (1, 1, 1), Fr(0), 21), (3, (1, 1, -1), Fr(1, 2), 22), (3, (1, -1, -1), Fr(1, 2), 23),
            (4, (1, -1, -1, -1), Fr(1, 2), 24)]
    for n, eta, lamv, seed in runs:
        t1 = time.time()
        rng = random.Random(seed)
        g0, lin, quad, cub = make_base(n, eta, rng, ctx, eps=eps, with_lin=False, with_cub=True)
        gj = FJ.base_from_coeffs(g0, lin, quad, cub, n, ctx)
        base = FJ.base_geometry(gj, n, ctx)
        if n == 3:
            gj_nc = FJ.base_from_coeffs(g0, lin, quad, {}, n, ctx)
            base_nc = FJ.base_geometry(gj_nc, n, ctx)
        y0 = rand_y0(n, eta, rng, ctx, g0)
        yr = [[g0[i][j] * 2 for j in range(n)] for i in range(n)]
        sig = (eta.count(1), eta.count(-1))

        def poly(sv):
            ex = sp.cancel(sv.as_expr())
            num, den = sp.fraction(ex)
            assert not den.has(E)
            pl = sp.Poly(num, E)
            return [sp.Rational(pl.coeff_monomial(E ** k)) / den for k in range(3)], pl.degree()
        s_off = scalar_Y(gj, base, n, y0, lamv, 2, 3, 'Z1', ctx)
        (c0, c1, c2), deg = poly(s_off)
        Y0i = FJ.finv(y0, ctx)
        trR = sum((Y0i[a_][b_] * base['Ric'][b_][a_] for a_ in range(n) for b_ in range(n)), ctx.zero)
        trR1 = sp.Rational(sp.cancel(trR.as_expr() / E))
        ok = deg <= 2 and c0 == sigma1(n, sp.Rational(lamv)) / 3 and c1 == trR1 / 2 and c2 != 0
        s_ray = scalar_Y(gj, base, n, yr, lamv, 2, 3, 'Z1', ctx)
        (r0, r1, r2), _ = poly(s_ray)
        Rs1 = sp.Rational(sp.cancel(base['Rs'].as_expr() / E))
        okr = r0 == sigma1(n, sp.Rational(lamv)) / 3 and r1 == Rs1 / (2 * 2) and r2 == 0
        extra = ''
        if n == 3:
            s_nc = scalar_Y(gj_nc, base_nc, n, y0, lamv, 2, 3, 'Z1', ctx)
            s_cH = scalar_Y(gj, base, n, y0, lamv, 4, 3, 'Z1', ctx)
            (d0, d1, d2), _ = poly(s_cH)
            okx = s_nc == s_off and d0 == c0 and d1 == c1 / 2 and d2 == c2 / 4
            ok &= okx
            extra = (f"; the same with the cubic terms of g removed (3-jet irrelevant): {s_nc == s_off}; at c_H = 4:"
                     f" linear coefficient halves, quadratic quarters: {d1 == c1 / 2 and d2 == c2 / 4}")
        check(f"n = {n}, base {sig}, lambda = {lamv}, c_H = 2, c_V = 3, g = g0 + eps (quad + cub), y0 off the section:"
              f" s_Y^Z1 = {c0} + ({c1}) eps + ({c2}) eps^2 with c0 = sigma_1/c_V, c1 = tr(y0^-1 Ric_1)/c_H, c2 != 0{extra}"
              f"  [{time.time() - t1:.0f} s]", ok)
        check(f"n = {n}, base {sig}: on the scaling ray y0 = 2 g(x0) the quadratic term vanishes: s_Y^Z1 = sigma_1/c_V +"
              f" R/(2 c_H) eps exactly (F4(b) for the section 2g)", okr)


# ====================================================================== O4
def o4():
    section("O4. The scaling ray y = c g(x0) (n = 4, general coordinates): s_Y and Ricci blocks of Z1 and Z2")
    ctx = FJ.Ctx('q')
    n = 4
    for eta, cval, variants, seed in (((1, -1, -1, -1), 2, ('Z1', 'Z2'), 41), ((1, 1, -1, -1), -2, ('Z1',), 42)):
        rng = random.Random(seed)
        g0, lin, quad, cub = make_base(n, eta, rng, ctx)
        gj = FJ.base_from_coeffs(g0, lin, quad, cub, n, ctx)
        base = FJ.base_geometry(gj, n, ctx)
        lamv = Fr(1, 2)
        yc = [[g0[i][j] * cval for j in range(n)] for i in range(n)]
        Rs = tofr(base['Rs'])
        Ricg = sp.Matrix(n, n, lambda a_, b_: tofr(base['Ric'][a_][b_]))
        sig = (eta.count(1), eta.count(-1))
        for var in variants:
            t1 = time.time()
            s, RF, G0 = scalar_Y(gj, base, n, yc, lamv, 2, 3, var, ctx, blocks=True)
            N = n * (n + 1) // 2
            hv0 = all(RF[a_][n + B] == 0 for a_ in range(n) for B in range(N))
            HH = sp.Matrix(n, n, lambda a_, b_: tofr(RF[a_][b_]))
            if var == 'Z1':
                s_pred = Rs / (cval * 2) + sigma1(4, sp.Rational(lamv)) / 3
                HH_pred = Ricg + sp.Rational(cval * 2, 4 * 3) / (4 * sp.Rational(lamv) - 1) * sp.diag(*eta)
            else:
                s_pred = Rs / 2 + sF(4) / 3
                HH_pred = Ricg
            check(f"base {sig}, {var}, lambda = 1/2, c_H = 2, c_V = 3, y = {cval} g(x0), random cubic base metric in"
                  f" general coordinates (Gamma(x0) != 0): s_Y = {'R/(c c_H) + sigma_1/c_V' if var == 'Z1' else 'R/c_H + s_F'};"
                  f" Ric_HV = 0; Ric_HH = {'Ric_g + c c_H g/(4 c_V (4 lambda - 1))' if var == 'Z1' else 'Ric_g'}"
                  f"  [{time.time() - t1:.0f} s]", tofr(s) == s_pred and hv0 and HH == HH_pred)


# ====================================================================== O5
def curv_tensor_space(n):
    """basis of algebraic curvature tensors (dicts over (k,l,i,j)), from symmetric maps on Lambda^2 + Bianchi"""
    pairs = list(itertools.combinations(range(n), 2))
    P2 = len(pairs)
    params = [(A, B) for A in range(P2) for B in range(A, P2)]

    def tensor(vec):
        Rt = {}
        for (A, B), v in zip(params, vec):
            if v == 0:
                continue
            for (I, J) in ((A, B), (B, A)) if A != B else ((A, B),):
                (k, l), (i, j) = pairs[I], pairs[J]
                for (kk, ll, sk) in ((k, l, 1), (l, k, -1)):
                    for (ii, jj, si) in ((i, j, 1), (j, i, -1)):
                        Rt[(kk, ll, ii, jj)] = Rt.get((kk, ll, ii, jj), 0) + sk * si * v
        return Rt
    # Bianchi: R_{kl ij} + R_{k ij l} + R_{k j l i} = 0
    units = [tensor([1 if q == p else 0 for q in range(len(params))]) for p in range(len(params))]
    rows = []
    for k, l, i, j in itertools.product(range(n), repeat=4):
        row = [Rt.get((k, l, i, j), 0) + Rt.get((k, i, j, l), 0) + Rt.get((k, j, l, i), 0) for Rt in units]
        if any(row):
            rows.append(row)
    ns = sp.Matrix(rows).nullspace()
    return [tensor(list(v)) for v in ns]


def ip(R1, R2):
    return sum(v * R2.get(k, 0) for k, v in R1.items())


def ricci_of(Rt, n):
    return sp.Matrix(n, n, lambda b, d: sum(Rt.get((a_, b, a_, d), 0) for a_ in range(n)))


def kn_tensor(h, g, n):
    Rt = {}
    for a_, b, c, d in itertools.product(range(n), repeat=4):
        v = h[a_, c] * g[b, d] + h[b, d] * g[a_, c] - h[a_, d] * g[b, c] - h[b, c] * g[a_, d]
        if v != 0:
            Rt[(a_, b, c, d)] = v
    return Rt


def o5():
    section("O5. Lemma O9 (Euclidean Z2): |A|^2 in eigen-coordinates; its O(4)-average is proportional to |Rm|^2")
    # (a) brute force n = 3: s_Y^Z2 = R/c_H + s_F - |A|^2 at diagonal y0, g0 = I, Gamma(0) = 0
    ctx = FJ.Ctx('q')
    n = 3
    rng = random.Random(51)
    g0, lin, quad, cub = make_base(n, (1, 1, 1), rng, ctx, with_lin=False)
    gj = FJ.base_from_coeffs(g0, lin, quad, cub, n, ctx)
    base = FJ.base_geometry(gj, n, ctx)
    lv = [sp.Rational(1, 2), sp.Integer(2), sp.Integer(3)]
    y0 = [[ctx.c(Fr(int(sp.Rational(lv[i]).p), int(sp.Rational(lv[i]).q))) if i == j else ctx.zero for j in range(n)]
          for i in range(n)]
    s = tofr(scalar_Y(gj, base, n, y0, Fr(0), 2, 3, 'Z2', ctx))
    Rl = {k: tofr(v) for k, v in base['R'].items()}
    f = lambda k, l: lv[k] / lv[l] + lv[l] / lv[k] - 2  # noqa: E731
    A2 = sp.Rational(3, 4 * 4) * sum(Rl[(k, l, i, j)] ** 2 * f(k, l) for k, l, i, j in itertools.product(range(n), repeat=4))
    pred = tofr(base['Rs']) / 2 + sF(3) / 3 - A2
    check("n = 3, base (3,0), Z2, lambda = 0, c_H = 2, c_V = 3, y0 = diag(1/2, 2, 3): brute-force s_Y = R/c_H + s_F -"
          " |A|^2 with |A|^2 = c_V/(4 c_H^2) sum_{ijkl} R_klij^2 (u_k/u_l + u_l/u_k - 2)", s == pred)
    # (b) Schur traces n = 4
    n = 4
    basis = curv_tensor_space(n)
    Ric = [ricci_of(Rt, n) for Rt in basis]
    W = sp.Matrix([[Rc[b, d] for b in range(n) for d in range(n)] for Rc in Ric]).T.nullspace()
    Wb = [{k: sum(vv[i] * basis[i].get(k, 0) for i in range(len(basis))) for k in
           set().union(*[set(B_.keys()) for B_ in basis])} for vv in W]
    gI = sp.eye(n)
    Hb = []
    for i in range(n):
        for j in range(i, n):
            h = sp.zeros(n, n)
            h[i, j] = h[j, i] = 1
            if i == j:
                if i == n - 1:
                    continue
                h[n - 1, n - 1] = -1
            Hb.append(kn_tensor(h, gI, n))
    Sb = [kn_tensor(gI, gI, n)]
    okdim = len(basis) == 20 and len(Wb) == 10 and len(Hb) == 9
    rng = random.Random(52)
    okall = True
    msg = []
    for _ in range(3):
        lv = [sp.Rational(rng.randint(1, 9), rng.randint(1, 4)) for _ in range(n)]
        fkl = {(k, l): lv[k] / lv[l] + lv[l] / lv[k] - 2 for k in range(n) for l in range(n)}
        Sf = sum(fkl[(k, l)] for k in range(n) for l in range(n) if k != l)
        trs = []

        def Q(R1, R2):
            return sum(v * R2.get(kk, 0) * fkl[(kk[0], kk[1])] for kk, v in R1.items())
        for V in (Wb, Hb, Sb):
            Gm = sp.Matrix(len(V), len(V), lambda i, j: ip(V[i], V[j]))
            Qm = sp.Matrix(len(V), len(V), lambda i, j: Q(V[i], V[j]))
            trs.append((Gm.inv() * Qm).trace() / len(V))
        tu = sum(lv)
        tui = sum(1 / v for v in lv)
        okall &= all(tr == Sf / 12 for tr in trs) and Sf == 2 * (tu * tui - n ** 2)
        msg.append(str(Sf / 12))
    check(f"n = 4, Euclidean: the space of algebraic curvature tensors has dim 20 = 10 (Weyl) + 9 (Ric_0) + 1 ({okdim});"
          f" for Q_u(R) = sum R_klij^2 f_kl at three random diagonal u, tr(Q|_V)/dim V is the same on all three"
          f" components (= sum_(k!=l) f_kl / 12 = {', '.join(msg)}), so avg_O(4) Q_u = (sum f/12) |Rm|^2, and"
          f" sum_(k!=l) f_kl = 2 (tr u tr u^-1 - n^2)", okdim and okall)


# ====================================================================== O6
def so_basis(eta):
    n = len(eta)
    out = []
    for i in range(n):
        for j in range(i + 1, n):
            X = sp.zeros(n, n)
            X[i, j] = eta[i]
            X[j, i] = -eta[j]
            out.append((X, eta[i] * eta[j] < 0))
    return out


def lie_on_sym(X, n):
    idx, E = sym_basis(n)
    cols = []
    for Eb in E:
        Yb = X.T * Eb + Eb * X
        cols.append([Yb[i, j] for (i, j) in idx])
    return sp.Matrix(cols).T


def o6():
    section("O6. Theorem O6 inputs: boosts on Sym^2, joint kernel R g, regular orbits, definite case, n = 2 Weyl Jacobian")
    for eta in ((1, -1), (1, 1, -1), (1, -1, -1), (1, -1, -1, -1), (1, 1, 1, -1), (1, 1, -1, -1), (1, 1, -1, -1, -1)):
        n = len(eta)
        etam = sp.diag(*eta)
        okso = True
        mats_b = []
        okeig = True
        for X, boost in so_basis(eta):
            okso &= X.T * etam + etam * X == sp.zeros(n, n)
            Mx = lie_on_sym(X, n)
            if boost:
                mats_b.append(Mx)
                ev = Mx.eigenvals()
                okeig &= all(sp.im(e) == 0 for e in ev) and set(ev) <= {-2, -1, 0, 1, 2} and Mx.is_diagonalizable()
        K = sp.Matrix.vstack(*mats_b).nullspace()
        idx, _ = sym_basis(n)
        gvec = sp.Matrix([etam[i, j] for (i, j) in idx])
        okker = len(K) == 1 and sp.Matrix.hstack(K[0], gvec).rank() == 1
        yreg = sp.diag(*[eta[i] * (i + 1 + (i == n - 1)) for i in range(n)])
        imgs = sp.Matrix.hstack(*[sp.Matrix([(X.T * yreg + yreg * X)[i, j] for (i, j) in idx]) for X, _ in so_basis(eta)])
        okreg = imgs.rank() == n * (n - 1) // 2
        sig = (eta.count(1), eta.count(-1))
        check(f"base {sig}: {len(mats_b)} boosts act on Sym^2 diagonalisably with real eigenvalues in {{0, +-1, +-2}}; their"
              f" joint kernel is R g ({okker}); the so(g)-orbit map at a regular y has rank dim so = {n * (n - 1) // 2}"
              f" (finite stabiliser): {okreg}", okso and okeig and okker and okreg)
    eta = (1, 1, 1, 1)
    oke = True
    for X, boost in so_basis(eta):
        ev = lie_on_sym(X, 4).eigenvals()
        oke &= (not boost) and all(e == 0 or sp.re(e) == 0 for e in ev)
    check("base (4,0): every generator of so(4) acts on Sym^2 with purely imaginary or zero eigenvalues (no boost; K"
          " compact)", oke)
    ytest = sp.diag(sp.Matrix([[-10, 11], [11, 10]]), -sp.eye(2))
    etam = sp.diag(1, -1, -1, -1)
    U = etam * ytest
    evs = list(sp.Matrix([[U[0, 0], U[0, 1]], [U[1, 0], U[1, 1]]]).eigenvals())
    inertia = (sum(m for e, m in ytest.eigenvals().items() if e > 0), sum(m for e, m in ytest.eigenvals().items() if e < 0))
    check(f"base (1,3): y = [[-10, 11], [11, 10]] + (-I_2) has signature {inertia}, tr(g^-1 y) = {U.trace()} and"
          f" g^-1 y has the complex eigenvalues {sorted(evs, key=str)}: tr(g^-1 (c y)) -> -oo, so e^(-beta tr u) is not"
          f" integrable on F_x", inertia == (1, 3) and U.trace() == -18 and all(sp.im(e) != 0 for e in evs))
    l1, l2, th = sp.symbols('l1 l2 theta', positive=True)
    Rm = sp.Matrix([[sp.cos(th), -sp.sin(th)], [sp.sin(th), sp.cos(th)]])
    Yp = Rm * sp.diag(l1, l2) * Rm.T
    Jac = sp.Matrix([[sp.diff(Yp[0, 0], v), sp.diff(Yp[0, 1], v), sp.diff(Yp[1, 1], v)] for v in (l1, l2, th)])
    jd = sp.simplify(Jac.det())
    a1, a2 = sp.symbols('a1 a2', real=True)
    dens = sp.simplify((sp.exp(a1) - sp.exp(a2)) * sp.exp(-(a1 + a2) / 2) - 2 * sp.sinh((a1 - a2) / 2))
    check(f"n = 2 polar coordinates y = R(theta) diag(l1, l2) R(theta)^T: |Jacobian| = |{jd}|, so the invariant measure"
          f" |det y|^(-3/2) dy = |l1 - l2| (l1 l2)^(-3/2) dl dtheta = 2 |sinh((a1 - a2)/2)| da dtheta (l = e^a)",
          sp.simplify(jd ** 2 - (l1 - l2) ** 2) == 0 and sp.simplify(dens.rewrite(sp.exp)) == 0)


# ====================================================================== O7
def o7(inv):
    section("O7. The gimel-Euler-Lagrange equation on the invariant sector (constant t, u): Lambda for ray weights")
    n = 4
    s1 = sigma1(4, lam)
    LA_Z2 = sp.simplify(-cH * sF(4) / cV / 2)
    LA_Z1 = sp.simplify(-cH * s1 / (2 * cV))
    LP_Z1 = sp.simplify(LA_Z1 + cH / (4 * cV * (n * lam - 1)))
    check(f"F5: Lambda_A(Z2) = {LA_Z2}; Lambda_A(Z1) = {sp.factor(LA_Z1)}; Lambda_P(Z1) = {sp.factor(LP_Z1)}; values at"
          f" lambda = 0, 1/2, 1 (units c_H/c_V): A {[LA_Z1.subs({lam: v, cH: 1, cV: 1}) for v in (0, sp.Rational(1, 2), 1)]},"
          f" P {[LP_Z1.subs({lam: v, cH: 1, cV: 1}) for v in (0, sp.Rational(1, 2), 1)]}",
          LA_Z2 == sp.Rational(9, 2) * cH / cV
          and sp.simplify(LA_Z1 - sp.Rational(5, 4) * cH / cV * (5 - 18 * lam) / (1 - 4 * lam)) == 0
          and sp.simplify(LP_Z1 - sp.Rational(3, 2) * cH / cV * (4 - 15 * lam) / (1 - 4 * lam)) == 0)
    Rsym, M1, M2, m0, eps = sp.symbols('R M_1 M_2 m_0 epsilon')
    out = {}
    for chi in (1, -1):
        d = inv[chi]
        a, Pt = d['a'], d['Pt']
        # Z2: on the ray s_Y = R/c_H + s_F, rho = c_H^2; nu = phi(s) ds, m0 = int phi
        L2 = cH ** 2 * m0 * (6 * a * (Rsym / cH + sF(4) / cV) + 14 * Pt)
        al2, be2 = sp.diff(L2, Rsym), sp.expand(L2 - Rsym * sp.diff(L2, Rsym))
        Lam2 = sp.simplify(-be2 / (2 * al2))
        ok2 = sp.simplify(Lam2 - (LA_Z2 - sp.Rational(7, 6) * cH * Pt / a)) == 0 and not Lam2.has(m0)
        # Z1: s_Y = e^-s R/c_H + sigma_1/c_V, rho = c_H^2 e^{2s}; M_k = int phi e^{ks}
        L1 = cH ** 2 * (6 * a * (M1 * Rsym / cH + M2 * s1 / cV) + 14 * Pt * M2)
        al1, be1 = sp.diff(L1, Rsym), sp.expand(L1 - Rsym * sp.diff(L1, Rsym))
        Lam1 = sp.simplify(-be1 / (2 * al1))
        ok1 = sp.simplify(Lam1 - M2 / M1 * (LA_Z1 - sp.Rational(7, 6) * cH * Pt / a)) == 0
        # leading order in T
        lead2 = sp.limit(Lam2.subs({t: eps * t, u: eps * u}, simultaneous=True), eps, 0)
        lead1 = sp.limit(Lam1.subs({t: eps * t, u: eps * u}, simultaneous=True), eps, 0)
        okl = sp.simplify(lead2 - LA_Z2) == 0 and sp.simplify(lead1 - M2 / M1 * LA_Z1) == 0
        ser2 = sp.series(Lam2.subs({t: eps * t, u: eps * u}, simultaneous=True), eps, 0, 2).removeO()
        corr2 = sp.factor(sp.expand(ser2 - LA_Z2).coeff(eps))
        check(f"{SG[chi]}, ray weight nu = phi(s) ds on y = e^s g: Z2: G + Lambda g = 0 with Lambda = Lambda_A(Z2) -"
              f" (7 c_H/6) P~/a, independent of phi; Z1: Lambda = (M_2/M_1) (Lambda_A(Z1) - (7 c_H/6) P~/a),"
              f" M_k = int phi e^(ks) ds; as T -> 0: Lambda -> Lambda_A (Z2), (M_2/M_1) Lambda_A (Z1); first-order"
              f" correction (Z2) = ({corr2}) |T|", ok2 and ok1 and okl)
        # BG identification in Z2
        D = d['D']
        LamBG = LA_Z2 - cH * (14 * Pt - D) / (11 * a)
        # derive: main + (1/2) <dg, t E1 + chi u E13>_HH on the ray, E_HH = x_i (Ric - s_Y c_H g/2) - c_i c_H g
        Gs, gs = sp.symbols('G g_')    # tensor placeholders: Ric - R g/2 = G; we track coefficients of G and g
        main = -6 * a * cH * Gs + cH ** 2 * sp.Rational(1, 2) * (6 * a * sF(4) / cV + 14 * Pt) * gs
        Ehh = a * (Gs - sp.Rational(1, 2) * cH * sF(4) / cV * gs) - cH * D * gs
        tot = sp.expand(main + cH ** 2 * sp.Rational(1, 2) / cH * Ehh)
        LamBG_der = sp.simplify(tot.coeff(gs) / tot.coeff(Gs))
        leadBG = sp.limit(LamBG.subs({t: eps * t, u: eps * u}, simultaneous=True), eps, 0)
        serBG = sp.series(LamBG.subs({t: eps * t, u: eps * u}, simultaneous=True), eps, 0, 2).removeO()
        corrBG = sp.factor(sp.expand(serBG - LA_Z2).coeff(eps))
        check(f"{SG[chi]}, Z2 with the Bourguignon-Gauduchon identification (varpi fixed): Lambda' = Lambda_A(Z2) -"
              f" c_H (14 P~ - D)/(11 a), D = t c1 + chi u c13; equal to Lambda at leading order (Lambda_A), first-order"
              f" correction ({corrBG}) |T|; closed forms Lambda = Lambda_A - c_H (t^2 + chi u^2)(7 kappa/(12 a) + 364/3),"
              f" Lambda' = Lambda_A - c_H (t^2 + chi u^2)(6 kappa/(11 a) + 104)",
              sp.simplify(LamBG_der - LamBG) == 0 and sp.simplify(leadBG - LA_Z2) == 0
              and sp.simplify(corrBG - corr2) != 0
              and sp.simplify(Lam2 - (LA_Z2 - cH * (t ** 2 + chi * u ** 2) * (7 * kap / (12 * a) + sp.Rational(364, 3)))) == 0
              and sp.simplify(LamBG - (LA_Z2 - cH * (t ** 2 + chi * u ** 2) * (6 * kap / (11 * a) + 104))) == 0)
        # comparison with (P) of I7(a) at y = 0
        c1 = d['c1']
        LamP4 = sp.Symbol('Lambda_P') - cH * c1 / x
        corrP = sp.factor(sp.expand(sp.series(LamP4.subs({t: eps * t, u: eps * u}, simultaneous=True), eps, 0, 2)
                                    .removeO()).coeff(eps))
        check(f"{SG[chi]}, y = 0 (Shiab#, W_A): (P) of I7(a) gives Lambda_4 = Lambda_P - c_H c1/x (first order"
              f" {corrP}); the gimel-equation gives Lambda_A + ({sp.factor(corr2.subs(y, 0))}) at first order",
              sp.simplify(corrP + cH * kap * t / x) == 0)
        # the branch t = -kappa/(208 x), y = 0
        Pty = Pt.subs(y, 0)
        br = sp.factor(Pty.subs(t, -kap / (208 * x)))
        check(f"{SG[chi]}, y = 0: on the W_A branch t = -kappa/(208 x) of I6, P~ = {br}, so the gimel-equation (ray weight)"
              f" has Lambda = Lambda_A exactly (Z2), (M_2/M_1) Lambda_A (Z1), for every u", br == 0)
        out[chi] = dict(Lam2=Lam2, LamBG=LamBG, Lam1=Lam1)
    # constant-mode system (Remark O14), Z2, y = 0
    for chi in (1, -1):
        d = inv[chi]
        a, Pt, c1, c13 = d['a'], d['Pt'], d['c1'], d['c13']
        sFv = sF(4) / cV
        Lam = LA_Z2 - sp.Rational(7, 6) * cH * Pt / a
        Et = sp.simplify((6 * x * (4 * Lam / cH + sFv) + 14 * c1).subs(y, 0))
        Eu = sp.factor((14 * c13).subs(y, 0))
        u2 = sp.solve(sp.numer(sp.together(Et.subs(t, -kap / (208 * x)))), u)
        u2v = sp.simplify(u2[0] ** 2) if u2 else None
        t2 = sp.solve(sp.numer(sp.together(Et.subs(u, 0))), t)
        t2v = sorted({sp.simplify(v ** 2) for v in t2}, key=str)
        u2pred = kap ** 2 / (43264 * x ** 2) - sp.Rational(3, 728) * sFv if chi == -1 else \
            sp.Rational(3, 728) * sFv - kap ** 2 / (43264 * x ** 2)
        okc = sp.simplify(u2v - u2pred) == 0
        okc &= t2v == [sp.Rational(3, 728) * sFv]
        LamU0 = sp.factor(Lam.subs({y: 0, u: 0}))
        check(f"{SG[chi]}, constant-mode system (Z2, ray weight, y = 0): u-equation {Eu} = 0; branch t = -kappa/(208 x):"
              f" u^2 = {u2v}, Lambda = Lambda_A(Z2); branch u = 0: t^2 = {t2v[0]}, Lambda = {LamU0}", okc)
    # (14,0) real form of Lambda
    xr, yr, sr, rr = sp.symbols("x' y' s r", real=True)
    Lr = sp.factor(sp.cancel(sp.expand((out[1]['Lam2'] - LA_Z2).subs(
        {x: -sp.I * xr, y: -sp.I * yr, t: sp.I * sr, u: sp.I * rr}, simultaneous=True))))
    check(f"(14,0) real form (x = -i x', y = -i y', t = i s, u = i r): Lambda - Lambda_A(Z2) = {Lr} is real",
          not Lr.has(sp.I))
    # Corollary O15: on-shell (S2) at the round-10 solutions (N12(c): Z1, (14,0), lambda = 1/5, flat base, T_e = 0)
    d = inv[1]
    l5 = sp.Rational(1, 5)
    LA15 = sp.simplify(LA_Z1.subs(lam, l5))
    GY = sp.Rational(15, 2) / cV                       # G(g_aleph) = (15/(2 c_V)) g_aleph (F7, N2)
    ups = sp.expand(d['c1'] - x * GY).subs({y: 0, u: 0})          # x G = c1 g; c13 = u (kappa + 208 x t) = 0 on u = 0
    LamO1 = sp.together((LA15 - sp.Rational(7, 6) * cH * d['Pt'] / d['a']).subs({y: 0, u: 0}))
    sol = sp.solve([ups, sp.numer(LamO1)], [t, kap], dict=True)
    xr_ = sp.Symbol("x'", positive=True)
    okO15 = len(sol) == 2
    for so in sol:
        tt, kk = so[t], so[kap]
        okO15 &= sp.simplify((kk ** 2 * cV).subs(x, -sp.I * xr_) - 12480 * xr_ ** 2) == 0      # c_V kappa^2 = 12480 x'^2
        okO15 &= sp.simplify((tt / sp.I + kk / (416 * x * sp.I)).subs(x, -sp.I * xr_)) == 0    # s = -kappa/(416 x')
    # det M = kappa^2 (kappa + 224 x' s)(kappa + 544 x' s) on u = 0 (J4(d) proof, N12), at s = -kappa/(416 x')
    s_star = -kap / (416 * xr_)
    dM = sp.factor((kap ** 2 * (kap + 224 * xr_ * s_star) * (kap + 544 * xr_ * s_star)))
    okO15 &= sp.simplify(LA15 - sp.Rational(35, 4) * cH / cV) == 0 and dM != 0
    check(f"(14,0), Z1, lambda = 1/5, flat base (N12(c)): Lambda_A(Z1) = {LA15}; the Upsilon = 0 condition x G = c1 g"
          f" (G = 15/(2 c_V) g, u = 0) together with Lambda_O(Z1) = 0 (the gimel-equation at a flat base, every admissible"
          f" weight) has the solutions c_V kappa^2 = 12480 x'^2, s = -kappa/(416 x') (t = i s, x = -i x'); det M = {dM} != 0",
          okO15)
    # Euclidean Z2 with a smooth natural weight (Theorem O10(b)): reduced Lagrangian with the O(4)-averaged |A|^2
    Fp, Ric0sq, GB = sp.symbols("F' |Ric_0|^2 E_GB")
    d = inv[1]
    a, Pt = d['a'], d['Pt']
    Lw = cH ** 2 * (6 * a * (m0 * Rsym / cH - cV / (24 * cH ** 2) * Fp * (GB + 4 * Ric0sq)) + m0 * (6 * a * sF(4) / cV
                                                                                               + 14 * Pt))
    al, ga = sp.diff(Lw, Rsym), sp.diff(Lw, Ric0sq)
    be = sp.expand(Lw - Rsym * al - Ric0sq * ga - GB * sp.diff(Lw, GB))
    Lamw = sp.simplify(-be / (2 * al))
    kw = sp.simplify(-ga / al)
    check(f"(14,0), Z2, smooth natural weight w (m_0 = int w, F' = int (tr u tr u^-1 - 16) w): reduced Lagrangian"
          f" alpha R + beta + gamma |Ric_0|^2 + (Gauss-Bonnet); EL: G + Lambda g + kappa_w E_|Ric_0|^2 = 0 with"
          f" Lambda = Lambda_A(Z2) - (7 c_H/6) P~/a (as for ray weights) and kappa_w = {kw}",
          sp.simplify(Lamw - out[1]['Lam2']) == 0 and sp.simplify(kw - cV * Fp / (6 * cH * m0)) == 0)


# ====================================================================== O8
def o8():
    section("O8. Proposition O13: a weight in the action changes the varpi-equation")
    # (a) toy: L = 1/2 w b_ij eps^{abc} T^i_a d_b T^j_c on R^3
    X = sp.symbols('X0:3')
    w = sp.Function('w')(*X)
    Tf = [[sp.Function(f'T{i}{a_}')(*X) for a_ in range(3)] for i in range(2)]
    b = sp.Matrix([[2, 1], [1, -3]])
    L = sp.Rational(1, 2) * w * sum(b[i, j] * sp.LeviCivita(a_, bb, c) * Tf[i][a_] * sp.diff(Tf[j][c], X[bb])
                                    for i in range(2) for j in range(2) for a_ in range(3) for bb in range(3)
                                    for c in range(3))
    ok = True
    for i in range(2):
        for a_ in range(3):
            f = Tf[i][a_]
            el = sp.diff(L, f) - sum(sp.diff(sp.diff(L, sp.diff(f, X[k])), X[k]) for k in range(3))
            pred = sum(b[i, j] * sp.LeviCivita(a_, bb, c) * (w * sp.diff(Tf[j][c], X[bb]) +
                                                            sp.Rational(1, 2) * sp.diff(w, X[bb]) * Tf[j][c])
                       for j in range(2) for bb in range(3) for c in range(3))
            ok &= sp.simplify(sp.expand(el - pred)) == 0
    check("toy (R^3, two components, symmetric b): the Euler-Lagrange operator of (1/2) w b(T ^ dT) is w b(*dT) +"
          " (1/2) b(*(dw ^ T)), i.e. w Upsilon + (1/2) Shiab(dw ^ T) in the GU notation", ok)
    # (b) Psi(T ^ xi ^ T) = 0 for all 48 types; (c) grade-2 part of Shiab(xi ^ T)
    N = IR.N
    for eta in (IR.EUC, IR.E77):
        chi = eta[7]
        Al = IR.Alg(eta)
        random.seed(81 + chi)
        okz = True
        nzc = 0
        for _ in range(2):
            T = [IR.rand_elt(list(range(15)), 8) for _ in range(N)]
            al = [IR.rand_elt(list(range(15)), 8) for _ in range(N)]
            xi = [random.randint(-2, 2) for _ in range(N)]
            zt = {}
            for a_, b_ in itertools.combinations(range(N), 2):
                v = IR.add(IR.scale(T[b_], xi[a_]), IR.scale(T[a_], -xi[b_]))
                if v:
                    zt[(a_, b_)] = v
            pz = Al.psi_form(zt)
            okz &= all(Al.pair_form(pz[ty], T) == 0 for ty in pz)
            nzc += sum(1 for ty in pz if Al.pair_form(pz[ty], al) != 0)
        check(f"{SG[chi]}: Psi_t(T ^ (xi ^ T)) = 0 for all 48 basis types of S_sa (two random dense T, random xi), while"
              f" Psi_t(alpha ^ (xi ^ T)) != 0 for {nzc} (type, sample) pairs: the symbol is antisymmetric, so"
              f" (1/2) w Psi(T ^ d_0 T) = (1/2) Psi(T~ ^ d_0 T~), T~ = w^(1/2) T", okz and nzc > 0)
        Phi1 = [Al.gamma(a_) for a_ in range(N)]
        Phi1G = [Al.gammaG(a_) for a_ in range(N)]
        okg = True
        for c in (0, 5, 9):
            def zeta_of(T):
                return {(min(a_, c), max(a_, c)): IR.scale(T[a_], 1 if c < a_ else -1) for a_ in range(N) if a_ != c and T[a_]}
            ref = IR.combine(Al.psi_form(zeta_of(Phi1)), IR.OPS['W_A'][1])
            ref = [{k: IR.unit_sp(IR.OPS['W_A'][0]) * QR(v) for k, v in Da.items() if IR.pc(k) == 2}
                   for Da in ref]
            for (xv, yv, zv, wv) in ((1, 0, 1, 0), (2, -1, 3, 1), (0, 1, 0, 2)):
                for tv, uv in ((1, 0), (2, 3), (1, chi)):
                    T = [IR.add(IR.scale(Phi1[a_], tv), IR.scale(Phi1G[a_], uv)) for a_ in range(N)]
                    pt = Al.psi_form(zeta_of(T))
                    tot = [dict() for _ in range(N)]
                    for op, cf in zip(IR.OPN, (xv, yv, zv, wv)):
                        D = IR.combine(pt, IR.OPS[op][1])
                        for a_ in range(N):
                            for k, v in D[a_].items():
                                if IR.pc(k) == 2:
                                    tot[a_][k] = tot[a_].get(k, 0) + cf * IR.unit_sp(IR.OPS[op][0]) * QR(v)
                    av = xv * tv + chi * yv * uv
                    okg &= all({k: sp.expand(v) for k, v in tot[a_].items() if sp.expand(v) != 0} ==
                               {k: sp.expand(av * v) for k, v in ref[a_].items() if sp.expand(av * v) != 0}
                               for a_ in range(N))
                    okg &= any(ref[a_] for a_ in range(N))
        check(f"{SG[chi]}: the grade-2 part of Shiab(xi ^ (t Phi^1 + u Phi^1 Gamma)) equals a(t, u) W_A(xi ^ Phi^1)_2,"
              f" a = x t + chi y u (members (1,0,1,0), (2,-1,3,1), (0,1,0,2); xi = e^0, e^5, e^9); W_A(xi ^ Phi^1)_2 != 0",
              okg)


# ====================================================================== O9
def o9():
    section("O9. Proposition O3(b): the gimel-derivative of pi' at fixed (epsilon, varpi)")
    e = sp.Matrix(2, 2, sp.symbols('e0:4'))
    vp = sp.Matrix(2, 2, sp.symbols('v0:4'))
    de = sp.Matrix(2, 2, sp.symbols('d0:4'))
    A0 = sp.Matrix(2, 2, sp.symbols('A0:4'))
    dA = sp.Matrix(2, 2, sp.symbols('B0:4'))
    s_ = sp.Symbol('s_')
    ei = e.inv()

    def pip(A):
        return e * vp * ei - (de + A * e - e * A) * ei
    der = (pip(A0 + s_ * dA).diff(s_)).applyfunc(sp.simplify)
    pred = -dA + e * dA * ei
    ok = (der - pred).applyfunc(sp.simplify) == sp.zeros(2, 2)
    at1 = (pred.subs(dict(zip(e, (1, 0, 0, 1))))).applyfunc(sp.simplify) == sp.zeros(2, 2)
    check("pi'(eps, varpi) = eps varpi eps^-1 - (d_0 eps) eps^-1 with d_0 eps = d eps + A_0 eps - eps A_0: d/ds pi' along"
          " A_0 + s dA_0 is -(1 - Ad_eps) dA_0 (2x2 symbolic), which vanishes at eps = 1", ok and at1)


def main():
    inv = o0()
    o1(inv)
    o2()
    o3()
    o4()
    o5()
    o6()
    o7(inv)
    o8()
    o9()
    P(f"\nTOTAL: {NCHK[0]} checks, {NTRUE[0]} True")
    P(f"runtime {time.time() - T0:.0f} s")


if __name__ == '__main__':
    main()
