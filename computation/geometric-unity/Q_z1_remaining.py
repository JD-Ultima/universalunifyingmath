"""GU round 13 (Q): the members of the family that round 10 (N) leaves open, in the Zorro metric Z1.

Formalization: formalization/geometric-unity/Q-z1-remaining-members.md.  Builds on N (N1, N2, N4, N6, N7, N9-N12), M (M1-M4,
M7, M8), L (L1, L3, L5, L8) and F (F4, F6, F7).  Imports N_z1_neighbourhood.py, M_remaining.py (and through them
L_nonradial.py, J_slaved_sector.py, F_jets.py, ...), all unmodified.  Exact arithmetic throughout: Python Fractions (also as
coefficients of N's truncated power series), sympy rationals.  Fixed seeds.  One core.  Every check is asserted.

Notation (as in the note): member p = (x, y, z, w), eps = +1 (14,0), -1 (7,7); a := w/x (M's c), A := a^2, tau := x t/kappa;
c := 4 cV (1 - 4 lambda), b := (1 - 5 lambda)/(cV (1 - 4 lambda)), beta_S := -5/(4 cV), B := 4 b c = 16 (1 - 5 lambda);
D1 = (1 + eps A) + 544 (1 - eps A) tau, N2 = numerator of M4's mu', L = (11 + eps A) + 2464 (1 - eps A) tau,
Q2 = 1 + 768 tau + (121856 - 147456 eps A) tau^2 (= det M/kappa^4 on u = 0).

Sections (labels = results of the note)
  Q0  inputs: the metric-free statements of round 9 (M's sections M0-M4, M6, M7, i.e. L's L0, the symmetries M1, the
      cone M2, the channel M3, mu' of M4, the x = +-y identities of M7, the (x,y) = 0 identities of M8), run silently.
  Q1  Lemma Q1(b): Z1 with (-cH, -cV) is -Z1(cH, cV) (2-jets off the section on a curved (1,3) base), Ric unchanged,
      signature (7,7) for both; c, b, beta_S change sign, B does not.
  Q2  Lemmas Q2, Q5: the N1(a) pairing c_V G_lambda(Z, sym (nabla_X Y)^V) = -(cH/2) Z(X,Y) used in Q2; Einstein constant
      -1/c in dimension 14 gives G = (6/c) g, c1 = 6x/c (= 15x/(2cV) at lambda = 1/5).
  Q3  Lemma Q6: the degree-1 system of N7 with mu_1 replaced by k = mu'(r0) (N's star_system with mu1 = k) on constant
      curvature and flat bases: solvable iff k c sigma^2 = 4b (4 k b c - 1), with the predicted solution.
  Q4  Lemma Q7: the persistence relations along Sigma; differentiating (C) with (B) gives (D); the linear cases.
  Q5  Proposition Q8: mu' and mu mu'' of the twisted cone from M's reduced tables; the exact certificate that the
      persistence conditions (C), (D) hold identically in tau only at the two linear cases, where rho = 0.
  Q6  Lemma Q9: H, q, q' in closed form; H and q have no common zero where det M != 0 (resultants).
  Q7  Proposition Q11: lambda = 1/5, second fibre derivatives of T~ = Ric + g/c at p (N's truncated series): T~ vanishes to
      order 1 in the HH and VV blocks and to order 2 over a flat base; the traced condition (O2) equals alpha'' |W|^2
      with alpha'' = -3 cV/(2 cH^2) (two Weyl tensors, three (cH, cV), with and without nabla R, homogeneity).
  Q8  Propositions Q12, Q13: the flat radial ODE for any mu; the x = +-y line, lambda = 23/100.
  Q9  Proposition Q15: existence of the flat-base constant solutions for every (14,0) member; det M at the roots.
Run (from computation/geometric-unity): python3 Q_z1_remaining.py > Q_z1_remaining.out 2>&1
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
import F_jets as FJ  # noqa: E402
import N_z1_neighbourhood as N  # noqa: E402
import M_remaining as M  # noqa: E402
import L_nonradial as L  # noqa: E402

T0 = time.time()
NCHK = [0]
n = 4
IDX = N.IDX
FC = N.FC
x, y, z, w, t, u, kap = L.x, L.y, L.z, L.w, L.t, L.u, L.kap
tau, A_, B_, rho = sp.symbols('tau A B rho')
SG = {1: '(14,0)', -1: '(7,7)'}


def P(*a):
    print(*a)
    sys.stdout.flush()


def check(label, ok):
    NCHK[0] += 1
    P(f"  [{NCHK[0]:3d}] {label}: {bool(ok)}")
    assert ok, label


def section(s):
    P(f"\n=== {s}   [{time.time() - T0:.0f} s]")


def polys(eps):
    D1 = (1 + eps * A_) + 544 * (1 - eps * A_) * tau
    N2 = (7 - 3 * eps * A_) + (7616 - 3008 * eps * A_ - 768 * A_ ** 2) * tau \
        + (2071552 - 3226624 * eps * A_ + 1155072 * A_ ** 2) * tau ** 2
    Lp = (11 + eps * A_) + 2464 * (1 - eps * A_) * tau
    Q2 = 1 + 768 * tau + (121856 - 147456 * eps * A_) * tau ** 2
    Nq = (7 + 19 * eps * A_) + (7616 - 2304 * eps * A_ - 4992 * A_ ** 2) * tau \
        + (2071552 - 3902464 * eps * A_ + 1830912 * A_ ** 2) * tau ** 2
    Nqp = (17 + 87 * eps * A_) + (18496 - 4288 * eps * A_ - 19968 * A_ ** 2) * tau \
        + (5030912 - 10732544 * eps * A_ + 5701632 * A_ ** 2) * tau ** 2
    return D1, N2, Lp, Q2, Nq, Nqp


def signature(Mx):
    """exact signature of a real symmetric rational matrix (Descartes' rule on the characteristic polynomial, all roots real)"""
    cp = sp.Matrix(Mx).charpoly(sp.Symbol('l')).all_coeffs()

    def changes(cs):
        cs = [c_ for c_ in cs if c_ != 0]
        return sum(1 for i in range(len(cs) - 1) if (cs[i] > 0) != (cs[i + 1] > 0))
    npos = changes(cp)
    neg = [c_ * (-1) ** (len(cp) - 1 - i) for i, c_ in enumerate(cp)]
    return npos, changes(neg)


# ====================================================================== Q0
def q0():
    section("Q0. Inputs: the metric-free statements of round 9 (M's sections M0-M4, M6, M7), run silently")
    buf = io.StringIO()
    n0 = M.NCHK[0]
    with contextlib.redirect_stdout(buf):
        M.m0()
        M.m1()
        M.m2()
        out = M.m3()
        M.m4(out)
        M.m6()
        M.m7()
    nl = M.NCHK[0] - n0
    lines = [ln for ln in buf.getvalue().splitlines() if ln.strip().startswith('[')]
    ok = len(lines) == nl and all(ln.rstrip().endswith(': True') for ln in lines)
    check(f"M.m0-m4, m6, m7 pass all their {nl} asserted checks (L's J1/J2 inputs; M1(a) rotation/boost and M1(b) swap of the"
          f" reduced system; the cone orbits and varpi; M3's channel of (x,0,0,w); M4's mu'; M7's x = +-y identities; M8's"
          f" theta-structure). None of them involves the metric g_aleph", ok)
    return out


# ====================================================================== Q1
def q1():
    section("Q1. Lemma Q1(b): the swap metric -g_aleph is Z1 with (-cH, -cV); Ric unchanged; (7,7) preserved")
    m = 14
    eta, lam = (1, -1, -1, -1), Fr(1, 2)
    rng = random.Random(111)
    gj = N.rand_general_base(rng, eta)
    base = FJ.base_geometry(gj, n, FC)
    y0 = [[Fr(eta[i]) * (i == j) + (Fr(rng.randint(-2, 2), 5) if i <= j else 0) for j in range(n)] for i in range(n)]
    for i in range(n):
        for j in range(i):
            y0[i][j] = y0[j][i]
    res = {}
    for cH, cV in ((Fr(2), Fr(3)), (Fr(-2), Fr(-3))):
        G, _, _, _ = FJ.total_metric(gj, base, n, y0, lam, cH, cV, 'Z1', FC)
        g0, dG, ddG = FJ.metric_point_data(G, m, FC)
        Ric, _, _ = FJ.ricci_at_point(G, m, FC)
        c, b, bS = N.consts(lam, cV)
        E = [Fr(0)] * n + [y0[i][j] for (i, j) in IDX]
        Ef = [sum(g0[a][bb] * E[bb] for bb in range(m)) for a in range(m)]
        okE = all(sum(Ric[a][bb] * E[bb] for bb in range(m)) == -Ef[a] / c for a in range(m))
        Gs, _, _, _ = FJ.total_metric(gj, base, n, [[Fr(eta[i]) * (i == j) for j in range(n)] for i in range(n)],
                                      lam, cH, cV, 'Z1', FC)
        gs = [[FJ.jval(Gs[a][bb], FC) for bb in range(m)] for a in range(m)]
        res[(cH, cV)] = dict(g0=g0, dG=dG, ddG=ddG, Ric=Ric, okE=okE, c=c, b=b, bS=bS, sig=signature(gs))
    r1, r2 = res[(Fr(2), Fr(3))], res[(Fr(-2), Fr(-3))]
    neg = lambda T: [neg(v) for v in T] if isinstance(T, list) else -T  # noqa: E731
    oklin = r2['g0'] == neg(r1['g0']) and r2['dG'] == neg(r1['dG']) and r2['ddG'] == neg(r1['ddG'])
    okR = r2['Ric'] == r1['Ric']
    okc = r2['c'] == -r1['c'] and r2['b'] == -r1['b'] and r2['bS'] == -r1['bS'] and 4 * r2['b'] * r2['c'] == 4 * r1['b'] * r1['c']
    check(f"curved (1,3) base in general coordinates, lambda = 1/2, y0 != g(x0): the 2-jet of Z1(-2,-3) is minus that of Z1(2,3)"
          f" ({oklin}); Ric equal ({okR}); Ric(E, .) = -E^flat/c holds for both with c -> -c ({r1['okE'] and r2['okE']}); c, b,"
          f" beta_S change sign, B = 4bc = {4 * r1['b'] * r1['c']} does not ({okc}); signature at the section point (n+, n-) ="
          f" {r1['sig']}, {r2['sig']}", oklin and okR and okc and r1['okE'] and r2['okE'] and r1['sig'] == (7, 7)
          and r2['sig'] == (7, 7))


# ====================================================================== Q2
def q2():
    section("Q2. Lemmas Q2, Q5: the N1(a) pairing used in Q2; the Einstein constant -1/c on int Z")
    lam = sp.Symbol('lambda')
    mu = lam / (1 - 4 * lam)
    rng = random.Random(121)
    ok = True
    for eta in ((1, 1, 1, 1), (1, -1, -1, -1), (1, 1, -1, -1)):
        Yv = sp.diag(*eta) + sp.Matrix(4, 4, lambda i, j: 0)
        Ym = sp.Matrix(4, 4, lambda i, j: sp.Rational(rng.randint(-2, 2), 7) if i <= j else 0)
        Yv = Yv + Ym + Ym.T - sp.diag(*[Ym[i, i] for i in range(4)])
        Yi = Yv.inv()
        Zm = sp.Matrix(4, 4, lambda i, j: sp.Rational(rng.randint(-3, 3)) if i <= j else 0)
        Zm = Zm + Zm.T - sp.diag(*[Zm[i, i] for i in range(4)])
        X = sp.Matrix([rng.randint(-2, 2) for _ in range(4)])
        Yw = sp.Matrix([rng.randint(-2, 2) for _ in range(4)])
        Xf, Yf = Yv * X, Yv * Yw
        sym = (Xf * Yf.T + Yf * Xf.T) / 2 + mu * (X.T * Yv * Yw)[0] * Yv
        G = (Yi * Zm * Yi * sym).trace() - lam * (Yi * Zm).trace() * (Yi * sym).trace()
        ok &= sp.simplify(G - (X.T * Zm * Yw)[0]) == 0
    check("for symmetric y (three signatures), symmetric Z, vectors X, Y and symbolic lambda: G_lambda(Z, X^flat . Y^flat +"
          " mu y(X,Y) y) = Z(X,Y), mu = lambda/(1 - 4 lambda); with N1(a), cV G_lambda(Z, sym (nabla_X Y)^V) = -(cH/2) Z(X,Y)", ok)
    c_, cV_, x_ = sp.symbols('c c_V x')
    Ric_coef = -1 / c_
    s = 14 * Ric_coef
    Gc = sp.simplify(Ric_coef - s / 2)
    lam5 = sp.Rational(1, 5)
    c5 = 4 * cV_ * (1 - 4 * lam5)
    check(f"dimension 14, Ric = -(1/c) g (N1(d)): G = {Gc} g, so -xG + c1 g = 0 gives c1 = 6x/c; at lambda = 1/5,"
          f" 6x/c = {sp.simplify(6 * x_ / c5)}", sp.simplify(Gc - 6 / c_) == 0 and sp.simplify(6 * x_ / c5 - 15 * x_ / (2 * cV_)) == 0)


# ====================================================================== Q3
def q3():
    section("Q3. Lemma Q6: the degree-1 system (*) with mu_1 -> k = mu'(r0) (N's star_system(D, k)), rhs / mu(r0)")
    cH, cV = Fr(2), Fr(3)
    oks = []

    def run(eta, lam, sig, ks, expect):
        c, b, bS = N.consts(lam, cV)
        K = sig * cH / 12
        R = N.radd((K, N.R1(eta)))
        D = N.point_data(N.normal_base(eta, R, None), lam, cH, cV)
        sc = N.scal_of(R, eta) / cH
        out = {}
        for kk in ks:
            rows, rhs = N.star_system(D, kk)
            r = N.rank_test(rows, rhs)
            good = None
            if r[0] == r[1]:
                x_, ok_, _ = FJ.solve_fit([[rw[j] for rw in rows] for j in range(28)], rhs, FC)
                E = N.E_vec(eta)
                good = ok_ and x_[:14] == [(2 * b - sc / 2) * e for e in E] and x_[14:] == [-kk * bS * (2 * b - sc / 2) * e for e in E]
            out[str(kk)] = (r, good)
        okk = all((out[str(kk)][0][0] == out[str(kk)][0][1]) == e_ and (not e_ or out[str(kk)][1]) for kk, e_ in zip(ks, expect))
        okk &= sc == sig and all(v[0][0] == 28 for v in out.values())
        return okk, out, (c, b)
    lam = Fr(1, 7)
    c, b, _ = N.consts(lam, cV)
    kp = -4 * b / (c * (1 - 16 * b * b))
    ok, out, _ = run((1, 1, 1, 1), lam, Fr(1), [kp, Fr(0), Fr(-7, 18)], [True, False, False])
    check(f"(4,0), lambda = 1/7, constant curvature sigma_g = 1: predicted k = -4b/(c(sigma^2 - 16 b^2)) = {kp}; (rank, augmented"
          f" rank, predicted solution v' = (2b - sigma/2) E, n' = -k beta_S v') = {out}: solvable exactly at the predicted k", ok)
    ok, out, _ = run((1, 1, 1, 1), lam, 4 * b, [Fr(0), Fr(-7, 18), Fr(5)], [False, False, False])
    check(f"(4,0), lambda = 1/7, sigma_g = 4b = {4 * b} (where 2b - sigma/2 = 0): not solvable for any sampled k, as (c) requires"
          f" (k c (sigma^2 - 16 b^2) = 0 != -4b): {out}", ok)
    kf = 1 / (16 * (1 - 5 * lam))
    ok, out, _ = run((1, 1, 1, 1), lam, Fr(0), [kf, Fr(-7, 18)], [True, False])
    check(f"(4,0), lambda = 1/7, flat base: solvable iff 4 k b c = 1, k = 1/B = {kf}: {out}", ok)
    ok, out, _ = run((1, 1, 1, 1), Fr(1, 5), Fr(0), [Fr(0), Fr(-7, 18), Fr(1)], [True, True, True])
    ok2, out2, _ = run((1, 1, 1, 1), Fr(1, 5), Fr(1), [Fr(0), Fr(-7, 18)], [True, False])
    check(f"(4,0), lambda = 1/5 (b = 0): flat base solvable for every sampled k {out}; sigma_g = 1 only at k = 0 {out2}", ok and ok2)
    lam = Fr(1, 2)
    c, b, _ = N.consts(lam, cV)
    kp = -4 * b / (c * (Fr(1, 4) - 16 * b * b))
    ok, out, _ = run((1, -1, -1, -1), lam, Fr(1, 2), [kp, Fr(-5, 12)], [True, False])
    check(f"(1,3), lambda = 1/2 ((7,7) side), sigma_g = 1/2: predicted k = {kp}: {out}", ok)


# ====================================================================== Q4
def q4():
    section("Q4. Lemma Q7: persistence along Sigma, theta -> q_theta = (x0, e^theta g(x0)) an integral curve of E")
    th, sg, b, c, k = sp.symbols('theta sigma b c k')
    rh = sp.Function('rho')(th)
    m0, m1, m2 = sp.symbols('mu mu1 mu2')        # mu, mu', mu'' at r(theta)
    Bs = 4 * b * c
    # (C): c mu' (e^{-2 theta} sigma^2 - 16 b^2) = -4b with e^{-theta} sigma = 4 b rho  <=>  B mu' (1 - rho^2) = 1
    C1 = c * m1 * (16 * b ** 2 * rh ** 2 - 16 * b ** 2) + 4 * b
    okC = sp.expand(C1 + 4 * b * (Bs * m1 * (1 - rh ** 2) - 1)) == 0
    # (B): dr/dtheta = c v_E = c mu (2b - e^{-theta} sigma/2) = (B/2) mu (1 - rho)
    okB = sp.expand(c * m0 * (2 * b - 4 * b * rh / 2) - Bs / 2 * m0 * (1 - rh)) == 0
    # d/dtheta [B mu'(r)(1 - rho^2)] with rho' = -rho, r' = (B/2) mu (1 - rho), and (1 - rho^2) = 1/(B mu')
    rr = sp.Symbol('rr')
    deriv = Bs * m2 * (Bs / 2 * m0 * (1 - rr)) * (1 - rr ** 2) + Bs * m1 * (-2 * rr * (-rr))
    deriv = deriv.subs(rr ** 2, 1 - 1 / (Bs * m1))
    Phi = m0 * m2 / m1 ** 2
    target = Bs * m1 / 2 * (Phi * (1 - rr) + 4 * (1 - 1 / (Bs * m1)))
    okD = sp.simplify(sp.expand(deriv) - sp.expand(target)) == 0
    check("with e^(-theta) sigma_g = 4 b rho: (C) <=> B mu'(1 - rho^2) = 1 and (B) <=> dr/dtheta = (B/2) mu (1 - rho)"
          f" ({okC}, {okB}); differentiating (C) along Sigma with rho' = -rho and (B) gives (B mu'/2)[Phi (1 - rho) + 4 rho^2],"
          f" Phi = mu mu''/mu'^2, so (D): Phi (1 - rho) + 4 rho^2 = 0 ({okD})", okC and okB and okD)
    # single point (c) for k = mu_1 reproduces N7(c)
    l_, cv_ = sp.symbols('lambda c_V')
    cs = 4 * cv_ * (1 - 4 * l_)
    bs = (1 - 5 * l_) / (cv_ * (1 - 4 * l_))
    m1n = sp.Rational(-7, 18)
    s2 = sp.solve(k * c * sg ** 2 - 4 * b * (4 * k * b * c - 1), sg ** 2)[0]
    okN7 = sp.simplify(s2.subs({k: m1n, b: bs, c: cs}) - 10 * (1 - 5 * l_) * (13 - 56 * l_) / (7 * cv_ ** 2 * (1 - 4 * l_) ** 2)) == 0
    lams = {}
    for m1v in (sp.Rational(-7, 18), sp.Rational(-5, 12), sp.Rational(-3, 289)):
        lams[str(m1v)] = sp.solve(4 * m1v * 4 * (1 - 5 * l_) - 1, l_)[0]
    okl = lams == {'-7/18': sp.Rational(13, 56), '-5/12': sp.Rational(23, 100), '-3/289': sp.Rational(337, 240)}
    check(f"k c sigma^2 = 4b(4kbc - 1) at k = mu_1 = -7/18 is N7(c) ({okN7}); linear mu = mu_1 r (Phi = 0, so (c) forces rho = 0):"
          f" lambda in {{1/5}} u {{4 mu_1 b c = 1}} = {{1/5, {lams}}} for mu_1 = -7/18 (Shiab#, W_A: N7(d)), -5/12 (x = +-y"
          f" line, M7(b)), -3/289 (a^2 = 11/6); 337/240 > 1/4 is not on the (14,0) side, where a^2 = 11/6 lives", okN7 and okl)


# ====================================================================== Q5
def q5(out):
    section("Q5. Proposition Q8: the twisted cone (x,0,0,ax), a != 0: mu', mu mu'' and the persistence certificate")
    ok_forms = True
    store = {}
    for eps in (1, -1):
        r = out[eps]
        D0 = r['D'].subs(u, 0)
        H0 = sp.cancel(r['HN'][0, 0].subs(u, 0) / D0)
        q0 = sp.cancel(r['QN'][('E1', 't', 't')].subs(u, 0) / D0 ** 2)
        mup = sp.cancel(x * (q0 - sp.diff(H0, t)) / H0 ** 2)
        mumupp = sp.cancel(x * sp.diff(mup, t) / H0)          # mu mu'' = mu (d mu'/dt)/F' = x (d mu'/dt)/H
        sub = {t: tau * kap / x, w: sp.sqrt(A_) * x}
        m1 = sp.cancel(mup.subs(sub))
        mm = sp.cancel(mumupp.subs(sub))
        D1, N2, Lp, Q2, _, _ = polys(eps)
        ok_forms &= sp.cancel(m1 + N2 / (18 * D1 ** 2)) == 0
        ok_forms &= sp.cancel(mm - A_ * (6 * A_ - 11 * eps) * Lp * Q2 / (81 * D1 ** 4)) == 0
        store[eps] = (m1, mm)
    check("from M's reduced tables on u = 0 (r = F(t), F''/F' = q/H, mu = x F'/H): mu' = -N2/(18 D1^2) (M4) and"
          " mu mu'' = x (d mu'/dt)/H = A (6A - 11 eps) L Q2/(81 D1^4), both functions of (tau, A) only; hence"
          " Phi = mu mu''/mu'^2 = 4 A (6A - 11 eps) L Q2/N2^2, and mu' is constant on an interval iff A (6A - 11 eps) = 0", ok_forms)
    for eps in (1, -1):
        t1 = time.time()
        D1, N2, Lp, Q2, _, _ = polys(eps)
        m1 = -N2 / (18 * D1 ** 2)
        Phi = 4 * A_ * (6 * A_ - 11 * eps) * Lp * Q2 / N2 ** 2
        X = 1 - 1 / (B_ * m1)                                   # rho^2 from (C)
        Rt = sp.together((Phi + 4 * X) ** 2 - X * Phi ** 2)      # (D) with rho = 1 + 4X/Phi, squared against (C)
        num = sp.expand(sp.numer(Rt))
        degt = sp.Poly(num, tau).degree()
        es = [sp.expand(num.subs(tau, v)) for v in (0, 1, -1, 2)]
        g = sp.factor(sp.gcd(sp.resultant(es[0], es[1], B_), sp.resultant(es[0], es[2], B_)))
        okg = sp.simplify(g / (A_ ** 4 * (3 * A_ - 7 * eps) ** 2 * (6 * A_ - 11 * eps) ** 4)).is_number
        cands = sorted(set(sp.Poly(g, A_).all_roots()))
        found, info = [], []
        for Av in cands:
            gB = sp.Poly(es[0].subs(A_, Av), B_)
            for e_ in es[1:]:
                gB = sp.gcd(gB, sp.Poly(e_.subs(A_, Av), B_))
            Bs_ = sorted(set(gB.all_roots())) if gB.degree() > 0 else []
            info.append((Av, gB.degree(), Bs_))
            for Bv in Bs_:
                if sp.expand(num.subs({A_: Av, B_: Bv})) == 0:
                    Xz = sp.cancel(X.subs({A_: Av, B_: Bv})) == 0
                    mconst = sp.cancel(Bv * m1.subs(A_, Av) - 1) == 0
                    found.append((Av, Bv, Xz, mconst))
        want = {1: [(0, sp.Rational(-18, 7)), (sp.Rational(11, 6), sp.Rational(-289, 3))],
                -1: [(sp.Rational(-11, 6), sp.Rational(-289, 3)), (0, sp.Rational(-18, 7))]}[eps]
        okc = [(f_[0], f_[1]) for f_ in found] == want and all(f_[2] and f_[3] for f_ in found)
        okc &= [i_[2] for i_ in info if i_[0] == sp.Rational(7, 3) * eps] == [[]]
        store[eps] = found
        check(f"{SG[eps]}: R~ := (Phi + 4X)^2 - X Phi^2, X = rho^2 = 1 - 1/(B mu'), has numerator of degree {degt} in tau; the gcd of"
              f" Res_B(R~(0), R~(1)) and Res_B(R~(0), R~(-1)) is {g}; at each root A the common B-roots of R~(0), R~(1), R~(-1),"
              f" R~(2) are (A, B-degree of their gcd, roots) = {info}; R~ vanishes identically in tau exactly at (A, B) ="
              f" {[(f_[0], f_[1]) for f_ in found]}, where X = 0 (rho = 0) and B mu' = 1 (the linear cases: lambda = 13/56 for"
              f" a = 0, 337/240 for a^2 = 11 eps/6)  [{time.time() - t1:.0f} s]", okg and okc)
    pos = [f_ for eps in (1, -1) for f_ in store[eps] if f_[0] > 0]
    lam337 = sp.solve(16 * (1 - 5 * sp.Symbol('l')) + sp.Rational(289, 3), sp.Symbol('l'))
    check(f"hence for real a != 0 (A > 0) the only case is {[(f_[0], f_[1]) for f_ in pos]} in (14,0), i.e. B = 16(1 - 5 lambda)"
          f" = -289/3, lambda = {lam337}, which is > 1/4 and so not on the (14,0) side: neither Q7(b) (mu' = 1/B constant) nor"
          f" Q7(c) (R~ = 0 with rho != 0) holds on any interval for an admissible lambda != 1/5",
          len(pos) == 1 and pos[0][0] == sp.Rational(11, 6) and lam337 == [sp.Rational(337, 240)] and lam337[0] > sp.Rational(1, 4))


# ====================================================================== Q6
def q6(out):
    section("Q6. Lemma Q9: H, q, q' of the twisted cone in closed form; zeros of H")
    ok = True
    for eps in (1, -1):
        r = out[eps]
        D0 = r['D'].subs(u, 0)
        H0 = sp.cancel(r['HN'][0, 0].subs(u, 0) / D0)
        q0 = sp.cancel(r['QN'][('E1', 't', 't')].subs(u, 0) / D0 ** 2)
        qp0 = sp.cancel(r['QpN'][('E1', 't', 't')].subs(u, 0) / D0 ** 2)
        sub = {t: tau * kap / x, w: sp.sqrt(A_) * x}
        D1, N2, Lp, Q2, Nq, Nqp = polys(eps)
        ok &= sp.cancel(H0.subs(sub) - 576 * x ** 2 * D1 / (kap * Q2)) == 0
        ok &= sp.cancel(q0.subs(sub) + 36864 * x ** 3 * Nq / (kap ** 2 * Q2 ** 2)) == 0
        ok &= sp.cancel(qp0.subs(sub) - 4608 * x ** 3 * Nqp / (kap ** 2 * Q2 ** 2)) == 0
        ok &= sp.expand(D0.subs(sub) - kap ** 4 * Q2) == 0
    check("on u = 0: H = 576 x^2 D1/(kappa Q2), q = -36864 x^3 Nq/(kappa^2 Q2^2), q' = 4608 x^3 Nq'/(kappa^2 Q2^2),"
          " det M = kappa^4 Q2 (both signatures)", ok)
    res = {}
    okr = True
    for eps in (1, -1):
        D1, N2, Lp, Q2, Nq, Nqp = polys(eps)
        rq = sp.factor(sp.resultant(D1, Nq, tau))
        rQ = sp.factor(sp.resultant(D1, Q2, tau))
        res[eps] = (rq, rQ)
        okr &= sp.cancel(rq / (A_ * (A_ - eps) * (6 * A_ - 11 * eps) ** 2)).is_number
        okr &= sp.cancel(rQ / (A_ * (6 * A_ - 11 * eps) ** 2)).is_number
    check(f"Res_tau(D1, Nq) and Res_tau(D1, Q2): {res}. For real a != 0: in (7,7) none vanishes; in (14,0) they vanish only at"
          f" A = 1 (where D1 = 2 is constant: H has no zero) and A = 11/6 (where the zero of H is a zero of det M). So H and q"
          f" have no common zero where det M != 0", okr)


# ====================================================================== Q7
def fibre_series_g(pb, lam, cH, cV, K):
    """as N.fibre_series, also returning the metric g_aleph at (0, g(0) + h) as truncated series"""
    gj = N.base_jets_at(pb, FC, [0] * n)
    base = FJ.base_geometry(gj, n, FC)
    A = N.Alg(10, K)
    ctx = N.FCtx(A)
    y0 = N.fibre_point(ctx, pb[0], A)
    G, m, _, Psi0 = FJ.total_metric(gj, base, n, y0, ctx.c(lam), ctx.c(cH), ctx.c(cV), 'Z1', ctx)
    assert all(v == 0 for row in Psi0 for v in row)
    Ric, _, _ = FJ.ricci_at_point(G, m, ctx)
    gs = [[FJ.jval(G[a][b], ctx) for b in range(m)] for a in range(m)]
    return A, ctx, Ric, gs


def o2_pairing(W, cH, cV, eta=(1, 1, 1, 1), cubic=None):
    lam = Fr(1, 5)
    m = 14
    pb = N.normal_base(eta, W, None)
    if cubic is not None:
        pb = (pb[0], pb[1], cubic)
    A, ctx, Ric, gs = fibre_series_g(pb, lam, cH, cV, 2)
    c = 4 * cV * (1 - 4 * lam)
    Tt = [[Ric[a][b] + gs[a][b] * (1 / c) for b in range(m)] for a in range(m)]
    blk = lambda a, b: (a < n) == (b < n)  # noqa: E731
    o0 = all(not Tt[a][b].part(0) for a in range(m) for b in range(m))
    o1 = all(not Tt[a][b].part(1) for a in range(m) for b in range(m) if blk(a, b))
    o2zero = all(not Tt[a][b].part(2) for a in range(m) for b in range(m))
    E = N.E_vec(eta)
    V0 = [sp.Matrix(N.symmat(v)) for v in N.traceless_V0(eta)]
    et = sp.diag(*eta)
    eti = et.inv()
    G = sp.Matrix(9, 9, lambda i, j: (eti * V0[i] * eti * V0[j]).trace())
    Gi = G.inv()
    dual = [sum((Gi[i, j] * V0[j] for j in range(9)), sp.zeros(4, 4)) for i in range(9)]
    fr = L.fr

    def co(U):
        return [fr(U[i, j]) for (i, j) in IDX]

    def d2(a, b, cu, cw):
        e = Tt[a][b]
        tot = Fr(0)
        for k in range(10):
            if not cu[k]:
                continue
            for l_ in range(10):
                if not cw[l_]:
                    continue
                mm = [0] * 10
                mm[k] += 1
                mm[l_] += 1
                tot += e.coef(tuple(mm)) * (2 if k == l_ else 1) * cu[k] * cw[l_]
        return tot
    S2 = Fr(0)
    for k in range(9):
        cu, cw = co(V0[k]), co(dual[k])
        trH = sum(Fr(eta[X]) * d2(X, X, cu, cw) for X in range(n))
        EE = sum(E[a] * E[b] * d2(a, b, cu, cw) for a in range(m) if E[a] for b in range(m) if E[b])
        S2 += trH / (4 * cH) - EE / c
    # the traceless HH part and its only invariant pairing, with P(U,W) = (U eta^-1 W + W eta^-1 U)_0
    S0 = Fr(0)
    Ef = [[fr(v) for v in row] for row in eti.tolist()]
    for k in range(9):
        for l_ in range(9):
            cu, cw = co(V0[k]), co(V0[l_])
            Qm = [[d2(X, Y, cu, cw) for Y in range(n)] for X in range(n)]
            trq = sum(Ef[X][X] * Qm[X][X] for X in range(n))
            Qm = [[Qm[X][Y] - trq / 4 * Fr(eta[X]) * (X == Y) for Y in range(n)] for X in range(n)]
            Pm = (dual[k] * eti * dual[l_] + dual[l_] * eti * dual[k]) / 2
            Pm = Pm - (eti * Pm).trace() / 4 * et
            Pf = [[fr(Pm[X, Y]) for Y in range(n)] for X in range(n)]
            S0 += sum(Ef[X][X] * Ef[Y][Y] * Qm[X][Y] * Pf[X][Y] for X in range(n) for Y in range(n))
    return dict(S2=S2, S0=S0, o0=o0, o1=o1, o2zero=o2zero)


def q7():
    section("Q7. Proposition Q11 (lambda = 1/5, zeros of H): the traced K-free condition (O2) = alpha'' |W|^2")
    eta = (1, 1, 1, 1)
    zeroW = {k_: Fr(0) for k_ in itertools.product(range(4), repeat=4)}
    alphas = {}
    for (cH, cV), seeds in (((Fr(2), Fr(3)), (71, 72)), ((Fr(1), Fr(1)), (73,)), ((Fr(3), Fr(1, 2)), (74,))):
        t1 = time.time()
        fl = o2_pairing(zeroW, cH, cV)
        check(f"(cH, cV) = ({cH}, {cV}), flat base: T~ = Ric + g/c vanishes to order 2 in the fibre (Einstein flat model, N2)"
              f" ({fl['o0'] and fl['o2zero']} ; all orders <= 2), so (O2) = 0 there ({fl['S2']})  [{time.time() - t1:.0f} s]",
              fl['o0'] and fl['o1'] and fl['o2zero'] and fl['S2'] == 0)
        for seed in seeds:
            t1 = time.time()
            Wt = N.weyl_of(N.rand_R(random.Random(seed), eta), eta)
            W2 = N.norm2(Wt, eta)
            r1 = o2_pairing(Wt, cH, cV)
            r2 = o2_pairing({k_: 2 * v for k_, v in Wt.items()}, cH, cV)
            al = r1['S2'] / W2
            alphas.setdefault((cH, cV), []).append(al)
            check(f"(cH, cV) = ({cH}, {cV}), random Weyl W (seed {seed}, |W|^2 = {W2}): T~ = 0 at p and its HH, VV blocks vanish to"
                  f" order 1 along the fibre ({r1['o0'] and r1['o1']}); sum_k [tr_H d_k d^k T~_HH/(4 cH) - d_k d^k T~(E,E)/c]"
                  f" = {r1['S2']} = alpha'' |W|^2 with alpha'' = {al}; for 2W it is 4 times as large ({r2['S2'] == 4 * r1['S2']}:"
                  f" no linear part); the traceless HH pairing with (U eta^-1 W + W eta^-1 U)_0 is {r1['S0']}"
                  f"  [{time.time() - t1:.0f} s]", r1['o0'] and r1['o1'] and r2['S2'] == 4 * r1['S2'] and al != 0 and r1['S0'] == 0)
    t1 = time.time()
    Wt = N.weyl_of(N.rand_R(random.Random(71), eta), eta)
    rc = o2_pairing(Wt, Fr(2), Fr(3), cubic=N.rand_cubic(random.Random(175)))
    W2 = N.norm2(Wt, eta)
    check(f"(cH, cV) = (2, 3), Weyl W (seed 71) plus random cubic terms (nabla R(x0) != 0): T~(p) = 0 and its HH, VV blocks"
          f" vanish to order 1 ({rc['o0'] and rc['o1']}); (O2) = {rc['S2']}, the same as without cubic terms"
          f" ({rc['S2'] == alphas[(Fr(2), Fr(3))][0] * W2}: no nabla R dependence, as the x -> -x parity requires)"
          f"  [{time.time() - t1:.0f} s]", rc['o0'] and rc['o1'] and rc['S2'] == alphas[(Fr(2), Fr(3))][0] * W2)
    fit = all(a_ == Fr(-3) * cv / (2 * ch ** 2) for (ch, cv), v in alphas.items() for a_ in v)
    check(f"alpha'' is the same for the two Weyl tensors at (2, 3) ({[str(a_) for a_ in alphas[(Fr(2), Fr(3))]]}; Schur: an O(4)-"
          f"invariant quadratic form on Weyl tensors is a multiple of |W|^2), and alpha'' = -3 cV/(2 cH^2) for all three"
          f" parameter sets ({fit}); never 0, so for a definite base (O2) forces W = 0",
          alphas[(Fr(2), Fr(3))][0] == alphas[(Fr(2), Fr(3))][1] and fit)


# ====================================================================== Q8
def q8():
    section("Q8. Propositions Q12, Q13: flat base for any mu; the x = +-y line (lambda = 23/100)")
    b, c, mu0, mu1 = sp.symbols('b c mu mu1')
    # (*) on the flat model, r = f(s): V_0 block nu = -mu beta_S; H block f'/(2c) = mu b; E block f'' = mu b c.
    # With f' = 2bc mu(f): f'' = 2bc mu'(f) f' = 4 b^2 c^2 mu' mu
    fpp = 2 * b * c * mu1 * (2 * b * c * mu0)
    cond = sp.factor(fpp - b * c * mu0)
    ok = sp.expand(cond - b * c * mu0 * (4 * b * c * mu1 - 1)) == 0
    check(f"flat base, r = f(s): f' = 2bc mu(f) and f'' = bc mu(f) are compatible iff {cond} = 0, i.e. b mu (4 b c mu'(f) - 1) = 0:"
          f" at lambda = 1/5 (b = 0) f' = 0; for b != 0 mu' = 1/(4bc) on the interval f(s)", ok)
    l_ = sp.Symbol('lambda')
    m1 = sp.Rational(-5, 12)
    lam = sp.solve(4 * m1 * 4 * (1 - 5 * l_) - 1, l_)
    check(f"x = +-y, ell != 0, (7,7): on the invariant line sigma = -kappa/208 the equation is linear with mu_1 = -5/12 (M7(b),"
          f" checked in Q0); Q7(d): lambda in {{1/5, {lam[0]}}}, both < 1/4, so none on the (7,7) side",
          lam == [sp.Rational(23, 100)] and sp.Rational(23, 100) < sp.Rational(1, 4))


# ====================================================================== Q9
def q9():
    section("Q9. Proposition Q15: flat-base constant solutions for every (14,0) member (lambda = 1/5, T_e = 0)")
    okD = True
    for eps in (1, -1):
        S, Mm, Lt, Bt = M.subs_sim(L.tables(eps), {y: 0, u: 0})
        D = sp.expand(Mm.det())
        okD &= sp.expand(D - kap ** 2 * (kap ** 2 + 768 * kap * t * x + 121856 * t ** 2 * x ** 2 - eps * 147456 * w ** 2 * t ** 2)) == 0
        okD &= not D.has(z)
        c1, c13 = M.cfun(eps, {y: 0})
        okD &= sp.expand(c13 - u * (kap + 208 * x * t)) == 0 and sp.expand(c1 - kap * t - x * (312 * t ** 2 + eps * 104 * u ** 2)) == 0
    check("members (x, 0, z, w) (every member with (x,y) != 0 is of this form up to M1(a)): det M on u = 0 is"
          " kappa^2 (kappa^2 + 768 kappa t x + 121856 t^2 x^2 - 147456 eps w^2 t^2), independent of z; c1 = kappa t + x(312 t^2"
          " +- 104 u^2), c13 = u (kappa + 208 x t)", okD)
    xp, s, r, cv, k = sp.symbols('xp s r c_V kappa', positive=True)
    I = sp.I
    sub = {x: -I * xp, t: I * s, u: I * r, kap: k}
    c1, c13 = M.cfun(1, {y: 0})
    e1 = sp.expand(xp * (c1.subs(sub) - 15 * (-I * xp) / (2 * cv)) / I)
    e13 = sp.expand(c13.subs(sub) / I)
    T, R = sp.symbols('T R')
    okE = sp.expand(e1 - (312 * T ** 2 + k * T + 104 * R ** 2 + sp.Rational(15, 2) * xp ** 2 / cv).subs({T: xp * s, R: xp * r})) == 0
    okE &= sp.expand(e13 - r * (k + 208 * xp * s)) == 0
    disc = sp.factor(sp.discriminant(312 * T ** 2 + k * T + sp.Rational(15, 2) * xp ** 2 / cv, T))
    br = sp.simplify((312 * T ** 2 + k * T).subs(T, -k / 208))
    check(f"(14,0) real variables (x = -i x', t = i s, u = i r, T = x' s, R = x' r): c1 = 15x/(2cV), c13 = 0 is"
          f" r (kappa + 208 T) = 0 and 312 T^2 + kappa T + 104 R^2 + 15 x'^2/(2cV) = 0 ({okE}); on r = 0 the discriminant is"
          f" {disc} (>= 0 iff cV kappa^2 >= 9360 x'^2); on T = -kappa/208 the left side is {br} + 104 R^2 + 15x'^2/(2cV) > 0"
          f" (cV > 0): no solution", okE and sp.simplify(disc - (k ** 2 - 9360 * xp ** 2 / cv)) == 0 and br == k ** 2 / 416)
    Aa = sp.Symbol('A', nonnegative=True)
    p1 = 312 * T ** 2 + k * T + sp.Rational(15, 2) * xp ** 2 / cv
    dM = k ** 2 + 768 * k * T + (121856 - 147456 * Aa) * T ** 2
    res = sp.factor(sp.resultant(p1, dM, T))
    om = sp.Symbol('omega', positive=True)
    rom = sp.factor(sp.numer(sp.together(res.subs(cv, 9360 * xp ** 2 / (om * k ** 2)))))
    Afree = sp.Symbol('A')
    prop = sp.solve(sp.Eq(sp.Rational(312, 1) * 768, 121856 - 147456 * Afree), Afree)
    lead = sp.factor(sp.Poly(rom, om).LC())
    cst = sp.factor(sp.Poly(rom, om).nth(0))
    ratio = sp.cancel(cst / (-(3504384 * Aa + 485199)))
    at0 = sp.solve(rom.subs(Aa, 0), om)
    check(f"det M at the roots, a := w'/x' (A = a^2): Res_T(312T^2 + kappa T + 15x'^2/(2cV), det M/kappa^2) in"
          f" omega = 9360 x'^2/(cV kappa^2) has numerator {rom} (degree <= 2 in omega, leading coefficient {lead}, constant term"
          f" {cst} < 0), so for each a at most two values of cV kappa^2/x'^2 put a root on det M = 0; both roots at once would"
          f" need the two quadratics proportional, i.e. A = {prop} < 0. At a = 0: omega = {at0} (cV kappa^2 = 277440/29, N12(c))",
          len(prop) == 1 and prop[0] < 0 and ratio.is_positive and sp.Poly(rom, om).degree() == 2
          and at0 == [sp.Rational(1131, 1156)])


def main():
    out = q0()
    q1()
    q2()
    q3()
    q4()
    q5(out)
    q6(out)
    q7()
    q8()
    q9()
    P(f"\n{NCHK[0]} checks, all True (each asserted).")
    P(f"runtime {time.time() - T0:.0f} s")


if __name__ == '__main__':
    main()
