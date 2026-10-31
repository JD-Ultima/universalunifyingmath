"""Referee cross-check of GU round 13 (formalization/geometric-unity/Q-z1-remaining-members.md), REVIEW10-GU.

Independent construction.  Nothing is imported from Q_z1_remaining.py, N_z1_neighbourhood.py, M_remaining.py,
L_nonradial.py, J_slaved_sector.py, I_elliptic_reduction.py or F_jets.py.  Repository imports (all referee scripts,
unmodified):
  - xcheck_round8.py (and through it xcheck_round6.py, xcheck_round5.py): the J1/J2 tables S, M, l, b, c1, c13 from the
    round-5/6 referee Clifford construction (own Clifford algebra, Shiab templates from their definitions);
  - xcheck_round10.py: the REVIEW8 multi-dual numbers MD, truncated jets, the textbook Ricci routine (controls S^3, H^3 in
    its own X0), the Zorro metric Z1 = c_H y(dx,dx) + c_V G_lambda|_y(Dy, Dy) as 2-jets at an arbitrary point of
    Y = Met(R^4), and its random-curvature / Weyl-base helpers.
Written here: the assembly of H, q, q' of the twisted cone from the tables; mu', mu mu'' (by a generic reparametrisation
F with F'' = (q/H) F'); an own persistence certificate (Groebner basis of ALL tau-coefficients of the reduced numerator);
the resultants of Q9 and Q15; the algebra of Q6, Q7, Q9(a), Q11 (O2), Q12, Q13; the (O2) contraction for alpha'' from exact
second fibre derivatives of T~ = Ric + g/c (nested multi-dual shift along V_0, and a second route by single-direction
second derivatives over a G_0-orthogonal basis); Q1(b) by brute force off the section; an own exact inertia routine.
Sections
  X0  inputs: the referee tables (xcheck_round8 X0 -> xcheck_round6 X0-X2)
  X1  Lemma Q1(b): Z1(-cH, -cV) = -Z1(cH, cV) off the section, Ric and G unchanged, N1(d), signature (7,7)
  X2  Lemmas Q2, Q5 and Propositions Q3, Q4: the G_lambda pairing; the F' of Q3; c1 - s c13 for the null members;
      M3's channel for (x, 0, 0, w) and c1 on y = 0
  X3  Proposition Q8: H, q, q', det M, mu', mu mu'' from the tables; where mu' is constant; the persistence certificate
  X4  Lemma Q9: (a) the block algebra; (b) the resultants
  X5  Lemmas Q6, Q7, Propositions Q11 (step 4), Q12: the component algebra (sympy)
  X6  Proposition Q11: alpha'' from the exact second fibre derivatives of T~ at lambda = 1/5 (base (4,0))
  X7  Proposition Q13: mu_1 = -5/12 on the x = +-y line and lambda in {1/5, 23/100}
  X8  Proposition Q15: existence for every (14,0) member without the rotation M1(a); det M at the roots
Exact arithmetic throughout (Fractions, sympy rationals).  Fixed seeds random.Random(20261131..20261136).  One core.
Every check is asserted and printed.  Checks marked "referee"/"observation" test statements of the review (REVIEW10-GU
R13-k); for those, True means the referee's claim was confirmed.
Run (from computation/geometric-unity): python3 xcheck_round13.py > xcheck_round13.out 2>&1
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
import xcheck_round8 as X8  # noqa: E402
import xcheck_round10 as X10  # noqa: E402

T0 = time.time()
RES = []
n = 4
x, y, z, w, kap, t, u = X8.x, X8.y, X8.z, X8.w, X8.kap, X8.t, X8.u
KS = X8.KS
tau, A_, B_ = sp.symbols('tau A B')
LAB = {1: '(14,0)', -1: '(7,7)'}
MD = X10.MD


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


def polys(ep):
    """the note's polynomials (Q.0 notation), A = a^2"""
    D1 = (1 + ep * A_) + 544 * (1 - ep * A_) * tau
    N2 = (7 - 3 * ep * A_) + (7616 - 3008 * ep * A_ - 768 * A_ ** 2) * tau \
        + (2071552 - 3226624 * ep * A_ + 1155072 * A_ ** 2) * tau ** 2
    L1 = (11 + ep * A_) + 2464 * (1 - ep * A_) * tau
    Q2 = 1 + 768 * tau + (121856 - 147456 * ep * A_) * tau ** 2
    Nq = (7 + 19 * ep * A_) + (7616 - 2304 * ep * A_ - 4992 * A_ ** 2) * tau \
        + (2071552 - 3902464 * ep * A_ + 1830912 * A_ ** 2) * tau ** 2
    return D1, N2, L1, Q2, Nq


def inertia(Mx):
    """(n+, n-, n0) of a symmetric Fraction matrix: Schur complements (Sylvester's law of inertia)"""
    M = [list(r) for r in Mx]
    m0 = len(M)
    pos = neg = 0
    while M:
        m = len(M)
        p = next((i for i in range(m) if M[i][i] != 0), None)
        if p is None:
            pr = next(((i, j) for i in range(m) for j in range(m) if i != j and M[i][j] != 0), None)
            if pr is None:
                break
            i, j = pr
            for k in range(m):                       # congruence e_i -> e_i + e_j
                M[i][k] += M[j][k]
            for k in range(m):
                M[k][i] += M[k][j]
            p = i
        piv = M[p][p]
        if piv > 0:
            pos += 1
        else:
            neg += 1
        idx = [k for k in range(m) if k != p]
        M = [[M[a][b] - M[a][p] * M[p][b] / piv for b in idx] for a in idx]
    return pos, neg, m0 - pos - neg


# ============================================================================================ X0 inputs
def x0():
    section("X0  Inputs: S, M, l, b, c1, c13 from the referee Clifford construction (xcheck_round8 X0 -> xcheck_round6 X0-X2)")
    t1 = time.time()
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        tabs = X8.x0()
    ok = len(X8.RES) == 1 and X8.RES[0][1] and ': False' not in buf.getvalue()
    check(f"xcheck_round8.x0() rerun as a library (38 internal checks of xcheck_round6: J1 on all 14 covectors, J2 on all 196 "
          f"tensors and sampled pairs; tables equal J's printed tables), both signatures {tick(t1)}", ok)
    return tabs


def cone_data(T, ep):
    """own assembly for p = (x, 0, 0, w): Sigma = -M^-1 S, E^(1) and E^(13) coefficients (J2 'Substituting the slaving')"""
    sub = {y: 0, z: 0}
    S = T['S'].subs(sub)
    M = T['M'].subs(sub)
    D = sp.expand(M.det())
    Sig = (-M.inv() * S).applyfunc(sp.cancel)
    dS = {0: Sig.diff(t), 1: Sig.diff(u)}
    out = {}
    for nm in ('E1', 'E13'):
        l_ = [sp.sympify(T['L'][(nm, k)]).subs(sub) for k in KS]
        H = [sp.cancel(sum(l_[i] * Sig[i, a] for i in range(4))) for a in range(2)]
        Q, Qp = {}, {}
        for a in range(2):
            for b in range(2):
                v = sum(l_[i] * dS[b][i, a] for i in range(4))
                vp = -v
                for (ik, k), (il, ll) in itertools.combinations_with_replacement(list(enumerate(KS)), 2):
                    b1, b2, b3 = (sp.sympify(e_).subs(sub) for e_ in T['B'][(nm, k, ll)])
                    wkl = sp.Rational(1, 2) if k == ll else 1
                    v += wkl * (b1 * Sig[ik, a] * Sig[il, b] + b2 * Sig[il, a] * Sig[ik, b])
                    vp += wkl * b3 * Sig[ik, a] * Sig[il, b]
                Q[(a, b)] = v
                Qp[(a, b)] = vp
        Qs = {(a, b): sp.cancel((Q[(a, b)] + Q[(b, a)]) / 2) for a in range(2) for b in range(2)}
        out[nm] = dict(H=H, Q=Qs, Qp={k: sp.cancel(v) for k, v in Qp.items()}, W=sp.cancel((Q[(0, 1)] - Q[(1, 0)]) / 2))
    return dict(D=D, Sig=Sig, E1=out['E1'], E13=out['E13'], c1=sp.expand(T['c1'].subs(sub)), c13=sp.expand(T['c13'].subs(sub)))


# ============================================================================================ X1 Q1(b)
def x1():
    section("X1  Lemma Q1(b): the swap metric -g_aleph is Z1 with (-cH, -cV), by brute force off the section")
    rng = random.Random(20261131)
    for eta, lam, (cH, cV) in (((1, -1, -1, -1), Fr(1, 2), (Fr(3), Fr(2))), ((1, 1, -1, -1), Fr(2, 3), (Fr(1), Fr(5)))):
        t1 = time.time()
        base = X10.rand_base(rng, eta)             # g = eta + random linear, quadratic, cubic terms: Gamma(0) != 0
        xP = [Fr(0)] * n
        gx = X10.gmat(eta)
        Q = [[Fr(int(i == j)) + Fr(rng.randint(-1, 1), rng.randint(3, 5)) for j in range(n)] for i in range(n)]
        yP = X10.mmul(X10.mT(Q), X10.mmul(gx, Q))   # y0 = Q^T g(0) Q != g(0), same signature
        res = {}
        for s in (1, -1):
            M, Psi0, _, _ = X10.zorro_jets(base, xP, yP, lam, s * cH, s * cV)
            g, dg, ddg = X10.two_jet(X10.S14, M)
            Ric = X10.ricci(g, dg, ddg)
            c = 4 * s * cV * (1 - 4 * lam)
            E = [Fr(0)] * n + X10.vec_of(yP)
            okE = all(sum(Ric[a][b] * E[b] for b in range(14)) == -sum(g[a][b] * E[b] for b in range(14)) / c for a in range(14))
            gi = X10.minv(g)
            sc = sum(gi[a][b] * Ric[b][a] for a in range(14) for b in range(14))
            G = [[Ric[a][b] - sc / 2 * g[a][b] for b in range(14)] for a in range(14)]
            Ms, _, _, _ = X10.zorro_jets(base, xP, gx, lam, s * cH, s * cV)
            gs = [[Ms[a][b].get(0, Fr(0)) for b in range(14)] for a in range(14)]
            res[s] = dict(g=g, dg=dg, ddg=ddg, Ric=Ric, G=G, okE=okE, sig=inertia(gs), sigoff=inertia(g),
                          Gam=any(v != 0 for row in Psi0 for v in row))
        r1, r2 = res[1], res[-1]
        neg = lambda T_: [neg(v) for v in T_] if isinstance(T_, list) else -T_  # noqa: E731
        oklin = r2['g'] == neg(r1['g']) and r2['dg'] == neg(r1['dg']) and r2['ddg'] == neg(r1['ddg'])
        okR = r2['Ric'] == r1['Ric'] and r2['G'] == r1['G']
        c1_, b1_, bS1 = X10.consts(lam, cV)
        c2_, b2_, bS2 = X10.consts(lam, -cV)
        okc = c2_ == -c1_ and b2_ == -b1_ and bS2 == -bS1 and 4 * b1_ * c1_ == 4 * b2_ * c2_
        sigs = (eta.count(1), eta.count(-1))
        check(f"base {sigs} in general coordinates (random polynomial of degree 3, Gamma(0) != 0: {r1['Gam']}), lambda = {lam}, "
              f"(cH, cV) = ({cH}, {cV}), at a point (0, y0), y0 = Q^T g(0) Q != g(0), off the section: the 2-jet of Z1(-cH, -cV) is minus that of Z1(cH, cV) ({oklin}); Ric and the "
              f"Einstein tensor G coincide ({okR}); Ric(E, .) = -E^flat/c for both, with c -> -c ({r1['okE'] and r2['okE']}); c, b, "
              f"beta_S change sign, B = 4bc = {4 * b1_ * c1_} does not ({okc}); inertia at the section point (0, g(0)) = "
              f"{r1['sig']}, {r2['sig']}, at the off-section point {r1['sigoff']}, {r2['sigoff']} {tick(t1)}",
              oklin and okR and okc and r1['okE'] and r2['okE'] and r1['sig'][:2] == (7, 7) and r2['sig'][:2] == (7, 7)
              and r1['sig'][2] == 0 and r1['Gam'])


# ============================================================================================ X2 Q2-Q5
def x2(tabs):
    section("X2  Lemma Q2 (G_lambda pairing), Proposition Q3 (F'), Proposition Q4 (c1 - s c13), Lemma Q5 (M3's channel)")
    lam = sp.Symbol('lambda')
    mu = lam / (1 - 4 * lam)
    rng = random.Random(20261132)
    ok = True
    for eta in ((1, 1, 1, 1), (1, -1, -1, -1), (1, 1, -1, -1), (-1, 1, 1, 1)):
        while True:
            P = sp.Matrix(4, 4, lambda i, j: int(i == j) + sp.Rational(rng.randint(-1, 1), rng.randint(2, 5)))
            if P.det() != 0:
                break
        Y = P.T * sp.diag(*eta) * P
        Yi = Y.inv()
        Zr = sp.Matrix(4, 4, lambda i, j: sp.Rational(rng.randint(-4, 4), rng.randint(1, 3)))
        Z = Zr + Zr.T
        X_ = sp.Matrix([rng.randint(-3, 3) for _ in range(4)])
        W_ = sp.Matrix([rng.randint(-3, 3) for _ in range(4)])
        Xf, Wf = Y * X_, Y * W_
        sym = (Xf * Wf.T + Wf * Xf.T) / 2 + mu * (X_.T * Y * W_)[0, 0] * Y
        G = (Yi * Z * Yi * sym).trace() - lam * (Yi * Z).trace() * (Yi * sym).trace()
        ok &= sp.cancel(G - (X_.T * Z * W_)[0, 0]) == 0
    check("Q2: for y = P^T eta P in four signatures, symmetric Z, vectors X, Y, symbolic lambda: G_lambda|_y(Z, X^flat . Y^flat "
          "+ mu y(X,Y) y) = Z(X,Y) with mu = lambda/(1 - 4 lambda) (so cV G_lambda(Z, sym (nabla_X Y)^V) = -(cH/2) Z(X,Y), N1(a))",
          ok)
    th, ep_ = sp.symbols('theta epsilon')
    Dm = kap ** 2 * (kap ** 2 - 147456 * ep_ * th ** 2)
    ratio = sp.cancel((ep_ * 184025088 * kap ** 5 * th / Dm ** 2) / (576 * kap ** 3 / Dm))
    Fp = (kap ** 2 - 147456 * ep_ * th ** 2) ** sp.Rational(-13, 12)
    okF = sp.simplify(sp.diff(Fp, th) / Fp - ratio) == 0
    check(f"Q3: with M8's Xi_0 coefficients, F''/F' must be {ratio}; F' = |kappa^2 - 147456 eps theta^2|^(-13/12) gives it, so "
          f"r = F(theta) has traceless Hessian 0 (concircular) where d theta != 0", okF)
    T7 = tabs[-1]
    oks = True
    for s in (1, -1):
        a = sp.Symbol('a')
        sub = {y: s * x, z: s * a * x, w: a * x}
        c1, c13 = sp.expand(T7['c1'].subs(sub)), sp.expand(T7['c13'].subs(sub))
        sg = x * (t - s * u)
        oks &= sp.expand(x * (c1 - s * c13) - sg * (kap + 208 * sg)) == 0
    check("Q4: (7,7) null members (x, s x, s a x, a x), s = +-1, from the referee c1, c13: x (c1 - s c13) = sigma (kappa + 208 "
          "sigma), sigma = t x + eps u y = x (t - s u) (the G-terms cancel since x_13 = s x_1)", oks)
    cds = {}
    for ep in (1, -1):
        t1 = time.time()
        cd = cone_data(tabs[ep], ep)
        cds[ep] = cd
        Sig = cd['Sig']
        ok = all(sp.simplify(Sig[k, 1]) == 0 for k in range(4))
        e13 = cd['E13']
        ok &= e13['H'] == [0, 0] and cd['E1']['H'][1] == 0
        ok &= all(e13['Q'][ab] == 0 and e13['Qp'][ab] == 0 for ab in ((0, 1), (1, 0), (1, 1)))
        ok &= sp.cancel(e13['Qp'][(0, 0)] + e13['Q'][(0, 0)]) == 0 and e13['W'] == 0
        ok &= sp.expand(cd['c13'] - u * (kap + 208 * x * t)) == 0
        ok &= sp.expand(cd['c1'] - kap * t - x * (312 * t ** 2 + ep * 104 * u ** 2)) == 0
        check(f"{LAB[ep]} Q5 inputs (M3) for (x, 0, 0, w), own assembly from the referee tables: Sigma has no du column; H_13 = 0, "
              f"H_1u = 0; E^(13) = u (kappa + 208 x t) g + Q^(13)_tt (dt dt - |dt|^2 g); c1 = kappa t + x (312 t^2 "
              f"{'+' if ep == 1 else '-'} 104 u^2) {tick(t1)}", ok)
    return cds


# ============================================================================================ X3 Q8
def x3(cds):
    section("X3  Proposition Q8: mu', mu mu'' of x (W_A + a W'_C Ad_Gamma) from the referee tables; own persistence certificate")
    m1s = {}
    for ep in (1, -1):
        t1 = time.time()
        cd = cds[ep]
        D0 = cd['D'].subs(u, 0)
        H0 = sp.cancel(cd['E1']['H'][0].subs(u, 0))
        q0 = sp.cancel(cd['E1']['Q'][(0, 0)].subs(u, 0))
        qp0 = sp.cancel(cd['E1']['Qp'][(0, 0)].subs(u, 0))
        D1, N2, L1, Q2, Nq = polys(ep)
        sub = {t: tau * kap / x, w: sp.sqrt(A_) * x}
        okf = sp.cancel(H0.subs(sub) - 576 * x ** 2 * D1 / (kap * Q2)) == 0
        okf &= sp.cancel(q0.subs(sub) + 36864 * x ** 3 * Nq / (kap ** 2 * Q2 ** 2)) == 0
        okf &= sp.expand(D0.subs(sub) - kap ** 4 * Q2) == 0
        Fp = sp.Function('Fp')(t)                      # F'(t), with F'' = (q/H) F'
        mu = x * Fp / H0
        rule = {sp.Derivative(Fp, t): q0 / H0 * Fp}
        mup = sp.cancel(sp.diff(mu, t).subs(rule) / Fp)                  # mu' = (d mu/dt)/F'
        mupp = sp.diff(mup, t) / Fp                                      # mu'' = (d mu'/dt)/F'
        mm = sp.cancel(mu * mupp)
        okm = not mm.has(Fp) and sp.cancel(mup - x * (q0 - sp.diff(H0, t)) / H0 ** 2) == 0
        m1 = sp.cancel(mup.subs(sub))
        mmt = sp.cancel(mm.subs(sub))
        okm &= sp.cancel(m1 + N2 / (18 * D1 ** 2)) == 0
        okm &= sp.cancel(mmt - A_ * (6 * A_ - 11 * ep) * L1 * Q2 / (81 * D1 ** 4)) == 0
        m1s[ep] = m1
        check(f"{LAB[ep]} on u = 0: H = 576 x^2 D1/(kappa Q2), q = -36864 x^3 Nq/(kappa^2 Q2^2), det M = kappa^4 Q2 ({okf}); with a "
              f"generic F (F'' = (q/H) F'), mu = x F'/H: mu' = x (q - H')/H^2 = -N2/(18 D1^2) and mu mu'' = A (6A - 11 eps) L1 Q2/"
              f"(81 D1^4), free of F ({okm}) {tick(t1)}", okf and okm)
    okc = True
    gg = {}
    for ep in (1, -1):
        D1, N2, L1, Q2, Nq = polys(ep)
        num = sp.Poly(sp.numer(sp.together(sp.diff(N2 / D1 ** 2, tau))), tau)
        g = sp.Poly(num.all_coeffs()[0], A_)
        for cf in num.all_coeffs()[1:]:
            g = sp.gcd(g, sp.Poly(cf, A_))
        gg[ep] = sp.factor(g.as_expr() / g.LC())
        okc &= sp.expand(gg[ep] - A_ * (A_ - sp.Rational(11, 6) * ep)) == 0
    check(f"Q8(a): the gcd over Q[A] of the tau-coefficients of the numerator of d mu'/d tau is {gg[1]} ((14,0)), {gg[-1]} ((7,7)): "
          f"mu' is constant on an interval iff A (6A - 11 eps) = 0", okc)
    store = {}
    oka = True
    for ep in (1, -1):
        t1 = time.time()
        D1, N2, L1, Q2, Nq = polys(ep)
        m1 = -N2 / (18 * D1 ** 2)
        Phi = 4 * A_ * (6 * A_ - 11 * ep) * L1 * Q2 / N2 ** 2
        X = 1 - 1 / (B_ * m1)
        Rtog = sp.together((Phi + 4 * X) ** 2 - X * Phi ** 2)       # the form used by the note (sympy together)
        Rt = sp.cancel(Rtog)                                       # reduced form
        num, den = sp.fraction(Rt)
        P = sp.Poly(sp.expand(num), tau)
        cfs = [sp.expand(c_) for c_ in P.all_coeffs()]
        gb = sp.groebner(cfs, B_, A_, order='lex')
        uni = [g_ for g_ in gb.exprs if not g_.has(B_)]
        assert len(uni) == 1
        pts = []
        for a_ in sorted(set(sp.Poly(uni[0], A_).all_roots())):
            gB = None
            for g_ in gb.exprs:
                if g_.has(B_):
                    pB = sp.Poly(sp.expand(g_.subs(A_, a_)), B_)
                    gB = pB if gB is None else sp.gcd(gB, pB)
            for b_ in sorted(set(gB.all_roots())):
                pts.append((a_, b_))
        want = sorted({(sp.Integer(0), sp.Rational(-18, 7)), (sp.Rational(11, 6) * ep, sp.Rational(-289, 3))})
        Xz = all(sp.cancel(X.subs({A_: a_, B_: b_})) == 0 for a_, b_ in pts)
        store[ep] = (pts, P.degree())
        check(f"{LAB[ep]} referee certificate: R~ = (Phi + 4X)^2 - X Phi^2, X = 1 - 1/(B mu'), reduced numerator of degree "
              f"{P.degree()} in tau (denominator {sp.factor(den)}); lex Groebner basis (B > A) of ALL its {len(cfs)} tau-"
              f"coefficients: {[sp.factor(g_) for g_ in gb.exprs]}; common zeros in C^2 (A, B) = {pts}, and X = 0 there "
              f"({Xz}) {tick(t1)}", sorted(pts) == want and Xz and P.degree() == 10)
        nA = sp.expand(sp.numer(Rtog))
        qq, rr = sp.div(sp.Poly(nA, tau, A_, B_), sp.Poly(sp.expand(num), tau, A_, B_))
        oka &= rr.is_zero and sp.Poly(nA, tau).degree() == 12
        qe = qq.as_expr()
        oka &= sp.cancel(qe / N2).is_number and sp.expand(N2.subs(tau, 0) - (7 - 3 * ep * A_)) == 0
    check("observation (R13-4): the un-cancelled numerator of sympy's together(), which the note uses, has degree 12 = (constant) x "
          "N2 x (reduced numerator); N2(0) = 7 - 3 eps A, so the factor (3A - 7 eps)^2 of the note's gcd and its candidate A = 7 "
          "eps/3 are artefacts of the common factor N2, not of R~", oka)
    lam_ = sp.Symbol('l')
    real_pos = [(a_, b_) for ep in (1, -1) for a_, b_ in store[ep][0] if a_ > 0]
    lv = [sp.solve(16 * (1 - 5 * lam_) - b_, lam_)[0] for _, b_ in real_pos]
    l0 = sp.solve(16 * (1 - 5 * lam_) + sp.Rational(18, 7), lam_)[0]
    check(f"hence for real a != 0 (A > 0): only (A, B) = {real_pos} ((14,0)), i.e. lambda = {lv}, > 1/4, not on the (14,0) side; "
          f"(A, B) = (0, -18/7) is W_A at lambda = {l0} (N7(d)). With mu' = 1/B constant only at A (6A - 11 eps) = 0, neither Q7(b) "
          f"nor Q7(c) holds on an interval for real a != 0 and admissible lambda != 1/5",
          real_pos == [(sp.Rational(11, 6), sp.Rational(-289, 3))] and lv == [sp.Rational(337, 240)] and l0 == sp.Rational(13, 56))


# ============================================================================================ X4 Q9
def x4():
    section("X4  Lemma Q9: (a) the blocks of -x Ric + alpha g + q dt dt = 0 at p in s(X); (b) common zeros of H and q")
    xx, cV, lam, cH, al, q, Et = sp.symbols('x c_V lambda c_H alpha q Et')
    c = 4 * cV * (1 - 4 * lam)
    b = (1 - 5 * lam) / (cV * (1 - 4 * lam))
    bS = sp.Rational(-5, 4) / cV
    al0 = sp.solve(-xx * bS + al, al)[0]                       # V_0 block, traceless part forces q theta0 theta0 = 0
    EE = sp.solve(xx + al0 * c + q * Et ** 2, q * Et ** 2)     # (E, E): -x Ric(E,E) + alpha |E|^2 + q (Et)^2, Ric(E,E) = -1
    HH = sp.simplify(-(-xx * (-cH / c) + al0 * cH))            # HH: -x (Ric_g - cH/c g) + alpha cH g + q dt_H dt_H
    ok = sp.simplify(1 + c * bS + b * c) == 0 and sp.simplify(EE[0] - xx * b * c) == 0 and sp.simplify(HH - xx * cH * b) == 0
    check("Q9(a): with N6(a) (Ric(E,E) = -1, Ric|V0 = beta_S g, Ric_HH = Ric_g - (cH/c) g, Ric_HV = 0): the V0 block gives alpha = "
          "x beta_S; then (E,E) gives q (Et)^2 = x b c, and HH gives x Ric_g + x cH b g = q dt_H dt_H (1 + c beta_S = -b c)", ok)
    res = {}
    okr = True
    for ep in (1, -1):
        D1, N2, L1, Q2, Nq = polys(ep)
        rq = sp.factor(sp.resultant(D1, Nq, tau))
        rQ = sp.factor(sp.resultant(D1, Q2, tau))
        res[ep] = (rq, rQ)
        okr &= sp.cancel(rq / (A_ * (A_ - ep) * (6 * A_ - 11 * ep) ** 2)).is_number
        okr &= sp.cancel(rQ / (A_ * (6 * A_ - 11 * ep) ** 2)).is_number
    D1, N2, L1, Q2, Nq = polys(1)
    okr &= D1.subs(A_, 1) == 2
    tH = sp.solve(D1.subs(A_, sp.Rational(11, 6)), tau)[0]
    okr &= Q2.subs({A_: sp.Rational(11, 6), tau: tH}) == 0
    check(f"Q9(b): Res_tau(D1, Nq), Res_tau(D1, Q2) = {res}; for A > 0 they vanish only in (14,0) at A = 1 (D1 = 2: H has no zero) "
          f"and A = 11/6 (the zero tau_H = {tH} of H is a zero of det M): H and q have no common zero where det M != 0", okr)


# ============================================================================================ X5 algebra
def x5():
    section("X5  Lemmas Q6, Q7; Proposition Q11 step 4; Proposition Q12 (sympy, from the component equations)")
    cH, cV, lam, k, mu0, sg, vE, nE = sp.symbols('c_H c_V lambda k mu_0 sigma v_E n_E')
    c = 4 * cV * (1 - 4 * lam)
    b = (1 - 5 * lam) / (cV * (1 - 4 * lam))
    bS = sp.Rational(-5, 4) / cV
    scal = sg * cH
    # step 3 (pair (X,U), output H, paired with Z):  -(cH/4) v_E U(X,Z) = mu0 [ -(cH b/2) U(X,Z) + (1/2) Ric_g(X, U~ Z) ]
    vEs = sp.solve(-cH / 4 * vE - mu0 * (-cH * b / 2 + scal / 8), vE)[0]
    # step 4 (pair (E,U), output V0): -k v_E c beta_S - n_E c = 0
    nEs = sp.solve(-k * vE * c * bS - nE * c, nE)[0]
    # step 5 (pair (X,E), output H): -(1/4) v_E + k c v_E (sigma/4 - 1/c) + c n_E = mu0 sigma/8
    st5 = -vE / 4 + k * c * vE * (sg / 4 - 1 / c) + c * nE - mu0 * sg / 8
    st5 = sp.simplify(st5.subs(nE, nEs).subs(vE, vEs) / mu0)
    target = k * c * sg ** 2 - 4 * b * (4 * k * b * c - 1)
    ok6 = sp.simplify(vEs - mu0 * (2 * b - sg / 2)) == 0 and sp.simplify(nEs + k * bS * vE) == 0
    ok6 &= sp.simplify(st5 * (-8) - target) == 0
    check("Q6: from the component equations with gradient coefficient k = mu'(r0) and right side mu0 = mu(r0) (each re-derived by "
          "hand from N1(e), N6(b)-(d), (f)): v_E = mu0 (2b - sigma_g/2), n_E = -k beta_S v_E, and step 5 is -(1/8) [k c sigma_g^2 "
          "- 4b (4kbc - 1)] = 0, i.e. Q6(c)", ok6)
    rr, rh, Bs = sp.symbols('r rho B')
    ps = sp.symbols('p0:4')
    mu = sum(ps[i] * rr ** i for i in range(4))                # a generic cubic mu(r)
    mu1, mu2 = sp.diff(mu, rr), sp.diff(mu, rr, 2)
    f = Bs * mu1 * (1 - rh ** 2) - 1                            # (C')
    Lf = sp.diff(f, rr) * Bs / 2 * mu * (1 - rh) + sp.diff(f, rh) * (-rh)   # d/dtheta along (B) and rho' = -rho
    target = Bs * mu1 / 2 * (mu * mu2 / mu1 ** 2 * (1 - rh) + 4 * rh ** 2)
    ok7 = sp.cancel(Lf - target - Bs * mu * mu2 * (1 - rh) / (2 * mu1) * f) == 0
    bb, cc, m1_ = sp.symbols('b c m1')
    ok7 &= sp.expand(cc * m1_ * (16 * bb ** 2 * rh ** 2 - 16 * bb ** 2) + 4 * bb
                     + 4 * bb * (4 * bb * cc * m1_ * (1 - rh ** 2) - 1)) == 0
    ok7 &= sp.expand(cc * mu * (2 * bb - 4 * bb * rh / 2) - 4 * bb * cc / 2 * mu * (1 - rh)) == 0
    check("Q7: with e^-theta sigma_g = 4 b rho, (C) c mu'(e^-2theta sigma^2 - 16 b^2) = -4b is -4b [B mu'(1 - rho^2) - 1] (B = 4bc) and "
          "(B) is dr/dtheta = (B/2) mu (1 - rho); for a generic cubic mu, the derivative of f = B mu'(1 - rho^2) - 1 along theta "
          "(r' = (B/2) mu (1 - rho), rho' = -rho) equals (B mu'/2)[Phi (1 - rho) + 4 rho^2] + (B mu mu''(1 - rho)/(2 mu')) f, "
          "Phi = mu mu''/mu'^2: so f = 0 on an interval gives (D)", ok7)
    Hh, Ee, cH_, xx = sp.symbols('P_H P_E c_H x')     # Hh = tr_g d d T~|_HH, Ee = d d T~(E,E)
    s2, beta, cs = sp.symbols('s2 beta c')
    eqH = -xx * (Hh / cH_ - 2 * s2) + 4 * beta          # trace over H w.r.t. g_aleph|_H = cH g of the HH block (dim 4)
    eqE = -xx * (Ee - s2 / 2 * cs) + beta * cs          # (E, E) component, g_aleph(E, E) = c
    bsol = sp.solve(eqE, beta)[0]
    o2 = sp.simplify(eqH.subs(beta, bsol) / (-4 * xx))
    ok11 = sp.simplify(o2 - (Hh / (4 * cH_) - Ee / cs)) == 0
    check(f"Q11 step 4: eliminating beta_AB and d_A d_B tr T~ between the H-trace of the HH block, -x[P_HH - s'' g/2] + beta g, "
          f"and the (E,E) component gives exactly (O2): tr_g P_HH/(4 cH) - P(E,E)/c = {o2} = 0", ok11)
    mu_, mup = sp.symbols('mu mu1')
    f1 = 2 * bb * cc * mu_
    f2 = 2 * bb * cc * mup * f1                          # f'' = 2bc mu'(f) f'
    cond = sp.factor(f2 - bb * cc * mu_)
    ok12 = sp.expand(cond - bb * cc * mu_ * (4 * bb * cc * mup - 1)) == 0
    check(f"Q12: on the flat model, r = f(s): the H and E blocks give f' = 2bc mu(f), f'' = bc mu(f); compatible iff {cond} = 0",
          ok12)


# ============================================================================================ X6 Q11 alpha''
def o2_nested(base, eta, cH, cV, lam=Fr(1, 5)):
    """second fibre derivatives of T~ = Ric + g/c at p = (0, eta) along V_0 by a nested multi-dual shift; returns (O2) contracted
    with G_0^-1 on V_0 (G_0 = tr(eta^-1 . eta^-1 .)), also with g_aleph^-1 on V_0, and the order-0/1 data"""
    g0 = X10.gmat(eta)
    V0 = X10.V0_basis(g0)
    k9 = range(9)
    yP = [[None] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            A = X10.AOF[(i, j)]
            inner = MD(g0[i][j], {jj: V0[jj][A] for jj in k9 if V0[jj][A]})
            yP[i][j] = MD(inner, {ii: V0[ii][A] for ii in k9 if V0[ii][A]})
    M, Psi0, _, _ = X10.zorro_jets(base, [Fr(0)] * n, yP, lam, cH, cV)
    g, dg, ddg = X10.two_jet(X10.S14, M)
    Ric = X10.ricci(g, dg, ddg)
    c = 4 * cV * (1 - 4 * lam)
    Tt = [[Ric[a][b] + g[a][b] * (1 / c) for b in range(14)] for a in range(14)]

    def d2(xv, i, j):
        if isinstance(xv, MD) and i in xv.d and isinstance(xv.d[i], MD):
            return xv.d[i].d.get(j, Fr(0))
        return Fr(0)

    def d1(xv, i):
        if isinstance(xv, MD) and i in xv.d:
            v = xv.d[i]
            return v.a if isinstance(v, MD) else v
        return Fr(0)
    o0 = all(X10.base_val(Tt[a][b]) == 0 for a in range(14) for b in range(14))
    o1 = all(d1(Tt[a][b], i) == 0 for a in range(14) for b in range(14) if (a < n) == (b < n) for i in k9)
    hv1 = any(d1(Tt[a][b], i) != 0 for a in range(n) for b in range(n, 14) for i in k9)
    g0i = X10.minv(g0)
    GS = [[X10.mtr(X10.mmul(X10.mmul(g0i, X10.mat_of(V0[i])), X10.mmul(g0i, X10.mat_of(V0[j])))) for j in k9] for i in k9]
    GSi = X10.minv(GS)
    E = [Fr(0)] * n + X10.vec_of(g0)
    S2 = Fr(0)
    for i in k9:
        for j in k9:
            if GSi[i][j]:
                trH = sum(Fr(eta[m_]) * d2(Tt[m_][m_], i, j) for m_ in range(n))
                EE = sum(E[a] * E[b] * d2(Tt[a][b], i, j) for a in range(14) if E[a] for b in range(14) if E[b])
                S2 += GSi[i][j] * (trH / (4 * cH) - EE / c)
    return dict(S2=S2, S2al=S2 / cV, o0=o0, o1=o1, hv1=hv1, flatPsi=all(v == 0 for row in Psi0 for v in row))


def o2_diag(base, eta, cH, cV, lam=Fr(1, 5)):
    """second route: sum_k d^2_{U_k} [tr_g T~_HH/(4cH) - T~(E,E)/c] / G_0(U_k, U_k) over a G_0-orthogonal basis of V_0 (eta = 1)"""
    basis = []
    for (i, j) in X10.IDX:
        if i != j:
            Mx = [[Fr(0)] * n for _ in range(n)]
            Mx[i][j] = Mx[j][i] = Fr(1)
            basis.append(Mx)
    for dg_ in ((1, -1, 0, 0), (1, 1, -2, 0), (1, 1, 1, -3)):
        basis.append([[Fr(dg_[i]) if i == j else Fr(0) for j in range(n)] for i in range(n)])
    g0 = X10.gmat(eta)
    c = 4 * cV * (1 - 4 * lam)
    E = [Fr(0)] * n + X10.vec_of(g0)
    tot = Fr(0)
    for U in basis:
        nrm = X10.mtr(X10.mmul(U, U))
        yP = [[MD(MD(g0[i][j], {0: U[i][j]} if U[i][j] else {}), {0: U[i][j]} if U[i][j] else {}) for j in range(n)]
              for i in range(n)]
        M, _, _, _ = X10.zorro_jets(base, [Fr(0)] * n, yP, lam, cH, cV)
        g, dg, ddg = X10.two_jet(X10.S14, M)
        Ric = X10.ricci(g, dg, ddg)
        Tt = [[Ric[a][b] + g[a][b] * (1 / c) for b in range(14)] for a in range(14)]

        def dd(xv):
            return xv.d[0].d.get(0, Fr(0)) if isinstance(xv, MD) and 0 in xv.d and isinstance(xv.d[0], MD) else Fr(0)
        trH = sum(Fr(eta[m_]) * dd(Tt[m_][m_]) for m_ in range(n))
        EE = sum(E[a] * E[b] * dd(Tt[a][b]) for a in range(14) if E[a] for b in range(14) if E[b])
        tot += (trH / (4 * cH) - EE / c) / nrm
    return tot


def x6():
    section("X6  Proposition Q11: (O2) = alpha'' |W|^2 at lambda = 1/5, base (4,0) (REVIEW8 jets and Ricci routine, own contraction)")
    eta = (1, 1, 1, 1)
    W1 = X10.weyl_from(X10.rand_curv(random.Random(20261133), eta), eta)
    W2 = X10.weyl_from(X10.rand_curv(random.Random(20261134), eta), eta)

    def nrm2(W):
        return sum(Fr(eta[k[0]] * eta[k[1]] * eta[k[2]] * eta[k[3]]) * W[k] ** 2 for k in W)

    def ricflat(W):
        return all(sum(Fr(eta[a]) * W[(a, b, a, d)] for a in range(n)) == 0 for b in range(n) for d in range(n))
    zero = {k: Fr(0) for k in W1}
    t1 = time.time()
    fl = o2_nested(X10.weyl_base(zero, eta), eta, Fr(2), Fr(3))
    check(f"flat base, (cH, cV) = (2, 3): T~(p) = 0, its HH, VV blocks have zero first fibre derivatives, and (O2) = {fl['S2']} "
          f"(the constant term of the quadratic function of W is 0, N2) {tick(t1)}", fl['o0'] and fl['o1'] and fl['S2'] == 0)
    vals, a2s = {}, {}
    for tag, W, cub, (cH, cV) in (('W1', W1, None, (Fr(2), Fr(3))), ('W1+cubic', W1, 20261135, (Fr(2), Fr(3))),
                                  ('2 W1', {k: 2 * v for k, v in W1.items()}, None, (Fr(2), Fr(3))),
                                  ('W2+cubic', W2, 20261136, (Fr(2), Fr(3))),
                                  ('W2+cubic', W2, 20261136, (Fr(1), Fr(1))),
                                  ('W1+cubic', W1, 20261135, (Fr(1, 2), Fr(5)))):
        t1 = time.time()
        base = X10.weyl_base(W, eta, random.Random(cub) if cub else None)
        Bg = X10.base_geometry(base)
        okB = all(Bg['Rl'][k] == W[k] for k in W) and ricflat(W)
        r = o2_nested(base, eta, cH, cV)
        a2 = r['S2'] / nrm2(W)
        pred = Fr(-3) * cV / (2 * cH ** 2)
        vals[(tag, cH, cV)] = r['S2']
        a2s[(tag, cH, cV)] = a2
        check(f"{tag}, (cH, cV) = ({cH}, {cV}): R_g(0) = W Weyl ({okB}), cubic terms {'yes' if cub else 'no'}; T~(p) = 0 and the "
              f"HH, VV blocks of d T~ vanish at p ({r['o0']}, {r['o1']}); HV block of d T~ at p nonzero: {r['hv1']}; |W|^2 = "
              f"{nrm2(W)}; (O2) contracted with G_0^-1 = {r['S2']} = alpha'' |W|^2, alpha'' = {a2} (note: -3 cV/(2 cH^2) = {pred}); "
              f"contracted with g_aleph^-1 on V_0: alpha''/cV = {a2 / cV} {tick(t1)}",
              okB and r['o0'] and r['o1'] and r['flatPsi'] and a2 == pred and r['hv1'] == bool(cub))
    check(f"alpha'' is unchanged by the cubic terms (nabla R(0) != 0; W1 at (2,3): {vals[('W1+cubic', 2, 3)]} = {vals[('W1', 2, 3)]}), "
          f"and quadratic: 2 W1 gives 4 times W1 ({vals[('2 W1', 2, 3)] == 4 * vals[('W1', 2, 3)]})",
          vals[('W1+cubic', 2, 3)] == vals[('W1', 2, 3)] and vals[('2 W1', 2, 3)] == 4 * vals[('W1', 2, 3)])
    t1 = time.time()
    dv = o2_diag(X10.weyl_base(W1, eta, random.Random(20261135)), eta, Fr(2), Fr(3))
    check(f"second route (nine single-direction second derivatives over a G_0-orthogonal basis of V_0): (O2) = {dv}, equal to the "
          f"nested-shift value {vals[('W1+cubic', 2, 3)]} {tick(t1)}", dv == vals[('W1+cubic', 2, 3)])
    ratios = sorted({cV / cH for (_, cH, cV) in a2s})
    okm = all(v == Fr(-3) * cV / (2 * cH ** 2) for (_, cH, cV), v in a2s.items()) and len(ratios) >= 3
    check(f"referee (R13-3): alpha'' = -3 cV/(2 cH^2) at three ratios cV/cH = {[str(q_) for q_ in ratios]} (G_0-contraction); "
          f"homothety (degree -1) and the fibre scaling N1(b) (degree -2 in cH) make alpha'' a monomial C cV/cH^2, so these fix C and "
          f"alpha'' != 0 for all cH, cV != 0; contracted with g_aleph^-1 on V_0 it is -3/(2 cH^2)", okm)


# ============================================================================================ X7 Q13
def x7():
    section("X7  Proposition Q13: the x = +-y line of M7(b) and the allowed lambda")
    xx, ll, kk, w_, uu = sp.symbols('x ell kappa w u', positive=True)
    zeta = kk * w_ + 208 * ll * uu
    Hc = 312 * xx * ll / zeta                       # Hessian coefficient on sigma = -kappa/208 (M7(b), xcheck_round9 Z4)
    qc = -105456 * xx * ll ** 2 / zeta ** 2         # dt-free gradient coefficient (du (x) du)
    Fp = zeta ** sp.Rational(-13, 8)
    okF = sp.simplify(sp.diff(Fp, uu) / Fp - qc / Hc) == 0
    F = sp.integrate(Fp, uu)
    mu = xx * Fp / Hc
    m1 = sp.simplify(mu / F)
    lam = sp.Symbol('lambda')
    lv = {str(mv): sp.solve(16 * mv * (1 - 5 * lam) - 1, lam)[0]
          for mv in (sp.Rational(-5, 12), sp.Rational(-7, 18), sp.Rational(-3, 289))}
    ok = okF and m1 == sp.Rational(-5, 12) and lv['-5/12'] == sp.Rational(23, 100) and sp.Rational(23, 100) < sp.Rational(1, 4)
    ok &= lv['-7/18'] == sp.Rational(13, 56) and lv['-3/289'] == sp.Rational(337, 240)
    check(f"from M7(b)'s line coefficients (Hessian 312 x ell/zeta, gradient -105456 x ell^2/zeta^2): F' = zeta^(-13/8) solves "
          f"F''/F' = q/H ({okF}), F = {sp.simplify(F)}, mu = x F'/H = {m1} F, linear; Q7(d) (sigma_g = 0, B mu_1 = 1, B = 16(1 - 5 "
          f"lambda)): lambda in {{1/5, {lv['-5/12']}}}, both < 1/4 (not (7,7)); the same rule gives {lv} for mu_1 = -5/12, -7/18, -3/289",
          ok)


# ============================================================================================ X8 Q15
def x8(tabs):
    section("X8  Proposition Q15: flat-base constant solutions for every (14,0) member, without M1(a); det M at the roots")
    T = tabs[1]
    c1, c13 = sp.expand(T['c1']), sp.expand(T['c13'])
    sg = t * x + u * y
    dl = t * y - u * x
    okc = sp.expand(y * c1 - x * c13 - dl * (kap + 208 * sg)) == 0
    okc &= sp.expand(x * c1 + y * c13 - (kap * sg + 312 * sg ** 2 + 104 * dl ** 2)) == 0
    xp, yp, zp, wp, s, r, cV = sp.symbols("x' y' z' w' s r c_V", real=True)
    k = sp.Symbol('kappa', positive=True)
    re = {x: -sp.I * xp, y: -sp.I * yp, z: -sp.I * zp, w: -sp.I * wp, t: sp.I * s, u: sp.I * r, kap: k}
    sgr = sp.expand(sg.subs(re))
    dlr = sp.expand(dl.subs(re))
    Nr = sp.expand((x ** 2 + y ** 2).subs(re))
    okc &= sp.expand(sgr - (xp * s + yp * r)) == 0 and sp.expand(dlr - (yp * s - xp * r)) == 0
    okc &= sp.expand(Nr + xp ** 2 + yp ** 2) == 0
    check("(14,0) referee c1, c13 for a general member (x, y, z, w): y c1 - x c13 = delta (kappa + 208 sigma) and x c1 + y c13 = "
          "kappa sigma + 312 sigma^2 + 104 delta^2 (sigma = t x + u y, delta = t y - u x); in the real form sigma = x' s + y' r, delta "
          "= y' s - x' r are real and N = x^2 + y^2 = -(x'^2 + y'^2)", okc)
    S_, D_, K2 = sp.symbols('sigma delta K2', real=True)
    Kp = sp.Symbol('K', positive=True)
    # I6 with G = 15/(2 cV) g: c_i = 15 x_i/(2 cV).  Combinations: delta (kappa + 208 sigma) = 0 and
    # kappa sigma + 312 sigma^2 + 104 delta^2 = 15 N/(2 cV) = -15 (x'^2 + y'^2)/(2 cV)
    e2 = k * S_ + 312 * S_ ** 2 + 104 * D_ ** 2 + sp.Rational(15, 2) * Kp / cV       # Kp = x'^2 + y'^2
    disc = sp.factor(sp.discriminant(e2.subs(D_, 0), S_))
    br = sp.expand(e2.subs(S_, -k / 208))
    okE = sp.simplify(disc - (k ** 2 - 9360 * Kp / cV)) == 0 and sp.expand(br - (k ** 2 / 416 + 104 * D_ ** 2 + sp.Rational(15, 2) * Kp / cV)) == 0
    check(f"Q15 without the rotation: c_i = 15 x_i/(2 cV) iff delta (kappa + 208 sigma) = 0 and kappa sigma + 312 sigma^2 + 104 delta^2 "
          f"+ 15 (x'^2 + y'^2)/(2 cV) = 0; on delta = 0 the discriminant is {disc} (>= 0 iff cV kappa^2 >= 9360 (x'^2 + y'^2), cV > 0); "
          f"on sigma = -kappa/208 the left side is {br} > 0; so the solutions are (t, u) || (x, y) (u = 0 only after M1(a))", okE)
    M = T['M']
    th_ = t * w - u * z
    sub_line = {t: S_ * x / (x ** 2 + y ** 2), u: S_ * y / (x ** 2 + y ** 2)}
    Dl = sp.factor(sp.cancel(M.det().subs(sub_line)))
    Afull = sp.cancel(((x * w - y * z) / (x ** 2 + y ** 2)) ** 2)
    pred = kap ** 2 * ((kap + 224 * S_) * (kap + 544 * S_) - 147456 * Afull * S_ ** 2)
    okd = sp.cancel(Dl - pred) == 0
    okd &= sp.cancel(Afull.subs(re) - ((xp * wp - yp * zp) / (xp ** 2 + yp ** 2)) ** 2) == 0
    check("det M on the line (t, u) = sigma (x, y)/N, general (14,0) member: kappa^2 [(kappa + 224 sigma)(kappa + 544 sigma) - "
          "147456 A sigma^2] with A = ((x w - y z)/N)^2 = ((x' w' - y' z')/(x'^2 + y'^2))^2 (= a^2 for the rotated member y = 0)", okd)
    om, Aa = sp.symbols('omega A')
    sg1 = sp.Symbol('sigma1')                                    # sigma = kappa sigma1 (homogeneity)
    p1 = 312 * sg1 ** 2 + sg1 + om / 1248                        # 15 (x'^2 + y'^2)/(2 cV) = omega kappa^2/1248
    dM = (1 + 224 * sg1) * (1 + 544 * sg1) - 147456 * Aa * sg1 ** 2
    P = sp.Poly(sp.expand(sp.resultant(p1, dM, sg1)), om)
    lc, c0 = sp.factor(P.LC()), sp.factor(P.nth(0))
    at1 = sp.factor(P.as_expr().subs(om, 1))
    prop = sp.solve(sp.Eq(312 * 768, 121856 - 147456 * Aa), Aa)
    at0 = [v for v in sp.solve(P.as_expr().subs(Aa, 0), om) if v > 0]
    rl, r0, r1_ = sp.cancel(lc / (144 * Aa - 119) ** 2), sp.cancel(c0 / (2304 * Aa + 319)), sp.cancel(at1 / (576 * Aa - 125) ** 2)
    okr = P.degree() == 2 and rl.is_number and rl > 0 and r0.is_number and r0 < 0 and r1_.is_number and r1_ > 0
    okr &= prop == [sp.Rational(-115, 144)] and at0 == [sp.Rational(1131, 1156)]
    okr &= sp.Rational(9360) / at0[0] == sp.Rational(277440, 29)
    check(f"Res_sigma(312 sigma^2 + kappa sigma + omega kappa^2/1248, det M/kappa^2) (kappa = 1 by homogeneity), omega = 9360 "
          f"(x'^2 + y'^2)/(cV kappa^2): degree 2 in omega, leading coefficient {lc} >= 0, constant term {c0} < 0, value at omega = 1 "
          f"{at1} >= 0. Hence for every real a EXACTLY ONE omega_* in (0, 1] (one value cV kappa^2/(x'^2 + y'^2) = 9360/omega_* "
          f">= 9360) puts a root on det M = 0; both roots at once would need A = {prop}; a = 0: omega_* = {at0} (cV kappa^2 = "
          f"277440/29, N12(c))", okr)
    Ae = sp.Rational(125, 576)
    roots = sp.solve(p1.subs(om, 1), sg1)
    dz = sp.expand(dM.subs({Aa: Ae, sg1: sp.Rational(-1, 624)}))
    check(f"referee (R13-1): at the threshold omega = 1 (cV kappa^2 = 9360 (x'^2 + y'^2)) the quadratic has the double root sigma = "
          f"kappa {roots}, and for A = a^2 = 125/576 (a real, e.g. a = 5 sqrt5/24) det M vanishes there ({dz}): so Q15's 'det M != 0 "
          f"at one of the two roots for every a' fails for that a at the threshold, where no root is admissible under (H1)",
          roots == [sp.Rational(-1, 624)] and dz == 0)


def main():
    tabs = x0()
    x1()
    cds = x2(tabs)
    x3(cds)
    x4()
    x5()
    x6()
    x7()
    x8(tabs)
    nt = sum(1 for _, ok in RES if ok)
    say(f"\nSUMMARY: {len(RES)} checks, {nt} True, {len(RES) - nt} False")
    say(f"runtime {time.time() - T0:.0f} s")


if __name__ == '__main__':
    main()
