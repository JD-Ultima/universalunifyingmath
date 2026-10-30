"""GU round 9 (M): the members of the family that round 8 (L) leaves open, under the hypotheses (H) of L.0.

Formalization: formalization/geometric-unity/M-z2-remaining-members.md.  Builds on L (L1-L9), J (J1, J2, Cor. J2.1) and
J5d (L0, L1).  Imports L_nonradial.py (and through it J_slaved_sector.py, I_elliptic_reduction.py, F_jets.py, ...),
unmodified.  Exact arithmetic throughout (Fractions, sympy rationals / Gaussian rationals); fixed seeds; one core.
Every check is asserted.

Notation (as in L): member p = (x, y, z, w), eps = +1 (14,0), -1 (7,7); ell = xz + eps yw, N = x^2 + eps y^2,
delta = ty - ux, sigma = tx + eps uy, theta = tw - uz.  All reduced coefficients are handled as polynomial numerators over
det M (H D, Q D^2, Q' D^2, varpi D^2), as in L's section L4.

Sections (labels = results of the note)
  M0  inputs: L's own section L0 (J1's S, M, c1, c13 recomputed with J's j0/j1; the J2 tables spot-checked against J's
      evaluator), run silently and asserted.
  M1  Lemma M1 (symmetries of the reduced system): the rotation (14,0) / boost (7,7) of (x,y), (z,w), (t,u) (symbolic
      one-parameter group); the (7,7) swap p -> -Kp, (t,u) -> (u,t), g_aleph -> -g_aleph.
  M2  Proposition M2 (the ell = 0 cone): orbit representatives; the Lambda^2 coefficients varpi on the cone.
  M3  Proposition M3 (the first-order channel of (x,0,0,w)): Sigma has no du column, H13 = 0, E13 = u (kappa + 208 x t) g
      + Q13_tt (dt dt - |dt|^2 g); the mirror statement for (0,y,z,0) in (7,7).
  M4  Proposition M4 (the single equation on u = 0): H, q, q', mu'(tau); mu' = -7/18 at w = 0; zeros of H and mu'.
  M5  Lemma M5 (fibre obstruction for nonlinear mu): pairings (L7's constants), Groebner certificate that the order-1
      condition (C1)_0 at mu' = 0 has no Ricci-flat solution, brute-force Taylor check with a random Weyl tensor.
  M6  Theorem M7 (N = 0 in (7,7)): (a) ell = 0: E1 - sgn E13 = (c1 - sgn c13) g, T_e depends on d sigma only;
      (b) ell != 0: R' = 0, K = 0, the invariant line sigma = -kappa/208, the single equation there, mu = -(5/12) r.
  M7  Proposition M8 ((x, y) = 0): everything depends on theta; the derivative-free combination; the single equation and
      its trace identity.
  M8  Proposition M9 (off the section): the first integrability condition of the fibre equation, relaxed in mu', is
      inconsistent at sample points y0 = P g P^T (random and constant curvature, four base signatures); consistent on s(X).
Run (from computation/geometric-unity): python3 M_remaining.py > M_remaining.out 2>&1   (~6 min on one core)
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
from sympy.polys.matrices import DomainMatrix  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import L_nonradial as L  # noqa: E402

T0 = time.time()
NCHK = [0]
KS = L.KS
x, y, z, w = L.x, L.y, L.z, L.w
t, u, kap = L.t, L.u, L.kap
I = sp.I
SG = {1: '(14,0)', -1: '(7,7)'}
EX = sp.expand
cc, m_ = sp.symbols('c m')


def P(*a):
    print(*a)
    sys.stdout.flush()


def check(label, ok):
    NCHK[0] += 1
    P(f"  [{NCHK[0]:3d}] {label}: {ok}")
    assert ok, label


def section(s):
    P(f"\n=== {s}   [{time.time() - T0:.0f} s]")


def subs_sim(tab, sub):
    """simultaneous substitution in the tables (L.subs_tables substitutes sequentially)"""
    S, M, Lt, B = tab
    f = lambda e: sp.sympify(e).subs(sub, simultaneous=True)  # noqa: E731
    return (S.applyfunc(f), M.applyfunc(f), {k: f(v) for k, v in Lt.items()},
            {k: (tuple(f(a) for a in v) if isinstance(v, tuple) else f(v)) for k, v in B.items()})


def redpoly(tab):
    """polynomial numerators of the reduced coefficients: D = det M, SN = -adj(M) S (Sigma = SN/D), HN = l_E SN (H = HN/D),
    QN, QpN (Q = QN/D^2, Q' = QpN/D^2) and WN (varpi = WN/D^2); same formulas as L.terms / J.terms"""
    S, M, Lt, B = tab
    D = EX(M.det())
    SN = (-M.adjugate() * S).applyfunc(EX)
    dD = {'t': EX(sp.diff(D, t)), 'u': EX(sp.diff(D, u))}
    dS2 = {f: (SN.diff({'t': t, 'u': u}[f]) * D - SN * dD[f]).applyfunc(EX) for f in 'tu'}
    lE = sp.Matrix([[Lt[(c, k)] for k in KS] for c in ('E1', 'E13')])
    HN = (lE * SN).applyfunc(EX)
    QN, QpN, WN = {}, {}, {}
    for nm in ('E1', 'E13'):
        for a, fa in enumerate('tu'):
            for b_, fb in enumerate('tu'):
                v = sum(Lt[(nm, k)] * dS2[fb][ik, a] for ik, k in enumerate(KS))
                vp = -v
                for (ik, k), (il, l) in itertools.combinations_with_replacement(list(enumerate(KS)), 2):
                    b1, b2, b3 = B[(nm, k, l)]
                    wkl = sp.Rational(1, 2) if k == l else 1
                    v += wkl * (b1 * SN[ik, a] * SN[il, b_] + b2 * SN[il, a] * SN[ik, b_])
                    vp += wkl * b3 * SN[ik, a] * SN[il, b_]
                QN[(nm, fa, fb)], QpN[(nm, fa, fb)] = EX(v), EX(vp)
        WN[nm] = EX((QN[(nm, 't', 'u')] - QN[(nm, 'u', 't')]) / 2)
    for nm in ('w3', 'w11'):
        v = sum(Lt[(nm, k)] * (dS2['t'][ik, 1] - dS2['u'][ik, 0]) for ik, k in enumerate(KS))
        for (ik, k), (il, l) in itertools.combinations(list(enumerate(KS)), 2):
            v += B[(nm, k, l)] * (SN[ik, 0] * SN[il, 1] - SN[ik, 1] * SN[il, 0])
        WN[nm] = EX(v)
    return dict(D=D, SN=SN, HN=HN, QN=QN, QpN=QpN, WN=WN)


def cfun(eps, sub=None):
    """c1, c13 of I6 (as in L)"""
    c1 = kap * t + x * (312 * t ** 2 + eps * 104 * u ** 2) + eps * 208 * y * t * u
    c13 = kap * u + 208 * x * t * u + y * (104 * t ** 2 + eps * 312 * u ** 2)
    if sub:
        c1, c13 = c1.subs(sub, simultaneous=True), c13.subs(sub, simultaneous=True)
    return EX(c1), EX(c13)


def iszero(e):
    return EX(e) == 0


# ====================================================================== M0
def m0():
    section("M0. Inputs: L's section L0 (J1 S, M, c1, c13 via J's j0/j1; J2 tables against J's evaluator), run silently")
    buf = io.StringIO()
    n0 = L.NCHK[0]
    with contextlib.redirect_stdout(buf):
        L.l0()
    nl = L.NCHK[0] - n0
    okl = nl == 8 and buf.getvalue().count(': True') == 8 and ': False' not in buf.getvalue()
    check(f"L.l0() passes all its {nl} asserted checks (both signatures): the transcribed S, M, c1, c13, l and b tables used"
          f" below are those of J", okl)


# ====================================================================== M1
def m1():
    section("M1. Lemma M1: symmetries of the reduced system (J1-J2 tables)")
    for eps in (1, -1):
        S, M, Lt, B = L.tables(eps)
        den = 1 + eps * m_ ** 2
        cs, sn = (1 - eps * m_ ** 2) / den, 2 * m_ / den           # cs^2 + eps sn^2 = 1
        R = sp.Matrix([[cs, -eps * sn], [sn, cs]])
        okR = sp.simplify(R.T * sp.diag(1, eps) * R - sp.diag(1, eps)) == sp.zeros(2, 2) and sp.simplify(R.det() - 1) == 0

        def rot(a, b):
            v = R * sp.Matrix([a, b])
            return v[0], v[1]
        (xn, yn), (zn, wn), (tn, un) = rot(x, y), rot(z, w), rot(t, u)
        sub = {x: xn, y: yn, z: zn, w: wn, t: tn, u: un}
        f = lambda e: sp.together(sp.sympify(e).subs(sub, simultaneous=True))  # noqa: E731
        zero = lambda e: sp.cancel(e) == 0  # noqa: E731
        okM = all(zero(f(M[i, j]) - M[i, j]) for i in range(4) for j in range(4))
        Sn = S.applyfunc(f)
        okS = all(zero(v) for v in (Sn * R - S))
        okT = True
        for k in KS:
            for a_, b_ in (('E1', 'E13'), ('w3', 'w11')):
                vv = sp.Matrix([Lt[(a_, k)], Lt[(b_, k)]])
                okT &= all(zero(e) for e in (vv.applyfunc(f) - R * vv))
        for k, l in itertools.combinations_with_replacement(KS, 2):
            for j in range(3):
                vv = sp.Matrix([B[('E1', k, l)][j], B[('E13', k, l)][j]])
                okT &= all(zero(e) for e in (vv.applyfunc(f) - R * vv))
            vv = sp.Matrix([B[('w3', k, l)], B[('w11', k, l)]])
            okT &= all(zero(e) for e in (vv.applyfunc(f) - R * vv))
        c1, c13 = cfun(eps)
        vv = sp.Matrix([c1, c13])
        okc = all(zero(e) for e in (vv.applyfunc(f) - R * vv))
        vv = sp.Matrix([x, y])
        okx = all(zero(e) for e in (vv.applyfunc(f) - R * vv))
        grp = 'SO(2)' if eps == 1 else 'SO+(1,1)'
        check(f"{SG[eps]} R = R(m) in {grp} (symbolic m; R^T diag(1,eps) R = diag(1,eps)) acting on (x,y), (z,w), (t,u):"
              f" M(Rp, Rtau) = M(p, tau), S(Rp) R = S(p) (so T_e is invariant), and the columns (l_E1, l_E13), (l_w3, l_w11),"
              f" (b_E1, b_E13), (b_w3, b_w11), (c1, c13), (x_1, x_13) all transform by R", okR and okM and okS and okT and okc and okx)
    # (7,7) swap: p -> -K p, (t,u) -> (u,t), g -> -g
    t1 = time.time()
    eps = -1
    tab = L.tables(eps)
    r0 = redpoly(tab)
    r1 = redpoly(subs_sim(tab, {x: -y, y: -x, z: -w, w: -z}))
    f = lambda e: EX(e.subs({t: u, u: t}, simultaneous=True))  # noqa: E731
    sw = {'t': 'u', 'u': 't'}
    ok = iszero(f(r1['D']) - r0['D'])
    for a_, b_, ia, ib in (('E1', 'E13', 0, 1), ('E13', 'E1', 1, 0)):
        for j in (0, 1):
            ok &= iszero(f(r1['HN'][ia, j]) + r0['HN'][ib, 1 - j])
        for f1 in 'tu':
            for g1 in 'tu':
                ok &= iszero(f(r1['QN'][(a_, f1, g1)]) + r0['QN'][(b_, sw[f1], sw[g1])])
                ok &= iszero(f(r1['QpN'][(a_, f1, g1)]) + r0['QpN'][(b_, sw[f1], sw[g1])])
    okW = all(iszero(f(r1['WN'][a_]) - r0['WN'][b_]) for a_, b_ in (('E1', 'E13'), ('E13', 'E1'), ('w3', 'w11'), ('w11', 'w3')))
    c1, c13 = cfun(eps)
    sub = {x: -y, y: -x, z: -w, w: -z, t: u, u: t}
    okc = iszero(c1.subs(sub, simultaneous=True) - c13) and iszero(c13.subs(sub, simultaneous=True) - c1)
    check(f"(7,7) swap p' = -K p = (-y, -x, -w, -z), (t', u') = (u, t): det M invariant; H'_(1,phi') = -H_(13,K phi),"
          f" H'_(13,phi') = -H_(1,K phi), and the same for Q, Q' (indices swapped); c'_1 = c_13, c'_13 = c_1, x'_1 = -x_13;"
          f" varpi'_(1,13,3,11) = varpi_(13,1,11,3).  Hence E'^(1)[-g] = -E^(13)[g], E'^(13)[-g] = -E^(1)[g]"
          f"  [{time.time() - t1:.0f} s]", ok and okW and okc)


# ====================================================================== M2
def m2():
    section("M2. Proposition M2: the ell = 0 cone; orbit representatives; varpi on the cone")
    # orbit maps (exact, rational group elements)
    okm = True
    # (14,0): (x,y) = (3,4), (z,w) = c (-y, x); rotation with (cos, sin) = (3/5, -4/5) sends (3,4) to (5,0)
    R = sp.Matrix([[sp.Rational(3, 5), sp.Rational(4, 5)], [-sp.Rational(4, 5), sp.Rational(3, 5)]])
    okm &= R * sp.Matrix([3, 4]) == sp.Matrix([5, 0]) and R * sp.Matrix([-4 * cc, 3 * cc]) == sp.Matrix([0, 5 * cc])
    # (7,7) N > 0: (5,3), (z,w) = c (y, x); boost (cosh, sinh) = (5/4, -3/4) sends (5,3) to (4,0)
    Bm = sp.Matrix([[sp.Rational(5, 4), -sp.Rational(3, 4)], [-sp.Rational(3, 4), sp.Rational(5, 4)]])
    okm &= Bm * sp.Matrix([5, 3]) == sp.Matrix([4, 0]) and Bm * sp.Matrix([3 * cc, 5 * cc]) == sp.Matrix([0, 4 * cc])
    # (7,7) N < 0: (3,5) -> (0,4); (z,w) = c (5,3) -> (4c, 0)
    okm &= Bm * sp.Matrix([3, 5]) == sp.Matrix([0, 4]) and Bm * sp.Matrix([5 * cc, 3 * cc]) == sp.Matrix([4 * cc, 0])
    # the cone: ell = 0 and (x,y) != 0  <=>  (z,w) = c J (x,y), J = (-y, x) in (14,0), K = (y, x) in (7,7)
    okc = iszero((x * (-cc * y) + y * (cc * x))) and iszero(x * (cc * y) - y * (cc * x))
    # on the cone, xw - yz = c N
    okc &= iszero(x * (cc * x) - y * (-cc * y) - cc * (x ** 2 + y ** 2)) and iszero(x * (cc * x) - y * (cc * y) - cc * (x ** 2 - y ** 2))
    # N = 0 in (7,7) and ell = 0: (z,w) = c K (x, sgn x) = c sgn (x, sgn x), i.e. (z,w) parallel to (x,y)
    for sg in (1, -1):
        okc &= sp.Matrix([sg * x, x]) == sg * sp.Matrix([x, sg * x])
    check("orbit representatives: (14,0) (3,4,-4c,3c) ~ (5,0,0,5c); (7,7) (5,3,3c,5c) ~ (4,0,0,4c) and (3,5,5c,3c) ~ (0,4,4c,0)"
          " (exact rational rotation/boost); on the cone (z,w) = c J(x,y) [(14,0)], c K(x,y) [(7,7)], xw - yz = c N; for"
          " N = 0 in (7,7) the cone members have (z,w) parallel to (x,y)", okm and okc)
    for eps in (1, -1):
        r = redpoly(subs_sim(L.tables(eps), {y: 0, z: 0}))
        W = r['WN']
        q13 = r['QN'][('E13', 't', 't')]
        okW = iszero(W['E13']) and iszero(W['E1'] + eps * q13) and sp.cancel(q13 / u).subs(u, 0) != 0
        okW &= sp.Poly(sp.cancel(q13 / u), t, u, x, w, kap).as_expr().subs({u: 0, t: 0}) != 0
        check(f"{SG[eps]} (x,0,0,w): varpi^(1) = {'-' if eps == 1 else '+'}Q13_tt = (nonzero polynomial) u x^2 kappa/(det M)^2,"
              f" a nonzero rational function, so L1 applies on the whole cone component (x,y) != 0, N != 0 (also via M1)", okW)
    # N = 0, ell = 0 (7,7) and (x,y) = 0: all varpi vanish identically
    okz = True
    for sg in (1, -1):
        r = redpoly(subs_sim(L.tables(-1), {y: sg * x, z: sg * cc * x, w: cc * x}))
        okz &= all(iszero(v) for v in r['WN'].values())
    for eps in (1, -1):
        r = redpoly(subs_sim(L.tables(eps), {x: 0, y: 0}))
        okz &= all(iszero(v) for v in r['WN'].values())
    check("(7,7) N = 0, ell = 0 [(x, +-x, +-cx, cx)] and (x,y) = 0 [both signatures]: all four varpi vanish identically, so L1"
          " does not apply there (it is not needed: M7(a), M8)", okz)


# ====================================================================== M3
def m3():
    section("M3. Proposition M3: the first-order channel of (x, 0, 0, w) (and of (0, y, z, 0) in (7,7))")
    out = {}
    for eps in (1, -1):
        t1 = time.time()
        r = redpoly(subs_sim(L.tables(eps), {y: 0, z: 0}))
        out[eps] = r
        ok = all(iszero(r['SN'][i, 1]) for i in range(4)) and iszero(r['HN'][1, 0]) and iszero(r['HN'][1, 1]) \
            and iszero(r['HN'][0, 1])
        ok &= iszero(r['QpN'][('E13', 't', 't')] + r['QN'][('E13', 't', 't')])
        ok &= all(iszero(r['QN'][('E13', a, b)]) and iszero(r['QpN'][('E13', a, b)]) for a, b in (('t', 'u'), ('u', 't'), ('u', 'u')))
        ok &= iszero(r['WN']['E13'])
        _, c13 = cfun(eps, {y: 0, z: 0})
        ok &= iszero(c13 - u * (kap + 208 * x * t))
        check(f"{SG[eps]} (x,0,0,w), x, w, kappa symbolic: Sigma has no du column; H_13t = H_13u = H_1u = 0; varpi^(13) = 0;"
              f" Q13_tu = Q13_ut = Q13_uu = 0 = Q'13 (same); Q'13_tt = -Q13_tt; c13 = u (kappa + 208 x t).  So"
              f" E^(13) = u (kappa + 208 x t) g + Q13_tt (dt (x) dt - |dt|^2 g)  [{time.time() - t1:.0f} s]", ok)
    # mirror: (0, y, z, 0) in (7,7)
    r = redpoly(subs_sim(L.tables(-1), {x: 0, w: 0}))
    ok = all(iszero(r['SN'][i, 0]) for i in range(4)) and iszero(r['HN'][0, 0]) and iszero(r['HN'][0, 1]) and iszero(r['HN'][1, 0])
    ok &= iszero(r['QpN'][('E1', 'u', 'u')] + r['QN'][('E1', 'u', 'u')])
    ok &= all(iszero(r['QN'][('E1', a, b)]) and iszero(r['QpN'][('E1', a, b)]) for a, b in (('t', 'u'), ('u', 't'), ('t', 't')))
    ok &= iszero(r['WN']['E1'])
    c1, _ = cfun(-1, {x: 0, w: 0})
    ok &= iszero(c1 - t * (kap - 208 * y * u))
    check("(7,7) (0,y,z,0): Sigma has no dt column, H_1 = 0, H_13t = 0, E^(1) = t (kappa - 208 y u) g + Q1_uu (du (x) du -"
          " |du|^2 g) (the swap image of the statement above)", ok)
    return out


# ====================================================================== M4
def m4(out):
    section("M4. Proposition M4: the single equation on u = 0 for (x, 0, 0, w); mu'(tau), tau = x t/kappa, c = w/x")
    tau = sp.Symbol('tau')
    for eps in (1, -1):
        r = out[eps]
        D0 = r['D'].subs(u, 0)
        H0 = sp.cancel(r['HN'][0, 0].subs(u, 0) / D0)
        q0 = sp.cancel(r['QN'][('E1', 't', 't')].subs(u, 0) / D0 ** 2)
        qp0 = sp.cancel(r['QpN'][('E1', 't', 't')].subs(u, 0) / D0 ** 2)
        c1, _ = cfun(eps, {y: 0, z: 0, u: 0})
        Dref = kap ** 2 * (kap ** 2 + 768 * kap * t * x + 121856 * t ** 2 * x ** 2 - eps * 147456 * t ** 2 * w ** 2)
        Href = 576 * (kap * (x ** 2 + eps * w ** 2) + 544 * t * x * (x ** 2 - eps * w ** 2)) / (Dref / kap ** 2)
        ok = iszero(D0 - Dref) and sp.cancel(H0 - Href) == 0 and iszero(c1 - (kap * t + 312 * x * t ** 2))
        mup = sp.factor(sp.cancel(x * (q0 - sp.diff(H0, t)) / H0 ** 2))
        mu_tau = sp.factor(sp.cancel(mup.subs({t: tau * kap / x, w: cc * x})))
        num = 7 - eps * 3 * cc ** 2 + (7616 - eps * 3008 * cc ** 2 - 768 * cc ** 4) * tau \
            + (2071552 - eps * 3226624 * cc ** 2 + 1155072 * cc ** 4) * tau ** 2
        den = 18 * ((1 + eps * cc ** 2) + 544 * (1 - eps * cc ** 2) * tau) ** 2
        ok &= sp.cancel(mu_tau + num / den) == 0
        ok &= sp.cancel(mu_tau.subs(cc, 0) + sp.Rational(7, 18)) == 0
        A_, B_, C_ = sp.Poly(num, tau).all_coeffs()
        disc = sp.factor(B_ ** 2 - 4 * A_ * C_)
        ok &= iszero(disc - 16384 * cc ** 2 * (cc ** 2 + eps * 35) * (6 * cc ** 2 - eps * 11) ** 2)
        ok &= sp.Poly(num, tau, cc).as_expr() != 0
        check(f"{SG[eps]} on u = 0: det M = {sp.factor(Dref)}, H = 576 (kappa (x^2 {'+' if eps == 1 else '-'} w^2) + 544 t x (x^2"
              f" {'-' if eps == 1 else '+'} w^2))/(..), c1 = kappa t + 312 x t^2; the reparametrised coefficient"
              f" mu' = x (q - dH/dt)/H^2 = -[{num}] / [18 ((1 {'+' if eps == 1 else '-'} c^2) + 544 (1 {'-' if eps == 1 else '+'} c^2) tau)^2]"
              f" (a function of tau = x t/kappa and c = w/x only); mu' = -7/18 at"
              f" c = 0 (L8 step 5); discriminant of the numerator in tau = {disc}", ok)
    okz = True
    for eps in (1, -1):
        hn = (1 + eps * cc ** 2) + 544 * (1 - eps * cc ** 2) * tau
        okz &= sp.cancel(hn.subs(tau, (1 + eps * cc ** 2) / (544 * (eps * cc ** 2 - 1)))) == 0
    check("hence mu' is not identically 0 for any c; in (7,7) it has no real zero for c^2 < 35, in (14,0) it has real zeros for"
          " every c != 0 (both are handled by Lemma M5, not by avoiding them); H vanishes exactly at"
          " tau_H = (1 + eps c^2)/(544 (eps c^2 - 1)) (no zero iff eps c^2 = 1, i.e. c^2 = 1 in (14,0)); that value is handled"
          " by the V-block argument of Theorem M6", okz)


# ====================================================================== M5
def weyl_param(eta):
    """the generic Ricci-flat algebraic curvature tensor (10 parameters) in base eta"""
    n = 4
    pairs = [(a, b) for a in range(n) for b in range(a + 1, n)]
    syms = {(i, j): sp.Symbol(f'R_{i}{j}') for i in range(6) for j in range(i, 6)}
    bian = {syms[(2, 3)]: -syms[(0, 5)] + syms[(1, 4)]}

    def Rs0(a, b, c, d):
        if a == b or c == d:
            return sp.Integer(0)
        sg = 1
        if a > b:
            a, b, sg = b, a, -sg
        if c > d:
            c, d, sg = d, c, -sg
        i, j = pairs.index((a, b)), pairs.index((c, d))
        if j < i:
            i, j = j, i
        return sg * syms[(i, j)].subs(bian)
    free = [s for k, s in syms.items() if k != (2, 3)]
    ric = [EX(sum(eta[a] * Rs0(a, b, a, d) for a in range(n))) for b in range(n) for d in range(b, n)]
    sol = sp.solve(ric, free, dict=True)[0]
    wfree = [s for s in free if s not in sol]
    cache = {}

    def Rs(a, b, c, d):
        key = (a, b, c, d)
        if key not in cache:
            cache[key] = EX(sp.sympify(Rs0(a, b, c, d)).subs(sol))
        return cache[key]
    okr = all(iszero(sum(eta[a] * Rs(a, b, a, d) for a in range(n))) for b in range(n) for d in range(n))
    return Rs, wfree, okr


def m5():
    section("M5. Lemma M5: the fibre obstruction for nonlinear mu (order 1 at mu' = 0), and the Ricci-flat certificate")
    # pairings (L7's constants) -> at mu' = 0, <(C1)_0, c1> = 18 + 9 (c_V/c_H)^2 |R|^2
    for eta in ((1, 1, 1, 1), (1, -1, -1, -1), (-1, 1, 1, 1), (1, 1, -1, -1)):
        t1 = time.time()
        Pr = L.Pairing(eta)
        p01, p11, p00 = Pr.pair(Pr.c0, Pr.c1), Pr.pair(Pr.c1, Pr.c1), Pr.pair(Pr.c0, Pr.c0)
        Kb = L.okbasis(eta)
        okl = all(Pr.quad(K1, K2, Pr.c1) == -9 * L.mtr2(K1, K2) for K1, K2 in itertools.product(Kb, repeat=2))
        check(f"base {(eta.count(1), eta.count(-1))}: <c0, c1> = {p01}, <c1, c1> = {p11}, <c0, c0> = {p00}, <q(K (x) K'), c1> ="
              f" -9 tr(K K') on all 36 basis elements; with sum eps_i eps_j tr(R_ij^2) = -|R|^2 (L, output L6) this gives"
              f" <(C1)_0, c1> = 18 + 9 (c_V/c_H)^2 |R|^2, so (C1)_0 forces |R|^2 = -2 c_H^2/c_V^2 (impossible for a definite"
              f" base)  [{time.time() - t1:.0f} s]", Pr.ok and p01 == 18 and p11 == 144 and p00 == Fr(27, 2) and okl)
    # Groebner certificate: (C1)_0 together with Ric = 0 has no (complex) solution
    for eta in ((1, -1, -1, -1), (-1, 1, 1, 1), (1, 1, -1, -1)):
        t1 = time.time()
        Rs, wfree, okr = weyl_param(eta)
        T = L.c1_tensor(eta, L.curv_B(eta, Rs), 0)
        Pb = L.pbasis(eta)
        eqs = []
        for X in Pb:
            for iy, Yv in enumerate(Pb):
                for iz, Zv in enumerate(Pb):
                    if iz <= iy:
                        continue
                    for e in T(X, Yv, Zv):
                        e = EX(e)
                        if e != 0:
                            eqs.append(e)
        polys = [sp.Poly(e, *wfree) for e in eqs]
        mons = sorted({mm for p in polys for mm in p.monoms()})
        rows = [[sp.QQ.from_sympy(sp.sympify(dict(zip(p.monoms(), p.coeffs())).get(mm, 0))) for mm in mons] for p in polys]
        Dm = DomainMatrix(rows, (len(rows), len(mons)), sp.QQ)
        rk = Dm.rank()
        Mx = sp.Matrix([[sp.Rational(v) for v in r] for r in rows])
        _, piv = Mx.T.rref()
        sel = [eqs[i] for i in piv]
        G = sp.groebner(sel, *wfree, order='grevlex')
        check(f"base {(eta.count(1), eta.count(-1))}: (C1)_0 for the generic Ricci-flat curvature tensor (10 Weyl parameters,"
              f" Ric = 0 checked: {okr}; {len(eqs)} scalar equations, rank {rk} = {len(sel)} independent ones): Groebner basis"
              f" over Q = {list(G.exprs)}, i.e. no solution even over C (scale (c_V/c_H)^2 absorbed into R)"
              f"  [{time.time() - t1:.0f} s]", okr and rk == len(sel) and list(G.exprs) == [1])
    # brute-force control: Taylor system of Hess_F w = Bt (mu constant, mu' = 0) at y = eta, random Weyl tensor with |W|^2 < 0
    eta = (1, -1, -1, -1)
    Rs, wfree, _ = weyl_param(eta)
    t1 = time.time()
    for seed in itertools.count(80):
        rng = random.Random(seed)
        vals = {s: sp.Rational(rng.randint(-3, 3), rng.randint(1, 2)) for s in wfree}
        Wl = lambda a, b, c, d: L.fr(Rs(a, b, c, d).subs(vals))  # noqa: E731
        W2 = sum(eta[a] * eta[b] * eta[c] * eta[d] * Wl(a, b, c, d) ** 2 for a, b, c, d in itertools.product(range(4), repeat=4))
        if W2 < 0:
            break
    af = Fr(-2) / W2
    fd = L.fibre_data(eta, Wl, maxorder=2, afac=af)
    unk, rowkeys, cols, _, rhs = L.hessian_system(fd, 2, 4, mu0=Fr(1), mu1=Fr(0))
    r1, r2 = L.dm_from_cols(len(rowkeys), cols).rank(), L.dm_from_cols(len(rowkeys), cols, rhs).rank()
    check(f"brute force, base (1, 3): random Ricci-flat W (seed {seed}, |W|^2 = {W2} < 0), A-term factor (c_V/c_H)^2 = -2/|W|^2"
          f" = {af} > 0 (so the c1-pairing condition holds), mu' = 0: the Taylor system for Hess_F w = Bt to order eps^2"
          f" ({len(cols)} unknowns) is {'consistent' if r1 == r2 else 'inconsistent'} (ranks {r1}, {r2}), as the certificate"
          f" predicts  [{time.time() - t1:.0f} s]", r1 != r2)
    t1 = time.time()
    R1 = lambda a, b, c, d: Fr((eta[a] if a == c else 0) * (eta[b] if b == d else 0) - (eta[a] if a == d else 0) * (eta[b] if b == c else 0))  # noqa: E731
    fd = L.fibre_data(eta, R1, maxorder=2, afac=Fr(-1, 12))
    unk, rowkeys, cols, _, rhs = L.hessian_system(fd, 2, 4, mu0=Fr(1), mu1=Fr(0))
    r1, r2 = L.dm_from_cols(len(rowkeys), cols).rank(), L.dm_from_cols(len(rowkeys), cols, rhs).rank()
    check(f"control (as in L): R = R_1 (not Ricci-flat), A-term factor -1/12 (so |R|^2 (c_V/c_H)^2 = -2), mu' = 0: consistent"
          f" (ranks {r1}, {r2}); the solvable case exists, so the Ricci-flat inconsistency is not an artefact"
          f"  [{time.time() - t1:.0f} s]", r1 == r2)


# ====================================================================== M6
def m6():
    section("M6. Theorem M7: N = 0 in (7,7), y = sgn x (a) ell = 0; (b) ell != 0: the invariant line sigma = -kappa/208")
    eps = -1
    for sg in (1, -1):
        t1 = time.time()
        r = redpoly(subs_sim(L.tables(eps), {y: sg * x, z: sg * cc * x, w: cc * x}))
        SN = r['SN']
        ok = all(iszero(SN[i, 1] + sg * SN[i, 0]) for i in range(4))
        ok &= all(iszero(r['HN'][0, j] - sg * r['HN'][1, j]) for j in range(2))
        ok &= all(iszero(r['QN'][('E1', f1, g1)] - sg * r['QN'][('E13', f1, g1)]) and
                  iszero(r['QpN'][('E1', f1, g1)] - sg * r['QpN'][('E13', f1, g1)]) for f1 in 'tu' for g1 in 'tu')
        c1, c13 = cfun(eps, {y: sg * x})
        sig = t * x - u * sg * x
        ok &= iszero(x * (c1 - sg * c13) - sig * (kap + 208 * sig))
        check(f"(a) (x, {'+' if sg == 1 else '-'}x, {'+' if sg == 1 else '-'}cx, cx) [ell = 0 = N]: Sigma (dt, du) = Sigma_t (dt"
              f" {'-' if sg == 1 else '+'} du) = Sigma_t d sigma / x; the derivative coefficients of E^(1) and {'' if sg == 1 else '-'}E^(13)"
              f" coincide; so E^(1) - sgn E^(13) = (c1 - sgn c13) g with x (c1 - sgn c13) = sigma (kappa + 208 sigma),"
              f" sigma = tx - uy  [{time.time() - t1:.0f} s]", ok)
    for sg in (1, -1):
        t1 = time.time()
        r = redpoly(subs_sim(L.tables(eps), {y: sg * x}))
        D, HN = r['D'], r['HN']
        hN = (sg * HN[0, :] - HN[1, :]).applyfunc(EX)            # (y, -x) H / x
        ell = x * z - sg * x * w
        vt, vu = sp.cancel(sg * x * hN[1] / (576 * ell)), sp.cancel(-sg * x * hN[0] / (576 * ell))
        okp = sp.denom(vt) == 1 and sp.denom(vu) == 1
        V = sp.Matrix([vt, vu])
        Qpp = sp.Matrix(2, 2, lambda i, j: sg * r['QN'][('E1', 'tu'[i], 'tu'[j])] - r['QN'][('E13', 'tu'[i], 'tu'[j])])
        Qppp = sp.Matrix(2, 2, lambda i, j: sg * r['QpN'][('E1', 'tu'[i], 'tu'[j])] - r['QpN'][('E13', 'tu'[i], 'tu'[j])])
        Dv = sp.Matrix([[sp.diff(vt, t), sp.diff(vt, u)], [sp.diff(vu, t), sp.diff(vu, u)]])
        acc = (Dv * V).applyfunc(EX)
        Rp = EX((hN * acc)[0] * D + (V.T * Qpp * V)[0])
        K = EX((V.T * (Qpp + Qppp) * V)[0])
        okhv = iszero(hN[0] * vt + hN[1] * vu)
        c1, c13 = cfun(eps, {y: sg * x})
        de, sig = sg * x * t - u * x, x * t - sg * x * u
        okc = iszero(sg * x * c1 - x * c13 - de * (kap + 208 * sig))
        inv = sp.factor(EX(vt * sp.diff(de, t) + vu * sp.diff(de, u)))
        okinv = sp.cancel(inv / (kap * ell * (kap + 16 * sig) * (kap + 208 * sig))) == -1
        check(f"(b) (x, {'+' if sg == 1 else '-'}x, z, w), ell = x (z {'-' if sg == 1 else '+'} w) != 0: kernel field v of the"
              f" curvature-free row h = (y,-x) H polynomial ({okp}), h.v = 0; along v the traceless residual R' = 0 AND the"
              f" trace coefficient K = v^T (Q'' + Q''') v = 0 identically, so E' = yE^(1) - xE^(13) = (y c1 - x c13) g ="
              f" delta (kappa + 208 sigma) g; v.grad(delta) = -kappa ell (kappa + 16 sigma)"
              f"(kappa + 208 sigma): the line sigma = -kappa/208 is v-invariant, delta = sgn sigma = 0 is not"
              f"  [{time.time() - t1:.0f} s]", okp and okhv and Rp == 0 and K == 0 and okc and okinv)
        # the single equation on the line sigma = -kappa/208, parametrised by u (t = sgn u - kappa/(208 x), dt = sgn du)
        line = {t: sg * u - kap / (208 * x)}
        Dl = sp.factor(D.subs(line))
        zeta = kap * w + 208 * u * ell
        okD = sp.cancel(Dl - 576 * kap ** 2 * zeta ** 2 / (169 * x ** 2)) == 0
        s_ = (sg, 1)
        co = {}
        for i, ch in ((0, 'E1'), (1, 'E13')):
            hh = sp.cancel((sg * HN[i, 0] + HN[i, 1]).subs(line) / D.subs(line))
            qq = sp.cancel(sum(r['QN'][(ch, a, b)] * s_[ia] * s_[ib] for ia, a in enumerate('tu') for ib, b in enumerate('tu')).subs(line)
                           / D.subs(line) ** 2)
            qp = sp.cancel(sum(r['QpN'][(ch, a, b)] * s_[ia] * s_[ib] for ia, a in enumerate('tu') for ib, b in enumerate('tu')).subs(line)
                           / D.subs(line) ** 2)
            co[ch] = (hh, qq, qp)
        hh, qq, qp = co['E1']
        okl = sp.cancel(hh - 312 * x * ell / zeta) == 0 and sp.cancel(qq + 105456 * x * ell ** 2 / zeta ** 2) == 0 \
            and sp.cancel(qp - 52728 * x * ell ** 2 / zeta ** 2) == 0
        okl &= all(sp.cancel(a - sg * b) == 0 for a, b in zip(co['E1'], co['E13']))
        okl &= sp.cancel(c1.subs(line) - sg * c13.subs(line)) == 0
        mup = sp.cancel(x * (qq - sp.diff(hh, u)) / hh ** 2)
        zu = sp.diff(zeta, u)
        # F' = |zeta|^(-13/8): F''/F' = q/h; mu = x F'/h = -(5/12) F with F = |zeta|^(-5/8) / (-(5/8) zeta_u)
        okF = sp.cancel(qq / hh + sp.Rational(13, 8) * zu / zeta) == 0 and sp.cancel(x * zeta / hh - sp.Rational(2, 3) * zeta ** 2 / zu) == 0
        okmu = sp.cancel(-sp.Rational(5, 8) * x * zu / (hh * zeta) + sp.Rational(5, 12)) == 0
        check(f"   on the line: det M = 576 kappa^2 zeta^2/(169 x^2), zeta = kappa w + 208 u ell; E^(13) = sgn E^(1) there, with"
              f" Hessian coefficient 312 x ell/zeta, gradient coefficients -105456 x ell^2/zeta^2, 52728 x ell^2/zeta^2;"
              f" mu' = {mup}; with F' = |zeta|^(-13/8), r = F: (Hess r)_0 = mu (Ric)_0 with mu = -(5/12) r exactly", okD and okl
              and mup == sp.Rational(-5, 12) and okF and okmu)


# ====================================================================== M7
def m7():
    section("M7. Proposition M8: (x, y) = 0; everything depends on theta = tw - uz")
    th = sp.Symbol('theta')
    G, Hs, dd, gg, B2 = sp.symbols('Gs Hess dthdth g dth2')   # formal tensors: Hess theta, dtheta dtheta, g, |dtheta|^2
    for eps in (1, -1):
        r = redpoly(subs_sim(L.tables(eps), {x: 0, y: 0}))
        D = r['D']
        okD = iszero(D - kap ** 2 * (kap ** 2 - eps * 147456 * (t * w - u * z) ** 2))
        SN = r['SN']
        ok = all(iszero(SN[i, 0] * z + SN[i, 1] * w) for i in range(4))        # columns proportional to (w, -z)
        rowc = {'E1': eps * w, 'E13': -z}
        for i, ch in ((0, 'E1'), (1, 'E13')):
            ok &= iszero(r['HN'][i, 0] - 576 * kap ** 3 * rowc[ch] * w) and iszero(r['HN'][i, 1] + 576 * kap ** 3 * rowc[ch] * z)
            for a, sa in (('t', w), ('u', -z)):
                for b, sb in (('t', w), ('u', -z)):
                    ok &= iszero(r['QN'][(ch, a, b)] - eps * 184025088 * kap ** 5 * rowc[ch] * (t * w - u * z) * sa * sb)
                    ok &= iszero(r['QpN'][(ch, a, b)] + eps * 92012544 * kap ** 5 * rowc[ch] * (t * w - u * z) * sa * sb)
        c1, c13 = cfun(eps, {x: 0, y: 0})
        ok &= iszero(z * c1 + eps * w * c13 - kap * (z * t + eps * w * u))
        okf = iszero(z * rowc['E1'] + eps * w * rowc['E13'])
        # single equation: (t,u) = theta (eps w, -z)/(z^2 + eps w^2) on z t + eps w u = 0
        ok1 = iszero((z * (eps * w) + eps * w * (-z))) and iszero((eps * w) * w - (-z) * z - (z ** 2 + eps * w ** 2))
        # trace identity: substitute Hess theta = -phi dth dth, box theta = -phi |dth|^2 with phi = F''/F'
        Dth = kap ** 2 * (kap ** 2 - eps * 147456 * th ** 2)
        chi, psi = 576 * kap ** 3 / Dth, eps * 184025088 * kap ** 5 * th / Dth ** 2
        phi = psi / chi
        okphi = sp.cancel(phi - eps * 319488 * th / (kap ** 2 - eps * 147456 * th ** 2)) == 0 and \
            sp.cancel(phi + sp.Rational(13, 12) * sp.diff(sp.log(kap ** 2 - eps * 147456 * th ** 2), th)) == 0
        Xi_g = kap * th / (z ** 2 + eps * w ** 2) + chi * phi * B2 - psi * B2 / 2      # g-coefficient after the substitution
        Xi_d = -chi * phi + psi                                                        # dtheta dtheta coefficient
        okt = sp.cancel(Xi_d) == 0 and sp.cancel(Xi_g - th * (kap / (z ** 2 + eps * w ** 2)
                                                             + eps * 92012544 * kap ** 5 * B2 / Dth ** 2)) == 0
        check(f"{SG[eps]} (0,0,z,w): det M = kappa^2 (kappa^2 {'-' if eps == 1 else '+'} 147456 theta^2); Sigma (dt,du) depends"
              f" on d theta only; E^(1) = kappa t g + eps w Xi_0, E^(13) = kappa u g - z Xi_0 with Xi_0 = (576 kappa^3/D)(Hess theta"
              f" - box theta g) + eps (184025088 kappa^5 theta/D^2)(dtheta dtheta - |dtheta|^2 g/2); z E^(1) + eps w E^(13) ="
              f" kappa (z t + eps w u) g; on z t + eps w u = 0 one equation Xi = kappa theta/(z^2 + eps w^2) g + Xi_0 = 0, whose"
              f" traceless part is killed by F' = |kappa^2 - eps 147456 theta^2|^(-13/12), and whose trace is then"
              f" theta (kappa/(z^2 + eps w^2) + eps 92012544 kappa^5 |dtheta|^2/D^2) g", okD and ok and okf and ok1 and okphi and okt)


# ====================================================================== M8
def s_taylor(fd, D):
    """Taylor coefficients (degrees 1..D) of s = (1/4) log det(eta^-1 y) at y = eta + delta, in the monomial basis of
    hessian_system: (1/4) sum_k (-1)^(k+1) tr(M^k)/k, M = eta^-1 delta"""
    Pp, idx = fd['Pp'], fd['idx']
    n = 4
    dl = [[{} for _ in range(n)] for _ in range(n)]
    for A, (i, j) in enumerate(idx):
        dl[i][j] = Pp.add(dl[i][j], Pp.var(A))
        if i != j:
            dl[j][i] = Pp.add(dl[j][i], Pp.var(A))
    eta = fd['eta']
    Mm = [[Pp.scale(dl[i][j], Fr(eta[i])) for j in range(n)] for i in range(n)]
    out, Pk = {}, Mm
    for k in range(1, D + 1):
        tr = Pp.ssum([Pk[i][i] for i in range(n)])
        out = Pp.add(out, Pp.scale(tr, Fr((-1) ** (k + 1), 4 * k)))
        Pk = [[Pp.trunc(Pp.ssum([Pp.mul(Pk[i][l], Mm[l][j], D) for l in range(n)]), D) for j in range(n)] for i in range(n)]
    return out


def m8():
    section("M8. Proposition M9: the first integrability condition of the fibre equation OFF the section (any mu, mu(w) != 0)")
    rng = random.Random(5)
    for eta in ((1, -1, -1, -1), (-1, 1, 1, 1), (1, 1, 1, 1), (1, 1, -1, -1)):
        t1 = time.time()
        R1 = lambda a, b, c, d, eta=eta: Fr((eta[a] if a == c else 0) * (eta[b] if b == d else 0)  # noqa: E731
                                            - (eta[a] if a == d else 0) * (eta[b] if b == c else 0))
        samples = [('random R #1', L.rand_curv(random.Random(90 + sum(eta)))), ('random R #2', L.rand_curv(random.Random(91 + sum(eta)))),
                   ('R = R_1', R1)]
        res, okk = [], True
        for lab, Rl in samples + [('control P = 1 (section)', L.rand_curv(random.Random(74)))]:
            Pm = sp.eye(4) if lab.startswith('control') else sp.Matrix(4, 4, lambda i, j: (1 if i == j else 0) + sp.Rational(rng.randint(-2, 2), 5))
            Pi = Pm.T.inv()
            Bm = {(c, d): sp.Matrix(4, 4, lambda a, b: eta[a] * sp.Rational(Rl(a, b, c, d))) for c in range(4) for d in range(4)}
            Bp = {k: Pm.T * v * Pi for k, v in Bm.items()}       # R_ij at y0 = P eta P^T, pulled back by y -> P y P^T
            Rp = lambda a, b, c, d, Bp=Bp: L.fr(Bp[(c, d)][a, b] * eta[a])  # noqa: E731
            lamF = Fr(1, 2) if -1 in eta else Fr(0)            # the side each real form requires (A4, A5), as in L's L2
            fd = L.fibre_data(eta, Rp, lam=lamF, maxorder=1, afac=Fr(1))
            fd['eta'] = eta
            unk, rowkeys, c0, _, rhs = L.hessian_system(fd, 1, 3, mu0=Fr(1), mu1=Fr(0))
            _, _, c1, _, _ = L.hessian_system(fd, 1, 3, mu0=Fr(1), mu1=Fr(1))
            deg1 = [j for j, mn in enumerate(unk) if sum(mn) == 1]
            rest = [j for j, mn in enumerate(unk) if sum(mn) > 1]
            okk &= all(c0[j] == c1[j] for j in rest)             # mu' enters only through the degree-1 columns
            dcols = [{i: c1[j].get(i, 0) - c0[j].get(i, 0) for i in set(c0[j]) | set(c1[j]) if c1[j].get(i, 0) != c0[j].get(i, 0)}
                     for j in deg1]
            cols = [c0[j] for j in rest] + [c0[j] for j in deg1] + dcols
            r1 = L.dm_from_cols(len(rowkeys), cols).rank()
            r2 = L.dm_from_cols(len(rowkeys), cols, rhs).rank()
            # structural kernel vector: the Taylor coefficients of s (Hess_F s = 0), with v' = 0
            st = s_taylor(fd, 3)
            kv = {j: st.get(mn, Fr(0)) for j, mn in enumerate(unk)}
            colsk = [c0[j] for j in range(len(unk))]
            Av = {}
            for j, col in enumerate(colsk):
                if kv[j]:
                    for i, val in col.items():
                        Av[i] = Av.get(i, 0) + val * kv[j]
            okk &= all(v_ == 0 for v_ in Av.values()) and any(kv.values())
            res.append((lab, r1, r2))
            okk &= (r1 != r2) if not lab.startswith('control') else (r1 == r2)
        check(f"base {(eta.count(1), eta.count(-1))}, lambda = {Fr(1, 2) if -1 in eta else 0}: at y0 = P g(x0) P^T (random rational P) the order-eps Taylor system of"
              f" Hess_F w = (mu0 + mu1 (w - w0)) Bt with mu0 = 1, RELAXED by replacing mu1 dw by a free covector, has"
              f" (rank, augmented rank) = {[(a, b) for _, a, b in res]} for {[l for l, _, _ in res]}: inconsistent at the three"
              f" off-section samples, consistent at the section (control); the Taylor coefficients of s span a kernel of the"
              f" coefficient matrix (Hess_F s = 0), so its rank is <= #unknowns - 1 identically  [{time.time() - t1:.0f} s]", okk)


def main():
    m0()
    m1()
    m2()
    out = m3()
    m4(out)
    m5()
    m6()
    m7()
    m8()
    P(f"\n{NCHK[0]} checks, all True (each asserted).")
    P(f"runtime {time.time() - T0:.0f} s")


if __name__ == '__main__':
    main()
