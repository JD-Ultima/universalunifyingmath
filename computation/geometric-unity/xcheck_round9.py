"""Referee cross-check of GU round 9 (formalization/geometric-unity/M-z2-remaining-members.md), REVIEW7-GU.

Imports only xcheck_round8.py (and, at run time, through it xcheck_round6.py and xcheck_round5.py), unmodified: the
round-5/6 referee construction (own Clifford algebra, Shiab templates from their definitions, the bilinear 2B, the
vector-sector patterns P_k, own grade extraction) supplies the J1/J2 tables S, M, l, b, c1, c13; xcheck_round8 supplies
its own bookkeeping `assemble`, the truncated polynomials `TP`, `rand_curv`, `endos`, `pq_bases` and the order-eps^2
Taylor system `fibre_system` at y = eta.  Nothing from M_remaining.py, L_nonradial.py, J_slaved_sector.py,
I_elliptic_reduction.py, G_symbols.py or F_jets.py is imported.  Written here:
  Z1  Lemma M1: the rotation/boost on the tables at exact rational group elements (symbolic in the member, kappa, t, u),
      and END-TO-END on the assembled coefficients (H, Q, Q', varpi, c) at random exact points; the (7,7) swap
      p -> -K p, (t, u) -> (u, t) end to end, with the sign bookkeeping of g -> -g;
  Z2  Proposition M2(b), M3: varpi on (x, 0, 0, w), the vanishing pattern and the form of E^(13); the mirror (0, y, z, 0);
  Z3  Proposition M4: H, q, q' on u = 0, mu'(tau, c) in closed form, the discriminant, the leading coefficient, that the
      numerator never vanishes identically (gcd over Q[c]), and WHERE mu' is constant in tau: besides c = 0 also at
      c^2 = 11/6 in (14,0) (referee), where H = 1632 x^2/(kappa + 928 t x), mu' = -3/289 and mu is linear in r;
  Z4  Theorem M7 by the DIRECT route along curves (A_k = a_k(rho) d rho): (a) for N = 0 = ell; (b) the traceless AND the
      trace coefficient of E' along the kernel field, y c1 - x c13 = delta (kappa + 208 sigma), v.grad delta (sign with
      L5(a)'s normalisation of v), the line sigma = -kappa/208, F, mu = -(5/12) r;
  Z5  Proposition M8 by the direct route (A_k = a_k(theta) d theta), and the trace identity;
  Z6  Lemma M5 algebra: own (C1) tensors and pairings; q is injective on Sym^2(so(eta)); B_R = -4 R kg^-1 R^T for the
      curvature operator R, so (C1)_0 <=> R^2 = -(c_H/c_V)^2/12 Id; Weyl operators commute with the Hodge star and are
      traceless with tr(R *) = 0 (hence no Ricci-flat solution, even over C), the 36 entries of R_W^2 + Id/12 have rank 12;
      Lambda^2 has signature (2,4) in base (2,2) (no real solution there); an explicit real Lorentzian, non-Ricci-flat
      solution R_L of (C1)_0 at (c_V/c_H)^2 = 1/12;
  Z7  Lemma M5 by brute force (xcheck_round8.fibre_system, order eps^2 at y = eta, mu' = 0): own Weyl tensors in (1,3) and
      (2,2) with the pairing condition met are inconsistent; R_L at (c_V/c_H)^2 = 1/12 is CONSISTENT (observation);
  Z8  Proposition M9 by an own construction: the relaxed order-eps system expanded DIRECTLY at y0 = P eta P^T (not via
      the pull-back y -> P y P^T used by the note), its ranks, the s-kernel, lambda-independence, R = 0 and R_L.
Exact arithmetic only (Fractions, Python/numpy integers, sympy rationals and Gaussian rationals); fixed seeds; one core.
Every check is asserted and printed.  Some checks test the referee's own statements (labelled "observation" or
"referee"); for those, True means the referee's claim was confirmed.  Sections Z0-Z8 are cited as xcheck9 section Zk.
Run (from computation/geometric-unity): python3 xcheck_round9.py > xcheck_round9.out 2>&1   (~3.5 min, 47 checks)
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
from sympy.combinatorics import Permutation  # noqa: E402
from sympy.polys.matrices import DomainMatrix  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import xcheck_round8 as X8  # noqa: E402

T0 = time.time()
RES = []
I = sp.I
x, y, z, w, kap, t, u = X8.x, X8.y, X8.z, X8.w, X8.kap, X8.t, X8.u
cc, tau, th, sgm = sp.symbols('c tau theta sigma')
KS = X8.KS
LAB = {1: '(14,0)', -1: '(7,7)'}
EX = sp.expand
PAIRS = [(a, b) for a in range(4) for b in range(a + 1, 4)]


def say(*a):
    print(*a)
    sys.stdout.flush()


def check(label, ok):
    ok = bool(ok)
    RES.append((label, ok))
    say(f"  [{len(RES):3d}] {label}: {ok}")
    assert ok, label


def section(s):
    say(f"\n=== {s}")


def tick(t1):
    return f"[{time.time() - t1:.0f} s]"


def zero(e):
    return EX(e) == 0


def gp(e, gens):
    return sp.Poly(EX(e), *gens, domain='QQ_I')


# ============================================================================================ helpers on the tables
def subT(T, sub):
    f = lambda e: sp.sympify(e).subs(sub, simultaneous=True)  # noqa: E731
    return {'S': T['S'].applyfunc(f), 'M': T['M'].applyfunc(f), 'c1': EX(f(T['c1'])), 'c13': EX(f(T['c13'])),
            'L': {k: f(v) for k, v in T['L'].items()},
            'B': {k: (tuple(f(a) for a in v) if isinstance(v, tuple) else f(v)) for k, v in T['B'].items()}}


def numerators(Tm):
    """D = det M, NS = -adj(M) S (Sigma = NS/D), and the assembled numerators H D, Q D^2, Q' D^2, varpi D^2
    (xcheck_round8.assemble applied to NS and dN = dNS D - NS dD)"""
    M, S = Tm['M'], Tm['S']
    D = EX(M.det())
    NS = (-M.adjugate() * S).applyfunc(EX)
    dNt = (NS.diff(t) * D - NS * sp.diff(D, t)).applyfunc(EX)
    dNu = (NS.diff(u) * D - NS * sp.diff(D, u)).applyfunc(EX)
    H, Qs, Qp, W = X8.assemble(Tm, NS, dNt, dNu)
    H = {k: [EX(v) for v in vv] for k, vv in H.items()}
    Qs = {k: {ab: EX(v) for ab, v in vv.items()} for k, vv in Qs.items()}
    Qp = {k: {ab: EX(v) for ab, v in vv.items()} for k, vv in Qp.items()}
    W = {k: EX(v) for k, v in W.items()}
    return D, NS, H, Qs, Qp, W


def coeffs_at(Tm, t0, u0):
    """the assembled reduced coefficients (H, Q_sym, Q', varpi, c1, c13) at the exact point (t0, u0), for a table Tm in
    which only t, u are symbolic"""
    M, S = Tm['M'], Tm['S']
    D = EX(M.det())
    NS = (-M.adjugate() * S).applyfunc(EX)
    pt = {t: t0, u: u0}
    D0 = D.subs(pt)
    Sg = NS.subs(pt) / D0
    dSt = (NS.diff(t).subs(pt) * D0 - NS.subs(pt) * sp.diff(D, t).subs(pt)) / D0 ** 2
    dSu = (NS.diff(u).subs(pt) * D0 - NS.subs(pt) * sp.diff(D, u).subs(pt)) / D0 ** 2
    H, Qs, Qp, W = X8.assemble(Tm, Sg, dSt, dSu)
    sim = lambda e: sp.expand(sp.sympify(e))  # noqa: E731
    Hm = sp.Matrix(2, 2, lambda i, a: sim(H[('E1', 'E13')[i]][a]))
    Qm = {nm: sp.Matrix(2, 2, lambda a, b: sim(Qs[nm][(a, b)])) for nm in ('E1', 'E13')}
    Qpm = {nm: sp.Matrix(2, 2, lambda a, b: sim(Qp[nm][(a, b)])) for nm in ('E1', 'E13')}
    Wd = {k: sim(v) for k, v in W.items()}
    return dict(D=D0, H=Hm, Q=Qm, Qp=Qpm, W=Wd, c=(sim(Tm['c1'].subs(pt)), sim(Tm['c13'].subs(pt))))


def rq(rng, lo=-5, hi=5, den=4):
    while True:
        v = sp.Rational(rng.randint(lo, hi), rng.randint(1, den))
        if v != 0:
            return v


# ============================================================================================ Z0 inputs
def z0():
    section("Z0  Inputs: S, M, l, b, c1, c13 from the referee construction (xcheck_round8 X0 -> xcheck_round6 X0-X2, own "
            "Clifford algebra); not from M_remaining.py, L_nonradial.py or J_slaved_sector.py")
    t1 = time.time()
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        tabs = X8.x0()
    ok = len(X8.RES) == 1 and X8.RES[0][1] and ': True' in buf.getvalue() and ': False' not in buf.getvalue()
    check(f"xcheck_round8.x0() rerun as a library: its check (38 internal checks of xcheck_round6: J1 on all 14 covectors, "
          f"J2 on all 196 tensors and sampled pairs, tables equal J's printed tables) passes, both signatures {tick(t1)}", ok)
    for ep in (1, -1):
        say(f"     {LAB[ep]}: c1 = {tabs[ep]['c1']};  c13 = {tabs[ep]['c13']}")
    return tabs


# ============================================================================================ Z1 Lemma M1
def z1(tabs):
    section("Z1  Lemma M1: the rotation (14,0) / boost (7,7) and the (7,7) swap, on the tables and end to end")
    for ep in (1, -1):
        t1 = time.time()
        T = tabs[ep]
        grp = ((sp.Rational(3, 5), sp.Rational(4, 5)), (sp.Rational(-5, 13), sp.Rational(12, 13))) if ep == 1 else \
            ((sp.Rational(5, 4), sp.Rational(3, 4)), (sp.Rational(13, 5), sp.Rational(-12, 5)))
        ok = True
        for c_, s_ in grp:
            R = sp.Matrix([[c_, -ep * s_], [s_, c_]])
            ok &= R.T * sp.diag(1, ep) * R == sp.diag(1, ep) and R.det() == 1
            rot = lambda a, b: tuple(R * sp.Matrix([a, b]))  # noqa: E731
            (xn, yn), (zn, wn), (tn, un) = rot(x, y), rot(z, w), rot(t, u)
            sub = {x: xn, y: yn, z: zn, w: wn, t: tn, u: un}
            f = lambda e: EX(sp.sympify(e).subs(sub, simultaneous=True))  # noqa: E731
            ok &= (T['M'].applyfunc(f) - T['M']).applyfunc(EX) == sp.zeros(4, 4)
            ok &= (T['S'].applyfunc(f) * R - T['S']).applyfunc(EX) == sp.zeros(4, 2)
            for k in KS:
                for a_, b_ in (('E1', 'E13'), ('w3', 'w11')):
                    vv = sp.Matrix([T['L'][(a_, k)], T['L'][(b_, k)]])
                    ok &= (vv.applyfunc(f) - R * vv).applyfunc(EX) == sp.zeros(2, 1)
            for k, l in itertools.combinations_with_replacement(KS, 2):
                for j in range(3):
                    vv = sp.Matrix([T['B'][('E1', k, l)][j], T['B'][('E13', k, l)][j]])
                    ok &= (vv.applyfunc(f) - R * vv).applyfunc(EX) == sp.zeros(2, 1)
            for k, l in itertools.combinations(KS, 2):
                vv = sp.Matrix([T['B'][('w3', k, l)], T['B'][('w11', k, l)]])
                ok &= (vv.applyfunc(f) - R * vv).applyfunc(EX) == sp.zeros(2, 1)
            vv = sp.Matrix([T['c1'], T['c13']])
            ok &= (vv.applyfunc(f) - R * vv).applyfunc(EX) == sp.zeros(2, 1)
        check(f"{LAB[ep]} M1(a) on the referee tables at two exact group elements {[(str(a), str(b)) for a, b in grp]} "
              f"(symbolic x, y, z, w, kappa, t, u): M(Rp, R tau) = M(p, tau), S(Rp) R = S(p), and the columns of l (E1/E13, "
              f"w3/w11), b and (c1, c13) transform by R {tick(t1)}", ok)
        # end to end: the assembled coefficients at a random exact point
        t1 = time.time()
        rng = random.Random(20261131 + ep)
        okE = True
        for _ in range(2):
            p = [rq(rng) for _ in range(4)]
            kv = sp.Rational(rng.randint(2, 9), rng.randint(1, 3))
            t0, u0 = sp.Rational(rng.randint(-3, 3), 50), sp.Rational(rng.randint(-3, 3), 50)
            c_, s_ = grp[0]
            R = sp.Matrix([[c_, -ep * s_], [s_, c_]])
            pr = list(R * sp.Matrix(p[:2])) + list(R * sp.Matrix(p[2:]))
            tr_ = R * sp.Matrix([t0, u0])
            A = coeffs_at(subT(T, {x: p[0], y: p[1], z: p[2], w: p[3], kap: kv}), t0, u0)
            B = coeffs_at(subT(T, {x: pr[0], y: pr[1], z: pr[2], w: pr[3], kap: kv}), tr_[0], tr_[1])
            okE &= A['D'] == B['D'] and A['D'] != 0
            okE &= (B['H'] * R - R * A['H']).applyfunc(sp.simplify) == sp.zeros(2, 2)
            for key in ('Q', 'Qp'):
                for i, nm in enumerate(('E1', 'E13')):
                    lhs = R.T * B[key][nm] * R
                    rhs = R[i, 0] * A[key]['E1'] + R[i, 1] * A[key]['E13']
                    okE &= (lhs - rhs).applyfunc(sp.simplify) == sp.zeros(2, 2)
            for a_, b_ in (('E1', 'E13'), ('w3', 'w11')):
                okE &= (sp.Matrix([B['W'][a_], B['W'][b_]]) - R * sp.Matrix([A['W'][a_], A['W'][b_]])).applyfunc(sp.simplify) \
                    == sp.zeros(2, 1)
            okE &= (sp.Matrix(B['c']) - R * sp.Matrix(A['c'])).applyfunc(sp.simplify) == sp.zeros(2, 1)
            okE &= any(v != 0 for v in A['W'].values()) and A['H'] != sp.zeros(2, 2)
        check(f"{LAB[ep]} M1(a) end to end (referee's own assembly of H, Q_sym, Q', varpi from Sigma = -M^-1 S, two random "
              f"exact members and points): det M equal, H' R = R H, R^T Q'^(i) R = sum_j R_ij Q^(j) (Q_sym and Q'), "
              f"(varpi1, varpi13)' = R (varpi1, varpi13), (varpi3, varpi11)' = R (varpi3, varpi11), c' = R c; so "
              f"(t, u) solves the reduced system for p iff R(t, u) does for Rp {tick(t1)}", okE)
    # (7,7) swap
    t1 = time.time()
    T = tabs[-1]
    Ds = EX(T['M'].det())
    okD = zero(Ds.subs({x: -y, y: -x, z: -w, w: -z, t: u, u: t}, simultaneous=True) - Ds)
    oks = zero(T['c1'].subs({x: -y, y: -x, z: -w, w: -z, t: u, u: t}, simultaneous=True) - T['c13']) and \
        zero(T['c13'].subs({x: -y, y: -x, z: -w, w: -z, t: u, u: t}, simultaneous=True) - T['c1'])
    rng = random.Random(20261133)
    okE = True
    sw = lambda Mx: sp.Matrix(2, 2, lambda a, b: Mx[1 - a, 1 - b])  # noqa: E731
    for _ in range(3):
        p = [rq(rng) for _ in range(4)]
        kv = sp.Rational(rng.randint(2, 9), rng.randint(1, 3))
        t0, u0 = sp.Rational(rng.randint(-3, 3), 50), sp.Rational(rng.randint(-3, 3), 50)
        A = coeffs_at(subT(T, {x: p[0], y: p[1], z: p[2], w: p[3], kap: kv}), t0, u0)
        B = coeffs_at(subT(T, {x: -p[1], y: -p[0], z: -p[3], w: -p[2], kap: kv}), u0, t0)
        okE &= A['D'] == B['D']
        okE &= all(B['H'][0, a] == -A['H'][1, 1 - a] and B['H'][1, a] == -A['H'][0, 1 - a] for a in range(2))
        for key in ('Q', 'Qp'):
            okE &= B[key]['E1'] == -sw(A[key]['E13']) and B[key]['E13'] == -sw(A[key]['E1'])
        okE &= B['W']['E1'] == A['W']['E13'] and B['W']['E13'] == A['W']['E1'] and B['W']['w3'] == A['W']['w11'] and \
            B['W']['w11'] == A['W']['w3']
        okE &= B['c'] == (A['c'][1], A['c'][0])
        # formal bookkeeping: E^(i) as coefficient vectors on the basis
        # [G, g, Hs_t, Hs_u, dt dt, dt du, du dt, du du, |dt|^2 g, <dt,du> g, |du|^2 g, dt^du]; g -> -g flips only c_i g
        xs = (p[0], p[1])

        def vec(C, i, xi, gsign, swap):
            nm = ('E1', 'E13')[i]
            Hr = [C['H'][i, 0], C['H'][i, 1]]
            Q, Qp_ = C['Q'][nm], C['Qp'][nm]
            v = [-xi, gsign * C['c'][i], Hr[0], Hr[1], Q[0, 0], Q[0, 1], Q[1, 0], Q[1, 1], Qp_[0, 0], Qp_[0, 1] + Qp_[1, 0],
                 Qp_[1, 1], C['W'][nm]]
            if swap:          # rewrite in the (t, u) basis of p when the primed scalars are (t', u') = (u, t)
                v = [v[0], v[1], v[3], v[2], v[7], v[6], v[5], v[4], v[10], v[9], v[8], -v[11]]
            return v
        xp = (-p[1], -p[0])
        okE &= vec(B, 0, xp[0], -1, True) == [-e for e in vec(A, 1, xs[1], 1, False)]
        okE &= vec(B, 1, xp[1], -1, True) == [-e for e in vec(A, 0, xs[0], 1, False)]
    check(f"(7,7) M1(b) swap p' = -K p, (t', u') = (u, t): det M invariant (symbolic), c'_1 = c_13, c'_13 = c_1 (symbolic); "
          f"end to end at three random exact points: H'_(1,phi') = -H_(13,K phi), H'_(13,phi') = -H_(1,K phi), Q, Q' likewise, "
          f"varpi'_(1,13,3,11) = varpi_(13,1,11,3); and, as coefficient vectors on (G, g, Hess-terms, gradient terms, "
          f"dt^du) with g -> -g flipping only c_i g: E'^(1)[-g] = -E^(13)[g], E'^(13)[-g] = -E^(1)[g] {tick(t1)}",
          okD and oks and okE)


# ============================================================================================ Z2 M2(b), M3
def z2(tabs):
    section("Z2  Proposition M2(b) and M3: (x, 0, 0, w) symbolic in x, w, kappa, t, u (and the mirror (0, y, z, 0) in (7,7))")
    out = {}
    for ep in (1, -1):
        t1 = time.time()
        Tm = subT(tabs[ep], {y: 0, z: 0})
        D, NS, H, Qs, Qp, W = numerators(Tm)
        out[ep] = (Tm, D, NS, H, Qs, Qp, W)
        ok = all(NS[k, 1] == 0 for k in range(4))
        ok &= H['E13'] == [0, 0] and H['E1'][1] == 0 and H['E1'][0] != 0
        ok &= all(Qs['E13'][ab] == 0 and Qp['E13'][ab] == 0 for ab in ((0, 1), (1, 0), (1, 1)))
        ok &= zero(Qp['E13'][(0, 0)] + Qs['E13'][(0, 0)]) and W['E13'] == 0
        ok &= zero(Tm['c13'] - u * (kap + 208 * x * t))
        q13 = Qs['E13'][(0, 0)]
        okw = zero(W['E1'] + ep * q13)
        qq, rr = sp.div(gp(q13, (x, w, kap, t, u)), gp(kap * u * x ** 2, (x, w, kap, t, u)))
        okw &= rr.is_zero and not qq.is_zero
        check(f"{LAB[ep]} M3: Sigma has no du column; H_13 = 0, H_1u = 0; Q^(13) has only the dt dt entry and Q'^(13)_tt = "
              f"-Q^(13)_tt; varpi^(13) = 0; c13 = u (kappa + 208 x t); so E^(13) = u (kappa + 208 x t) g + Q^(13)_tt "
              f"(dt dt - |dt|^2 g) {tick(t1)}", ok)
        check(f"{LAB[ep]} M2(b): varpi^(1) = {'-' if ep == 1 else '+'}Q^(13)_tt, and Q^(13)_tt (det M)^2 = kappa u x^2 P(x, w, "
              f"kappa, t, u) with P a nonzero polynomial", okw)
    # mirror in (7,7)
    t1 = time.time()
    Tm = subT(tabs[-1], {x: 0, w: 0})
    D, NS, H, Qs, Qp, W = numerators(Tm)
    ok = all(NS[k, 0] == 0 for k in range(4)) and H['E1'] == [0, 0] and H['E13'][0] == 0
    ok &= all(Qs['E1'][ab] == 0 and Qp['E1'][ab] == 0 for ab in ((0, 0), (0, 1), (1, 0)))
    ok &= zero(Qp['E1'][(1, 1)] + Qs['E1'][(1, 1)]) and W['E1'] == 0 and zero(Tm['c1'] - t * (kap - 208 * y * u))
    check(f"(7,7) mirror (0, y, z, 0): Sigma has no dt column, H_1 = 0, H_13t = 0, varpi^(1) = 0, E^(1) = t (kappa - 208 y u) g "
          f"+ Q^(1)_uu (du du - |du|^2 g) {tick(t1)}", ok)
    return out


# ============================================================================================ Z3 M4
def z3(out):
    section("Z3  Proposition M4: the single equation on u = 0 for (x, 0, 0, w); mu'(tau, c), tau = x t/kappa, c = w/x")
    for ep in (1, -1):
        t1 = time.time()
        Tm, D, NS, H, Qs, Qp, W = out[ep]
        D0 = D.subs(u, 0)
        H0 = sp.cancel(H['E1'][0].subs(u, 0) / D0)
        q0 = sp.cancel(Qs['E1'][(0, 0)].subs(u, 0) / D0 ** 2)
        qp0 = sp.cancel(Qp['E1'][(0, 0)].subs(u, 0) / D0 ** 2)
        d1 = kap ** 2 + 768 * kap * t * x + 121856 * x ** 2 * t ** 2 - ep * 147456 * w ** 2 * t ** 2
        ok = zero(D0 - kap ** 2 * d1)
        ok &= sp.cancel(H0 - 576 * (kap * (x ** 2 + ep * w ** 2) + 544 * t * x * (x ** 2 - ep * w ** 2)) / d1) == 0
        ok &= zero(Tm['c1'].subs(u, 0) - (kap * t + 312 * x * t ** 2))
        mup = sp.cancel(x * (q0 - sp.diff(H0, t)) / H0 ** 2)
        mt = sp.factor(sp.cancel(mup.subs({t: tau * kap / x, w: cc * x})))
        num = 7 - ep * 3 * cc ** 2 + (7616 - ep * 3008 * cc ** 2 - 768 * cc ** 4) * tau + \
            (2071552 - ep * 3226624 * cc ** 2 + 1155072 * cc ** 4) * tau ** 2
        den = 18 * ((1 + ep * cc ** 2) + 544 * (1 - ep * cc ** 2) * tau) ** 2
        ok &= sp.cancel(mt + num / den) == 0 and sp.cancel(mt.subs(cc, 0) + sp.Rational(7, 18)) == 0
        A2, A1, A0 = sp.Poly(num, tau).all_coeffs()
        disc = sp.factor(A1 ** 2 - 4 * A2 * A0)
        okd = zero(disc - 16384 * cc ** 2 * (cc ** 2 + ep * 35) * (6 * cc ** 2 - ep * 11) ** 2)
        okd &= zero(A2 - 1024 * (cc ** 2 - ep) * (1128 * cc ** 2 - ep * 2023))
        g3 = sp.gcd(sp.gcd(sp.Poly(A0, cc), sp.Poly(A1, cc)), sp.Poly(A2, cc))
        okd &= g3.degree() == 0
        # where is mu' constant in tau?  numerator of d mu'/d tau: gcd over Q[c] of its tau-coefficients
        dn = sp.Poly(sp.numer(sp.together(sp.diff(num / den, tau))), tau)
        gg = sp.Poly(dn.all_coeffs()[0], cc)
        for cf in dn.all_coeffs()[1:]:
            gg = sp.gcd(gg, sp.Poly(cf, cc))
        okg = sp.factor(gg.as_expr() / gg.LC()) == sp.factor(cc ** 2 * (cc ** 2 - ep * sp.Rational(11, 6)))
        tH = (1 + ep * cc ** 2) / (544 * (ep * cc ** 2 - 1))
        okd &= sp.cancel(((1 + ep * cc ** 2) + 544 * (1 - ep * cc ** 2) * tau).subs(tau, tH)) == 0
        hn = (1 + ep * cc ** 2) + 544 * (1 - ep * cc ** 2) * tau
        dm = 1 + 768 * tau + (121856 - 147456 * ep * cc ** 2) * tau ** 2
        rs = sp.factor(sp.resultant(hn, dm, tau))
        okr = zero(rs - (-ep) * 4096 * cc ** 2 * (6 * cc ** 2 - ep * 11) ** 2)
        check(f"{LAB[ep]} M4: discriminant of N in tau = 16384 c^2 (c^2 {'+' if ep == 1 else '-'} 35)(6 c^2 "
              f"{'-' if ep == 1 else '+'} 11)^2, leading coefficient 1024 (c^2 {'-' if ep == 1 else '+'} 1)(1128 c^2 "
              f"{'-' if ep == 1 else '+'} 2023); the three tau-coefficients of N have no common root c in C (gcd 1), so mu' "
              f"is not identically 0 for any c; the numerator of H vanishes at tau_H = (1 + eps c^2)/(544 (eps c^2 - 1)) {tick(t1)}",
              ok and okd)
        tail = ("(14,0) at c^2 = 11/6, against M4's 'c != 0: mu' is a nonconstant rational function of tau'" if ep == 1
                else "(14,0); in (7,7) never")
        check(f"{LAB[ep]} M4, referee: d mu'/d tau vanishes identically exactly for c^2 (c^2 - eps 11/6) = 0, and the numerator "
              f"of H and det M|_(u=0) have a common root exactly there (resultant {rs}); so for REAL c != 0 mu' is constant in "
              f"tau only in {tail}", okg and okr)
        if ep == 1:
            t2 = time.time()
            c6 = sp.sqrt(sp.Rational(66)) / 6                      # c^2 = 11/6
            dd = kap + 928 * t * x
            D6 = sp.factor(D0.subs(w, c6 * x))
            H6 = sp.factor(sp.cancel(H0.subs(w, c6 * x)))
            q6 = sp.factor(sp.cancel(q0.subs(w, c6 * x)))
            m6 = sp.simplify(sp.cancel(mup.subs(w, c6 * x)))
            ok6 = zero(D6 - kap ** 2 * (kap - 160 * t * x) * dd) and sp.cancel(H6 - 1632 * x ** 2 / dd) == 0
            ok6 &= sp.cancel(q6 + 1542144 * x ** 3 / dd ** 2) == 0 and m6 == sp.Rational(-3, 289)
            # F' = |d|^(-502/493): F''/F' = q/H; F = |d|^(-9/493)/(-(9/493) 928 x) (nonzero), mu = x F'/H = -(3/289) F
            ok6 &= sp.cancel(q6 / H6 + sp.Rational(502, 493) * sp.diff(dd, t) / dd) == 0
            FpF = -sp.Rational(9, 493) * sp.diff(dd, t) / dd                       # F'/F
            ok6 &= sp.cancel(x / H6 * FpF + sp.Rational(3, 289)) == 0
            check(f"(14,0) referee, the member x (W_A + c W'_C Ad_Gamma) with c^2 = 11/6 (from the referee tables, not the closed "
                  f"form): det M|_(u=0) = kappa^2 (kappa - 160 t x)(kappa + 928 t x), H = 1632 x^2/(kappa + 928 t x) (no zero "
                  f"where det M != 0), q = -1542144 x^3/(kappa + 928 t x)^2, mu' = {m6}; with F' = |kappa + 928 t x|^(-502/493), "
                  f"r = F(t) is proportional to |kappa + 928 t x|^(-9/493), defined and nonzero on all of U, and mu = x F'/H = "
                  f"-(3/289) r EXACTLY: the remaining equation is linear, so L6 and L9 (linear case of L7) apply as in M7(b) "
                  f"step 6: no solution on any open U for this member {tick(t2)}", ok6)
        say(f"     {LAB[ep]} on u = 0 (for the record; the note says q, q' are printed in output section M4, which they are "
            f"not): q = {sp.factor(q0)};  q' = {sp.factor(qp0)}")


# ============================================================================================ Z4 M7
def direct(Tm, D, NS, vt, vu, lam, beta12, beta3):
    """direct route along an integral curve of v = (vt, vu): a_k = (NS v)_k/D; returns numerators (times D^2) of the
    traceless (d rho d rho) coefficient and of the |d rho|^2 g coefficient of sum_k lam_k (N_k^T - tr N_k g) + quadratic"""
    nk = [EX(NS[k, 0] * vt + NS[k, 1] * vu) for k in range(4)]

    def vg(e):
        return EX(vt * sp.diff(e, t) + vu * sp.diff(e, u))
    vD = vg(D)
    hv = EX(sum(lam[k] * nk[k] for k in range(4)))
    lin = sum(lam[k] * (vg(nk[k]) * D - nk[k] * vD) for k in range(4) if lam[k] != 0)
    q12, q3 = 0, 0
    for (ik, k), (il, ll) in itertools.combinations_with_replacement(list(enumerate(KS)), 2):
        wkl = sp.Rational(1, 2) if k == ll else 1
        if beta12[(k, ll)] != 0:
            q12 += wkl * beta12[(k, ll)] * nk[ik] * nk[il]
        if beta3[(k, ll)] != 0:
            q3 += wkl * beta3[(k, ll)] * nk[ik] * nk[il]
    return hv, EX(lin + q12), EX(-lin + q3)


def z4(tabs):
    section("Z4  Theorem M7 (N = 0 in (7,7), y = sgn x) by the direct route along curves")
    ep = -1
    T = tabs[ep]
    for sg in (1, -1):
        t1 = time.time()
        Tm = subT(T, {y: sg * x, z: sg * cc * x, w: cc * x})
        D, NS, H, Qs, Qp, W = numerators(Tm)
        ok = all(zero(NS[k, 1] + sg * NS[k, 0]) for k in range(4))
        ok &= all(zero(H['E1'][a] - sg * H['E13'][a]) for a in range(2))
        ok &= all(zero(Qs['E1'][ab] - sg * Qs['E13'][ab]) and zero(Qp['E1'][ab] - sg * Qp['E13'][ab]) for ab in Qs['E1'])
        sig = x * (t - sg * u)
        ok &= zero(x * (Tm['c1'] - sg * Tm['c13']) - sig * (kap + 208 * sig))
        ok &= all(v == 0 for v in W.values())
        check(f"(7,7) M7(a) (x, {'+' if sg == 1 else '-'}x, {'+' if sg == 1 else '-'}cx, cx): Sigma_u = -sgn Sigma_t (T_e = "
              f"Sigma_t d sigma/x), H, Q, Q' of E^(1) equal sgn times those of E^(13), x (c1 - sgn c13) = sigma (kappa + 208 "
              f"sigma) with sigma = x (t - sgn u), all four varpi = 0 {tick(t1)}", ok)
    for sg in (1, -1):
        t1 = time.time()
        Tm = subT(T, {y: sg * x})
        D, NS, H, Qs, Qp, W = numerators(Tm)
        ell = x * z - sg * x * w
        gens = (x, z, w, kap, t, u)
        hN = [EX(sg * x * H['E1'][a] - x * H['E13'][a]) for a in range(2)]       # h D, h = (y, -x) H
        vt_, rt = sp.div(gp(hN[1], gens), gp(576 * ell, gens))
        vu_, ru = sp.div(gp(-hN[0], gens), gp(576 * ell, gens))
        okp = rt.is_zero and ru.is_zero
        vt, vu = vt_.as_expr(), vu_.as_expr()
        okv = not (vt_.is_zero and vu_.is_zero)
        l1, l13 = X8.lrow(Tm, 'E1'), X8.lrow(Tm, 'E13')
        lam = [EX(sg * x * l1[k] - x * l13[k]) for k in range(4)]
        b12, b3 = {}, {}
        for k, ll in itertools.combinations_with_replacement(KS, 2):
            p1, p13 = Tm['B'][('E1', k, ll)], Tm['B'][('E13', k, ll)]
            b12[(k, ll)] = EX(sg * x * (p1[0] + p1[1]) - x * (p13[0] + p13[1]))
            b3[(k, ll)] = EX(sg * x * p1[2] - x * p13[2])
        hv, Rtl, Ktr = direct(Tm, D, NS, vt, vu, lam, b12, b3)
        de, sig = sg * x * t - u * x, x * t - sg * x * u
        okc = zero(sg * x * Tm['c1'] - x * Tm['c13'] - de * (kap + 208 * sig)) and zero(de - sg * sig)
        inv = EX(vt * sp.diff(de, t) + vu * sp.diff(de, u))
        ratio = sp.cancel(inv / (kap * ell * (kap + 16 * sig) * (kap + 208 * sig)))
        check(f"(7,7) M7(b) (x, {'+' if sg == 1 else '-'}x, z, w), ell = x (z {'-' if sg == 1 else '+'} w): v = det M (h_u, -h_t)/"
              f"(576 ell) (L5(a)'s normalisation) is a nonzero polynomial ({okp and okv}); h.v = 0; by the DIRECT route along "
              f"integral curves of v the (d rho d rho)_0 coefficient R' and the |d rho|^2 g coefficient of E' both vanish "
              f"identically, so E' = (y c1 - x c13) g = delta (kappa + 208 sigma) g with delta = sgn sigma {tick(t1)}",
              okp and okv and hv == 0 and Rtl == 0 and Ktr == 0 and okc)
        check(f"(7,7) M7(b), sgn = {sg:+d}, referee: with L5(a)'s v, v.grad delta = {ratio} kappa ell (kappa + 16 sigma)(kappa + "
              f"208 sigma), i.e. -sgn kappa ell (..)(..); so the note's '-kappa ell (..)(..)' holds for sgn = +1 and has the "
              f"opposite sign for sgn = -1 (zero sets, hence the argument, unchanged)", ratio == -sg)
        # the line sigma = -kappa/208: t = sgn u - kappa/(208 x), dt = sgn du; direct route with A_k = a_k(u) du
        t1 = time.time()
        line = {t: sg * u - kap / (208 * x)}
        Dl = sp.factor(D.subs(line))
        zeta = kap * w + 208 * u * ell
        okD = sp.cancel(Dl - 576 * kap ** 2 * zeta ** 2 / (169 * x ** 2)) == 0
        a = [sp.cancel((sg * NS[k, 0] + NS[k, 1]).subs(line) / Dl) for k in range(4)]
        ap = [sp.cancel(sp.diff(a[k], u)) for k in range(4)]
        co = {}
        for nm in ('E1', 'E13'):
            lr = X8.lrow(Tm, nm)
            hc = sp.cancel(sum(lr[k] * a[k] for k in range(4)))
            qc = sum(lr[k] * ap[k] for k in range(4))
            qpc = -sum(lr[k] * ap[k] for k in range(4))
            for (ik, k), (il, ll) in itertools.combinations_with_replacement(list(enumerate(KS)), 2):
                b1, b2, b3_ = Tm['B'][(nm, k, ll)]
                wkl = sp.Rational(1, 2) if k == ll else 1
                qc += wkl * (b1 + b2) * a[ik] * a[il]
                qpc += wkl * b3_ * a[ik] * a[il]
            co[nm] = (hc, sp.cancel(qc), sp.cancel(qpc))
        hc, qc, qpc = co['E1']
        okl = sp.cancel(hc - 312 * x * ell / zeta) == 0 and sp.cancel(qc + 105456 * x * ell ** 2 / zeta ** 2) == 0 and \
            sp.cancel(qpc - 52728 * x * ell ** 2 / zeta ** 2) == 0
        okl &= all(sp.cancel(p_ - sg * q_) == 0 for p_, q_ in zip(co['E13'], co['E1']))
        okl &= zero(Tm['c13'].subs(line) - sg * Tm['c1'].subs(line))
        zu = sp.diff(zeta, u)
        okF = sp.cancel(qc / hc + sp.Rational(13, 8) * zu / zeta) == 0               # F' = |zeta|^(-13/8)
        okF &= sp.cancel(x / hc * (-sp.Rational(5, 8) * zu / zeta) + sp.Rational(5, 12)) == 0   # mu/F, F'/F = -(5/8) zeta_u/zeta
        mup = sp.cancel(x * (qc - sp.diff(hc, u)) / hc ** 2)
        check(f"(7,7) M7(b) on the line sigma = -kappa/208 (direct route, A_k = a_k(u) du): det M = 576 kappa^2 zeta^2/(169 x^2), "
              f"zeta = kappa w + 208 ell u; E^(1) has Hessian coefficient 312 x ell/zeta, gradient coefficients "
              f"-105456 x ell^2/zeta^2 and 52728 x ell^2/zeta^2; E^(13) = sgn E^(1) (coefficients and c's); F' = |zeta|^(-13/8) "
              f"solves F''/F' = q/h, mu = x F'/h = -(5/12) F, mu' = {mup} {tick(t1)}", okD and okl and okF and mup == sp.Rational(-5, 12))


# ============================================================================================ Z5 M8
def z5(tabs):
    section("Z5  Proposition M8 ((x, y) = 0) by the direct route, A_k = a_k(theta) d theta, theta = t w - u z")
    for ep in (1, -1):
        t1 = time.time()
        Tm = subT(tabs[ep], {x: 0, y: 0})
        M, S = Tm['M'], Tm['S']
        D = EX(M.det())
        NS = (-M.adjugate() * S).applyfunc(EX)
        okD = zero(D - kap ** 2 * (kap ** 2 - ep * 147456 * (t * w - u * z) ** 2))
        ok = all(zero(NS[k, 0] * z + NS[k, 1] * w) for k in range(4))
        sub_th = {t: (th + u * z) / w}
        a = [sp.cancel((NS[k, 0] / (w * D)).subs(sub_th)) for k in range(4)]
        ok &= all(sp.diff(ak, u) == 0 for ak in a)
        Dth = kap ** 2 * (kap ** 2 - ep * 147456 * th ** 2)
        ap = [sp.cancel(sp.diff(ak, th)) for ak in a]
        rowc = {'E1': ep * w, 'E13': -z}
        for nm in ('E1', 'E13'):
            lr = X8.lrow(Tm, nm)
            hc = sp.cancel(sum(lr[k] * a[k] for k in range(4)))
            qc = sum(lr[k] * ap[k] for k in range(4))
            qpc = -sum(lr[k] * ap[k] for k in range(4))
            for (ik, k), (il, ll) in itertools.combinations_with_replacement(list(enumerate(KS)), 2):
                b1, b2, b3 = Tm['B'][(nm, k, ll)]
                wkl = sp.Rational(1, 2) if k == ll else 1
                qc += wkl * (b1 + b2) * a[ik] * a[il]
                qpc += wkl * b3 * a[ik] * a[il]
            ok &= sp.cancel(hc - 576 * kap ** 3 * rowc[nm] / Dth) == 0
            ok &= sp.cancel(qc - ep * 184025088 * kap ** 5 * th * rowc[nm] / Dth ** 2) == 0
            ok &= sp.cancel(qpc + ep * 92012544 * kap ** 5 * th * rowc[nm] / Dth ** 2) == 0
        ok &= zero(Tm['c1'] - kap * t) and zero(Tm['c13'] - kap * u)
        # derivative-free combination and the reduced equation
        okc = zero(z * rowc['E1'] + ep * w * rowc['E13']) and zero(z * Tm['c1'] + ep * w * Tm['c13'] - kap * (z * t + ep * w * u))
        tu = (th * ep * w / (z ** 2 + ep * w ** 2), -th * z / (z ** 2 + ep * w ** 2))
        okc &= zero(sp.cancel(z * tu[0] + ep * w * tu[1])) and sp.cancel(tu[0] * w - tu[1] * z - th) == 0
        chi, psi = 576 * kap ** 3 / Dth, ep * 184025088 * kap ** 5 * th / Dth ** 2
        phi = sp.cancel(psi / chi)
        okc &= sp.cancel(phi + sp.Rational(13, 12) * sp.diff(sp.log(kap ** 2 - ep * 147456 * th ** 2), th)) == 0
        B2, C_ = sp.symbols('B2 C')
        # with Hess theta = -phi dth dth (nabla dr = 0): Xi = kappa theta/(z^2 + eps w^2) g + chi (Hess - box g) + psi (dd - B2 g/2)
        Xi_g = kap * th / (z ** 2 + ep * w ** 2) + chi * phi * B2 - psi * B2 / 2
        Xi_d = -chi * phi + psi
        okc &= sp.cancel(Xi_d) == 0 and sp.cancel(Xi_g - th * (kap / (z ** 2 + ep * w ** 2) + ep * 92012544 * kap ** 5 * B2 / Dth ** 2)) == 0
        # |dth|^2 = C |Delta|^(13/6): the bracket becomes kappa/(z^2 + eps w^2) + eps 92012544 kappa C |Delta|^(1/6)
        Dl = sp.Symbol('Delta', positive=True)
        okc &= sp.simplify(92012544 * kap ** 5 * C_ * Dl ** sp.Rational(13, 6) / (kap ** 2 * Dl) ** 2 - 92012544 * kap * C_ * Dl ** sp.Rational(1, 6)) == 0
        check(f"{LAB[ep]} M8 (referee tables, direct route): det M = kappa^2 (kappa^2 {'-' if ep == 1 else '+'} 147456 theta^2); "
              f"the columns of Sigma are proportional to (w, -z) and a_k = Sigma_kt/w depends on theta only; E^(1) = kappa t g "
              f"+ eps w Xi_0, E^(13) = kappa u g - z Xi_0 with Xi_0 = (576 kappa^3/D)(Hess theta - box theta g) + eps (184025088 "
              f"kappa^5 theta/D^2)(dtheta dtheta - |dtheta|^2 g/2) {tick(t1)}", okD and ok)
        check(f"{LAB[ep]} M8: z E^(1) + eps w E^(13) = kappa (z t + eps w u) g; on z t + eps w u = 0, (t, u) = theta (eps w, -z)/"
              f"(z^2 + eps w^2); F''/F' = psi/chi = -(13/12) d/dtheta log|kappa^2 - eps 147456 theta^2|; with nabla dr = 0, "
              f"Xi = theta (kappa/(z^2 + eps w^2) + eps 92012544 kappa^5 |dtheta|^2/D^2) g, and with |dtheta|^2 = C |Delta|^(13/6) "
              f"the bracket is kappa/(z^2 + eps w^2) + eps 92012544 kappa C |Delta|^(1/6)", okc)


# ============================================================================================ Z6 M5 algebra
def comm(A, B):
    return np.einsum('...ij,...jk->...ik', A, B) - np.einsum('...ij,...jk->...ik', B, A)


def c0t(Xs):            # 4 c0(X_a, X_b, X_c) = -[[X_b, X_c], X_a]
    BC = comm(Xs[None, :, None], Xs[None, None, :])
    return -comm(BC, Xs[:, None, None])


def c1t(Xs):            # c1 = tr(X_a X_b) X_c - tr(X_a X_c) X_b
    g = np.einsum('aij,bji->ab', Xs, Xs)
    return np.einsum('ab,cij->abcij', g, Xs) - np.einsum('ac,bij->abcij', g, Xs)


def qt(Xs, K, Kp):      # 4 q(K (x) K'), 4 d(X, Y) Z = tr(Z [K, X]) [K', Y] + tr(Z [K, Y]) [K', X]
    P = comm(K[None], Xs)
    Q = comm(Kp[None], Xs)
    tauv = np.einsum('cij,aji->ca', Xs, P)
    dd = np.einsum('ca,bij->abcij', tauv, Q) + np.einsum('cb,aij->abcij', tauv, Q)
    return dd - np.transpose(dd, (0, 2, 1, 3, 4))


def dual(Xb):
    G = [[Fr(int(np.trace(Xb[a] @ Xb[b]))) for b in range(9)] for a in range(9)]
    Gi = X8.finv(G)
    den = 1
    for r in Gi:
        for v in r:
            den = den * v.denominator // np.gcd(den, v.denominator)
    return np.array([sum(int(Gi[a][b] * den) * Xb[b] for b in range(9)) for a in range(9)]), den


def kcoords(B, eta):
    return [Fr(B[i][j]) * eta[i] for (i, j) in PAIRS]      # K_(ij) has entries [i][j] = eta_i, [j][i] = -eta_j


def curv_op(R, eta):
    """matrix (K basis) of the curvature operator so(eta) -> so(eta), E_cd -> R(e_c, e_d); K_(ij) = eta_i eta_j E_ij"""
    Bc = X8.endos(R, eta)
    cols = [[eta[i] * eta[j] * v for v in kcoords(Bc[(i, j)], eta)] for (i, j) in PAIRS]
    return sp.Matrix(6, 6, lambda r, c: sp.Rational(cols[c][r]))


def BR_coeff(R, eta):
    """coefficient matrix (K basis) of B_R = sum_cd eps_c eps_d R_cd (x) R_cd"""
    Bc = X8.endos(R, eta)
    b = sp.zeros(6, 6)
    for c in range(4):
        for d in range(4):
            v = sp.Matrix([sp.Rational(e) for e in kcoords(Bc[(c, d)], eta)])
            b += eta[c] * eta[d] * v * v.T
    return b


def R_from_pairs(Rp):
    R = [[[[Fr(0)] * 4 for _ in range(4)] for _ in range(4)] for _ in range(4)]
    for I_, (a, b) in enumerate(PAIRS):
        for J_, (c, d) in enumerate(PAIRS):
            v = Fr(Rp[I_][J_])
            R[a][b][c][d], R[b][a][c][d], R[a][b][d][c], R[b][a][d][c] = v, -v, -v, v
    return R


def norm2(R, eta):
    return sum(eta[a] * eta[b] * eta[c] * eta[d] * R[a][b][c][d] ** 2 for a, b, c, d in itertools.product(range(4), repeat=4))


def ricci(R, eta):
    return [[sum(eta[a] * R[a][b][a][d] for a in range(4)) for d in range(4)] for b in range(4)]


def weyl(R, eta):
    Ric = ricci(R, eta)
    sc = sum(eta[b] * Ric[b][b] for b in range(4))
    g = lambda a, b: Fr(eta[a]) if a == b else Fr(0)  # noqa: E731
    W = [[[[R[a][b][c][d] - Fr(1, 2) * (Ric[a][c] * g(b, d) + Ric[b][d] * g(a, c) - Ric[a][d] * g(b, c) - Ric[b][c] * g(a, d))
            + sc / 6 * (g(a, c) * g(b, d) - g(a, d) * g(b, c)) for d in range(4)] for c in range(4)] for b in range(4)]
         for a in range(4)]
    return W


def bianchi_ok(R):
    return all(R[a][b][c][d] + R[a][c][d][b] + R[a][d][b][c] == 0 for a, b, c, d in itertools.product(range(4), repeat=4))


def R_L():
    """referee's example (base (1,3) or (3,1)): R_{01,12} = R_{02,23} = R_{03,31} = 1, all other pair components 0"""
    Rp = [[0] * 6 for _ in range(6)]
    for (I_, J_, v) in ((0, 3, 1), (1, 5, 1), (2, 4, -1)):          # (01,12), (02,23), (03,13) = -(03,31)
        Rp[I_][J_] = Rp[J_][I_] = v
    return R_from_pairs(Rp)


def qB(Xb, Kb, b):
    out = np.zeros((9, 9, 9, 4, 4), dtype=object)
    for a in range(6):
        for c in range(6):
            if b[a, c] != 0:
                out = out + qt(Xb, Kb[a], Kb[c]).astype(object) * Fr(int(sp.numer(b[a, c])), int(sp.denom(b[a, c])))
    return out


def z6():
    section("Z6  Lemma M5, the algebra of (C1)_0: own (C1) tensors; (C1)_0 <=> R^2 = -(c_H/c_V)^2/12 Id; Weyl tensors; "
            "an explicit Lorentzian solution")
    SIG = ((1, 1, 1, 1), (1, -1, -1, -1), (-1, 1, 1, 1), (1, 1, -1, -1))
    for eta in SIG:
        t1 = time.time()
        sigl = f"({eta.count(1)},{eta.count(-1)})"
        Xb, Kb = X8.pq_bases(eta)
        Xd, den = dual(Xb)
        C0, C1 = c0t(Xb), c1t(Xb)
        C1d = c1t(Xd)
        pair01 = Fr(int(np.einsum('abcij,abcji->', C0, C1d)), 4 * den ** 3)
        okq = True
        for K in Kb:
            for Kp in Kb:
                okq &= Fr(int(np.einsum('abcij,abcji->', qt(Xb, K, Kp), C1d)), 4 * den ** 3) == -9 * int(np.trace(K @ Kp))
        # q injective on Sym^2(so)
        vecs = [(qt(Xb, Kb[a], Kb[b]) + qt(Xb, Kb[b], Kb[a])).flatten() for a in range(6) for b in range(a, 6)]
        Mt = np.array(vecs)
        Mt = Mt[:, np.nonzero(np.any(Mt != 0, axis=0))[0]]
        rk = DomainMatrix([[sp.QQ(int(v)) for v in r] for r in Mt.tolist()], Mt.shape, sp.QQ).rank()
        kg = sp.Matrix(6, 6, lambda a, b: int(np.trace(Kb[a] @ Kb[b])))
        kgi = kg.inv()
        cas = sp.Matrix(kgi)
        okcas = all(v == 0 for v in (qB(Xb, Kb, cas) + 3 * C0.astype(object)).flatten())
        check(f"base {sigl}: own (C1) tensors: <c0, c1> = {pair01}, <q(K (x) K'), c1> = -9 tr(KK') on all 36 basis pairs; "
              f"q restricted to Sym^2(so(eta)) has rank {rk} = dim Sym^2(so(eta)) (injective); q(Cas) = -3 c0 {tick(t1)}",
              pair01 == 18 and okq and rk == 21 and okcas)
        t1 = time.time()
        rng = random.Random(20261134 + sum(eta))
        R = X8.rand_curv(rng, eta)
        Rm = curv_op(R, eta)
        b = BR_coeff(R, eta)
        ok = kgi * Rm.T * kg == Rm                                   # self-adjoint w.r.t. the trace form (pair symmetry)
        ok &= b == -4 * Rm * kgi * Rm.T and b == -4 * Rm * Rm * kgi
        n2 = norm2(R, eta)
        Bc = X8.endos(R, eta)
        sumtr = sum(eta[i] * eta[j] * X8.ftr(X8.fmul(Bc[(i, j)], Bc[(i, j)])) for i in range(4) for j in range(4))
        ok &= sumtr == -n2 and sp.Rational(n2) == 4 * (Rm * Rm).trace()
        pq = sum(b[a, c] * (-9) * int(np.trace(Kb[a] @ Kb[c])) for a in range(6) for c in range(6))
        ok &= sp.Rational(pq) == 9 * sp.Rational(n2)
        f_ = Fr(rng.randint(1, 9), rng.randint(1, 9))
        lhs = C0.astype(object) * Fr(1, 4) + qB(Xb, Kb, b) * (f_ / 4)
        rhs = qB(Xb, Kb, b * sp.Rational(f_.numerator, f_.denominator) - cas / 3) * Fr(1, 4)
        ok &= all(v == 0 for v in (lhs - rhs).flatten())
        check(f"base {sigl}, random curvature tensor R (|R|^2 = {n2}): the curvature operator is self-adjoint for the trace "
              f"form; B_R = sum eps_c eps_d R_cd (x) R_cd = -4 R kg^-1 R^T = -4 R^2 Cas; sum eps eps tr R_ij^2 = -|R|^2 = -4 tr "
              f"R^2; <q(B_R), c1> = 9 |R|^2 (so <(C1)_0, c1> = 18 + 9 f |R|^2, f = (c_V/c_H)^2: M5(b)); and c0 + f q(B_R) = "
              f"q(f B_R - Cas/3) (f = {f_}). With q injective: (C1)_0 <=> f B_R = Cas/3 <=> R^2 = -Id/(12 f) (referee) "
              f"{tick(t1)}", ok)
    # Weyl tensors, symbolically
    for eta in SIG:
        t1 = time.time()
        sigl = f"({eta.count(1)},{eta.count(-1)})"
        syms = {(I_, J_): sp.Symbol(f'r{I_}{J_}') for I_ in range(6) for J_ in range(I_, 6)}
        bian = {syms[(2, 3)]: syms[(1, 4)] - syms[(0, 5)]}        # R_{01,23} - R_{02,13} + R_{03,12} = 0
        Rp = [[syms[(min(I_, J_), max(I_, J_))].subs(bian) for J_ in range(6)] for I_ in range(6)]
        free = [s for k, s in syms.items() if k != (2, 3)]

        def Rs(a, b, c, d):
            if a == b or c == d:
                return sp.Integer(0)
            sg = 1
            if a > b:
                a, b, sg = b, a, -sg
            if c > d:
                c, d, sg = d, c, -sg
            return sg * Rp[PAIRS.index((a, b))][PAIRS.index((c, d))]
        okb = all(EX(Rs(a, b, c, d) + Rs(a, c, d, b) + Rs(a, d, b, c)) == 0 for a, b, c, d in itertools.product(range(4), repeat=4))
        g = lambda a, b: eta[a] if a == b else 0  # noqa: E731
        Ric = [[EX(sum(eta[a] * Rs(a, b, a, d) for a in range(4))) for d in range(4)] for b in range(4)]
        sc = EX(sum(eta[b] * Ric[b][b] for b in range(4)))

        def Ws(a, b, c, d):
            return EX(Rs(a, b, c, d) - sp.Rational(1, 2) * (Ric[a][c] * g(b, d) + Ric[b][d] * g(a, c) - Ric[a][d] * g(b, c)
                                                            - Ric[b][c] * g(a, d)) + sc / 6 * (g(a, c) * g(b, d) - g(a, d) * g(b, c)))
        okric = all(EX(sum(eta[a] * Ws(a, b, a, d) for a in range(4))) == 0 for b in range(4) for d in range(4))
        h = [eta[a] * eta[b] for (a, b) in PAIRS]
        RW = sp.Matrix(6, 6, lambda I_, J_: Ws(*PAIRS[I_], *PAIRS[J_]) * h[J_])
        St = sp.Matrix(6, 6, lambda I_, J_: Permutation(list(PAIRS[I_]) + list(PAIRS[J_])).signature()
                       if len(set(PAIRS[I_]) | set(PAIRS[J_])) == 4 else 0) * sp.diag(*h)
        s2 = St * St
        e2 = 1 if eta.count(-1) % 2 == 0 else -1
        okh = s2 == e2 * sp.eye(6)
        okh &= (RW * St - St * RW).applyfunc(EX) == sp.zeros(6, 6) and EX(RW.trace()) == 0 and EX((RW * St).trace()) == 0
        E = (RW * RW + sp.eye(6) / 12).applyfunc(EX)
        polys = [sp.Poly(e, *free) for e in E if e != 0]
        mons = sorted({m for p in polys for m in p.monoms()})
        rows = [[sp.QQ.from_sympy(dict(zip(p.monoms(), p.coeffs())).get(m, sp.Integer(0))) for m in mons] for p in polys]
        rk = DomainMatrix(rows, (len(rows), len(mons)), sp.QQ).rank()
        wdim = sp.Matrix([[sp.diff(Ws(a, b, c, d), s) for s in free] for a, b, c, d in itertools.product(range(4), repeat=4)]).rank()
        npos = sum(1 for v in h if v > 0)
        check(f"base {sigl}, generic Weyl tensor W (Weyl projection of the generic 20-parameter curvature tensor; {wdim} "
              f"independent components; Ric(W) = 0 {okric}): on Lambda^2 (metric signature ({npos},{6 - npos})), *^2 = {e2:+d}, "
              f"[R_W, *] = 0, tr R_W = 0, tr(R_W *) = 0 (Bianchi). Hence R_W is traceless on each 3-dim eigenspace of * (over C), "
              f"so R_W^2 = k Id, k != 0, is impossible (diagonalisable, eigenvalues +-sqrt(k), three of them cannot sum to 0): "
              f"no Ricci-flat solution of (C1)_0 even over C (M5(c), referee's proof). The 36 entries of R_W^2 + Id/12 have "
              f"rank {rk} (the note's 'rank 12') {tick(t1)}", okb and okric and okh and rk == 12 and wdim == 10)
    # real solutions: none in (2,2) (Lambda^2 of signature (2,4)), an explicit one in (1,3) and (3,1)
    t1 = time.time()
    RL = R_L()
    okL = bianchi_ok(RL)
    for eta in ((1, -1, -1, -1), (-1, 1, 1, 1)):
        Xb, Kb = X8.pq_bases(eta)
        Rm = curv_op(RL, eta)
        okL &= Rm * Rm == -sp.eye(6) and norm2(RL, eta) == -24 and any(v != 0 for r in ricci(RL, eta) for v in r)
        res = c0t(Xb).astype(object) + qB(Xb, Kb, BR_coeff(RL, eta)) * Fr(1, 12)
        okL &= all(v == 0 for v in res.flatten())
        res2 = c0t(Xb).astype(object) + qB(Xb, Kb, BR_coeff(RL, eta)) * Fr(1, 6)
        okL &= any(v != 0 for v in res2.flatten())
    check(f"observation (referee): R_L with R_(01,12) = R_(02,23) = R_(03,31) = 1 (and all symmetries; Bianchi holds) is a real "
          f"algebraic curvature tensor in bases (1,3) and (3,1) with curvature operator R_L^2 = -Id, Ric(R_L) != 0, |R_L|^2 = "
          f"-24 = -2 c_H^2/c_V^2 at (c_V/c_H)^2 = 1/12; (C1)_0 holds EXACTLY for (c_V/c_H)^2 = 1/12 > 0 (and fails at 1/6). "
          f"In base (2,2), Lambda^2 has signature (2,4), so no real self-adjoint R with R^2 = -k Id (k > 0) exists (R would "
          f"map a 4-dim negative subspace onto a positive one); in (4,0) none either. So real non-Ricci-flat solutions of (C1)_0 "
          f"exist exactly in the Lorentzian bases {tick(t1)}", okL)


# ============================================================================================ Z7 M5 brute force
def z7():
    section("Z7  Lemma M5 by brute force: xcheck_round8.fibre_system (order eps^2 at y = eta, Hess_F w = (1 + mu'(w - 1)) "
            "(Ric_F + s (A., A.)), c_H = c_V = 1, so f = (c_V/c_H)^2 = s)")
    for eta, seed in (((1, -1, -1, -1), 20261135), ((1, 1, -1, -1), 20261136)):
        t1 = time.time()
        sigl = f"({eta.count(1)},{eta.count(-1)})"
        rng = random.Random(seed)
        while True:
            W = weyl(X8.rand_curv(rng, eta), eta)
            n2 = norm2(W, eta)
            if n2 < 0:
                break
        f_ = Fr(-2) / n2
        okW = all(v == 0 for r in ricci(W, eta) for v in r) and bianchi_ok(W)
        rA, rAb, nu, nr = X8.fibre_system(eta, W, Fr(1), Fr(1), f_, Fr(0))
        # own linear algebra on (C1)_0 itself: c0 + f q(B_W) = 0 has no solution f (c0 not proportional to q(B_W))
        Xb, Kb = X8.pq_bases(eta)
        qw = qB(Xb, Kb, BR_coeff(W, eta)).flatten()
        c0v = c0t(Xb).astype(object).flatten()
        idx = next(i for i, v in enumerate(qw) if v != 0)
        lam_ = -c0v[idx] / qw[idx]
        nolin = any(c0v[i] + lam_ * qw[i] != 0 for i in range(len(qw)))
        check(f"base {sigl}: own Weyl tensor W (Weyl projection of a random curvature tensor; Ric = 0 {okW}), |W|^2 = {n2} < 0, "
              f"f = -2/|W|^2 = {f_} so that <(C1)_0, c1> = 0: c0 + f' q(B_W) = 0 has no solution f' (own linear algebra), and "
              f"the Taylor system with mu' = 0 ({nr} equations, {nu} unknowns) is inconsistent (ranks {rA}, {rAb}), as M5(c) "
              f"predicts {tick(t1)}", okW and nolin and rA != rAb)
    t1 = time.time()
    RL = R_L()
    runs = [((1, -1, -1, -1), Fr(1, 12), Fr(0), True), ((1, -1, -1, -1), Fr(1, 12), Fr(-7, 18), False),
            ((1, -1, -1, -1), Fr(1, 6), Fr(0), False), ((-1, 1, 1, 1), Fr(1, 12), Fr(0), True)]
    res = []
    ok = True
    for eta, f_, mu1, exp in runs:
        rA, rAb, nu, nr = X8.fibre_system(eta, RL, Fr(1), Fr(1), f_, mu1)
        res.append((f"({eta.count(1)},{eta.count(-1)})", str(f_), str(mu1), rA, rAb))
        ok &= (rA == rAb) == exp
    check(f"observation (referee): R_L (real, Lorentzian, not Ricci-flat): the order-eps^2 Taylor system is CONSISTENT for "
          f"f = 1/12, mu' = 0, and inconsistent for mu' = -7/18 or f = 1/6, i.e. exactly as (C1) predicts; (base, f, mu', rank, "
          f"augmented rank) = {res}. So in Lorentzian bases M5 alone does not exclude a fibre solution at s(X) with mu'(w(p)) "
          f"= 0 and physical constants; M6 step 6 (Ric_g = 0 on pi(O)) is needed there {tick(t1)}", ok)


# ============================================================================================ Z8 M9
def relaxed_system(eta, R, y0, cA, lam):
    """Own construction of the relaxed order-eps system of Hess_F w = mu(w) (Ric_F + (A., A.)) DIRECTLY at the fibre point y0
    (coordinates Y_A, y = y0 + sum Y_A E_A): rows = degree <= 1 Taylor coefficients of each (A, B) component, unknowns = the
    coefficients of w of degree 1..3 and a free covector xi replacing mu' dw(y0) in mu' (w - w0) Bt; mu(w(y0)) = 1.
    Gamma(h, k) = -(1/2)(h y^-1 k + k y^-1 h), Ric_F = -G_1/4, (A U, A W) = cA sum eps_i eps_j G_lam(U, R_ij.y) G_lam(W, R_ij.y),
    R.y = -(R^T y + y R), cA = c_V^2/(4 c_H^2).  Returns rows, rhs, the s-kernel vector test and the unknown count."""
    T = X8.TP(10, 1)
    y0i = X8.finv(y0)
    Y0i = X8.pm_const(T, y0i)
    Dl = [[T.var(X8.SYMIDX.index((min(i, j), max(i, j)))) for j in range(4)] for i in range(4)]
    yi = X8.pm_add(T, Y0i, X8.pm_mul(T, Y0i, X8.pm_mul(T, Dl, Y0i)), -1)
    Yv = X8.pm_add(T, X8.pm_const(T, y0), Dl)
    EA = [X8.pm_const(T, X8.Emat(A)) for A in range(10)]
    QA = [X8.pm_mul(T, yi, EA[A]) for A in range(10)]
    trQ = [X8.pm_tr(T, QA[A]) for A in range(10)]
    Gam = {}
    for A in range(10):
        for B in range(10):
            Mx = X8.pm_scl(T, X8.pm_add(T, X8.pm_mul(T, EA[A], QA[B]), X8.pm_mul(T, EA[B], QA[A])), Fr(-1, 2))
            Gam[(A, B)] = [Mx[i][j] for (i, j) in X8.SYMIDX]
    Bc = X8.endos(R, eta)
    fA = {}
    for ij, Bm in Bc.items():
        By = X8.pm_scl(T, X8.pm_add(T, X8.pm_mul(T, X8.pm_const(T, X8.fT(Bm)), Yv), X8.pm_mul(T, Yv, X8.pm_const(T, Bm))), -1)
        yBy = X8.pm_mul(T, yi, By)
        trB = X8.pm_tr(T, yBy)
        for A in range(10):
            fA[(ij, A)] = T.add(X8.pm_trprod(T, QA[A], yBy), T.mul(trQ[A], trB), -lam)
    Bt = {}
    for A in range(10):
        for B in range(A, 10):
            ric = T.scl(T.add(X8.pm_trprod(T, QA[A], QA[B]), T.mul(trQ[A], trQ[B]), Fr(-1, 4)), -1)
            at = T.ssum([T.scl(T.mul(fA[((i, j), A)], fA[((i, j), B)]), cA * eta[i] * eta[j]) for i in range(4) for j in range(4)])
            Bt[(A, B)] = T.add(ric, at)
    W3 = X8.TP(10, 3)
    unk = [k for k in range(len(W3.mons)) if W3.dg[k] >= 1]
    uix = {k: i for i, k in enumerate(unk)}
    nw = len(unk)
    rows, rhs = [], []
    for A in range(10):
        for B in range(A, 10):
            for kk, n_ in enumerate(T.mons):
                row = {}
                m = list(n_)
                m[A] += 1
                m[B] += 1
                j = uix[W3.ix[tuple(m)]]
                row[j] = row.get(j, 0) + Fr(m[A] * (m[A] - 1) if A == B else m[A] * m[B])
                for C in range(10):
                    for k1, gv in Gam[(A, B)][C].items():
                        n2 = [p - q for p, q in zip(n_, T.mons[k1])]
                        if min(n2) < 0:
                            continue
                        n2[C] += 1
                        j = uix[W3.ix[tuple(n2)]]
                        row[j] = row.get(j, 0) - gv * n2[C]
                if sum(n_) == 1:
                    b0 = Bt[(A, B)].get(0, Fr(0))
                    if b0:
                        C = n_.index(1)
                        row[nw + C] = row.get(nw + C, 0) - b0
                rows.append({j: v for j, v in row.items() if v})
                rhs.append(Bt[(A, B)].get(kk, Fr(0)))
    # s = (1/4) log det(y0^-1 y): its Taylor coefficients (degrees 1..3) with xi = 0 are a kernel vector
    Dl3 = [[W3.var(X8.SYMIDX.index((min(i, j), max(i, j)))) for j in range(4)] for i in range(4)]
    Mm = X8.pm_mul(W3, X8.pm_const(W3, y0i), Dl3)
    s, Pk = {}, Mm
    for k in range(1, 4):
        s = W3.add(s, W3.scl(X8.pm_tr(W3, Pk), Fr((-1) ** (k + 1), 4 * k)))
        Pk = X8.pm_mul(W3, Pk, Mm)
    kv = {uix[k]: v for k, v in s.items() if k in uix}
    okker = bool(kv) and all(sum(row.get(j, 0) * v for j, v in kv.items()) == 0 for row in rows)
    return rows, rhs, okker, nw + 10


def ranks(rows, rhs, nunk):
    dA, dAb = {}, {}
    for r, row in enumerate(rows):
        ent = {j: sp.QQ(v.numerator, v.denominator) for j, v in row.items()}
        if ent:
            dA[r] = dict(ent)
        if rhs[r]:
            ent[nunk] = sp.QQ(rhs[r].numerator, rhs[r].denominator)
        if ent:
            dAb[r] = ent
    return DomainMatrix(dA, (len(rows), nunk), sp.QQ).rank(), DomainMatrix(dAb, (len(rows), nunk + 1), sp.QQ).rank()


def z8():
    section("Z8  Proposition M9 by an own construction: the relaxed order-eps system expanded directly at y0 = P eta P^T "
            "(the note uses the pull-back y -> P y P^T to eta)")
    eta = (1, -1, -1, -1)
    etaM = [[Fr(eta[i]) if i == j else Fr(0) for j in range(4)] for i in range(4)]
    rng = random.Random(20261137)

    def randP():
        while True:
            Pm = [[Fr(rng.randint(-2, 2), 3) + (1 if i == j else 0) for j in range(4)] for i in range(4)]
            if sp.Matrix(Pm).det() != 0:
                return Pm
    t1 = time.time()
    R = X8.rand_curv(rng, eta)
    P1 = randP()
    y0 = X8.fmul(X8.fmul(P1, etaM), X8.fT(P1))
    rows, rhs, okk, nunk = relaxed_system(eta, R, y0, Fr(1, 4), Fr(1, 2))
    r1 = ranks(rows, rhs, nunk)
    rows2, rhs2, _, _ = relaxed_system(eta, R, y0, Fr(1, 4), Fr(2))
    okl = rows2 == rows and rhs2 == rhs
    check(f"base (1,3), random curvature tensor, random P (y0 = P eta P^T != eta), c_V/c_H = 1: relaxed system ({len(rows)} "
          f"equations, {nunk} unknowns) has (rank, augmented rank) = {r1}: inconsistent; the Taylor coefficients of s give a "
          f"kernel vector ({okk}); the system is identical for lambda = 1/2 and lambda = 2 (the lambda-term of G_lam(U, R.y) "
          f"carries tr(y^-1 R.y) = -2 tr R = 0, and Gamma, Ric_F are lambda-free), so M9 holds for every lambda {tick(t1)}",
          r1 == (294, 295) and okk and okl)
    t1 = time.time()
    rows, rhs, okk, nunk = relaxed_system(eta, R, etaM, Fr(1, 4), Fr(1, 2))
    r0 = ranks(rows, rhs, nunk)
    zeroR = [[[[Fr(0)] * 4 for _ in range(4)] for _ in range(4)] for _ in range(4)]
    rows, rhs, okk2, nunk = relaxed_system(eta, zeroR, y0, Fr(1, 4), Fr(1, 2))
    rz = ranks(rows, rhs, nunk)
    rows, rhs, okk3, nunk = relaxed_system(eta, R_L(), y0, Fr(1, 48), Fr(1, 2))
    rl = ranks(rows, rhs, nunk)
    check(f"base (1,3) controls: the same R at the section y0 = eta: {r0} (consistent, as L7 requires at order eps); R = 0 at "
          f"the off-section y0: {rz} (consistent: 0 lies in the exceptional set Z, as M9 says); the referee's R_L at "
          f"(c_V/c_H)^2 = 1/12, off-section: {rl} (inconsistent: R_L is not in Z although it passes M5 at s(X)) {tick(t1)}",
          r0 == (294, 294) and rz[0] == rz[1] and rl == (294, 295) and okk and okk2 and okk3)
    t1 = time.time()
    eta = (1, 1, -1, -1)
    etaM = [[Fr(eta[i]) if i == j else Fr(0) for j in range(4)] for i in range(4)]
    R = X8.rand_curv(rng, eta)
    P1 = randP()
    y0 = X8.fmul(X8.fmul(P1, etaM), X8.fT(P1))
    rows, rhs, okk, nunk = relaxed_system(eta, R, y0, Fr(1, 4), Fr(1, 2))
    r2 = ranks(rows, rhs, nunk)
    check(f"base (2,2), random curvature tensor, random P: relaxed system (rank, augmented rank) = {r2} (inconsistent), "
          f"s-kernel ({okk}) {tick(t1)}", r2 == (294, 295) and okk)


# ============================================================================================ main
def main():
    say("Referee cross-check of GU round 9 (M-z2-remaining-members.md), REVIEW7-GU.  Exact arithmetic; fixed seeds; one core.")
    tabs = z0()
    z1(tabs)
    out = z2(tabs)
    z3(out)
    z4(tabs)
    z5(tabs)
    z6()
    z7()
    z8()
    nf = [lab for lab, ok in RES if not ok]
    say(f"\nSUMMARY: {len(RES)} checks, {len(RES) - len(nf)} True, {len(nf)} False")
    say(f"runtime {time.time() - T0:.0f} s")


if __name__ == '__main__':
    main()
