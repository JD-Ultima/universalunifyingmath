"""GU round 8 (L): Upsilon = 0 on a neighbourhood of s(X) beyond the radial Z2 ansatz (non-radial t, u).

Formalization: formalization/geometric-unity/L-nonradial-z2.md.  Builds on J (J1, J2, Cor. J2.1, J5), J5d (L0, L1), F (F1,
F3, L-F1, L-F5) and K.  Imports J_slaved_sector.py (and through it I_elliptic_reduction.py, G_symbols.py, ...) and
F_jets.py, both unmodified.  Exact arithmetic throughout (Fractions, sympy rationals / Gaussian rationals); fixed seeds;
one core.  Every check is asserted.

Sections (labels = results of the note)
  L0  inputs: J1 (S, M) recomputed with J's own j0/j1; the J2 tables (l, b) transcribed and spot-checked against J's
      evaluator on sample basis tensors and pairs; then everything else uses the transcription.
  L1  Lambda^2 step (Lemma L1): exact varpi for Shiab#, W_A and the family x W_A + z W_C Ad_Gamma; leading varpi of every
      member; some varpi is a nonzero rational function, so dt ^ du = 0 on every open solution set.
  L2  Lemma L2 (Z2 has no open set with Ric in span(g, ds (x) ds) etc.): its inputs checked by brute force -- the Ricci
      tensor of Z2 at points OFF the section (F_jets, polynomial base metrics in general coordinates) against the O'Neill
      blocks Ric_F + (AU,AW) (VV) and Ric_g - 2(A_X,A_Y) (HH); Ric(E, .) = 0 and |E|^2 for the Euler field.
  L3  Lemma L3 (fibre lemma: a concircular function on a fibre R x SL(4)/SO(p,q) is affine in s with psi = 0): the
      algebraic core (the so(eta)-commutant in p is 0, [p,p] = so(eta)) and a Taylor-jet check of the Hessian equation.
  L4  Proposition L5 (curvature-free channel): y E1 - x E13 has no curvature term; its Hessian row h, the kernel field v, the
      residual R' = -1916928 (x^2 +- y^2)(xz +- yw)^2 (ty - ux) for EVERY member (symbolic); invariance of the line
      ty = ux under v iff xw = yz; the family x W_A + z W_C Ad_Gamma; W_A's E13; random exact test of q = mu^2 R'.
  L5  Theorem L8, step 5 (the residual line): for (z, w) parallel to (x, y) the odd equation on (t, u) = sigma (x, y) is x_i times a
      single equation, with mu' = -7/18 after the reparametrisation r = F(sigma) (symbolic in x, y, c, kappa).
  L6  Lemma L7 (fibre obstruction at a point of s(X)): the order-eps^2 integrability condition (C1); the pairing
      certificate <(C1), 4 c0 - 3 c1> = -360 mu' (identities on all of so (x) so); span ranks; brute-force Taylor
      solution of the fibre equation (inconsistent for mu' = -7/18, controls consistent where predicted).
Run (from computation/geometric-unity): python3 L_nonradial.py > L_nonradial.out 2>&1   (~6.5 min on one core)
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
import J_slaved_sector as J  # noqa: E402
import F_jets as FJ  # noqa: E402

T0 = time.time()
NCHK = [0]
KS = (0, 2, 12, 14)
x, y, z, w = J.xs, J.ys, J.zs, J.ws
t, u, kap = J.t_, J.u_, J.kap
I = sp.I
SG = {1: '(14,0)', -1: '(7,7)'}


def P(*a):
    print(*a)
    sys.stdout.flush()


def check(label, ok):
    NCHK[0] += 1
    P(f"  [{NCHK[0]:3d}] {label}: {ok}")
    assert ok, label


def section(s):
    P(f"\n=== {s}   [{time.time() - T0:.0f} s]")


def fr(v):
    """int / Fraction / sympy rational -> Fraction"""
    if isinstance(v, Fr):
        return v
    if isinstance(v, int):
        return Fr(v)
    v = sp.Rational(v)
    return Fr(int(v.p), int(v.q))


def QQ(v):
    if isinstance(v, Fr):
        return sp.QQ(v.numerator, v.denominator)
    return sp.QQ.from_sympy(sp.sympify(v))


# ====================================================================== tables (transcription of J1/J2)
def tables(eps):
    """S, M (J1) and the J2 tables l (derivative map) and b (quadratic map), eps = +1 (14,0), -1 (7,7)"""
    S = sp.Matrix([[0, 0], [-24 * x, -24 * eps * y], [24 * I * w, -24 * I * z], [0, 0]])
    sig, th, de = t * x + eps * u * y, t * w - u * z, t * y - u * x
    M = sp.Matrix([[kap, 0, 208 * eps * de, 0], [0, kap + 224 * sig, 384 * I * eps * th, 16 * eps * de],
                   [16 * de, -384 * I * th, kap + 544 * sig, 0], [0, 208 * de, 0, kap]])
    L = {('E1', 0): 0, ('E13', 0): 0, ('w3', 0): -x, ('w11', 0): -y,
         ('E1', 2): 24 * x, ('E13', 2): 24 * y, ('w3', 2): 11 * I * z, ('w11', 2): 11 * I * w,
         ('E1', 12): 24 * I * w * eps, ('E13', 12): -24 * I * z, ('w3', 12): 11 * y * eps, ('w11', 12): -11 * x,
         ('E1', 14): 0, ('E13', 14): 0, ('w3', 14): -I * w * eps, ('w11', 14): I * z}
    L = {k: sp.sympify(v) for k, v in L.items()}
    B = {}
    for k, l in itertools.combinations_with_replacement(KS, 2):
        B[('E1', k, l)] = (0, 0, 0)
        B[('E13', k, l)] = (0, 0, 0)
        B[('w3', k, l)] = 0
        B[('w11', k, l)] = 0
    r = sp.Rational
    B[('w3', 0, 2)], B[('w11', 0, 2)] = -r(22, 3) * x, -r(22, 3) * y
    B[('E1', 0, 12)], B[('E13', 0, 12)] = (0, -16 * y * eps, 16 * y * eps), (0, 16 * x, -16 * x)
    B[('w3', 0, 14)], B[('w11', 0, 14)] = r(2, 3) * y * eps, -r(2, 3) * x
    B[('E1', 2, 2)], B[('E13', 2, 2)] = (-224 * x, -224 * x, -176 * x), (-224 * y, -224 * y, -176 * y)
    B[('E1', 2, 12)] = (-32 * I * w * eps, -32 * I * w * eps, -352 * I * w * eps)
    B[('E13', 2, 12)] = (32 * I * z, 32 * I * z, 352 * I * z)
    B[('w3', 2, 12)], B[('w11', 2, 12)] = -r(266, 3) * y * eps, r(266, 3) * x
    B[('E1', 2, 14)], B[('E13', 2, 14)] = (16 * y * eps, 0, -16 * y * eps), (-16 * x, 0, 16 * x)
    B[('E1', 12, 12)] = (-160 * x * eps, -160 * x * eps, 528 * x * eps)
    B[('E13', 12, 12)] = (-160 * y * eps, -160 * y * eps, 528 * y * eps)
    B[('w3', 12, 14)], B[('w11', 12, 14)] = r(22, 3) * x * eps, r(22, 3) * y * eps
    B = {k: (tuple(sp.sympify(a) for a in v) if isinstance(v, tuple) else sp.sympify(v)) for k, v in B.items()}
    return S, M, L, B


def subs_tables(tab, sub):
    S, M, L, B = tab
    return (S.subs(sub), M.subs(sub), {k: v.subs(sub) for k, v in L.items()},
            {k: (tuple(a.subs(sub) for a in v) if isinstance(v, tuple) else v.subs(sub)) for k, v in B.items()})


def terms(Sig, dSig, L, B):
    """as J_slaved_sector.terms: Q (coefficient of dphi (x) dpsi), Q' (of <dphi, dpsi> g) in E1, E13, and the Lambda^2
    coefficients varpi (of dt ^ du) in E1_anti, E13_anti, w3, w11, for A_k = Sigma_k,t dt + Sigma_k,u du"""
    ph = ('t', 'u')
    Q, Qp, W = {}, {}, {}
    for nm in ('E1', 'E13'):
        for a, fa in enumerate(ph):
            for b_, fb in enumerate(ph):
                v = sum(L[(nm, k)] * dSig[fb][ik, a] for ik, k in enumerate(KS))
                vp = -v
                for (ik, k), (il, l) in itertools.combinations_with_replacement(list(enumerate(KS)), 2):
                    b1, b2, b3 = B[(nm, k, l)]
                    wkl = sp.Rational(1, 2) if k == l else 1
                    v += wkl * (b1 * Sig[ik, a] * Sig[il, b_] + b2 * Sig[il, a] * Sig[ik, b_])
                    vp += wkl * b3 * Sig[ik, a] * Sig[il, b_]
                Q[(nm, fa, fb)], Qp[(nm, fa, fb)] = v, vp
        W[nm] = (Q[(nm, 't', 'u')] - Q[(nm, 'u', 't')]) / 2
    for nm in ('w3', 'w11'):
        v = sum(L[(nm, k)] * (dSig['t'][ik, 1] - dSig['u'][ik, 0]) for ik, k in enumerate(KS))
        for (ik, k), (il, l) in itertools.combinations(list(enumerate(KS)), 2):
            v += B[(nm, k, l)] * (Sig[ik, 0] * Sig[il, 1] - Sig[ik, 1] * Sig[il, 0])
        W[nm] = v
    return Q, Qp, W


def reduced(tab):
    """exact Sigma, dSigma, H, Q, Q', varpi for a (numeric or partly symbolic) member"""
    S, M, L, B = tab
    detM = sp.factor(M.det())
    Sig = (-M.adjugate() * S / M.det()).applyfunc(sp.factor)
    dSig = {'t': Sig.diff(t).applyfunc(sp.factor), 'u': Sig.diff(u).applyfunc(sp.factor)}
    lE = sp.Matrix([[L[(c, k)] for k in KS] for c in ('E1', 'E13')])
    H = (lE * Sig).applyfunc(sp.factor)
    Q, Qp, W = terms(Sig, dSig, L, B)
    return dict(detM=detM, Sig=Sig, dSig=dSig, H=H, Q=Q, Qp=Qp, W=W)


# ====================================================================== L0
def l0():
    section("L0. Inputs: J1 recomputed with J's j0/j1; J2 tables transcribed and spot-checked against J's evaluator")
    for eta in (J.EUC, J.E77):
        eps = eta[7]
        C = J.Ctx(eta)
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            c1, c13 = J.j0(C)
            d1 = J.j1(C)
        S, M, L, B = tables(eps)
        check(f"{SG[eps]} J1: S and M of J's j1 equal the transcription",
              sp.expand(d1['S'] - S) == sp.zeros(4, 2) and sp.expand(d1['M'] - M) == sp.zeros(4, 4))
        ref1 = kap * t + x * (312 * t ** 2 + eps * 104 * u ** 2) + eps * 208 * y * t * u
        ref13 = kap * u + 208 * x * t * u + y * (104 * t ** 2 + eps * 312 * u ** 2)
        check(f"{SG[eps]} I6/J0: c1, c13 of J's j0", sp.expand(c1 - ref1) == 0 and sp.expand(c13 - ref13) == 0)
        Al = C.Al
        # derivative map: E = l (N^T - tr N g), w = l dA, read off at (p,q) = (0,1); pattern checked at (0,1), (2,2), (3,9)
        okd = True
        for k in KS:
            data = {}
            for (p, q) in ((0, 1), (2, 2), (3, 9)):
                Nm = [[0] * 14 for _ in range(14)]
                Nm[p][q] = 1
                po = J.per_op(Al.psi_form(C.dzeta(k, Nm)))
                for op in J.OPN:
                    (E1, E13, w3, w11), okx = C.extract(C.hat(po[op]))
                    okd &= okx
                    data[(op, p, q)] = (E1, E13, w3, w11)
            for ic, nm in enumerate(('E1', 'E13')):
                coef = {op: data[(op, 0, 1)][ic][1][0] for op in J.OPN}
                okd &= sp.expand(J.famsum(coef) - L[(nm, k)]) == 0
                for op in J.OPN:
                    for (p, q) in ((0, 1), (2, 2), (3, 9)):
                        Ea = data[(op, p, q)][ic]
                        al = coef[op]
                        okd &= all(Ea[a][b] == al * ((a == q and b == p) - (p == q) * eta[p] * (a == b) * eta[a])
                                   for a in range(14) for b in range(14))
            for ic, nm in ((2, 'w3'), (3, 'w11')):
                coef = {op: data[(op, 0, 1)][ic][0][1] for op in J.OPN}
                okd &= sp.expand(J.famsum(coef) - L[(nm, k)]) == 0
        check(f"{SG[eps]} J2 derivative table l: J's evaluator on basis tensors (0,1), (2,2), (3,9), k = 0, 2, 12, 14,"
              f" four operators, equals the transcription (and the pattern l (N^T - tr N g))", okd)
        okq = True
        for k, l in itertools.combinations_with_replacement(KS, 2):
            data = {}
            for (p, q) in ((0, 1), (0, 0)):
                X, Yv = C.P(k, C.e(p)), C.P(l, C.e(q))
                po = J.per_op(J.q2B(Al, X, Yv, {'x': Al.quad3(X), 'y': Al.quad3(Yv)}))
                for op in J.OPN:
                    (E1, E13, w3, w11), okx = C.extract(C.hat(po[op], 3))
                    okq &= okx
                    data[(op, p, q)] = (E1, E13, w3, w11)
            for ic, nm in enumerate(('E1', 'E13')):
                for j, getter in enumerate((lambda d: d[(0, 1)][ic][0][1], lambda d: d[(0, 1)][ic][1][0],
                                            lambda d: d[(0, 0)][ic][1][1] * eta[1] * eta[0])):
                    coef = {op: getter({pq: data[(op,) + pq] for pq in ((0, 1), (0, 0))}) for op in J.OPN}
                    okq &= sp.expand(J.famsum(coef) - B[(nm, k, l)][j]) == 0
            for ic, nm in ((2, 'w3'), (3, 'w11')):
                coef = {op: data[(op, 0, 1)][ic][0][1] for op in J.OPN}
                okq &= sp.expand(J.famsum(coef) - B[(nm, k, l)]) == 0
        check(f"{SG[eps]} J2 quadratic table b: J's evaluator on the pairs (e^0, e^1), (e^0, e^0) for all ten (k, l),"
              f" four operators, equals the transcription", okq)


# ====================================================================== L1
def l1():
    section("L1. Lemma L1 (Lambda^2 step): some varpi is a nonzero rational function of (t, u), so dt ^ du = 0 on every open"
            " solution set")
    out = {}
    for eps in (1, -1):
        tab = tables(eps)
        for nm, sub in (('Shiab#', {x: 1, y: 0, z: 1, w: 0}), ('W_A', {x: 1, y: 0, z: 0, w: 0}),
                        ('x W_A + z W_C AdG', {y: 0, w: 0})):
            r = reduced(subs_tables(tab, sub))
            keys = ('E13',) if nm.startswith('x W_A') else ('E1', 'E13', 'w3', 'w11')
            W = {k: sp.factor(sp.cancel(r['W'][k])) for k in keys}
            out[(eps, nm)] = (r, W)
            P(f"  {SG[eps]} {nm}: varpi({', '.join(keys)}) = {[W[k] for k in keys]}")
        dM = kap ** 4 + 768 * kap ** 3 * t + 121856 * kap ** 2 * t ** 2 - eps * 154112 * kap ** 2 * u ** 2 \
            - eps * 2555904 * kap * t * u ** 2 + 11075584 * u ** 4
        ref13 = -9216 * kap ** 2 * (-5 * kap ** 4 + 1600 * kap ** 3 * t + 609280 * kap ** 2 * t ** 2
                                    + eps * 930304 * kap ** 2 * u ** 2 + eps * 56016896 * kap * t * u ** 2
                                    - 587005952 * u ** 4) / dM ** 2
        check(f"{SG[eps]} Shiab#: exact varpi(E13) equals J's Cor. J2.1(e) formula, a nonzero rational function",
              sp.cancel(out[(eps, 'Shiab#')][1]['E13'] - ref13) == 0)
        refA = eps * 1916928 * kap * u / (kap ** 2 + 224 * kap * t - eps * 3328 * u ** 2) ** 2
        WA = out[(eps, 'W_A')][1]
        check(f"{SG[eps]} W_A: exact varpi(E1) = {'+' if eps == 1 else '-'}1916928 kappa u/(...)^2 (J), the other three 0",
              sp.cancel(WA['E1'] - refA) == 0 and all(WA[k] == 0 for k in ('E13', 'w3', 'w11')))
        Wf = out[(eps, 'x W_A + z W_C AdG')][1]
        check(f"{SG[eps]} family x W_A + z W_C AdG: varpi(E13) = x z^2 (nonzero polynomial)/(det M)^2, so nonzero for xz != 0",
              sp.cancel(Wf['E13'].subs({z: 1, x: 1}) - ref13) == 0 and sp.factor(Wf['E13']).has(z))
        # W_A-type members x W_A: varpi(E1) proportional to u
        WAx = reduced(subs_tables(tab, {y: 0, z: 0, w: 0}))['W']
        wa1 = sp.factor(sp.cancel(WAx['E1']))
        check(f"{SG[eps]} members x W_A: varpi(E1) = {wa1} (nonzero multiple of u), the other three 0",
              sp.cancel(wa1 / u).subs(u, 0) != 0 and all(sp.cancel(WAx[k]) == 0 for k in ('E13', 'w3', 'w11')))
        # every member: leading order (Neumann series, as J)
        S, M, L, B = tab
        M1, M13 = M.diff(t), M.diff(u)
        Ks = t * M1 + u * M13
        Sig_lo = -S / kap + Ks * S / kap ** 2
        _, _, Wl = terms(Sig_lo, {'t': M1 * S / kap ** 2, 'u': M13 * S / kap ** 2}, L, B)
        e_ = sp.Symbol('e_')
        lead = {k: sp.factor(sp.expand(sp.expand(v).subs(kap, 1 / e_)).coeff(e_, 2)) / kap ** 2 for k, v in Wl.items()}
        ell = x * z + eps * y * w
        okl = sp.expand(lead['w3'] * kap ** 2 - (-88320 * I * y * ell * (1 if eps == 1 else -1))) == 0 and \
            sp.expand(lead['w11'] * kap ** 2 - 88320 * I * x * ell) == 0
        check(f"{SG[eps]} every member, leading order: varpi(w3) = {'-' if eps == 1 else '+'}88320 i y (xz {'+' if eps == 1 else '-'} yw)/kappa^2,"
              f" varpi(w11) = +88320 i x (...)/kappa^2 (J, Cor. J2.1(e)): some varpi is nonzero whenever (x,y) != 0 and xz +- yw != 0",
              okl)
    return out


# ====================================================================== L2
def l2():
    section("L2. Lemma L2, inputs off the section: brute-force Ricci of Z2 (F_jets) at (x0, y0), y0 != g(x0), polynomial"
            " base metric in general coordinates, against the O'Neill blocks")
    ctx = FJ.Ctx('q')
    n = 4
    for eta, seed in (((1, 1, 1, 1), 11), ((1, -1, -1, -1), 12), ((1, 1, -1, -1), 13)):
        t1 = time.time()
        rng = random.Random(seed)

        def rsym(s=2):
            Mx = [[None] * n for _ in range(n)]
            for i in range(n):
                for j in range(i, n):
                    v = ctx.c(Fr(rng.randint(-s, s), rng.randint(1, 3)))
                    Mx[i][j] = v
                    Mx[j][i] = v
            return Mx
        g0 = [[ctx.c(eta[i]) if i == j else ctx.zero for j in range(n)] for i in range(n)]
        lin = [rsym() for _ in range(n)]
        quad = {(a, b): rsym() for a in range(n) for b in range(a, n)}
        cub = {(a, b, c): rsym() for a in range(n) for b in range(a, n) for c in range(b, n)}
        gj = FJ.base_from_coeffs(g0, lin, quad, cub, n, ctx)
        base = FJ.base_geometry(gj, n, ctx)
        lamF = Fr(1, 2) if -1 in eta else Fr(0)          # the side each real form requires (A4, A5)
        lam = ctx.c(lamF)
        cH, cV = ctx.c(Fr(2)), ctx.c(Fr(3))
        y0 = [[g0[i][j] + (ctx.c(Fr(rng.randint(-2, 2), 3)) if i <= j else 0) for j in range(n)] for i in range(n)]
        for i in range(n):
            for j in range(i):
                y0[i][j] = y0[j][i]
        G, m, idx, Psi0 = FJ.total_metric(gj, base, n, y0, lam, cH, cV, 'Z2', ctx)
        Ric, G0, G0i = FJ.ricci_at_point(G, m, ctx)
        RF = FJ.frame_blocks(Ric, n, idx, Psi0, ctx)
        N = len(idx)
        # exact rational copies for the predictions
        tofr = lambda v: Fr(int(sp.fraction(v.as_expr())[0]), int(sp.fraction(v.as_expr())[1]))  # noqa: E731
        Rup = {(c, a, i, j): eta[c] * tofr(base['R'][(c, a, i, j)]) for c, a, i, j in itertools.product(range(n), repeat=4)}
        Y0 = sp.Matrix(n, n, lambda i, j: sp.Rational(tofr(y0[i][j])))
        Y0i = Y0.inv()
        Ev = []
        for (i, j) in idx:
            Mx = sp.zeros(n, n)
            Mx[i, j] = 1
            Mx[j, i] = 1
            Ev.append(Mx)

        def Bdot(i, j):
            Bm = sp.Matrix(n, n, lambda c, a: Rup[(c, a, i, j)])     # B^c_a = R(d_i, d_j)^c_a
            return -(Bm.T * Y0 + Y0 * Bm)
        Kd = {(i, j): Bdot(i, j) for i in range(n) for j in range(n)}
        YE = [Y0i * E for E in Ev]
        YK = {k: Y0i * v for k, v in Kd.items()}
        cHr, cVr, lr = sp.Integer(2), sp.Integer(3), sp.Rational(lamF)
        Gl = lambda A_, B_, l_: (A_ * B_).trace() - l_ * A_.trace() * B_.trace()  # noqa: E731
        GK = {(A, k): (YE[A] * YK[k]).trace() for A in range(N) for k in Kd}
        okVV = True
        for A in range(N):
            for Bb in range(N):
                pred = -Gl(YE[A], YE[Bb], sp.Rational(1, 4))
                pred += cVr ** 2 / (4 * cHr ** 2) * sum(eta[i] * eta[j] * GK[(A, (i, j))] * GK[(Bb, (i, j))]
                                                        for i in range(n) for j in range(n))
                okVV &= sp.Rational(tofr(RF[n + A][n + Bb])) == pred
        okHH = True
        for a in range(n):
            for b in range(n):
                Phi = sum(eta[i] * (YK[(a, i)] * YK[(b, i)]).trace() for i in range(n))
                pred = sp.Rational(tofr(base['Ric'][a][b])) - cVr / (2 * cHr) * Phi
                okHH &= sp.Rational(tofr(RF[a][b])) == pred
        Evec = [tofr(y0[i][j]) for (i, j) in idx]
        okE = all(sum((tofr(RF[p][n + Bb]) * Evec[Bb] for Bb in range(N)), Fr(0)) == 0 for p in range(m))
        EE = cVr * Gl(Y0i * Y0, Y0i * Y0, lr)
        okHV = any(RF[a][n + Bb] != 0 for a in range(n) for Bb in range(N))
        sig = (eta.count(1), eta.count(-1))
        check(f"base {sig}, lambda = {lamF}, c_H = 2, c_V = 3, y0 off the section: VV block of Ric(Z2) = Ric_F + (AU, AW),"
              f" (AU,AW) = c_V^2/(4 c_H^2) sum_ij eps_i eps_j G(U, R_ij.y) G(W, R_ij.y)", okVV)
        check(f"base {sig}: HH block = Ric_g - 2(A_X, A_Y), (A_X, A_Y) = c_V/(4 c_H) Phi_y(X, Y) (J5d L0(a))", okHH)
        check(f"base {sig}: Euler field E = y d_y: Ric(E, .) = 0 and g(E, E) = 4 c_V (1 - 4 lambda) = {EE}; HV block"
              f" nonzero off the section (as expected): {okHV}  [{time.time() - t1:.0f} s]",
              okE and EE == 4 * cVr * (1 - 4 * lr) and okHV)


# ====================================================================== fibre jets (truncated polynomials)
class Poly:
    """truncated polynomials in NV variables: dict exponent-tuple -> Fraction"""

    def __init__(self, NV):
        self.NV = NV
        self.z0 = tuple([0] * NV)

    def add(self, a, b, s=1):
        o = dict(a)
        for k, v in b.items():
            o[k] = o.get(k, 0) + s * v
            if o[k] == 0:
                del o[k]
        return o

    def mul(self, a, b, D):
        o = {}
        for ka, va in a.items():
            da = sum(ka)
            for kb, vb in b.items():
                if da + sum(kb) > D:
                    continue
                k = tuple(p + q for p, q in zip(ka, kb))
                o[k] = o.get(k, 0) + va * vb
        return {k: v for k, v in o.items() if v != 0}

    def scale(self, a, s):
        return {k: v * s for k, v in a.items()} if s != 0 else {}

    def deriv(self, a, i):
        o = {}
        for k, v in a.items():
            if k[i] > 0:
                kk = list(k)
                kk[i] -= 1
                o[tuple(kk)] = o.get(tuple(kk), 0) + v * k[i]
        return o

    def trunc(self, a, D):
        return {k: v for k, v in a.items() if sum(k) <= D}

    def const(self, c):
        return {self.z0: Fr(c)} if c != 0 else {}

    def var(self, A):
        e = [0] * self.NV
        e[A] = 1
        return {tuple(e): Fr(1)}

    def ssum(self, lst):
        o = {}
        for p in lst:
            o = self.add(o, p)
        return o

    def monos(self, deg):
        out = []
        for comb in itertools.combinations_with_replacement(range(self.NV), deg):
            e = [0] * self.NV
            for c in comb:
                e[c] += 1
            out.append(tuple(e))
        return out


def fibre_data(eta, Rlow, lam=Fr(1, 2), maxorder=2, afac=Fr(1)):
    """the fibre F = Met(R^4) at y = eta + delta (10 coordinates delta_A), metric G_lambda (c_V = 1), to order maxorder:
    Christoffel symbols, Ric_F = -G_{1/4} (F1) and the O'Neill term (AU, AW) of the base point with R = R_g(x0)
    (c_V/c_H absorbed into R; afac multiplies the A-term, for controls)"""
    n = 4
    idx = [(i, j) for i in range(n) for j in range(i, n)]
    NV = len(idx)
    Pp = Poly(NV)
    Do = maxorder + 1
    d = [[{} for _ in range(n)] for _ in range(n)]
    for A, (i, j) in enumerate(idx):
        d[i][j] = Pp.add(d[i][j], Pp.var(A))
        if i != j:
            d[j][i] = Pp.add(d[j][i], Pp.var(A))
    Y = [[Pp.add(d[i][j], Pp.const(eta[i]) if i == j else {}) for j in range(n)] for i in range(n)]

    def mm(A_, B_, D):
        return [[Pp.trunc(Pp.ssum([Pp.mul(A_[i][k], B_[k][j], D) for k in range(n)]), D) for j in range(n)]
                for i in range(n)]
    Ed = [[Pp.scale(d[i][j], Fr(eta[i])) for j in range(n)] for i in range(n)]
    Yi = [[Pp.const(eta[i]) if i == j else {} for j in range(n)] for i in range(n)]
    term = [[dict(Yi[i][j]) for j in range(n)] for i in range(n)]
    for _ in range(Do):
        term = mm(Ed, term, Do)
        term = [[Pp.scale(term[i][j], -1) for j in range(n)] for i in range(n)]
        Yi = [[Pp.add(Yi[i][j], term[i][j]) for j in range(n)] for i in range(n)]
    Ev = []
    for (i, j) in idx:
        Mx = [[{} for _ in range(n)] for _ in range(n)]
        Mx[i][j] = Pp.const(1)
        Mx[j][i] = Pp.const(1)
        Ev.append(Mx)
    tr = lambda Mx: Pp.ssum([Mx[i][i] for i in range(n)])  # noqa: E731
    YiE = [mm(Yi, Ev[A], Do) for A in range(NV)]
    trYiE = [tr(YiE[A]) for A in range(NV)]

    def G(a, b, l_, D):
        return Pp.add(tr(mm(YiE[a], YiE[b], D)), Pp.mul(trYiE[a], trYiE[b], D), -l_)
    Gm = [[G(a, b, lam, Do) for b in range(NV)] for a in range(NV)]
    G0 = sp.Matrix(NV, NV, lambda a, b: sp.Rational(Gm[a][b].get(Pp.z0, 0)))
    G0i = G0.inv()
    G0iF = [[Fr(int(sp.fraction(G0i[a, b])[0]), int(sp.fraction(G0i[a, b])[1])) for b in range(NV)] for a in range(NV)]
    dGm = [[Pp.add(Gm[a][b], Pp.const(Gm[a][b].get(Pp.z0, 0)), -1) for b in range(NV)] for a in range(NV)]

    def mN(A_, B_, D):
        return [[Pp.ssum([Pp.mul(A_[i][k], B_[k][j], D) for k in range(NV) if A_[i][k] and B_[k][j]])
                 for j in range(NV)] for i in range(NV)]
    G0iP = [[Pp.const(G0iF[a][b]) for b in range(NV)] for a in range(NV)]
    Gi = [[dict(G0iP[a][b]) for b in range(NV)] for a in range(NV)]
    term = [[dict(G0iP[a][b]) for b in range(NV)] for a in range(NV)]
    Xm = mN(G0iP, dGm, maxorder)
    for _ in range(maxorder):
        term = mN(Xm, term, maxorder)
        term = [[Pp.scale(term[i][j], -1) for j in range(NV)] for i in range(NV)]
        Gi = [[Pp.add(Gi[i][j], term[i][j]) for j in range(NV)] for i in range(NV)]
    dG = [[[Pp.deriv(Gm[a][b], c) for b in range(NV)] for a in range(NV)] for c in range(NV)]
    Gam1 = [[[Pp.scale(Pp.add(Pp.add(dG[a][c][b], dG[b][c][a]), dG[c][a][b], -1), Fr(1, 2)) for b in range(NV)]
             for a in range(NV)] for c in range(NV)]
    Gam = [[[Pp.trunc(Pp.ssum([Pp.mul(Gi[e][c], Gam1[c][a][b], maxorder) for c in range(NV)]), maxorder)
             for b in range(NV)] for a in range(NV)] for e in range(NV)]
    RicF = [[Pp.scale(G(a, b, Fr(1, 4), maxorder), -1) for b in range(NV)] for a in range(NV)]
    Bm = {(c, dd): [[Fr(eta[a]) * fr(Rlow(a, b, c, dd)) for b in range(n)] for a in range(n)]
          for c in range(n) for dd in range(n)}

    def Bdot(Bx):
        return [[Pp.scale(Pp.ssum([Pp.add(Pp.scale(Y[c][b], Bx[c][a]), Pp.scale(Y[a][c], Bx[c][b])) for c in range(n)]), -1)
                 for b in range(n)] for a in range(n)]
    K = {k: Bdot(v) for k, v in Bm.items()}
    GK = {(A, k): tr(mm(YiE[A], mm(Yi, K[k], maxorder), maxorder)) for A in range(NV) for k in K}
    Aterm = [[Pp.scale(Pp.ssum([Pp.scale(Pp.mul(GK[(a, (i, j))], GK[(b, (i, j))], maxorder), Fr(eta[i] * eta[j]))
                                for i in range(n) for j in range(n)]), Fr(1, 4) * afac)
              for b in range(NV)] for a in range(NV)]
    Bt = [[Pp.add(RicF[a][b], Aterm[a][b]) for b in range(NV)] for a in range(NV)]
    return dict(Pp=Pp, NV=NV, Gam=Gam, Gm=Gm, Bt=Bt, idx=idx)


def hessian_system(fd, maxorder, D, extra_unknowns=None, mu0=None, mu1=None):
    """linear system for the Taylor coefficients (degree 1..D) of w at y = eta in
         Hess_F w - (mu0 + mu1 (w - w0)) Bt = 0      (Bt = Ric_F + (A.,A.); w0 = 1)
       to order maxorder; or, with extra_unknowns = 'psi', Hess_F w - psi G = 0 (psi a polynomial of degree <= maxorder)."""
    Pp, NV, Gam, Bt, Gm = fd['Pp'], fd['NV'], fd['Gam'], fd['Bt'], fd['Gm']
    unk = [m for dd in range(1, D + 1) for m in Pp.monos(dd)]
    rowkeys = [(a, b, beta) for a in range(NV) for b in range(a, NV) for dd in range(maxorder + 1) for beta in Pp.monos(dd)]
    rk = {k: i for i, k in enumerate(rowkeys)}
    cols = []
    for m in unk:
        phi = {m: Fr(1)}
        col = {}
        for a in range(NV):
            for b in range(a, NV):
                e = Pp.deriv(Pp.deriv(phi, a), b)
                for c in range(NV):
                    if Gam[c][a][b]:
                        e = Pp.add(e, Pp.mul(Gam[c][a][b], Pp.deriv(phi, c), maxorder), -1)
                if mu1 is not None and mu1 != 0:
                    e = Pp.add(e, Pp.scale(Pp.mul(phi, Bt[a][b], maxorder), mu1), -1)
                for k, v in Pp.trunc(e, maxorder).items():
                    col[rk[(a, b, k)]] = v
        cols.append(col)
    psi_cols = []
    if extra_unknowns == 'psi':
        for dd in range(maxorder + 1):
            for m in Pp.monos(dd):
                col = {}
                for a in range(NV):
                    for b in range(a, NV):
                        e = Pp.scale(Pp.mul({m: Fr(1)}, Pp.trunc(Gm[a][b], maxorder), maxorder), -1)
                        for k, v in e.items():
                            col[rk[(a, b, k)]] = v
                psi_cols.append((m, col))
    rhs = {}
    if mu0 is not None:
        for a in range(NV):
            for b in range(a, NV):
                for k, v in Pp.trunc(Pp.scale(Bt[a][b], mu0), maxorder).items():
                    rhs[rk[(a, b, k)]] = v
    return unk, rowkeys, cols, psi_cols, rhs


def dm_from_cols(nrows, cols, rhs=None):
    used = sorted({i for c in cols for i in c} | (set(rhs) if rhs else set()))
    pos = {i: r for r, i in enumerate(used)}
    C = len(cols) + (1 if rhs is not None else 0)
    rows = [[sp.QQ(0)] * C for _ in used]
    for j, col in enumerate(cols):
        for i, v in col.items():
            rows[pos[i]][j] = QQ(v)
    if rhs is not None:
        for i, v in rhs.items():
            rows[pos[i]][C - 1] = QQ(v)
    return DomainMatrix(rows, (len(used), C), sp.QQ)


# ====================================================================== L3
def pbasis(eta):
    n = 4
    E = sp.diag(*eta)
    out = []
    for a in range(n):
        for b in range(a, n):
            h = sp.zeros(n, n)
            h[a, b] = 1
            h[b, a] = 1
            H = E.inv() * h
            out.append(H - H.trace() / n * sp.eye(n))
    sel, Mx = [], []
    for H in out:
        if sp.Matrix(Mx + [list(H)]).rank() > len(Mx):
            Mx.append(list(H))
            sel.append(H)
    return sel


def kbasis(eta):
    n = 4
    E = sp.diag(*eta)
    out = []
    for a in range(n):
        for b in range(a + 1, n):
            k = sp.zeros(n, n)
            k[a, b] = 1
            k[b, a] = -1
            out.append(E.inv() * k)       # eta-skew: K^T eta = -eta K
    return out


def br(A, C):
    return A * C - C * A


def l3():
    section("L3. Lemma L3 (fibre lemma): concircular functions on the fibre R x SL(4)/SO(p,q) are affine in s, psi = 0")
    for eta in ((1, 1, 1, 1), (1, -1, -1, -1), (1, 1, -1, -1)):
        Pb, Kb = pbasis(eta), kbasis(eta)
        E = sp.diag(*eta)
        okk = all((Kx.T * E + E * Kx) == sp.zeros(4, 4) for Kx in Kb)
        brs = [list(br(X, Y)) for X, Y in itertools.combinations(Pb, 2)]
        rk = sp.Matrix(brs).rank()
        vs = sp.symbols('v0:9')
        V = sum((vs[i] * Pb[i] for i in range(9)), sp.zeros(4, 4))
        eqs = [e for X, Y in itertools.combinations(Pb, 2) for e in br(br(X, Y), V)]
        sol = sp.solve(eqs, vs, dict=True)
        okc = len(sol) == 1 and all(v == 0 for v in sol[0].values()) and len(sol[0]) == 9
        check(f"base {(eta.count(1), eta.count(-1))}: [p, p] spans so(eta) (rank {rk} = 6; so(eta) checked), and"
              f" {{v in p : [[X, Y], v] = 0 for all X, Y in p}} = 0, i.e. R^S(X, Y) v = 0 for all X, Y forces v = 0",
              okk and rk == 6 and okc)
    # Taylor-jet check of the concircular equation Hess_F w = psi g_F at y = eta (unknown w of degree <= 4, psi <= 2)
    for eta in ((1, 1, 1, 1), (1, -1, -1, -1)):
        t1 = time.time()
        fd = fibre_data(eta, lambda *a: 0, maxorder=2)
        unk, rowkeys, cols, psi_cols, _ = hessian_system(fd, 2, 4, extra_unknowns='psi')
        allc = cols + [c for _, c in psi_cols]
        Dm = dm_from_cols(len(rowkeys), allc)
        Nsp = Dm.nullspace().to_Matrix()
        nw = len(cols)
        psi_idx = {m: nw + i for i, (m, _) in enumerate(psi_cols)}
        Pp = fd['Pp']
        lowpsi = [psi_idx[m] for dd in (0, 1) for m in Pp.monos(dd)]
        ok_psi = all(Nsp[r, c] == 0 for r in range(Nsp.shape[0]) for c in lowpsi)
        # degree-1 part of w must be proportional to ds = (1/4) tr(eta^-1 delta)
        idx = fd['idx']
        ds = [Fr(eta[i], 4) if i == j else Fr(0) for (i, j) in idx]
        lin_idx = [unk.index(m) for m in Pp.monos(1)]
        okw = True
        dsr = [sp.Rational(v.numerator, v.denominator) for v in ds]
        for r in range(Nsp.shape[0]):
            vec = [sp.Rational(Nsp[r, c]) for c in lin_idx]
            okw &= all(vec[a] * dsr[b] - vec[b] * dsr[a] == 0 for a in range(10) for b in range(10))
        check(f"base {(eta.count(1), eta.count(-1))}: jets at y = eta of solutions of Hess_F w = psi G_lambda (to order 2,"
              f" {len(allc)} unknowns): psi(eta) = 0, d psi(eta) = 0 and dw(eta) is a multiple of ds"
              f" (solution space dim {Nsp.shape[0]})  [{time.time() - t1:.0f} s]", ok_psi and okw)


# ====================================================================== L4
def l4():
    section("L4. Proposition L5: the curvature-free channel y E1 - x E13; residual R'; invariant line; the family; W_A")
    s_ = sp.Symbol('sigma')
    for eps in (1, -1):
        t1 = time.time()
        S, M, L, B = tables(eps)
        Mk = M.subs(kap, 1)
        D = sp.expand(Mk.det())
        SN = (-Mk.adjugate() * S).applyfunc(sp.expand)
        lE = sp.Matrix([[L[(c, k)] for k in KS] for c in ('E1', 'E13')])
        HN = (lE * SN).applyfunc(sp.expand)
        hN = (y * HN[0, :] - x * HN[1, :]).applyfunc(sp.expand)
        ell = x * z + eps * y * w
        g0 = 576 * ell
        vt, vu = sp.cancel(hN[1] / g0), sp.cancel(-hN[0] / g0)
        okpoly = sp.denom(vt) == 1 and sp.denom(vu) == 1
        vv = [vt, vu]
        vgrad = lambda f: vv[0] * sp.diff(f, t) + vv[1] * sp.diff(f, u)  # noqa: E731
        lc = {k: y * L[('E1', k)] - x * L[('E13', k)] for k in KS}
        SvN = [sp.expand(SN[ik, 0] * vt + SN[ik, 1] * vu) for ik in range(4)]
        term1 = sp.expand(hN[0] * sp.expand(vgrad(vt)) + hN[1] * sp.expand(vgrad(vu)))
        vgD = sp.expand(vgrad(D))
        term2 = 0
        for ik, k in enumerate(KS):
            if lc[k] != 0:
                a = sp.expand(vt * vgrad(SN[ik, 0]) + vu * vgrad(SN[ik, 1]))
                term2 += lc[k] * (a * D - SvN[ik] * vgD)
        term3 = 0
        for (ik, k), (il, l) in itertools.combinations_with_replacement(list(enumerate(KS)), 2):
            b1 = [y * p - x * q for p, q in zip(B[('E1', k, l)], B[('E13', k, l)])]
            wkl = sp.Rational(1, 2) if k == l else 1
            if b1[0] + b1[1] != 0:
                term3 += wkl * (b1[0] + b1[1]) * SvN[ik] * SvN[il]
        total = sp.expand(term1 * D + sp.expand(term2) + sp.expand(term3))
        qq, rr = sp.div(sp.Poly(total, t, u, x, y, z, w), sp.Poly(sp.expand(D ** 2), t, u, x, y, z, w))
        N_ = x ** 2 + eps * y ** 2
        de = t * y - u * x
        okR = rr.is_zero and sp.expand(qq.as_expr() + 1916928 * N_ * ell ** 2 * de) == 0
        check(f"{SG[eps]} every member (kappa = 1): curvature-free row h = (y, -x) H, kernel field v = (h_u, -h_t) det M /"
              f" (576 (xz {'+' if eps == 1 else '-'} yw)) is polynomial ({okpoly}); the residual R' := h.(Dv)v + v^T Q'_sym v"
              f" = -1916928 (x^2 {'+' if eps == 1 else '-'} y^2)(xz {'+' if eps == 1 else '-'} yw)^2 (ty - ux) exactly"
              f"  [{time.time() - t1:.0f} s]", okpoly and okR)
        inv = sp.factor(sp.expand((y * vt - x * vu).subs({t: s_ * x, u: s_ * y})))
        ref = -(w * x - y * z) * (160 * s_ * N_ - 1)
        check(f"{SG[eps]} on the line (t, u) = sigma (x, y): v.grad(ty - ux) = -(xw - yz)(160 sigma (x^2"
              f" {'+' if eps == 1 else '-'} y^2) - 1): the line is v-invariant iff xw = yz", sp.expand(inv - ref) == 0)
        # the family x W_A + z W_C AdG  (y = w = 0), kappa symbolic
        r = reduced(subs_tables(tables(eps), {y: 0, w: 0}))
        H = r['H']
        h13 = H[1, :]
        ref_h = [-221184 * kap ** 2 * u * x * z ** 2, 576 * kap * z ** 2 * (kap ** 2 + 224 * kap * t * x - eps * 3328 * u ** 2 * x ** 2)]
        okh = all(sp.cancel(h13[i] * r['detM'] - ref_h[i]) == 0 for i in range(2))
        vf = [kap ** 2 + 224 * kap * t * x - eps * 3328 * u ** 2 * x ** 2, 384 * kap * u * x]
        okvf = sp.cancel(h13[0] * vf[0] + h13[1] * vf[1]) == 0
        Dv = sp.Matrix([[sp.diff(vf[i], t), sp.diff(vf[i], u)] for i in range(2)])
        Qs = sp.Matrix(2, 2, lambda i, j: (r['Q'][('E13', 'tu'[i], 'tu'[j])] + r['Q'][('E13', 'tu'[j], 'tu'[i])]) / 2)
        vfm = sp.Matrix(vf)
        Rf = sp.factor(sp.cancel((h13 * (Dv * vfm))[0] + (vfm.T * Qs * vfm)[0]))
        check(f"{SG[eps]} family x W_A + z W_C AdG (kappa symbolic): H13 = z^2 (..)/det M, kernel field v = (kappa^2 + 224 kappa"
              f" t x {'-' if eps == 1 else '+'} 3328 u^2 x^2, 384 kappa u x), R = {Rf}; v is tangent to u = 0",
              okh and okvf and sp.expand(Rf + 1916928 * kap * u * x ** 4) == 0 and vf[1].subs(u, 0) == 0)
        sol = sp.solve([ref_h[0], ref_h[1]], [t, u], dict=True)
        okz = all(sp.factor(r['detM'].subs(s)) == 0 for s in sol if s)
        check(f"{SG[eps]} family, z != 0: H13 = 0 only at u = 0, kappa + 224 t x = 0, where det M = 0", okz)
        # W_A-type members x W_A: H13 = 0, Sigma has no u-column, E13 = c13 g - 16 x Sigma_C Sigma_B (dt dt - |dt|^2 g)
        ra = reduced(subs_tables(tables(eps), {y: 0, z: 0, w: 0}))
        Sg = ra['Sig']
        okA = ra['H'][1, 0] == 0 and ra['H'][1, 1] == 0 and all(Sg[i, 1] == 0 for i in range(4))
        Q13 = sp.cancel(ra['Q'][('E13', 't', 't')] + 16 * x * Sg[1, 0] * Sg[3, 0])
        Qp13 = sp.cancel(ra['Qp'][('E13', 't', 't')] - 16 * x * Sg[1, 0] * Sg[3, 0])
        okA &= Q13 == 0 and Qp13 == 0 and all(sp.cancel(ra['Q'][('E13', a, b)]) == 0 for a, b in (('t', 'u'), ('u', 't'), ('u', 'u')))
        SB = sp.factor(Sg[3, 0])
        SC = sp.factor(Sg[1, 0])
        okA &= sp.cancel(SB / u).subs(u, 0) != 0 and SC.subs(u, 0) != 0
        check(f"{SG[eps]} members x W_A: H13 = 0, Sigma = (0, Sigma_C, 0, Sigma_B) dt with Sigma_C = {SC},"
              f" Sigma_B = {SB}; E13 = c13 g - 16 x Sigma_C Sigma_B (dt (x) dt - |dt|^2 g)", okA)
        # random exact test of the reduction q = mu^2 R along an integral curve (gamma' = mu v, gamma'' = mu' v + mu^2 (Dv) v)
        rng = random.Random(50 + eps)
        mem = {x: 2, y: -1, z: 3, w: 1, kap: 1}
        rm = reduced(subs_tables(tables(eps), mem))
        hm = (mem[y] * rm['H'][0, :] - mem[x] * rm['H'][1, :])
        Qm = sp.Matrix(2, 2, lambda i, j: mem[y] * (rm['Q'][('E1', 'tu'[i], 'tu'[j])]) - mem[x] * (rm['Q'][('E13', 'tu'[i], 'tu'[j])]))
        vtm, vum = vt.subs(mem), vu.subs(mem)
        Rm = -1916928 * (mem[x] ** 2 + eps * mem[y] ** 2) * (mem[x] * mem[z] + eps * mem[y] * mem[w]) ** 2 * (t * mem[y] - u * mem[x])
        okq = True
        for _ in range(3):
            pt = {t: sp.Rational(rng.randint(-9, 9), 1000), u: sp.Rational(rng.randint(-9, 9), 1000)}
            mu, mup = sp.Rational(rng.randint(1, 9), 7), sp.Rational(rng.randint(-9, 9), 5)
            v0 = sp.Matrix([vtm.subs(pt), vum.subs(pt)])
            Dv0 = sp.Matrix([[sp.diff(vtm, t).subs(pt), sp.diff(vtm, u).subs(pt)], [sp.diff(vum, t).subs(pt), sp.diff(vum, u).subs(pt)]])
            g1 = mu * v0
            g2 = mup * v0 + mu ** 2 * Dv0 * v0
            hq = (hm.subs(pt) * g1)[0]
            qhat = (hm.subs(pt) * g2)[0] + (g1.T * Qm.subs(pt) * g1)[0]
            okq &= sp.simplify(hq) == 0 and sp.simplify(qhat - mu ** 2 * Rm.subs(pt)) == 0
        check(f"{SG[eps]} member (2, -1, 3, 1), three random exact points and curve jets gamma' = mu v, gamma'' = mu' v +"
              f" mu^2 (Dv) v: h.gamma' = 0 and the traceless dr (x) dr coefficient of y E1 - x E13 equals mu^2 R'", okq)


# ====================================================================== L5
def l5():
    section("L5. Theorem L8, step 5: the residual line (t, u) = sigma (x, y) for (z, w) = c (x, y): one equation, mu' = -7/18")
    s_, c_ = sp.symbols('sigma c')
    for eps in (1, -1):
        t1 = time.time()
        tab = subs_tables(tables(eps), {z: c_ * x, w: c_ * y})
        S, M, L, B = tab
        Dm = sp.expand(M.det())
        SN = (-M.adjugate() * S).applyfunc(sp.expand)
        Sig = SN / Dm
        dSig = {'t': (SN.diff(t) * Dm - SN * sp.diff(Dm, t)) / Dm ** 2, 'u': (SN.diff(u) * Dm - SN * sp.diff(Dm, u)) / Dm ** 2}
        line = {t: s_ * x, u: s_ * y}
        Sl = Sig.subs(line).applyfunc(sp.factor)
        dSl = {k: v.subs(line).applyfunc(sp.factor) for k, v in dSig.items()}
        lE = sp.Matrix([[L[(cc, k)] for k in KS] for cc in ('E1', 'E13')])
        Hl = (lE * Sl).applyfunc(sp.factor)
        Q, Qp, _ = terms(Sl, dSl, L, B)
        vec = (x, y)
        N_ = x ** 2 + eps * y ** 2
        out = {}
        for ic, (ch, xi) in enumerate((('E1', x), ('E13', y))):
            Ht = sp.factor(Hl[ic, 0] * x + Hl[ic, 1] * y)
            Qt = sp.factor(sum(Q[(ch, a, b)] * vec[i] * vec[j] for i, a in enumerate('tu') for j, b in enumerate('tu')))
            Qpt = sp.factor(sum(Qp[(ch, a, b)] * vec[i] * vec[j] for i, a in enumerate('tu') for j, b in enumerate('tu')))
            out[ch] = (Ht, Qt, Qpt)
            ok = sp.cancel(Ht - 576 * xi * kap * N_ / (kap ** 2 + 224 * kap * s_ * N_)) == 0 and \
                sp.cancel(Qt + 258048 * xi * N_ ** 2 / (kap + 224 * s_ * N_) ** 2) == 0 and \
                sp.cancel(Qpt - 78336 * xi * N_ ** 2 / (kap + 224 * s_ * N_) ** 2) == 0
            check(f"{SG[eps]} {ch} on the line: Hessian coefficient 576 x_i N/(kappa + 224 sigma N), gradient coefficients"
                  f" -258048 x_i N^2/(..)^2 and 78336 x_i N^2/(..)^2, N = x^2 {'+' if eps == 1 else '-'} y^2 (independent of c)", ok)
        c1 = kap * t + x * (312 * t ** 2 + eps * 104 * u ** 2) + eps * 208 * y * t * u
        c13 = kap * u + 208 * x * t * u + y * (104 * t ** 2 + eps * 312 * u ** 2)
        okc = sp.expand(c1.subs(line) - x * (kap * s_ + 312 * N_ * s_ ** 2)) == 0 and \
            sp.expand(c13.subs(line) - y * (kap * s_ + 312 * N_ * s_ ** 2)) == 0
        check(f"{SG[eps]} zeroth order on the line: c_i = x_i (kappa sigma + 312 N sigma^2); so E_i = x_i E~ with ONE tensor"
              f" equation E~, and the curvature-free combination vanishes identically on the line", okc)
        dline = sp.factor(Dm.subs(line))
        check(f"{SG[eps]} det M on the line = {dline}: nonzero there iff (kappa + 224 N sigma)(kappa + 544 N sigma) != 0, so"
              f" r below is smooth wherever det M != 0",
              sp.expand(Dm.subs(line) - kap ** 2 * (kap + 224 * N_ * s_) * (kap + 544 * N_ * s_)) == 0)
        Ht, Qt, _ = out['E1']
        mup = sp.factor(x * (Qt - sp.diff(Ht, s_)) / Ht ** 2)
        Fp = 1 / (kap + 224 * s_ * N_) ** 2
        okF = sp.cancel(sp.diff(Fp, s_) / Fp - Qt / Ht) == 0
        F = -1 / (224 * N_ * (kap + 224 * s_ * N_))
        mu = sp.cancel(x * Fp / Ht)
        check(f"{SG[eps]} reparametrisation r = F(sigma) = -1/(224 N (kappa + 224 N sigma)) (F''/F' = Q~/H~): traceless part"
              f" (Hess r)_0 = mu (Ric)_0 with mu = x F'/H~ = -(7/18) r exactly; mu' = {mup}  [{time.time() - t1:.0f} s]",
              okF and sp.cancel(sp.diff(F, s_) - Fp) == 0 and sp.cancel(mu + sp.Rational(7, 18) * F) == 0
              and mup == sp.Rational(-7, 18))
    # the y = 0 case u = 0: Shiab#, W_A, the family (u = 0 is the line)
    for eps in (1, -1):
        cf = {}
        for nm, sub in (('Shiab#', {x: 1, y: 0, z: 1, w: 0}), ('W_A', {x: 1, y: 0, z: 0, w: 0})):
            r = reduced(subs_tables(tables(eps), sub))
            cf[nm] = [sp.factor(sp.cancel(sp.sympify(e).subs(u, 0))) for e in
                      (r['H'][0, 0], r['H'][1, 0], r['Q'][('E1', 't', 't')], r['Qp'][('E1', 't', 't')],
                       r['Q'][('E13', 't', 't')], r['Qp'][('E13', 't', 't')], r['Sig'][1, 0])]
        check(f"{SG[eps]} at u = 0 (du = 0): Shiab# and W_A have the same coefficients (H1t, H13t, Q1tt, Q'1tt, Q13tt, Q'13tt,"
              f" Sigma_Ct) = {cf['W_A']}", all(sp.cancel(a - b) == 0 for a, b in zip(cf['Shiab#'], cf['W_A']))
              and cf['W_A'][1] == 0 and cf['W_A'][4] == 0 and cf['W_A'][5] == 0)


# ====================================================================== L6
# --- plain-Fraction 4 x 4 matrices for the pairings
def mmul(A, B):
    return [[sum(A[i][k] * B[k][j] for k in range(4)) for j in range(4)] for i in range(4)]


def madd(A, B, s=1):
    return [[A[i][j] + s * B[i][j] for j in range(4)] for i in range(4)]


def msc(A, s):
    return [[A[i][j] * s for j in range(4)] for i in range(4)]


def mbr(A, B):
    return madd(mmul(A, B), mmul(B, A), -1)


def mtr2(A, B):
    return sum(A[i][k] * B[k][i] for i in range(4) for k in range(4))


def obasis(eta):
    """G_0-orthogonal basis X_a of p (eta-self-adjoint traceless 4 x 4), with G_aa = tr(X_a X_a)"""
    n = 4
    out = []
    for a in range(n):
        for b in range(a + 1, n):
            X = [[Fr(0)] * n for _ in range(n)]
            X[a][b] = Fr(eta[a])
            X[b][a] = Fr(eta[b])
            out.append(X)
    for dvec in ((1, -1, 0, 0), (1, 1, -2, 0), (1, 1, 1, -3)):
        out.append([[Fr(dvec[i]) if i == j else Fr(0) for j in range(n)] for i in range(n)])
    G = [mtr2(X, X) for X in out]
    ok = all(mtr2(out[i], out[j]) == 0 for i in range(9) for j in range(9) if i != j) and all(g != 0 for g in G)
    E = [[Fr(eta[i]) if i == j else Fr(0) for j in range(n)] for i in range(n)]
    ok &= all(mmul(mmul(E, list(map(list, zip(*X)))), E) == X for X in out)      # eta X^T eta = X (self-adjoint)
    ok &= all(sum(X[i][i] for i in range(n)) == 0 for X in out)
    return out, G, ok


def okbasis(eta):
    n = 4
    out = []
    for a in range(n):
        for b in range(a + 1, n):
            K = [[Fr(0)] * n for _ in range(n)]
            K[a][b] = Fr(eta[a])
            K[b][a] = Fr(-eta[b])
            out.append(K)                       # K = eta^-1 (e_ab - e_ba): K^T eta = -eta K
    return out


class Pairing:
    """<E, F> = sum_abc tr(E(X_a, X_b, X_c) F(X^a, X^b, X^c)) with the dual basis X^a = X_a / G_aa"""

    def __init__(self, eta):
        self.X, self.G, self.ok = obasis(eta)
        X = self.X
        self.trip = list(itertools.product(range(9), repeat=3))
        self.w = {abc: 1 / (self.G[abc[0]] * self.G[abc[1]] * self.G[abc[2]]) for abc in self.trip}
        self.c0 = {(a, b, c): msc(mbr(mbr(X[b], X[c]), X[a]), Fr(-1, 4)) for a, b, c in self.trip}
        self.c1 = {(a, b, c): madd(msc(X[c], mtr2(X[a], X[b])), msc(X[b], mtr2(X[a], X[c])), -1) for a, b, c in self.trip}

    def pair(self, E, F):
        return sum(self.w[k] * mtr2(E[k], F[k]) for k in self.trip)

    def quad(self, K1, K2, F):
        """<q(K1 (x) K2), F>, q(X, Y, Z) = d2A(X, Y)Z - d2A(X, Z)Y, d2A(X, Y)Z = 1/4 [tr(Z [K1, X]) [K2, Y] + tr(Z [K1, Y]) [K2, X]]"""
        X = self.X
        A1 = [mbr(K1, Xa) for Xa in X]
        A2 = [mbr(K2, Xa) for Xa in X]
        tau = [[mtr2(X[c], A1[a]) for a in range(9)] for c in range(9)]      # tau[c][a] = tr(X_c [K1, X_a])
        s = Fr(0)
        for (a, b, c) in self.trip:
            Fm = F[(a, b, c)]
            t2 = [mtr2(A2[d], Fm) for d in (a, b, c)]
            s += self.w[(a, b, c)] * (tau[c][a] * t2[1] + tau[c][b] * t2[0] - tau[b][a] * t2[2] - tau[b][c] * t2[0])
        return s / 4


def c1_tensor(eta, Bfun, mu1):
    """(C1)(X, Y, Z) = -1/4 [[Y, Z], X] + mu1 (tr(XY) Z - tr(XZ) Y) + d2A(X, Y)Z - d2A(X, Z)Y (sympy), with
    d2A(X, Y)Z = 1/4 sum_s coef_s [tr(Z [K_s, X]) [K'_s, Y] + tr(Z [K_s, Y]) [K'_s, X]] for B = sum_s coef_s K_s (x) K'_s"""
    def d2A(X, Y, Z):
        s = sp.zeros(4, 4)
        for coef, K1, K2 in Bfun:
            s += coef * ((Z * br(K1, X)).trace() * br(K2, Y) + (Z * br(K1, Y)).trace() * br(K2, X))
        return s / 4
    return lambda X, Y, Z: -br(br(Y, Z), X) / 4 + mu1 * ((X * Y).trace() * Z - (X * Z).trace() * Y) + d2A(X, Y, Z) - d2A(X, Z, Y)


def curv_B(eta, Rlow):
    n = 4
    Bm = {(c, d): sp.Matrix(n, n, lambda a, b: eta[a] * Rlow(a, b, c, d)) for c in range(n) for d in range(n)}
    return [(eta[i] * eta[j], Bm[(i, j)], Bm[(i, j)]) for i in range(n) for j in range(n)]


def rand_curv(rng, scale=3):
    n = 4
    pairs = [(a, b) for a in range(n) for b in range(a + 1, n)]
    Sx = {}
    for i in range(6):
        for j in range(i, 6):
            Sx[(i, j)] = Fr(rng.randint(-scale, scale), rng.randint(1, 2))
    Sx[(2, 3)] = -Sx[(0, 5)] + Sx[(1, 4)]     # first Bianchi: R0123 - R0213 + R0312 = 0

    def Rl(a, b, c, d):
        if a == b or c == d:
            return Fr(0)
        sg = 1
        if a > b:
            a, b, sg = b, a, -sg
        if c > d:
            c, d, sg = d, c, -sg
        i, j = pairs.index((a, b)), pairs.index((c, d))
        if j < i:
            i, j = j, i
        return sg * Sx[(i, j)]
    return Rl


def l6():
    section("L6. Lemma L7: the fibre obstruction at a point of s(X) (order eps^2 of Hess_F w = mu(w)(Ric_F + (A.,A.)))")
    for eta in ((1, 1, 1, 1), (1, -1, -1, -1), (-1, 1, 1, 1), (1, 1, -1, -1)):
        t1 = time.time()
        sig = (eta.count(1), eta.count(-1))
        Pr = Pairing(eta)
        p00, p01, p11 = Pr.pair(Pr.c0, Pr.c0), Pr.pair(Pr.c0, Pr.c1), Pr.pair(Pr.c1, Pr.c1)
        Kb = okbasis(eta)
        okl = True
        for K1, K2 in itertools.product(Kb, repeat=2):
            tk = mtr2(K1, K2)
            okl &= Pr.quad(K1, K2, Pr.c1) == -9 * tk and Pr.quad(K1, K2, Pr.c0) == Fr(-27, 4) * tk
        check(f"base {sig}: orthogonal basis of p ({Pr.ok}); <c0, c0> = {p00}, <c0, c1> = {p01}, <c1, c1> = {p11}; for all 36"
              f" basis elements K (x) K' of so (x) so: <q, c1> = -9 tr(K K'), <q, c0> = -(27/4) tr(K K'); hence"
              f" <q(B), 4 c0 - 3 c1> = 0 for every B in so (x) so and <(C1), 4 c0 - 3 c1> = -360 mu'  [{time.time() - t1:.0f} s]",
              Pr.ok and (p00, p01, p11) == (Fr(27, 2), 18, 144) and okl and 4 * p00 - 3 * p01 == 0 and 4 * p01 - 3 * p11 == -360)
        Rl = rand_curv(random.Random(60 + sum(eta)))
        E = [[Fr(eta[i]) if i == j else Fr(0) for j in range(4)] for i in range(4)]
        Bm = {(c, d): [[eta[a] * Rl(a, b, c, d) for b in range(4)] for a in range(4)] for c in range(4) for d in range(4)}
        R2 = sum(eta[a] * eta[b] * eta[c] * eta[d] * Rl(a, b, c, d) ** 2 for a, b, c, d in itertools.product(range(4), repeat=4))
        sB = sum(eta[i] * eta[j] * mtr2(Bm[(i, j)], Bm[(i, j)]) for i in range(4) for j in range(4))
        okso = all(mmul(list(map(list, zip(*B_))), E) == msc(mmul(E, B_), -1) for B_ in Bm.values())
        check(f"base {sig}: B_ij = R(e_i, e_j) lie in so(eta) ({okso}); sum_ij eps_i eps_j tr(B_ij^2) = -|R|^2 (random R:"
              f" {sB} = -({R2})), so <q, c1> = 9 |R|^2 and <q, c0> = (27/4) |R|^2 for every curvature tensor", okso and sB == -R2)
    # span certificate: the coefficient vectors of (C1) over the quadratic monomials of a generic curvature tensor
    for eta in ((1, 1, 1, 1), (1, -1, -1, -1), (1, 1, -1, -1)):
        t1 = time.time()
        n = 4
        pairs = [(a, b) for a in range(n) for b in range(a + 1, n)]
        syms = {}
        for i in range(6):
            for j in range(i, 6):
                syms[(i, j)] = sp.Symbol(f'R_{i}{j}')
        bian = {syms[(2, 3)]: -syms[(0, 5)] + syms[(1, 4)]}
        free = [s for k, s in syms.items() if k != (2, 3)]

        def Rs(a, b, c, d):
            if a == b or c == d:
                return 0
            sg = 1
            if a > b:
                a, b, sg = b, a, -sg
            if c > d:
                c, d, sg = d, c, -sg
            i, j = pairs.index((a, b)), pairs.index((c, d))
            if j < i:
                i, j = j, i
            return sg * syms[(i, j)].subs(bian)
        m1 = sp.Symbol('m1')
        T = c1_tensor(eta, curv_B(eta, Rs), m1)
        Pb = pbasis(eta)
        eqs = []
        for X in Pb:
            for iy, Yv in enumerate(Pb):
                for iz, Zv in enumerate(Pb):
                    if iz <= iy:
                        continue
                    for e in T(X, Yv, Zv):
                        e = sp.expand(e)
                        if e != 0:
                            eqs.append(e)
        one = tuple([0] * len(free))
        polys = [sp.Poly(e, *free) for e in eqs]
        mons = sorted({m for p in polys for m in p.monoms()} - {one})
        rows, cc0, cc1 = [], [], []
        for p in polys:
            dct = dict(zip(p.monoms(), p.coeffs()))
            rows.append([dct.get(m, 0) for m in mons])
            c = sp.sympify(dct.get(one, 0))
            cc0.append(c.subs(m1, 0))
            cc1.append(c.coeff(m1))
        Dq = DomainMatrix([[QQ(v) for v in r] for r in rows], (len(rows), len(mons)), sp.QQ)
        D0 = DomainMatrix([[QQ(v) for v in r] + [QQ(a)] for r, a in zip(rows, cc0)], (len(rows), len(mons) + 1), sp.QQ)
        D01 = DomainMatrix([[QQ(v) for v in r] + [QQ(a), QQ(b)] for r, a, b in zip(rows, cc0, cc1)], (len(rows), len(mons) + 2), sp.QQ)
        rq, r0, r01 = Dq.rank(), D0.rank(), D01.rank()
        check(f"base {(eta.count(1), eta.count(-1))}: (C1) for the generic curvature tensor (20 components, {len(eqs)} scalar"
              f" equations): rank of the quadratic parts {rq}; with c0 added {r0}; with c1 added {r01}: c0 lies in the span"
              f" and c1 does not, so (C1) has a solution R only if mu' = 0  [{time.time() - t1:.0f} s]",
              rq == 21 and r0 == 21 and r01 == 22)
    # brute force: Taylor solution of the fibre equation at y = eta (independent of the (C1) derivation)
    R1 = lambda eta: (lambda a, b, c, d: Fr((eta[a] if a == c else 0) * (eta[b] if b == d else 0)  # noqa: E731
                                            - (eta[a] if a == d else 0) * (eta[b] if b == c else 0)))
    for eta, lab, Rl, mu1, af, expect in (
            ((1, -1, -1, -1), 'random R', rand_curv(random.Random(71)), Fr(-7, 18), Fr(1), False),
            ((1, 1, 1, 1), 'random R', rand_curv(random.Random(72)), Fr(-7, 18), Fr(1), False),
            ((1, 1, -1, -1), 'random R', rand_curv(random.Random(73)), Fr(-7, 18), Fr(1), False),
            ((1, -1, -1, -1), "control: R = R_1, A-term times -1/12, mu' = 0", R1((1, -1, -1, -1)), Fr(0), Fr(-1, 12), True),
            ((1, -1, -1, -1), "control: R = R_1, A-term times -1/12, mu' = -7/18", R1((1, -1, -1, -1)), Fr(-7, 18), Fr(-1, 12), False)):
        t1 = time.time()
        fd = fibre_data(eta, Rl, maxorder=2, afac=af)
        unk, rowkeys, cols, _, rhs = hessian_system(fd, 2, 4, mu0=Fr(1), mu1=mu1)
        r1, r2 = dm_from_cols(len(rowkeys), cols).rank(), dm_from_cols(len(rowkeys), cols, rhs).rank()
        check(f"brute force, base {(eta.count(1), eta.count(-1))}, {lab}, mu(w) = 1 + mu' (w - 1), w(eta) = 1: the Taylor system"
              f" for w to order eps^2 ({len(cols)} unknowns) is {'consistent' if r1 == r2 else 'inconsistent'} (ranks {r1}, {r2}),"
              f" as (C1) predicts  [{time.time() - t1:.0f} s]", (r1 == r2) == expect)
    # order eps: consistent, and every solution has dw(eta) = 0
    t1 = time.time()
    eta = (1, -1, -1, -1)
    fd = fibre_data(eta, rand_curv(random.Random(74)), maxorder=1)
    unk, rowkeys, cols, _, rhs = hessian_system(fd, 1, 3, mu0=Fr(1), mu1=Fr(-7, 18))
    Da = dm_from_cols(len(rowkeys), cols, rhs)
    ra = Da.rank()
    okc = dm_from_cols(len(rowkeys), cols).rank() == ra
    Ma = Da.to_Matrix()
    okz = True
    for m in fd['Pp'].monos(1):
        e = sp.zeros(1, Ma.shape[1])
        e[0, unk.index(m)] = 1
        okz &= DomainMatrix.from_Matrix(Ma.col_join(e)).convert_to(sp.QQ).rank() == ra
    check(f"brute force, base (1, 3), random R, mu' = -7/18: the order-eps system is consistent, and adding dw(eta) = 0"
          f" changes no rank, i.e. every solution has vanishing vertical gradient at p  [{time.time() - t1:.0f} s]", okc and okz)


def main():
    l0()
    l1()
    l2()
    l3()
    l4()
    l5()
    l6()
    P(f"\n{NCHK[0]} checks, all True (each asserted).")
    P(f"runtime {time.time() - T0:.0f} s")


if __name__ == '__main__':
    main()
