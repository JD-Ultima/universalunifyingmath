"""GU round 6 (J): the slaved even sector.  Does the elliptic twist (z, w) change the reduced equations once T_e != 0?

Formalization: formalization/geometric-unity/J-slaved-sector.md.  Builds on I (I2: parity split, I2.2: slaving
T_e = -L_{T_o}^{-1} Shiab(d_0 T_o); I3: twist identity; I6: invariant sector), H (H3 coordinates, H7 family), F (F2-F6
closed forms along s(X)).  Imports I_elliptic_reduction.py (and through it G_symbols.py, H_ellipticity.py,
shiab_checks.py), unmodified.

Setting (as in I).  d = 14, eta = (+^14) and (+^7,-^7); g = Cl(eta) (x) C; unitary gauge, A = nabla^0 + T,
F_A = F_0 + d_0 T + T ^ T.  Upsilon(alpha) = sum_a Tr(D_a alpha_a); we work with the 1-form  Y_a = eta_aa D_a
(Upsilon = *Y), so that the kappa-term is Y = kappa T and the I6 algebraic part is
   Y_odd = -[x G - c1 g]_ab e^a (x) gamma^b - [y G - c13 g]_ab e^a (x) gamma^b Gamma.
Family (H7, code normalisation): Shiab = x W_A + y W'_A + z W_C.AdG + w W'_C.AdG; Shiab# = (1, 0, 1, 0).
Ansatz (the "vector sector").  T_o = t Phi^1 + u Phi^1 Gamma (t, u functions; I6) and T_e = sum_k P_k(A_k) with
   P_0(v)_a = v_a 1,  P_14(v)_a = v_a Gamma,  P_2(v)_a = sum_{b != a} v_b gamma_a gamma^b,  P_12(v) = P_2(v) Gamma,
A_k covector fields (k in {0, 2, 12, 14}; names A, C, D, B).  These are the four copies of the vector representation in
Lambda^1 (x) g_ev.
Sections
  J0  re-derivation of c1, c13 (I6) with the same evaluator.
  J1  closure of the even equation on the vector sector: Shiab(v ^ Phi^1), Shiab(v ^ Phi^1 Gamma) and 2B(P_k(v), Phi^1),
      2B(P_k(v), Phi^1 Gamma) are exact combinations of P_j(v) for every basis covector v (hence all v): source matrix S,
      slaving block M(t,u) = kappa 1 + t M1 + u M13, det M, slaved T_e = Sigma (dt, du), Sigma = -M^-1 S.
  J2  the odd equation on the sector: the derivative map nabla A_k -> Shiab(d_0 P_k(A_k)) on all 196 basis tensors and the
      quadratic map (A_k, A_l) -> 2B(P_k(A_k), P_l(A_l)) on all basis pairs, exact tables; assembled reduced equation;
      kinetic matrix H (exact and to O(1/kappa)), its determinant; Lambda^2-type equations.
  J3  end-to-end check on random 2-jets (fixed seed, exact Gaussian rationals, numeric family members): the slaved T_e
      solves the even equation exactly, and the directly evaluated odd part equals the assembled formula.
  J4  reduction along s(X) (F3-F6 closed forms, n = 4): (P) and (E) for basic t, u; (E) for the isotropy-invariant 2-jet
      ansatz (including the radial sub-case); balance identity; exact sample solutions in (14,0) and (7,7).
  J4b (after review REVIEW4-GU, R6-1) the 2-jet ansatz with a normal gradient a != 0: closed form for every a, the (14,0)
      bounds for Shiab# in Z1 and Z2, the referee's Z1 configuration, an exact solution attaining the Z1 bound.
  J5  the constant-curvature leaf test of Theorem J5(d), bases (4,0), (1,3), (2,2).
  Also added after the review: exact det H for every member (Jacobi), the exact image condition, exact varpi for W_A (J2).
Exact arithmetic throughout (integers, Fractions, Gaussian rationals, sympy).  Seeds fixed.  One core.
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

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import I_elliptic_reduction as I  # noqa: E402
from G_symbols import FULL, N  # noqa: E402

T0 = time.time()
random.seed(20261030)
pc = lambda v: bin(v).count('1')  # noqa: E731
EUC, E77, SIG = I.EUC, I.E77, I.SIG
OPN = I.OPN
UNIT = {op: I.unit_sp(I.OPS[op][0]) for op in OPN}
xs, ys, zs, ws = I.xs, I.ys, I.zs, I.ws
FAM = I.FAMC
t_, u_, kap = sp.symbols('t u kappa')
KS = (0, 2, 12, 14)
KNM = {0: 'A', 2: 'C', 12: 'D', 14: 'B'}
CH = ('E1', 'E13', 'w3', 'w11')
FAST = os.environ.get('J_FAST') == '1'      # debugging only: orbit representatives instead of all basis tensors/pairs
FAST_REPS = {1: {(0, 0), (0, 1), (1, 0)}, -1: {(0, 0), (7, 7), (0, 1), (1, 0), (7, 8), (8, 7), (0, 7), (7, 0)}}


def P_(*a):
    print(*a)
    sys.stdout.flush()


class GI:
    """exact Gaussian integers on Python ints (end-to-end check of J3)"""
    __slots__ = ('r', 'i')

    def __init__(s, r=0, i=0):
        s.r, s.i = r, i

    def __add__(s, o):
        return GI(s.r + o.r, s.i + o.i) if isinstance(o, GI) else GI(s.r + o, s.i)
    __radd__ = __add__

    def __neg__(s):
        return GI(-s.r, -s.i)

    def __sub__(s, o):
        return s + (-o)

    def __mul__(s, o):
        if isinstance(o, GI):
            return GI(s.r * o.r - s.i * o.i, s.r * o.i + s.i * o.r)
        return GI(s.r * o, s.i * o)
    __rmul__ = __mul__

    def __eq__(s, o):
        return (s.r == o.r and s.i == o.i) if isinstance(o, GI) else (s.i == 0 and s.r == o)

    def __ne__(s, o):
        return not s.__eq__(o)

    def __hash__(s):
        return hash((s.r, s.i))


def gi(v):
    """sympy Gaussian integer -> GI (integrality asserted)"""
    v = sp.expand_complex(sp.sympify(v))
    re, im = sp.re(v), sp.im(v)
    assert re.is_integer and im.is_integer, v
    return GI(int(re), int(im))


class Ctx:
    """per-signature algebra: patterns, Hodge-index raising, odd-pattern extraction"""

    def __init__(self, eta):
        self.eta = eta
        self.Al = I.Alg(eta)
        self.G = self.Al.Gam()
        self._g3 = {}
        self.Phi1 = [self.Al.gamma(a) for a in range(N)]
        self.Phi1G = [self.Al.gammaG(a) for a in range(N)]

    def e(self, p):
        return [1 if i == p else 0 for i in range(N)]

    def g3(self, a, b, c):
        k = (a, b, c)
        if k not in self._g3:
            X = self.Al.mul(self.Al.mul(self.Al.gamma(a), self.Al.gamma(b)), self.Al.gamma(c))
            (m, s), = X.items()
            self._g3[k] = (m, s)
        return self._g3[k]

    def P(self, k, v):
        """P_k(v) as a g-valued 1-form (list of 14 Clifford elements)"""
        Al, eta = self.Al, self.eta
        out = []
        for a in range(N):
            X = {}
            if k in (0, 14):
                if v[a] != 0:
                    X = {0: v[a]} if k == 0 else {FULL: v[a]}
            else:
                for b in range(N):
                    if b != a and v[b] != 0:
                        X = I.add(X, I.scale(Al.mul(Al.gamma(a), Al.gamma(b)), v[b] * eta[b]))
                if k == 12:
                    X = Al.mul(X, self.G)
            out.append(X)
        return out

    def hat(self, D, div=1):
        return [{k: (Fr(v) if isinstance(v, int) else v) * self.eta[a] / div for k, v in D[a].items()}
                for a in range(N)]

    def fitV(self, Y, c):
        """Y (1-form) = sum_k s_k P_k(e^c) exactly?  returns ({k: s_k}, ok)"""
        s = {}
        for k in KS:
            pat = self.P(k, self.e(c))
            a0 = next(a for a in range(N) if pat[a])
            m0 = next(iter(pat[a0]))
            s[k] = Y[a0].get(m0, 0) / pat[a0][m0]
        rec = [dict() for _ in range(N)]
        for k in KS:
            if s[k] != 0:
                pat = self.P(k, self.e(c))
                for a in range(N):
                    I.acc_into(rec[a], pat[a], s[k])
        return s, all(I.add(rec[a], Y[a], -1) == {} for a in range(N))

    def extract(self, Y):
        """odd 1-form Y -> (E1, E13, w3, w11):  Y_a = sum_b E1_ab gamma^b + sum_b E13_ab gamma^b Gamma
        + sum_{b,c} w3_bc gamma_a gamma^b gamma^c + sum_{b,c} w11_bc gamma_a gamma^b gamma^c Gamma; exact reconstruction"""
        Al, eta = self.Al, self.eta
        E1 = [[0] * N for _ in range(N)]
        E13 = [[0] * N for _ in range(N)]
        w3 = [[0] * N for _ in range(N)]
        w11 = [[0] * N for _ in range(N)]
        for a in range(N):
            for m, v in Y[a].items():
                g = pc(m)
                if g == 1:
                    b = m.bit_length() - 1
                    E1[a][b] = v * eta[b]
                elif g == 13:
                    b = (FULL ^ m).bit_length() - 1
                    (_, sg), = Al.gammaG(b).items()
                    E13[a][b] = v * sg * eta[b]
        for a in range(N):
            for b, c in itertools.combinations([i for i in range(N) if i != a], 2):
                m, sg = self.g3(a, b, c)
                sg *= eta[b] * eta[c]
                v = Y[a].get(m, 0)
                if v != 0 and w3[b][c] == 0:
                    w3[b][c] = v / (2 * sg)
                    w3[c][b] = -w3[b][c]
                s2 = Al.cl.sign(m, FULL) * sg
                v2 = Y[a].get(m ^ FULL, 0)
                if v2 != 0 and w11[b][c] == 0:
                    w11[b][c] = v2 / (2 * s2)
                    w11[c][b] = -w11[b][c]
        rec = self.build(E1, E13, w3, w11)
        ok = all(I.add(rec[a], Y[a], -1) == {} for a in range(N))
        return (E1, E13, w3, w11), ok

    def build(self, E1, E13, w3, w11):
        Al, eta = self.Al, self.eta
        R = [dict() for _ in range(N)]
        for a in range(N):
            for b in range(N):
                if E1[a][b] != 0:
                    I.acc_into(R[a], I.scale(Al.gamma(b), E1[a][b] * eta[b]))
                if E13[a][b] != 0:
                    I.acc_into(R[a], I.scale(Al.gammaG(b), E13[a][b] * eta[b]))
                for c in range(N):
                    if a in (b, c) or b == c:
                        continue
                    m, sg = self.g3(a, b, c)
                    if w3[b][c] != 0:
                        I.acc_into(R[a], {m: w3[b][c] * eta[b] * eta[c] * sg})
                    if w11[b][c] != 0:
                        I.acc_into(R[a], {m ^ FULL: w11[b][c] * eta[b] * eta[c] * sg * Al.cl.sign(m, FULL)})
        return R

    def dzeta(self, k, Nmat):
        """2-form d_0 P_k(A) for nabla_c A_b = Nmat[c][b]:  zeta_(c,a) = P_k(N_c.)_a - P_k(N_a.)_c  (c < a)"""
        rows = {c: self.P(k, Nmat[c]) for c in range(N) if any(v != 0 for v in Nmat[c])}
        z = {}
        for c, a in itertools.combinations(range(N), 2):
            X = {}
            if c in rows and rows[c][a]:
                X = I.add(X, rows[c][a])
            if a in rows and rows[a][c]:
                X = I.add(X, rows[a][c], -1)
            if X:
                z[(c, a)] = X
        return z

    def vwedge(self, v, P1):
        """2-form v ^ P1 (P1 a g-valued 1-form, v a covector)"""
        z = {}
        for c, a in itertools.combinations(range(N), 2):
            X = {}
            if v[c] != 0 and P1[a]:
                X = I.add(X, I.scale(P1[a], v[c]))
            if v[a] != 0 and P1[c]:
                X = I.add(X, I.scale(P1[c], -v[a]))
            if X:
                z[(c, a)] = X
        return z


def per_op(pt):
    return {op: I.combine(pt, I.OPS[op][1]) for op in OPN}


def q2B(Al, X, Y, cache=None):
    """per-type dual of 3 * 2B(X, Y) = 3Q(X+Y) - 3Q(X) - 3Q(Y)"""
    qs = Al.quad3([I.add(X[a], Y[a]) for a in range(N)])
    qx = Al.quad3(X) if cache is None or 'x' not in cache else cache['x']
    qy = Al.quad3(Y) if cache is None or 'y' not in cache else cache['y']
    z = [dict() for _ in range(N)]
    return {t: [I.add(I.add(qs.get(t, z)[a], qx.get(t, z)[a], -1), qy.get(t, z)[a], -1) for a in range(N)]
            for t in set(qs) | set(qx) | set(qy)}


def famsum(d):
    """{op: rational} -> sympy linear form in x, y, z, w (with the operator units)"""
    return sp.expand(sum(UNIT[op] * FAM[op] * sp.Rational(Fr(d[op]).numerator, Fr(d[op]).denominator)
                         for op in OPN if d.get(op, 0) != 0))


# ====================================================================== J0
def j0(C):
    Al, eta = C.Al, C.eta
    both = [I.add(C.Phi1[a], C.Phi1G[a]) for a in range(N)]
    vals = {}
    for lab, T in (('tt', C.Phi1), ('uu', C.Phi1G), ('sum', both)):
        po = per_op(Al.quad3(T))
        for op in OPN:
            D = po[op]
            c1 = Fr(D[0].get(1, 0) * eta[0], 3)
            (k0, s0), = Al.gammaG(0).items()
            c13 = Fr(D[0].get(k0, 0) * s0 * eta[0], 3)
            vals[(lab, op)] = (c1, c13)
    c1, c13 = kap * t_, kap * u_
    for op in OPN:
        tt, uu, sm = vals[('tt', op)], vals[('uu', op)], vals[('sum', op)]
        r = lambda q: sp.Rational(q.numerator, q.denominator)  # noqa: E731
        c1 += UNIT[op] * FAM[op] * (r(tt[0]) * t_ ** 2 + r(uu[0]) * u_ ** 2 + r(sm[0] - tt[0] - uu[0]) * t_ * u_)
        c13 += UNIT[op] * FAM[op] * (r(tt[1]) * t_ ** 2 + r(uu[1]) * u_ ** 2 + r(sm[1] - tt[1] - uu[1]) * t_ * u_)
    c1, c13 = sp.expand(c1), sp.expand(c13)
    sg = eta[7]
    ref1 = kap * t_ + xs * (312 * t_ ** 2 + sg * 104 * u_ ** 2) + sg * 208 * ys * t_ * u_
    ref13 = kap * u_ + 208 * xs * t_ * u_ + ys * (104 * t_ ** 2 + sg * 312 * u_ ** 2)
    ok = sp.expand(c1 - ref1) == 0 and sp.expand(c13 - ref13) == 0
    P_(f"J0 {SIG[eta[7]]}: c1 = {c1};  c13 = {c13}  (= I6: {ok})")
    assert ok
    return c1, c13


# ====================================================================== J1
def j1(C):
    Al, eta = C.Al, C.eta
    P_(f"\n=== J1 {SIG[eta[7]]}. Closure of the even equation on the vector sector (exact, every basis covector e^c) ===")
    src = {}
    ok_all = True
    for lab, P1 in (('t', C.Phi1), ('u', C.Phi1G)):
        ref = None
        for c in range(N):
            po = per_op(Al.psi_form(C.vwedge(C.e(c), P1)))
            fits = {}
            for op in OPN:
                s, ok = C.fitV(C.hat(po[op]), c)
                ok_all &= ok
                fits[op] = s
            ok_all &= ref is None or fits == ref
            ref = fits
        src[lab] = ref
    S = sp.Matrix([[famsum({op: src[lab][op][k] for op in OPN}) for lab in ('t', 'u')] for k in KS])
    P_(f"  source Y_ev[Shiab(dt ^ Phi^1 + du ^ Phi^1 Gamma)] = sum_k P_k(S_k,t dt + S_k,u du), rows k = (0, 2, 12, 14):"
       f" S = {S.tolist()}")
    blocks = {}
    for lab, To in (('t', C.Phi1), ('u', C.Phi1G)):
        qT = Al.quad3(To)
        ref = None
        for c in range(N):
            fits = {}
            for k in KS:
                X = C.P(k, C.e(c))
                po = per_op(q2B(Al, X, To, {'y': qT}))
                for op in OPN:
                    s, ok = C.fitV(C.hat(po[op], 3), c)
                    ok_all &= ok
                    fits[(op, k)] = s
            ok_all &= ref is None or fits == ref
            ref = fits
        blocks[lab] = ref
    M1 = sp.Matrix([[famsum({op: blocks['t'][(op, k)][j] for op in OPN}) for k in KS] for j in KS])
    M13 = sp.Matrix([[famsum({op: blocks['u'][(op, k)][j] for op in OPN}) for k in KS] for j in KS])
    M = kap * sp.eye(4) + t_ * M1 + u_ * M13
    P_(f"  Y_ev[kappa T_e + 2B(T_e, T_o)] = sum_j P_j((M A)_j), A = (A_0, A_2, A_12, A_14):  M = kappa 1 + t M1 + u M13,")
    P_(f"   M = {M.tolist()}")
    P_(f"  every output is exactly of this form for all 14 basis covectors, both maps, all four operators: {ok_all}")
    assert ok_all
    sgn = eta[7]
    sig_, th_, de_ = t_ * xs + sgn * u_ * ys, t_ * ws - u_ * zs, t_ * ys - u_ * xs
    detM = sp.expand(M.det())
    P_(f"  det M = {sp.factor(detM)}")
    # compact form of M with sigma = t x + eps u y, theta = t w - u z, delta = t y - u x  (eps = +1 (14,0), -1 (7,7))
    s1, s2, s3 = sp.symbols('sigma theta delta')
    Mc = sp.Matrix([[kap, 0, 208 * s3 * sgn, 0],
                    [0, kap + 224 * s1, 384 * sp.I * s2 * sgn, 16 * s3 * sgn],
                    [16 * s3, -384 * sp.I * s2, kap + 544 * s1, 0],
                    [0, 208 * s3, 0, kap]])
    okc = sp.expand(Mc.subs({s1: sig_, s2: th_, s3: de_}) - M) == sp.zeros(4, 4)
    dc = sp.factor(sp.expand(Mc.det()))
    P_(f"  M = [[k, 0, 208 e d, 0], [0, k + 224 s, 384 i e th, 16 e d], [16 d, -384 i th, k + 544 s, 0], [0, 208 d, 0, k]]"
       f" with s = tx + e uy, th = tw - uz, d = ty - ux, e = {sgn:+d}: {okc};  det M = {dc}")
    assert okc
    rho_ = kap ** 2 - 3328 * sgn * de_ ** 2
    okd = sp.expand(detM - ((rho_ + 224 * kap * sig_) * (rho_ + 544 * kap * sig_) - sgn * (384 * kap * th_) ** 2)) == 0
    P_(f"  det M = (r + 224 k s)(r + 544 k s) - e (384 k th)^2 with r = k^2 - 3328 e d^2: {okd}")
    assert okd
    Sig = -M.adjugate() * S / detM
    lo = -S / kap
    P_("  slaved T_e: (A_0, A_2, A_12, A_14) = Sigma (dt, du)^T, Sigma = -M^-1 S;  to leading order Sigma = -S/kappa:"
       f" A_2 = C = {sp.factor(lo[1, 0])} dt + {sp.factor(lo[1, 1])} du,  A_12 = D = {sp.factor(lo[2, 0])} dt + "
       f"{sp.factor(lo[2, 1])} du;  A_0 = A and A_14 = B are O(kappa^-2)")
    for nm, sub in (('Shiab#', {xs: 1, zs: 1, ys: 0, ws: 0}), ('W_A', {xs: 1, zs: 0, ys: 0, ws: 0})):
        Ss = sp.simplify(Sig.subs(sub))
        P_(f"  {nm}: Sigma (exact) = {[[sp.factor(Ss[i, j]) for j in range(2)] for i in range(4)]}")
    zw_in = any({zs, ws} & set(sp.expand(Sig[i, j] * detM).free_symbols) for i in range(4) for j in range(2))
    P_(f"  z, w occur in the slaved T_e: {zw_in} (through the D-row of S, already at O(1/kappa), and through M)")
    return dict(S=S, M=M, M1=M1, M13=M13, detM=detM, Sig=Sig)


# ====================================================================== J2
def j2(C, d1):
    Al, eta = C.Al, C.eta
    P_(f"\n=== J2 {SIG[eta[7]]}. The odd equation on the vector sector (exact) ===")
    # ---- derivative map
    lin = {}
    ok_all = True
    for k in KS:
        data = {}
        for p in range(N):
            for q in range(N):
                if FAST and (p, q) not in FAST_REPS[eta[7]]:
                    continue
                Nm = [[0] * N for _ in range(N)]
                Nm[p][q] = 1
                po = per_op(Al.psi_form(C.dzeta(k, Nm)))
                for op in OPN:
                    data[(op, p, q)] = C.extract(C.hat(po[op]))
        for op in OPN:
            (E1, E13, w3, w11), _ = data[(op, 0, 0)]
            (F1, F13, v3, v11), _ = data[(op, 0, 1)]
            res = {}
            for nm, Ea, Fa in (('E1', E1, F1), ('E13', E13, F13)):
                res[nm] = (Fa[0][1], Fa[1][0], Ea[1][1] * eta[0] * eta[1])     # (alpha, alpha', alpha'')
            res['w3'], res['w11'] = v3[0][1], v11[0][1]
            for p in range(N):
                for q in range(N):
                    if (op, p, q) not in data:
                        continue
                    (E1, E13, w3, w11), okx = data[(op, p, q)]
                    ok_all &= okx
                    for nm, Ea in (('E1', E1), ('E13', E13)):
                        al0, al1, al2 = res[nm]
                        ok_all &= all(Ea[a][b] == al0 * (a == p and b == q) + al1 * (a == q and b == p) +
                                      al2 * (eta[p] if p == q else 0) * (eta[a] if a == b else 0)
                                      for a in range(N) for b in range(N))
                    for nm, wa in (('w3', w3), ('w11', w11)):
                        ok_all &= all(wa[a][b] == res[nm] * ((a == p and b == q) - (a == q and b == p))
                                      for a in range(N) for b in range(N))
            lin[(op, k)] = res
    struct = all(lin[(op, k)][nm][0] == 0 and lin[(op, k)][nm][2] == -lin[(op, k)][nm][1]
                 for op in OPN for k in KS for nm in ('E1', 'E13'))
    L = {}
    for k in KS:
        for nm in ('E1', 'E13'):
            L[(nm, k)] = famsum({op: lin[(op, k)][nm][1] for op in OPN})
        for nm in ('w3', 'w11'):
            L[(nm, k)] = famsum({op: lin[(op, k)][nm] for op in OPN})
    P_("  derivative map, N_k = nabla A_k (N_ab = nabla_a A_b):  Y_odd[Shiab(d_0 P_k(A_k))] has E1, E13 = l_k (N_k^T - tr N_k g)"
       " (alpha = 0, alpha'' = -alpha' for all k and all four operators: " + str(struct) + "),  w3, w11 = l_k dA_k:")
    for k in KS:
        P_(f"    k = {k:2d} ({KNM[k]}): " + ";  ".join(f"{nm}: {L[(nm, k)]}" for nm in CH))
    P_(f"  checked exactly on all 196 basis tensors e^p (x) e^q, k = 0, 2, 12, 14, four operators: {ok_all}")
    assert ok_all and struct
    # ---- quadratic map
    bil = {}
    okq = True
    singles = {}
    for k in KS:
        for p in range(N):
            if FAST and p not in (0, 1, 7, 8):
                continue
            singles[(k, p)] = Al.quad3(C.P(k, C.e(p)))
    for k, l in itertools.combinations_with_replacement(KS, 2):
        data = {}
        pairs = [(p, q) for p in range(N) for q in range(N) if (k != l or p <= q)]
        if FAST:
            pairs = [pq for pq in pairs if pq in FAST_REPS[eta[7]]]
        for p, q in pairs:
            X, Y = C.P(k, C.e(p)), C.P(l, C.e(q))
            po = per_op(q2B(Al, X, Y, {'x': singles[(k, p)], 'y': singles[(l, q)]}))
            for op in OPN:
                data[(op, p, q)] = C.extract(C.hat(po[op], 3))
        for op in OPN:
            (E1, E13, w3, w11), _ = data[(op, 0, 1)]
            (F1, F13, v3, v11), _ = data[(op, 0, 0)]
            res = {}
            for nm, Ea, Fa in (('E1', E1, F1), ('E13', E13, F13)):
                res[nm] = (Ea[0][1], Ea[1][0], Fa[1][1] * eta[1] * eta[0])      # (b1, b2, b3)
            res['w3'], res['w11'] = w3[0][1], w11[0][1]
            for p, q in pairs:
                (E1, E13, w3, w11), okx = data[(op, p, q)]
                okq &= okx
                for nm, Ea in (('E1', E1), ('E13', E13)):
                    b1, b2, b3 = res[nm]
                    okq &= all(Ea[a][b] == b1 * (a == p and b == q) + b2 * (a == q and b == p) +
                               b3 * (eta[p] if p == q else 0) * (eta[a] if a == b else 0)
                               for a in range(N) for b in range(N))
                for nm, wa in (('w3', w3), ('w11', w11)):
                    okq &= all(wa[a][b] == res[nm] * ((a == p and b == q) - (a == q and b == p))
                               for a in range(N) for b in range(N))
            bil[(op, k, l)] = res
    B = {}
    for k, l in itertools.combinations_with_replacement(KS, 2):
        for nm in ('E1', 'E13'):
            B[(nm, k, l)] = tuple(famsum({op: bil[(op, k, l)][nm][j] for op in OPN}) for j in range(3))
        for nm in ('w3', 'w11'):
            B[(nm, k, l)] = famsum({op: bil[(op, k, l)][nm] for op in OPN})
    P_("  quadratic map: Y_odd[2B(P_k(A), P_l(A'))] has E = b1 A(x)A' + b2 A'(x)A + b3 <A,A'> g,  w = b4 A ^ A' (nonzero entries):")
    for k, l in itertools.combinations_with_replacement(KS, 2):
        ent = [f"{nm}: {B[(nm, k, l)]}" for nm in ('E1', 'E13') if any(v != 0 for v in B[(nm, k, l)])]
        ent += [f"{nm}: {B[(nm, k, l)]}" for nm in ('w3', 'w11') if B[(nm, k, l)] != 0]
        P_(f"    ({KNM[k]}, {KNM[l]}): " + ("; ".join(ent) if ent else "0"))
    twist_only = all(bil[(op, k, l)] == {'E1': (0, 0, 0), 'E13': (0, 0, 0), 'w3': 0, 'w11': 0}
                     for op in ('W_C.AdG', "W'_C.AdG") for k, l in itertools.combinations_with_replacement(KS, 2)
                     if (k, l) != (2, 12))
    P_(f"  checked exactly on every basis pair (e^p, e^q) (p <= q on the diagonal blocks), four operators: {okq};"
       f"  the twisted pair enters Q(T_e) only through the (C, D) cross term: {twist_only}")
    assert okq and twist_only
    # ---- assembly
    S, M, M1, M13, Sig = d1['S'], d1['M'], d1['M1'], d1['M13'], d1['Sig']
    lE = sp.Matrix([[L[(nm, k)] for k in KS] for nm in ('E1', 'E13')])
    H = lE * Sig
    H0 = -lE * S / kap
    P_("  assembled (A_k = Sigma_k,t dt + Sigma_k,u du):  E^(i) = -x_i G + c_i g + sum_phi H_i,phi (nabla^2 phi - Box phi g)"
       " + sum Q_i,phi,psi dphi(x)dpsi + sum Q'_i,phi,psi <dphi,dpsi> g,  (x_1, x_13) = (x, y);"
       "  H = l_E Sigma (2 x 2, rows E1, E13, columns t, u)")
    P_(f"  H to leading order: H = -l_E S / kappa = {[[sp.factor(H0[i, j]) for j in range(2)] for i in range(2)]},"
       f"  det = {sp.factor(H0.det())}")
    KK = sp.Matrix([[1, 0], [0, eta[7]]]) * H0          # channel-13 row times its Gram sign in (7,7)
    P_(f"  symmetric (Gram-weighted) form diag(1, {eta[7]:+d}) H_0 = {[[sp.factor(KK[i, j]) for j in range(2)] for i in range(2)]},"
       f" symmetric: {sp.expand(KK - KK.T) == sp.zeros(2, 2)}")
    if eta[7] == 1:
        xr, yr, zr, wr = sp.symbols("x' y' z' w'", real=True)
        Kr = sp.expand(-KK.subs({xs: -sp.I * xr, ys: -sp.I * yr, zs: -sp.I * zr, ws: -sp.I * wr}))
        P_(f"  (14,0) u(128)-real form (x = -i x', ..., t = i s, u = i r; equation multiplied by -i): Hessian coefficients"
           f" -K with K = {Kr.tolist()} = (576/kappa) Gram((x', w'), (y', -z')), det K = {sp.factor(Kr.det())}")
    for nm, sub in (('Shiab#', {xs: 1, zs: 1, ys: 0, ws: 0}), ('W_A', {xs: 1, zs: 0, ys: 0, ws: 0})):
        Hs = sp.simplify(H.subs(sub))
        P_(f"  {nm}: H (exact) = {[[sp.factor(Hs[i, j]) for j in range(2)] for i in range(2)]};  det H = {sp.factor(Hs.det())}")
    # gradient-quadratic coefficients and Lambda^2 coefficients, to leading order (Neumann series of M^-1)
    K_ = t_ * M1 + u_ * M13
    Sig_lo = -S / kap + K_ * S / kap ** 2
    Q, Qp, W = terms(Sig_lo, {'t': M1 * S / kap ** 2, 'u': M13 * S / kap ** 2}, L, B)
    eps = sp.symbols('epsilon')

    def lead(expr):              # coefficient of kappa^-2 (all these start at that order)
        e = sp.expand(sp.expand(expr).subs(kap, 1 / eps))
        assert all(sp.Poly(e, eps).coeff_monomial(eps ** j) == 0 for j in (0, 1))
        return sp.factor(sp.Poly(e, eps).coeff_monomial(eps ** 2) / kap ** 2)
    P_("  gradient-quadratic terms (leading order, O(1/kappa^2); exact ones are used in J3):")
    for nm in ('E1', 'E13'):
        P_(f"    {nm}: Q_tt = {lead(Q[(nm, 't', 't')])}, Q_uu = {lead(Q[(nm, 'u', 'u')])}, "
           f"(Q_tu + Q_ut)/2 = {lead((Q[(nm, 't', 'u')] + Q[(nm, 'u', 't')]) / 2)};  Q'_tt = {lead(Qp[(nm, 't', 't')])}, "
           f"Q'_uu = {lead(Qp[(nm, 'u', 'u')])}, (Q'_tu + Q'_ut)/2 = {lead((Qp[(nm, 't', 'u')] + Qp[(nm, 'u', 't')]) / 2)}")
    P_("  Lambda^2-type equations  varpi dt ^ du = 0  (antisymmetric part of E1, E13; w3; w11), leading order:")
    for nm in CH:
        P_(f"    {nm}: varpi = {lead(W[nm])};  at Shiab#: {lead(W[nm].subs({xs: 1, zs: 1, ys: 0, ws: 0}))};"
           f"  at W_A: {lead(W[nm].subs({xs: 1, zs: 0, ys: 0, ws: 0}))}")
    sh = {xs: 1, zs: 1, ys: 0, ws: 0}
    Ssh = sp.simplify(Sig.subs(sh))
    _, _, Wsh = terms(Ssh, {'t': sp.diff(Ssh, t_), 'u': sp.diff(Ssh, u_)}, L, B)
    w13 = sp.factor(sp.simplify(Wsh['E13'].subs(sh)))
    P_(f"  Shiab#, exact: varpi(E13) = {w13}  (nonzero at T_o = 0, so dt ^ du = 0 wherever Upsilon = 0 holds on an open set"
       f" and this does not vanish)")
    # exact det H for every member (review R6-2): l_E vanishes on the A, B columns and S on the A, B rows, so
    # H = -l_CD (M^-1)_CD S_CD; Jacobi's complementary-minor identity det((M^-1)_CD) = det(M_AB)/det M = kappa^2/det M
    sgn = eta[7]
    detM = d1['detM']
    CDi, ABi = [1, 2], [0, 3]
    offb = all(lE[i, j] == 0 for i in range(2) for j in ABi) and all(S[j, k] == 0 for j in ABi for k in range(2))
    jac = sp.expand(M.adjugate().extract(CDi, CDi).det() - kap ** 2 * detM) == 0
    prodLS = sp.factor(lE.extract([0, 1], CDi).det() * S.extract(CDi, [0, 1]).det())
    okp = sp.expand(prodLS - 576 ** 2 * (xs * zs + sgn * ys * ws) ** 2) == 0
    P_(f"  exact, every member: l_E = 0 on the A, B columns and S = 0 on the A, B rows: {offb}; det adj(M)_CD = kappa^2 det M"
       f" (Jacobi): {jac}; det l_CD det S_CD = {prodLS};  hence det H = 576^2 (xz {'+' if sgn == 1 else '-'} yw)^2 kappa^2 / det M"
       f" exactly: {offb and jac and okp}")
    assert offb and jac and okp
    colC_ok = sp.expand(lE[0, 1] * ys - lE[1, 1] * xs) == 0
    detD = sp.factor(sp.Matrix([[lE[0, 2], xs], [lE[1, 2], ys]]).det())
    P_(f"  exact image condition: column C of l_E is proportional to (x, y): {colC_ok}; det[column D, (x, y)] = {detD};"
       f" so im H(t, u) is contained in span(x, y) whenever xz {'+' if sgn == 1 else '-'} yw = 0, for every (t, u)")
    assert colC_ok and sp.expand(detD / (xs * zs + sgn * ys * ws)).free_symbols <= set()
    swa = {xs: 1, zs: 0, ys: 0, ws: 0}
    Swa = sp.simplify(Sig.subs(swa))
    _, _, Wwa = terms(Swa, {'t': sp.diff(Swa, t_), 'u': sp.diff(Swa, u_)}, L, B)
    Wwa = {nm: sp.factor(sp.simplify(sp.sympify(v).subs(swa))) for nm, v in Wwa.items()}
    ref = (1 if sgn == 1 else -1) * 1916928 * kap * u_ / (kap ** 2 + 224 * kap * t_ - sgn * 3328 * u_ ** 2) ** 2
    okw = sp.simplify(Wwa['E1'] - ref) == 0 and all(Wwa[nm] == 0 for nm in ('E13', 'w3', 'w11'))
    P_(f"  W_A, exact varpi: E1 = {Wwa['E1']}; E13 = {Wwa['E13']}; w3 = {Wwa['w3']}; w11 = {Wwa['w11']}  (O(kappa^-3), nonzero"
       f" for u != 0: the antisymmetric part of E^(1) is varpi dt ^ du pointwise, so Upsilon = 0 at a point with u != 0"
       f" forces dt ^ du = 0 there): {okw}")
    assert okw
    return dict(L=L, B=B, H=H, H0=H0, lE=lE)


def terms(Sig, dSig, L, B):
    """gradient-quadratic coefficients Q (of dphi (x) dpsi), Q' (of <dphi, dpsi> g) in E1, E13 and the Lambda^2
    coefficients varpi (of dt ^ du) in E1_anti, E13_anti, w3, w11, for A_k = Sigma_k,t dt + Sigma_k,u du"""
    ph = ('t', 'u')
    Q, Qp, W = {}, {}, {}
    for nm in ('E1', 'E13'):
        for a, fa in enumerate(ph):
            for b, fb in enumerate(ph):
                v = sum(L[(nm, k)] * dSig[fb][ik, a] for ik, k in enumerate(KS))
                vp = -v
                for (ik, k), (il, l) in itertools.combinations_with_replacement(list(enumerate(KS)), 2):
                    b1, b2, b3 = B[(nm, k, l)]
                    wkl = sp.Rational(1, 2) if k == l else 1
                    v += wkl * (b1 * Sig[ik, a] * Sig[il, b] + b2 * Sig[il, a] * Sig[ik, b])
                    vp += wkl * b3 * Sig[ik, a] * Sig[il, b]
                Q[(nm, fa, fb)], Qp[(nm, fa, fb)] = v, vp
        W[nm] = (Q[(nm, 't', 'u')] - Q[(nm, 'u', 't')]) / 2
    for nm in ('w3', 'w11'):
        v = sum(L[(nm, k)] * (dSig['t'][ik, 1] - dSig['u'][ik, 0]) for ik, k in enumerate(KS))
        for (ik, k), (il, l) in itertools.combinations(list(enumerate(KS)), 2):
            v += B[(nm, k, l)] * (Sig[ik, 0] * Sig[il, 1] - Sig[ik, 1] * Sig[il, 0])
        W[nm] = v
    return Q, Qp, W


# ====================================================================== J3
def j3(C, d1, d2):
    """random 2-jet, numeric member; everything scaled to Gaussian integers (exact, Python ints)"""
    Al, eta = C.Al, C.eta
    P_(f"\n=== J3 {SIG[eta[7]]}. End-to-end check on a random 2-jet (exact; scaled Gaussian integers) ===")
    t1 = time.time()
    rng = random.Random(20261030 + eta[7])
    num = {xs: 2, ys: -1, zs: 3, ws: 1, kap: 5, t_: sp.Rational(1, 7), u_: sp.Rational(-1, 3)}
    do = 21                                                    # t = 3/do, u = -7/do
    Sig = d1['Sig']
    SigN = Sig.subs(num).applyfunc(sp.expand_complex)
    dSigN = {f: sp.diff(Sig, v).subs(num).applyfunc(sp.expand_complex) for f, v in (('t', t_), ('u', u_))}

    def den(Mx):
        return sp.ilcm(*[sp.fraction(c)[1] for e in Mx for c in (sp.re(e), sp.im(e))])
    D1 = int(den(SigN))
    D2 = int(sp.ilcm(D1, den(dSigN['t']), den(dSigN['u'])))
    dphi = {f: [rng.randint(-2, 2) for _ in range(N)] for f in ('t', 'u')}
    Hes = {}
    for f in ('t', 'u'):
        Hm = [[0] * N for _ in range(N)]
        for a in range(N):
            for b in range(a, N):
                if rng.random() < 0.3:
                    Hm[a][b] = Hm[b][a] = rng.randint(-2, 2)
        Hes[f] = Hm
    At = [[gi(D1 * (SigN[ik, 0] * dphi['t'][c] + SigN[ik, 1] * dphi['u'][c])) for c in range(N)] for ik in range(4)]
    Te = [dict() for _ in range(N)]                            # D1 * T_e
    for ik, k in enumerate(KS):
        Pk = C.P(k, At[ik])
        for a in range(N):
            I.acc_into(Te[a], Pk[a])
    Nt = []                                                    # D2 * nabla A_k
    for ik in range(4):
        Nm = [[None] * N for _ in range(N)]
        for c in range(N):
            for b in range(N):
                v = 0
                for j, f in enumerate(('t', 'u')):
                    v += SigN[ik, j] * Hes[f][c][b]
                    for g_ in ('t', 'u'):
                        v += dSigN[g_][ik, j] * dphi[g_][c] * dphi[f][b]
                Nm[c][b] = gi(D2 * v)
        Nt.append(Nm)
    cf = {op: gi(UNIT[op] * num[FAM[op]]) for op in OPN}

    def fam_out(pt):
        po = per_op(pt)
        out = [dict() for _ in range(N)]
        for op in OPN:
            for a in range(N):
                I.acc_into(out[a], {m: cf[op] * v * eta[a] for m, v in po[op][a].items()})
        return out
    To = [I.add(I.scale(C.Phi1[a], 3), I.scale(C.Phi1G[a], -7)) for a in range(N)]      # do * T_o
    z1, z2 = C.vwedge(dphi['t'], C.Phi1), C.vwedge(dphi['u'], C.Phi1G)
    zs_ = {k: I.add(z1.get(k, {}), z2.get(k, {})) for k in set(z1) | set(z2)}
    Ysrc = fam_out(Al.psi_form({k: v for k, v in zs_.items() if v}))
    Yb = fam_out(q2B(Al, Te, To))                               # 3 * 2B(D1 T_e, do T_o)
    res_ev = [I.add(I.add(I.scale(Ysrc[a], 3 * D1 * do), I.scale(Te[a], 3 * do * num[kap])), Yb[a]) for a in range(N)]
    ok_ev = all(r == {} for r in res_ev)
    P_(f"  numeric member (x, y, z, w, kappa, t, u) = (2, -1, 3, 1, 5, 1/7, -1/3), random dt, du, Hess t, Hess u (seed"
       f" {20261030 + eta[7]}): T_e = Sigma (dt, du) solves Y_ev = Shiab(d_0 T_o) + kappa T_e + 2B(T_e, T_o) = 0 exactly:"
       f" {ok_ev}  ({sum(len(x) for x in Te)} monomials in T_e)")
    assert ok_ev
    # odd part, derivative and quadratic pieces separately (they scale differently)
    zt = {}
    for ik, k in enumerate(KS):
        for key, v in C.dzeta(k, Nt[ik]).items():
            zt[key] = I.add(zt.get(key, {}), v)
    Yd = fam_out(Al.psi_form({k: v for k, v in zt.items() if v}))      # D2 * Y[Shiab(d_0 T_e)]
    Yq = fam_out(Al.quad3(Te))                                          # 3 D1^2 * Y[Q(T_e)]
    odd_only = all(pc(m) % 2 == 1 for Y in (Yd, Yq) for a in range(N) for m in Y[a])
    L, B = d2['L'], d2['B']
    L0 = {k_: 0 for k_ in L}
    B0 = {k_: ((0, 0, 0) if k_[0] in ('E1', 'E13') else 0) for k_ in B}
    HN = (d2['lE'] * SigN).subs(num).applyfunc(sp.expand_complex)
    box = {f: sum(eta[a] * Hes[f][a][a] for a in range(N)) for f in ('t', 'u')}
    ip = {(f, g_): sum(eta[a] * dphi[f][a] * dphi[g_][a] for a in range(N)) for f in ('t', 'u') for g_ in ('t', 'u')}
    oks = []
    for lab, (LL, BB), scale_, Ydir, withH in (('Shiab(d_0 T_e)', (L, B0), D2, Yd, True),
                                                ('Q(T_e)', (L0, B), 3 * D1 ** 2, Yq, False)):
        Qs, Qps, Ws = terms(SigN, dSigN, LL, BB)
        Qs, Qps, Ws = ({k_: sp.expand_complex(sp.sympify(v).subs(num)) for k_, v in X.items()} for X in (Qs, Qps, Ws))
        pred = {nm: [[0] * N for _ in range(N)] for nm in CH}
        for i, nm in enumerate(('E1', 'E13')):
            for a in range(N):
                for b in range(N):
                    v = 0
                    for j, f in enumerate(('t', 'u')):
                        if withH:
                            v += HN[i, j] * (Hes[f][a][b] - (box[f] * eta[a] if a == b else 0))
                        for g_ in ('t', 'u'):
                            v += Qs[(nm, f, g_)] * dphi[f][a] * dphi[g_][b]
                            if a == b:
                                v += Qps[(nm, f, g_)] * ip[(f, g_)] * eta[a]
                    pred[nm][a][b] = gi(scale_ * v)
        for nm in ('w3', 'w11'):
            for a in range(N):
                for b in range(N):
                    pred[nm][a][b] = gi(scale_ * Ws[nm] * (dphi['t'][a] * dphi['u'][b] - dphi['t'][b] * dphi['u'][a]))
        Yp = C.build(pred['E1'], pred['E13'], pred['w3'], pred['w11'])
        oks.append(all(I.add(Yp[a], Ydir[a], -1) == {} for a in range(N)))
    P_(f"  direct Y_odd[Shiab(d_0 T_e)] and Y_odd[Q(T_e)] are g_odd-valued: {odd_only}; equal to the assembled formula"
       f" (H and the derivative parts of Q, Q', varpi; the quadratic parts of Q, Q', varpi) exactly: {oks}"
       f"  ({time.time() - t1:.0f} s)")
    assert odd_only and all(oks)


# ====================================================================== J4
def j4(res):
    P_("\n=== J4. Reduction along s(X), n = 4 (closed forms of F2-F6; exact sympy) ===")
    lam, cH, cV, Lg = sp.symbols('lambda c_H c_V Lambda_g')
    n, Nf = 4, 10
    mu = lam / (1 - n * lam)
    trG0 = (Nf + n * mu) / cV                  # tr_{g_aleph} G_0 on V
    trTT = n / (cV * (1 - n * lam))            # tr_{g_aleph} tau(x)tau = |tau|^2
    geo = {}
    for Z in ('Z1', 'Z2'):
        if Z == 'Z1':
            e1 = cH / (4 * cV * (n * lam - 1))
            sig = -sp.Rational(23, 2) - (1 + lam) / (1 - 4 * lam)
            ric_V = (-sp.Rational(5, 4), sp.Rational(1, 4))          # -(n+1)/4 G_{1/(n+1)}
        else:
            e1 = 0
            sig = -sp.Integer(9)
            ric_V = (-sp.Integer(1), sp.Rational(1, 4))              # -(n/4) G_{1/n}
        sY = 4 * Lg / cH + sig / cV
        # G^Y in the basis (h = g^Y_HH = c_H g, G_0, tau tau) for an Einstein base Ric(g) = Lambda_g g
        G = (sp.simplify((Lg + e1) / cH - sY / 2), sp.simplify(ric_V[0] - sY / 2 * cV), sp.simplify(ric_V[1] + sY / 2 * cV * lam))
        LP = sp.simplify(e1 - cH * sig / (2 * cV))
        geo[Z] = dict(G=G, LP=LP, beta=(1 / (2 * cV * (1 - n * lam)) if Z == 'Z1' else 0))
        P_(f"  {Z}: on s(X), Ric(g) = Lambda_g g:  G^Y = ({G[0]}) h + ({G[1]}) G_0 + ({G[2]}) tau(x)tau  (h = c_H g on H),"
           f"  Lambda_P = {sp.factor(LP)};  Hess of a function with dt = a tau, t const on s(X): HH = "
           f"({geo[Z]['beta']}) a h")
    # ---- J4a/b: basic t, u
    P_("  (B) basic t = pi^* tau, u = pi^* upsilon.  On s(X): Hess_HH = Hess_g, Hess_VV = 0 (fibres totally geodesic, F3/F4(a)),"
       " Hess(X, U) = 0 (Z2: A = 0) or g_aleph(II(X, grad phi), U) (Z1, F4(a)); dphi horizontal.")
    P_("    (P) keeps the HH and HV blocks: Z2 -> the 4D system  -x_i (G(g) + Lambda_P g) + c_H c_i g + sum H_i,phi (nabla^2 tau_phi"
       " - Box tau_phi g) + sum Q dtau dtau + sum Q' <dtau,dtau> g = 0 (i = 1, 13) and varpi dtau ^ dupsilon = 0 (coefficients"
       " of J2, unchanged); Z1 -> in addition HV:  II(X, H_i,t grad tau + H_i,u grad upsilon) = 0, i.e. H (dtau, dupsilon)^T = 0.")
    P_("    (E) adds the VV block, -x_i G^Y_VV + (scalar) g_VV = 0:  needs Ric^Y_VV proportional to g_VV, i.e. Z1 with"
       " lambda = 1/5 (F6) when (x, y) != 0; and HV as in (P)-Z1.")
    # ---- J4c: isotropy-invariant 2-jet ansatz
    P_("  (J) isotropy-invariant 2-jet along s(X): t, u constant on s(X), dphi = a_phi tau (a const),"
       " Hess phi = beta_phi h + p_phi G_0 + q_phi tau(x)tau with beta_phi = beta a_phi; G(g) traceless part must vanish"
       " (g Einstein) when (x, y) != 0.")
    out = {}
    for sgn, d in res.items():
        L, B = d['L'], d['B']
        c1, c13 = d['c']
        # generic symbols for the slaving matrix Sigma(t, u) and its derivatives: every identity below is then an
        # identity for any Sigma; the actual Sigma of J1 is substituted only in the sample solutions
        Sig = sp.Matrix(4, 2, lambda i, j: sp.Symbol(f"S{KS[i]}_{'tu'[j]}"))
        dS = {f: sp.Matrix(4, 2, lambda i, j: sp.Symbol(f"dS{f}{KS[i]}_{'tu'[j]}")) for f in ('t', 'u')}
        H = d['lE'] * Sig
        out[(sgn, 'sym')] = (Sig, dS)
        for Z in ('Z1', 'Z2'):
            G, LP, beta = geo[Z]['G'], geo[Z]['LP'], geo[Z]['beta']
            at, au, pt, pu, qt, qu = sp.symbols('a_t a_u p_t p_u q_t q_u')
            a = {'t': at, 'u': au}
            p = {'t': pt, 'u': pu}
            q = {'t': qt, 'u': qu}
            al = [Sig[ik, 0] * at + Sig[ik, 1] * au for ik in range(4)]
            comps = {}
            for i, (nm, xi, ci) in enumerate((('E1', xs, c1), ('E13', ys, c13))):
                Nk = []
                for ik in range(4):
                    hcoef = sum(Sig[ik, j] * beta * a[f] for j, f in enumerate(('t', 'u')))
                    g0 = sum(Sig[ik, j] * p[f] for j, f in enumerate(('t', 'u')))
                    tt = sum(Sig[ik, j] * q[f] for j, f in enumerate(('t', 'u'))) + \
                        sum(dS[g_][ik, j] * a[g_] * a[f] for j, f in enumerate(('t', 'u')) for g_ in ('t', 'u'))
                    tr = 4 * hcoef + trG0 * g0 + trTT * tt
                    Nk.append((hcoef, g0, tt, tr))
                gal = (1, cV, -lam * cV)              # g_aleph in the basis (h, G_0, tau tau)
                e = [-xi * G[j] + ci * gal[j] for j in range(3)]
                for ik, k in enumerate(KS):
                    lk = L[(nm, k)]
                    hcoef, g0, tt, tr = Nk[ik]
                    e[0] += lk * (hcoef - tr * gal[0])
                    e[1] += lk * (g0 - tr * gal[1])
                    e[2] += lk * (tt - tr * gal[2])
                for (ik, k), (il, l) in itertools.combinations_with_replacement(list(enumerate(KS)), 2):
                    b1, b2, b3 = B[(nm, k, l)]
                    wkl = sp.Rational(1, 2) if k == l else 1
                    e[2] += wkl * (b1 + b2) * al[ik] * al[il]
                    for j in range(3):
                        e[j] += wkl * b3 * al[ik] * al[il] * trTT * gal[j]
                comps[nm] = e
            out[(sgn, Z)] = comps
            # balance identity: [h] - [G_0]/c_V = -x_i (G_h - G_G0/c_V) + sum H_i,phi (beta_phi - p_phi/c_V)
            for i, nm in enumerate(('E1', 'E13')):
                xi = (xs, ys)[i]
                lhs = comps[nm][0] - comps[nm][1] / cV
                rhs = -xi * (G[0] - G[1] / cV) + sum(H[i, j] * (beta * a[f] - p[f] / cV) for j, f in enumerate(('t', 'u')))
                assert sp.simplify(sp.expand(lhs - rhs)) == 0
            rho = sp.simplify(G[0] - G[1] / cV - Lg / cH)
            P_(f"  {SIG[sgn]} {Z}: balance identity [h] - [G_0]/c_V of E^(i):  x_i (Lambda_g/c_H + rho) = sum_phi H_i,phi"
               f" (beta a_phi - p_phi/c_V) with rho = {sp.factor(rho)}  (polynomial identity for generic Sigma, i = 1, 13)")
        # a = 0 branch: closed form (generic Sigma)
        for Z in ('Z1', 'Z2'):
            comps = out[(sgn, Z)]
            sub0 = {sp.Symbol('a_t'): 0, sp.Symbol('a_u'): 0}
            G = geo[Z]['G']
            Th = sp.simplify(trG0 * G[1] + trTT * G[2] - 9 * G[0] + 0 * Lg)   # tr_V G^Y_VV - 9 (G^Y)_h
            # claim: on a = 0 the six equations hold iff c_i = x_i * Theta and (P_i, Q_i) as below
            ok = True
            for i, (nm, xi, ci) in enumerate((('E1', xs, c1), ('E13', ys, c13))):
                e = [sp.expand(comps[nm][j].subs(sub0)) for j in range(3)]
                Pi = sum(H[i, j] * sp.Symbol('p_' + f) for j, f in enumerate(('t', 'u')))
                Qi = sum(H[i, j] * sp.Symbol('q_' + f) for j, f in enumerate(('t', 'u')))
                Bi = trG0 * Pi + trTT * Qi
                ok &= sp.simplify(e[0] - (-xi * G[0] + ci - Bi)) == 0
                ok &= sp.simplify(e[1] - (-xi * G[1] + cV * (ci - Bi) + Pi)) == 0
                ok &= sp.simplify(e[2] - (-xi * G[2] - lam * cV * (ci - Bi) + Qi)) == 0
            Th0 = sp.simplify(Th.subs(Lg, 0))
            slope = sp.simplify(sp.diff(Th, Lg))
            P_(f"  {SIG[sgn]} {Z}, branch a = 0 (T_e = 0 on s(X), nabla T_e != 0):  E^(i) = [-x_i G_h + c_i - B_i] h +"
               f" [-x_i G_G0 + c_V (c_i - B_i) + P_i] G_0 + [-x_i G_tt - lambda c_V (c_i - B_i) + Q_i] tau tau with"
               f" P_i = sum H_i,phi p_phi, Q_i = sum H_i,phi q_phi, B_i = tr: {ok};  eliminating P_i, Q_i, B_i:"
               f"  c_i(t, u) = x_i Theta,  Theta = tr_V G^Y_VV - 9 G^Y_h = {Th0} + ({slope}) Lambda_g")
            assert ok
            out[(sgn, Z, 'Theta')] = (Th0, slope)
    return out, geo


def j4_samples(res, out, geo):
    """exact sample solutions of (E) at s(X) in ansatz (J), a = 0 branch, Shiab# and W_A, both signatures, Z1 and Z2.
    (14,0): u(128)-real form x = -i x', z = -i z', t = i s, u = i r (s, r, kappa real); then Lambda_g must be real and
    p_phi, q_phi (Hessian coefficients of t = i s) in i R.  (7,7): real code normalisation."""
    P_("\n  exact sample solutions of (E) at s(X) in ansatz (J), branch a = 0 (all six components, Einstein base):")
    lam, cH, cV, Lg = sp.symbols('lambda c_H c_V Lambda_g')
    pt, pu, qt, qu = sp.symbols('p_t p_u q_t q_u')
    for sgn, d in res.items():
        c1, c13 = d['c']
        unit = -sp.I if sgn == 1 else 1
        tun = sp.I if sgn == 1 else 1
        lv = 0 if sgn == 1 else sp.Rational(1, 2)
        for opn, sub in (('Shiab#', {xs: unit, zs: unit, ys: 0, ws: 0}), ('W_A', {xs: unit, zs: 0, ys: 0, ws: 0})):
            for Z in ('Z1', 'Z2'):
                for tv in (sp.Rational(1, 1000), sp.Rational(-1, 1000)):
                    num = dict(sub)
                    num.update({kap: 1, cH: 1, cV: 1, lam: lv, u_: 0, t_: tun * tv})
                    Th0, slope = out[(sgn, Z, 'Theta')]
                    c1n, c13n = sp.expand_complex(c1.subs(num)), sp.expand_complex(c13.subs(num))
                    assert c13n == 0
                    Lgv = sp.expand_complex((c1n / sub[xs] - Th0.subs(num)) / slope.subs(num))
                    Sg, dSg = out[(sgn, 'sym')]
                    SigN = d['Sig'].subs(num).applyfunc(sp.expand_complex)
                    rep = {Sg[i, j]: SigN[i, j] for i in range(4) for j in range(2)}
                    rep.update({dSg[f][i, j]: 0 for f in ('t', 'u') for i in range(4) for j in range(2)})  # multiply a = 0
                    eqs = [sp.expand(sp.together(out[(sgn, Z)][nm][j].subs({sp.Symbol('a_t'): 0, sp.Symbol('a_u'): 0})
                                                  .subs(rep).subs(num).subs(Lg, Lgv))) for nm in ('E1', 'E13') for j in range(3)]
                    sol = sp.solve(eqs, [pt, pu, qt, qu], dict=True)
                    sol = [{k_: sp.expand_complex(v) for k_, v in sol[0].items()}] if sol else sol
                    ok = bool(sol) and all(sp.simplify(e.subs(sol[0])) == 0 for e in eqs)
                    real = sp.im(Lgv) == 0 and all(sp.im(sp.expand_complex(v / tun)) == 0 for v in sol[0].values())
                    Hn = d['H'].subs(num).applyfunc(sp.expand_complex)
                    P_(f"   {SIG[sgn]} {opn:6s} {Z} lambda = {lv}, c_H = c_V = kappa = 1, u = 0, "
                       f"{'t = i*' if sgn == 1 else 't = '}{tv}:  Lambda_g = {Lgv};  det H = {sp.expand_complex(Hn.det())};"
                       f"  {sol[0]};  all six = 0: {ok};  real in the real form: {real}")
                    assert ok and real


def j4b(res, out, geo):
    """ansatz (J) with a normal gradient a != 0 (review R6-1): general-a closed form J4(c'), the (14,0) bounds for
    Shiab# in Z1 and Z2, the referee's Z1 configuration, and an exact solution attaining the Z1 bound"""
    P_("\n=== J4b. Ansatz (J) with a != 0: closed form for every a, (14,0) bounds, exact Z1 configurations ===")
    lam, cH, cV, Lg = sp.symbols('lambda c_H c_V Lambda_g')
    at, au, pt, pu, qt, qu = sp.symbols('a_t a_u p_t p_u q_t q_u')
    n, Nf = 4, 10
    mu = lam / (1 - n * lam)
    trG0 = (Nf + n * mu) / cV
    trTT = n / (cV * (1 - n * lam))
    for sgn, d in res.items():
        Sg, dSg = out[(sgn, 'sym')]
        L, B = d['L'], d['B']
        c1, c13 = d['c']
        H = d['lE'] * Sg
        al = [Sg[ik, 0] * at + Sg[ik, 1] * au for ik in range(4)]
        Ha = H * sp.Matrix([at, au])
        for Z in ('Z1', 'Z2'):
            G, beta = geo[Z]['G'], geo[Z]['beta']
            Th = trG0 * G[1] + trTT * G[2] - 9 * G[0]
            ok = True
            for i, (nm, xi, ci) in enumerate((('E1', xs, c1), ('E13', ys, c13))):
                e = out[(sgn, Z)][nm]
                Bi = 0
                for (ik, k), (il, l) in itertools.combinations_with_replacement(list(enumerate(KS)), 2):
                    b1, b2, b3 = B[(nm, k, l)]
                    Bi += (sp.Rational(1, 2) if k == l else 1) * (b1 + b2 + b3) * al[ik] * al[il]
                lhs = ci - xi * Th - 13 * beta * Ha[i] + trTT * Bi
                rhs = -9 * e[0] + trG0 * e[1] + trTT * e[2]
                ok &= sp.expand(sp.together(lhs - rhs)) == 0
                # (P_i, Q_i) enter e_G0, e_tt with the scalar matrix [[1 - c_V trG0, -c_V trTT], [lam c_V trG0, 1 + lam c_V trTT]]
                for j, f in enumerate(('t', 'u')):
                    pf, qf = sp.Symbol('p_' + f), sp.Symbol('q_' + f)
                    ok &= sp.simplify(sp.diff(e[1], pf) - (1 - cV * trG0) * H[i, j]) == 0
                    ok &= sp.simplify(sp.diff(e[1], qf) + cV * trTT * H[i, j]) == 0
                    ok &= sp.simplify(sp.diff(e[2], pf) - lam * cV * trG0 * H[i, j]) == 0
                    ok &= sp.simplify(sp.diff(e[2], qf) - (1 + lam * cV * trTT) * H[i, j]) == 0
            dPQ = sp.simplify((1 - cV * trG0) * (1 + lam * cV * trTT) + lam * cV ** 2 * trG0 * trTT)
            P_(f"  {SIG[sgn]} {Z}: -9 e_h + tr G_0 e_G0 + |tau|^2 e_tt = c_i - x_i Theta - 13 beta (H a)_i + |tau|^2 B_i(alpha),"
               f" B_i = sum_(k<=l) w_kl (b1+b2+b3)^(i)_kl alpha_k alpha_l, alpha_k = sum_phi Sigma_k,phi a_phi (identity for generic"
               f" Sigma, dSigma; dSigma drops out); the other two equations fix (P_i, Q_i) through a scalar 2x2 matrix of"
               f" determinant {dPQ}: {ok}")
            assert ok and dPQ == -9
    # ---- (14,0) Shiab#: reality, completed square, bounds
    d = res[1]
    Sig, c1 = d['Sig'], d['c'][0]
    s_, r_, as_, ar_ = sp.symbols('s r a_s a_r', real=True)
    rf = {xs: -sp.I, zs: -sp.I, ys: 0, ws: 0}
    rng = random.Random(20261030)
    okr = True
    for _ in range(4):
        kv = rng.choice((1, 2, 3))
        num = dict(rf)
        num.update({kap: kv, t_: sp.I * sp.Rational(rng.randint(-9, 9) * kv, 1000),
                    u_: sp.I * sp.Rational(rng.randint(-9, 9) * kv, 1000)})
        SigN = Sig.subs(num).applyfunc(sp.expand_complex)
        aN = [sp.I * sp.Rational(rng.randint(-9, 9), 100), sp.I * sp.Rational(rng.randint(-9, 9), 100)]
        alN = [sp.expand_complex(SigN[ik, 0] * aN[0] + SigN[ik, 1] * aN[1]) for ik in range(4)]
        okr &= sp.re(alN[0]) == 0 and sp.im(alN[1]) == 0 and sp.re(alN[2]) == 0 and sp.im(alN[3]) == 0
    P_(f"  (14,0) u(128) real form, Shiab#: alpha_A in iR, alpha_C in R, alpha_D in iR, alpha_B in R at 4 random exact points"
       f" (T_e in u(128), as the slaving preserves the real form): {okr}")
    assert okr
    aC, aD = sp.symbols('alpha_C alpha_D')
    alg = {0: sp.Symbol('alpha_A'), 2: aC, 12: aD, 14: sp.Symbol('alpha_B')}
    B1 = 0
    for k, l in itertools.combinations_with_replacement(KS, 2):
        b1, b2, b3 = d['B'][('E1', k, l)]
        B1 += (sp.Rational(1, 2) if k == l else 1) * (b1 + b2 + b3) * alg[k] * alg[l]
    B1s = sp.expand(B1.subs({ys: 0, ws: 0}))
    Ha1 = sp.expand(sum(d['L'][('E1', k)] * alg[k] for k in KS).subs({ys: 0, ws: 0}))
    okB = sp.expand(B1s - xs * (-312 * aC ** 2 + 104 * aD ** 2)) == 0 and sp.expand(Ha1 - 24 * xs * aC) == 0
    e1 = sp.expand(-sp.I * c1.subs({xs: -sp.I, ys: 0, t_: sp.I * s_, u_: sp.I * r_}))
    oke = sp.expand(e1 - (kap * s_ + 312 * s_ ** 2 + 104 * r_ ** 2)) == 0
    bt, tt2, dlt = sp.symbols('beta tau2 delta', real=True)
    Thx = -e1 - 13 * bt * 24 * aC + tt2 * (-312 * aC ** 2 - 104 * dlt ** 2)          # alpha_D = i delta
    sq = -312 * tt2 * (aC + bt / (2 * tt2)) ** 2 + 78 * bt ** 2 / tt2 - 104 * tt2 * dlt ** 2 - e1
    oksq = sp.simplify(Thx - sq) == 0
    lam, cH, cV = sp.symbols('lambda c_H c_V')
    b_Z1, tau2 = 1 / (2 * cV * (1 - 4 * lam)), 4 / (cV * (1 - 4 * lam))
    oksh = sp.simplify(78 * b_Z1 ** 2 / tau2 - sp.Rational(39, 8) / (cV * (1 - 4 * lam))) == 0
    P_(f"  Shiab# (y = w = 0, any z): B_1 = x(-312 alpha_C^2 + 104 alpha_D^2) and (H a)_1 = 24 x alpha_C: {okB};"
       f"  e1 = kappa s + 312 s^2 + 104 r^2 >= -kappa^2/1248 for all real s, r: {oke}")
    P_(f"  hence, with c1/x = -e1 and alpha_D = i delta:  Theta = -e1 - 312 beta alpha_C - 312 |tau|^2 alpha_C^2 - 104 |tau|^2"
       f" delta^2 = -e1 + 78 beta^2/|tau|^2 - 312|tau|^2 (alpha_C + beta/(2|tau|^2))^2 - 104|tau|^2 delta^2: {oksq};"
       f"  78 beta^2/|tau|^2 = 39/(8 c_V (1 - 4 lambda)) in Z1 ({oksh}), 0 in Z2")
    assert okB and oke and oksq and oksh
    P_("  => for every a, every (t, u) and c_H > 0 (definite (14,0), lambda < 1/4):  Lambda_g = (c_H/11)(Theta_0 - Theta) >="
       " (c_H/11)(Theta_0 - kappa^2/1248) in Z2 and >= (c_H/11)(Theta_0 - kappa^2/1248 - 39/(8 c_V (1 - 4 lambda))) in Z1")
    # ---- exact Z1 configurations (lambda = 0, c_H = c_V = kappa = 1, real form): the referee's one, and one at the bound
    Lg = sp.symbols('Lambda_g')
    Sg, dSg = out[(1, 'sym')]
    base = dict(rf)
    base.update({kap: 1, cH: 1, cV: 1, lam: 0, u_: 0})
    Th0 = out[(1, 'Z1', 'Theta')][0].subs(base)
    bound = sp.Rational(1, 11) * (Th0 - sp.Rational(1, 1248) - sp.Rational(39, 8))
    a0bound = sp.Rational(1, 11) * (Th0 - sp.Rational(1, 1248))

    def check(tv, atv, given=None):
        num = dict(base)
        num[t_] = tv
        SigN = Sig.subs(num).applyfunc(sp.expand_complex)
        dSN = {f: sp.diff(Sig, v).subs(num).applyfunc(sp.expand_complex) for f, v in (('t', t_), ('u', u_))}
        rep = {Sg[i, j]: SigN[i, j] for i in range(4) for j in range(2)}
        rep.update({dSg[f][i, j]: dSN[f][i, j] for f in ('t', 'u') for i in range(4) for j in range(2)})
        rep.update({at: atv, au: 0, pu: 0, qu: 0})
        eqs = [sp.expand_complex(sp.together(out[(1, 'Z1')][nm][j].subs(rep).subs(num))) for nm in ('E1', 'E13')
               for j in range(3)]
        if given is None:
            sol = sp.solve(eqs, [pt, qt, Lg], dict=True)
            assert len(sol) == 1
            sol = {k_: sp.expand_complex(v) for k_, v in sol[0].items()}
        else:
            sol = given
        res_ = [sp.simplify(e.subs(sol)) for e in eqs]
        detMn = sp.expand_complex(res[1]['detM'].subs(num))
        real = sp.im(sol[Lg]) == 0 and sp.re(sol[pt]) == 0 and sp.re(sol[qt]) == 0 and sp.re(atv) == 0 and sp.re(tv) == 0
        return sol, all(r == 0 for r in res_) and detMn != 0, real, detMn
    given = {pt: sp.Rational(-24174663, 11000000000) * sp.I, qt: sp.Rational(17, 4000) * sp.I,
             Lg: sp.Rational(-984211, 1375000)}
    sol, ok6, real, dMn = check(sp.I / 1000, sp.Rational(-51, 16000) * sp.I, given)
    P_(f"  referee's configuration (R6-1): Z1, lambda = 0, c_H = c_V = kappa = 1, Shiab#, t = i/1000, u = 0, a_t = -51i/16000,"
       f" a_u = 0, p_t = -24174663i/11000000000, q_t = 17i/4000, p_u = q_u = 0, Lambda_g = -984211/1375000 (det M = {dMn}):  all six"
       f" components vanish: {ok6}; real in the real form: {real};  a = 0 lower bound {a0bound} = {float(a0bound):.6f},"
       f" Lambda_g = {float(given[Lg]):.6f} is below it: {given[Lg] < a0bound}")
    assert ok6 and real and given[Lg] < a0bound
    # attaining the Z1 bound: s = -kappa/624, alpha_C = -beta/(2|tau|^2) = -1/16 (lambda = 0, c_V = 1), alpha_C = 24 a_s/(kappa + 224 s)
    sv = sp.Rational(-1, 624)
    asv = sp.Rational(-1, 16) * (1 + 224 * sv) / 24
    sol, ok6, real, dMn2 = check(sp.I * sv, sp.I * asv)
    P_(f"  bound attained: t = i(-1/624), u = 0, a_t = i({asv}), a_u = 0:  unique (p_t, q_t, Lambda_g) = ({sol[pt]}, {sol[qt]},"
       f" {sol[Lg]}) (det M = {dMn2}); all six vanish: {ok6}; real: {real}; Lambda_g = (c_H/11)(Theta_0 - 1/1248 - 39/8) = {bound}:"
       f" {sp.simplify(sol[Lg] - bound) == 0}")
    assert ok6 and real and sp.simplify(sol[Lg] - bound) == 0
    P_("  (14,0) range of Lambda_g for Shiab# under (J): exactly [(c_H/11)(Theta_0 - kappa^2/1248), oo) in Z2 (every a) and exactly"
       " [(c_H/11)(Theta_0 - kappa^2/1248 - 39/(8 c_V (1 - 4 lambda))), oo) in Z1 (the slice u = a_u = 0 is connected,"
       " has det M != 0 for s >= -kappa/624, attains the minimum and is unbounded above)")
    rad = sp.solve(sp.Eq(-cH / cV, sp.Rational(1, 11) * cH * (-sp.Rational(9, 2) / cV - kap ** 2 / 1248)), kap ** 2)
    P_(f"  radial Z2 in (14,0) (J4(e): Lambda_g = -c_H/c_V) meets the Z2 bound iff kappa^2 >= {rad[0]}, i.e. c_V kappa^2 >= 8112")
    assert sp.simplify(rad[0] - 8112 / cV) == 0


def j5_leaf_test():
    """indefinite case of Theorem J5: the order-eps^2 leaf condition sum_i eps_i G_0(R(X,E_i).h, R(Y,E_i).h) = 0 for a base
    of constant curvature k != 0 (exact, sympy)"""
    P_("\n=== J5. Radial ansatz on a neighbourhood, Z2: order-eps^2 leaf condition for constant-curvature bases"
       " (bases (4,0), (1,3) and (2,2); the last two give (7,7) for lambda > 1/4, A5) ===")
    n = 4
    for sig in ((1, 1, 1, 1), (1, -1, -1, -1), (1, 1, -1, -1)):
        g = sp.diag(*sig)
        gi = g.inv()

        def Rend(a, b):   # R(e_a, e_b) Z = g(e_b, Z) e_a - g(e_a, Z) e_b  (k = 1)
            return sp.Matrix(n, n, lambda c, dd: (1 if c == a else 0) * g[b, dd] - (1 if c == b else 0) * g[a, dd])
        hs = sp.symbols('h0:10')
        h = sp.zeros(n, n)
        it = iter(hs)
        for i in range(n):
            for j in range(i, n):
                h[i, j] = h[j, i] = next(it)
        h = h.subs(hs[0], sp.solve((gi * h).trace(), hs[0])[0])          # g-traceless h (tangent to the leaf)
        act = lambda B, X: -(B.T * X + X * B)  # noqa: E731
        G0 = lambda A_, B_: (gi * A_ * gi * B_).trace()  # noqa: E731
        Fs = [sp.expand(sum(sig[i] * G0(act(Rend(a, i), h), act(Rend(b, i), h)) for i in range(n)))
              for a, b in ((0, 0), (0, 1), (1, 1))]
        P_(f"  base signature {sig}: F(e0,e0; h) = {sp.factor(Fs[0])};  nonzero for (e0,e0), (e0,e1), (e1,e1):"
           f" {[F != 0 for F in Fs]}")
        assert all(F != 0 for F in Fs)


def j4_radial(res, out, geo):
    P_("\n  radial sub-case of (J) (t = f(s), u = h(s), s = log|det(g^-1 y)|/4 near s(X): p_phi = 0); geometry only, both"
       " signatures:")
    lam, cH, cV, Lg = sp.symbols('lambda c_H c_V Lambda_g')
    for Z in ('Z1', 'Z2'):
        G = geo[Z]['G']
        rho = sp.simplify(G[0] - G[1] / cV - Lg / cH)
        assert sp.diff(rho, Lg) == 0
        if Z == 'Z2':
            assert sp.simplify(rho - 1 / cV) == 0
            P_(f"   Z2: rho = {rho}; beta = 0 and p = 0, so the balance identity gives x_i (Lambda_g/c_H + {rho}) = 0:"
               f"  Lambda_g = -c_H/c_V for every member with (x, y) != 0 and every (z, w)")
        else:
            assert sp.simplify(rho - (5 * lam - 1) / (cV * (4 * lam - 1))) == 0
            P_(f"   Z1: rho = {sp.factor(rho)}; x_i (Lambda_g/c_H + rho) = beta sum_phi H_i,phi a_phi, beta = {geo[Z]['beta']}:"
               f" the normal gradient a shifts Lambda_g; at a = 0, Lambda_g = {sp.factor(-cH * rho)}"
               f" (= 0 at lambda = 1/5, the F6 value: {sp.simplify((-cH * rho).subs(lam, sp.Rational(1, 5))) == 0})")


def main():
    if FAST:
        P_("*** J_FAST=1: debugging mode, orbit representatives only; not a verification run ***")
    res = {}
    for eta in (EUC, E77):
        C = Ctx(eta)
        c = j0(C)
        d1 = j1(C)
        d2 = j2(C, d1)
        j3(C, d1, d2)
        res[eta[7]] = dict(c=c, S=d1['S'], Sig=d1['Sig'], H=d2['H'], L=d2['L'], B=d2['B'], lE=d2['lE'], detM=d1['detM'])
        P_(f"  [{SIG[eta[7]]} done, {time.time() - T0:.0f} s]")
    out, geo = j4(res)
    j4_samples(res, out, geo)
    j4_radial(res, out, geo)
    j4b(res, out, geo)
    j5_leaf_test()
    P_(f"\nruntime {time.time() - T0:.0f} s")


if __name__ == '__main__':
    main()
