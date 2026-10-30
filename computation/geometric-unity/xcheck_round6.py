"""Referee cross-check of GU round 6 (formalization/geometric-unity/J-slaved-sector.md), independent construction.

Imports only xcheck_round5.py, unmodified: the round-5 referee's construction (own Clifford algebra, the four Shiab
templates built from their definitions, Q from p(alpha, Q(T)) = (1/3) Psi(alpha^T^T) + (1/3) Psi(T^[T^alpha])).
Nothing from J_slaved_sector.py, I_elliptic_reduction.py, G_symbols.py, H_ellipticity.py or shiab_checks.py.
New here:
  - the vector-sector patterns P_0, P_2, P_12, P_14 and an exact pattern fit;
  - a direct bilinear 2B (dual of 3*2B(X,Y) = Psi(a^(X^Y+Y^X)) + Psi(X^[Y^a]) + Psi(Y^[X^a])), checked against
    polarisation of the round-5 Q;
  - an own grade extraction of g_odd-valued 1-forms (grades 1, 13, 3, 11) with exact reconstruction;
  - own bookkeeping of the six equations of J4 from the F3/F4 inputs (with generic Sigma, l and b symbols), and their
    reduction for a nonzero normal gradient a (X4b): c_i = x_i Theta + 13 beta (H a)_i - |tau|^2 B_i(alpha), with an exact
    (14,0) Z1 solution below the half-line that J4(d) states as exact;
  - own constant-curvature leaf test for J5(d), including base signature (2,2).
Conventions (J.0, round 5): Upsilon = *Y, Y_a = eta_a D_a where Psi(alpha ^ zeta) = sum_a Tr(alpha_a D_a); the operator
is PHASE[op] times the real template; family Shiab = x W_A + y W'_A + z W_C.AdG + w W'_C.AdG; the kappa-term is Y = kappa T.
Exact arithmetic (ints, Fractions, sympy).  Seeds fixed.  Single core.  Every check prints True/False; a False is
recorded, the run continues.  Sampling (stated where used): the bilinear map of J2 is evaluated on sampled basis pairs
(orbit representatives of the signed-permutation symmetry plus extra pairs), not on all pairs.
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
import xcheck_round5 as R5  # noqa: E402

T0 = time.time()
random.seed(20261030)
N, FULL = R5.N, R5.FULL
OPS, PHASE = R5.OPS, R5.PHASE
I = sp.I
x, y, z, w, kap, t, u = sp.symbols('x y z w kappa t u')
FAM = dict(zip(OPS, (x, y, z, w)))
KS = (0, 2, 12, 14)
KN = {0: 'A', 2: 'C', 12: 'D', 14: 'B'}
CH = ('E1', 'E13', 'w3', 'w11')
pc = R5.pc
RES = []
SH = {x: 1, y: 0, z: 1, w: 0}          # Shiab# (code normalisation)
WA = {x: 1, y: 0, z: 0, w: 0}


def say(*a):
    print(*a)
    sys.stdout.flush()


def check(label, ok):
    ok = bool(ok)
    RES.append((label, ok))
    say(f"  [{ok}] {label}")
    return ok


def e(p):
    return [1 if i == p else 0 for i in range(N)]


def rat(v):
    v = Fr(v)
    return sp.Rational(v.numerator, v.denominator)


def fam(d):
    """{op: rational template coefficient} -> family coefficient (sympy)"""
    return sp.expand(sum(FAM[op] * PHASE[op] * rat(d.get(op, 0)) for op in OPS))


def zero(expr):
    return sp.simplify(sp.expand(expr)) == 0


def mzero(Mx):
    return all(zero(v) for v in Mx)


# ------------------------------------------------------------------------------------------------ algebra context
class Ctx:
    def __init__(self, eta):
        self.eta = eta
        self.S = R5.Shiab(eta)
        self.cl = self.S.cl
        self.G = {FULL: 1}
        self.Phi1 = [{1 << a: 1} for a in range(N)]
        self.Phi1G = [self.cl.mul({1 << a: 1}, self.G) for a in range(N)]
        self.eps = eta[7]                      # +1 in (14,0), -1 in (7,7)
        self.lab = '(14,0)' if eta[7] == 1 else '(7,7)'
        self._g3 = {}

    def P(self, k, v):
        """P_0(v)_a = v_a 1, P_14(v)_a = v_a Gamma, P_2(v)_a = sum_{b!=a} v_b gamma_a gamma^b, P_12 = P_2 Gamma"""
        cl, eta = self.cl, self.eta
        res = []
        for a in range(N):
            X = {}
            if k in (0, 14):
                if v[a]:
                    X = {0 if k == 0 else FULL: v[a]}
            else:
                for b in range(N):
                    if b != a and v[b]:
                        X = R5.add(X, R5.scl(cl.mul({1 << a: 1}, {1 << b: eta[b]}), v[b]))
                if k == 12 and X:
                    X = cl.mul(X, self.G)
            res.append(X)
        return res

    def Yform(self, D, div=1):
        return [{m: Fr(v * self.eta[a], div) for m, v in D[a].items()} for a in range(N)]

    def g3(self, a, b, c):
        key = (a, b, c)
        if key not in self._g3:
            X = self.cl.mul(self.cl.mul({1 << a: 1}, {1 << b: 1}), {1 << c: 1})
            (m, s), = X.items()
            self._g3[key] = (m, s)
        return self._g3[key]

    # ---- 2B: dual of 3 * 2B(X, Y) (template of op)
    def second(self, op, T1, T2):
        S, cl = self.S, self.cl
        D = [dict() for _ in range(N)]
        for b, c in itertools.combinations(range(N), 2):
            for a in range(N):
                if a in (b, c) or not T1[a]:
                    continue
                Pm = S.Mt(op, (a, b, c), T1[a])
                if not Pm:
                    continue
                if T2[b]:
                    R5.acc(D[c], R5.add(cl.mul(Pm, T2[b]), cl.mul(T2[b], Pm), -1))
                if T2[c]:
                    R5.acc(D[b], R5.add(cl.mul(Pm, T2[c]), cl.mul(T2[c], Pm), -1), -1)
        return D

    def B2(self, op, X, Y):
        z1, z2 = self.S.wedge(X, Y), self.S.wedge(Y, X)
        zz = {}
        for k in set(z1) | set(z2):
            v = R5.add(z1.get(k, {}), z2.get(k, {}))
            if v:
                zz[k] = v
        D = self.S.dual(op, zz)
        for Dx in (self.second(op, X, Y), self.second(op, Y, X)):
            for a in range(N):
                R5.acc(D[a], Dx[a])
        return D

    # ---- fit Y = sum_k s_k P_k(e^c)
    def fitV(self, Yf, c):
        coef, ok, supp = {}, True, set()
        for k in KS:
            pat = self.P(k, e(c))
            ratios = set()
            for a in range(N):
                for m, pv in pat[a].items():
                    supp.add((a, m))
                    ratios.add(Fr(Yf[a].get(m, 0)) / pv)
            ok &= len(ratios) == 1
            coef[k] = next(iter(ratios))
        ok &= all((a, m) in supp for a in range(N) for m in Yf[a])
        return coef, ok

    # ---- odd 1-form -> (E1, E13, W3, W11); exact reconstruction flag
    def decomp(self, Yf):
        eta, cl = self.eta, self.cl
        E1 = [[0] * N for _ in range(N)]
        E13 = [[0] * N for _ in range(N)]
        W3, W11 = {}, {}
        ok = True
        for a in range(N):
            for m, v in Yf[a].items():
                g = pc(m)
                if g == 1:
                    b = m.bit_length() - 1
                    E1[a][b] = v * eta[b]
                elif g == 13:
                    b = (FULL ^ m).bit_length() - 1
                    E13[a][b] = v * eta[b] * cl.sign(1 << b, FULL)
                elif g == 3 and (m >> a) & 1:
                    b, c = [i for i in range(N) if (m >> i) & 1 and i != a]
                    _, s = self.g3(a, b, c)
                    W3.setdefault((b, c), Fr(v) / (2 * eta[b] * eta[c] * s))
                elif g == 11 and (FULL ^ m) >> a & 1:
                    tri = FULL ^ m
                    b, c = [i for i in range(N) if (tri >> i) & 1 and i != a]
                    mm, s = self.g3(a, b, c)
                    W11.setdefault((b, c), Fr(v) / (2 * eta[b] * eta[c] * s * cl.sign(mm, FULL)))
                else:
                    ok = False
        R = [dict() for _ in range(N)]
        for a in range(N):
            for b in range(N):
                if E1[a][b]:
                    R5.acc(R[a], {1 << b: E1[a][b] * eta[b]})
                if E13[a][b]:
                    R5.acc(R[a], {FULL ^ (1 << b): E13[a][b] * eta[b] * cl.sign(1 << b, FULL)})
            for (b, c), val in W3.items():
                if a not in (b, c):
                    m, s = self.g3(a, b, c)
                    R5.acc(R[a], {m: 2 * val * eta[b] * eta[c] * s})
            for (b, c), val in W11.items():
                if a not in (b, c):
                    m, s = self.g3(a, b, c)
                    R5.acc(R[a], {FULL ^ m: 2 * val * eta[b] * eta[c] * s * cl.sign(m, FULL)})
        ok &= all(R5.add(R[a], Yf[a], -1) == {} for a in range(N))
        return (E1, E13, W3, W11), ok


def vwedge(v, P1):
    zt = {}
    for b, c in itertools.combinations(range(N), 2):
        X = {}
        if v[b] and P1[c]:
            X = R5.add(X, R5.scl(P1[c], v[b]))
        if v[c] and P1[b]:
            X = R5.add(X, R5.scl(P1[b], -v[c]))
        if X:
            zt[(b, c)] = X
    return zt


# ================================================================================================ X0: base
def x0(C):
    say(f"\n######## {C.lab} ########")
    say("X0  base normalisation: c1, c13 of I6 from the round-5 Q (Y_a = c1 gamma_a + c13 gamma_a Gamma)")
    c1, c13 = kap * t, kap * u
    sg = C.cl.sign(1, FULL)
    for op in OPS:
        vals = {}
        for tv, uv in ((1, 0), (0, 1), (1, 1)):
            T = [R5.add(R5.scl(C.Phi1[a], tv), R5.scl(C.Phi1G[a], uv)) for a in range(N)]
            Yq = C.Yform(C.S.Q3(op, T), 3)
            vals[(tv, uv)] = (Yq[0].get(1, 0), Yq[0].get(FULL ^ 1, 0) * sg)
        for i in range(2):
            A_, C_ = rat(vals[(1, 0)][i]), rat(vals[(0, 1)][i])
            B_ = rat(vals[(1, 1)][i]) - A_ - C_
            q = FAM[op] * PHASE[op] * (A_ * t ** 2 + B_ * t * u + C_ * u ** 2)
            if i == 0:
                c1 += q
            else:
                c13 += q
    c1, c13 = sp.expand(c1), sp.expand(c13)
    ep = C.eps
    ref1 = kap * t + x * (312 * t ** 2 + ep * 104 * u ** 2) + ep * 208 * y * t * u
    ref13 = kap * u + 208 * x * t * u + y * (104 * t ** 2 + ep * 312 * u ** 2)
    say(f"     c1 = {c1};  c13 = {c13}")
    check(f"{C.lab} X0 c1, c13 equal I6 (whole family)", zero(c1 - ref1) and zero(c13 - ref13))
    # 2B against polarisation of the round-5 Q, one sample per op
    ok = True
    X = C.P(2, [random.randint(-2, 2) for _ in range(N)])
    Yv = C.P(12, [random.randint(-2, 2) for _ in range(N)])
    for op in OPS:
        Xs = [R5.add(X[a], Yv[a]) for a in range(N)]
        q, qx, qy = C.S.Q3(op, Xs), C.S.Q3(op, X), C.S.Q3(op, Yv)
        pol = [R5.add(R5.add(q[a], qx[a], -1), qy[a], -1) for a in range(N)]
        D = C.B2(op, X, Yv)
        ok &= all(R5.add(pol[a], D[a], -1) == {} for a in range(N))
    check(f"{C.lab} X0 direct bilinear 2B = polarisation of Q (random P_2, P_12 inputs, all four operators)", ok)
    return c1, c13


# ================================================================================================ X1: J1
def x1(C):
    say(f"X1  J1 closure of the even equation on the vector sector (all 14 basis covectors, both maps, four operators)")
    ok_fit, ok_same, ok_even = True, True, True
    src = {}
    for lab, P1 in (('t', C.Phi1), ('u', C.Phi1G)):
        ref = None
        for c in range(N):
            zeta = vwedge(e(c), P1)
            fits = {}
            for op in OPS:
                Yf = C.Yform(C.S.dual(op, zeta))
                ok_even &= all(pc(m) % 2 == 0 for a in range(N) for m in Yf[a])
                s, ok = C.fitV(Yf, c)
                ok_fit &= ok
                fits[op] = s
            ok_same &= ref is None or fits == ref
            ref = ref or fits
        src[lab] = ref
    S = sp.Matrix([[fam({op: src[lab][op][k] for op in OPS}) for lab in 'tu'] for k in KS])
    blk = {}
    for lab, To in (('t', C.Phi1), ('u', C.Phi1G)):
        ref = None
        for c in range(N):
            fits = {}
            for k in KS:
                X = C.P(k, e(c))
                for op in OPS:
                    Yf = C.Yform(C.B2(op, X, To), 3)
                    ok_even &= all(pc(m) % 2 == 0 for a in range(N) for m in Yf[a])
                    s, ok = C.fitV(Yf, c)
                    ok_fit &= ok
                    fits[(op, k)] = s
            ok_same &= ref is None or fits == ref
            ref = ref or fits
        blk[lab] = ref
    M1 = sp.Matrix([[fam({op: blk['t'][(op, k)][j] for op in OPS}) for k in KS] for j in KS])
    M13 = sp.Matrix([[fam({op: blk['u'][(op, k)][j] for op in OPS}) for k in KS] for j in KS])
    M = kap * sp.eye(4) + t * M1 + u * M13
    say(f"     S = {S.tolist()}")
    say(f"     M = {M.tolist()}")
    check(f"{C.lab} J1 outputs are g_ev-valued and exactly of the form sum_j s_j P_j(e^c) for every c", ok_fit and ok_even)
    check(f"{C.lab} J1 same coefficients for all 14 covectors (hence all v by linearity)", ok_same)
    ep = C.eps
    S_wu = sp.Matrix([[0, 0], [-24 * x, -24 * ep * y], [24 * I * w, -24 * I * z], [0, 0]])
    check(f"{C.lab} J1(a) S equals the write-up's matrix", mzero(S - S_wu))
    sig, th, de = t * x + ep * u * y, t * w - u * z, t * y - u * x
    M_wu = sp.Matrix([[kap, 0, 208 * ep * de, 0],
                      [0, kap + 224 * sig, 384 * I * ep * th, 16 * ep * de],
                      [16 * de, -384 * I * th, kap + 544 * sig, 0],
                      [0, 208 * de, 0, kap]])
    check(f"{C.lab} J1(b) M equals the write-up's compact form (sigma, theta, delta, epsilon = {ep:+d})", mzero(M - M_wu))
    detM = sp.expand(M.det())
    rho = kap ** 2 - 3328 * ep * de ** 2
    check(f"{C.lab} J1(b) det M = (rho + 224 k s)(rho + 544 k s) - eps (384 k theta)^2",
          zero(detM - ((rho + 224 * kap * sig) * (rho + 544 * kap * sig) - ep * (384 * kap * th) ** 2)))
    # z, w only through theta; A, B couple to C, D only through delta
    hs_, ht_, hd_ = sp.symbols('sigma_ theta_ delta_')
    Mc = sp.Matrix([
        [kap, 0, 208 * ep * hd_, 0], [0, kap + 224 * hs_, 384 * I * ep * ht_, 16 * ep * hd_],
        [16 * hd_, -384 * I * ht_, kap + 544 * hs_, 0], [0, 208 * hd_, 0, kap]])
    okth = mzero(Mc.subs({hs_: sig, ht_: th, hd_: de}) - M) and not ({z, w} & (sig.free_symbols | de.free_symbols))
    okAB = all(Mc[i, j].free_symbols <= {hd_} for i in (0, 3) for j in (1, 2)) and \
        all(Mc[i, j].free_symbols <= {hd_} for i in (1, 2) for j in (0, 3))
    check(f"{C.lab} J1(b) z, w enter M only through theta; A, B couple to C, D only through delta", okth and okAB)
    # Sigma = -M^-1 S
    Sig = -M.adjugate() * S / detM
    lo = -S / kap
    check(f"{C.lab} J1(d) leading order C = 24/k (x dt + eps' y du), D = 24i/k (-w dt + z du) [eps' = {ep:+d}]",
          zero(lo[1, 0] - 24 * x / kap) and zero(lo[1, 1] - 24 * ep * y / kap) and
          zero(lo[2, 0] + 24 * I * w / kap) and zero(lo[2, 1] - 24 * I * z / kap) and lo[0, 0] == 0 and lo[3, 1] == 0)
    SigA = Sig.subs(WA).applyfunc(sp.factor)
    dA = kap ** 2 + 224 * kap * t - 3328 * ep * u ** 2
    okWA = zero(SigA[1, 0] - 24 * kap / dA) and zero(SigA[3, 0] - 4992 * u / dA) and \
        all(zero(SigA[i, j]) for i in range(4) for j in range(2) if (i, j) not in ((1, 0), (3, 0)))
    check(f"{C.lab} J1(d) W_A: Sigma has only Sigma_Ct = 24k/(k^2+224kt-3328 eps u^2), Sigma_Bt = 4992u/(...)", okWA)
    SigS = Sig.subs(SH).applyfunc(sp.factor)
    numpt = {kap: 5, t: sp.Rational(1, 7), u: sp.Rational(-1, 3)}
    check(f"{C.lab} J1(d) Shiab#: all four components of T_e nonzero at a point with u != 0",
          all(any(sp.simplify(SigS[i, j].subs(numpt)) != 0 for j in range(2)) for i in range(4)))
    # reality of the slaved field in the real form
    if ep == 1:
        s_, r_ = sp.Rational(1, 7), sp.Rational(-1, 3)
        okr = True
        for mem in ({x: -I, y: 0, z: -I, w: 0}, {x: -2 * I, y: I, z: -3 * I, w: -I}):
            SN = Sig.subs(mem).subs({kap: 5, t: I * s_, u: I * r_}).applyfunc(sp.expand_complex)
            # dt = i ds: coefficient of ds is i Sigma; required: A in iR, C in R, D in iR, B in R
            ph = [sp.expand_complex(I * SN[i, j]) for i in range(4) for j in range(2)]
            okr &= all(sp.re(ph[2 * i + j]) == 0 for i in (0, 2) for j in range(2))
            okr &= all(sp.im(ph[2 * i + j]) == 0 for i in (1, 3) for j in range(2))
        check("(14,0) slaved T_e lies in u(128) (A, D imaginary, C, B real coefficients) for real (s, r), two members", okr)
    else:
        okr = True
        for mem in (SH, {x: 2, y: -1, z: 3, w: 1}):
            SN = Sig.subs(mem).subs({kap: 5, t: sp.Rational(1, 7), u: sp.Rational(-1, 3)}).applyfunc(sp.expand_complex)
            okr &= all(sp.re(SN[i, j]) == 0 for i in (0, 2) for j in range(2))
            okr &= all(sp.im(SN[i, j]) == 0 for i in (1, 3) for j in range(2))
        check("(7,7) slaved T_e lies in u(64,64) (A, D imaginary, C, B real coefficients) for real (t, u), two members", okr)
    return dict(S=S, M=M, M1=M1, M13=M13, detM=detM, Sig=Sig)


# ================================================================================================ X2: J2
def dzeta(C, k, p, q):
    """d_0 P_k(A) for nabla_a A_b = delta_ap delta_bq:  zeta_(c,a) = P_k(N_c.)_a - P_k(N_a.)_c"""
    Pq = C.P(k, e(q))
    zt = {}
    for c, a in itertools.combinations(range(N), 2):
        X = {}
        if c == p and Pq[a]:
            X = R5.add(X, Pq[a])
        if a == p and Pq[c]:
            X = R5.add(X, Pq[c], -1)
        if X:
            zt[(c, a)] = X
    return zt


WU_L = {  # write-up table of derivative coefficients (E1, E13, w3, w11)
    1: {0: (0, 0, -x, -y), 2: (24 * x, 24 * y, 11 * I * z, 11 * I * w), 12: (24 * I * w, -24 * I * z, 11 * y, -11 * x),
        14: (0, 0, -I * w, I * z)},
    -1: {0: (0, 0, -x, -y), 2: (24 * x, 24 * y, 11 * I * z, 11 * I * w), 12: (-24 * I * w, -24 * I * z, -11 * y, -11 * x),
         14: (0, 0, I * w, I * z)}}
R3 = sp.Rational(1, 3)
WU_B = {  # write-up table of quadratic coefficients: {(k,l): {'E1': (b1,b2,b3), 'E13': (...), 'w3': b4, 'w11': b4}}
    1: {(0, 2): {'w3': -22 * x * R3, 'w11': -22 * y * R3},
        (0, 12): {'E1': (0, -16 * y, 16 * y), 'E13': (0, 16 * x, -16 * x)},
        (0, 14): {'w3': 2 * y * R3, 'w11': -2 * x * R3},
        (2, 2): {'E1': (-224 * x, -224 * x, -176 * x), 'E13': (-224 * y, -224 * y, -176 * y)},
        (2, 12): {'E1': (-32 * I * w, -32 * I * w, -352 * I * w), 'E13': (32 * I * z, 32 * I * z, 352 * I * z),
                  'w3': -266 * y * R3, 'w11': 266 * x * R3},
        (2, 14): {'E1': (16 * y, 0, -16 * y), 'E13': (-16 * x, 0, 16 * x)},
        (12, 12): {'E1': (-160 * x, -160 * x, 528 * x), 'E13': (-160 * y, -160 * y, 528 * y)},
        (12, 14): {'w3': 22 * x * R3, 'w11': 22 * y * R3}},
    -1: {(0, 2): {'w3': -22 * x * R3, 'w11': -22 * y * R3},
         (0, 12): {'E1': (0, 16 * y, -16 * y), 'E13': (0, 16 * x, -16 * x)},
         (0, 14): {'w3': -2 * y * R3, 'w11': -2 * x * R3},
         (2, 2): {'E1': (-224 * x, -224 * x, -176 * x), 'E13': (-224 * y, -224 * y, -176 * y)},
         (2, 12): {'E1': (32 * I * w, 32 * I * w, 352 * I * w), 'E13': (32 * I * z, 32 * I * z, 352 * I * z),
                   'w3': 266 * y * R3, 'w11': 266 * x * R3},
         (2, 14): {'E1': (-16 * y, 0, 16 * y), 'E13': (-16 * x, 0, 16 * x)},
         (12, 12): {'E1': (160 * x, 160 * x, -528 * x), 'E13': (160 * y, 160 * y, -528 * y)},
         (12, 14): {'w3': -22 * x * R3, 'w11': -22 * y * R3}}}

PAIRS = {1: [(0, 1), (0, 0), (1, 0), (5, 5), (3, 9), (9, 3), (13, 2), (4, 11)],
         -1: [(0, 1), (0, 0), (1, 0), (7, 7), (8, 8), (0, 7), (7, 0), (7, 8), (8, 7), (10, 12), (2, 11), (12, 4)]}


def x2(C):
    eta, ep = C.eta, C.eps
    say(f"X2  J2 derivative map on all 196 basis tensors e^p (x) e^q, k = 0, 2, 12, 14, four operators")
    lin = {}
    ok_rec, ok_struct, ok_odd = True, True, True
    for k in KS:
        data = {}
        for p in range(N):
            for q in range(N):
                zt = dzeta(C, k, p, q)
                for op in OPS:
                    Yf = C.Yform(C.S.dual(op, zt))
                    ok_odd &= all(pc(m) % 2 == 1 for a in range(N) for m in Yf[a])
                    dec, ok = C.decomp(Yf)
                    ok_rec &= ok
                    data[(op, p, q)] = dec
        for op in OPS:
            E1, E13, W3, W11 = data[(op, 0, 1)]
            l = {'E1': E1[1][0], 'E13': E13[1][0], 'w3': W3.get((0, 1), 0), 'w11': W11.get((0, 1), 0)}
            lin[(op, k)] = l
            for p in range(N):
                for q in range(N):
                    E1, E13, W3, W11 = data[(op, p, q)]
                    for nm, Ea in (('E1', E1), ('E13', E13)):
                        ok_struct &= all(Ea[a][b] == l[nm] * ((a == q and b == p) -
                                                            (p == q) * eta[p] * (a == b) * eta[a])
                                         for a in range(N) for b in range(N))
                    for nm, Wd in (('w3', W3), ('w11', W11)):
                        ok_struct &= all(Wd.get((b, c), 0) == l[nm] * ((b == p and c == q) - (b == q and c == p))
                                         for b, c in itertools.combinations(range(N), 2))
    L = {(nm, k): fam({op: lin[(op, k)][nm] for op in OPS}) for k in KS for nm in CH}
    for k in KS:
        say(f"     k = {k:2d} ({KN[k]}): " + ";  ".join(f"{nm}: {L[(nm, k)]}" for nm in CH))
    check(f"{C.lab} J2 derivative outputs are g_odd-valued, only grades 1, 13, 3, 11 (exact reconstruction)",
          ok_rec and ok_odd)
    check(f"{C.lab} J2 derivative map is exactly l (N^T - tr N g) in E1, E13 and l dA in w3, w11 (all 196 tensors)",
          ok_struct)
    check(f"{C.lab} J2 derivative table equals the write-up's table",
          all(zero(L[(nm, k)] - WU_L[ep][k][i]) for k in KS for i, nm in enumerate(CH)))
    # ---- quadratic map on sampled pairs
    say(f"X2b J2 quadratic map 2B(P_k(e^p), P_l(e^q)) on sampled basis pairs {PAIRS[ep]} (p <= q when k = l)")
    bil = {}
    okq_rec, okq_struct, okq_odd = True, True, True
    for k, l in itertools.combinations_with_replacement(KS, 2):
        pairs = [pq for pq in PAIRS[ep] if k != l or pq[0] <= pq[1]]
        data = {}
        for p, q in pairs:
            X, Yv = C.P(k, e(p)), C.P(l, e(q))
            for op in OPS:
                Yf = C.Yform(C.B2(op, X, Yv), 3)
                okq_odd &= all(pc(m) % 2 == 1 for a in range(N) for m in Yf[a])
                dec, ok = C.decomp(Yf)
                okq_rec &= ok
                data[(op, p, q)] = dec
        for op in OPS:
            E1o, E13o, W3o, W11o = data[(op, 0, 1)]
            E1d, E13d, _, _ = data[(op, 0, 0)]
            b = {'E1': (E1o[0][1], E1o[1][0], E1d[1][1] * eta[0] * eta[1]),
                 'E13': (E13o[0][1], E13o[1][0], E13d[1][1] * eta[0] * eta[1]),
                 'w3': W3o.get((0, 1), 0), 'w11': W11o.get((0, 1), 0)}
            bil[(op, k, l)] = b
            for p, q in pairs:
                E1, E13, W3, W11 = data[(op, p, q)]
                for nm, Ea in (('E1', E1), ('E13', E13)):
                    b1, b2, b3 = b[nm]
                    okq_struct &= all(Ea[a][c] == b1 * (a == p and c == q) + b2 * (a == q and c == p) +
                                      b3 * (p == q) * eta[p] * (a == c) * eta[a] for a in range(N) for c in range(N))
                for nm, Wd in (('w3', W3), ('w11', W11)):
                    okq_struct &= all(Wd.get((c1_, c2_), 0) == b[nm] * ((c1_ == p and c2_ == q) - (c1_ == q and c2_ == p))
                                      for c1_, c2_ in itertools.combinations(range(N), 2))
    B = {}
    for k, l in itertools.combinations_with_replacement(KS, 2):
        for nm in ('E1', 'E13'):
            B[(nm, k, l)] = tuple(fam({op: bil[(op, k, l)][nm][j] for op in OPS}) for j in range(3))
        for nm in ('w3', 'w11'):
            B[(nm, k, l)] = fam({op: bil[(op, k, l)][nm] for op in OPS})
    for k, l in itertools.combinations_with_replacement(KS, 2):
        ent = [f"{nm}: {B[(nm, k, l)]}" for nm in ('E1', 'E13') if any(v != 0 for v in B[(nm, k, l)])]
        ent += [f"{nm}: {B[(nm, k, l)]}" for nm in ('w3', 'w11') if B[(nm, k, l)] != 0]
        say(f"     ({KN[k]}, {KN[l]}): " + ("; ".join(ent) if ent else "0"))
    check(f"{C.lab} J2 quadratic outputs are g_odd-valued, only grades 1, 13, 3, 11 (sampled pairs)", okq_rec and okq_odd)
    check(f"{C.lab} J2 quadratic map has the form b1 A(x)A' + b2 A'(x)A + b3 <A,A'> g, b4 A^A' (sampled pairs)", okq_struct)
    okt = True
    for k, l in itertools.combinations_with_replacement(KS, 2):
        wu = WU_B[ep].get((k, l), {})
        for nm in ('E1', 'E13'):
            ref = wu.get(nm, (0, 0, 0))
            okt &= all(zero(B[(nm, k, l)][j] - ref[j]) for j in range(3))
        for nm in ('w3', 'w11'):
            okt &= zero(B[(nm, k, l)] - wu.get(nm, 0))
    check(f"{C.lab} J2 quadratic table equals the write-up's table (all 10 pairs (k,l))", okt)
    tw = all(all(v == 0 for v in bil[(op, k, l)]['E1'] + bil[(op, k, l)]['E13']) and bil[(op, k, l)]['w3'] == 0
             and bil[(op, k, l)]['w11'] == 0
             for op in ('W_C.AdG', "W'_C.AdG") for k, l in itertools.combinations_with_replacement(KS, 2) if (k, l) != (2, 12))
    check(f"{C.lab} J2(d) the twisted pair enters Q(T_e) only through the (C, D) cross term", tw)
    return L, B


# ================================================================================================ X3: Cor. J2.1
def coeffs(Sig, dSt, dSu, L, B):
    """kinetic H, gradient Q, Q' and Lambda^2 coefficients varpi (own bookkeeping)"""
    lE = sp.Matrix([[L[(nm, k)] for k in KS] for nm in ('E1', 'E13')])
    H = lE * Sig
    dS = {'t': dSt, 'u': dSu}
    ph = 'tu'
    Q, Qp, W = {}, {}, {}
    for nm in ('E1', 'E13'):
        for a, fa in enumerate(ph):
            for b, fb in enumerate(ph):
                v = sum(L[(nm, k)] * dS[fb][i, a] for i, k in enumerate(KS))      # l d_psi Sigma_k,phi dphi(x)dpsi
                vp = -v
                for (i, k), (j, l) in itertools.combinations_with_replacement(list(enumerate(KS)), 2):
                    b1, b2, b3 = B[(nm, k, l)]
                    wkl = sp.Rational(1, 2) if k == l else 1
                    v += wkl * (b1 * Sig[i, a] * Sig[j, b] + b2 * Sig[j, a] * Sig[i, b])
                    vp += wkl * b3 * Sig[i, a] * Sig[j, b]
                Q[(nm, fa, fb)], Qp[(nm, fa, fb)] = v, vp
        # antisymmetric part: coefficient of dt^du
        v = sum(L[(nm, k)] * (dSu[i, 0] - dSt[i, 1]) for i, k in enumerate(KS)) / 2
        for (i, k), (j, l) in itertools.combinations(list(enumerate(KS)), 2):
            b1, b2, _ = B[(nm, k, l)]
            v += (b1 - b2) * (Sig[i, 0] * Sig[j, 1] - Sig[i, 1] * Sig[j, 0]) / 2
        W[nm] = v
    for nm in ('w3', 'w11'):
        v = sum(L[(nm, k)] * (dSt[i, 1] - dSu[i, 0]) for i, k in enumerate(KS))
        for (i, k), (j, l) in itertools.combinations(list(enumerate(KS)), 2):
            v += B[(nm, k, l)] * (Sig[i, 0] * Sig[j, 1] - Sig[i, 1] * Sig[j, 0])
        W[nm] = v
    return lE, H, Q, Qp, W


def lead2(expr):
    """coefficient of kappa^-2 (asserting nothing of lower order)"""
    eps_ = sp.Symbol('eps_')
    ex = sp.expand(sp.expand(expr).subs(kap, 1 / eps_))
    P = sp.Poly(ex, eps_)
    return sp.factor(P.coeff_monomial(eps_ ** 2) / kap ** 2), all(P.coeff_monomial(eps_ ** j) == 0 for j in (0, 1))


def x3(C, d1, L, B):
    ep = C.eps
    say("X3  Cor. J2.1 (own assembly from the X1/X2 tables)")
    S, M, M1, M13, detM, Sig = d1['S'], d1['M'], d1['M1'], d1['M13'], d1['detM'], d1['Sig']
    lE = sp.Matrix([[L[(nm, k)] for k in KS] for nm in ('E1', 'E13')])
    H0 = -lE * S / kap
    if ep == 1:
        H0_wu = 576 / kap * sp.Matrix([[x ** 2 + w ** 2, x * y - w * z], [x * y - w * z, y ** 2 + z ** 2]])
    else:
        H0_wu = 576 / kap * sp.Matrix([[x ** 2 - w ** 2, w * z - x * y], [x * y - w * z, z ** 2 - y ** 2]])
    check(f"{C.lab} J2.1(b) leading-order H = -l_E S / kappa equals the write-up's matrix", mzero(H0 - H0_wu))
    check(f"{C.lab} J2.1(b) det H_0 = 576^2 (xz {'+' if ep == 1 else '-'} yw)^2 / kappa^2",
          zero(H0.det() - 576 ** 2 * (x * z + ep * y * w) ** 2 / kap ** 2))
    # exact: det H is proportional to (xz +- yw)^2 for every (t, u) (l_E and S vanish off the C, D block)
    lCD, SCD = lE[:, 1:3], S[1:3, :]
    check(f"{C.lab} l_E and S vanish off the (C, D) block, so det H = det(l_CD) det((M^-1)_CD) det(S_CD) exactly",
          all(lE[i, j] == 0 for i in range(2) for j in (0, 3)) and all(S[i, j] == 0 for i in (0, 3) for j in range(2)))
    check(f"{C.lab} det l_CD * det S_CD = 576^2 (xz {'+' if ep == 1 else '-'} yw)^2",
          zero(lCD.det() * SCD.det() - 576 ** 2 * (x * z + ep * y * w) ** 2))
    # exact det H for the whole family (Jacobi: det (M^-1)_CD = det M_(A,B) / det M = kappa^2 / det M)
    rngH = random.Random(11)
    okdet = True
    for _ in range(3):
        pt_ = {v: sp.Rational(rngH.randint(-7, 7), rngH.randint(1, 4)) for v in (x, y, z, w, t, u)}
        pt_[kap] = sp.Rational(rngH.randint(1, 9), rngH.randint(1, 3))
        Mn, Sn, ln = M.subs(pt_), S.subs(pt_), lE.subs(pt_)
        Hn = -ln * Mn.inv() * Sn
        okdet &= zero(Hn.det() - (576 ** 2 * (x * z + ep * y * w) ** 2 * kap ** 2).subs(pt_) / Mn.det())
    check(f"{C.lab} exact, whole family: det H = 576^2 (xz {'+' if ep == 1 else '-'} yw)^2 kappa^2 / det M "
          f"(3 random rational points; the write-up states this only for Shiab# and at leading order)", okdet)
    if ep == 1:
        Ksym = sp.diag(1, 1) * H0
        xr, yr, zr, wr = sp.symbols("x' y' z' w'", real=True)
        Kr = sp.expand(-Ksym.subs({x: -I * xr, y: -I * yr, z: -I * zr, w: -I * wr}))
        Kwu = 576 / kap * sp.Matrix([[xr ** 2 + wr ** 2, xr * yr - wr * zr], [xr * yr - wr * zr, yr ** 2 + zr ** 2]])
        check("(14,0) J2.1(b) real form: Hessian coefficient -K with K = (576/kappa) Gram((x',w'),(y',-z'))", mzero(Kr - Kwu))
        check("(14,0) J2.1(b) det K = 576^2 (x'z' + y'w')^2 / kappa^2",
              zero(Kr.det() - 576 ** 2 * (xr * zr + yr * wr) ** 2 / kap ** 2))
        # K is (576/kappa) x Gram: positive semidefinite only for kappa > 0
        Kn = Kr.subs({xr: 1, yr: 0, zr: 1, wr: 0, kap: -1})
        check("(14,0) remark: for kappa < 0, K is negative definite (so 'positive semidefinite' needs kappa > 0)",
              all(v < 0 for v in Kn.eigenvals()))
    else:
        KK = sp.diag(1, -1) * H0
        check("(7,7) J2.1(b) diag(1,-1) H_0 symmetric with det = -576^2 (xz - yw)^2/kappa^2 (indefinite)",
              mzero(KK - KK.T) and zero(KK.det() + 576 ** 2 * (x * z - y * w) ** 2 / kap ** 2))
    # exact H for Shiab# and W_A
    SigS = Sig.subs(SH).applyfunc(sp.factor)
    HS = (lE.subs(SH) * SigS).applyfunc(sp.factor)
    dMS = sp.factor(detM.subs(SH))
    ue = -ep
    HS_wu = 576 * kap / dMS * sp.Matrix([[kap ** 2 + 544 * kap * t - 3328 * ep * u ** 2, -ep * 384 * kap * u],
                                         [-384 * kap * u, kap ** 2 + 224 * kap * t - 3328 * ep * u ** 2]])
    check(f"{C.lab} J2.1(c) exact H for Shiab# equals (576 k/det M)[[k^2+544kt-3328 eps u^2, -eps 384ku],"
          f"[-384ku, k^2+224kt-3328 eps u^2]]", mzero(HS - HS_wu))
    check(f"{C.lab} J2.1(c) det H# = 576^2 kappa^2 / det M#", zero(HS.det() - 576 ** 2 * kap ** 2 / dMS))
    HA = (lE.subs(WA) * Sig.subs(WA)).applyfunc(sp.factor)
    HA_wu = 576 * kap / (kap ** 2 + 224 * kap * t - 3328 * ep * u ** 2) * sp.Matrix([[1, 0], [0, 0]])
    check(f"{C.lab} J2.1(c) H for W_A = 576k/(k^2+224kt-3328 eps u^2) diag(1, 0) (rank 1)", mzero(HA - HA_wu))
    _ = ue
    # leading-order Lambda^2 coefficients (Neumann series)
    K_ = t * M1 + u * M13
    Sig_lo = -S / kap + K_ * S / kap ** 2
    _, _, Qlo, Qplo, Wlo = coeffs(Sig_lo, M1 * S / kap ** 2, M13 * S / kap ** 2, L, B)
    fz = x * z + ep * y * w
    if ep == 1:
        wu = {'E1': -46080 * w * fz, 'E13': 46080 * z * fz, 'w3': -88320 * I * y * fz, 'w11': 88320 * I * x * fz}
    else:
        wu = {'E1': 46080 * w * fz, 'E13': 46080 * z * fz, 'w3': 88320 * I * y * fz, 'w11': 88320 * I * x * fz}
    okw, okord = True, True
    for nm in CH:
        lv, okl = lead2(Wlo[nm])
        okord &= okl
        okw &= zero(lv - wu[nm] / kap ** 2)
        say(f"     varpi[{nm}] leading = {lv}")
    check(f"{C.lab} J2.1(e) leading Lambda^2 coefficients equal the write-up's (all O(kappa^-2))", okw and okord)
    lq, okl = lead2(Qlo[('E1', 't', 't')])
    check(f"{C.lab} J2.1(d) gradient terms O(kappa^-2) and z, w-dependent (e.g. leading Q^(1)_tt = {lq})",
          okl and bool({z, w} & lq.free_symbols))
    # exact varpi for Shiab# and W_A
    dSt, dSu = SigS.diff(t), SigS.diff(u)
    _, _, _, _, WS = coeffs(SigS, dSt, dSu, {k_: v.subs(SH) for k_, v in L.items()},
                            {k_: (tuple(c.subs(SH) for c in v) if isinstance(v, tuple) else v.subs(SH)) for k_, v in B.items()})
    w13 = sp.factor(sp.simplify(WS['E13']))
    ref13 = -9216 * kap ** 2 * (-5 * kap ** 4 + 1600 * kap ** 3 * t + 609280 * kap ** 2 * t ** 2 + ep * 930304 * kap ** 2 * u ** 2
                                + ep * 56016896 * kap * t * u ** 2 - 587005952 * u ** 4) / dMS ** 2
    say(f"     Shiab# exact varpi[E13] = {w13}")
    check(f"{C.lab} J2.1(e) exact varpi^(13) for Shiab# equals the printed formula; = 46080/kappa^2 at T_o = 0",
          zero(w13 - ref13) and zero(w13.subs({t: 0, u: 0}) - 46080 / kap ** 2))
    SigA = Sig.subs(WA).applyfunc(sp.factor)
    _, _, _, _, WA_ = coeffs(SigA, SigA.diff(t), SigA.diff(u), {k_: v.subs(WA) for k_, v in L.items()},
                             {k_: (tuple(c.subs(WA) for c in v) if isinstance(v, tuple) else v.subs(WA)) for k_, v in B.items()})
    wa1 = sp.factor(sp.simplify(WA_['E1']))
    say(f"     W_A exact varpi: E1 = {wa1}; E13 = {sp.factor(sp.simplify(WA_['E13']))}; w3 = "
        f"{sp.factor(sp.simplify(WA_['w3']))}; w11 = {sp.factor(sp.simplify(WA_['w11']))}")
    check(f"{C.lab} observation: for W_A the exact varpi^(1) is nonzero for u != 0 (O(kappa^-3)), so dt^du = 0 is "
          f"forced for W_A too wherever u != 0", not zero(wa1) and zero(wa1.subs(u, 0)))
    return dict(lE=lE, H0=H0, SigS=SigS, dMS=dMS)


# ================================================================================================ X4: J4
lam, cH, cV, Lg = sp.symbols('lambda c_H c_V Lambda_g')
n4 = 4


def vertical_traces():
    """tr_{g_V} G_0 and tr_{g_V} tau(x)tau from an explicit 10 x 10 Gram computation at g = 1 (n = 4)"""
    basis = []
    for i in range(n4):
        for j in range(i, n4):
            Mx = sp.zeros(n4, n4)
            Mx[i, j] = Mx[j, i] = 1
            basis.append(Mx)
    G0 = sp.Matrix(10, 10, lambda a, b: (basis[a] * basis[b]).trace())
    TT = sp.Matrix(10, 10, lambda a, b: basis[a].trace() * basis[b].trace())
    gV = cV * (G0 - lam * TT)
    gVi = gV.inv()
    return sp.simplify((gVi * G0).trace()), sp.simplify((gVi * TT).trace())


def geometry(trG0, trTT):
    """G^Y on s(X) for an Einstein base Ric(g) = Lambda_g g in the basis (h, G_0, tau tau), from F3 / F4 inputs"""
    geo = {}
    for Z in ('Z1', 'Z2'):
        if Z == 'Z2':
            e1 = 0
            ricV = (-1, sp.Rational(1, 4))                       # -(n/4) G_{1/n}       (F1.4 / F3)
            sigma_F = -9                                         # s_Y = R/c_H - 9/c_V  (F3)
            beta = 0                                             # II = 0
        else:
            e1 = cH / (4 * cV * (n4 * lam - 1))                  # F4(c)
            ricV = (-sp.Rational(5, 4), sp.Rational(1, 4))       # -(5/4) G_{1/5}       (F4(c))
            sigma_F = -sp.Rational(23, 2) - (1 + lam) / (1 - 4 * lam)      # F4(b)
            mu = lam / (1 - n4 * lam)
            beta = (1 + n4 * mu) / (2 * cV)                      # Hess_HH = -a tau(II) = beta a h   (F4(a))
        ric = ((Lg + e1) / cH, ricV[0], ricV[1])
        s_tr = n4 * ric[0] + trG0 * ric[1] + trTT * ric[2]      # trace of the F Ricci blocks
        s_F = 4 * Lg / cH + sigma_F / cV
        gal = (1, cV, -lam * cV)
        G = tuple(sp.simplify(ric[j] - s_F / 2 * gal[j]) for j in range(3))
        geo[Z] = dict(G=G, beta=sp.simplify(beta), consistent=zero(s_tr - s_F), gal=gal)
    return geo


def six(G, gal, beta, trG0, trTT, xi, ci, lrow, Brow, Sig, dS, a, p, q):
    """(h, G_0, tau tau) coefficients of E^(i) at a point of s(X) under (J)"""
    ph = ('t', 'u')
    alpha = [sum(Sig[i, j] * a[f] for j, f in enumerate(ph)) for i in range(4)]
    E = [-xi * G[j] + ci * gal[j] for j in range(3)]
    for i, k in enumerate(KS):
        nh = sum(Sig[i, j] * beta * a[f] for j, f in enumerate(ph))
        nG = sum(Sig[i, j] * p[f] for j, f in enumerate(ph))
        nT = sum(Sig[i, j] * q[f] for j, f in enumerate(ph)) + \
            sum(dS[g_][i, j] * a[g_] * a[f] for j, f in enumerate(ph) for g_ in ph)
        trN = n4 * nh + trG0 * nG + trTT * nT
        for j, nv in enumerate((nh, nG, nT)):
            E[j] += lrow[k] * (nv - trN * gal[j])
    for (i, k), (j, l) in itertools.combinations_with_replacement(list(enumerate(KS)), 2):
        b1, b2, b3 = Brow[(k, l)]
        wkl = sp.Rational(1, 2) if k == l else 1
        E[2] += wkl * (b1 + b2) * alpha[i] * alpha[j]
        for jj in range(3):
            E[jj] += wkl * b3 * alpha[i] * alpha[j] * trTT * gal[jj]
    return E


def x4(res):
    say("\n######## J4 (reduction along s(X), n = 4; own bookkeeping; F3/F4 inputs taken as given) ########")
    trG0, trTT = vertical_traces()
    check("vertical traces: tr G_0 = (9 + 1/(1-4 lam))/c_V, |tau|^2 = 4/(c_V(1-4 lam)) (explicit 10x10 Gram)",
          zero(trG0 - (9 + 1 / (1 - 4 * lam)) / cV) and zero(trTT - 4 / (cV * (1 - 4 * lam))))
    geo = geometry(trG0, trTT)
    check("F inputs consistent: trace of the F3/F4(c) Ricci blocks equals s_Y of F3/F4(b), in Z1 and Z2",
          geo['Z1']['consistent'] and geo['Z2']['consistent'])
    for Z in ('Z1', 'Z2'):
        say(f"     {Z}: G^Y = ({sp.factor(geo[Z]['G'][0])}) h + ({sp.factor(geo[Z]['G'][1])}) G_0 + "
            f"({sp.factor(geo[Z]['G'][2])}) tau tau;  beta = {geo[Z]['beta']}")
    check("Z1 beta = 1/(2 c_V (1 - 4 lam)) from II of F4(a); Z2 beta = 0",
          zero(geo['Z1']['beta'] - 1 / (2 * cV * (1 - 4 * lam))) and geo['Z2']['beta'] == 0)
    # generic symbols: Sigma, dSigma, l, b, x_i, c_i
    Sg = sp.Matrix(4, 2, lambda i, j: sp.Symbol(f"S{KS[i]}{'tu'[j]}"))
    dS = {f: sp.Matrix(4, 2, lambda i, j: sp.Symbol(f"dS{f}{KS[i]}{'tu'[j]}")) for f in 'tu'}
    lg = {k: sp.Symbol(f"l{k}") for k in KS}
    Bg = {(k, l): tuple(sp.Symbol(f"b{k}_{l}_{j}") for j in range(3))
          for k, l in itertools.combinations_with_replacement(KS, 2)}
    xi, ci = sp.symbols('x_i c_i')
    at, au, pt, pu, qt, qu = sp.symbols('a_t a_u p_t p_u q_t q_u')
    a, p, q = {'t': at, 'u': au}, {'t': pt, 'u': pu}, {'t': qt, 'u': qu}
    Hrow = [sum(lg[k] * Sg[i, j] for i, k in enumerate(KS)) for j in range(2)]
    okb, okc = True, True
    rhos = {}
    Theta = {}
    for Z in ('Z1', 'Z2'):
        G, gal, beta = geo[Z]['G'], geo[Z]['gal'], geo[Z]['beta']
        E = six(G, gal, beta, trG0, trTT, xi, ci, lg, Bg, Sg, dS, a, p, q)
        lhs = E[0] - E[1] / cV
        rho = sp.simplify(G[0] - G[1] / cV - Lg / cH)
        rhos[Z] = rho
        rhs = -xi * (Lg / cH + rho) + sum(Hrow[j] * (beta * a[f] - p[f] / cV) for j, f in enumerate('tu'))
        okb &= zero(lhs - rhs) and not sp.diff(rho, Lg)
        # a = 0 branch: linear in (c_i, P_i, Q_i); solve and compare with the stated closed form
        E0 = [sp.expand(v.subs({at: 0, au: 0})) for v in E]
        Pi, Qi = sp.symbols('P_i Q_i')
        E0s = [sp.expand(v.subs({pt: Pi / Hrow[0], pu: 0, qt: Qi / Hrow[0], qu: 0})) for v in E0]
        sol = sp.solve(E0s, [ci, Pi, Qi], dict=True)
        Th = sp.simplify(trG0 * G[1] + trTT * G[2] - 9 * G[0])
        Theta[Z] = Th
        okc &= len(sol) == 1 and zero(sol[0][ci] - xi * Th) and zero(sol[0][Pi] - xi * (G[1] - cV * G[0])) and \
            zero(sol[0][Qi] - xi * (G[2] + lam * cV * G[0]))
        say(f"     {Z}: rho = {sp.factor(rho)};  Theta = {sp.factor(Th.subs(Lg, 0))} + ({sp.simplify(sp.diff(Th, Lg))}) Lambda_g")
    check("J4(b) balance identity [h] - [G_0]/c_V: x_i (Lambda_g/c_H + rho) = sum H_i,phi (beta a_phi - p_phi/c_V), "
          "identity in generic Sigma, dSigma, l and b (Z1, Z2)", okb)
    check("J4(b) rho_Z2 = 1/c_V, rho_Z1 = (5 lam - 1)/(c_V (4 lam - 1))",
          zero(rhos['Z2'] - 1 / cV) and zero(rhos['Z1'] - (5 * lam - 1) / (cV * (4 * lam - 1))))
    check("J4(c) a = 0 branch: c_i = x_i Theta, sum H p = x_i (G_G0 - c_V G_h), sum H q = x_i (G_tt + lam c_V G_h) "
          "is the unique solution (generic Sigma, l, b)", okc)
    T0Z2 = -sp.Rational(9, 2) / cV
    T0Z1 = 3 * (2 - 15 * lam) / (2 * cV * (4 * lam - 1))
    check("J4(c) Theta = Theta_0 - (11/c_H) Lambda_g with Theta_0(Z2) = -9/(2c_V), Theta_0(Z1) = 3(2-15lam)/(2c_V(4lam-1)); "
          "Theta_0(Z1, 1/5) = 15/(2c_V)",
          zero(Theta['Z2'] - (T0Z2 - 11 * Lg / cH)) and zero(Theta['Z1'] - (T0Z1 - 11 * Lg / cH)) and
          zero(T0Z1.subs(lam, sp.Rational(1, 5)) - sp.Rational(15, 2) / cV))
    G6 = [gv.subs({lam: sp.Rational(1, 5), Lg: 0}) for gv in geo['Z1']['G']]
    gal6 = [1, cV, -cV / 5]
    check("J4(c) at the F6 point (Z1, lam = 1/5, Lambda_g = 0): G^Y = (15/(2c_V)) g_aleph and the p, q right-hand sides "
          "vanish", all(zero(G6[j] - sp.Rational(15, 2) / cV * gal6[j]) for j in range(3)) and
          zero(G6[1] - cV * G6[0]) and zero(G6[2] + cV / 5 * G6[0]))
    check("J4(e) radial Z1 at a = 0: Lambda_g = c_H (1 - 5 lam)/(c_V (4 lam - 1)), zero exactly at lam = 1/5",
          zero(-cH * rhos['Z1'] - cH * (1 - 5 * lam) / (cV * (4 * lam - 1))) and
          sp.solve(sp.numer(sp.together(-cH * rhos['Z1'])), lam) == [sp.Rational(1, 5)])

    # ---------------------------------------------------------------- J4(c) image condition, exactly
    for ep, d in res.items():
        lE = d['lE']
        fz = x * z + ep * y * w
        okim = True
        for j in range(4):
            dt_ = sp.expand(sp.Matrix([[lE[0, j], x], [lE[1, j], y]]).det())
            okim &= dt_ == 0 or sp.fraction(sp.cancel(dt_ / fz))[1] == 1
        check(f"{d['lab']} J4(c) image condition exact: every column of l_E is (x,y) up to a multiple of (xz {'+' if ep == 1 else '-'} yw), "
              f"so im H(t,u) is contained in span(x,y) whenever det H_0 = 0 (not only at leading order)", okim)

    # ---------------------------------------------------------------- J4(d) ranges and det M on the branches
    for ep, d in res.items():
        c1, c13, detM = d['c1'], d['c13'], d['detM']
        c1S, c13S, dMS = c1.subs(SH), sp.factor(c13.subs(SH)), sp.factor(detM.subs(SH))
        if ep == -1:
            ok = zero(c13S - u * (kap + 208 * t))
            ok &= zero(dMS.subs(u, 0) - kap ** 2 * (kap + 224 * t) * (kap + 544 * t))
            vv = sp.Symbol('v', nonnegative=True)
            dbr = sp.expand(dMS.subs(t, -kap / 208))
            ok &= zero(dbr - (vv ** 2 + sp.Rational(554, 13) * kap ** 2 * vv + sp.Rational(21, 169) * kap ** 4).subs(vv, 3328 * u ** 2))
            ok &= zero(c1S.subs(u, 0) - (kap * t + 312 * t ** 2))
            ok &= zero(c1S.subs(t, -kap / 208) - (kap ** 2 / 416 - 104 * u ** 2))
            check("(7,7) J4(d) Shiab#: curve u(k+208t) = 0; det M = k^2(k+224t)(k+544t) on u = 0, = v^2+554/13 k^2 v+21/169 k^4 "
                  "> 0 on t = -k/208; c1 covers [-k^2/1248, oo) (t >= -k/624, det M != 0) and (-oo, k^2/416]: every real "
                  "Lambda_g", ok and -sp.Rational(1, 1248) < sp.Rational(1, 416))
        else:
            s_, r_ = sp.symbols('s r', real=True)
            sub = {x: -I, y: 0, z: -I, w: 0, t: I * s_, u: I * r_}
            e1 = sp.expand(-I * c1.subs(sub))
            e13 = sp.factor(sp.expand(-I * c13.subs(sub)))
            dMr = sp.factor(sp.expand(detM.subs(sub)))
            ok = zero(e1 - (kap * s_ + 312 * s_ ** 2 + 104 * r_ ** 2)) and zero(e13 - r_ * (kap + 208 * s_))
            ok &= zero(dMr.subs(r_, 0) - kap ** 2 * (kap + 224 * s_) * (kap + 544 * s_))
            ok &= zero(sp.expand(e1.subs(s_, -kap / 208)) - (kap ** 2 / 416 + 104 * r_ ** 2))
            check("(14,0) J4(d) Shiab# real form: e1 = ks + 312s^2 + 104r^2 on r(k+208s) = 0; e1 >= -k^2/1248; det M = "
                  "k^2(k+224s)(k+544s) on r = 0 (nonzero for s >= -k/624, k > 0): a = 0 branch gives exactly "
                  "[c_H/11 (Theta_0 - k^2/1248), oo)", ok)
    return geo, trG0, trTT


def x4b(res, geo, trG0, trTT):
    """J4 with a != 0 (nonzero normal gradient, still ansatz (J)): closed form, and the scope of the (14,0) half-line"""
    say("\nX4b J4 beyond the a = 0 branch (ansatz (J) with a_phi != 0)")
    # (i) closed form: eliminating p, q (H invertible), the six equations are c_i = x_i Theta + 13 beta (H a)_i - |tau|^2 B_i(alpha),
    #     B_i(alpha) = sum_{k<=l} w_kl (b1 + b2 + b3)^(i)_kl alpha_k alpha_l, alpha_k = sum_phi Sigma_k,phi a_phi (dSigma drops out)
    rng = random.Random(7)

    def R():
        return sp.Rational(rng.randint(-9, 9), rng.randint(1, 5))
    ci_, pt_, qt_ = sp.symbols('c_i p_t q_t')
    okcf = True
    for Z in ('Z1', 'Z2'):
        for _trial in range(4):
            num = {lam: sp.Rational(rng.randint(-5, 0), 7), cH: R() + 11, cV: R() + 13, Lg: R()}
            G = [g.subs(num) for g in geo[Z]['G']]
            gal = [sp.sympify(g).subs(num) for g in geo[Z]['gal']]
            beta = sp.sympify(geo[Z]['beta']).subs(num)
            tg, tt = trG0.subs(num), trTT.subs(num)
            Sg = sp.Matrix(4, 2, lambda i, j: R())
            dS = {f: sp.Matrix(4, 2, lambda i, j: R()) for f in 'tu'}
            lg = {k: R() for k in KS}
            Bg = {(k, l): (R(), R(), R()) for k, l in itertools.combinations_with_replacement(KS, 2)}
            xi, a = R(), {'t': R(), 'u': R()}
            E = six(G, gal, beta, tg, tt, xi, ci_, lg, Bg, Sg, dS, a, {'t': pt_, 'u': R()}, {'t': qt_, 'u': R()})
            sol = sp.solve(E, [ci_, pt_, qt_], dict=True)
            H = [sum(lg[k] * Sg[i, j] for i, k in enumerate(KS)) for j in range(2)]
            al = [Sg[i, 0] * a['t'] + Sg[i, 1] * a['u'] for i in range(4)]
            Bf = sum((sp.Rational(1, 2) if k == l else 1) * sum(Bg[(k, l)]) * al[i] * al[j]
                     for (i, k), (j, l) in itertools.combinations_with_replacement(list(enumerate(KS)), 2))
            Th = tg * G[1] + tt * G[2] - 9 * G[0]
            okcf &= len(sol) == 1 and sol[0][ci_] == xi * Th + 13 * beta * (H[0] * a['t'] + H[1] * a['u']) - tt * Bf
    check("J4, a != 0: the six equations reduce to c_i = x_i Theta + 13 beta (H a)_i - |tau|^2 B_i(alpha) (dSigma drops out); "
          "exact identity test at 4 random rational points per metric (Z1, Z2), generic Sigma, dSigma, l, b", okcf)
    # (ii) Z2, Shiab#, (14,0): B_1(alpha)/x = -312 alpha_C^2 + 104 alpha_D^2, alpha_C real, alpha_D imaginary in u(128)
    d = res[1]
    Bsum = {kl: sp.expand(sum(d['B'][('E1',) + kl]).subs({y: 0, w: 0})) for kl in itertools.combinations_with_replacement(KS, 2)}
    okB = all(zero(v) for kl, v in Bsum.items() if kl not in ((2, 2), (12, 12))) and \
        zero(Bsum[(2, 2)] / 2 + 312 * x) and zero(Bsum[(12, 12)] / 2 - 104 * x)
    check("(14,0) Z2, x W_A + z W_C.AdG (incl. Shiab#): Theta = -e1 - 312|tau|^2 alpha_C^2 - 104|tau|^2 |alpha_D|^2 <= -e1 for every a "
          "(alpha_C real, alpha_D imaginary, X1 reality check), so the half-line of J4(d) is exact in Z2 also for a != 0 "
          "(an argument the write-up does not give)", okB)
    # (ii') with a != 0 the cosmological constant at s(X) depends on (z, w): Z2, T_o = 0 on s(X) (Sigma = -S/kappa), a = (0, 1)
    Sig0 = d['Sig'].subs({t: 0, u: 0})
    al0 = [Sig0[i, 1] for i in range(4)]
    Bf0 = sp.expand(sum((sp.Rational(1, 2) if k == l else 1) * sum(d['B'][('E1', k, l)]) * al0[i] * al0[j]
                        for (i, k), (j, l) in itertools.combinations_with_replacement(list(enumerate(KS)), 2)))
    say(f"     Z2, (14,0), T_o = 0 on s(X), a = (0, 1): Theta - c1/x = |tau|^2 B_1/x with B_1 = {sp.factor(Bf0)}")
    check("J4 with a != 0: the effective cosmological constant at s(X) depends on (z, w) (here through B_1; in Z1 also "
          "through 13 beta (H a)), so 'never on the cosmological term' holds only for the zeroth-order term c_i g",
          bool({z, w} & Bf0.free_symbols))
    # (iii) Z1, Shiab#, (14,0): exact real solution of all six equations below the a = 0 half-line
    sub = {x: -I, y: 0, z: -I, w: 0, kap: 1}
    num = {cH: 1, cV: 1, lam: 0}
    G = [g.subs(num) for g in geo['Z1']['G']]
    gal = [sp.sympify(g).subs(num) for g in geo['Z1']['gal']]
    beta = sp.sympify(geo['Z1']['beta']).subs(num)
    tg, tt = trG0.subs(num), trTT.subs(num)
    sv = sp.Rational(1, 1000)
    SigN = d['Sig'].subs(sub).subs({t: I * sv, u: 0}).applyfunc(sp.simplify)
    dSt = d['Sig'].subs(sub).diff(t).subs({t: I * sv, u: 0}).applyfunc(sp.simplify)
    dSu = d['Sig'].subs(sub).diff(u).subs({t: I * sv, u: 0}).applyfunc(sp.simplify)
    vv = -beta / (48 * tt)                                  # optimal v = a_s/(kappa + 224 s) on this slice
    a_s = vv * (1 + 224 * sv)
    a = {'t': I * a_s, 'u': 0}
    c1 = d['c1'].subs(sub).subs({t: I * sv, u: 0})
    c13 = d['c13'].subs(sub).subs({t: I * sv, u: 0})
    lrow1 = {k: sp.sympify(d['L'][('E1', k)]).subs(sub) for k in KS}
    lrow13 = {k: sp.sympify(d['L'][('E13', k)]).subs(sub) for k in KS}
    B1 = {kl: tuple(sp.sympify(c).subs(sub) for c in d['B'][('E1',) + kl]) for kl in itertools.combinations_with_replacement(KS, 2)}
    B13 = {kl: tuple(sp.sympify(c).subs(sub) for c in d['B'][('E13',) + kl]) for kl in itertools.combinations_with_replacement(KS, 2)}
    T0Z1 = (3 * (2 - 15 * lam) / (2 * cV * (4 * lam - 1))).subs(num)
    H1 = [sum(lrow1[k] * SigN[i, j] for i, k in enumerate(KS)) for j in range(2)]
    al = [SigN[i, 0] * a['t'] for i in range(4)]
    Bf1 = sum((sp.Rational(1, 2) if k == l else 1) * sum(B1[(k, l)]) * al[i] * al[j]
              for (i, k), (j, l) in itertools.combinations_with_replacement(list(enumerate(KS)), 2))
    Th = sp.simplify((c1 - 13 * beta * H1[0] * a['t'] + tt * Bf1) / (-I))
    Lval = sp.nsimplify(sp.simplify((T0Z1 - Th) / 11))
    pt2, qt2 = sp.symbols('p_t q_t')
    E1 = six([g.subs(Lg, Lval) for g in G], gal, beta, tg, tt, -I, c1, lrow1, B1, SigN, {'t': dSt, 'u': dSu}, a,
             {'t': pt2, 'u': 0}, {'t': qt2, 'u': 0})
    E13 = six([g.subs(Lg, Lval) for g in G], gal, beta, tg, tt, 0, c13, lrow13, B13, SigN, {'t': dSt, 'u': dSu}, a,
              {'t': pt2, 'u': 0}, {'t': qt2, 'u': 0})
    sol = sp.solve([sp.expand(v) for v in E1[:2]], [pt2, qt2], dict=True)
    ok = len(sol) == 1
    if ok:
        allz = all(zero(sp.expand(v.subs(sol[0]))) for v in E1 + E13)
        real = sp.im(Lval) == 0 and all(sp.re(sp.expand_complex(v)) == 0 for v in sol[0].values())
        bound0 = (T0Z1 - sp.Rational(1, 1248)) / 11
        dM = sp.simplify(d['detM'].subs(sub).subs({t: I * sv, u: 0}))
        say(f"     Z1, lambda = 0, c_H = c_V = kappa = 1, Shiab#, t = i/1000, u = 0, a_t = i*({a_s}), a_u = 0: Lambda_g = {Lval} "
            f"= {float(Lval):.6f};  p_t = {sp.expand_complex(sol[0][pt2])}, q_t = {sp.expand_complex(sol[0][qt2])};  det M = {dM};"
            f"  a = 0 lower bound {bound0} = {float(bound0):.6f}")
        ok = allz and real and Lval < bound0 and dM != 0
    lowZ1 = sp.simplify(78 * geo['Z1']['beta'] ** 2 / trTT)
    check("(14,0) J4(d) 'In (14,0) it takes exactly the values in [c_H/11 (Theta_0 - k^2/1248), oo)' fails in Z1 for ansatz (J): "
          "exact real solution of all six equations (a != 0, det M != 0) with Lambda_g below that bound", ok)
    check("(14,0) Z1, slice u = a_u = 0: the a != 0 range extends down to c_H/11 (Theta_0 - k^2/1248 - 39/(8 c_V (1 - 4 lam))) "
          "(78 beta^2/|tau|^2 = 39/(8 c_V (1 - 4 lam)))", zero(lowZ1 - sp.Rational(39, 8) / (cV * (1 - 4 * lam))))


# ================================================================================================ X5: J5
def x5():
    say("\n######## J5 (radial ansatz, Z2) ########")
    # (a) Euler field algebra
    lam_ = sp.Symbol('lambda')
    yv = sp.Matrix([[3, 1, 0, 2], [1, 2, 1, 0], [0, 1, 4, 1], [2, 0, 1, 5]])
    kv = sp.Matrix([[1, -2, 0, 1], [-2, 0, 3, 1], [0, 3, -1, 2], [1, 1, 2, 4]])
    Gl = lambda A_, B_: (yv.inv() * A_ * yv.inv() * B_).trace() - lam_ * (yv.inv() * A_).trace() * (yv.inv() * B_).trace()  # noqa: E731
    check("J5(a) G_lam|_y(y, k) = (1 - 4 lam) tr(y^-1 k) and G_lam|_y(y, y) = 4(1 - 4 lam), so |ds|^2 = 1/(4 c_V (1 - 4 lam))",
          zero(Gl(yv, kv) - (1 - 4 * lam_) * (yv.inv() * kv).trace()) and zero(Gl(yv, yv) - 4 * (1 - 4 * lam_)))
    # (c) B^T y + y B = 0 for y = E_ii forces B = 0
    Bs = sp.Matrix(4, 4, lambda i, j: sp.Symbol(f"B{i}{j}"))
    eqs = []
    for i in range(4):
        Eii = sp.zeros(4, 4)
        Eii[i, i] = 1
        eqs += list(Bs.T * Eii + Eii * Bs)
    sol = sp.solve(eqs, list(Bs), dict=True)
    check("J5(c) B^T y + y B = 0 for y = e_i (x) e_i (i = 1..4) forces B = 0", len(sol) == 1 and all(v == 0 for v in sol[0].values()))
    # (d) order-eps^2 leaf condition for constant curvature bases, signatures (4,0), (1,3), (2,2)
    out_ = {}
    for sig in ((1, 1, 1, 1), (1, -1, -1, -1), (1, 1, -1, -1)):
        g = sp.diag(*sig)
        gi = g.inv()
        hs = sp.symbols('h0:10')
        h = sp.zeros(4, 4)
        it = iter(hs)
        for i in range(4):
            for j in range(i, 4):
                h[i, j] = h[j, i] = next(it)
        h = h.subs(hs[0], sp.solve((gi * h).trace(), hs[0])[0])

        def R(a_, b_):     # R(e_a, e_b) Z = g(e_b, Z) e_a - g(e_a, Z) e_b  (curvature +1)
            Mx = sp.zeros(4, 4)
            for c in range(4):
                for dd in range(4):
                    Mx[c, dd] = (c == a_) * g[b_, dd] - (c == b_) * g[a_, dd]
            return Mx

        def act(Bm, X):
            return -(Bm.T * X + X * Bm)
        trless = all(zero((gi * act(R(a_, b_), h)).trace()) for a_ in range(4) for b_ in range(4))
        F = {}
        for a_ in range(4):
            for b_ in range(a_, 4):
                F[(a_, b_)] = sp.expand(sum(sig[i] * (gi * act(R(a_, i), h) * gi * act(R(b_, i), h)).trace() for i in range(4)))
        nz = [ab for ab, v in F.items() if v != 0]
        out_[sig] = (trless, nz, F[(0, 0)])
        say(f"     base {sig}: R.h traceless: {trless};  F(e_a, e_b; h) not identically 0 for (a,b) in {nz};  F(e0,e0) = {sp.factor(F[(0, 0)])}")
    check("J5(d) R.h is g-traceless (the lambda-term drops), all three signatures", all(v[0] for v in out_.values()))
    check("J5(d) leaf condition fails for constant curvature k != 0, Riemannian base", bool(out_[(1, 1, 1, 1)][1]))
    check("J5(d) leaf condition fails for constant curvature k != 0, Lorentzian base", bool(out_[(1, -1, -1, -1)][1]))
    check("J5(d) leaf condition fails for constant curvature k != 0, base (2,2) (also gives (7,7), A5; not in the write-up)",
          bool(out_[(1, 1, -1, -1)][1]))


# ================================================================================================ main
def main():
    res = {}
    for eta in (R5.EUC, R5.E77):
        C = Ctx(eta)
        t1 = time.time()
        c1, c13 = x0(C)
        d1 = x1(C)
        L, B = x2(C)
        d3 = x3(C, d1, L, B)
        res[C.eps] = dict(lab=C.lab, c1=c1, c13=c13, detM=d1['detM'], Sig=d1['Sig'], L=L, B=B, lE=d3['lE'])
        say(f"  [{C.lab} done, {time.time() - t1:.0f} s]")
    geo, trG0, trTT = x4(res)
    x4b(res, geo, trG0, trTT)
    x5()
    nf = [lab for lab, ok in RES if not ok]
    say(f"\nSUMMARY: {len(RES)} checks, {len(RES) - len(nf)} True, {len(nf)} False")
    for lab in nf:
        say(f"  False: {lab}")
    say(f"\nruntime {time.time() - T0:.0f} s")


if __name__ == '__main__':
    main()
