"""Referee cross-check of GU round 5 (formalization/geometric-unity/I-elliptic-reduction.md), independent construction.

Nothing from the repository is imported.  Own Clifford algebra (bitmask monomials E_K, own sign rule
E_A E_B = (-1)^{#{(i,j): i in A, j in B, i > j}} prod_{i in A&B} eta_i E_{A^B}), own Shiab operators built from the
definitions of D.1 / C (R-S):
    W_s xi = [Phi11 ^ xi]_s,  Phi11 = *Phi^3,  Phi^3 = sum e^abc (x) i gamma_abc  (Gamma-partner: i gamma_abc Gamma),
    [a ^ b]_C = ab - ba,  [a ^ b]_A = i (ab + ba)   (the draft's i^((1 +- 1)/2), C-gauge-shiab.md R-S),
    W_C o Ad_Gamma: precompose with Ad_Gamma = grade involution.
Self-adjoint form (D2): Psi(alpha ^ zeta) = b(alpha ^ Shiab zeta) = sum_a sum_{b<c} Tr(alpha_a M_abc(zeta_bc)),
    M_abc(Y) = [i g_abc, Y]_s (with Ad_Gamma for the twisted pair),  g_abc = gamma^a gamma^b gamma^c,  Tr 1 = 1.
So W_A = -{g, .}, W'_A = -{g Gamma, .}, W_C Ad_G = i [g, Ad_G .], W'_C Ad_G = i [g Gamma, Ad_G .]  (real templates
below carry these phases separately).  The 48 basis types O_(j,l), E_(j,l) of H3 are rebuilt from their definition
Tr((Gamma) q(omega) (X_j Y_l + Y_j X_l)).  Inner product <alpha, beta> = sum_a eta_aa Tr(alpha_a beta_a) = p(alpha, *beta).
True field equation (D5/G2): Upsilon = Shiab F_0 + Shiab(d_0 T) + kappa *T + Q(T),
    p(alpha, Q(T)) = (1/3) Psi(alpha ^ T ^ T) + (1/3) Psi(T ^ [T ^ alpha]),  [T ^ alpha]_bc = [T_b, a_c] - [T_c, a_b].
Exact integer / rational / Gaussian arithmetic.  Seeds fixed.  Single core.
"""
import itertools
import random
import sys
import time
from fractions import Fraction
from functools import lru_cache

import sympy as sp

T0 = time.time()
N = 14
FULL = (1 << N) - 1
EUC = (1,) * 14
E77 = (1,) * 7 + (-1,) * 7
SIG = {EUC: '(14,0)', E77: '(7,7)'}


def P(*a):
    print(*a)
    sys.stdout.flush()


def pc(x):
    return bin(x).count('1')


# ----------------------------------------------------------------------------------------------- Clifford algebra
class Cl:
    def __init__(self, eta):
        self.eta = eta

        @lru_cache(maxsize=1 << 21)
        def sign(A, B):
            s = 0
            b = B
            while b:
                low = b & -b
                j = low.bit_length() - 1
                s += pc(A >> (j + 1))
                b ^= low
            v = -1 if s & 1 else 1
            c = A & B
            while c:
                low = c & -c
                v *= eta[low.bit_length() - 1]
                c ^= low
            return v
        self.sign = sign

    def mul(self, X, Y):
        out = {}
        sg = self.sign
        for A, a in X.items():
            for B, b in Y.items():
                k = A ^ B
                out[k] = out.get(k, 0) + sg(A, B) * a * b
        return {k: v for k, v in out.items() if v != 0}

    def tr(self, X, Y):              # normalised Tr(XY)
        return sum(a * Y[K] * self.sign(K, K) for K, a in X.items() if K in Y)

    def gam_up(self, a):             # gamma^a = eta_a gamma_a
        return {1 << a: self.eta[a]}

    def g3(self, a, b, c):           # gamma^a gamma^b gamma^c
        return self.mul(self.mul(self.gam_up(a), self.gam_up(b)), self.gam_up(c))


def add(X, Y, s=1):
    out = dict(X)
    for k, v in Y.items():
        out[k] = out.get(k, 0) + s * v
    return {k: v for k, v in out.items() if v != 0}


def acc(D, Y, s=1):
    for k, v in Y.items():
        D[k] = D.get(k, 0) + s * v
        if D[k] == 0:
            del D[k]


def scl(X, c):
    return {k: c * v for k, v in X.items() if c * v != 0}


def hat(X):                          # Ad_Gamma = grade involution (d = 14 even)
    return {k: (-v if pc(k) & 1 else v) for k, v in X.items()}


def grade_part(X, pred):
    return {k: v for k, v in X.items() if pred(pc(k))}


def form_zero(D):
    return not any(D)


# ----------------------------------------------------------------------------------------------- Shiab templates
OPS = ('W_A', "W'_A", 'W_C.AdG', "W'_C.AdG")
PHASE = {'W_A': -1, "W'_A": -1, 'W_C.AdG': sp.I, "W'_C.AdG": sp.I}   # real template -> actual operator


class Shiab:
    def __init__(self, eta):
        self.eta = eta
        self.cl = Cl(eta)
        self.G = {FULL: 1}
        self.g = {}
        for a, b, c in itertools.permutations(range(N), 3):
            self.g[(a, b, c)] = self.cl.g3(a, b, c)
        self.gG = {k: self.cl.mul(v, self.G) for k, v in self.g.items()}

    def M(self, op, abc, Y):
        """real template M_abc(Y) (operator = PHASE[op] * template)"""
        cl = self.cl
        if op == 'W_A':
            g = self.g[abc]
            return add(cl.mul(g, Y), cl.mul(Y, g))
        if op == "W'_A":
            g = self.gG[abc]
            return add(cl.mul(g, Y), cl.mul(Y, g))
        if op == 'W_C.AdG':
            g, Yh = self.g[abc], hat(Y)
            return add(cl.mul(g, Yh), cl.mul(Yh, g), -1)
        if op == "W'_C.AdG":
            g, Yh = self.gG[abc], hat(Y)
            return add(cl.mul(g, Yh), cl.mul(Yh, g), -1)
        if op == 'W_C':                                           # untwisted commutator (skew), for the lemma
            g = self.g[abc]
            return add(cl.mul(g, Y), cl.mul(Y, g), -1)
        raise ValueError(op)

    def Mt(self, op, abc, X):
        """transpose w.r.t. Tr: Tr(X M(Y)) = Tr(Mt(X) Y), computed from the definition (not assumed = M)"""
        cl = self.cl
        if op in ('W_A', "W'_A"):
            g = self.g[abc] if op == 'W_A' else self.gG[abc]
            return add(cl.mul(X, g), cl.mul(g, X))
        if op in ('W_C.AdG', "W'_C.AdG"):
            g = self.g[abc] if op == 'W_C.AdG' else self.gG[abc]
            return hat(add(cl.mul(X, g), cl.mul(g, X), -1))     # Tr(X[g,Yh]) = Tr((Xg - gX) Yh) = Tr(hat(Xg - gX) Y)
        if op == 'W_C':
            g = self.g[abc]
            return add(cl.mul(X, g), cl.mul(g, X), -1)
        raise ValueError(op)

    # ---- dual of zeta -> Psi(. ^ zeta): list D_a with Psi(alpha ^ zeta) = sum_a Tr(alpha_a D_a)
    def dual(self, op, zeta):
        D = [dict() for _ in range(N)]
        for (b, c), Y in zeta.items():
            if not Y:
                continue
            for a in range(N):
                if a in (b, c):
                    continue
                acc(D[a], self.M(op, (a, b, c), Y))
        return D

    def psi(self, op, alpha, zeta):
        D = self.dual(op, zeta)
        return sum(self.cl.tr(alpha[a], D[a]) for a in range(N))

    def wedge(self, S, T):            # (S ^ T)_bc = S_b T_c - S_c T_b, b < c
        out = {}
        for b, c in itertools.combinations(range(N), 2):
            v = add(self.cl.mul(S[b], T[c]), self.cl.mul(S[c], T[b]), -1)
            if v:
                out[(b, c)] = v
        return out

    def cubic(self, op, T):           # C(T) = Psi(T ^ T ^ T)
        return self.psi(op, T, self.wedge(T, T))

    def Q3(self, op, T):
        """dual of 3 Q(T) = dual of alpha -> Psi(alpha ^ T ^ T) + Psi(T ^ [T ^ alpha])"""
        cl = self.cl
        D = self.dual(op, self.wedge(T, T))
        for b, c in itertools.combinations(range(N), 2):
            for a in range(N):
                if a in (b, c) or not T[a]:
                    continue
                Pm = self.Mt(op, (a, b, c), T[a])
                if not Pm:
                    continue
                if T[b]:
                    acc(D[c], add(cl.mul(Pm, T[b]), cl.mul(T[b], Pm), -1))
                if T[c]:
                    acc(D[b], add(cl.mul(Pm, T[c]), cl.mul(T[c], Pm), -1), -1)
        return D

    def pair(self, alpha, D):
        return sum(self.cl.tr(alpha[a], D[a]) for a in range(N))


# ---- the 48 basis types of H3, from their definition; all types at once
class Types(Shiab):
    def typed_M(self, abc, Y):
        """{t: M_t,abc(Y)} with Tr(X M_t(Y)) = O/E_t(abc; X, Y), O_(j,l) = Tr(g (X_j Y_l + Y_j X_l)),
        E_(j,l) = Tr(Gamma g (X_j Y_l + Y_j X_l))  (j < l); each monomial of Y feeds one O- and one E-type."""
        cl = self.cl
        g = self.g[abc]
        (A, sa), = g.items()
        Gg = cl.mul(self.G, g)
        (B, sb), = Gg.items()
        out = {}
        for K, y in Y.items():
            k = pc(K)
            for kind, Mm, s in (('O', A, sa), ('E', B, sb)):
                m = pc(K ^ Mm)
                if m < k:                    # X_j with j = m, Y_l with l = k:  P_j(Y_l g)
                    t, val = (kind, m, k), cl.sign(K, Mm) * s * y
                else:                        # X_l with l = m, Y_j with j = k:  P_l(g Y_j)
                    t, val = (kind, k, m), cl.sign(Mm, K) * s * y
                d = out.setdefault(t, {})
                d[K ^ Mm] = d.get(K ^ Mm, 0) + val
        return out

    def typed_dual(self, zeta):
        D = {}
        for (b, c), Y in zeta.items():
            for a in range(N):
                if a in (b, c):
                    continue
                for t, v in self.typed_M((a, b, c), Y).items():
                    acc(D.setdefault(t, [dict() for _ in range(N)])[a], v)
        return D

    def typed_Q3(self, T):
        cl = self.cl
        D = self.typed_dual(self.wedge(T, T))
        for b, c in itertools.combinations(range(N), 2):
            for a in range(N):
                if a in (b, c) or not T[a]:
                    continue
                for t, Pm in self.typed_M((a, b, c), T[a]).items():   # types are Tr-symmetric (checked below)
                    Pm = {k: v for k, v in Pm.items() if v}
                    if not Pm:
                        continue
                    Dt = D.setdefault(t, [dict() for _ in range(N)])
                    if T[b]:
                        acc(Dt[c], add(cl.mul(Pm, T[b]), cl.mul(T[b], Pm), -1))
                    if T[c]:
                        acc(Dt[b], add(cl.mul(Pm, T[c]), cl.mul(T[c], Pm), -1), -1)
        return D


# ----------------------------------------------------------------------------------------------- random data
def rmask(grades):
    return sum(1 << i for i in random.sample(range(N), random.choice(grades)))


def relt(grades, nterm, r=2):
    X = {}
    for _ in range(nterm):
        K = rmask(grades)
        X[K] = X.get(K, 0) + random.randint(-r, r)
    return {k: v for k, v in X.items() if v}


ODD = [1, 3, 5, 7, 9, 11, 13]
EVEN = [0, 2, 4, 6, 8, 10, 12, 14]
ALL = list(range(15))


def rsym(r=2, dens=0.4):
    M = [[0] * N for _ in range(N)]
    for a in range(N):
        for b in range(a, N):
            if random.random() < dens:
                M[a][b] = M[b][a] = random.randint(-r, r)
    return M


def KN(h, k):
    R = {}
    for a, b, c, d in itertools.product(range(N), repeat=4):
        v = h[a][c] * k[b][d] + h[b][d] * k[a][c] - h[a][d] * k[b][c] - h[b][c] * k[a][d]
        if v:
            R[(a, b, c, d)] = v
    return R


def einstein(R, eta):
    Ric = [[Fraction(0)] * N for _ in range(N)]
    for (a, b, c, d), v in R.items():
        if a == c:
            Ric[b][d] += eta[a] * v
    s = sum(eta[b] * Ric[b][b] for b in range(N))
    return [[Ric[a][b] - (s / 2 if a == b else 0) * eta[a] for b in range(N)] for a in range(N)], s


def curv_form(S, R):                 # 4 F_R = sum_{a<b} e^ab (x) sum_{c,d} R_abcd gamma^c gamma^d
    z = {}
    for (a, b, c, d), v in R.items():
        if a < b and c != d:
            acc(z.setdefault((a, b), {}), scl(S.cl.mul(S.cl.gam_up(c), S.cl.gam_up(d)), v))
    return {k: v for k, v in z.items() if v}


def einstein_dual(S, G, right=None, fac=1):
    """dual of -*(G_ab e^a (x) gamma^b [right]):  D_a = -eta_a sum_b G_ab gamma^b [right]"""
    D = []
    for a in range(N):
        X = {}
        for b in range(N):
            if G[a][b]:
                Y = S.cl.gam_up(b) if right is None else S.cl.mul(S.cl.gam_up(b), right)
                acc(X, scl(Y, -S.eta[a] * G[a][b] * fac))
        D.append(X)
    return D


def eqD(D1, D2):
    return all(add(D1[a], D2[a], -1) == {} for a in range(N))


# ======================================================================================================== checks
def main():
    random.seed(20261030)
    res = {}
    for eta in (EUC, E77):
        lab = SIG[eta]
        S = Shiab(eta)
        Ty = Types(eta)
        cl = S.cl
        P(f"\n######## {lab} ########")
        # ---- X0: algebra sanity: Gamma^2, Ad_Gamma = grade involution, gamma_a^2 = eta_a
        G2 = cl.mul(S.G, S.G)
        Gi = scl(S.G, G2[0]) if 0 in G2 else None
        okAd = all(cl.mul(cl.mul(S.G, {K: 1}), Gi) == hat({K: 1}) for K in (rmask(ALL) for _ in range(40)))
        P(f"X0  Gamma^2 = {G2[0]}; Gamma E_K Gamma^-1 = (-1)^|K| E_K on 40 random monomials: {okAd}")
        assert okAd
        # ---- X1: self-adjointness of the four templates (Tr(X M(Y)) symmetric), skewness of untwisted W_C
        ok_sym, ok_skew = True, True
        for _ in range(6):
            abc = tuple(random.sample(range(N), 3))
            X, Y = relt(ALL, 6), relt(ALL, 6)
            for op in OPS:
                ok_sym &= cl.tr(X, S.M(op, abc, Y)) == cl.tr(Y, S.M(op, abc, X))
                ok_sym &= add(S.Mt(op, abc, X), S.M(op, abc, X), -1) == {}
            ok_skew &= cl.tr(X, S.M('W_C', abc, Y)) == -cl.tr(Y, S.M('W_C', abc, X))
        P(f"X1  four templates Tr-symmetric (Mt = M) on random data: {ok_sym}; untwisted W_C Tr-skew: {ok_skew}")
        assert ok_sym and ok_skew
        # ---- X2: Einstein normalisation (D3.4) on a random algebraic curvature tensor; twisted pair kills F_R
        scal = 0
        while scal == 0:                              # Scal != 0 so that the trace part of G is tested too
            h1, k1, h2, k2 = rsym(), rsym(), rsym(), rsym()
            R = KN(h1, k1)
            for key, v in KN(h2, k2).items():
                R[key] = R.get(key, 0) + v
            R = {k: v for k, v in R.items() if v}
            G, scal = einstein(R, eta)
        F4 = curv_form(S, R)
        oks = []
        for op in OPS:
            Dop = [scl(x, 1) for x in S.dual(op, F4)]
            if op == 'W_A':
                oks.append(eqD([scl(x, -1) for x in Dop], einstein_dual(S, G, fac=4)))
            elif op == "W'_A":
                oks.append(eqD([scl(x, -1) for x in Dop], einstein_dual(S, G, S.G, fac=4)))
            else:
                oks.append(form_zero(Dop))
        P(f"X2  random algebraic R (Scal = {scal}): W_A F_R = -*(G e^a(x)gamma^b): {oks[0]};  W'_A F_R = -*(G e^a(x)"
          f"gamma^b Gamma): {oks[1]};  W_C.AdG F_R = 0: {oks[2]};  W'_C.AdG F_R = 0: {oks[3]}   [universal form factor"
          f" K = 1 in the convention Psi(alpha^zeta) = sum Tr(alpha_a [i g_abc, zeta_bc]_s)]")
        assert all(oks)
        # ---- X3: coordinates of the four operators in the 48-type basis (each monomial triple meets one type)
        coords = {op: {} for op in OPS}
        seen = set()
        okc = True
        for it in range(6000):
            a, b, c = random.sample(range(N), 3)
            J = rmask(ALL)
            L = J ^ (1 << a | 1 << b | 1 << c) ^ (FULL if it % 2 else 0)     # witness of an O- (even it) or E-type
            tv = Ty.typed_M((a, b, c), {L: 1})
            vals = {t: Ty.cl.tr({J: 1}, v) for t, v in tv.items()}
            vals = {t: v for t, v in vals.items() if v}
            if not vals:
                continue
            assert len(vals) == 1
            (t, vt), = vals.items()
            seen.add(t)
            for op in OPS:
                vo = cl.tr({J: 1}, S.M(op, (a, b, c), {L: 1})) * PHASE[op]
                r = sp.nsimplify(vo) / vt
                if t in coords[op]:
                    okc &= coords[op][t] == r
                else:
                    coords[op][t] = r
        nz = {op: {t: v for t, v in coords[op].items() if v != 0} for op in OPS}
        P(f"X3  {len(seen)} basis types met by 6000 random monomial triples; ratios consistent: {okc}")
        for op in OPS:
            P(f"     {op:9s}: " + ", ".join(f"{k}{j},{l}: {v}" for (k, j, l), v in sorted(nz[op].items(),
                                                                                         key=lambda x: (x[0][0], x[0][1]))))
        assert okc and len(seen) == 48
        # ---- X3b: real forms from the coordinates.  Basis of the real form: X_K = i^eps(|K|) gamma_K with
        # (14,0), u(128) (gamma_a Hermitian): gamma_K anti-Hermitian iff |K| = 2,3 mod 4  -> eps(k) = 1 for k = 0,1 mod 4;
        # (7,7), u(64,64) with gamma_a in u(64,64) (I0b): gamma_K^# = (-1)^(k(k+1)/2) gamma_K  -> eps(k) = 1 for k = 0,3 mod 4.
        # Psi = sum p_t Psi_t is real on the real form iff p_t in i^(eps_j + eps_l) R (Psi_t real on monomials).
        if eta == EUC:
            eps = lambda k: 1 if k % 4 in (0, 1) else 0       # noqa: E731
            scale_, nm_ = -sp.I, "-i x (code normalisation), i.e. Phi^E"
        else:
            eps = lambda k: 1 if k % 4 in (0, 3) else 0       # noqa: E731
            scale_, nm_ = 1, "code normalisation"
        okr = True
        for op in OPS:
            for (k_, j, l), v in nz[op].items():
                w = sp.expand(scale_ * v / sp.I ** (eps(j) + eps(l)))
                okr &= sp.im(w) == 0
        P(f"X3b real form {'u(128)' if eta == EUC else 'u(64,64)'}: all four operators are real for {nm_} "
          f"(p_t in i^(eps_j+eps_l) R on every block): {okr}")
        assert okr
        # ---- X4: gradient identity for my Q (3Q = dC), dense random T
        T = [relt(ALL, 5) for _ in range(N)]
        al = [relt(ALL, 5) for _ in range(N)]
        okg = True
        for op in OPS:
            Cp = S.cubic(op, [add(T[i], al[i]) for i in range(N)])
            Cm = S.cubic(op, [add(T[i], al[i], -1) for i in range(N)])
            dC = S.pair(al, S.Q3(op, T))
            okg &= Cp - Cm == 2 * dC + 2 * S.cubic(op, al)
        P(f"X4  C(T+a) - C(T-a) = 2 dC_T(a) + 2 C(a) with dC_T(a) = <a, 3Q(T)> (all four, dense random T): {okg}")
        assert okg
        # ---- X5: parity split (I2) and the vanishing of B(T_e,T_o)_odd, all 48 types and the four operators
        Te = [relt(EVEN, 4) for _ in range(N)]
        To = [relt(ODD, 4) for _ in range(N)]
        Tt = [add(Te[i], To[i]) for i in range(N)]
        qe, qo, qt = Ty.typed_Q3(Te), Ty.typed_Q3(To), Ty.typed_Q3(Tt)
        zero = [dict() for _ in range(N)]
        ok_odd = all(form_zero([grade_part(x, lambda g: g % 2 == 0) for x in q.get(t, zero)])
                     for q in (qe, qo) for t in seen)
        Bt = {t: [add(add(qt.get(t, zero)[a], qe.get(t, zero)[a], -1), qo.get(t, zero)[a], -1) for a in range(N)]
              for t in seen}                                   # 2 B_t(T_e, T_o) (times 3)
        B_odd_zero = all(form_zero([grade_part(x, lambda g: g % 2 == 1) for x in Bt[t]]) for t in seen)
        B_ev_nz = sum(1 for t in seen if not form_zero(Bt[t]))
        P(f"X5  all 48 types: Q_t(T_e), Q_t(T_o) have zero g_ev part: {ok_odd};  2B_t(T_e,T_o) has ZERO g_odd part for "
          f"every type: {B_odd_zero};  types with 2B_t(T_e,T_o) != 0 (necessarily g_ev): {B_ev_nz}")
        assert ok_odd and B_odd_zero
        # ---- X6: twist identity (I3) and the untwisted lemma
        P("X6  twist identity:")
        lemma = S.pair(Tt, S.dual('W_C', S.wedge(Tt, Tt)))
        P(f"     untwisted lemma p(T, W_C(T^T)) = {lemma} (random mixed T)")
        assert lemma == 0
        for op in ('W_C.AdG', "W'_C.AdG"):
            for _try in range(10):                    # draw until the identity is tested at a nonzero value
                Te2 = [relt(EVEN, 4) for _ in range(N)]
                To2 = [relt(ODD, 4) for _ in range(N)]
                Tt2 = [add(Te2[i], To2[i]) for i in range(N)]
                lhs = S.cubic(op, Tt2)
                rhs = 2 * S.psi(op, To2, S.wedge(Te2, Te2))
                if lhs != 0 or rhs != 0:
                    break
            Qo = S.Q3(op, To2)
            Qe = S.Q3(op, Te2)
            Wee = S.dual(op, S.wedge(Te2, Te2))
            Woo = S.dual(op, S.wedge(To2, To2))
            ok_e = eqD(Qe, [scl(x, 2) for x in Wee])
            P(f"     {op:9s}: C(T) = {lhs}, 2 Psi(T_o^T_e^T_e) = {rhs}, equal: {lhs == rhs};  Q(T_o) = 0: {form_zero(Qo)};"
              f"  Q(T_e) = (2/3) Shiab(T_e^T_e): {ok_e};  Shiab(T_o^T_o) != 0 (survives in printed (9.7)): {not form_zero(Woo)}")
            assert lhs == rhs and lhs != 0 and form_zero(Qo) and ok_e and not form_zero(Woo)
        # ---- X7: subspace {Psi in S_sa : Q(T) = 0 for all odd T}: rank mod p of the sampled system
        prime = 2 ** 61 - 1
        cols = {t: {} for t in seen}
        for smp in range(3):
            Ts = [relt(ODD, 3) for _ in range(N)]
            for t, D in Ty.typed_Q3(Ts).items():
                for a in range(N):
                    for k, v in D[a].items():
                        cols[t][(smp, a, k)] = v % prime
        tl = sorted(seen)
        vecs = [dict(cols[t]) for t in tl]
        rank, piv = 0, []
        basis = []
        for v in vecs:
            v = dict(v)
            for (pk, bv) in basis:
                if pk in v and v[pk]:
                    f = v[pk] * pow(bv[pk], -1, prime) % prime
                    for kk, x in bv.items():
                        v[kk] = (v.get(kk, 0) - f * x) % prime
                    v = {kk: x for kk, x in v.items() if x}
            if v:
                pk = next(iter(v))
                basis.append((pk, v))
        rank = len(basis)
        nrows = len({r for c in cols.values() for r in c})
        named = {}
        for nm, coeff in (('O0,3', {('O', 0, 3): 1}), ('E0,11', {('E', 0, 11): 1}),
                          ('W_C.AdG', nz['W_C.AdG']), ("W'_C.AdG", nz["W'_C.AdG"])):
            Ts = [relt(ODD, 4) for _ in range(N)]
            Dq = Ty.typed_Q3(Ts)
            tot = [dict() for _ in range(N)]
            for t, cf in coeff.items():
                for a in range(N):
                    acc(tot[a], scl(Dq.get(t, zero)[a], cf))
            named[nm] = form_zero(tot)
        P(f"X7  sampled system ({nrows} rows x 48): rank mod (2^61-1) = {rank} => dim{{Psi: Q(T_o) = 0 for all odd T}}"
          f" <= {48 - rank} (rank over Q >= rank mod p);  the four named elements give Q(T_o) = 0 on a fresh sample:"
          f" {named}")
        res[lab + 'rank'] = rank
        # ---- X8: invariant sector T = t Phi^1 + u Phi^1 Gamma, exact c1, c13 per operator
        Phi1 = [{1 << a: 1} for a in range(N)]
        Phi1G = [cl.mul({1 << a: 1}, S.G) for a in range(N)]
        tt, uu, kap = sp.symbols('t u kappa')
        xs, ys, zs, ws = sp.symbols('x y z w')
        FAM = dict(zip(OPS, (xs, ys, zs, ws)))
        c1, c13 = kap * tt, kap * uu
        okinv = True
        for op in OPS:
            data = {}
            for (tv, uv) in ((1, 0), (0, 1), (1, 1), (2, 1), (1, -2)):
                Tv = [add(scl(Phi1[a], tv), scl(Phi1G[a], uv)) for a in range(N)]
                D = S.Q3(op, Tv)                                    # dual of 3Q; Q = *beta with beta_a = eta_a D_a
                # beta_a = q1 gamma_a + q13 gamma_a Gamma for every a
                (kG, sG), = Phi1G[0].items()
                q1 = Fraction(eta[0] * D[0].get(1, 0), 3)
                q13 = Fraction(eta[0] * D[0].get(kG, 0) * sG, 3) if sG else 0
                for a in range(N):
                    exp = add(scl(Phi1[a], q1 * eta[a]), scl(Phi1G[a], q13 * eta[a]))
                    okinv &= add({k: Fraction(v, 3) for k, v in D[a].items()}, exp, -1) == {}
                data[(tv, uv)] = (q1, q13)
            A1, A2, A3, B1, B2, B3 = sp.symbols('A1:4 B1:4')
            e1 = [sp.Eq(A1 * tv ** 2 + A2 * tv * uv + A3 * uv ** 2, sp.Rational(v[0].numerator, v[0].denominator))
                  for (tv, uv), v in data.items()]
            e2 = [sp.Eq(B1 * tv ** 2 + B2 * tv * uv + B3 * uv ** 2, sp.Rational(v[1].numerator, v[1].denominator))
                  for (tv, uv), v in data.items()]
            s1, s2 = sp.solve(e1, [A1, A2, A3], dict=True), sp.solve(e2, [B1, B2, B3], dict=True)
            assert len(s1) == 1 and len(s2) == 1
            ph = PHASE[op]
            c1 += FAM[op] * ph * (s1[0][A1] * tt ** 2 + s1[0][A2] * tt * uu + s1[0][A3] * uu ** 2)
            c13 += FAM[op] * ph * (s2[0][B1] * tt ** 2 + s2[0][B2] * tt * uu + s2[0][B3] * uu ** 2)
        c1, c13 = sp.expand(c1), sp.expand(c13)
        sg = 1 if eta == EUC else -1
        c1_note = kap * tt + xs * (312 * tt ** 2 + sg * 104 * uu ** 2) + sg * 208 * ys * tt * uu
        c13_note = kap * uu + 208 * xs * tt * uu + ys * (104 * tt ** 2 + sg * 312 * uu ** 2)
        P(f"X8  invariant sector (every T-term of the invariant form, all a: {okinv}):  c1 = {c1};  c13 = {c13};"
          f"  equal to I6: {sp.expand(c1 - c1_note) == 0 and sp.expand(c13 - c13_note) == 0}")
        assert okinv and sp.expand(c1 - c1_note) == 0 and sp.expand(c13 - c13_note) == 0
        res[lab] = (c1, c13)
        # printed (9.7): Upsilon^pr = Shiab F_A + kappa *T, T-part kappa *T + Shiab(T^T) (no -h/3)
        p1, p13 = kap * tt, kap * uu
        for op in OPS:
            data = {}
            for (tv, uv) in ((1, 0), (0, 1), (1, 1), (2, 1), (1, -2)):
                Tv = [add(scl(Phi1[a], tv), scl(Phi1G[a], uv)) for a in range(N)]
                D = S.dual(op, S.wedge(Tv, Tv))
                (kG, sG), = Phi1G[0].items()
                data[(tv, uv)] = (sp.Integer(eta[0] * D[0].get(1, 0)), sp.Integer(eta[0] * D[0].get(kG, 0) * sG))
            A1, A2, A3, B1, B2, B3 = sp.symbols('A1:4 B1:4')
            s1 = sp.solve([sp.Eq(A1 * a_ ** 2 + A2 * a_ * b_ + A3 * b_ ** 2, v[0]) for (a_, b_), v in data.items()],
                          [A1, A2, A3], dict=True)
            s2 = sp.solve([sp.Eq(B1 * a_ ** 2 + B2 * a_ * b_ + B3 * b_ ** 2, v[1]) for (a_, b_), v in data.items()],
                          [B1, B2, B3], dict=True)
            assert len(s1) == 1 and len(s2) == 1
            p1 += FAM[op] * PHASE[op] * (s1[0][A1] * tt ** 2 + s1[0][A2] * tt * uu + s1[0][A3] * uu ** 2)
            p13 += FAM[op] * PHASE[op] * (s2[0][B1] * tt ** 2 + s2[0][B2] * tt * uu + s2[0][B3] * uu ** 2)
        res[lab + 'printed'] = (sp.expand(p1), sp.expand(p13))
        P(f"X8p printed (9.7) in the same sector: c1^pr = {sp.expand(p1)};  c13^pr = {sp.expand(p13)}")
        # ---- X9: constraint system on the Lambda^1 and [2,1]_0 constituents at v = e^0 (= symbol of Shiab o d)
        def zeta_from(alpha):          # v ^ alpha with v = e^0: (e^0 ^ alpha)_{0c} = alpha_c
            return {(0, c): alpha[c] for c in range(1, N) if alpha[c]}
        tests = {'Lambda^1': (Phi1, Phi1G),
                 '[2,1]_0': ([{1 << 2: 1} if a == 2 else ({1 << 3: -1} if a == 3 else {}) for a in range(N)],
                             [cl.mul({1 << 2: 1}, S.G) if a == 2 else (cl.mul({1 << 3: -1}, S.G) if a == 3 else {})
                              for a in range(N)])}
        for nm, (inA, inB) in tests.items():
            outs = {}
            for col, alpha in enumerate((inA, inB)):
                z = zeta_from(alpha)
                for op in OPS:
                    D = S.dual(op, z)
                    for a in range(N):
                        for k, v in D[a].items():
                            outs.setdefault(pc(k), []).append((op, col, a, k, v))
            grades = sorted(outs)
            rows = []
            for g in grades:
                # group by (op, col): each must be proportional to one reference pattern in grade g
                pats = {}
                for op, col, a, k, v in outs[g]:
                    pats.setdefault((op, col), {})[(a, k)] = v
                ref = next(iter(pats.values()))
                k0 = next(iter(ref))
                row = [0, 0]
                for (op, col), pat in pats.items():
                    r = Fraction(pat.get(k0, 0), ref[k0])
                    assert set(pat) == set(ref) and all(Fraction(pat[q]) == r * ref[q] for q in ref)
                    row[col] += FAM[op] * PHASE[op] * sp.Rational(r.numerator, r.denominator)
                rows.append(row)
            Mx = sp.Matrix(rows)
            det = sp.factor(Mx.det()) if Mx.shape == (2, 2) else None
            P(f"X9  constraint Shiab(v ^ (a A + b B)) on the {nm} constituent (v = e^0): output grades {grades};  "
              f"det = {det}")
            res[(lab, nm)] = det
        # Lambda^3 constituent: zeta_ab = sum_c eps_abc gamma^c on {0,1,2} (and its Gamma-partner)
        L3 = {(0, 1): cl.gam_up(2), (0, 2): scl(cl.gam_up(1), -1), (1, 2): cl.gam_up(0)}
        L3G = {k: cl.mul(v, S.G) for k, v in L3.items()}
        outs = {}
        for col, z in enumerate((L3, L3G)):
            for op in OPS:
                D = S.dual(op, z)
                for a in range(N):
                    for k, v in D[a].items():
                        outs.setdefault(pc(k), {}).setdefault((op, col), {})[(a, k)] = v
        rows = []
        for g in sorted(outs):
            pats = outs[g]
            ref = next(iter(pats.values()))
            k0 = next(iter(ref))
            row = [0, 0]
            for (op, col), pat in pats.items():
                r = Fraction(pat.get(k0, 0), ref[k0])
                assert set(pat) == set(ref) and all(Fraction(pat[q]) == r * ref[q] for q in ref)
                row[col] += FAM[op] * PHASE[op] * sp.Rational(r.numerator, r.denominator)
            rows.append(row)
        Mx = sp.Matrix(rows)
        minors = sorted({sp.factor(Mx[[i, j], :].det()) for i, j in itertools.combinations(range(Mx.shape[0]), 2)}
                        - {0}, key=str)
        P(f"X9  Lambda^3 constituent: output grades {sorted(outs)}; nonzero 2x2 minors (up to a common constant): {minors}")
        # ---- X10: symmetric soldering sector (I5(a)) exact on two random symmetric theta
        ok5 = True
        for _ in range(2):
            th = rsym(r=3, dens=0.5)
            Tth = []
            for a in range(N):
                X = {}
                for b in range(N):
                    if th[a][b]:
                        acc(X, scl(cl.gam_up(b), th[a][b]))
                Tth.append(X)
            Rth = {}
            for a, b, c, d in itertools.product(range(N), repeat=4):
                v = 4 * (th[a][c] * th[b][d] - th[a][d] * th[b][c])
                if v:
                    Rth[(a, b, c, d)] = v
            Gth, _ = einstein(Rth, eta)
            DA = [scl(x, Fraction(-1, 3)) for x in S.Q3('W_A', Tth)]          # Q_{W_A} dual (phase -1)
            DpA = [scl(x, Fraction(-1, 3)) for x in S.Q3("W'_A", Tth)]
            ok5 &= eqD(DA, einstein_dual(S, Gth))                             # Q_{W_A} = -*(G(R_theta) gamma)
            ok5 &= eqD(DpA, einstein_dual(S, Gth, S.G, fac=Fraction(1, 3)))   # Q_{W'_A} = -(1/3)*(G(R_theta) gamma Gamma)
            ok5 &= form_zero(S.Q3('W_C.AdG', Tth)) and form_zero(S.Q3("W'_C.AdG", Tth))
        P(f"X10 symmetric theta (2 samples): Q_W_A = -*(G(R_theta) e(x)gamma), Q_W'_A = -(1/3)*(G(R_theta) e(x)gamma Gamma),"
          f" twisted Q = 0, with R_theta = 4(theta_ac theta_bd - theta_ad theta_bc): {ok5}")
        assert ok5
    # ---- X11: reduction algebra (I7) from the X8 coefficients and the F closed forms
    P("\n######## reduction algebra (I7) ########")
    tt, uu, kap = sp.symbols('t u kappa')
    xs, ys, zs, ws = sp.symbols('x y z w')
    cH, cV = sp.symbols('c_H c_V', positive=True)
    s_, r_ = sp.symbols('s r', real=True)
    xr, yr = sp.symbols("x' y'", real=True)
    c1, c13 = res['(14,0)']
    sub = {xs: -sp.I * xr, ys: -sp.I * yr, tt: sp.I * s_, uu: sp.I * r_}
    e1, e13 = sp.expand(-sp.I * c1.subs(sub)), sp.expand(-sp.I * c13.subs(sub))
    sh1, sh13 = e1.subs({xr: 1, yr: 0}), e13.subs({xr: 1, yr: 0})
    P(f"(14,0) u(128)-real (x = -i x', t = i s, u = i r): x'G + e1 g = 0, y'G + e13 g = 0 with e1 = {e1}, e13 = {e13};"
      f" Shiab#: e1 = {sh1}, e13 = {sp.factor(sh13)}")
    # (P): Lambda_4 = Lambda_P + c_H e1 on e13 = 0; minimum over the two branches
    br_r0 = sh1.subs(r_, 0)
    smin = sp.solve(sp.diff(br_r0, s_), s_)[0]
    P(f"  (P) branch r = 0: e1 = {br_r0}, min at s = {smin}: {sp.simplify(br_r0.subs(s_, smin))};  branch s = -kappa/208:"
      f" e1 = {sp.expand(sh1.subs(s_, -kap / 208))}  => Lambda_4 in [Lambda_P - c_H kappa^2/1248, oo) for x' = 1")
    xg = sp.symbols('X', positive=True)
    gen = sp.expand((kap * s_ + xg * 312 * s_ ** 2) / xg)
    P(f"  general x' = X > 0, y' = 0: branch r = 0 minimum of e1/x' = {sp.simplify(gen.subs(s_, sp.solve(sp.diff(gen, s_), s_)[0]))}"
      f"  (so the stated bound c_H kappa^2/1248 is the x' = 1 value)")
    # (E): F6 case G = (15/(2 c_V)) g_aleph
    eqE = sp.expand(sh1 + sp.Rational(15, 2) / cV)
    disc = sp.discriminant(eqE.subs(r_, 0), s_)
    P(f"  (E) [F6: Z1, lambda = 1/5, Ric(g) = 0, G = 15/(2c_V) g_aleph]: {eqE} = 0 on {sp.factor(sh13)} = 0;"
      f" branch r = 0 discriminant = {sp.factor(disc)} (>= 0 iff c_V kappa^2 >= 9360); branch s = -kappa/208:"
      f" {sp.expand(eqE.subs(s_, -kap / 208))} > 0")
    lam = sp.symbols('lambda')
    LP_Z1 = sp.Rational(3, 2) * cH / cV * (4 - 15 * lam) / (1 - 4 * lam)          # F5 table, n = 4
    L4_E = sp.simplify(LP_Z1.subs(lam, sp.Rational(1, 5)) + cH * (-sp.Rational(15, 2) / cV))
    P(f"  configurations solving Upsilon = 0 AT s(X) in this sector (= (E); the only case Cor. I2.1 can reach there):"
      f" e1 = -15/(2c_V), Lambda_P(Z1, 1/5) = {LP_Z1.subs(lam, sp.Rational(1, 5))}, so their (P) value is Lambda_4 = {L4_E}"
      f" (fixed, = 0), not free")
    c1, c13 = res['(7,7)']
    sh1, sh13 = c1.subs({xs: 1, ys: 0}), c13.subs({xs: 1, ys: 0})
    P(f"(7,7): Shiab#: c1 = {sh1}, c13 = {sp.factor(sh13)};  branch u = 0: Lambda_4 - Lambda_P = -c_H({sh1.subs(uu, 0)})"
      f" <= c_H kappa^2/1248;  branch t = -kappa/208: -c_H({sp.expand(sh1.subs(tt, -kap / 208))}) >= -c_H kappa^2/416;"
      f"  union = R")
    q1, q13 = res['(7,7)printed']
    q1s, q13s = q1.subs({xs: 1, ys: 0}), q13.subs({xs: 1, ys: 0})
    P(f"(7,7) printed (9.7), Shiab# (= W_A here): c1^pr = {q1s}, c13^pr = {q13s};  for kappa != 0 the curve forces u = 0,"
      f" so Lambda_4 - Lambda_P = -c_H({q1s.subs(uu, 0)}) <= c_H kappa^2/1248: NOT every real value (the second branch"
      f" t = -kappa/208 comes from the -h/3 term, absent in the draft)")
    P(f"(7,7) members with x = y = 0 (no Einstein output): c1 = {c1.subs({xs: 0, ys: 0})}, c13 = "
      f"{c13.subs({xs: 0, ys: 0})}; Shiab F_0 = 0, so T = 0 solves Upsilon = 0 identically (the (E) 'no solution' needs (x,y) != 0)")
    P(f"\nX9 summary: (14,0) Lambda^1 det = {res[('(14,0)', 'Lambda^1')]}, [2,1]_0 det = {res[('(14,0)', '[2,1]_0')]};"
      f" (7,7) Lambda^1 det = {res[('(7,7)', 'Lambda^1')]}, [2,1]_0 det = {res[('(7,7)', '[2,1]_0')]}")
    P(f"\nruntime {time.time() - T0:.0f} s")


if __name__ == '__main__':
    main()
