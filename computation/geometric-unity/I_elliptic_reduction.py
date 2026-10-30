"""GU round 5 (I): what the elliptic twist W_C o Ad_Gamma does to Upsilon = 0, and to its reduction along a section.

Formalization: formalization/geometric-unity/I-elliptic-reduction.md.  Builds on D (D2: S_sa = {Shiab = Psi ^ .};
D5: true Upsilon), G (G1/G2: unitary gauge, Upsilon = Shiab F_A + kappa *T - h(T)/3), H (H3: 48-type basis of S_sa and
the coordinates of W_A, W'_A, W_C o Ad_Gamma, W'_C o Ad_Gamma; H7: the elliptic family), F (F3-F6: Ricci blocks of
g_aleph along s(X) in the readings Z1/Z2).  Imports G_symbols.py, H_ellipticity.py, shiab_checks.py (unmodified).

Setting.  d = 14, eta = (+^14) [Euclidean, where Shiab# is elliptic] and (+^7,-^7) [draft]; g = Cl(eta) (x) C (complex
form of u(128) / u(64,64)); Clifford monomials E_K; Gamma = E_full; Ad_Gamma = (-1)^grade.  Unitary gauge (G1):
T = S in Omega^1(g), A = nabla^0 + T, F_A = F_0 + d_0 T + T ^ T, T ^ T = {(a,b): T_a T_b - T_b T_a}.
For Shiab in S_sa the field equation (D5/G2) is the gradient of
   L(T) = Psi(T ^ F_0) + 1/2 Psi(T ^ d_0 T) + 1/3 Psi(T ^ T ^ T) + kappa/2 <T, T>,
   Upsilon(alpha) = Psi(alpha ^ F_0) + Psi(alpha ^ d_0 T) + Q(T)(alpha) + kappa <alpha, T>,
   Q(T)(alpha) = 1/3 Psi(alpha ^ T ^ T) + 1/3 Psi(T ^ [T ^ alpha])          ( = Shiab(T^T) - h(T)/3 of D5 )
with Psi(alpha ^ zeta) = b(alpha ^ Shiab zeta), b = Tr/128.  Upsilon is stored as its dual D: Upsilon(alpha) = sum_a
Tr(D_a alpha_a); then Upsilon = *(sum_a e^a (x) eta^aa D_a).  Family (H7, code normalisation of G_symbols, Phi^3 built
from i gamma_abc):  Shiab = x W_A + y W'_A + z W_C Ad_Gamma + w W'_C Ad_Gamma;  Shiab# = W_A + W_C Ad_Gamma.
Every quantity is computed per basis type t (48 types), hence for all of S_sa at once.
Sections
  I0  validation: per-type Shiab = operator code (G_symbols / H_ellipticity) on random inputs; dual of
      alpha -> Psi(T^[T^alpha]) = direct evaluation; gradient identity C(T+a) - C(T-a) = 2 dC_T(a) + 2 C(a);
      W_A(F_R) = -*(G e(x)gamma), W'_A(F_R) = -*(G e(x)gamma Gamma), twisted operators kill F_R; (7,7): the four
      operators (code normalisation) preserve u(64,64) (explicit matrices of clifford77.py).
  I1  grade selection rules of the four operators; Spin(14)-constituent table on input grades 0,1,2,12,13,14.
  I2  parity: Q_t(T_e) and Q_t(T_o) are odd and B_t(T_e,T_o) is even for all 48 types (=> Upsilon_ev = Shiab(d_0 T_o)
      + kappa *T_e + 2B(T_e,T_o)_ev; no cross term in Upsilon_odd).
  I3  twist identity C_CG(T) = 2 Psi_CG(T_o ^ T_e ^ T_e); Q_CG(T_o) = 0; Q_CG(T_e) = (2/3) W_CG(T_e^T_e);
      the subspace of S_sa with Q(T_o) = 0 for all odd T (exact null space).
  I4  soldering sector T = theta_ab e^a (x) gamma^b + theta'_ab e^a (x) gamma^b Gamma: exact Q for symmetric theta;
      constraint equations Shiab(d_0 T) = 0 per constituent: determinants.
  I5  Spin(14)-invariant sector T = t Phi^1 + u Phi^1 Gamma: exact Upsilon for the whole family.
  I6  reduction along s(X) (Z1, Z2; readings (P), (E)) in both real forms (sympy, from F3-F6 closed forms).
Exact arithmetic throughout (Python integers, Fractions, sympy rationals / Gaussian rationals).  Seeds fixed.
"""
import itertools
import os
import random
import sys
import time
from fractions import Fraction

for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')          # 128x128 products in I0b: threading overhead dominates
import numpy as np
import sympy as sp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import G_symbols as GS  # noqa: E402
import H_ellipticity as H  # noqa: E402
from G_symbols import Cliff, mask, N, FULL  # noqa: E402
from shiab_checks import random_curvature, kn  # noqa: E402

T0 = time.time()
random.seed(20210401)
pc = lambda v: bin(v).count('1')  # noqa: E731
EUC = [1] * 14
E77 = [1] * 7 + [-1] * 7
SIG = {1: '(14,0)', -1: '(7,7)'}
TYPES, IDX, NT, NAME = H.TYPES, H.IDX, H.NT, H.NAME
ODD = [1, 3, 5, 7, 9, 11, 13]
EVEN = [0, 2, 4, 6, 8, 10, 12, 14]


def P(*a):
    print(*a)
    sys.stdout.flush()


def perm_sign(seq):
    s = 1
    for i in range(len(seq)):
        for j in range(i + 1, len(seq)):
            if seq[i] > seq[j]:
                s = -s
    return s


def add(X, Y, s=1):
    out = dict(X)
    for k, v in Y.items():
        out[k] = out.get(k, 0) + s * v
    return {k: v for k, v in out.items() if v != 0}


def acc_into(acc, Y, s=1):
    for k, v in Y.items():
        acc[k] = acc.get(k, 0) + s * v
        if acc[k] == 0:
            del acc[k]


def scale(X, c):
    return {k: c * v for k, v in X.items() if c * v != 0}


class Alg:
    """Clifford elements as dicts mask -> coefficient; Psi duals per basis type."""

    def __init__(self, eta):
        self.eta = list(eta)
        self.cl = Cliff(eta)
        self.P = H.Psi(eta)
        self._ss = {}

    def ss(self, J):                      # E_J E_J = ss(J) * 1
        v = self._ss.get(J)
        if v is None:
            v = self._ss[J] = self.cl.sign(J, J)
        return v

    def mul(self, X, Y):
        out = {}
        sg = self.cl.sign
        for J, a in X.items():
            for K, b in Y.items():
                k = J ^ K
                out[k] = out.get(k, 0) + sg(J, K) * a * b
        return {k: v for k, v in out.items() if v != 0}

    def comm(self, X, Y):
        return add(self.mul(X, Y), self.mul(Y, X), -1)

    def tr(self, X, Y):                   # normalised trace Tr(XY), Tr 1 = 1
        return sum(a * Y[J] * self.ss(J) for J, a in X.items() if J in Y)

    def gamma(self, a):                   # gamma_a
        return {1 << a: 1}

    def gammaG(self, a):                  # gamma_a Gamma
        return {FULL ^ (1 << a): self.cl.sign(1 << a, FULL)}

    def Gam(self):
        return {FULL: 1}

    def psi_dual(self, A, Y):
        """{t: D_t}: Psi_t(e^A; X, Y) = Tr(D_t X) for all X (A a sorted 3-form mask)"""
        out = {}
        ev = self.P.eval
        for L, c in Y.items():
            for J in (A ^ L, A ^ L ^ FULL):
                for t, v in ev(A, J, L).items():
                    d = out.setdefault(t, {})
                    d[J] = d.get(J, 0) + c * v * self.ss(J)
        return {t: {k: v for k, v in d.items() if v != 0} for t, d in out.items()}

    def psi_form(self, zeta):
        """dual of beta -> Psi(beta ^ zeta) (zeta: {(a,b), a<b: elt}): {t: [D_0..D_13]}"""
        out = {}
        for (a, b), Y in zeta.items():
            if not Y:
                continue
            for c in range(N):
                if c in (a, b):
                    continue
                s = perm_sign((c, a, b))
                for t, D in self.psi_dual(mask((c, a, b)), Y).items():
                    slot = out.setdefault(t, [dict() for _ in range(N)])
                    acc_into(slot[c], D, s)
        return out

    def psi_htilde(self, T):
        """dual of alpha -> Psi(T ^ [T ^ alpha]), [T^alpha]_jk = [T_j, a_k] - [T_k, a_j]: {t: [D_a]}"""
        out = {}
        for i in range(N):
            if not T[i]:
                continue
            for j, k in itertools.combinations(range(N), 2):
                if i in (j, k):
                    continue
                s = perm_sign((i, j, k))
                for t, Dp in self.psi_dual(mask((i, j, k)), T[i]).items():
                    slot = out.setdefault(t, [dict() for _ in range(N)])
                    if T[j]:
                        acc_into(slot[k], self.comm(Dp, T[j]), s)
                    if T[k]:
                        acc_into(slot[j], self.comm(Dp, T[k]), -s)
        return out

    def wedge11(self, P1, Q1):
        out = {}
        for a, b in itertools.combinations(range(N), 2):
            v = add(self.mul(P1[a], Q1[b]), self.mul(P1[b], Q1[a]), -1)
            if v:
                out[(a, b)] = v
        return out

    def bracket11(self, P1, Q1):
        out = {}
        for a, b in itertools.combinations(range(N), 2):
            v = add(self.comm(P1[a], Q1[b]), self.comm(P1[b], Q1[a]), -1)
            if v:
                out[(a, b)] = v
        return out

    def pair_form(self, D, alpha):
        return sum(self.tr(D[a], alpha[a]) for a in range(N))

    def quad3(self, T):
        """per type: dual of 3 Q(T) = Psi(alpha^T^T) + Psi(T^[T^alpha]) = dC_T (C(T) = Psi(T^T^T))"""
        p1, p2 = self.psi_form(self.wedge11(T, T)), self.psi_htilde(T)
        z = [dict() for _ in range(N)]
        return {t: [add(p1.get(t, z)[a], p2.get(t, z)[a]) for a in range(N)] for t in set(p1) | set(p2)}


def combine(per_type, coeffs):
    out = [dict() for _ in range(N)]
    for t, c in coeffs.items():
        if c != 0 and t in per_type:
            for a in range(N):
                acc_into(out[a], per_type[t][a], c)
    return out


def grades_of(D):
    return sorted({pc(K) for Da in D for K in Da})


def parity_part(D, p):
    return [{k: v for k, v in Da.items() if pc(k) % 2 == p} for Da in D]


def _vec(d):
    return {IDX[(n[0], int(n[1:].split(',')[0]), int(n[1:].split(',')[1]))]: v for n, v in d.items()}


# coordinates from H_ellipticity.out (H1), stored as unit * integer sign pattern; validated in I0
OPS = {
    'W_A': (-2, _vec({'O0,3': 1, 'O1,2': 1, 'O2,5': 1, 'O3,4': 1, 'O4,7': 1, 'O5,6': 1, 'O6,9': 1, 'O7,8': 1,
                      'O8,11': 1, 'O9,10': 1, 'O10,13': 1, 'O11,12': 1})),
    "W'_A": (2, _vec({'E0,11': 1, 'E1,10': 1, 'E2,9': 1, 'E2,13': 1, 'E3,8': 1, 'E3,12': 1, 'E4,7': 1, 'E4,11': 1,
                      'E5,6': 1, 'E5,10': 1, 'E6,9': 1, 'E7,8': 1})),
    'W_C.AdG': (2j, _vec({'O1,4': -1, 'O2,3': 1, 'O3,6': -1, 'O4,5': 1, 'O5,8': -1, 'O6,7': 1, 'O7,10': -1,
                          'O8,9': 1, 'O9,12': -1, 'O10,11': 1, 'O11,14': -1, 'O12,13': 1})),
    "W'_C.AdG": (2j, _vec({'E1,12': 1, 'E2,11': -1, 'E3,10': 1, 'E3,14': 1, 'E4,9': -1, 'E4,13': -1, 'E5,8': 1,
                           'E5,12': 1, 'E6,7': -1, 'E6,11': -1, 'E7,10': 1, 'E8,9': -1})),
}
OPN = list(OPS)
xs, ys, zs, ws = sp.symbols('x y z w')
FAMC = dict(zip(OPN, (xs, ys, zs, ws)))


def unit_sp(u):
    return sp.Integer(int(u.real)) + sp.I * sp.Integer(int(u.imag))


def rand_elt(grades, nterm, r=2):
    X = {}
    for _ in range(nterm):
        K = mask(tuple(random.sample(range(N), random.choice(grades))))
        X[K] = X.get(K, 0) + random.randint(-r, r)
    return {k: v for k, v in X.items() if v}


def z3_form(Al, Z, right=None):
    """zeta_ab = sum_c Z[a][b][c] gamma^c (a<b) (optionally times `right` on the right)"""
    z = {}
    for a, b in itertools.combinations(range(N), 2):
        X = {}
        for c in range(N):
            v = int(Z[a, b, c])
            if v:
                X[1 << c] = X.get(1 << c, 0) + v * Al.eta[c]
        X = {k: v for k, v in X.items() if v}
        if X and right is not None:
            X = Al.mul(X, right)
        if X:
            z[(a, b)] = X
    return z


def z4_form(Al, Z, right=None, half=False):
    """zeta_ab = sum_{c<d} c0 Z[a][b][c][d] gamma^c gamma^d, c0 = 1/2 if half (spin curvature F_R) else 1"""
    z = {}
    for a, b in itertools.combinations(range(N), 2):
        X = {}
        for c, d in itertools.combinations(range(N), 2):
            v = Z[a][b][c][d]
            if v:
                k = mask((c, d))
                X[k] = X.get(k, 0) + (Fraction(int(v), 2) if half else int(v)) * Al.eta[c] * Al.eta[d] * \
                    Al.cl.sign(1 << c, 1 << d)
        X = {k: v for k, v in X.items() if v}
        if X and right is not None:
            X = Al.mul(X, right)
        if X:
            z[(a, b)] = X
    return z


def T_from(Al, th, thp=None):
    """T = theta_ab e^a (x) gamma^b + theta'_ab e^a (x) gamma^b Gamma"""
    T = []
    for a in range(N):
        X = {}
        for b in range(N):
            if th[a][b]:
                X = add(X, scale(Al.gamma(b), th[a][b] * Al.eta[b]))
            if thp is not None and thp[a][b]:
                X = add(X, scale(Al.gammaG(b), thp[a][b] * Al.eta[b]))
        T.append(X)
    return T


def rsym(r=2, dens=0.4):
    M = [[0] * N for _ in range(N)]
    for a in range(N):
        for b in range(a, N):
            if random.random() < dens:
                M[a][b] = M[b][a] = random.randint(-r, r)
    return M


def rasym(r=2, dens=0.3):
    M = [[0] * N for _ in range(N)]
    for a in range(N):
        for b in range(a + 1, N):
            if random.random() < dens:
                v = random.randint(-r, r)
                M[a][b], M[b][a] = v, -v
    return M


# ====================================================================== I0
def i0_validate():
    P("=== I0. Validation of the per-type evaluator (exact) ===")
    for eta in (EUC, E77):
        Al = Alg(eta)
        okall = True
        for nm, base, tw in (('W_A', 'WA', False), ("W'_A", 'WpA', False), ('W_C.AdG', 'WC', True),
                             ("W'_C.AdG", 'WpC', True)):
            cl, _, ops = H.shiab_ops_ext(eta, base)
            unit, sv = OPS[nm]
            ok = True
            for _ in range(25):
                a, b = random.sample(range(N), 2)
                K = mask(tuple(random.sample(range(N), random.randint(0, 14))))
                zeta = {(min(a, b), max(a, b)): {K: 1 if b < a else -1}}
                D = combine(Al.psi_form(zeta), sv)
                sgn = (-1) ** pc(K) if tw else 1
                for m in range(N):
                    ref = GS.apply_op(cl, ops[(b, a)][m], K) if m in ops[(b, a)] else {}
                    ref = {k: sgn * complex(v) for k, v in ref.items() if v != 0}
                    got = {k: complex(unit * (-1) ** m * v) for k, v in D[m].items()}
                    ok &= ref == got
            okall &= ok
            P(f"  {SIG[eta[7]]} {nm:9s}: Psi-basis operator = operator code on 25 random inputs e^b^e^a (x) E_K "
              f"(13-form coefficient X_m = (-1)^m D_m): {ok}")
        assert okall
        T = [rand_elt(list(range(15)), 12) for _ in range(N)]
        al = [rand_elt(list(range(15)), 12) for _ in range(N)]
        pz = Al.psi_form(Al.bracket11(T, al))
        direct = {t: Al.pair_form(pz[t], T) for t in pz}
        ph = Al.psi_htilde(T)
        via = {t: Al.pair_form(ph[t], al) for t in ph}
        ok1 = all(direct.get(t, 0) == via.get(t, 0) for t in range(NT))

        def C(X):
            p = Al.psi_form(Al.wedge11(X, X))
            return {t: Al.pair_form(p[t], X) for t in p}
        cp = C([add(T[a], al[a]) for a in range(N)])
        cm = C([add(T[a], al[a], -1) for a in range(N)])
        ca = C(al)
        q = Al.quad3(T)
        dC = {t: Al.pair_form(q[t], al) for t in q}
        ok2 = all(cp.get(t, 0) - cm.get(t, 0) == 2 * dC.get(t, 0) + 2 * ca.get(t, 0) for t in range(NT))
        nz = sum(1 for t in range(NT) if dC.get(t, 0) != 0)
        P(f"  {SIG[eta[7]]} dual of alpha -> Psi(T^[T^alpha]) = direct value (random dense T, alpha; all 48 types): {ok1}; "
          f"gradient identity C(T+a) - C(T-a) = 2 dC_T(a) + 2 C(a), C = Psi(T^T^T): {ok2} ({nz} types nonzero)")
        assert ok1 and ok2 and nz > 30
        rng = np.random.default_rng(1)
        R = random_curvature(N, rng, terms=1).astype(np.int64)
        etam = np.diag(eta).astype(np.int64)
        Ric = np.einsum('ac,abcd->bd', etam, R)
        s = int(np.einsum('bd,bd->', etam, Ric))
        G = [[Fraction(int(Ric[a, b])) - Fraction(s, 2) * (eta[a] if a == b else 0) for b in range(N)] for a in range(N)]
        pt = Al.psi_form(z4_form(Al, R, half=True))
        oks = []
        for nm in OPN:
            unit, sv = OPS[nm]
            D = combine(pt, sv)
            if nm in ('W_A', "W'_A"):
                gf = Al.gamma if nm == 'W_A' else Al.gammaG
                ok = True
                for a in range(N):
                    exp = {}
                    for b in range(N):
                        for k, v in gf(b).items():
                            exp[k] = exp.get(k, 0) - eta[a] * G[a][b] * eta[b] * v
                    ok &= {k: complex(unit * v) for k, v in D[a].items()} == {k: complex(v) for k, v in exp.items() if v}
                oks.append(ok)
            else:
                oks.append(not any(D))
        P(f"  {SIG[eta[7]]} Levi-Civita spin curvature F_R (random algebraic R): W_A F_R = -*(G e^a (x) gamma^b): {oks[0]}; "
          f"W'_A F_R = -*(G e^a (x) gamma^b Gamma): {oks[1]}; W_C.AdG F_R = 0: {oks[2]}; W'_C.AdG F_R = 0: {oks[3]}")
        assert all(oks)


def i0_real77():
    """I0b: in (7,7) the four operators in code normalisation map u(64,64)-valued 2-forms to u(64,64)-valued 13-forms
    (explicit real Cl(7,7) of clifford77.py, exact Gaussian-integer arithmetic), so real (x, y, z, w) is the real form."""
    from clifford77 import clifford_pq
    from shiab_checks import Ext, comm, acomm
    rng = np.random.default_rng(20210401)
    pos, neg = clifford_pq(7)
    gam = [g_.astype(complex) for g_ in pos + neg]
    Id = np.eye(128, dtype=complex)

    def mono(idx):
        M = Id.copy()
        for i in idx:
            M = M @ gam[i]
        return M
    Gm, Bf = mono(range(14)), mono(range(7, 14))
    Bfi = np.round(np.linalg.inv(Bf).real).astype(complex)
    Gi = np.round(np.linalg.inv(Gm).real).astype(complex)
    assert np.array_equal(Bfi @ Bf, Id) and np.array_equal(Gi @ Gm, Id)

    def in_u(X):
        return np.array_equal(X.conj().T @ Bf + Bf @ X, 0 * Id)

    def rand_u():
        M = rng.integers(-1, 2, (128, 128)) + 1j * rng.integers(-1, 2, (128, 128))
        return M - Bfi @ M.conj().T @ Bf
    E = Ext(E77)
    Phi11 = E.star({abc: 1j * mono(abc) for abc in itertools.combinations(range(N), 3)})
    Phi11p = E.star({abc: 1j * mono(abc) @ Gm for abc in itertools.combinations(range(N), 3)})
    ok = all(in_u(g_) and in_u(g_ @ Gm) for g_ in gam)
    for trial in range(3):
        X = rand_u()
        assert in_u(X)
        xi = {(0, 1 + trial): X}
        adg = {k: Gm @ v @ Gi for k, v in xi.items()}
        for Pp in (Phi11, Phi11p):
            for br, inp in ((acomm, xi), (comm, adg)):
                out = E.wedge(Pp, inp, br)
                ok &= bool(out) and all(in_u(v) for v in out.values())
    P(f"  (7,7) real form: gamma_a, gamma_a Gamma in u(64,64), and W_A, W'_A, W_C.AdG, W'_C.AdG (code normalisation) map "
      f"random u(64,64)-valued 2-forms into u(64,64) (3 samples, exact): {ok}")
    assert ok


# ====================================================================== I1
def constituent_vectors(Al):
    eta = np.diag(Al.eta).astype(np.int64)
    n = N
    R = random_curvature(n, np.random.default_rng(1), terms=1).astype(np.int64)
    Ric = np.einsum('ac,abcd->bd', eta, R)
    s = int(np.einsum('bd,bd->', eta, Ric))
    W = 312 * R - kn(26 * Ric - s * eta, eta).astype(np.int64)            # 312 * Weyl(R)
    h = np.zeros((n, n), dtype=np.int64)
    h[0, 1] = h[1, 0] = 1
    beta = np.zeros((n, n), dtype=np.int64)
    beta[0, 1], beta[1, 0] = 1, -1
    L4 = np.zeros((n,) * 4, dtype=np.int64)
    for p in itertools.permutations(range(4)):
        L4[p] = perm_sign(p)
    Hk = np.zeros((n,) * 4, dtype=np.int64)
    for (a, b, c, d), v in (((0, 1, 2, 3), 1), ((2, 3, 0, 1), -1)):
        Hk[a, b, c, d] += v
        Hk[b, a, c, d] -= v
        Hk[a, b, d, c] -= v
        Hk[b, a, d, c] += v
    two = [('Weyl', W), ('S^2_0', kn(h, eta).astype(np.int64)), ('trivial', kn(eta, eta).astype(np.int64)),
           ('Lambda^4', L4), ('Lambda^2', kn(beta, eta).astype(np.int64)), ('[2,1,1]_0', Hk)]
    L3 = np.zeros((n,) * 3, dtype=np.int64)
    for p in itertools.permutations(range(3)):
        L3[p] = perm_sign(p)
    v1 = np.zeros((n,) * 3, dtype=np.int64)
    for a in range(n):
        for b in range(n):
            for c in range(n):
                v1[a, b, c] = (a == 0) * eta[b, c] - (b == 0) * eta[a, c]
    h21 = np.zeros((n,) * 3, dtype=np.int64)
    h21[0, 2, 2], h21[2, 0, 2], h21[0, 3, 3], h21[3, 0, 3] = 1, -1, -1, 1
    one = [('Lambda^3', L3), ('Lambda^1', v1), ('[2,1]_0', h21)]
    G = Al.Gam()
    out = [(0, 'Lambda^2', {(0, 1): {0: 1}})]
    out += [(1, nm, z3_form(Al, Z)) for nm, Z in one]
    out += [(2, nm, z4_form(Al, Z)) for nm, Z in two]
    out += [(12, nm, z4_form(Al, Z, G)) for nm, Z in two]
    out += [(13, nm, z3_form(Al, Z, G)) for nm, Z in one]
    out.append((14, 'Lambda^2', {(0, 1): {FULL: 1}}))
    return out, one


def i1_selection():
    P("\n=== I1. Which components of a g-valued 2-form each operator sees (exact) ===")
    P("  (a) grade selection rules (from the H1 coordinates; a type (j,l) pairs input grade l with output grade j and vice versa):")
    for nm in OPN:
        _, sv = OPS[nm]
        rule = {k: sorted({TYPES[t][1] if TYPES[t][2] == k else TYPES[t][2] for t in sv if k in TYPES[t][1:]})
                for k in range(15)}
        P(f"    {nm:9s}: " + "; ".join(f"{k}->{v if v else '-'}" for k, v in rule.items()))
    P("  (b) Spin(14)-constituents of the input (grades 0,1,2 and their Gamma-partners 14,13,12); listed: output grades of"
      " the image of one test vector (each constituent has multiplicity 1, so 'nonzero' = injective on it, Schur):")
    for eta in (EUC, E77):
        Al = Alg(eta)
        vecs, _ = constituent_vectors(Al)
        P(f"   {SIG[eta[7]]}")
        for g, nm, zt in vecs:
            assert all(pc(K) == g for v in zt.values() for K in v)
            pt = Al.psi_form(zt)
            row = []
            for op in OPN:
                D = combine(pt, OPS[op][1])
                row.append(f"{op}: {grades_of(D) if any(D) else 0}")
            P(f"    input grade {g:2d} {nm:10s}  " + "   ".join(row))


# ====================================================================== I2, I3
def i2_i3_parity_twist():
    P("\n=== I2. Parity split of Upsilon (exact, all 48 types) ===")
    P("  Lemma H2: every type pairs an even with an odd grade, so Shiab maps g_ev-valued forms to g_odd-valued ones and back:"
      f" {all((j + l) % 2 for _, j, l in TYPES)}")
    b_odd0 = True
    for eta in (EUC, E77):
        Al = Alg(eta)
        Te = [rand_elt(EVEN, 7) for _ in range(N)]
        To = [rand_elt(ODD, 7) for _ in range(N)]
        T = [add(Te[a], To[a]) for a in range(N)]
        qe, qo, qt = Al.quad3(Te), Al.quad3(To), Al.quad3(T)
        z = [dict() for _ in range(N)]
        ok = all(not any(parity_part(q.get(t, z), 0)) for q in (qe, qo) for t in range(NT))
        cross_ev = sum(1 for t in range(NT) if any(parity_part(qt.get(t, z), 0)))
        b_odd0 = b_odd0 and all(not any(parity_part([add(add(qt.get(t, z)[a], qe.get(t, z)[a], -1), qo.get(t, z)[a], -1)
                                                     for a in range(N)], 1)) for t in range(NT))
        P(f"  {SIG[eta[7]]}: Q_t(T_e) and Q_t(T_o) have zero g_ev-part for all 48 types (random T_e, T_o): {ok}; "
          f"types with nonzero g_ev-part of Q_t(T_e + T_o) (= 2 B_t(T_e,T_o)_ev): {cross_ev}")
        assert ok and cross_ev > 0
        if eta is EUC:
            I3data = (Al, Te, To, T, qe, qo)
        else:
            I3data77 = (Al, Te, To, T, qe, qo)
    assert b_odd0
    P("  => Upsilon_ev = Shiab(d_0 T_o) + kappa *T_e + 2 B(T_e, T_o)_ev  (linear in T_e);"
      " Upsilon_odd = Shiab F_0 + Shiab(d_0 T_e) + kappa *T_o + Q(T_o) + Q(T_e), since 2 B_t(T_e, T_o) = Q_t(T_e + T_o)"
      f" - Q_t(T_e) - Q_t(T_o) has zero g_odd-part for all 48 types, both signatures: {b_odd0}")

    P("\n=== I3. The twisted pair: cubic term on odd/even fields (exact) ===")
    for Al, Te, To, T, qe, qo in (I3data, I3data77):
        for nm in ('W_C.AdG', "W'_C.AdG"):
            _, sv = OPS[nm]

            def cub(X, Y, Z):
                p = Al.psi_form(Al.wedge11(Y, Z))
                return sum(c * Al.pair_form(p.get(t, [{}] * N), X) for t, c in sv.items())
            lhs, rhs = cub(T, T, T), 2 * cub(To, Te, Te)
            Do, De = combine(qo, sv), combine(qe, sv)
            Wee = combine(Al.psi_form(Al.wedge11(Te, Te)), sv)
            Woo = combine(Al.psi_form(Al.wedge11(To, To)), sv)
            ok_e = all(De[a] == {k: 2 * v for k, v in Wee[a].items()} for a in range(N))
            P(f"  {SIG[Al.eta[7]]} {nm:9s}: C(T) = 2 Psi(T_o ^ T_e ^ T_e): {lhs == rhs} (value {lhs}); Q(T_o) = 0: "
              f"{not any(Do)}; Q(T_e) = (2/3) Shiab(T_e ^ T_e): {ok_e}; but Shiab(T_o ^ T_o) != 0: {any(Woo)} "
              f"(cancelled by -h/3)")
            assert lhs == rhs and lhs != 0 and not any(Do) and ok_e and any(Woo)
        rows = {}
        for smp in range(3):
            Ts = [rand_elt(ODD, 3) for _ in range(N)]
            for t, D in Al.quad3(Ts).items():
                for a in range(N):
                    for k, v in D[a].items():
                        rows.setdefault((smp, a, k), {})[t] = v
        M = sp.Matrix([[r.get(t, 0) for t in range(NT)] for r in rows.values()])
        ns = M.nullspace()
        desc = []
        for v in ns:
            sup = {NAME[i]: v[i] for i in range(NT) if v[i] != 0}
            hit = [nm for nm in OPS if {NAME[t]: c for t, c in OPS[nm][1].items()} == sup or
                   {NAME[t]: -c for t, c in OPS[nm][1].items()} == sup]
            desc.append(hit[0] if hit else '+'.join(sup))
        P(f"  {SIG[Al.eta[7]]}: {{Shiab in S_sa : Q(T) = 0 for the 3 random odd T}} has exact dimension {len(ns)}: "
          f"{desc}  ({M.shape[0]} x {NT} exact system)")
        assert len(ns) == 4


# ====================================================================== I4
def i4_soldering():
    P("\n=== I4. Soldering sector T = theta_ab e^a(x)gamma^b + theta'_ab e^a(x)gamma^b Gamma, T_e = 0 (exact) ===")
    a1, a2, a3, a4 = sp.symbols('a1:5')
    for eta in (EUC, E77):
        Al = Alg(eta)
        th, thp, be = rsym(), rsym(), rasym()
        cases = [('sym theta', T_from(Al, th)), ("sym theta, sym theta'", T_from(Al, th, thp)),
                 ('theta + antisym beta', T_from(Al, [[th[a][b] + be[a][b] for b in range(N)] for a in range(N)]))]
        for lab, T in cases:
            q = Al.quad3(T)
            P(f"  {SIG[eta[7]]} output grades of Q(T), {lab:22s}: " +
              "  ".join(f"{op}: {grades_of(combine(q, OPS[op][1]))}" for op in OPN))
        E = np.diag(eta).astype(object)
        eqs = {'W_A': [], "W'_A": []}
        for smp in range(2):
            thn = rsym(r=3, dens=0.5)
            Th = np.array(thn, dtype=object)
            th2 = Th.dot(E).dot(Th)
            tr1 = sum(E[a, a] * Th[a, a] for a in range(N))
            tr2 = sum(th2[a, a] * E[a, a] for a in range(N))
            q = Al.quad3(T_from(Al, thn))
            for op, gf in (('W_A', Al.gamma), ("W'_A", Al.gammaG)):
                unit, sv = OPS[op]
                D = combine(q, sv)
                basis = {next(iter(gf(b))) for b in range(N)}
                assert all(k in basis for Da in D for k in Da)
                for a in range(N):
                    for b in range(N):
                        (k, sg), = gf(b).items()
                        d = Fraction(D[a].get(k, 0), sg)
                        Mab = sp.Rational(d.numerator, d.denominator) * eta[a] * eta[b] * unit_sp(unit) / 3
                        eqs[op].append(sp.Eq(Mab, a1 * th2[a, b] + a2 * tr1 * Th[a, b] + (a3 * tr2 + a4 * tr1 ** 2) * E[a, b]))
        for op, chan in (('W_A', 'gamma'), ("W'_A", 'gamma Gamma')):
            sol = sp.solve(eqs[op], [a1, a2, a3, a4], dict=True)
            assert len(sol) == 1
            P(f"  {SIG[eta[7]]} {op}, symmetric theta: Q(T) = *(M_ab e^a (x) {chan}^b), M = a1 theta^2 + a2 tr(theta) theta + "
              f"(a3 tr(theta^2) + a4 tr(theta)^2) g: {sol[0]} ({len(eqs[op])} exact equations, 2 samples)")
        P("    [-G(R_theta) with R_theta = 4(theta_ac theta_bd - theta_ad theta_bc) has (a1..a4) = (4, -4, -2, 2): "
          "so Q_{W_A}(theta) = W_A(T^T) (h vanishes) and Q_{W'_A}(theta) = W'_A(T^T)/3]")
    P("  constraint (g_ev-part) equations Shiab(d_0 T) = 0 for T = theta e(x)gamma + theta' e(x)gamma Gamma: per Spin(14)-"
      "constituent rho of d^nabla theta, d^nabla theta' in Lambda^2(x)Lambda^1, input a v + b v Gamma (v in rho):")
    for eta in (EUC, E77):
        Al = Alg(eta)
        _, one = constituent_vectors(Al)
        for nm, Z in one:
            outs = {}
            for lab, right in (('a', None), ('b', Al.Gam())):
                pt = Al.psi_form(z3_form(Al, Z, right))
                for op in OPN:
                    unit, sv = OPS[op]
                    D = combine(pt, sv)
                    for g in grades_of(D):
                        outs.setdefault(g, []).append((op, lab, unit, [{k: v for k, v in Da.items() if pc(k) == g} for Da in D]))
            rows, gl = [], []
            for g, lst in sorted(outs.items()):
                ref = lst[0][3]
                a0 = next(a for a in range(N) if ref[a])
                k0 = next(iter(ref[a0]))
                row = [0, 0]
                for op, lab, unit, Dg in lst:
                    r = Fraction(Dg[a0].get(k0, 0), ref[a0][k0])
                    assert all({k: Fraction(v) for k, v in Dg[a].items()} == {k: r * v for k, v in ref[a].items()} for a in range(N))
                    row[0 if lab == 'a' else 1] += unit_sp(unit) * sp.Rational(r.numerator, r.denominator) * FAMC[op]
                rows.append(row)
                gl.append(g)
            M = sp.Matrix(rows)
            minors = sorted({sp.factor(M[[i, j], :].det()) for i, j in itertools.combinations(range(M.shape[0]), 2)} - {0},
                            key=str)
            P(f"   {SIG[eta[7]]} {nm:9s}: output grades {gl}; rows [a, b] = {[list(r_) for r_ in rows]}; nonzero 2x2 minors "
              f"{minors}; rank at Shiab# (x=z=1): {M.subs({xs: 1, zs: 1, ys: 0, ws: 0}).rank()}, at W_A: "
              f"{M.subs({xs: 1, zs: 0, ys: 0, ws: 0}).rank()}, at W'_A: {M.subs({xs: 0, zs: 0, ys: 1, ws: 0}).rank()}")


# ====================================================================== I5
def i5_invariant():
    P("\n=== I5. Spin(14)-invariant sector T = t Phi^1 + u Phi^1 Gamma (t, u functions; exact) ===")
    t, u, kap = sp.symbols('t u kappa')
    res = {}
    for eta in (EUC, E77):
        Al = Alg(eta)
        Phi1 = [Al.gamma(a) for a in range(N)]
        Phi1G = [Al.gammaG(a) for a in range(N)]
        both = [add(Phi1[a], Phi1G[a]) for a in range(N)]
        vals = {}
        for lab, T in (('tt', Phi1), ('uu', Phi1G), ('sum', both)):
            q = Al.quad3(T)
            for op in OPN:
                unit, sv = OPS[op]
                D = combine(q, sv)
                c1 = sp.Rational(D[0].get(1, 0) * eta[0])
                (k0, s0), = Al.gammaG(0).items()
                c13 = sp.Rational(D[0].get(k0, 0) * s0 * eta[0])
                for a in range(N):                  # D_a = eta_aa (c1 gamma_a + c13 gamma_a Gamma) exactly
                    exp = add(scale(Al.gamma(a), c1 * eta[a]), scale(Al.gammaG(a), c13 * eta[a]))
                    assert {k: sp.Rational(v) for k, v in D[a].items()} == {k: sp.Rational(v) for k, v in exp.items()}
                vals[(lab, op)] = (unit_sp(unit) * c1 / 3, unit_sp(unit) * c13 / 3)
        C1 = kap * t
        C13 = kap * u
        for op in OPN:
            tt, uu, sm = vals[('tt', op)], vals[('uu', op)], vals[('sum', op)]
            cr = [sm[i] - tt[i] - uu[i] for i in range(2)]
            C1 += FAMC[op] * sp.expand(tt[0] * t ** 2 + uu[0] * u ** 2 + cr[0] * t * u)
            C13 += FAMC[op] * sp.expand(tt[1] * t ** 2 + uu[1] * u ** 2 + cr[1] * t * u)
        C1, C13 = sp.expand(C1), sp.expand(C13)
        res[eta[7]] = (C1, C13)
        P(f"  {SIG[eta[7]]}: algebraic part (every T-term is Spin(14)-invariant, checked componentwise):  Upsilon = "
          f"-*([x G - c1 g]_ab e^a(x)gamma^b + [y G - c13 g]_ab e^a(x)gamma^b Gamma) + derivative terms,")
        P(f"           c1  = {C1}")
        P(f"           c13 = {C13}")
        P(f"           (z, w absent: {not ({zs, ws} & (C1.free_symbols | C13.free_symbols))})")
        assert not ({zs, ws} & (C1.free_symbols | C13.free_symbols))
        # derivative terms: d_0 T = dt ^ Phi^1 + du ^ Phi^1 Gamma (Phi^1, Gamma parallel); directions c = 0, 5, 9
        der = {}
        for lab, P1 in (('dt', Phi1), ('du', Phi1G)):
            for c in (0, 5, 9):
                zt = {}
                for a in range(N):
                    if a != c:
                        zt[(min(a, c), max(a, c))] = scale(P1[a], 1 if c < a else -1)
                pt = Al.psi_form(zt)
                for op in OPN:
                    unit, sv = OPS[op]
                    D = combine(pt, sv)
                    der.setdefault((lab, op), []).append(D)
        rows = {}
        allg = sorted({pc(k) for D_ in der.values() for D in D_ for Da in D for k in Da})
        for ic, c in enumerate((0, 5, 9)):
            for g in allg:
                ref = None
                row = [0, 0]
                for i, lab in enumerate(('dt', 'du')):
                    for op in OPN:
                        unit, _ = OPS[op]
                        Dg = [{k: v for k, v in Da.items() if pc(k) == g} for Da in der[(lab, op)][ic]]
                        if not any(Dg):
                            continue
                        if ref is None:
                            ref = Dg
                            a0 = next(a for a in range(N) if ref[a])
                            k0 = next(iter(ref[a0]))
                        rr = Fraction(Dg[a0].get(k0, 0), ref[a0][k0])
                        assert all({k: Fraction(v) for k, v in Dg[a].items()} == {k: rr * v for k, v in ref[a].items()}
                                   for a in range(N))       # same Spin(14)-pattern (Schur), exact
                        row[i] += unit_sp(unit) * sp.Rational(rr.numerator, rr.denominator) * FAMC[op]
                rows.setdefault(g, []).append(row)
        same = all(r_ == rows[g][0] for g in rows for r_ in rows[g])
        Mrows = [rows[g][0] for g in allg]
        M = sp.Matrix(Mrows)
        P(f"  {SIG[eta[7]]}: derivative terms Shiab(dt ^ Phi^1 + du ^ Phi^1 Gamma) occupy output grades {allg}; rows [dt, du] "
          f"(relative to a fixed pattern) = {Mrows}; identical for directions e^0, e^5, e^9: {same};  det = {sp.factor(M.det())}")
        assert same and allg == [2, 12]
        # reading (P) keeps only form indices b tangent to s(X): still every component of dt (du) is determined
        okP = True
        for lab, P1, op in (('dt', Phi1, 'W_A'), ('du', Phi1G, 'W_C.AdG')):
            pats = {}
            for c in range(N):
                zt = {(min(a, c), max(a, c)): scale(P1[a], 1 if c < a else -1) for a in range(N) if a != c}
                pats[c] = combine(Al.psi_form(zt), OPS[op][1])
            for b in (0, 1):
                mons = set()
                for c in range(N):
                    Db = pats[c][b]
                    okP &= (not Db) if c == b else (bool(Db) and not (set(Db) & mons))
                    mons |= set(Db)
        P(f"  {SIG[eta[7]]}: W_A(dt^Phi^1) and W_C.AdG(du^Phi^1 Gamma): the form-index-b component involves every d_c (c != b)"
          f" through pairwise distinct monomials (b = 0, 1), so two tangential indices b already force dt = 0 (du = 0): {okP}")
        assert okP
    return res


# ====================================================================== I6
def i6_reduction(res):
    P("\n=== I6. Reduction along s(X), n = 4, d = 14 (closed forms of F3-F6; exact sympy) ===")
    lam, cH, cV, kap = sp.symbols('lambda c_H c_V kappa')
    n = 4
    s1 = -sp.Rational(23, 2) - (1 + lam) / (1 - 4 * lam)                         # F4(b), n = 4
    LP_Z1 = sp.simplify(cH / (4 * cV * (n * lam - 1)) - cH * s1 / (2 * cV))         # from F4(b),(c): G^Y_HH = G(g) + LP g
    LP_Z2 = sp.simplify(cH * sp.Rational(n * (n - 1) * (n + 2), 16) / cV)
    ok = sp.simplify(LP_Z1 - sp.Rational(3, 2) * cH / cV * (4 - 15 * lam) / (1 - 4 * lam)) == 0
    P(f"  Lambda_P(Z1) = {sp.factor(LP_Z1)} (= F5 table: {ok});  Lambda_P(Z2) = {LP_Z2};  Lambda_Y (F6: Z1, lambda = 1/5, "
      f"Ric(g) = 0, G^Y = -Lambda_Y g_aleph on s(X)) = {-sp.Rational(12 * 5, 8) / cV}")
    assert ok
    LY = -sp.Rational(15, 2) / cV
    s_, r_ = sp.symbols('s r', real=True)
    xr, yr, zr, wr = sp.symbols("x' y' z' w'", real=True)
    for sig, C in res.items():
        C1, C13 = C
        if sig == 1:      # u(128): Shiab = -i(x'W_A + ...), T = i(s Phi^1 + r Phi^1 Gamma), kappa real
            sub = {xs: -sp.I * xr, ys: -sp.I * yr, zs: -sp.I * zr, ws: -sp.I * wr, sp.Symbol('t'): sp.I * s_,
                   sp.Symbol('u'): sp.I * r_}
            # x G - c1 g = -i (x' G - i c1 g): Upsilon_gamma = 0 <=> x' G + e1 g = 0 with e1 = -i c1
            e1 = sp.expand(-sp.I * C1.subs(sub))
            e13 = sp.expand(-sp.I * C13.subs(sub))
            P(f"  (14,0), u(128)-real (x = -i x' etc., t = i s, u = i r, s, r, x', y', kappa real):  Upsilon = 0 <=>"
              f"  x' G(g_aleph) + ({e1}) g_aleph = 0  and  y' G(g_aleph) + ({e13}) g_aleph = 0")
            ok_real = all(sp.im(cf) == 0 for cf in sp.Poly(e1, s_, r_, xr, yr, kap).coeffs() +
                          sp.Poly(e13, s_, r_, xr, yr, kap).coeffs())
            P(f"     (real polynomials: {ok_real})")
            assert ok_real
            sh1, sh13 = e1.subs({xr: 1, yr: 0}), e13.subs({xr: 1, yr: 0})
            P(f"   Shiab# (x' = z' = 1, y' = w' = 0) [algebraic terms identical for W_A: z' absent]:  G(g_aleph) = -({sh1}) g_aleph,"
              f"  0 = {sp.factor(sh13)}")
            br1 = sh1.subs(r_, 0)
            sq1 = sp.expand(312 * (s_ + kap / 624) ** 2 - kap ** 2 / 1248)
            br2 = sp.expand(sh1.subs(s_, -kap / 208))
            P(f"   (P): Lambda_4 = Lambda_P + c_H*({sh1}) on {{{sp.factor(sh13)} = 0}}:  branch r = 0: {br1} = "
              f"312 (s + kappa/624)^2 - kappa^2/1248 [{sp.expand(br1 - sq1) == 0}];  branch s = -kappa/208: {br2};  so for c_H > 0"
              f" Lambda_4 takes exactly the values in [Lambda_P - c_H kappa^2/1248, oo): not fixed by (lambda, c_H, c_V)")
            assert sp.expand(br1 - sq1) == 0 and sp.expand(br2 - (kap ** 2 / 416 + 104 * r_ ** 2)) == 0
            eqE = sp.expand(-LY + sh1)                         # G = -LY g on s(X) (Z1, lambda = 1/5, Ric(g) = 0)
            disc = sp.discriminant(eqE.subs(r_, 0), s_)
            P(f"   (E) on s(X) [possible only for Z1, lambda = 1/5, Ric(g) = 0 (F6); base (4,0), fiber (10,0) for c_V > 0]:  "
              f"{eqE} = 0 and {sp.factor(sh13)} = 0.  Branch r = 0: discriminant {sp.factor(disc)} >= 0 <=> c_V kappa^2 >= 9360"
              f" (c_V > 0).  Branch s = -kappa/208: {sp.expand(eqE.subs(s_, -kap / 208))} = 0, no real solution for c_V > 0.")
            assert sp.simplify(disc - (kap ** 2 - 9360 / cV)) == 0
            assert sp.expand(eqE.subs(s_, -kap / 208) - (15 / (2 * cV) + kap ** 2 / 416 + 104 * r_ ** 2)) == 0
        else:             # u(64,64): real x, y, z, w, t, u, kappa (code normalisation)
            tt, uu = sp.Symbol('t'), sp.Symbol('u')
            sh1, sh13 = C1.subs({xs: 1, ys: 0}), C13.subs({xs: 1, ys: 0})
            P(f"  (7,7), u(64,64)-real (code normalisation, real x,y,z,w,t,u,kappa):  Upsilon = 0 <=> x G = c1 g and y G = c13 g")
            P(f"   Shiab# (x = z = 1, y = w = 0) [algebraic terms identical for W_A]:  G(g_aleph) = ({sh1}) g_aleph,  0 = {sp.factor(sh13)}")
            br1 = sh1.subs(uu, 0)
            sq1 = sp.expand(312 * (tt + kap / 624) ** 2 - kap ** 2 / 1248)
            br2 = sp.expand(sh1.subs(tt, -kap / 208))
            P(f"   (P): Lambda_4 = Lambda_P - c_H*({sh1}) on {{{sp.factor(sh13)} = 0}}:  branch u = 0: c1 = {br1} = "
              f"312 (t + kappa/624)^2 - kappa^2/1248 [{sp.expand(br1 - sq1) == 0}], so Lambda_4 - Lambda_P takes every value"
              f" <= c_H kappa^2/1248;  branch t = -kappa/208: c1 = {br2}, so Lambda_4 - Lambda_P takes every value >="
              f" -c_H kappa^2/416;  union (c_H > 0): every real Lambda_4")
            assert sp.expand(br1 - sq1) == 0 and sp.expand(br2 - (kap ** 2 / 416 - 104 * uu ** 2)) == 0
            P("   (E) on s(X): impossible in (7,7) for every T in this sector and every (x,y,z,w) with (x,y) != 0: it needs"
              " g_aleph Einstein on s(X), which by F6 forces Z1 with lambda = 1/5, signature (8,6)/(10,4), not (7,7)")
    P("  Symmetric displacement theta (I4) for Shiab = x W_A + ... (x != 0, y = 0): Upsilon = 0 <=> d^nabla theta = 0 and"
      " G(g_aleph) + G(R_theta) = (kappa/x) theta; (P)-pull-back: G(g) + Lambda_P g + s^*G(R_theta) - (kappa/x) s^*theta = 0"
      " (write-up I7(c)); z, w absent.")


def main():
    i0_validate()
    i0_real77()
    i1_selection()
    i2_i3_parity_twist()
    i4_soldering()
    res = i5_invariant()
    i6_reduction(res)
    P(f"\nruntime {time.time() - T0:.0f} s")


if __name__ == '__main__':
    main()
