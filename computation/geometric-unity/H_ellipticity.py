"""GU round 4 (H): is the bosonic deformation complex elliptic in Euclidean signature for SOME self-adjoint
Shiab operator?  (GU-057, GU-072, GU-G16.)

Formalization: formalization/geometric-unity/H-ellipticity.md.  Builds on D (Theorem D2: S_sa = {Shiab = Psi ^ .},
Psi in Hom_K(Lambda^3, S^2 g)) and G (Prop. G10, Theorem G11: the symbol sequence at xi is exact iff
rank sigma(delta_2)(xi) = 14 dim g).  Reuses G_symbols.py (monomial Clifford algebra, exact sparse ranks, the
W_A / W'_A operator code), D_rep_counts.py (D_7 characters) and shiab_checks.py (curvature helpers); none is modified.

Setting.  d = 14, eta = (+^14) (and (+^7,-^7) for comparison); g = Cl(eta) (x) C = gl(128, C) with real forms
u(128) [Euclidean] / u(64,64) [draft]; Clifford monomials E_K, K a bitmask of {0..13}.
Basis of S_sa (H0).  For each (j, l), j < l, the D_7-invariants Lambda^3 -> Lambda^j (x) Lambda^l are
   O_(j,l)(w; X, Y) = Tr(q(w) X_j Y_l) + Tr(q(w) Y_j X_l)            (nonzero iff A xor J xor L = 0 is possible)
   E_(j,l)(w; X, Y) = Tr(Gamma q(w) X_j Y_l) + Tr(Gamma q(w) Y_j X_l)  (nonzero iff A xor J xor L = all is possible)
with X_j the grade-j part, q(e^A) = gamma^A = eta^A E_A the (Pin-equivariant) quantization, Gamma = E_full, Tr/128.
Symbol (H2).  At xi = e^0 the symbol of Shiab o d restricted to W (x) g, W = xi^perp, is the antisymmetric form
   beta(e^a X, e^c Y) = Psi(e^0 ^ e^a ^ e^c; X, Y);
exactness of the complex at xi  <=>  kappa != 0 and beta nondegenerate (Lemma H1 of the write-up).
Schur reduction (H3).  Under Spin(13) = Stab(xi), W (x) g = sum_rho m_rho rho with rho in {Lambda^0..Lambda^6 (m = 4, 8),
H_1..H_6 (hooks [2,1^(j-1)], m = 4)}; copies are labelled by Clifford grades; beta = sum_rho C_rho (x) <,>_rho,
and beta is nondegenerate iff det C_rho != 0 for all rho (C_rho = odd-copies x even-copies block).
Tests
  H0  basis: 48 types = D_7 character multiplicities (block by block); symmetry, parity, exact infinitesimal
      Spin(14)-invariance on random monomials; real structure for u(128).
  H1  W_A, W'_A (G_symbols code), W_C o Ad_Gamma, W'_C o Ad_Gamma: exact coordinates in the basis (overdetermined
      exact solve); W_C alone is not in S_sa; entrywise identity of the assembled symbol with G_symbols' matrix.
  H2  Schur blocks C_rho in (14,0): consistency (antisymmetry, parity, three test vectors, orthogonality across rho),
      sum m_rho dim rho = 212992; exact determinants and their factorization.
  H3  discriminant and exact full sparse ranks of beta (212992 x 212992) at W_A, W'_A, W_A + W'_A, a random point,
      points on single components; Schur-predicted ranks must agree.
  H4  Einstein output of every basis element on Levi-Civita curvature (exact rationals).
  H5  draft's family span(W_A, W'_A) and the twisted family span(W_A, W'_A, W_C Ad_Gamma, W'_C Ad_Gamma):
      determinants, sum-of-squares identities, an explicit u(128)-real elliptic operator, rank sigma(delta_2) with kappa.
  H6  so(64,64) / so(128): 12-dim S_sa, Schur blocks, maximal rank.
  H7  (7,7) at a non-null covector, for comparison with G11: pencil determinants and a generic point.
Exact arithmetic throughout: Python integers / sympy rationals and Gaussian rationals; sparse ranks by exact integer
elimination (G_symbols.sparse_rank).  Fixed seeds.  Runtime is dominated by sympy's factorisation of the five quartics
over Q(i) (H3); H_SKIP_QI=1 skips that step only (the committed .out is a full run).
"""
import itertools
import os
import sys
import time
from math import comb

import numpy as np
import scipy.sparse as sps
import sympy as sp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import G_symbols as GS  # noqa: E402
from G_symbols import Cliff, mask, N, FULL, DIMG, vsign, sparse_rank, apply_op, _KS  # noqa: E402

pc = lambda x: bin(x).count('1')  # noqa: E731
EUC = [1] * 14
E77 = [1] * 7 + [-1] * 7
WB = list(range(1, 14))          # W = xi^perp for xi = e^0
WFULL = FULL ^ 1
GSO = {1, 2, 5, 6, 9, 10, 13, 14}  # Clifford grades of so(64,64) / so(128)


# ====================================================================== H0. basis of S_sa = Hom_K(Lambda^3, S^2 g)
def make_types():
    out = []
    for kind in 'OE':
        for j in range(15):
            for l in range(j + 1, 15):
                ll = l if kind == 'O' else 14 - l          # E_(j,l) <-> O_(j,14-l) via Y -> Gamma Y
                if (3 + j + ll) % 2:
                    continue
                s_aj, s_al, s_jl = (3 + j - ll) // 2, (3 + ll - j) // 2, (j + ll - 3) // 2
                if min(s_aj, s_al, s_jl) < 0 or s_aj + s_al + s_jl > 14:
                    continue
                out.append((kind, j, l))
    return out


TYPES = make_types()
IDX = {t: i for i, t in enumerate(TYPES)}
NT = len(TYPES)
NAME = [f"{k}{j},{l}" for k, j, l in TYPES]


def eps_u(k):
    """gamma_K (|K| = k) is Hermitian in (14,0) iff k = 0,1 mod 4; the u(128) basis is X_K = i^eps(k) gamma_K."""
    return 1 if k % 4 in (0, 1) else 0


class Psi:
    """Basis evaluation Psi_t(e^A; E_J, E_L) (e^A a sorted monomial 3-form); returns {t: integer}."""

    def __init__(self, eta):
        self.eta = list(eta)
        self.cl = Cliff(eta)

    def etaM(self, M):
        s = 1
        for i in range(N):
            if M >> i & 1:
                s *= self.eta[i]
        return s

    def chain(self, masks):
        c, cur = 1, 0
        for m in masks:
            c *= self.cl.sign(cur, m)
            cur ^= m
        return c, cur

    def tr(self, kind, Z, X, Y):
        """Tr((Gamma) E_Z E_X E_Y)/128 for monomials (Gamma iff kind == 'E')"""
        c, x = self.chain(([FULL] if kind == 'E' else []) + [Z, X, Y])
        return c if x == 0 else 0

    def eval_cl(self, Z, X, Y):
        """Psi_t(Z; X, Y) with Z, X, Y Clifford monomials (Z the quantized 3-form itself): {t: value}"""
        gx, gy = pc(X), pc(Y)
        if gx == gy:
            return {}
        lo, hi = (X, Y) if gx < gy else (Y, X)
        key = (min(gx, gy), max(gx, gy))
        out = {}
        for kind in 'OE':
            t = (kind,) + key
            if t in IDX:
                v = self.tr(kind, Z, lo, hi)
                if v:
                    out[IDX[t]] = v
        return out

    def eval(self, A, J, L):
        eA = self.etaM(A)
        return {t: eA * v for t, v in self.eval_cl(A, J, L).items()}


def lin_eval(P, Zd, Xd, Yd):
    """Psi on linear combinations (dicts mask -> coefficient): vector over types"""
    out = [0] * NT
    for Z, cz in Zd.items():
        for X, cx in Xd.items():
            for Y, cy in Yd.items():
                for t, v in P.eval_cl(Z, X, Y).items():
                    out[t] += cz * cx * cy * v
    return out


def witness(t, perm):
    """monomials (A, J, L) with Psi_t(A; J, L) != 0: A = P u Q, J = P u R, L = Q u R (O-type, |P|,|Q|,|R| fixed by the
    grades), and L -> complement for E-type; indices relabelled by the permutation perm"""
    kind, j, l = t
    ll = l if kind == 'O' else 14 - l
    s_aj, s_al, s_jl = (3 + j - ll) // 2, (3 + ll - j) // 2, (j + ll - 3) // 2
    Pp, Qq, Rr = perm[:s_aj], perm[s_aj:s_aj + s_al], perm[s_aj + s_al:s_aj + s_al + s_jl]
    A, J, L = mask(Pp + Qq), mask(Pp + Rr), mask(Qq + Rr)
    return A, J, (L if kind == 'O' else FULL ^ L)


def comm_mono(cl, g, K):
    """[E_g, E_K] as dict"""
    s1, s2 = cl.sign(g, K), cl.sign(K, g)
    return {g ^ K: s1 - s2} if s1 != s2 else {}


def h0_basis(rng):
    print("=== H0. Basis of S_sa = Hom_Spin(14)(Lambda^3, S^2 g), g = gl(128) (real forms u(128), u(64,64)) ===")
    from D_rep_counts import Dm, hom
    G = Dm(7)
    L = {k: G.ext(k) for k in range(15)}
    w3 = (1, 1, 1, 0, 0, 0, 0)
    mult = {}
    for j in range(15):
        for l in range(j + 1, 15):
            mult[(j, l)] = hom({w3: 1}, G.decompose(L[j] * L[l]))
    diag = sum(hom({w3: 1}, G.decompose(L[j].sym2())) for j in range(15))
    ours = {}
    for k, j, l in TYPES:
        ours[(j, l)] = ours.get((j, l), 0) + 1
    agree = all(ours.get(b, 0) == m for b, m in mult.items()) and diag == 0
    print(f"  {NT} types ({sum(t[0] == 'O' for t in TYPES)} O-type, {sum(t[0] == 'E' for t in TYPES)} E-type); "
          f"D_7 characters: dim Hom(Lambda^3, Lambda^j (x) Lambda^l) summed = {sum(mult.values())}, "
          f"diagonal S^2 Lambda^j blocks = {diag}; block-by-block agreement with the types: {agree}")
    assert agree and NT == 48
    print("  types: " + ", ".join(f"p{i}={NAME[i]}" for i in range(NT)))
    print(f"  parity: every type pairs an even with an odd grade (j + l odd): {all((j + l) % 2 for _, j, l in TYPES)}")
    for eta, lab in ((EUC, '(14,0)'), (E77, '(7,7)')):
        P = Psi(eta)
        cl = P.cl
        ok_inv = ok_sym = ok_wit = True
        hits = set()
        for it in range(4000):
            a, b = sorted(int(x) for x in rng.choice(N, 2, replace=False))
            g = mask((a, b))
            t0 = it % NT
            if it < 3000:                       # a random witness of type t0 (relabelled indices)
                A, J, L_ = witness(TYPES[t0], [int(x) for x in rng.permutation(N)])
            else:                               # unconstrained random monomials
                A = mask(tuple(int(x) for x in rng.choice(N, 3, replace=False)))
                J = mask(tuple(int(x) for x in rng.choice(N, TYPES[t0][1], replace=False)))
                L_ = mask(tuple(int(x) for x in rng.choice(N, TYPES[t0][2], replace=False)))
            e1 = P.eval_cl(A, J, L_)
            if it < 3000:
                ok_wit &= e1.get(t0, 0) != 0
            hits |= set(e1)
            ok_sym &= e1 == P.eval_cl(A, L_, J)
            tot = [x + y + z for x, y, z in zip(lin_eval(P, comm_mono(cl, g, A), {J: 1}, {L_: 1}),
                                                lin_eval(P, {A: 1}, comm_mono(cl, g, J), {L_: 1}),
                                                lin_eval(P, {A: 1}, {J: 1}, comm_mono(cl, g, L_)))]
            ok_inv &= all(v == 0 for v in tot)
        print(f"  {lab}: 4000 samples (3000 on witnesses of each type, 1000 random): symmetric in (X, Y): {ok_sym}; exact "
              f"infinitesimal Spin(14)-invariance Psi([g,.],.,.) + Psi(.,[g,.],.) + Psi(.,.,[g,.]) = 0: {ok_inv}; "
              f"each type nonzero on its witnesses (distinct supports => independent): {ok_wit and len(hits) == NT}")
        assert ok_sym and ok_inv and ok_wit and len(hits) == NT
    # real structure: on X_J = i^eps gamma_J, Psi_t is real iff p_t in i^(eps_j + eps_l) R
    ph = [sp.I ** (eps_u(j) + eps_u(l)) for _, j, l in TYPES]
    print("  real form u(128) (Euclidean gamma_a Hermitian, X_K = i^eps(|K|) gamma_K anti-Hermitian): Psi = sum p_t Psi_t is "
          "real on u(128) iff p_t = i^(eps_j+eps_l) r_t, r_t real; in the basis X_K the matrix of beta at r equals the "
          "gamma-basis matrix at the integer point p'_t = (-1)^(eps_j+eps_l) r_t.  Phases i^(eps_j+eps_l) = i on "
          f"{sum(1 for z in ph if z == sp.I)} types, -1 on {sum(1 for z in ph if z == -1)}, 1 on {sum(1 for z in ph if z == 1)}")
    return ph


# ====================================================================== H1. operators of the draft family in the basis
def shiab_ops_ext(eta, which):
    """G_symbols.shiab_ops for WA/WpA (identical code path) plus the commutator versions WC/WpC."""
    if which in ('WA', 'WpA'):
        return GS.shiab_ops(eta, which)
    cl = Cliff(eta)
    E = GS.Ext(eta)
    prime = which == 'WpC'
    Phi3 = {abc: GS.CL({mask(abc) ^ (FULL if prime else 0): 1j * (cl.sign(mask(abc), FULL) if prime else 1)})
            for abc in itertools.combinations(range(N), 3)}
    Phi11 = E.star(Phi3)
    B = GS.br(cl, 'C')
    ops = {}
    for b in range(N):
        for a in range(N):
            if a == b:
                continue
            res = E.wedge(Phi11, {tuple(sorted((a, b))): (1 if b < a else -1) * GS.IDOP}, B)
            d = {}
            for I, op in res.items():
                m = next(i for i in range(N) if i not in I)
                op = GS.OP({k: v for k, v in op.d.items() if v != 0})
                if op.d:
                    d[m] = op
            ops[(b, a)] = d
    return cl, E, ops


def fit_operator(eta, which, twist, rng, nsamp=2000):
    """Exact coordinates of the operator in the Psi basis, from Psi(e^c ^ e^b ^ e^a; E_J, E_K) =
    (-1)^c Tr(E_J Shiab_c(e^b ^ e^a (x) E_K))/128 (Theorem D2) on nsamp random (c, b, a, K) and all J that can
    occur; an overdetermined exact linear system.  twist: precompose with Ad_Gamma (E_K -> (-1)^|K| E_K)."""
    cl, _, ops = shiab_ops_ext(eta, which)
    P = Psi(eta)
    rows, rhs = [], []
    for _ in range(nsamp):
        c, b, a = (int(x) for x in rng.choice(N, 3, replace=False))
        g = int(rng.integers(0, 15))
        K = mask(tuple(int(x) for x in rng.choice(N, g, replace=False)))
        A = mask((c, b, a))
        sgn = perm_sign((c, b, a))                    # e^c ^ e^b ^ e^a = sgn e^A
        img = apply_op(cl, ops[(b, a)][c], K) if c in ops[(b, a)] else {}
        tw = (-1) ** pc(K) if twist else 1
        for J in sorted(set(img) | {K ^ A, K ^ A ^ FULL}):
            val = complex((-1) ** c * img.get(J, 0) * cl.sign(J, J) * tw)
            row = [0] * NT
            for t, v in P.eval(A, J, K).items():
                row[t] = sgn * v
            rows.append(row)
            rhs.append(sp.Integer(round(val.real)) + sp.I * sp.Integer(round(val.imag)))
    M = sp.Matrix(rows)
    try:
        sol, params = M.gauss_jordan_solve(sp.Matrix(rhs))
    except ValueError:
        return None, len(rows)
    assert params.shape[0] == 0
    return [sp.nsimplify(sol[i]) for i in range(NT)], len(rows)


def perm_sign(seq):
    s = 1
    for x in range(len(seq)):
        for y in range(x + 1, len(seq)):
            if seq[x] > seq[y]:
                s = -s
    return s


def fmt_vec(v):
    return "{" + ", ".join(f"{NAME[i]}: {v[i]}" for i in range(NT) if v[i] != 0) + "}"


# ====================================================================== the bilinear form beta at xi = e^0
POP = np.bitwise_count(_KS).astype(np.int64)
TYPE_TAB = {k: -np.ones((15, 15), dtype=np.int64) for k in 'OE'}
for _t, _i in IDX.items():
    TYPE_TAB[_t[0]][_t[1], _t[2]] = _i
    TYPE_TAB[_t[0]][_t[2], _t[1]] = _i


def beta_entries(eta, p, grades=None):
    """Sparse integer matrix of beta on W (x) g: entry [(c-1) 2^14 + L, (a-1) 2^14 + K] = Psi(e^0^e^a^e^c; E_K, E_L);
    grades: restrict g to these Clifford grades (so(64,64))."""
    cl = Cliff(eta)
    p = np.asarray([int(x) for x in p], dtype=np.int64)
    keep = np.ones(DIMG, dtype=bool) if grades is None else np.isin(POP, sorted(grades))
    rows, cols, vals = [], [], []
    for a in WB:
        for c in WB:
            if a == c:
                continue
            A = 1 | (1 << a) | (1 << c)
            sgn = 1 if a < c else -1
            eA = eta[0] * eta[a] * eta[c]
            for kind in 'OE':
                L = _KS ^ A ^ (FULL if kind == 'E' else 0)
                tix = TYPE_TAB[kind][POP, POP[L]]
                ok = (tix >= 0) & keep & keep[L]
                K_, L_, t_ = _KS[ok], L[ok], tix[ok]
                lowK = POP[K_] < POP[L_]
                lo, hi = np.where(lowK, K_, L_), np.where(lowK, L_, K_)
                if kind == 'O':
                    cc = vsign(cl, np.full_like(lo, A), lo) * vsign(cl, hi, hi)
                else:
                    cc = cl.sign(FULL, A) * vsign(cl, np.full_like(lo, FULL ^ A), lo) * vsign(cl, hi, hi)
                v = sgn * eA * cc * p[t_]
                nz = v != 0
                rows.append((c - 1) * DIMG + L_[nz])
                cols.append((a - 1) * DIMG + K_[nz])
                vals.append(v[nz])
    return np.concatenate(rows), np.concatenate(cols), np.concatenate(vals)


def full_rank(eta, p, grades=None):
    r, c, v = beta_entries(eta, p, grades)
    return sparse_rank(r, c, v, 13 * DIMG, 13 * DIMG)[0]


# ====================================================================== Spin(13) Schur blocks
class Forms:
    """forms on W (bitmasks); vectors in W (x) Lambda W as dicts (a, mask) -> int"""

    def __init__(self, eta):
        self.eta = list(eta)

    def etaM(self, M):
        s = 1
        for i in range(N):
            if M >> i & 1:
                s *= self.eta[i]
        return s

    @staticmethod
    def wedge1(a, J):
        return (0, 0) if J >> a & 1 else ((-1) ** pc(J & ((1 << a) - 1)), J | (1 << a))

    @staticmethod
    def iota(a, J):
        return ((-1) ** pc(J & ((1 << a) - 1)), J ^ (1 << a)) if J >> a & 1 else (0, 0)

    def starW(self, J):
        Jc = WFULL ^ J
        seq = [i for i in WB if J >> i & 1] + [i for i in WB if Jc >> i & 1]
        s = 1
        for x in range(len(seq)):
            for y in range(x + 1, len(seq)):
                if seq[x] > seq[y]:
                    s = -s
        return self.etaM(J) * s, Jc

    def phi_plus(self, J):      # w -> sum_a eta^aa e^a (x) e^a ^ w
        out = {}
        for a in WB:
            s, M = self.wedge1(a, J)
            if s:
                out[(a, M)] = out.get((a, M), 0) + self.eta[a] * s
        return out

    def phi_minus(self, J):     # w -> sum_a e^a (x) iota_{e_a} w
        out = {}
        for a in WB:
            s, M = self.iota(a, J)
            if s:
                out[(a, M)] = out.get((a, M), 0) + s
        return out

    def star_vec(self, v):
        out = {}
        for (a, M), c in v.items():
            s, Mc = self.starW(M)
            out[(a, Mc)] = out.get((a, Mc), 0) + s * c
        return out

    def quantize(self, v, with0):
        """W (x) Lambda^m W -> W (x) g: e^M -> gamma^M = eta^M E_M, or gamma^0 gamma^M"""
        out = {}
        for (a, M), c in v.items():
            cc = c * self.etaM(M) * (self.eta[0] if with0 else 1)
            K = (M | 1) if with0 else M
            out[(a, K)] = out.get((a, K), 0) + cc
        return out


RHOS = [('L', i) for i in range(7)] + [('H', j) for j in range(1, 7)]


def rho_name(rho):
    return f"Lambda^{rho[1]}" if rho[0] == 'L' else f"H_{rho[1]}"


def dim_rho(rho):
    k, i = rho
    return comb(13, i) if k == 'L' else 13 * comb(13, i) - comb(13, i + 1) - comb(13, i - 1)


def copies(F, rho):
    k, i = rho
    if k == 'L':
        fm = [('+', F.phi_plus), ('-', F.phi_minus),
              ('*+', lambda J: F.star_vec(F.phi_plus(J))), ('*-', lambda J: F.star_vec(F.phi_minus(J)))]
        if i == 0:
            fm = [fm[0], fm[2]]
    else:
        fm = [('u', lambda u: u), ('*u', F.star_vec)]
    out = []
    for lab, f in fm:
        for w0 in (False, True):
            out.append((lab + ('g0' if w0 else ''), lambda om, f=f, w0=w0: F.quantize(f(om), w0)))
    return out


def test_vector(F, rho, which):
    k, i = rho
    if k == 'L':
        J = [mask(range(1, 1 + i)), mask(range(14 - i, 14)), mask(range(4, 4 + i))][which]
        return J, F.etaM(J)
    a1, a2 = [(1, 2), (5, 6), (8, 9)][which]
    rest = [x for x in WB if x not in (a1, a2)]
    S = mask(rest[:i - 1] if which < 2 else (rest[len(rest) - (i - 1):] if i > 1 else []))
    u = {(a1, S | (1 << a1)): F.wedge1(a1, S)[0], (a2, S | (1 << a2)): -F.wedge1(a2, S)[0]}
    wed, con = {}, {}
    for (a, M), c in u.items():
        s, M2 = F.wedge1(a, M)
        if s:
            wed[M2] = wed.get(M2, 0) + s * c
        s, M2 = F.iota(a, M)
        if s:
            con[M2] = con.get(M2, 0) + s * c * F.eta[a]
    assert all(v == 0 for v in wed.values()) and all(v == 0 for v in con.values())   # u lies in H_j
    return u, sum(F.eta[a] * F.etaM(M) * c * c for (a, M), c in u.items())


def beta_vec(P, v, w):
    """beta(v, w) as an integer vector over the 48 types"""
    res = [0] * NT
    for (a, K), x in v.items():
        for (c, L), y in w.items():
            if a == c:
                continue
            A = 1 | (1 << a) | (1 << c)
            s = 1 if a < c else -1
            for t, val in P.eval(A, K, L).items():
                res[t] += s * val * x * y
    return res


def schur_data(eta, rho, which=0):
    F, P = Forms(eta), Psi(eta)
    cps = copies(F, rho)
    om, nrm = test_vector(F, rho, which)
    vecs = [(lab, f(om)) for lab, f in cps]
    grades = [pc(next(iter(v))[1]) for _, v in vecs]
    assert all(len({pc(K) for (_, K) in v}) == 1 for _, v in vecs)
    M = [[beta_vec(P, v, w) for _, w in vecs] for _, v in vecs]
    return [lab for lab, _ in vecs], grades, M, nrm, vecs


PSYM = sp.symbols(f'p0:{NT}')


def schur_C(eta, rho, grades_keep=None, which=0):
    labs, grades, M, nrm, _ = schur_data(eta, rho, which)
    keep = [x for x in range(len(labs)) if grades_keep is None or grades[x] in grades_keep]
    odd = [x for x in keep if grades[x] % 2]
    ev = [x for x in keep if grades[x] % 2 == 0]
    C = sp.Matrix(len(odd), len(ev), lambda i, j: sum(sp.Rational(M[odd[i]][ev[j]][t], nrm) * PSYM[t] for t in range(NT)))
    return C, [f"{labs[x]}[{grades[x]}]" for x in odd], [f"{labs[x]}[{grades[x]}]" for x in ev]


def at(C, p):
    return C.subs(dict(zip(PSYM, p)))


def h2_schur(eta, show=True):
    if show:
        print("\n=== H2. Spin(13) = Stab(e^0) reduction in (14,0): W (x) g = sum_rho m_rho rho; beta = sum_rho C_rho (x) <,>_rho ===")
        print("  copies of rho in W (x) g: label[Clifford grade]; '+': w -> sum e^a (x) e^a ^ w, '-': w -> sum e^a (x) i_a w,"
              " '*': Hodge star on W, 'g0': times gamma^0; H_j: u = e^1(x)e^(1S) - e^2(x)e^(2S)")
    total, blocks, ok_all = 0, {}, True
    P = Psi(eta)
    for rho in RHOS:
        labs, grades, M, nrm, vecs = schur_data(eta, rho)
        m = len(labs)
        total += m * dim_rho(rho)
        anti = all(M[x][y][t] == -M[y][x][t] for x in range(m) for y in range(m) for t in range(NT))
        par = all(M[x][y][t] == 0 for x in range(m) for y in range(m) for t in range(NT) if (grades[x] - grades[y]) % 2 == 0)
        cons = True
        for w in (1, 2):
            _, _, M2, n2, vecs2 = schur_data(eta, rho, w)
            cons &= all(M[x][y][t] * n2 == M2[x][y][t] * nrm for x in range(m) for y in range(m) for t in range(NT))
            if not (rho[0] == 'L' and rho[1] == 0):   # different test vectors are orthogonal
                cons &= all(v == 0 for x in range(m) for y in range(m) for v in beta_vec(P, vecs[x][1], vecs2[y][1]))
        C, ro, co = schur_C(eta, rho)
        blocks[rho] = C
        ok_all &= anti and par and cons
        if show:
            print(f"  {rho_name(rho):9s} dim {dim_rho(rho):5d}, m = {m}: odd copies {ro}, even copies {co}; antisymmetric: {anti}; "
                  f"parity-odd: {par}; Schur-consistent on 3 test vectors (+ orthogonality): {cons}")
            print(f"            C = {[[str(C[i, j]) for j in range(C.shape[1])] for i in range(C.shape[0])]}")
    # orthogonality across different rho (test vectors 0)
    orth = True
    data = {rho: schur_data(eta, rho)[4] for rho in RHOS}
    for r1, r2 in itertools.combinations(RHOS, 2):
        for _, v in data[r1]:
            for _, w in data[r2]:
                if any(beta_vec(P, v, w)):
                    orth = False
    if show:
        print(f"  beta(copy of rho, copy of rho') = 0 for all rho != rho' (all copies, test vectors 0): {orth}")
        print(f"  sum_rho m_rho dim rho = {total} (= 13 * 2^14 = {13 * DIMG}: {total == 13 * DIMG})")
    assert ok_all and orth and total == 13 * DIMG
    return blocks


def predicted_rank(blocks, p):
    return sum(2 * at(C, p).rank() * dim_rho(rho) for rho, C in blocks.items() if C.shape[0] and C.shape[1])


# ====================================================================== main
def main():
    t0 = time.time()
    rng = np.random.default_rng(20210401)
    ph = h0_basis(rng)

    print("\n=== H1. Draft-family operators in the basis (exact; G_symbols operator code) ===")
    vec = {}
    for eta, lab in ((EUC, '(14,0)'), (E77, '(7,7)')):
        for which, tw in (('WA', False), ('WpA', False), ('WC', True), ('WpC', True), ('WC', False), ('WpC', False)):
            v, nrows = fit_operator(eta, which, tw, rng)
            nm = which + (' o Ad_Gamma' if tw else '')
            if v is None:
                print(f"  {lab} {nm:14s}: no solution in S_sa ({nrows} exact equations) -> not self-adjoint")
            else:
                vec[(lab, which, tw)] = v
                print(f"  {lab} {nm:14s}: {nrows} exact equations, unique solution {fmt_vec(v)}")
    WA, WpA = vec[('(14,0)', 'WA', False)], vec[('(14,0)', 'WpA', False)]
    WCg, WpCg = vec[('(14,0)', 'WC', True)], vec[('(14,0)', 'WpC', True)]
    # entrywise identity of the assembled symbol with G_symbols' matrix (xi = e^0, kappa = 0)
    for eta, lab in ((EUC, '(14,0)'), (E77, '(7,7)')):
        for which in ('WA', 'WpA'):
            cl, E, ops = GS.shiab_ops(eta, which)
            rr, cc, vv, nr, nc = GS.symbol_entries(ops, cl, E, [1] + [0] * 13, 0)
            Mg = sps.coo_matrix((vv, (rr, cc)), shape=(nr, nc)).tocsr()
            Mg.sum_duplicates()
            Mg.eliminate_zeros()
            rowc, rowJ = np.repeat(np.arange(14), DIMG), np.tile(_KS, 14)
            scal = (((-1) ** rowc) * vsign(cl, rowJ, rowJ)).astype(np.int64)
            Mg = (sps.diags(scal, dtype=np.int64) @ Mg).tocsr()
            z0 = Mg[:DIMG].nnz == 0 and Mg[:, :DIMG].nnz == 0
            Mg = Mg[DIMG:, DIMG:14 * DIMG]
            r, c, v = beta_entries(eta, vec[(lab, which, False)])
            Mb = sps.coo_matrix((v, (r, c)), shape=(13 * DIMG, 13 * DIMG)).tocsr()
            Mb.sum_duplicates()
            D = Mb - Mg
            D.eliminate_zeros()
            print(f"  {lab} {which}: symbol of G_symbols' operator at xi = e^0 (rows rescaled by (-1)^c Tr(E_J E_J)/128) "
                  f"vanishes on xi (x) g and lands in Ann(xi (x) g): {z0}; equals beta of the fitted Psi entrywise "
                  f"({Mb.nnz} nonzeros): {D.nnz == 0}")
            assert z0 and D.nnz == 0

    blocks = h2_schur(EUC)

    print("\n=== H3. Exact determinants, discriminant, and exact full ranks of beta (need 212992) in (14,0) ===")
    facs = {}
    for rho, C in blocks.items():
        d = sp.factor(C.det())
        print(f"  det C[{rho_name(rho)}] = {d}")
        for f, e in sp.factor_list(C.det())[1]:
            facs.setdefault(sp.srepr(sp.expand(f)), (f, rho))
    deg = sum(sp.Poly(f, *PSYM).total_degree() for f, _ in facs.values())
    print(f"  distinct irreducible factors over Q: {len(facs)}; Delta = their product has degree {deg}; "
          f"none is identically zero: {all(f != 0 for f, _ in facs.values())}")
    allv = set().union(*[f.free_symbols for f, _ in facs.values()])
    print(f"  variables occurring in Delta: {len(allv)} of {NT}")
    rr = sp.symbols(f'r0:{NT}')
    for f, rho in facs.values():
        fr = sp.expand(f.subs({PSYM[t]: ph[t] * rr[t] for t in range(NT)}))
        c0 = sp.Poly(fr, *rr).coeffs()[0]
        fr = sp.expand(fr / c0)                                  # divide by a unit (power of i times a rational)
        real = all(sp.im(c) == 0 for c in sp.Poly(fr, *rr).coeffs())
        signs = set()
        for _ in range(200):
            v = {rr[t]: int(x) for t, x in enumerate(rng.integers(-4, 5, NT))}
            signs.add(sp.sign(fr.subs(v)))
            if {1, -1} <= signs:
                break
        tq = time.time()
        irr = len(sp.factor_list(f, gaussian=True)[1]) == 1 if not os.environ.get('H_SKIP_QI') else True
        print(f"  factor from {rho_name(rho):9s}: degree {sp.Poly(f, *PSYM).total_degree()}, {len(f.free_symbols):2d} variables; "
              f"irreducible over Q(i): {irr} ({time.time() - tq:.0f} s); real polynomial on the u(128) form (after a unit): "
              f"{real}; takes both signs at real integer points: {({1, -1} <= signs)}")
        assert real and irr and {1, -1} <= signs
    tests = []
    pr = [int(x) for x in rng.integers(1, 6, NT) * rng.choice([-1, 1], NT)]
    tests += [("W_A", WA), ("W'_A", WpA), ("W_A + W'_A", [a + b for a, b in zip(WA, WpA)]), ("random integer p", pr)]
    q2 = list(pr)                        # point on Q_2 = p3 p21 + p28 p32 = 0 only
    q2[IDX[('O', 11, 12)]] = 1
    q2[IDX[('O', 2, 3)]] = -q2[IDX[('E', 2, 11)]] * q2[IDX[('E', 3, 12)]]
    tests.append(("random, on det C[H_2] = 0", q2))
    C3 = blocks[('L', 3)]
    i7 = IDX[('O', 4, 5)]
    d3 = sp.expand(C3.det())
    sub = {PSYM[t]: pr[t] for t in range(NT) if t != i7}
    Acoef = sp.diff(d3, PSYM[i7]).subs(sub)
    Bcoef = d3.subs({PSYM[i7]: 0}).subs(sub)
    assert sp.degree(d3, PSYM[i7]) == 1 and Acoef != 0
    d3pt = [int(Acoef) * x for x in pr]
    d3pt[i7] = int(-Bcoef)
    assert sp.expand(d3.subs(dict(zip(PSYM, d3pt)))) == 0
    tests.append(("random, on det C[Lambda^3] = 0", d3pt))
    for lab, p in tests:
        pred = predicted_rank(blocks, p)
        zero = [rho_name(r) for r, C in blocks.items() if at(C, p).det() == 0]
        fr = full_rank(EUC, p)
        print(f"  {lab:32s}: vanishing det C_rho: {zero if zero else 'none'}; Schur-predicted rank {pred}; "
              f"exact sparse rank {fr}; agree: {pred == fr}; nondegenerate: {fr == 13 * DIMG}")
        assert pred == fr
    print(f"  (G_symbols.out, (14,0): rank Shiab(xi^.) = 122850 for W_A, W'_A, W_A+W'_A; rank sigma(delta_2) = 139234 = "
          f"122850 + 16384, as Lemma H1 predicts)")

    print("\n=== H4. Einstein output on Levi-Civita curvature (exact rationals; Y of grade 1, 13 and 3, 5, 9, 11) ===")
    from shiab_checks import random_curvature, ricci
    from fractions import Fraction
    R = random_curvature(N, np.random.default_rng(1)).astype(np.int64)
    Ric, s = ricci(R.astype(float), np.eye(N))
    Ric, s = Ric.astype(np.int64), int(round(s))
    P = Psi(EUC)
    F2 = {}
    for a, b in itertools.combinations(range(N), 2):
        F2[(a, b)] = {mask((c, d)): int(R[a, b, c, d]) for c, d in itertools.combinations(range(N), 2) if R[a, b, c, d]}

    def out(Ymask, m):
        """vector over types of sum_{a<b} Psi(e^m ^ e^a ^ e^b; E_Y, 2 F_ab)"""
        res = [0] * NT
        for (a, b), Fab in F2.items():
            if m in (a, b):
                continue
            A = mask((m, a, b))
            sg = 1
            for x, y in ((m, a), (m, b), (a, b)):
                if x > y:
                    sg = -sg
            for L in (A ^ Ymask, A ^ Ymask ^ FULL):
                if L in Fab:
                    for t, v in P.eval(A, Ymask, L).items():
                        res[t] += sg * v * Fab[L]
        return res
    G = Ric - Fraction(1, 2) * s * np.eye(N, dtype=object)
    clE = Cliff(EUC)
    # Y = gamma_n (grade 1) and Y = gamma_n Gamma = sign(n, full) E_(full - n) (grade 13)
    for grade_lab, Yof, sY in (("1", lambda n: 1 << n, lambda n: 1),
                               ("13", lambda n: FULL ^ (1 << n), lambda n: clE.sign(1 << n, FULL))):
        K = {t: [[sY(n) * out(Yof(n), m)[t] for n in range(N)] for m in range(N)] for t in range(NT)}
        nz = [t for t in range(NT) if any(K[t][m][n] for m in range(N) for n in range(N))]
        for t in nz:
            m0, n0 = next((m, n) for m in range(N) for n in range(N) if G[m][n] != 0)
            lam = Fraction(K[t][m0][n0]) / G[m0][n0]
            exact = all(Fraction(K[t][m][n]) == lam * G[m][n] for m in range(N) for n in range(N))
            print(f"  output grade {grade_lab:2s}: only type {NAME[t]} contributes ({len(nz)} type(s)); "
                  f"Y = {'gamma_n' if grade_lab == '1' else 'gamma_n Gamma'}: K_mn = lambda (Ric - Scal/2 eta)_mn exactly: {exact}, "
                  f"lambda = {lam} per unit p_t")
            assert exact
    zero_other = True
    for g in (3, 5, 9, 11):
        for _ in range(3):
            Y = mask(tuple(int(x) for x in rng.choice(N, g, replace=False)))
            zero_other &= all(not any(out(Y, m)) for m in range(N))
    print(f"  outputs of grade 3, 5, 9, 11 vanish for all 48 types (3 random monomials each, all m): {zero_other}")
    assert zero_other

    print("\n=== H5. The draft's self-adjoint family span(W_A, W'_A) and the twisted family ===")
    x, y, z, w = sp.symbols('x y z w')
    pen = [x * a + y * b for a, b in zip(WA, WpA)]
    fam = [x * a + y * b + z * c + w * d for a, b, c, d in zip(WA, WpA, WCg, WpCg)]
    sos_ok = True
    for rho, C in blocks.items():
        dp = sp.factor(at(C, pen).det())
        df = sp.factor(at(C, fam).det())
        line = f"  {rho_name(rho):9s}: det on x W_A + y W'_A = {dp};  on + z W_C Ad_G + w W'_C Ad_G = {df}"
        if rho[0] == 'L' and rho[1] >= 2:
            e = sp.expand(at(C, fam).det())
            q4, q2_ = (x ** 2 + y ** 2 - z ** 2 - w ** 2) ** 2, (x * z + y * w) ** 2
            cq4 = sp.Poly(e, x, y, z, w).coeff_monomial(x ** 4)
            cq2 = sp.Poly(sp.expand(e - cq4 * q4), x, y, z, w).coeff_monomial(x ** 2 * z ** 2)
            ok = sp.expand(e - cq4 * q4 - cq2 * q2_) == 0
            sos_ok &= ok and cq4 * cq2 > 0
            line += f"\n             = {sp.factor(cq4)} (x^2+y^2-z^2-w^2)^2 + {sp.factor(cq2)} (xz+yw)^2 exactly: {ok}"
        print(line)
    print(f"  sum-of-squares identities hold with coefficients of one sign: {sos_ok}")
    assert sos_ok
    # u(128)-real coordinates: Phi^E = Sum e^abc (x) gamma_abc = -i Phi(code)  => coefficients (x,y,z,w) = -i (x',y',z',w')
    xr, yr, zr, wr = sp.symbols("x' y' z' w'", real=True)
    famr = [sp.expand(f.subs({x: -sp.I * xr, y: -sp.I * yr, z: -sp.I * zr, w: -sp.I * wr})) for f in fam]
    realok = True
    for f, ph_t in zip(famr, ph):
        r_ = sp.expand(f / ph_t)
        realok &= sp.im(r_) == 0
    print(f"  with the u(128)-valued Phi^E = sum e^abc (x) gamma_abc (= -i Phi of G_symbols), the operators "
          f"x' W_A + y' W'_A + z' W_C Ad_G + w' W'_C Ad_G (x',y',z',w' real) are u(128)-real: {realok}; "
          f"their determinants are, up to units, (x'z'+y'w') and a(x'^2+y'^2-z'^2-w'^2)^2 + b(x'z'+y'w')^2 with ab > 0")
    assert realok
    sharp = [sp.expand(f.subs({xr: 1, yr: 0, zr: 1, wr: 0})) for f in famr]   # W_A^E + W_C^E o Ad_Gamma
    p_sharp = [int(sp.expand(f / ph_t)) * int(sp.re(ph_t * ph_t)) for f, ph_t in zip(sharp, ph)]
    dets = {rho_name(r): at(C, sharp).det() for r, C in blocks.items()}
    fr = full_rank(EUC, p_sharp)
    print(f"  example  Shiab# = W_A^E + W_C^E o Ad_Gamma  (u(128)-real; Ad_Gamma X = Gamma X Gamma^-1 = grade involution):")
    print(f"     det C_rho = {dets}")
    print(f"     exact rank of beta in the u(128) basis = {fr} of {13 * DIMG}: nondegenerate {fr == 13 * DIMG}")
    # assemble sigma(delta_2)(e^0), kappa = 1, on Lambda^1 (x) g (+) g -> Lambda^13 (x) g (operator rows)
    cl = Cliff(EUC)
    r, c, v = beta_entries(EUC, p_sharp)
    rc = r // DIMG + 1
    v = v * ((-1) ** rc) * vsign(cl, r % DIMG, r % DIMG)              # operator matrix (rows: 13-form missing c)
    rows = np.concatenate([DIMG + r, _KS])
    cols = np.concatenate([DIMG + c, 14 * DIMG + _KS])
    Es = GS.Ext(EUC)
    vals = np.concatenate([v, np.full(DIMG, -1 * Es.star_coeff((0,))[1], dtype=np.int64)])
    rk2 = sparse_rank(rows, cols, vals, 14 * DIMG, 15 * DIMG)[0]
    print(f"     rank sigma(delta_2)(e^0), kappa = 1 (true or printed Upsilon, T = 0) = {rk2} of {14 * DIMG}: "
          f"symbol sequence exact: {rk2 == 14 * DIMG}")
    assert fr == 13 * DIMG and rk2 == 14 * DIMG
    print(f"     Einstein output: p(O1,2) = {sharp[IDX[('O', 1, 2)]]}, p(E2,13) = {sharp[IDX[('E', 2, 13)]]} "
          f"(same as W_A^E: Shiab#(F_R) = W_A^E(F_R), the exact Einstein tensor up to the factor i of the real form)")

    print("\n=== H6. so(64,64) / so(128): grades 1,2,5,6,9,10,13,14 (need 13 * 8128 = 105664) ===")
    so_t = [t for t in range(NT) if TYPES[t][1] in GSO and TYPES[t][2] in GSO]
    print(f"  types inside so: {len(so_t)} ({', '.join(NAME[t] for t in so_t)}) = dim S_sa(so(64,64)) = 12 (D1): {len(so_t) == 12}")
    gen = lower = predWA = 0
    pso = [pr[t] if t in so_t else 0 for t in range(NT)]
    for rho in RHOS:
        C, ro, co = schur_C(EUC, rho, GSO)
        rk = at(C, pso).rank() if C.shape[0] and C.shape[1] else 0
        mx = min(C.shape)
        gen += 2 * rk * dim_rho(rho)
        lower += abs(C.shape[0] - C.shape[1]) * dim_rho(rho)
        predWA += 2 * (at(C, WA).rank() if mx else 0) * dim_rho(rho)
        print(f"  {rho_name(rho):9s}: odd {ro} even {co}: C = {C.tolist()}; generic rank {rk} (max possible {mx})")
    rso = full_rank(EUC, pso, GSO)
    rwa = full_rank(EUC, WA, GSO)
    print(f"  corank >= {lower} from unequal odd/even multiplicities alone; Schur generic rank {gen}; exact sparse rank at "
          f"a random point {rso}, at W_A {rwa} (Schur {predWA}); maximal rank over all 12 dimensions = {gen}, i.e. deficit >= {13 * 8128 - gen} for every member")
    assert gen == rso and predWA == rwa and gen < 13 * 8128

    print("\n=== H7. (7,7), non-null xi = e^0 (comparison with G11; the null cone is excluded by G11.1) ===")
    b77 = h2_schur(E77, show=False)
    WA7, WpA7 = vec[('(7,7)', 'WA', False)], vec[('(7,7)', 'WpA', False)]
    pen7 = [x * a + y * b for a, b in zip(WA7, WpA7)]
    print("  det C_rho on x W_A + y W'_A: " + "; ".join(f"{rho_name(r)}: {sp.factor(at(C, pen7).det())}" for r, C in b77.items()))
    for lab, p in (("W_A + W'_A", [a + b for a, b in zip(WA7, WpA7)]), ("W_A - W'_A", [a - b for a, b in zip(WA7, WpA7)]),
                   ("W_A", WA7), ("random integer p", pr)):
        pred = predicted_rank(b77, p)
        extra = f"; exact sparse rank {full_rank(E77, p)}" if lab.startswith("random") else ""
        print(f"  {lab:18s}: Schur-predicted rank {pred}{extra}")
    print("  (G_symbols.out (7,7) non-null: 106496 for W_A +- W'_A, 122850 for W_A)")
    print(f"\nruntime {time.time() - t0:.0f} s")


if __name__ == '__main__':
    main()
