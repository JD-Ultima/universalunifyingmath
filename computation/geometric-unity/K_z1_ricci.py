"""GU round 7 (K): a complete proof of Theorem F4(c) (Ricci blocks of the Zorro metric Z1 along s(X)) for every n.
Write-up: formalization/geometric-unity/K-z1-ricci-general-n.md.  Record: metrics/geometric-unity/K-z1-ricci.yaml.

Setting (F, R-Z): X^n with metric g, Y = Met(X), coordinates z = (x^mu, y_A), A = (i <= j), Dy = dy - Psi dx with
Psi_{(ij) mu} = Gamma^k_{mu i} y_kj + Gamma^k_{mu j} y_ik (Levi-Civita of g induced on Sym^2 T*X),
    Z1 = cH y(dx,dx) + cV G_lambda|_y(Dy,Dy),   Z2 = cH g(dx,dx) + cV G_lambda|_y(Dy,Dy),   delta = Z1 - Z2,
G_lambda|_y(h,k) = tr(y^-1 h y^-1 k) - lambda tr(y^-1 h) tr(y^-1 k).  Curvature conventions of F (unit sphere Ric = (n-1)g).
The proof (K-z1-ricci-general-n.md) is: Ric(Z1) = Ric(Z2) + dRic, where at a point p of s(X)
  (K1) dRic is a universal expression in (h, nabla delta, nabla^2 delta) at p (h = Z2, nabla = nabla^h, delta(p) = 0);
  (K2) at p, nabla delta and nabla^2 delta have only the components
         (nabla_U delta)(X,Y) = cH U(X,Y),   (nabla^2_{U,W} delta)(X,Y) = (cH/2)(U g^-1 W + W g^-1 U)(X,Y);
  (K3) evaluating (K1) on these data gives, for every n, dRic_HH = cH/(4 cV (n lam - 1)) g, dRic_HV = 0,
       dRic_VV = -(1/4) G_0  (exact O(g)-invariant diagram calculus with a symbolic loop weight n);
  then F3 (Z2) gives F4(c).
Sections (all exact: Fractions or sympy rational-function fields; fixed seed; one core):
  K0  controls of the independent curvature routine (spheres / hyperbolic spaces).
  K1  the variation formula on random exact 2-jets, dimensions m = 2..5 (definite and indefinite).
  K2  the jets of delta and the Z2 facts used, on curved bases in non-normal coordinates (n = 2, 3, 4):
      nabla pi_* = 0 at p, nabla K at p, nabla delta and nabla^2 delta at p; fibre Christoffel symbols (lambda-free).
  K3  general-n algebra: Brauer-diagram evaluation of (K1) on the data of (K2) with loop weight n (symbolic);
      cross-check against explicit tensors for n = 1..4; F4(a), F4(b) and F6 consistency.
  K4  independent brute-force Ricci tensor of Z1 (and Z2) at a point of s(X): n = 2..6, two or three signatures each,
      random curved bases in non-normal coordinates; also the chain dRic(brute force) = (K1) on brute-force (K2) data.
  K5  the flat model of the F4(c) sketch at a point y0 != eta (homogeneity), n = 3, 4.
Run from computation/geometric-unity:  python3 K_z1_ricci.py > K_z1_ricci.out 2>&1
This script does not import F_jets.py or any other repository code.
"""
import itertools
import random
import sys
import time
from collections import defaultdict
from fractions import Fraction as Fr

import numpy as np
from sympy import QQ, Rational
from sympy.polys.fields import field

T0 = time.time()
RNG = random.Random(20261030)
CHECKS = []


def P(*a):
    print(*a)
    sys.stdout.flush()


def check(label, ok):
    ok = bool(ok)
    CHECKS.append((label, ok))
    P(f"  [{ok}] {label}")
    assert ok, label


# ====================================================================== exact scalars
class Num:
    """Exact scalar context: Python Fractions (K = None) or a sympy rational-function field K."""

    def __init__(self, K=None):
        self.K = K
        self.zero = self.c(0)
        self.one = self.c(1)

    def c(self, x, y=1):
        if self.K is not None and not isinstance(x, (int, Fr)):
            return x / y if y != 1 else x          # already a field element
        q = Fr(x) / Fr(y)
        if self.K is None:
            return q
        return self.K(Rational(q.numerator, q.denominator))


def oarr(shape, val=0):
    a = np.empty(shape, dtype=object)
    a.fill(val)
    return a


def mat_inv(M, C):
    """Exact inverse (Gauss-Jordan) of a square object array."""
    m = M.shape[0]
    A = [[M[i, j] for j in range(m)] + [C.one if i == j else C.zero for j in range(m)] for i in range(m)]
    for col in range(m):
        piv = next(r for r in range(col, m) if A[r][col] != 0)
        A[col], A[piv] = A[piv], A[col]
        p = A[col][col]
        A[col] = [x / p for x in A[col]]
        rc = A[col]
        nz = [j for j in range(2 * m) if rc[j] != 0]
        for r in range(m):
            if r != col and A[r][col] != 0:
                f = A[r][col]
                row = A[r]
                for j in nz:
                    row[j] = row[j] - f * rc[j]
    out = oarr((m, m))
    for i in range(m):
        for j in range(m):
            out[i, j] = A[i][m + j]
    return out


def is_zero(a):
    return all(x == 0 for x in np.asarray(a, dtype=object).flat)


def eq(a, b):
    return is_zero(np.asarray(a, dtype=object) - np.asarray(b, dtype=object))


def ein(spec, *ops):
    return np.einsum(spec, *ops, optimize=True)


# ====================================================================== Ricci tensor from a 2-jet (own routine)
def christoffel(g0, dg, C):
    """Gamma^a_{db} at the origin from g0 and dg[c,a,b] = d_c g_ab."""
    gi = mat_inv(g0, C)
    G1 = (np.einsum('deb->edb', dg) + np.einsum('bed->edb', dg) - dg) * C.c(1, 2)   # G1[e,d,b] = Gamma_{e,db}
    return gi, G1, ein('ae,edb->adb', gi, G1)


def ricci_2jet(g0, dg, ddg, C):
    """Ric_{db} = d_a Gamma^a_{db} - d_d Gamma^a_{ab} + Gamma^a_{ae} Gamma^e_{db} - Gamma^a_{de} Gamma^e_{ab}
       at the origin, from g0[a,b], dg[c,a,b] = d_c g_ab, ddg[c,d,a,b] = d_c d_d g_ab (any coordinates)."""
    gi, G1, Gam = christoffel(g0, dg, C)
    M = ein('ap,cpq->caq', gi, dg)                                   # M_c = g^-1 d_c g
    w = -ein('q,qe->e', np.einsum('aaq->q', M), gi)                  # w^e = d_a g^{ae}
    half = C.c(1, 2)
    t1 = ein('e,edb->db', w, G1) + (ein('ae,adeb->db', gi, ddg) + ein('ae,abed->db', gi, ddg)
                                    - ein('ae,aedb->db', gi, ddg)) * half
    t2 = (-ein('dpq,bqp->db', M, M) + ein('ae,dbae->db', gi, ddg)) * half
    t3 = ein('e,edb->db', np.einsum('aae->e', Gam), Gam) - ein('ade,eab->db', Gam, Gam)
    return t1 - t2 + t3, gi, Gam


# ====================================================================== random exact data
def rand_q(lo=-3, hi=3, dmax=3):
    return Fr(RNG.randint(lo, hi), RNG.randint(1, dmax))


def rand_sym(m, C, **kw):
    a = oarr((m, m))
    for i in range(m):
        for j in range(i, m):
            a[i, j] = a[j, i] = C.c(rand_q(**kw))
    return a


def rand_g0(eps, nondiag, C):
    m = len(eps)
    D = oarr((m, m), C.zero)
    for i in range(m):
        D[i, i] = C.c(eps[i])
    if not nondiag:
        return D
    U = oarr((m, m), C.zero)                                        # unimodular: det g0 = prod(eps)
    for i in range(m):
        U[i, i] = C.one
        for j in range(i + 1, m):
            U[i, j] = C.c(RNG.randint(-1, 1))
    return ein('ki,kl,lj->ij', U, D, U)


def rand_jets(m, C, sym_last=True):
    """Random exact symmetric first and second derivative arrays d1[c,a,b], d2[c,d,a,b] (Schwarz-symmetric)."""
    d1 = oarr((m, m, m))
    for c in range(m):
        d1[c] = rand_sym(m, C)
    d2 = oarr((m, m, m, m))
    for c in range(m):
        for d in range(c, m):
            s = rand_sym(m, C)
            d2[c, d] = s
            d2[d, c] = s
    return d1, d2


# ====================================================================== the variation formula (K1), as data
# Delta Ric_{bc} at a point where delta = 0; D1[a,b,c] = nabla_a delta_bc, D2[a,b,c,d] = nabla_a nabla_b delta_cd,
# T[d,b,c] = D1[b,d,c] + D1[c,d,b] - D1[d,b,c]; every pair listed in the third entry is contracted with h^-1.
VAR_TERMS = [
    (Fr(1, 2), [('D2', 'abdc')], ['ad']),
    (Fr(1, 2), [('D2', 'acbd')], ['ad']),
    (Fr(-1, 2), [('D2', 'adbc')], ['ad']),
    (Fr(-1, 2), [('D2', 'bcad')], ['ad']),
    (Fr(-1, 2), [('D1', 'apq'), ('T', 'dbc')], ['ap', 'dq']),
    (Fr(1, 2), [('D1', 'bpq'), ('D1', 'cad')], ['ap', 'dq']),
    (Fr(1, 4), [('T', 'pad'), ('T', 'qbc')], ['ap', 'dq']),
    (Fr(-1, 4), [('T', 'pbd'), ('T', 'qac')], ['ap', 'dq']),
]


def T_from_D1(D1):
    return np.einsum('bdc->dbc', D1) + np.einsum('cdb->dbc', D1) - D1


def var_numeric(hinv, D1, D2, C, terms=VAR_TERMS, out='bc'):
    """Evaluate the variation formula with explicit arrays (any basis)."""
    tens = dict(D1=D1, D2=D2, T=T_from_D1(D1))
    acc = None
    for coef, facs, pairs in terms:
        subs = [lab for _, lab in facs] + list(pairs)
        ops = [tens[nm] for nm, _ in facs] + [hinv] * len(pairs)
        val = ein(','.join(subs) + '->' + out, *ops) * C.c(coef)
        acc = val if acc is None else acc + val
    return acc


def nabla_delta(d1, d2, Gam, nh=None):
    """nabla delta and nabla^2 delta at a point where delta = 0, from d1 = d delta, d2 = dd delta and Gamma there:
       nabla_a nabla_b delta_cd = d_a d_b delta_cd - Gamma^e_bc d_a delta_ed - Gamma^e_bd d_a delta_ce
                                  - Gamma^e_ab d_e delta_cd - Gamma^e_ac d_b delta_ed - Gamma^e_ad d_b delta_ce.
       If nh is given, delta_cd is assumed to vanish identically unless c, d < nh (sparse evaluation)."""
    D1 = d1.copy()
    if nh is None:
        D2 = (d2 - ein('ebc,aed->abcd', Gam, d1) - ein('ebd,ace->abcd', Gam, d1) - ein('eab,ecd->abcd', Gam, d1)
              - ein('eac,bed->abcd', Gam, d1) - ein('ead,bce->abcd', Gam, d1))
        return D1, D2
    assert is_zero(d1[:, nh:, :]) and is_zero(d2[:, :, nh:, :])
    h1, Gh = d1[:, :nh, :nh], Gam[:nh]
    D2 = d2.copy()
    D2[:, :, :, :nh] = D2[:, :, :, :nh] - ein('ebc,aed->abcd', Gh, h1) - ein('eac,bed->abcd', Gh, h1)
    D2[:, :, :nh, :] = D2[:, :, :nh, :] - ein('ebd,ace->abcd', Gh, h1) - ein('ead,bce->abcd', Gh, h1)
    D2[:, :, :nh, :nh] = D2[:, :, :nh, :nh] - ein('eab,ecd->abcd', Gam, h1)
    return D1, D2


# ====================================================================== Y = Met(X): jets of the Zorro metrics
def sym_index(n):
    return [(i, j) for i in range(n) for j in range(i, n)]


def basis_E(n, C):
    idx = sym_index(n)
    E = oarr((len(idx), n, n), C.zero)
    for A, (i, j) in enumerate(idx):
        E[A, i, j] = C.one
        E[A, j, i] = C.one
    return idx, E


class Base:
    """Base metric g(x) = g0 + g1[r] x^r + (1/2) g2[r,s] x^r x^s + (1/6) g3[r,s,t] x^r x^s x^t (exact).
       Stores Gamma, d Gamma, dd Gamma at 0 (Gam[k,mu,i] = Gamma^k_{mu i}) and Ric_g(0)."""

    def __init__(self, g0, g1, g2, g3, C):
        n = g0.shape[0]
        self.n, self.C, self.g0, self.g1, self.g2, self.g3 = n, C, g0, g1, g2, g3
        gi = mat_inv(g0, C)
        half = C.c(1, 2)
        L0 = (np.einsum('ilj->lij', g1) + np.einsum('jli->lij', g1) - np.einsum('lij->lij', g1)) * half  # Gamma_{l,ij}
        L1 = (np.einsum('rilj->rlij', g2) + np.einsum('rjli->rlij', g2) - g2) * half                      # d_r Gamma_{l,ij}
        L2 = (np.einsum('rsilj->rslij', g3) + np.einsum('rsjli->rslij', g3) - g3) * half
        dgi = -ein('ka,rab,bl->rkl', gi, g1, gi)
        ddgi = (ein('ka,rab,bc,scd,dl->rskl', gi, g1, gi, g1, gi) + ein('ka,sab,bc,rcd,dl->rskl', gi, g1, gi, g1, gi)
                - ein('ka,rsab,bl->rskl', gi, g2, gi))
        self.Gam = ein('kl,lij->kij', gi, L0)
        self.dGam = ein('rkl,lij->rkij', dgi, L0) + ein('kl,rlij->rkij', gi, L1)
        self.ddGam = (ein('rskl,lij->rskij', ddgi, L0) + ein('rkl,slij->rskij', dgi, L1)
                      + ein('skl,rlij->rskij', dgi, L1) + ein('kl,rslij->rskij', gi, L2))
        self.Ric, self.gi, _ = ricci_2jet(g0, g1, g2, C)
        self.Rs = ein('ab,ab->', self.gi, self.Ric)
        # second route to Ric(0) from the Christoffel jets, and nabla Ric(0) (to show tests are not on parallel Ricci)
        G, dG, ddG = self.Gam, self.dGam, self.ddGam
        self.Ric_b = (np.einsum('aadb->db', dG) - np.einsum('daab->db', dG) + ein('aae,edb->db', G, G)
                      - ein('ade,eab->db', G, G))
        dRic = (np.einsum('raadb->rdb', ddG) - np.einsum('rdaab->rdb', ddG) + ein('raae,edb->rdb', dG, G)
                + ein('aae,redb->rdb', G, dG) - ein('rade,eab->rdb', dG, G) - ein('ade,reab->rdb', G, dG))
        self.nablaRic = dRic - ein('erd,eb->rdb', G, self.Ric) - ein('erb,de->rdb', G, self.Ric)


def random_base(n, eps, nondiag, C):
    g0 = rand_g0(eps, nondiag, C)
    g1 = oarr((n, n, n))
    for r in range(n):
        g1[r] = rand_sym(n, C)
    g2 = oarr((n, n, n, n))
    for r, s in itertools.combinations_with_replacement(range(n), 2):
        a = rand_sym(n, C)
        g2[r, s] = a
        g2[s, r] = a
    g3 = oarr((n, n, n, n, n))
    for r, s, t in itertools.combinations_with_replacement(range(n), 3):
        a = rand_sym(n, C)
        for perm in set(itertools.permutations((r, s, t))):
            g3[perm] = a
    return Base(g0, g1, g2, g3, C)


def flat_base(g0, C):
    n = g0.shape[0]
    return Base(g0, oarr((n, n, n), C.zero), oarr((n, n, n, n), C.zero), oarr((n, n, n, n, n), C.zero), C)


class YJets:
    """2-jets at z0 = (x = 0, y = y0) of the building blocks of the Zorro metrics, in coordinates (x^mu, y_A).
       Directions c = 0..n-1 are x^mu, c = n + A are y_A.  For a matrix-valued function F:
       F0 = F(z0), F1[c] = d_c F, F2[c,d] = d_c d_d F."""

    def __init__(self, base, y0=None, order=2):
        C, n = base.C, base.n
        self.base, self.C, self.n = base, C, n
        y0 = base.g0 if y0 is None else y0
        self.y0 = y0
        idx, E = basis_E(n, C)
        self.idx, self.E = idx, E
        N = len(idx)
        m = n + N
        self.N, self.m = N, m
        I, J = [i for i, _ in idx], [j for _, j in idx]

        def psi(Gam, y):        # Psi[A, mu] for Gamma^k_{mu i} = Gam[k, mu, i] and a symmetric y
            full = ein('kmi,kj->ijm', Gam, y) + ein('kmj,ik->ijm', Gam, y)
            return full[I, J, :]

        # Psi
        self.Psi0 = psi(base.Gam, y0)
        P1 = oarr((m, N, n), C.zero)
        P2 = oarr((m, m, N, n), C.zero)
        for r in range(n):
            P1[r] = psi(base.dGam[r], y0)
            for s in range(n):
                P2[r, s] = psi(base.ddGam[r, s], y0)
            for B in range(N):
                v = psi(base.dGam[r], E[B])
                P2[r, n + B] = v
                P2[n + B, r] = v
        for B in range(N):
            P1[n + B] = psi(base.Gam, E[B])
        self.Psi1, self.Psi2 = P1, P2
        # DeWitt pieces G0part(h,k) = tr(W h W k), TT(h,k) = tr(W h) tr(W k), W = y^-1; e_A = W E_A
        W = mat_inv(y0, C)
        e = ein('ij,Ajk->Aik', W, E)
        t = np.einsum('Aii->A', e)
        F2 = ein('Aij,Bji->AB', e, e)
        PP = ein('Cij,Djk->CDik', e, e)
        F3 = ein('CAij,Bji->CAB', PP, e)
        F4 = ein('CDij,ABji->CDAB', PP, PP) if order >= 2 else oarr((N, N, N, N), C.zero)
        d1 = -(F3 + np.einsum('ACB->CAB', F3))
        d2 = (F4 + np.einsum('DCAB->CDAB', F4) + np.einsum('CADB->CDAB', F4) + np.einsum('DACB->CDAB', F4)
              + np.einsum('ACDB->CDAB', F4) + np.einsum('ADCB->CDAB', F4))
        T0_ = ein('A,B->AB', t, t)
        T1 = -(ein('CA,B->CAB', F2, t) + ein('A,CB->CAB', t, F2))
        T2 = (ein('DCA,B->CDAB', F3, t) + ein('CDA,B->CDAB', F3, t) + ein('CA,DB->CDAB', F2, F2)
              + ein('DA,CB->CDAB', F2, F2) + ein('A,DCB->CDAB', t, F3) + ein('A,CDB->CDAB', t, F3))
        if order < 2:
            d2, T2 = F4, F4
        self.Gv = {}
        for key, (v0, v1, v2) in (('G0', (F2, d1, d2)), ('TT', (T0_, T1, T2))):
            G1 = oarr((m, N, N), C.zero)
            G2 = oarr((m, m, N, N), C.zero)
            G1[n:] = v1
            G2[n:, n:] = v2
            self.Gv[key] = (v0, G1, G2)
        # R = Gv Psi and S = Psi^T Gv Psi, for each DeWitt piece
        self.R, self.S = {}, {}
        for key, (G0_, G1, G2) in self.Gv.items():
            R0 = ein('AB,Bm->Am', G0_, self.Psi0)
            R1 = ein('cAB,Bm->cAm', G1, self.Psi0) + ein('AB,cBm->cAm', G0_, P1)
            S0 = ein('Am,An->mn', self.Psi0, R0)
            S1 = ein('cAm,An->cmn', P1, R0) + ein('Am,cAn->cmn', self.Psi0, R1)
            if order >= 2:
                R2 = (ein('cdAB,Bm->cdAm', G2, self.Psi0) + ein('cAB,dBm->cdAm', G1, P1)
                      + ein('dAB,cBm->cdAm', G1, P1) + ein('AB,cdBm->cdAm', G0_, P2))
                S2 = (ein('cdAm,An->cdmn', P2, R0) + ein('cAm,dAn->cdmn', P1, R1) + ein('dAm,cAn->cdmn', P1, R1)
                      + ein('Am,cdAn->cdmn', self.Psi0, R2))
            else:                                         # first-order mode: second derivatives not computed
                R2, S2 = oarr((m, m, N, n), C.zero), oarr((m, m, n, n), C.zero)
            self.R[key] = (R0, R1, R2)
            self.S[key] = (S0, S1, S2)
        # horizontal blocks: Z1 -> y (tautological), Z2 -> g(x)
        Y1 = oarr((m, n, n), C.zero)
        for B in range(N):
            Y1[n + B] = E[B]
        H1 = oarr((m, n, n), C.zero)
        H2 = oarr((m, m, n, n), C.zero)
        H1[:n] = base.g1
        H2[:n, :n] = base.g2
        self.YH = {'Z1': (y0, Y1, oarr((m, m, n, n), C.zero)), 'Z2': (base.g0, H1, H2)}

    def metric(self, variant, lam, cH, cV):
        """2-jet (G0, dG, ddG) of Z1 / Z2 at z0:  G_HH = cH YH + cV Psi^T Gv Psi,  G_HV = -cV (Gv Psi)^T,
           G_VV = cV Gv, with Gv = G0part - lam TT."""
        C, n, N, m = self.C, self.n, self.N, self.m
        out = []
        for k in range(3):
            Gv = self.Gv['G0'][k] - self.Gv['TT'][k] * lam
            R = self.R['G0'][k] - self.R['TT'][k] * lam
            S = self.S['G0'][k] - self.S['TT'][k] * lam
            YH = self.YH[variant][k]
            shape = [(), (m,), (m, m)][k]
            G = oarr(shape + (m, m), C.zero)
            G[..., :n, :n] = YH * cH + S * cV
            RT = np.swapaxes(R, -1, -2)
            G[..., :n, n:] = RT * (-cV)
            G[..., n:, :n] = R * (-cV)
            G[..., n:, n:] = Gv * cV
            out.append(G)
        return out

    def frame(self):
        """Adapted frame at z0: e_mu = d_mu + Psi0[A,mu] d_A (horizontal lifts), e_A = d_A.  Rows = frame vectors."""
        C, n, N, m = self.C, self.n, self.N, self.m
        F = oarr((m, m), C.zero)
        for i in range(m):
            F[i, i] = C.one
        F[:n, n:] = self.Psi0.T
        return F


def to_frame(T, F):
    """Components of a covariant tensor (array with k equal axes of size m) in the frame F."""
    k = T.ndim
    letters = 'abcdefgh'[:k]
    caps = 'ABCDEFGH'[:k]
    spec = ','.join(f'{caps[i]}{letters[i]}' for i in range(k)) + ',' + letters + '->' + caps
    return ein(spec, *([F] * k), T)


def dewitt(y, U, Wm, lam, C):
    yi = mat_inv(y, C)
    a, b = ein('ij,jk->ik', yi, U), ein('ij,jk->ik', yi, Wm)
    return ein('ij,ji->', a, b) - lam * np.einsum('ii->', a) * np.einsum('ii->', b)


# ====================================================================== K3 engine: O(g)-invariant strand diagrams
# A tensor on T_pY = H + V (H = T_xX, V = Sym^2 T*_xX) built from g alone is a linear combination of strand
# diagrams: an H-slot has one endpoint (s,0), a V-slot two endpoints (s,0), (s,1) (it eats a symmetric U_ij);
# a diagram is a perfect matching of the endpoints, each strand being g or g^-1.  Contracting two slots with
# h^-1 = (cH g)^-1 + (cV G_lam)^-1 adds strands; a closed loop is tr(id) = n.  The loop weight is a parameter, so
# the same code evaluates the formula for symbolic n (the Brauer calculus, valid for every n >= 1: we only use the
# map diagrams -> tensors, never its injectivity).
def dg_sym_V(tensor):
    """Average every component over the swaps of the two endpoints of each V-slot."""
    out = {}
    for types, comps in tensor.items():
        vs = [s for s, t in enumerate(types) if t == 'V']
        new = []
        for coef, strands in comps:
            for flips in itertools.product((0, 1), repeat=len(vs)):
                fl = {s: f for s, f in zip(vs, flips)}

                def mp(ep):
                    s, k = ep
                    return (s, k ^ fl.get(s, 0)) if s in fl else ep
                new.append((coef / 2 ** len(vs), tuple(tuple(sorted((mp(a), mp(b)))) for a, b in strands)))
        out[types] = new
    return out


def dg_perm(tensor, perm):
    """New tensor whose slot i is the old slot perm[i]."""
    inv = {old: new for new, old in enumerate(perm)}
    out = {}
    for types, comps in tensor.items():
        nt = tuple(types[perm[i]] for i in range(len(perm)))
        out[nt] = [(c, tuple(((inv[a[0]], a[1]), (inv[b[0]], b[1])) for a, b in st)) for c, st in comps]
    return out


def dg_lin(*pairs):
    out = defaultdict(list)
    for coef, tensor in pairs:
        for types, comps in tensor.items():
            out[types] += [(coef * c, st) for c, st in comps]
    return dict(out)


def dg_reduce(strands, free_eps):
    """Join strands; return (sorted pairing of the free endpoints, number of closed loops)."""
    adj = defaultdict(list)
    for i, (u, v) in enumerate(strands):
        adj[u].append((v, i))
        adj[v].append((u, i))
    used = [False] * len(strands)
    seen = set()
    pairs = []
    for s in sorted(free_eps):
        if s in seen:
            continue
        cur = s
        while True:
            nxt = [(w, i) for w, i in adj[cur] if not used[i]]
            if not nxt:
                break
            w, i = nxt[0]
            used[i] = True
            cur = w
        assert cur in free_eps and cur != s
        seen.update((s, cur))
        pairs.append(tuple(sorted((s, cur))))
    loops = 0
    for i in range(len(strands)):
        if used[i]:
            continue
        loops += 1
        start = strands[i][0]
        cur = start
        while True:
            nxt = [(w, j) for w, j in adj[cur] if not used[j]]
            if not nxt:
                break
            w, j = nxt[0]
            used[j] = True
            cur = w
        assert cur == start
    return tuple(sorted(pairs)), loops


def dg_evaluate(terms, tensors, free_types, C, loopw, lam, cH, cV):
    """Evaluate a list of terms (coef, [(tensor name, slot labels)], [h^-1 pairs]) on diagram tensors.
       free_types: {label: 'H' | 'V'}.  Returns {pairing of free endpoints: coefficient}."""
    mu = lam / (1 - loopw * lam)
    hinv = {'H': lambda a, b: [(1 / cH, [((a, 0), (b, 0))])],
            'V': lambda a, b: [(1 / (2 * cV), [((a, 0), (b, 0)), ((a, 1), (b, 1))]),
                               (1 / (2 * cV), [((a, 0), (b, 1)), ((a, 1), (b, 0))]),
                               (mu / cV, [((a, 0), (a, 1)), ((b, 0), (b, 1))])]}
    res = defaultdict(lambda: C.zero)
    free_eps = set()
    for lab, t in free_types.items():
        free_eps.update([(lab, 0)] if t == 'H' else [(lab, 0), (lab, 1)])
    for coef, facs, pairs in terms:
        choices = []
        for pr in pairs:
            fixed = [free_types[l] for l in pr if l in free_types]
            choices.append(fixed[:1] if fixed else ['H', 'V'])
        for assign in itertools.product(*choices):
            ty = dict(free_types)
            for pr, t in zip(pairs, assign):
                ty[pr[0]] = ty[pr[1]] = t
            lists = []
            ok = True
            for nm, labs in facs:
                key = tuple(ty[l] for l in labs)
                comps = tensors[nm].get(key)
                if not comps:
                    ok = False
                    break
                lists.append([(c, [((labs[a[0]], a[1]), (labs[b[0]], b[1])) for a, b in st]) for c, st in comps])
            if not ok:
                continue
            for pr, t in zip(pairs, assign):
                lists.append(hinv[t](pr[0], pr[1]))
            for combo in itertools.product(*lists):
                c = C.c(coef)
                strands = []
                for cc, st in combo:
                    c = c * cc
                    strands += st
                pairing, loops = dg_reduce(strands, free_eps)
                res[pairing] = res[pairing] + c * loopw ** loops
    return {k: v for k, v in res.items() if v != 0}


def dg_orbits(res, free_types):
    """Sum coefficients over the orbits of the V-slot swaps (a bilinear form on Sym^2 sees only the orbit sums)."""
    vs = [l for l, t in free_types.items() if t == 'V']
    out = defaultdict(int)
    for pairing, c in res.items():
        images = []
        for flips in itertools.product((0, 1), repeat=len(vs)):
            fl = dict(zip(vs, flips))

            def mp(ep):
                return (ep[0], ep[1] ^ fl[ep[0]]) if ep[0] in fl else ep
            images.append(tuple(sorted(tuple(sorted((mp(a), mp(b)))) for a, b in pairing)))
        out[min(images)] = out[min(images)] + c
    return {k: v for k, v in out.items() if v != 0}


def dg_data(C, cH):
    """Lemma K2 data at p as diagram tensors: D1 = nabla delta, D2 = nabla^2 delta."""
    h = C.c(1, 2)
    D1 = dg_sym_V({('V', 'H', 'H'): [(cH, (((0, 0), (1, 0)), ((0, 1), (2, 0))))]})
    D2 = dg_sym_V({('V', 'V', 'H', 'H'): [
        (cH * h, (((2, 0), (0, 0)), ((0, 1), (1, 0)), ((1, 1), (3, 0)))),     # (cH/2) (U g^-1 W)(X, Y)
        (cH * h, (((2, 0), (1, 0)), ((1, 1), (0, 0)), ((0, 1), (3, 0))))]})   # (cH/2) (W g^-1 U)(X, Y)
    T = dg_lin((1, dg_perm(D1, (1, 0, 2))), (1, dg_perm(D1, (1, 2, 0))), (-1, D1))
    return dict(D1=D1, D2=D2, T=T)


G0_REP = ((('b', 0), ('c', 0)), (('b', 1), ('c', 1)))      # orbit of G_0(U,W) = tr(g^-1 U g^-1 W)
TT_REP = ((('b', 0), ('b', 1)), (('c', 0), ('c', 1)))      # tr(g^-1 U) tr(g^-1 W)
HH_REP = ((('b', 0), ('c', 0)),)                            # g(X, Y)


# ====================================================================== main
def section(title):
    P("")
    P(title)


def ricci_frame_blocks(Ric, F, n):
    RF = to_frame(Ric, F)
    return RF[:n, :n], RF[:n, n:], RF[n:, n:]


def predicted_blocks(Y, variant, lam, cH, cV):
    """Closed forms on s(X): F3 (Z2) and F4(c) (Z1), in the adapted frame (e_mu, E_A)."""
    C, n, N, E, g0 = Y.C, Y.n, Y.N, Y.E, Y.base.g0
    nn = C.c(n)
    HH = Y.base.Ric.copy()
    VV = oarr((N, N))
    for A in range(N):
        for B in range(N):
            VV[A, B] = -nn / 4 * dewitt(g0, E[A], E[B], 1 / nn, C)
    if variant == 'Z1':
        HH = HH + g0 * (cH / (4 * cV * (nn * lam - 1)))
        for A in range(N):
            for B in range(N):
                VV[A, B] = VV[A, B] - C.c(1, 4) * dewitt(g0, E[A], E[B], C.zero, C)
    return HH, oarr((n, N), C.zero), VV


def frame_data(Y, C):
    """Explicit Lemma K2 data in the adapted frame at p (h_p, nabla delta, nabla^2 delta), symbolic or numeric."""
    n, N, m, E, g0 = Y.n, Y.N, Y.m, Y.E, Y.base.g0
    gi = mat_inv(g0, C)
    D1 = oarr((m, m, m), C.zero)
    D2 = oarr((m, m, m, m), C.zero)
    for A in range(N):
        D1[n + A, :n, :n] = E[A]
        for B in range(N):
            D2[n + A, n + B, :n, :n] = (ein('ij,jk,kl->il', E[A], gi, E[B]) + ein('ij,jk,kl->il', E[B], gi, E[A])) * C.c(1, 2)
    return D1, D2


if __name__ == '__main__':
    P("K_z1_ricci.py  (exact: Fractions and sympy rational-function fields; seed random.Random(20261030))")
    Q = Num()

    # ------------------------------------------------------------------ K0
    section("K0  controls of the curvature routine: 4|dx|^2/(1 +- |x|^2)^2 has Ric = +-(m-1) g at 0")
    for m in (2, 3, 4):
        for sgn in (1, -1):
            g0 = oarr((m, m), Q.zero)
            dd = oarr((m, m, m, m), Q.zero)
            for i in range(m):
                g0[i, i] = Q.c(4)
                for k in range(m):
                    dd[k, k, i, i] = Q.c(-16 * sgn)
            Ric, _, _ = ricci_2jet(g0, oarr((m, m, m), Q.zero), dd, Q)
            check(f"m={m}, {'sphere' if sgn > 0 else 'hyperbolic space'}: Ric == {sgn}*(m-1)*g", eq(Ric, sgn * (m - 1) * g0))

    # ------------------------------------------------------------------ K1
    section("K1  variation formula (Lemma K1): at a point where delta = 0, with nabla = nabla^h,\n"
            "    Ric(h+delta) - Ric(h) = (1/2) h^ad (nabla_a nabla_b delta_dc + nabla_a nabla_c delta_bd - nabla_a nabla_d delta_bc\n"
            "      - nabla_b nabla_c delta_ad) - (1/2) h^ap h^dq nabla_a delta_pq T_dbc + (1/2) h^ap h^dq nabla_b delta_pq nabla_c delta_ad\n"
            "      + (1/4) h^ap h^dq T_pad T_qbc - (1/4) h^ap h^dq T_pbd T_qac,   T_dbc = nabla_b delta_dc + nabla_c delta_db - nabla_d delta_bc.\n"
            "    Test: random exact 2-jets of h and delta (delta(0) = 0), both Riccis by the independent routine.")
    for m, eps in [(2, [1, 1]), (2, [1, -1]), (3, [1, 1, 1]), (3, [1, -1, -1]), (4, [1, 1, -1, -1]), (5, [1, -1, 1, 1, -1])]:
        h0 = rand_g0(eps, True, Q)
        h1, h2 = rand_jets(m, Q)
        e1, e2 = rand_jets(m, Q)
        Rh, hinv, Gam = ricci_2jet(h0, h1, h2, Q)
        Rhd, _, _ = ricci_2jet(h0, h1 + e1, h2 + e2, Q)
        D1, D2 = nabla_delta(e1, e2, Gam)
        pred = var_numeric(hinv, D1, D2, Q)
        lin = var_numeric(hinv, D1, D2, Q, terms=VAR_TERMS[:4])
        check(f"m={m} sig={(eps.count(1), eps.count(-1))}: Ric(h+delta)-Ric(h) == formula "
              f"(nabla^2 delta symmetric in its derivative slots: {eq(D2, np.einsum('abcd->bacd', D2))}; "
              f"quadratic part != 0: {not eq(pred, lin)})", eq(Rhd - Rh, pred))

    # ------------------------------------------------------------------ K2
    section("K2  Lemma K2 on curved bases (random 3-jets, non-normal coordinates), Z2 at p = (0, g(0)) in s(X):\n"
            "    (a) s_* = horizontal lift; (b) nabla pi_* = 0 at p, i.e. Gamma^{Z2,rho}_{ab}(p) = Gamma^{g,rho}_{ab}(0) for\n"
            "    base directions a, b and 0 if a or b is vertical (A = 0 on s(X), T = 0); (c) nabla K at p has only the VV part\n"
            "    (1/2)(U g^-1 W + W g^-1 U); (d) nabla delta, nabla^2 delta at p (adapted frame) are exactly those of K2;\n"
            "    (e) fibre Christoffel symbols Gamma^F(U,W) = -(1/2)(U y^-1 W + W y^-1 U), independent of lambda.")
    KS, lamS, cHS, cVS = field('lam,cH,cV', QQ)
    S = Num(KS)
    for n, eps in [(2, [1, -1]), (3, [1, 1, -1]), (4, [1, -1, -1, -1])]:
        t = time.time()
        base = random_base(n, eps, True, S)
        Y = YJets(base, order=1)
        N, m, E = Y.N, Y.m, Y.E
        G0, G1, _ = Y.metric('Z2', lamS, cHS, cVS)
        gi, _, GamY = christoffel(G0, G1, S)
        ok_a = all(Y.Psi0[A, mu] == base.g1[mu, i, j] for A, (i, j) in enumerate(Y.idx) for mu in range(n))
        ok_b = eq(GamY[:n, :n, :n], base.Gam) and is_zero(GamY[:n, n:, :]) and is_zero(GamY[:n, :, n:])
        # nabla K at p:  (nabla_a K)_b = nabla^E_a(K_b) - Gamma^c_{ab} K_c, K_mu = -(Gamma_mu^T y + y Gamma_mu), K_B = E_B
        y0 = Y.y0

        def psimat(Gam_, y):
            return ein('kmi,kj->mij', Gam_, y) + ein('kmj,ik->mij', Gam_, y)
        K0 = oarr((m, n, n), S.zero)
        K0[:n] = -psimat(base.Gam, y0)
        K0[n:] = E
        dK = oarr((m, m, n, n), S.zero)                      # dK[a, b] = d_a K_b
        for r in range(n):
            dK[r, :n] = -psimat(base.dGam[r], y0)
        for B in range(N):
            dK[n + B, :n] = -psimat(base.Gam, E[B])
        nK = dK - ein('cab,cij->abij', GamY, K0)
        for r in range(n):                                   # connection term of nabla^E along base directions
            Gr = base.Gam[:, r, :]
            nK[r] = nK[r] - ein('ki,bkj->bij', Gr, K0) - ein('bik,kj->bij', K0, Gr)
        F = Y.frame()
        nKf = ein('Aa,Bb,abij->ABij', F, F, nK)
        expK = oarr((m, m, n, n), S.zero)
        y0i = mat_inv(y0, S)
        for A in range(N):
            for B in range(N):
                expK[n + A, n + B] = (ein('ij,jk,kl->il', E[A], y0i, E[B]) + ein('ij,jk,kl->il', E[B], y0i, E[A])) * S.c(1, 2)
        ok_c = eq(nKf, expK)
        # delta = cH (y - g(x))(dx, dx): coordinate jets at p
        d1 = oarr((m, m, m), S.zero)
        d2 = oarr((m, m, m, m), S.zero)
        for r in range(n):
            d1[r, :n, :n] = base.g1[r] * (-cHS)
            for s_ in range(n):
                d2[r, s_, :n, :n] = base.g2[r, s_] * (-cHS)
        for B in range(N):
            d1[n + B, :n, :n] = E[B] * cHS
        D1, D2 = nabla_delta(d1, d2, GamY)
        ok_sparse = eq(D2, nabla_delta(d1, d2, GamY, nh=n)[1])          # the sparse evaluation used in K4
        D1f, D2f = to_frame(D1, F), to_frame(D2, F)
        X1, X2 = frame_data(Y, S)
        ok_d = eq(D1f, X1 * cHS) and eq(D2f, X2 * cHS) and ok_sparse
        # fibre Christoffels of cV G_lam at y0, as matrices sum_C Gamma^C_AB E_C
        Gv0 = Y.Gv['G0'][0] - Y.Gv['TT'][0] * lamS
        Gv1 = Y.Gv['G0'][1][n:] - Y.Gv['TT'][1][n:] * lamS
        _, _, GamF = christoffel(Gv0 * cVS, Gv1 * cVS, S)
        ok_e = all(eq(ein('C,Cij->ij', GamF[:, A, B], E),
                      -(ein('ij,jk,kl->il', E[A], y0i, E[B]) + ein('ij,jk,kl->il', E[B], y0i, E[A])) * S.c(1, 2))
                   for A in range(N) for B in range(N))
        P(f"  n={n} sig={(eps.count(1), eps.count(-1))}, g0 non-diagonal, Gamma_g(0) != 0: {not is_zero(base.Gam)}, "
          f"Ric_g(0) != 0: {not is_zero(base.Ric)}; coefficients in Q(lam, cH, cV)   [{time.time() - t:.1f}s]")
        check(f"n={n}: (a) s_* = lift {ok_a}; (b) nabla pi_* = 0 at p {ok_b}; (c) nabla K at p {ok_c}; "
              f"(d) nabla delta, nabla^2 delta at p {ok_d}; (e) fibre Christoffels {ok_e}", ok_a and ok_b and ok_c and ok_d and ok_e)

    # ------------------------------------------------------------------ K3
    section("K3  Lemma K3 (general n): the variation formula of K1 on the data of K2, by O(g)-invariant strand diagrams\n"
            "    with symbolic loop weight n; coefficients in Q(n, lam, cH, cV), mu = lam/(1 - n lam).")
    KB, nB, lamB, cHB, cVB = field('n,lam,cH,cV', QQ)
    B_ = Num(KB)
    data = dg_data(B_, cHB)
    hh = dg_orbits(dg_evaluate(VAR_TERMS, data, {'b': 'H', 'c': 'H'}, B_, nB, lamB, cHB, cVB), {'b': 'H', 'c': 'H'})
    vv = dg_orbits(dg_evaluate(VAR_TERMS, data, {'b': 'V', 'c': 'V'}, B_, nB, lamB, cHB, cVB), {'b': 'V', 'c': 'V'})
    hv = dg_evaluate(VAR_TERMS, data, {'b': 'H', 'c': 'V'}, B_, nB, lamB, cHB, cVB)
    P(f"  dRic_HH = ({hh.get(HH_REP, 0)}) g")
    P(f"  dRic_VV = ({vv.get(G0_REP, 0)}) G_0 + ({vv.get(TT_REP, 0)}) tr(g^-1 .) tr(g^-1 .)")
    P(f"  dRic_HV: {len(hv)} diagrams (an H-V bilinear form has 3 endpoints, so no perfect matching exists)")
    check("dRic_HH == cH/(4 cV (n lam - 1)) g for every n", set(hh) <= {HH_REP} and hh.get(HH_REP) == cHB / (4 * cVB * (nB * lamB - 1)))
    check("dRic_VV == -(1/4) G_0 for every n", set(vv) <= {G0_REP} and vv.get(G0_REP) == B_.c(-1, 4))
    check("dRic_HV == 0 for every n", hv == {})
    # term-by-term table (for the write-up)
    P("  term-by-term (HH coefficient of g; VV coefficient of G_0, of tr tr):")
    for k, term in enumerate(VAR_TERMS, 1):
        a = dg_orbits(dg_evaluate([term], data, {'b': 'H', 'c': 'H'}, B_, nB, lamB, cHB, cVB), {'b': 'H', 'c': 'H'})
        b = dg_orbits(dg_evaluate([term], data, {'b': 'V', 'c': 'V'}, B_, nB, lamB, cHB, cVB), {'b': 'V', 'c': 'V'})
        P(f"    term {k}: HH {a.get(HH_REP, 0)};  VV {b.get(G0_REP, 0)}, {b.get(TT_REP, 0)}")
    # F4(a): II = vertical part of the Christoffel difference C^A_{mu nu} = (1/2) h^{AB} T_{B mu nu}
    ii = dg_orbits(dg_evaluate([(Fr(1, 2), [('T', 'Bmn')], ['AB'])], data, {'A': 'V', 'm': 'H', 'n': 'H'},
                               B_, nB, lamB, cHB, cVB), {'A': 'V', 'm': 'H', 'n': 'H'})
    muB = lamB / (1 - nB * lamB)
    exp_ii = {((('A', 0), ('m', 0)), (('A', 1), ('n', 0))): -cHB / (2 * cVB),
              ((('A', 0), ('A', 1)), (('m', 0), ('n', 0))): -cHB / (2 * cVB) * muB}
    check("C^V(X,Y) = -(cH/(2 cV)) [X^flat . Y^flat + mu g(X,Y) g]  (= II of F4(a))", ii == exp_ii)
    # F4(b): trace of the Ricci blocks
    fib = dg_sym_V({('V', 'V'): [(-nB / 4, (((0, 0), (1, 0)), ((0, 1), (1, 1)))),
                                 (B_.c(1, 4), (((0, 0), (0, 1)), ((1, 0), (1, 1))))]})      # Ric_F = -(n/4) G_{1/n}
    tr_d = dg_evaluate([(c, f, p + ['bc']) for c, f, p in VAR_TERMS], data, {}, B_, nB, lamB, cHB, cVB).get((), 0)
    tr_f = dg_evaluate([(Fr(1), [('F', 'bc')], ['bc'])], dict(F=fib), {}, B_, nB, lamB, cHB, cVB).get((), 0)
    sigma1 = -nB * (nB - 1) * (nB + 2) / 8 - nB / 4 * ((1 + lamB) / (1 - nB * lamB) + (nB + 1) / 2)
    check("fibre trace h^-1 . Ric_F = s_F = -n(n-1)(n+2)/(8 cV)", tr_f == -nB * (nB - 1) * (nB + 2) / (8 * cVB))
    check("s_Y - R/cH = h^-1 . (Ric_F + dRic) == sigma_1/cV of F4(b)", tr_f + tr_d == sigma1 / cVB)
    check("-(n/4) G_{1/n} - (1/4) G_0 == -((n+1)/4) G_{1/(n+1)} (coefficients of G_0 and tr tr)",
          (-nB / 4 - B_.c(1, 4), nB / 4 * (1 / nB)) == (-(nB + 1) / 4, (nB + 1) / 4 * (1 / (nB + 1))))
    # cross-check the engine against explicit tensors at n = 1..4 (coefficients in Q(lam, cH, cV))
    P("  engine cross-check: explicit h_p, nabla delta, nabla^2 delta in the adapted frame, n = 1..4:")
    for n, eps in [(1, [-1]), (2, [1, -1]), (3, [1, 1, 1]), (4, [1, -1, -1, -1])]:
        t = time.time()
        Y = YJets(flat_base(rand_g0(eps, True, S), S))
        N, m, E, g0 = Y.N, Y.m, Y.E, Y.base.g0
        X1, X2 = frame_data(Y, S)
        h = oarr((m, m), S.zero)
        h[:n, :n] = g0 * cHS
        h[n:, n:] = (Y.Gv['G0'][0] - Y.Gv['TT'][0] * lamS) * cVS
        dR = var_numeric(mat_inv(h, S), X1 * cHS, X2 * cHS, S)
        data_n = dg_data(S, cHS)
        dh = dg_orbits(dg_evaluate(VAR_TERMS, data_n, {'b': 'H', 'c': 'H'}, S, S.c(n), lamS, cHS, cVS),
                       {'b': 'H', 'c': 'H'}).get(HH_REP, 0)
        dv = dg_orbits(dg_evaluate(VAR_TERMS, data_n, {'b': 'V', 'c': 'V'}, S, S.c(n), lamS, cHS, cVS), {'b': 'V', 'c': 'V'})
        exp = oarr((m, m), S.zero)
        exp[:n, :n] = g0 * dh
        for A in range(N):
            for B in range(N):
                exp[n + A, n + B] = (dv.get(G0_REP, 0) * dewitt(g0, E[A], E[B], S.zero, S)
                                     + dv.get(TT_REP, 0) * (dewitt(g0, E[A], E[B], S.zero, S) - dewitt(g0, E[A], E[B], S.one, S)))
        check(f"n={n} sig={(eps.count(1), eps.count(-1))}: explicit evaluation == diagram evaluation at loop weight {n}   "
              f"[{time.time() - t:.1f}s]", eq(dR, exp))
    # F6 (Z1 half) from the blocks: Ric = kappa g_aleph on s(X) iff lam = 1/(n+1) and Ric(g) = 0
    import sympy as sp
    nS, lS, hS, vS, kS = sp.symbols('n lam cH cV kappa')
    sol = sp.solve([-(nS + 1) / 4 - kS * vS, sp.Rational(1, 4) + kS * vS * lS], [kS, lS], dict=True)
    ric_g = sp.simplify(kS * hS - hS / (4 * vS * (nS * lS - 1))).subs(sol[0]) if sol else None
    check(f"VV block proportional to cV G_lam iff (kappa, lam) = {sol}; then Ric_g = ({sp.simplify(ric_g)}) g  (F6, Z1 half)",
          len(sol) == 1 and sp.simplify(sol[0][lS] - 1 / (nS + 1)) == 0 and sp.simplify(ric_g) == 0)

    # ------------------------------------------------------------------ K4
    section("K4  independent brute force: Ricci tensor of Z1 and Z2 at p = (0, g(0)) from the exact coordinate 2-jet of the\n"
            "    metric on Y^(n(n+3)/2); random curved base 3-jets, g0 non-diagonal, non-normal coordinates. Blocks in the adapted\n"
            "    frame vs F4(c) (Z1) and F3 (Z2); 'chain': Ric(Z1) - Ric(Z2) == formula K1 evaluated on the brute-force\n"
            "    nabla^{Z2} delta, nabla^{Z2} nabla^{Z2} delta (coordinates).  n = 2, 3: coefficients in Q(lam, cH, cV);\n"
            "    n = 4, 5, 6: exact rationals at the listed (lam, cH, cV).")
    PARAMS = [(Fr(1, 2), Fr(1), Fr(1)), (Fr(2, 7), Fr(-5, 3), Fr(3, 4)), (Fr(1), Fr(3, 2), Fr(-2, 3)),
              (Fr(-3, 5), Fr(2), Fr(5, 2)), (Fr(1, 5), Fr(7, 4), Fr(-1, 3)), (Fr(3, 8), Fr(-1, 2), Fr(-4, 5))]
    runs = [(2, [1, 1], 'sym'), (2, [1, -1], 'sym'), (3, [1, 1, 1], 'sym'), (3, [1, 1, -1], 'sym'),
            (4, [1, -1, -1, -1], 2), (4, [1, 1, -1, -1], 2), (4, [1, 1, 1, 1], 2),
            (5, [1, 1, 1, 1, 1], 2), (5, [1, 1, 1, -1, -1], 2), (6, [1, 1, 1, 1, 1, -1], 2), (6, [1, 1, 1, -1, -1, -1], 2)]
    pk = 0
    summary = []
    for n, eps, mode in runs:
        t = time.time()
        Cn = S if mode == 'sym' else Q
        base = random_base(n, eps, True, Cn)
        Y = YJets(base)
        N, m, E = Y.N, Y.m, Y.E
        F = Y.frame()
        nonpar = not is_zero(base.nablaRic)
        notein = not eq(base.Ric * base.g0[0, 0], base.g0 * base.Ric[0, 0])
        P(f"  n={n} sig={(eps.count(1), eps.count(-1))} (dim Y = {m}): Gamma_g(0) != 0 {not is_zero(base.Gam)}, "
          f"Ric_g(0) not prop. to g {notein}{' (n = 2: always prop.)' if n == 2 else ''}, nabla Ric_g(0) != 0 {nonpar}, Ric routes agree {eq(base.Ric, base.Ric_b)}"
          f"   [jets {time.time() - t:.1f}s]")
        check(f"n={n}: test base is generic (non-normal coordinates, nabla Ric != 0{', Ric not Einstein' if n > 2 else ''})",
              not is_zero(base.Gam) and (notein or n == 2) and nonpar and eq(base.Ric, base.Ric_b))
        if mode == 'sym':
            plist = [(lamS, cHS, cVS)]
        else:
            plist = [tuple(Q.c(x) for x in PARAMS[(pk + j) % len(PARAMS)]) for j in range(mode)]
            pk += mode
        # coordinate jets of delta = cH (y - g(x))(dx, dx) (cH factored out)
        d1 = oarr((m, m, m), Cn.zero)
        d2 = oarr((m, m, m, m), Cn.zero)
        for r in range(n):
            d1[r, :n, :n] = -base.g1[r]
            for s_ in range(n):
                d2[r, s_, :n, :n] = -base.g2[r, s_]
        for B in range(N):
            d1[n + B, :n, :n] = E[B]
        for lam, cH, cV in plist:
            t = time.time()
            res = {}
            for variant in ('Z2', 'Z1'):
                G0, G1, G2 = Y.metric(variant, lam, cH, cV)
                Ric, gi, Gam = ricci_2jet(G0, G1, G2, Cn)
                res[variant] = (Ric, gi, Gam)
                HH, HV, VV = ricci_frame_blocks(Ric, F, n)
                pHH, pHV, pVV = predicted_blocks(Y, variant, lam, cH, cV)
                ok = (eq(HH, pHH), eq(HV, pHV), eq(VV, pVV))
                lab = ('Ric_HH - Ric_g == cH/(4cV(n lam-1)) g' if variant == 'Z1' else 'Ric_HH == Ric_g') + f' {ok[0]}; ' \
                      + f'Ric_HV == 0 {ok[1]}; ' + ('Ric_VV == -((n+1)/4) G_1/(n+1)' if variant == 'Z1' else 'Ric_VV == -(n/4) G_1/n') \
                      + f' {ok[2]}'
                summary.append((n, tuple(eps), variant, all(ok)))
                check(f"n={n} sig={(eps.count(1), eps.count(-1))} {variant} (lam, cH, cV) = "
                      f"{'symbolic' if mode == 'sym' else (str(lam), str(cH), str(cV))}: {lab}", all(ok))
            RZ2, giZ2, GamZ2 = res['Z2']
            D1, D2 = nabla_delta(d1 * cH, d2 * cH, GamZ2, nh=n)
            chain = eq(res['Z1'][0] - RZ2, var_numeric(giZ2, D1, D2, Cn))
            check(f"      chain: Ric(Z1) - Ric(Z2) == K1 formula on brute-force nabla delta, nabla^2 delta   "
                  f"[{time.time() - t:.1f}s]", chain)

    # ------------------------------------------------------------------ K5
    section("K5  the flat model of the F4(c) sketch, cH y dx^2 + cV G_lam(dy, dy), at a point y0 off the section of eta:\n"
            "    Ric_{mu nu} = cH/(4 cV (n lam - 1)) y0, Ric_{mu A} = 0, Ric_AB = -(n/4) G_1/n|y0 - (1/4) G_0|y0 (coordinates).")
    for n, eps, Cn in [(3, [1, 1, -1], S), (4, [1, -1, -1, -1], Q)]:
        t = time.time()
        eta = rand_g0(eps, False, Cn)
        y0 = rand_g0([e * (k + 2) for k, e in enumerate(eps)], True, Cn)
        Y = YJets(flat_base(eta, Cn), y0=y0)
        lam, cH, cV = (lamS, cHS, cVS) if Cn is S else (Q.c(2, 9), Q.c(-3, 4), Q.c(5, 3))
        G0, G1, G2 = Y.metric('Z1', lam, cH, cV)
        Ric, _, _ = ricci_2jet(G0, G1, G2, Cn)
        N, E = Y.N, Y.E
        nn = Cn.c(n)
        pVV = oarr((N, N))
        for A in range(N):
            for B in range(N):
                pVV[A, B] = -nn / 4 * dewitt(y0, E[A], E[B], 1 / nn, Cn) - Cn.c(1, 4) * dewitt(y0, E[A], E[B], Cn.zero, Cn)
        ok = (eq(Ric[:n, :n], y0 * (cH / (4 * cV * (nn * lam - 1)))), is_zero(Ric[:n, n:]), eq(Ric[n:, n:], pVV))
        check(f"n={n} sig(y0)={(eps.count(1), eps.count(-1))}, y0 != eta: HH {ok[0]}, HV {ok[1]}, VV {ok[2]}"
              f"   [{time.time() - t:.1f}s]", all(ok))

    # ------------------------------------------------------------------ summary
    section("Summary")
    P(f"  K4 brute force: {sum(1 for s in summary if s[3])} of {len(summary)} (n, signature, variant, parameters) cases pass; "
      f"n covered: {sorted(set(s[0] for s in summary))}")
    P(f"  all {len(CHECKS)} checks pass: {all(ok for _, ok in CHECKS)}")
    P(f"\nTotal time {time.time() - T0:.1f}s")
