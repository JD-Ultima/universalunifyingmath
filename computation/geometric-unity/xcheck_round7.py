"""Referee cross-check of round 7 (REVIEW5-GU): K-z1-ricci-general-n.md (Theorem F4(c) for every n) and
J5d-indefinite-leaf.md (the indefinite case of Theorem J5(d)).

Imports NO repository code (not K_z1_ricci.py, not J5d_leaf.py, not F_*.py).  Everything is rebuilt here:
  * truncated multivariate Taylor polynomials ("jets") with exact Fractions, used to build the 2-jets of the Zorro
    metrics Z1, Z2 on Y = Met(X) at a point of s(X) from a random polynomial base metric;
  * a textbook Ricci routine (Christoffel symbols and their first derivatives from the metric 2-jet);
  * the referee's own form of Lemma K1 (re-derived by hand in REVIEW5-GU) and of the K2 data;
  * algebraic curvature tensors by Bianchi projection, Weyl tensors by Kulkarni-Nomizu, the forms F and Phi_y of
    J5d from their definitions (action B.y = -(B^T y + y B), G_0(h,k) = tr(g^-1 h g^-1 k)).
Sections
  X0  controls: unit-sphere metric; Sherman-Morrison facts (i)-(iii) of Lemma K3, n = 1..6.
  X1  Lemma K1 on random exact 2-jets (dimensions 4 (2,2) and 5 (3,2)).
  X2  K3 term table: the eight terms evaluated on the K2 data in the adapted frame, n = 1..4, rational parameters.
  X3  brute force: Ricci of Z1 and Z2 at p in s(X) on random curved polynomial bases, n = 3 ((2,1), (3,0)) and
      n = 4 ((1,3), (2,2)); K2 jets of delta from the brute-force Z2 connection; chain Ric(Z1)-Ric(Z2) = K1 formula.
  X4  J5d L0: Phi_{g+eps h} = eps^2 F(X,Y;h,h) + O(eps^3) (eps^0, eps^1 vanish; h traceless or not); lambda-term.
  X5  J5d L1: C1, C2 as defined in the write-up; their identities; the h <-> h' asymmetry of the polarised form
      (C2 != C2'); the symmetrised trace C2sym (the one the leaf condition actually constrains) and its closed form;
      Einstein case R = W + (Lambda/3) R_1; n = 4 in (4,0), (1,3), (2,2), plus n = 3 (1,2) and n = 5 (2,3).
  X6  J5d L3/L4: positivity identity at a random T; an explicit (2,2) solution R = n (x) n (own construction) checked
      against F = 0 (quadratic and polarised) and the full condition Phi_y = 0 at random y.
Exact arithmetic; fixed seed random.Random(20261107); one core.
Run from computation/geometric-unity:  python3 xcheck_round7.py > xcheck_round7.out 2>&1
"""
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

T0 = time.time()
RNG = random.Random(20261107)
RES = []


def say(*a):
    print(*a)
    sys.stdout.flush()


def check(label, ok):
    ok = bool(ok)
    RES.append((label, ok))
    say(f"  [{ok}] {label}")


def rq(lo=-3, hi=3, dmax=3):
    return Fr(RNG.randint(lo, hi), RNG.randint(1, dmax))


def oz(shape):
    a = np.empty(shape, dtype=object)
    a.fill(Fr(0))
    return a


def ein(spec, *ops):
    return np.einsum(spec, *ops, optimize=True)


def zero(a):
    return all(x == 0 for x in np.asarray(a, dtype=object).flat)


def inv(M):
    m = M.shape[0]
    A = [[Fr(M[i, j]) for j in range(m)] + [Fr(int(i == j)) for j in range(m)] for i in range(m)]
    for c in range(m):
        p = next(r for r in range(c, m) if A[r][c] != 0)
        A[c], A[p] = A[p], A[c]
        pv = A[c][c]
        A[c] = [x / pv for x in A[c]]
        for r in range(m):
            if r != c and A[r][c] != 0:
                f = A[r][c]
                A[r] = [x - f * y for x, y in zip(A[r], A[c])]
    out = oz((m, m))
    for i in range(m):
        for j in range(m):
            out[i, j] = A[i][m + j]
    return out


# ============================================================================ jets (truncated Taylor polynomials)
# a jet is a dict {sorted tuple of variable indices: Fraction}; the key () is the constant term
def jmul(a, b, D):
    out = {}
    for ka, va in a.items():
        la = len(ka)
        for kb, vb in b.items():
            if la + len(kb) > D:
                continue
            k = tuple(sorted(ka + kb))
            out[k] = out.get(k, 0) + va * vb
    return {k: v for k, v in out.items() if v != 0}


def jadd(*js):
    out = {}
    for j in js:
        for k, v in j.items():
            out[k] = out.get(k, 0) + v
    return {k: v for k, v in out.items() if v != 0}


def jsc(a, c):
    return {k: v * c for k, v in a.items()} if c != 0 else {}


def jtrunc(a, D):
    return {k: v for k, v in a.items() if len(k) <= D}


def jder(a, c):
    out = {}
    for k, v in a.items():
        m = k.count(c)
        if m:
            kk = list(k)
            kk.remove(c)
            kk = tuple(kk)
            out[kk] = out.get(kk, 0) + m * v
    return out


def mmul(A, B, D):
    n, m, p = len(A), len(B), len(B[0])
    return [[jadd(*[jmul(A[i][k], B[k][j], D) for k in range(m)]) for j in range(p)] for i in range(n)]


def minv_jet(M, D):
    """inverse of a matrix of jets, truncated at degree D: sum_k (-M0^-1 N)^k M0^-1"""
    n = len(M)
    M0 = oz((n, n))
    for i in range(n):
        for j in range(n):
            M0[i, j] = M[i][j].get((), Fr(0))
    I0 = inv(M0)
    I0j = [[({(): I0[i, j]} if I0[i, j] != 0 else {}) for j in range(n)] for i in range(n)]
    N = [[{k: v for k, v in M[i][j].items() if k != ()} for j in range(n)] for i in range(n)]
    K = mmul(I0j, N, D)
    K = [[jsc(x, -1) for x in row] for row in K]
    out, term = I0j, I0j
    for _ in range(D):
        term = mmul(K, term, D)
        out = [[jadd(out[i][j], term[i][j]) for j in range(n)] for i in range(n)]
    return out


def jets_to_arrays(G, m):
    """G: m x m matrix of jets in variables 0..m-1 -> (G0, G1[c,a,b] = d_c G_ab, G2[c,d,a,b])"""
    G0, G1, G2 = oz((m, m)), oz((m, m, m)), oz((m, m, m, m))
    for a in range(m):
        for b in range(m):
            jt = G[a][b]
            G0[a, b] = jt.get((), Fr(0))
            for c in range(m):
                G1[c, a, b] = jt.get((c,), Fr(0))
                for d in range(m):
                    k = (min(c, d), max(c, d))
                    G2[c, d, a, b] = jt.get(k, Fr(0)) * (2 if c == d else 1)
    return G0, G1, G2


# ============================================================================ Ricci from a 2-jet (textbook route)
def christ(G0, G1, G2):
    gi = inv(G0)
    L = (np.einsum('bdc->dbc', G1) + np.einsum('cdb->dbc', G1) - G1) * Fr(1, 2)        # L[d,b,c] = Gamma_{d,bc}
    Gam = ein('ad,dbc->abc', gi, L)
    dgi = -ein('ap,epq,qd->ead', gi, G1, gi)
    dL = (np.einsum('ebdc->edbc', G2) + np.einsum('ecdb->edbc', G2) - G2) * Fr(1, 2)  # dL[e,d,b,c] = d_e Gamma_{d,bc}
    dGam = ein('ead,dbc->eabc', dgi, L) + ein('ad,edbc->eabc', gi, dL)
    return gi, Gam, dGam


def ricci(G0, G1, G2):
    """Ric_bd = d_a Gamma^a_bd - d_d Gamma^a_ab + Gamma^a_ae Gamma^e_bd - Gamma^a_de Gamma^e_ab (unit sphere: (m-1) g)"""
    gi, Gam, dGam = christ(G0, G1, G2)
    Ric = (np.einsum('aabd->bd', dGam) - np.einsum('daab->bd', dGam) + ein('aae,ebd->bd', Gam, Gam)
           - ein('ade,eab->bd', Gam, Gam))
    return Ric, gi, Gam


# ============================================================================ the referee's Lemma K1
def k1_terms(hinv, D1, D2):
    """the eight terms of Lemma K1 (derived by hand in REVIEW5-GU), D1[a,b,c] = nabla_a delta_bc,
    D2[a,b,c,d] = nabla_a nabla_b delta_cd at a zero of delta; T[d,b,c] = D1[b,d,c] + D1[c,d,b] - D1[d,b,c]"""
    T = np.einsum('bdc->dbc', D1) + np.einsum('cdb->dbc', D1) - D1
    h = Fr(1, 2)
    q = Fr(1, 4)
    return [ein('ad,abdc->bc', hinv, D2) * h,
            ein('ad,acbd->bc', hinv, D2) * h,
            -ein('ad,adbc->bc', hinv, D2) * h,
            -ein('ad,bcad->bc', hinv, D2) * h,
            -ein('ap,dq,apq,dbc->bc', hinv, hinv, D1, T) * h,
            ein('ap,dq,bpq,cad->bc', hinv, hinv, D1, D1) * h,
            ein('ap,dq,pad,qbc->bc', hinv, hinv, T, T) * q,
            -ein('ap,dq,pbd,qac->bc', hinv, hinv, T, T) * q]


def nabla_delta(d1, d2, Gam):
    """nabla delta, nabla^2 delta at a point where delta = 0 (d1 = d delta, d2 = dd delta there)"""
    D1 = d1
    D2 = (d2 - ein('ebc,aed->abcd', Gam, d1) - ein('ebd,ace->abcd', Gam, d1)
          - ein('eab,ecd->abcd', Gam, D1) - ein('eac,bed->abcd', Gam, D1) - ein('ead,bce->abcd', Gam, D1))
    return D1, D2


def rand_g0(eta):
    m = len(eta)
    U = oz((m, m))
    for i in range(m):
        U[i, i] = Fr(1)
        for j in range(i + 1, m):
            U[i, j] = Fr(RNG.randint(-1, 1))
    Dm = oz((m, m))
    for i in range(m):
        Dm[i, i] = Fr(eta[i])
    return ein('ki,kl,lj->ij', U, Dm, U)


def rsym(m):
    a = oz((m, m))
    for i in range(m):
        for j in range(i, m):
            a[i, j] = a[j, i] = rq()
    return a


# ============================================================================ X0
def X0():
    say("\nX0  controls")
    m = 3
    G0, G1, G2 = oz((m, m)), oz((m, m, m)), oz((m, m, m, m))
    for i in range(m):                    # 4|dx|^2/(1+|x|^2)^2 = 4 - 8|x|^2 + ...
        G0[i, i] = Fr(4)
        for c in range(m):
            G2[c, c, i, i] = Fr(-16)
    Ric, _, _ = ricci(G0, G1, G2)
    check("unit sphere 4|dx|^2/(1+|x|^2)^2 at 0, m = 3: Ric == (m-1) g", zero(Ric - 2 * G0))
    ok = True
    for n in range(1, 7):
        lam = Fr(2, 7) if n != 3 else Fr(-3, 5)
        mu = lam / (1 - n * lam)
        g = rand_g0([1] * (n // 2) + [-1] * (n - n // 2)) if n > 1 else oz((1, 1)) + Fr(-1)
        gi = inv(g)
        idx = [(i, j) for i in range(n) for j in range(i, n)]
        E = []
        for (i, j) in idx:
            e = oz((n, n))
            e[i, j] = e[j, i] = Fr(1)
            E.append(e)
        N = len(E)
        Gl = oz((N, N))
        for A in range(N):
            for B in range(N):
                a, b = ein('ij,jk->ik', gi, E[A]), ein('ij,jk->ik', gi, E[B])
                Gl[A, B] = ein('ij,ji->', a, b) - lam * np.trace(a) * np.trace(b)
        Gi = inv(Gl)
        S = ein('AB,Aij,Bkl->ijkl', Gi, np.array(E), np.array(E))
        rhs = (ein('ik,jl->ijkl', g, g) + ein('il,jk->ijkl', g, g)) * Fr(1, 2) + ein('ij,kl->ijkl', g, g) * mu
        ok &= zero(S - rhs)
        ok &= zero(ein('AB,Aij,jk,Bkl->il', Gi, np.array(E), gi, np.array(E)) - g * (Fr(n + 1, 2) + mu))
        ok &= zero(ein('AB,Aij,ij,Bkl->kl', Gi, np.array(E), gi, np.array(E)) - g / (1 - n * lam))
    check("Lemma K3 (i)-(iii): G_lam^-1 = (1/2)(g g + g g) + mu g g, sum G^AB E_A g^-1 E_B = ((n+1)/2 + mu) g,"
          " sum G^AB tr(g^-1 E_A) E_B = g/(1 - n lam), exact, n = 1..6, non-diagonal indefinite g", ok)


# ============================================================================ X1
def X1():
    say("\nX1  Lemma K1 on random exact 2-jets (h, delta with delta(0) = 0)")
    for eta in ((1, 1, -1, -1), (1, 1, 1, -1, -1)):
        m = len(eta)
        h0 = rand_g0(eta)
        h1, h2, d1, d2 = oz((m, m, m)), oz((m, m, m, m)), oz((m, m, m)), oz((m, m, m, m))
        for c in range(m):
            h1[c], d1[c] = rsym(m), rsym(m)
            for d in range(c, m):
                s, t = rsym(m), rsym(m)
                h2[c, d] = h2[d, c] = s
                d2[c, d] = d2[d, c] = t
        Rh, hinv, Gam = ricci(h0, h1, h2)
        Rhd, _, _ = ricci(h0, h1 + d1, h2 + d2)
        D1, D2 = nabla_delta(d1, d2, Gam)
        terms = k1_terms(hinv, D1, D2)
        lhs = Rhd - Rh
        quad = sum(terms[4:], oz((m, m)))
        sym = zero(D2 - np.einsum('bacd->abcd', D2))
        check(f"m = {m}, signature ({eta.count(1)},{eta.count(-1)}): Ric(h+delta) - Ric(h) == sum of the eight terms"
              f" (quadratic part nonzero: {not zero(quad)}; nabla^2 delta symmetric in its derivative slots: {sym})",
              zero(lhs - sum(terms, oz((m, m)))) and not zero(quad) and sym)


# ============================================================================ K2 data in the adapted frame (abstract)
def k2_data(g, lam, cH, cV):
    n = g.shape[0]
    gi = inv(g)
    idx = [(i, j) for i in range(n) for j in range(i, n)]
    E = []
    for (i, j) in idx:
        e = oz((n, n))
        e[i, j] = e[j, i] = Fr(1)
        E.append(e)
    N = len(E)
    m = n + N
    h = oz((m, m))
    h[:n, :n] = g * cH
    for A in range(N):
        for B in range(N):
            a, b = ein('ij,jk->ik', gi, E[A]), ein('ij,jk->ik', gi, E[B])
            h[n + A, n + B] = (ein('ij,ji->', a, b) - lam * np.trace(a) * np.trace(b)) * cV
    D1, D2 = oz((m, m, m)), oz((m, m, m, m))
    for A in range(N):
        D1[n + A, :n, :n] = E[A] * cH
        for B in range(N):
            D2[n + A, n + B, :n, :n] = (ein('ij,jk,kl->il', E[A], gi, E[B]) + ein('ij,jk,kl->il', E[B], gi, E[A])) * (cH / 2)
    return h, D1, D2, E, gi


def G0form(E, gi):
    N = len(E)
    out, tt = oz((N, N)), oz((N, N))
    for A in range(N):
        for B in range(N):
            a, b = ein('ij,jk->ik', gi, E[A]), ein('ij,jk->ik', gi, E[B])
            out[A, B] = ein('ij,ji->', a, b)
            tt[A, B] = np.trace(a) * np.trace(b)
    return out, tt


# ============================================================================ X2
def X2():
    say("\nX2  K3 term table: the eight terms on the K2 data (adapted frame), every block")
    for n, eta, (lam, cH, cV) in ((1, (-1,), (Fr(1, 3), Fr(2), Fr(-3, 4))), (2, (1, -1), (Fr(1, 3), Fr(-5, 3), Fr(3, 4))),
                                  (3, (1, 1, -1), (Fr(1), Fr(3, 2), Fr(-2, 3))), (4, (1, -1, -1, -1), (Fr(-3, 5), Fr(2), Fr(5, 2)))):
        g = rand_g0(eta)
        h, D1, D2, E, gi = k2_data(g, lam, cH, cV)
        m = h.shape[0]
        terms = k1_terms(inv(h), D1, D2)
        mu = lam / (1 - n * lam)
        X = (cH / cV) * (Fr(n + 1, 2) + mu)
        G0f, _ = G0form(E, gi)
        N = len(E)
        expHH = [0, 0, -X / 2, 0, 0, 0, -cH / (4 * cV * (1 - n * lam)), X / 2]
        expVV = [0, 0, 0, Fr(-1, 2), 0, Fr(1, 2), 0, Fr(-1, 4)]
        ok = True
        for k in range(8):
            ok &= zero(terms[k][:n, :n] - g * expHH[k])
            ok &= zero(terms[k][n:, n:] - G0f * expVV[k])
            ok &= zero(terms[k][:n, n:]) and zero(terms[k][n:, :n])
        tot = sum(terms, oz((m, m)))
        ok &= zero(tot[:n, :n] - g * (cH / (4 * cV * (n * lam - 1)))) and zero(tot[n:, n:] + G0f / 4)
        check(f"n = {n}, sig ({eta.count(1)},{eta.count(-1)}), (lam, cH, cV) = ({lam}, {cH}, {cV}): every entry of the"
              f" K3 table (HH, VV), every HV block 0; totals cH/(4cV(n lam-1)) g and -(1/4) G_0", ok)


# ============================================================================ X3 brute force on Y = Met(X)
def X3():
    say("\nX3  brute force: Ricci of Z1, Z2 at p = (0, g(0)) from the exact 2-jet on Y; random polynomial bases")
    cases = [(3, (1, 1, -1), [(Fr(1, 2), Fr(1), Fr(1)), (Fr(-2, 7), Fr(-5, 3), Fr(3, 4))]),
             (3, (1, 1, 1), [(Fr(1), Fr(3, 2), Fr(-2, 3))]),
             (4, (1, -1, -1, -1), [(Fr(1, 2), Fr(1), Fr(1)), (Fr(3, 8), Fr(-1, 2), Fr(-4, 5))]),
             (4, (1, 1, -1, -1), [(Fr(-1, 3), Fr(7, 4), Fr(-1, 3))])]
    for n, eta, params in cases:
        t0 = time.time()
        idx = [(i, j) for i in range(n) for j in range(i, n)]
        N = len(idx)
        m = n + N
        Aof = {}
        for A, (i, j) in enumerate(idx):
            Aof[(i, j)] = Aof[(j, i)] = A
        # random polynomial base metric of degree 3 in x_0..x_{n-1}
        g0 = rand_g0(eta)
        gj = [[({(): g0[i, j]} if g0[i, j] != 0 else {}) for j in range(n)] for i in range(n)]
        for deg in (1, 2, 3):
            for key in itertools.combinations_with_replacement(range(n), deg):
                c = rsym(n)
                for i in range(n):
                    for j in range(n):
                        if c[i, j] != 0:
                            gj[i][j][key] = gj[i][j].get(key, 0) + c[i, j]
        g2j = [[jtrunc(x, 2) for x in row] for row in gj]
        gij = minv_jet(g2j, 2)
        dg = [[[jder(gj[i][j], r) for j in range(n)] for i in range(n)] for r in range(n)]   # degree 2
        Gam = [[[None] * n for _ in range(n)] for _ in range(n)]                          # Gam[k][mu][i]
        for k in range(n):
            for mu_ in range(n):
                for i in range(n):
                    acc = {}
                    for l_ in range(n):
                        s = jadd(dg[mu_][l_][i], dg[i][l_][mu_], jsc(dg[l_][mu_][i], -1))
                        acc = jadd(acc, jmul(gij[k][l_], s, 2))
                    Gam[k][mu_][i] = jsc(acc, Fr(1, 2))
        bG0, bG1, bG2 = jets_to_arrays(g2j, n)
        Ricg, _, _ = ricci(bG0, bG1, bG2)
        sg = ein('ab,ab->', inv(bG0), Ricg)
        einst = not zero(Ricg - bG0 * (sg / n))
        # y = g0 + Y, Y_ij = z_{n + A(i,j)}
        yj = [[jadd({(): g0[i, j]} if g0[i, j] != 0 else {}, {(n + Aof[(i, j)],): Fr(1)}) for j in range(n)] for i in range(n)]
        yij = minv_jet(yj, 2)
        Psi = [[None] * n for _ in range(N)]                                              # Psi[A][mu]
        for mu_ in range(n):
            for A, (i, j) in enumerate(idx):
                acc = {}
                for k in range(n):
                    acc = jadd(acc, jmul(Gam[k][mu_][i], yj[k][j], 2), jmul(Gam[k][mu_][j], yj[i][k], 2))
                Psi[A][mu_] = acc
        Emat = []
        for (i, j) in idx:
            e = [[{} for _ in range(n)] for _ in range(n)]
            e[i][j] = {(): Fr(1)}
            e[j][i] = {(): Fr(1)}
            Emat.append(e)
        PA = [mmul(yij, Emat[A], 2) for A in range(N)]
        trA = [jadd(*[PA[A][i][i] for i in range(n)]) for A in range(N)]
        G0p = [[None] * N for _ in range(N)]
        TTp = [[None] * N for _ in range(N)]
        for A in range(N):
            for B in range(A, N):
                G0p[A][B] = G0p[B][A] = jadd(*[jmul(PA[A][i][j], PA[B][j][i], 2) for i in range(n) for j in range(n)])
                TTp[A][B] = TTp[B][A] = jmul(trA[A], trA[B], 2)
        Psi0 = oz((N, n))
        for A in range(N):
            for mu_ in range(n):
                Psi0[A, mu_] = Psi[A][mu_].get((), Fr(0))
        Fr_ = oz((m, m))
        for i in range(m):
            Fr_[i, i] = Fr(1)
        Fr_[:n, n:] = Psi0.T
        E_arr = []
        for (i, j) in idx:
            e = oz((n, n))
            e[i, j] = e[j, i] = Fr(1)
            E_arr.append(e)
        G0f, TTf = G0form(E_arr, inv(g0))
        say(f"  n = {n}, base ({eta.count(1)},{eta.count(-1)}), dim Y = {m}: Gamma_g(0) != 0: "
            f"{any(Gam[k][a][b].get((), 0) != 0 for k in range(n) for a in range(n) for b in range(n))}, "
            f"Ric_g(0) != 0: {not zero(Ricg)}, Ric_g(0) not Einstein: {einst}   [jets {time.time() - t0:.1f}s]")
        for lam, cH, cV in params:
            t1 = time.time()
            Gv = [[jadd(G0p[A][B], jsc(TTp[A][B], -lam)) for B in range(N)] for A in range(N)]
            Q = [[jadd(*[jmul(Gv[A][B], Psi[B][nu], 2) for B in range(N)]) for nu in range(n)] for A in range(N)]
            res = {}
            for var in ('Z1', 'Z2'):
                Hb = yj if var == 'Z1' else g2j
                G = [[{} for _ in range(m)] for _ in range(m)]
                for mu_ in range(n):
                    for nu in range(n):
                        G[mu_][nu] = jadd(jsc(Hb[mu_][nu], cH), jsc(jadd(*[jmul(Psi[A][mu_], Q[A][nu], 2) for A in range(N)]), cV))
                    for A in range(N):
                        G[mu_][n + A] = G[n + A][mu_] = jsc(Q[A][mu_], -cV)
                for A in range(N):
                    for B in range(N):
                        G[n + A][n + B] = jsc(Gv[A][B], cV)
                res[var] = jets_to_arrays(G, m)
            Ric1, _, _ = ricci(*res['Z1'])
            Ric2, hinv2, Gam2 = ricci(*res['Z2'])
            R1f = ein('ia,ab,jb->ij', Fr_, Ric1, Fr_)
            R2f = ein('ia,ab,jb->ij', Fr_, Ric2, Fr_)
            ok1 = (zero(R1f[:n, :n] - Ricg - g0 * (cH / (4 * cV * (n * lam - 1)))) and zero(R1f[:n, n:])
                   and zero(R1f[n:, n:] + (G0f - TTf / (n + 1)) * Fr(n + 1, 4)))
            ok2 = zero(R2f[:n, :n] - Ricg) and zero(R2f[:n, n:]) and zero(R2f[n:, n:] + (G0f - TTf / n) * Fr(n, 4))
            check(f"n = {n} ({eta.count(1)},{eta.count(-1)}), (lam, cH, cV) = ({lam}, {cH}, {cV}): Ric(Z1) at p == F4(c)"
                  f" [HH = Ric_g + cH/(4cV(n lam-1)) g, HV = 0, VV = -((n+1)/4) G_1/(n+1)]", ok1)
            check(f"n = {n} ({eta.count(1)},{eta.count(-1)}), (lam, cH, cV) = ({lam}, {cH}, {cV}): Ric(Z2) at p == F3"
                  f" [HH = Ric_g, HV = 0, VV = -(n/4) G_1/n]", ok2)
            # K2 from the brute-force Z2 connection: delta = Z1 - Z2
            d0 = res['Z1'][0] - res['Z2'][0]
            d1 = res['Z1'][1] - res['Z2'][1]
            d2 = res['Z1'][2] - res['Z2'][2]
            D1, D2 = nabla_delta(d1, d2, Gam2)
            D1f = ein('ia,jb,kc,abc->ijk', Fr_, Fr_, Fr_, D1)
            D2f = ein('ia,jb,kc,ld,abcd->ijkl', Fr_, Fr_, Fr_, Fr_, D2)
            e1, e2 = oz((m, m, m)), oz((m, m, m, m))
            gi0 = inv(g0)
            for A in range(N):
                e1[n + A, :n, :n] = E_arr[A] * cH
                for B in range(N):
                    e2[n + A, n + B, :n, :n] = (ein('ij,jk,kl->il', E_arr[A], gi0, E_arr[B])
                                                + ein('ij,jk,kl->il', E_arr[B], gi0, E_arr[A])) * (cH / 2)
            ok3 = zero(d0) and zero(D1f - e1) and zero(D2f - e2)
            check(f"   K2: nabla^Z2 delta, nabla^Z2 nabla^Z2 delta at p (adapted frame) are exactly"
                  f" (nabla_U delta)(X,Y) = cH U(X,Y), (nabla^2_UW delta)(X,Y) = (cH/2)(U g^-1 W + W g^-1 U)(X,Y), all else 0", ok3)
            chain = sum(k1_terms(hinv2, D1, D2), oz((m, m)))
            check("   chain: Ric(Z1) - Ric(Z2) (brute force, coordinates) == referee's K1 formula on the brute-force jets"
                  f"   [{time.time() - t1:.1f}s]", zero(Ric1 - Ric2 - chain))


# ============================================================================ J5d helpers
def acurv(n, sig, scale=3):
    """random algebraic curvature tensor (dict of Fractions), convention R_abcd = g(R(e_c,e_d)e_b, e_a)"""
    pairs = [(a, b) for a in range(n) for b in range(a + 1, n)]
    P = {p: i for i, p in enumerate(pairs)}
    M = [[None] * len(pairs) for _ in pairs]
    for i in range(len(pairs)):
        for j in range(i, len(pairs)):
            M[i][j] = M[j][i] = Fr(RNG.randint(-scale, scale), RNG.randint(1, 2))

    def c(a, b, cc, d):
        if a == b or cc == d:
            return Fr(0)
        s = 1
        if a > b:
            a, b, s = b, a, -s
        if cc > d:
            cc, d, s = d, cc, -s
        return s * M[P[(a, b)]][P[(cc, d)]]
    R0 = {k: c(*k) for k in itertools.product(range(n), repeat=4)}
    R = {}
    for (a, b, cc, d), v in R0.items():
        R[a, b, cc, d] = v - (R0[a, b, cc, d] + R0[a, cc, d, b] + R0[a, d, b, cc]) / 3
    return R


def is_acurv(R, n):
    for a, b, c, d in itertools.product(range(n), repeat=4):
        v = R[a, b, c, d]
        if v != -R[b, a, c, d] or v != -R[a, b, d, c] or v != R[c, d, a, b]:
            return False
        if R[a, b, c, d] + R[a, c, d, b] + R[a, d, b, c] != 0:
            return False
    return True


def ric_of(R, n, sig):
    return {(b, d): sum(sig[a] * R[a, b, a, d] for a in range(n)) for b in range(n) for d in range(n)}


def norms(R, n, sig):
    Ric = ric_of(R, n, sig)
    R2 = sum(v * v * sig[k[0]] * sig[k[1]] * sig[k[2]] * sig[k[3]] for k, v in R.items())
    Ric2 = sum(v * v * sig[k[0]] * sig[k[1]] for k, v in Ric.items())
    return R2, Ric2


def kn(h, k, n):
    return {(a, b, c, d): h[a, c] * k[b, d] + h[b, d] * k[a, c] - h[a, d] * k[b, c] - h[b, c] * k[a, d]
            for a, b, c, d in itertools.product(range(n), repeat=4)}


def weyl_of(R, n, sig):
    g = {(a, b): Fr(sig[a]) if a == b else Fr(0) for a in range(n) for b in range(n)}
    Ric = ric_of(R, n, sig)
    s = sum(sig[a] * Ric[a, a] for a in range(n))
    Pm = {k: Ric[k] - s / (2 * (n - 1)) * g[k] for k in g}
    K = kn(Pm, g, n)
    return {k: R[k] - K[k] / (n - 2) for k in R}


def endos(R, n, sig):
    """R(e_c,e_d) as n x n Fraction matrices: (.)^p_q = eta_p R_pqcd"""
    return {(c, d): [[sig[p] * R[p, q, c, d] for q in range(n)] for p in range(n)] for c in range(n) for d in range(n)}


def mat(n, f):
    return [[f(i, j) for j in range(n)] for i in range(n)]


def mm(A, B):
    n = len(A)
    return [[sum(A[i][k] * B[k][j] for k in range(n)) for j in range(n)] for i in range(n)]


def madd(A, B, c=1):
    return [[A[i][j] + c * B[i][j] for j in range(len(A))] for i in range(len(A))]


def tr(A):
    return sum(A[i][i] for i in range(len(A)))


def tp(A):
    return [list(r) for r in zip(*A)]


def act(B, y):          # B.y = -(B^T y + y B)
    return [[-x for x in row] for row in madd(mm(tp(B), y), mm(y, B))]


def minv(A):
    n = len(A)
    M = oz((n, n))
    for i in range(n):
        for j in range(n):
            M[i, j] = A[i][j]
    I_ = inv(M)
    return [[I_[i, j] for j in range(n)] for i in range(n)]


def Rvec(Ed, X, i, n):
    return [[sum(X[a] * Ed[(a, i)][p][q] for a in range(n)) for q in range(n)] for p in range(n)]


def Fbil(Ed, n, sig, X, Y, h, hp):
    """F(X,Y;h,h') = sum_i eta_i G_0(R(X,e_i).h, R(Y,e_i).h'), G_0(h,k) = tr(g^-1 h g^-1 k)"""
    gi = mat(n, lambda i, j: Fr(sig[i]) if i == j else Fr(0))
    out = Fr(0)
    for i in range(n):
        B, C = Rvec(Ed, X, i, n), Rvec(Ed, Y, i, n)
        out += sig[i] * tr(mm(mm(gi, act(B, h)), mm(gi, act(C, hp))))
    return out


def Phi(Ed, n, sig, X, Y, y):
    yi = minv(y)
    out = Fr(0)
    for i in range(n):
        B, C = Rvec(Ed, X, i, n), Rvec(Ed, Y, i, n)
        out += sig[i] * tr(mm(mm(yi, act(B, y)), mm(yi, act(C, y))))
    return out


def tsym(n, p, q):
    return mat(n, lambda i, j: Fr(int((i, j) == (p, q)) + int((i, j) == (q, p)), 2))


def unit(n, a):
    return [Fr(int(i == a)) for i in range(n)]


def traces(R, n, sig):
    Ed = endos(R, n, sig)
    e = [unit(n, a) for a in range(n)]
    S = {(p, q): tsym(n, p, q) for p in range(n) for q in range(n)}
    C1 = sum(sig[a] * sig[p] * sig[q] * Fbil(Ed, n, sig, e[a], e[a], S[p, q], S[p, q])
             for a in range(n) for p in range(n) for q in range(n))
    C2 = sum(sig[q] * Fbil(Ed, n, sig, e[a], e[b], S[a, q], S[b, q]) for a in range(n) for b in range(n) for q in range(n))
    C2p = sum(sig[q] * Fbil(Ed, n, sig, e[a], e[b], S[b, q], S[a, q]) for a in range(n) for b in range(n) for q in range(n))

    def Qf(a, b, hh):                      # the quadratic form that the leaf condition constrains
        return Fbil(Ed, n, sig, e[a], e[b], hh, hh)
    C2pol = sum(sig[q] * (Qf(a, b, madd(S[a, q], S[b, q])) - Qf(a, b, S[a, q]) - Qf(a, b, S[b, q])) / 2
                for a in range(n) for b in range(n) for q in range(n))
    return C1, C2, C2p, C2pol


def _rs(n):
    a = [[Fr(0)] * n for _ in range(n)]
    for i in range(n):
        for j in range(i, n):
            a[i][j] = a[j][i] = rq()
    return a


# ============================================================================ X4
def X4():
    say("\nX4  J5d L0: expansion of the leaf condition at y = g + eps h (n = 4)")
    n = 4
    for sig in ((1, -1, -1, -1), (1, 1, -1, -1)):
        R = acurv(n, sig)
        Ed = endos(R, n, sig)
        g = mat(n, lambda i, j: Fr(sig[i]) if i == j else Fr(0))
        gi = g
        ok_lam = True
        for _ in range(3):
            y = _rs(n)
            while minv_safe(y) is None:
                y = _rs(n)
            K = [[Fr(0)] * n for _ in range(n)]
            for i in range(n):
                for j in range(i + 1, n):
                    K[i][j] = rq()
                    K[j][i] = -K[i][j]
            B = mm(gi, K)                             # B in so(g)
            ok_lam &= tr(mm(minv(y), act(B, y))) == 0
        check(f"({sig.count(1)},{sig.count(-1)}): tr(y^-1 (B.y)) = 0 for B in so(g) at 3 random y (the lambda-term of"
              " G_lambda drops at every y)", ok_lam)
        ok = True
        for traceless in (True, False):
            h = _rs(n)
            if traceless:
                t = sum(sig[i] * h[i][i] for i in range(n))
                h[0][0] -= t * sig[0]
            X = [rq() for _ in range(n)]
            Y = [rq() for _ in range(n)]
            # y^-1 = sum_k (-eps g^-1 h)^k g^-1 as an eps-polynomial to order 3; Phi to order 2
            Hm = mm(gi, h)
            terms = [gi]
            for k in range(1, 3):
                terms.append([[-x for x in row] for row in mm(Hm, terms[-1])])
            coef = []
            for order in range(3):
                tot = Fr(0)
                for i in range(n):
                    Bx, By = Rvec(Ed, X, i, n), Rvec(Ed, Y, i, n)
                    a0, a1 = act(Bx, g), act(Bx, h)
                    b0, b1 = act(By, g), act(By, h)
                    for i1 in range(3):
                        for j1 in range(2):
                            for i2 in range(3):
                                for j2 in range(2):
                                    if i1 + j1 + i2 + j2 != order:
                                        continue
                                    P1 = mm(terms[i1], a0 if j1 == 0 else a1)
                                    P2 = mm(terms[i2], b0 if j2 == 0 else b1)
                                    tot += sig[i] * tr(mm(P1, P2))
                coef.append(tot)
            ok &= coef[0] == 0 and coef[1] == 0 and coef[2] == Fbil(Ed, n, sig, X, Y, h, h)
        check(f"({sig.count(1)},{sig.count(-1)}): Phi_(g + eps h)(X,Y) = 0 + 0 eps + F(X,Y;h,h) eps^2 + O(eps^3), for a"
              " g-traceless and for a general h (random R, X, Y)", ok)


def minv_safe(y):
    n = len(y)
    M = sp.Matrix(n, n, lambda i, j: sp.Rational(y[i][j].numerator, y[i][j].denominator))
    return None if M.det() == 0 else True


# ============================================================================ X5
def X5():
    say("\nX5  J5d L1: the traces C1, C2 of the write-up, the h <-> h' asymmetry, and the symmetrised trace")
    cases = [(4, (1, 1, 1, 1)), (4, (1, -1, -1, -1)), (4, (1, 1, -1, -1)), (3, (1, -1, -1)), (5, (1, 1, -1, -1, -1))]
    for n, sig in cases:
        t0 = time.time()
        R = acurv(n, sig)
        assert is_acurv(R, n)
        R2, Ric2 = norms(R, n, sig)
        C1, C2, C2p, C2pol = traces(R, n, sig)
        lab = f"n = {n} ({sig.count(1)},{sig.count(-1)})"
        check(f"{lab}: write-up identities C1 = (n+2)|R|^2, C2 = 3/4|R|^2 + (n+3)/2|Ric|^2 (C2 with h != h' in the"
              f" polarised form)", C1 == (n + 2) * R2 and C2 == Fr(3, 4) * R2 + Fr(n + 3, 2) * Ric2)
        Ed = endos(R, n, sig)
        X = [rq() for _ in range(n)]
        Y = [rq() for _ in range(n)]
        h, hp = _rs(n), _rs(n)
        gi = mat(n, lambda i, j: Fr(sig[i]) if i == j else Fr(0))
        H, Hp = mm(gi, h), mm(gi, hp)
        comm = lambda A, B: madd(mm(A, B), mm(B, A), -1)  # noqa: E731
        anti = -sum(sig[i] * tr(mm(comm(Rvec(Ed, X, i, n), Rvec(Ed, Y, i, n)), comm(Hp, H))) for i in range(n))
        dF = Fbil(Ed, n, sig, X, Y, h, hp) - Fbil(Ed, n, sig, X, Y, hp, h)
        check(f"{lab}: observation: F(X,Y;h,h') - F(X,Y;h',h) = -sum_i eta_i tr([R(X,e_i),R(Y,e_i)][H',H]) (H = g^-1 h)"
              f" and is nonzero ({dF}); the leaf condition constrains only the h-symmetric part", dF == anti and dF != 0)
        check(f"{lab}: C2 != C2' (h, h' swapped: {C2} vs {C2p}), so F(X,Y;h,h) = 0 does not by itself give C2 = 0",
              C2 != C2p)
        check(f"{lab}: symmetrised trace C2sym := (C2 + C2')/2 equals the polarisation of the quadratic form"
              f" F(e_a,e_b;h,h) and C2sym = (n+8)/8 |R|^2 + (n+4)/4 |Ric|^2   [{time.time() - t0:.1f}s]",
              C2pol == (C2 + C2p) / 2 and C2pol == Fr(n + 8, 8) * R2 + Fr(n + 4, 4) * Ric2)
    say("  Hence under the leaf condition C1 = 0 and C2sym = 0, so |R|^2 = |Ric|^2 = 0 (and then C2 = 0 a posteriori);")
    say("  for n = 4: 4 C2sym - C1 = 8 |Ric|^2, and = 32 Lambda^2 for an Einstein tensor.")
    n = 4
    for sig in ((1, 1, 1, 1), (1, -1, -1, -1), (1, 1, -1, -1)):
        t0 = time.time()
        W = weyl_of(acurv(n, sig), n, sig)
        Lam = rq(-4, 4)
        while Lam == 0:
            Lam = rq(-4, 4)
        R1 = {(a, b, c, d): (sig[a] if a == c else 0) * (sig[b] if b == d else 0) - (sig[a] if a == d else 0) * (sig[b] if b == c else 0)
              for a, b, c, d in itertools.product(range(n), repeat=4)}
        RE = {k: W[k] + Lam / 3 * R1[k] for k in W}
        RicE = ric_of(RE, n, sig)
        einst = all(RicE[b, d] == (Lam * sig[b] if b == d else 0) for b in range(n) for d in range(n))
        W2, Wric2 = norms(W, n, sig)
        C1, C2, C2p, C2pol = traces(RE, n, sig)
        lab = f"({sig.count(1)},{sig.count(-1)}) Einstein R = W + (Lambda/3) R_1, Lambda = {Lam}, |W|^2 = {W2}"
        check(f"{lab}: W algebraic, Ric(W) = 0, Ric(R) = Lambda g; write-up: C1 = 6|W|^2 + 16 Lambda^2,"
              f" C2 = 3/4|W|^2 + 16 Lambda^2, 8 C2 - C1 = 112 Lambda^2",
              is_acurv(W, n) and Wric2 == 0 and all(v == 0 for v in ric_of(W, n, sig).values()) and einst
              and C1 == 6 * W2 + 16 * Lam ** 2 and C2 == Fr(3, 4) * W2 + 16 * Lam ** 2 and 8 * C2 - C1 == 112 * Lam ** 2)
        check(f"{lab}: corrected: C2sym = 3/2|W|^2 + 12 Lambda^2, 4 C2sym - C1 = 32 Lambda^2 != 0   [{time.time() - t0:.1f}s]",
              C2pol == Fr(3, 2) * W2 + 12 * Lam ** 2 and 4 * C2pol - C1 == 32 * Lam ** 2 and C2pol != C2)


# ============================================================================ X6
def X6():
    say("\nX6  J5d L3, L4: positivity identity; an explicit (2,2) solution R = n (x) n")
    n = 4
    for sig in ((1, -1, -1, -1), (1, 1, -1, -1)):
        R = acurv(n, sig)
        Ed = endos(R, n, sig)
        T = [rq() for _ in range(n)]
        gTT = sum(sig[i] * T[i] ** 2 for i in range(n))
        Tb = [sig[i] * T[i] for i in range(n)]
        hT = mat(n, lambda i, j: Tb[i] * Tb[j])
        lhs = Fbil(Ed, n, sig, T, T, hT, hT)
        rhs = 0
        for a in range(n):
            B = Rvec(Ed, T, a, n)
            v = [sum(B[p][q] * T[q] for q in range(n)) for p in range(n)]
            rhs += sig[a] * sum(sig[p] * v[p] ** 2 for p in range(n))
        check(f"({sig.count(1)},{sig.count(-1)}): F(T,T;T^b (x) T^b) = 2 g(T,T) sum_a eta_a g(R(T,e_a)T, R(T,e_a)T) at a random"
              f" T (g(T,T) = {gTT}); over T^perp this is L3's sum of squares when T^perp is definite", lhs == 2 * gTT * rhs)
    # own construction of a null anti-self-dual 2-form in (2,2)
    sig = (1, 1, -1, -1)
    pairs = [(a, b) for a in range(4) for b in range(a + 1, 4)]
    star = sp.zeros(6, 6)
    for P, (a, b) in enumerate(pairs):
        for Q, (c, d) in enumerate(pairs):
            # (*w)_ab = (1/2) eps_ab^cd w_cd = sum_{c<d} eps_abcd eta_c eta_d w_cd
            star[P, Q] = sp.LeviCivita(a, b, c, d) * sig[c] * sig[d]
    ip = sp.diag(*[sig[a] * sig[b] for a, b in pairs])
    asd = (star + sp.eye(6)).nullspace()
    found = None
    for coefs in itertools.product(range(-2, 3), repeat=len(asd)):
        if all(c == 0 for c in coefs):
            continue
        v = sum((c * w for c, w in zip(coefs, asd)), sp.zeros(6, 1))
        if (v.T * ip * v)[0, 0] == 0:
            found = v
            break
    nf = {}
    for P, (a, b) in enumerate(pairs):
        nf[a, b] = Fr(int(found[P].p), int(found[P].q))
        nf[b, a] = -nf[a, b]
    for a in range(4):
        nf[a, a] = Fr(0)
    c0 = Fr(5, 3)
    R = {(a, b, c, d): c0 * nf[a, b] * nf[c, d] for a, b, c, d in itertools.product(range(4), repeat=4)}
    Ed = endos(R, 4, sig)
    okF, okP = True, True
    for _ in range(4):
        X = [rq() for _ in range(4)]
        Y = [rq() for _ in range(4)]
        h, hp = _rs(4), _rs(4)
        okF &= Fbil(Ed, 4, sig, X, Y, h, h) == 0 and Fbil(Ed, 4, sig, X, Y, h, hp) == 0
        y = _rs(4)
        while minv_safe(y) is None:
            y = _rs(4)
        okP &= Phi(Ed, 4, sig, X, Y, y) == 0
    Ric = ric_of(R, 4, sig)
    check(f"(2,2): n = {list(found)} on e^a^e^b (a<b) is anti-self-dual (*n = -n) and null; R = (5/3) n (x) n is an"
          " algebraic curvature tensor with Ric = 0 (Lambda = 0)",
          star * found == -found and (found.T * ip * found)[0, 0] == 0 and is_acurv(R, 4) and all(v == 0 for v in Ric.values())
          and any(v != 0 for v in R.values()))
    check("(2,2): for this R, F(X,Y;h,h) = 0 and even the polarised F(X,Y;h,h') = 0 at 4 random (X,Y,h,h'), and the full"
          " condition Phi_y(X,Y) = 0 at 4 random nondegenerate y (L4)", okF and okP)
    Rc = {(a, b, c, d): (sig[a] if a == c else 0) * (sig[b] if b == d else 0) - (sig[a] if a == d else 0) * (sig[b] if b == c else 0)
          for a, b, c, d in itertools.product(range(4), repeat=4)}
    Edc = endos(Rc, 4, sig)
    e0 = unit(4, 0)
    hc = _rs(4)
    check("(2,2) control: constant curvature fails the leaf condition, F(e_0,e_0;h,h) != 0 for a random h",
          Fbil(Edc, 4, sig, e0, e0, hc, hc) != 0)


def main():
    say("xcheck_round7.py  (REVIEW5-GU; exact Fractions/sympy; seed random.Random(20261107); no repository imports)")
    X0()
    X1()
    X2()
    X3()
    X4()
    X5()
    X6()
    nf = [lab for lab, ok in RES if not ok]
    say(f"\nSUMMARY: {len(RES)} checks, {len(RES) - len(nf)} True, {len(nf)} False")
    for lab in nf:
        say(f"  FALSE: {lab}")
    say(f"runtime {time.time() - T0:.0f} s")


if __name__ == '__main__':
    main()
