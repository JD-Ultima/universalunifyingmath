"""GU round 2 (F): the Zorro metric on Y = Met(X^n) and the 14 -> 4 reduction (gaps GU-G06, GU-G14).
Write-up: formalization/geometric-unity/F-zorro-reduction.md.  Claims: GU-006, GU-021, GU-024, GU-054.

Reconstruction (R-Z): a metric g on X (the draft's gimel) -> its Levi-Civita connection (aleph) -> the induced
linear connection on Sym^2 T*X -> horizontal distribution H on Y (Dy = dy - (Gamma y + y Gamma) dx) ->
    Z1:  g_aleph = cH * y(pi_*, pi_*)  + cV * G_lambda|_y(Dy, Dy)      (tautological horizontal metric, (3.7))
    Z2:  g_aleph = cH * g(pi_*, pi_*)  + cV * G_lambda|_y(Dy, Dy)      (base metric pulled back; Riemannian submersion)
G_lambda(h,k) = tr(y^-1 h y^-1 k) - lambda tr(y^-1 h) tr(y^-1 k).  Section s = g : X -> Y.

Everything exact: coefficients in Q(lambda, cH, cV) (and, for n = 2, also the 3-jet of a conformal factor phi).
Metric 2-jets at a point of Y are built by truncated Taylor arithmetic (F_jets.py); base metrics are
polynomials with random rational coefficients in general (non-normal) coordinates.

  R0  controls: round sphere / hyperbolic space (n = 3): Ric = +-(n-1) g.
  R1  section: s* g_aleph = cH g; ds(T X) = H along s(X) (connection coefficients = dg); V perp H.
  R2  n = 2, fully general metric g = e^{2 phi} eta (phi symbolic 3-jet), Riemannian and Lorentzian:
      full Riemann tensor of Y on s(X); Gauss equation; vertical block = fiber curvature; mixed sum; blocks.
  R3  n = 3, 4, 5, random rational metrics, several signatures, both variants: second fundamental form of s(X),
      Ricci blocks HH/HV/VV and scalar curvature of g_aleph on s(X) vs closed forms (Theorems F2, F3).
  R4  control off the section (Z2): s_Y = R/cH + s_F - |A|^2 (O'Neill, Besse 9.37) at y0 != g(x0).
  R5  the 4D equations and lambda-dependence (n = 4): cosmological terms from (P) projection and (A) reduced action.
"""
import itertools
import random
import sys
import time
from fractions import Fraction
from F_jets import (Ctx, jzero, jconst, jadd, jmul, jscale, jtrunc, jval, jgrad, jhess, finv, fmatmul,
                    base_from_coeffs, base_geometry, total_metric, frame_blocks, ricci_at_point,
                    riemann_at_point, trace, sym_index, supp, metric_point_data, christoffel)

T0 = time.time()
random.seed(20210401)
ctx = Ctx('lam,cH,cV')
lam, cH, cV = ctx.gens


def P(*a):
    print(*a)
    sys.stdout.flush()


def F(x, y=1):
    return Fraction(x, y)


# ------------------------------------------------------------------ helpers
def basis_mats(n, idx, C):
    mats = []
    for (i, j) in idx:
        M = [[C.zero] * n for _ in range(n)]
        for (a, b) in supp((i, j)):
            M[a][b] = C.one
        mats.append(M)
    return mats


def tr(M, C):
    return sum((M[i][i] for i in range(len(M))), C.zero)


def G_lam(yi, h, k, lv, C):
    return tr(fmatmul(fmatmul(yi, h, C), fmatmul(yi, k, C), C), C) - lv * tr(fmatmul(yi, h, C), C) * tr(fmatmul(yi, k, C), C)


def vec_to_sym(v, n, idx, C):
    M = [[C.zero] * n for _ in range(n)]
    for A, (i, j) in enumerate(idx):
        M[i][j] = v[A]
        M[j][i] = v[A]
    return M


def closed_forms(n, variant, C, lamv, cHv, cVv):
    """Predicted values on s(X) (Theorems F2/F3)."""
    nn = C.c(F(n))
    mu = lamv / (1 - nn * lamv)
    sF = C.c(F(-n * (n - 1) * (n + 2), 8))           # times 1/cV
    if variant == 'Z2':
        return dict(gammaHH=C.zero, VV_extra=C.zero, sigma=sF, II_coef=C.zero, mu=mu)
    sigma = sF - nn / 4 * ((1 + lamv) / (1 - nn * lamv) + (nn + 1) / 2)
    return dict(gammaHH=cHv / (4 * cVv * (nn * lamv - 1)), VV_extra=C.c(F(-1, 4)), sigma=sigma,
                II_coef=-cHv / (2 * cVv), mu=mu)


def second_fundamental_form(G, m, n, idx, Psi0, gjets, C):
    """II(e_mu, e_nu) = vertical (frame) part of nabla_{E_mu} E_nu at z0, E_nu = s_* d_nu."""
    g0, dG, ddG = metric_point_data(G, m, C)
    g0i = finv(g0, C)
    _, Gam = christoffel(g0i, dG, m, C)
    N = len(idx)
    II = {}
    for mu in range(n):
        for nu in range(n):
            Emu = [C.one if a == mu else C.zero for a in range(n)] + [Psi0[A][mu] for A in range(N)]
            Enu = [C.one if a == nu else C.zero for a in range(n)] + [Psi0[A][nu] for A in range(N)]
            v = [C.zero] * m
            for A, (i, j) in enumerate(idx):
                v[n + A] = jhess(gjets[i][j], mu, nu, C)
            for c in range(m):
                v[c] += sum((Gam[c][a][b] * Emu[a] * Enu[b] for a in range(m) for b in range(m)
                             if Emu[a] != 0 and Enu[b] != 0), C.zero)
            II[(mu, nu)] = [v[n + A] - sum((Psi0[A][r] * v[r] for r in range(n)), C.zero) for A in range(N)]
    return II


def predicted_II(n, idx, g0, cf, C):
    out = {}
    for mu in range(n):
        for nu in range(n):
            vals = []
            for (i, j) in idx:
                sym = (g0[i][mu] * g0[j][nu] + g0[j][mu] * g0[i][nu]) / 2
                vals.append(cf['II_coef'] * (sym + cf['mu'] * g0[mu][nu] * g0[i][j]))
            out[(mu, nu)] = vals
    return out


def check_section(n, idx, Psi0, G, base, gjets, C):
    """R1: s*g_aleph = cH g (frame HH block), horizontal lift = s_* (Psi0 = dg), V perp H at z0."""
    ok_psi = all(Psi0[A][mu] == jgrad(gjets[i][j], mu, C) for A, (i, j) in enumerate(idx) for mu in range(n))
    G0 = [[jval(G[a][b], C) for b in range(len(G))] for a in range(len(G))]
    Gf = frame_blocks(G0, n, idx, Psi0, C)
    ok_pull = all(Gf[a][b] == cH_of(C) * base['g0'][a][b] for a in range(n) for b in range(n))
    ok_perp = all(Gf[a][n + B] == 0 for a in range(n) for B in range(len(idx)))
    return ok_psi, ok_pull, ok_perp, Gf


def cH_of(C):
    return C.gens[1]


# ------------------------------------------------------------------ R0 controls
P("F_zorro_reduction.py  (exact; coefficients in Q(lambda, cH, cV))\n")
P("R0  controls")
for sgn, name in [(1, 'sphere'), (-1, 'hyperbolic')]:
    n = 3
    G = [[jzero(2) for _ in range(n)] for _ in range(n)]
    for i in range(n):
        G[i][i][0][()] = ctx.c(4)
        for k in range(n):
            G[i][i][2][(k, k)] = ctx.c(-8 * sgn)
    Ric, g0, g0i = ricci_at_point(G, n, ctx)
    P(f"  {name} 4|dx|^2/(1 {'+' if sgn > 0 else '-'} |x|^2)^2, n=3: Ric == {sgn}*(n-1)*g: "
      f"{all(Ric[a][b] == sgn * 2 * g0[a][b] for a in range(n) for b in range(n))}")
P("")

# ------------------------------------------------------------------ R2  n = 2 fully general (symbolic conformal factor)
P("R2  n = 2: general metric g = e^{2 phi} eta near a point (isothermal / null coordinates), phi(0) = 0,")
P("    phi = p_i x^i + p_ij x^i x^j / 2 + p_ijk x^i x^j x^k / 6 with all p symbolic. Full Riemann of Y^5 on s(X).")
C2 = Ctx('lam,cH,cV,p0,p1,p00,p01,p11,p000,p001,p011,p111')
l2, h2, v2, p0, p1, p00, p01, p11, p000, p001, p011, p111 = C2.gens
pd1 = {(0,): p0, (1,): p1}
pd2 = {(0, 0): p00 / 2, (0, 1): p01, (1, 1): p11 / 2}
pd3 = {(0, 0, 0): p000 / 6, (0, 0, 1): p001 / 2, (0, 1, 1): p011 / 2, (1, 1, 1): p111 / 6}
phi = [dict(), pd1, pd2, pd3]
two_phi = jscale(phi, 2)
expo = jconst(C2.one, 3)
term = jconst(C2.one, 3)
for k in range(1, 4):
    term = jscale(jmul(term, two_phi, 3), C2.c(F(1, k)))
    expo = jadd(expo, term)
for eps in ([1, 1], [1, -1]):
    n = 2
    gj = [[jscale(expo, C2.c(eps[i])) if i == j else jzero(3) for j in range(n)] for i in range(n)]
    base = base_geometry(gj, n, C2)
    R_pred = -2 * (eps[0] * p00 + eps[1] * p11)       # R = -2 e^{-2phi} eta^{ij} d_i d_j phi at 0
    P(f"  eta = diag{tuple(eps)}:  R(g)(0) = {base['Rs']}   (== -2 eta^ij p_ij: {base['Rs'] == R_pred})")
    for variant in ('Z2', 'Z1'):
        t = time.time()
        G, m, idx, Psi0 = total_metric(gj, base, n, base['g0'], l2, h2, v2, variant, C2)
        ok_psi, ok_pull, ok_perp, Gf = check_section(n, idx, Psi0, G, base, gj, C2)
        R, g0Y, g0Yi = riemann_at_point(G, m, C2)
        Ef = [[C2.zero] * m for _ in range(m)]
        for mu in range(n):
            Ef[mu][mu] = C2.one
            for A in range(len(idx)):
                Ef[mu][n + A] = Psi0[A][mu]
        for A in range(len(idx)):
            Ef[n + A][n + A] = C2.one
        # frame components of Riemann (sequential contraction)
        T = R
        for slot in range(4):
            T2 = {}
            for key in itertools.product(range(m), repeat=4):
                acc = C2.zero
                for a in range(m):
                    if Ef[key[slot]][a] != 0:
                        kk = list(key)
                        kk[slot] = a
                        acc += Ef[key[slot]][a] * T[tuple(kk)]
                T2[key] = acc
            T = T2
        Rf = T
        gf = Gf
        gfi = finv(gf, C2)
        cf = closed_forms(n, variant, C2, l2, h2, v2)
        # second fundamental form
        II = second_fundamental_form(G, m, n, idx, Psi0, gj, C2)
        IIp = predicted_II(n, idx, base['g0'], cf, C2)
        ok_II = all(II[k] == IIp[k] for k in II)

        def ip(u, w):
            return sum((u[A] * gf[n + A][n + B] * w[B] for A in range(len(idx)) for B in range(len(idx))), C2.zero)
        # Gauss: R^Y_abcd = cH R^g_abcd - <II_ca, II_db> + <II_cb, II_da>   (tangential indices)
        ok_gauss = all(Rf[(a, b, c, d)] == h2 * base['R'][(a, b, c, d)] - ip(II[(c, a)], II[(d, b)]) + ip(II[(c, b)], II[(d, a)])
                       for a, b, c, d in itertools.product(range(n), repeat=4))
        # vertical block = fiber curvature  -(cV/4) tr([[H_c,H_d],H_b] H_a)
        y0i = finv(base['g0'], C2)
        Hm = [fmatmul(y0i, M, C2) for M in basis_mats(n, idx, C2)]

        def comm(a, b):
            ab, ba = fmatmul(a, b, C2), fmatmul(b, a, C2)
            return [[ab[i][j] - ba[i][j] for j in range(n)] for i in range(n)]
        N = len(idx)
        ok_vert = all(Rf[(n + a, n + b, n + c, n + d)] == -v2 / 4 * tr(fmatmul(comm(comm(Hm[c], Hm[d]), Hm[b]), Hm[a], C2), C2)
                      for a, b, c, d in itertools.product(range(N), repeat=4))
        s_mix = sum((gfi[mu][rho] * gfi[n + A][n + B] * Rf[(mu, n + A, rho, n + B)]
                     for mu in range(n) for rho in range(n) for A in range(N) for B in range(N)), C2.zero)
        nn = C2.c(F(n))
        s_mix_pred = C2.zero if variant == 'Z2' else -(nn * (nn + 1) / 2 + nn * l2 / (1 - nn * l2)) / (4 * v2)
        odd = sum(1 for key, v in Rf.items() if sum(1 for k in key if k >= n) % 2 == 1 and v != 0)
        Ricf = [[sum((gfi[a][c] * Rf[(c, b, a, d)] for a in range(m) for c in range(m) if gfi[a][c] != 0), C2.zero)
                 for d in range(m)] for b in range(m)]
        s = sum((gfi[a][b] * Ricf[a][b] for a in range(m) for b in range(m) if gfi[a][b] != 0), C2.zero)
        ok_s = (s == base['Rs'] / h2 + cf['sigma'] / v2)
        ok_HH = all(Ricf[a][b] == base['Ric'][a][b] + cf['gammaHH'] * base['g0'][a][b] for a in range(n) for b in range(n))
        ok_HV = all(Ricf[a][n + B] == 0 for a in range(n) for B in range(N))
        Ev = basis_mats(n, idx, C2)
        ok_VV = all(Ricf[n + A][n + B] == -F(n, 4) * G_lam(y0i, Ev[A], Ev[B], C2.c(F(1, n)), C2)
                    + cf['VV_extra'] * G_lam(y0i, Ev[A], Ev[B], C2.zero, C2) for A in range(N) for B in range(N))
        P(f"   {variant}: s_*=horizontal lift {ok_psi}; s*g_aleph = cH g {ok_pull}; V perp H {ok_perp}; "
          f"II == closed form {ok_II}; Gauss eq. {ok_gauss}; vertical block = fiber curvature {ok_vert}")
        P(f"        mixed scalar sum = {s_mix} (pred {s_mix_pred}: {s_mix == s_mix_pred}); "
          f"#nonzero frame components with odd # vertical indices = {odd}")
        P(f"        Ric_HH {ok_HH}, Ric_HV = 0 {ok_HV}, Ric_VV {ok_VV}, s_Y = R/cH + sigma/cV {ok_s}:  s_Y = {s}"
          f"   [{time.time() - t:.1f}s]")
P("")

# ------------------------------------------------------------------ R3  n = 3, 4, 5
P("R3  n = 3, 4, 5: random rational polynomial metrics (3-jets) in general coordinates; Ricci blocks on s(X)")


def rsym(n, s=3):
    M = [[None] * n for _ in range(n)]
    for i in range(n):
        for j in range(i, n):
            v = ctx.c(F(random.randint(-s, s), random.randint(1, 3)))
            M[i][j] = v
            M[j][i] = v
    return M


def random_metric(n, eps, nondiag):
    g0 = [[ctx.c(eps[i]) if i == j else ctx.zero for j in range(n)] for i in range(n)]
    if nondiag:
        Pm = [[ctx.c(F(random.randint(-1, 1))) + (3 if i == j else 0) for j in range(n)] for i in range(n)]
        PT = [list(r) for r in zip(*Pm)]
        g0 = fmatmul(fmatmul(PT, g0, ctx), Pm, ctx)
    lin = [rsym(n) for _ in range(n)]
    quad = {(a, b): rsym(n) for a in range(n) for b in range(a, n)}
    cub = {(a, b, c): rsym(n) for a in range(n) for b in range(a, n) for c in range(b, n)}
    return base_from_coeffs(g0, lin, quad, cub, n, ctx)


results = []
runs = [(3, [1, 1, 1], True), (3, [1, 1, -1], True),
        (4, [1, -1, -1, -1], True), (4, [1, 1, 1, -1], False), (4, [1, 1, 1, 1], True), (4, [1, 1, -1, -1], False),
        (5, [1, 1, 1, 1, 1], False)]
for n, eps, nondiag in runs:
    gj = random_metric(n, eps, nondiag)
    base = base_geometry(gj, n, ctx)
    sig = (eps.count(1), eps.count(-1))
    nab_nonzero = any(base['nablaRic'][r][a][b] != 0 for r in range(n) for a in range(n) for b in range(n))
    for variant in ('Z2', 'Z1'):
        t = time.time()
        G, m, idx, Psi0 = total_metric(gj, base, n, base['g0'], lam, cH, cV, variant, ctx)
        ok_psi, ok_pull, ok_perp, Gf = check_section(n, idx, Psi0, G, base, gj, ctx)
        cf = closed_forms(n, variant, ctx, lam, cH, cV)
        II = second_fundamental_form(G, m, n, idx, Psi0, gj, ctx)
        ok_II = all(II[k] == v for k, v in predicted_II(n, idx, base['g0'], cf, ctx).items())
        Ric, G0, G0i = ricci_at_point(G, m, ctx)
        s = trace(Ric, G0i, ctx)
        RF = frame_blocks(Ric, n, idx, Psi0, ctx)
        N = len(idx)
        g0 = base['g0']
        y0i = finv(g0, ctx)
        Ev = basis_mats(n, idx, ctx)
        dHH = [[RF[a][b] - base['Ric'][a][b] for b in range(n)] for a in range(n)]
        ok_HH = all(dHH[a][b] == cf['gammaHH'] * g0[a][b] for a in range(n) for b in range(n))
        ok_HV = all(RF[a][n + B] == 0 for a in range(n) for B in range(N))
        ok_VV = all(RF[n + A][n + B] == -F(n, 4) * G_lam(y0i, Ev[A], Ev[B], ctx.c(F(1, n)), ctx)
                    + cf['VV_extra'] * G_lam(y0i, Ev[A], Ev[B], ctx.zero, ctx) for A in range(N) for B in range(N))
        sres = s - base['Rs'] / cH
        ok_s = (sres == cf['sigma'] / cV)
        P(f"  n={n} sig={sig} g0 {'non-diag' if nondiag else 'diag'} {variant}: s_*=lift {ok_psi}, s*g=cH g {ok_pull}, "
          f"V perp H {ok_perp}, II {ok_II}; Ric_HH - Ric_g = {cf['gammaHH']}*g {ok_HH}; Ric_HV = 0 {ok_HV} "
          f"(nabla Ric != 0: {nab_nonzero}); Ric_VV {ok_VV}; s_Y - R/cH = {sres} {ok_s}   [{time.time() - t:.1f}s]")
        results.append((n, sig, variant, ok_psi and ok_pull and ok_perp and ok_II and ok_HH and ok_HV and ok_VV and ok_s))
P(f"  R3 all checks pass: {all(r[3] for r in results)}\n")

# ------------------------------------------------------------------ R4 O'Neill control off the section
P("R4  control (Z2) at a point y0 != g(x0) of the fiber: s_Y = R/cH + s_F - |A|^2 (O'Neill; Besse 9.37),")
P("    |A|^2 = sum_{ijkl} (cH g)^{ik} (cH g)^{jl} cV G_lambda(A_ij, A_kl),  A_ij = (1/2) R^{Sym2}(d_i,d_j) y0")
for n, eps in [(2, [1, -1]), (3, [1, 1, -1])]:
    gj = random_metric(n, eps, True)
    base = base_geometry(gj, n, ctx)
    g0 = base['g0']
    while True:
        y0 = [[g0[i][j] + ctx.c(F(random.randint(-2, 2), 3)) * (1 if i <= j else 0) for j in range(n)] for i in range(n)]
        for i in range(n):
            for j in range(i):
                y0[i][j] = y0[j][i]
        try:
            y0i = finv(y0, ctx)
            break
        except StopIteration:
            pass
    G, m, idx, Psi0 = total_metric(gj, base, n, y0, lam, cH, cV, 'Z2', ctx)
    Ric, G0, G0i = ricci_at_point(G, m, ctx)
    s = trace(Ric, G0i, ctx)
    gi = base['g0i']
    Rup = {(c, a, i, j): sum((gi[c][e] * base['R'][(e, a, i, j)] for e in range(n)), ctx.zero)
           for c, a, i, j in itertools.product(range(n), repeat=4)}
    Aij = {}
    for i in range(n):
        for j in range(n):
            Aij[(i, j)] = [[-(sum((Rup[(c, a, i, j)] * y0[c][b] + Rup[(c, b, i, j)] * y0[a][c] for c in range(n)), ctx.zero)) / 2
                            for b in range(n)] for a in range(n)]
    A2 = sum((gi[i][k] * gi[j][l] * G_lam(y0i, Aij[(i, j)], Aij[(k, l)], lam, ctx)
              for i, j, k, l in itertools.product(range(n), repeat=4)), ctx.zero) * cV / (cH * cH)
    sF = ctx.c(F(-n * (n - 1) * (n + 2), 8)) / cV
    P(f"  n={n}: s_Y(y0) == R/cH + s_F - |A|^2: {s == base['Rs'] / cH + sF - A2};  |A|^2 = {A2}")
P("")

# ------------------------------------------------------------------ R5 4D equations, lambda dependence
P("R5  4D equations on s(X) (n general, then n = 4).  sigma := cV * (s_Y - R/cH) on s(X).")
P("    (P) projection: s*(G^Y) = G(g) + Lambda_P g      (horizontal block of the Einstein tensor of g_aleph)")
P("    (A) reduced action  S[g] = int_X s*(s_Y) dvol(s* g_aleph) = cH^(n/2) int (R/cH + sigma/cV) dvol_g")
P("        => Euler-Lagrange: G(g) + Lambda_A g = 0,   Lambda_A = -cH sigma / (2 cV).")
for variant in ('Z2', 'Z1'):
    for n in (4,):
        cf = closed_forms(n, variant, ctx, lam, cH, cV)
        LA = -cH * cf['sigma'] / (2 * cV)
        LP = cf['gammaHH'] - cH * cf['sigma'] / (2 * cV)
        P(f"  {variant}, n={n}: sigma = {cf['sigma']}")
        P(f"         Lambda_A = {LA}")
        P(f"         Lambda_P = {LP}")
        P(f"         Lambda_P - Lambda_A = {LP - LA}")
        for lv in (F(0), F(1, 2), F(1), F(5, 18), F(4, 15)):
            cfv = closed_forms(n, variant, ctx, ctx.c(lv), cH, cV)
            la = -cH * cfv['sigma'] / (2 * cV)
            lp = cfv['gammaHH'] - cH * cfv['sigma'] / (2 * cV)
            P(f"           lambda = {str(lv):5s}: Lambda_A = {la},  Lambda_P = {lp}")
# zeros of Lambda_A, Lambda_P for Z1 as functions of lambda (n = 4)
import sympy as sp
L, ch, cv = sp.symbols('lam cH cV')
for nm, expr in [('Lambda_A', sp.Rational(5, 4) * ch / cv * (5 - 18 * L) / (1 - 4 * L)),
                 ('Lambda_P', sp.Rational(3, 2) * ch / cv * (4 - 15 * L) / (1 - 4 * L))]:
    cf = closed_forms(4, 'Z1', ctx, lam, cH, cV)
    mine = (-cH * cf['sigma'] / (2 * cV)) if nm == 'Lambda_A' else (cf['gammaHH'] - cH * cf['sigma'] / (2 * cV))
    same = sp.simplify(sp.sympify(str(mine), locals={'lam': L, 'cH': ch, 'cV': cv}) - expr) == 0
    P(f"  Z1 n=4 {nm} == {expr}: {same};  zero at lambda = {sp.solve(sp.numer(sp.together(expr)), L)}")
P("  General n (Z1):  Lambda_A = cH n/(16 cV) [n^2+2n-1 + 2(1+lam)/(1-n lam)],")
P("                   Lambda_P = cH/(16 cV) [n(n^2+2n-1) + (2n(1+lam)-4)/(1-n lam)];  Z2: Lambda_A = Lambda_P = cH n(n-1)(n+2)/(16 cV).")
ok_gen = True
for n in (2, 3, 4, 5, 6):
    nn = ctx.c(F(n))
    cf = closed_forms(n, 'Z1', ctx, lam, cH, cV)
    LA = -cH * cf['sigma'] / (2 * cV)
    LP = cf['gammaHH'] - cH * cf['sigma'] / (2 * cV)
    ok_gen &= (LA == cH * nn / (16 * cV) * (nn * nn + 2 * nn - 1 + 2 * (1 + lam) / (1 - nn * lam)))
    ok_gen &= (LP == cH / (16 * cV) * (nn * (nn * nn + 2 * nn - 1) + (2 * nn * (1 + lam) - 4) / (1 - nn * lam)))
P(f"  general-n Lambda formulas consistent with the closed forms (n = 2..6): {ok_gen}")
P("  Sign under the (7,7) choice (Theorem A4: lambda > 1/4, sign(cV cH) = +1 for base (1,3), -1 for (3,1)):")
P("    vacuum equation G + Lambda g = 0; physical cosmological constant Lambda_phys = -Lambda in (+---), +Lambda in (-+++).")
P("    Z2: Lambda_phys = -9|cH/cV|/2 < 0 (anti-de Sitter sign) in both conventions, for every lambda.")
for lv in (F(13, 50), F(27, 100), F(1, 2), F(1), F(10)):
    cfv = closed_forms(4, 'Z1', ctx, ctx.c(lv), ctx.one, ctx.one)
    la = -cfv['sigma'] / 2
    lp = cfv['gammaHH'] - cfv['sigma'] / 2
    P(f"    Z1, lambda = {str(lv):6s}: Lambda_A*(cV/cH) = {la}, Lambda_P*(cV/cH) = {lp}  "
      f"-> Lambda_phys sign in (+---) with cH = cV = 1 (A, P) = ({'-' if sp.sympify(str(la)) > 0 else '+'}, "
      f"{'-' if sp.sympify(str(lp)) > 0 else '+'})")
P(f"\nTotal time {time.time() - T0:.1f}s")
