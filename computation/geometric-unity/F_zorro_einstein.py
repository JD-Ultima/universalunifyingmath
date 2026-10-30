"""GU round 2 (F): Einstein condition for the Zorro metric along the section, and two auxiliary checks.
Write-up: formalization/geometric-unity/F-zorro-reduction.md (Theorem F4, Remarks).  Uses F_jets.py.

  E1  n = 1 sanity check: Y = Met(R) = R_{>0}; Z1 metric cH y dx^2 + cV (1-lambda) dy^2/y^2 is a hyperbolic
      plane of curvature K = -1/(4 cV (1-lambda)); compare with the closed forms of Theorem F3 at n = 1.
  E2  Theorem F4 (Z1): along s(X), Ric(g_aleph) = kappa g_aleph  <=>  lambda = 1/(n+1) and Ric(g) = 0,
      kappa = -(n+1)/(4 cV).  Test n = 4, base (1,3):
        (i)   Ricci-flat-at-the-point metric g = eta - W x x / 3 (W a random exact Weyl tensor), lambda = 1/5: Einstein?
        (ii)  same metric, lambda = 1/2 (draft-literal) and lambda = 1: Einstein?
        (iii) a random non-Ricci-flat metric, lambda = 1/5: Einstein?
        (iv)  Z2, lambda = 1/5 and 1/2: Einstein?
        (v)   off the section: (x0, y0) with y0 != g(x0), lambda = 1/5: Ric(g_aleph) proportional to g_aleph?
  E4  explanation of E2(v)-control: for flat aleph, Z1 at lambda = 1/(n+1), cH = 2 cV is locally isometric to the
      symmetric space SL(n+1,R)/SO(p',q') with metric cV tr((M^-1 dM)^2): with M = L^T diag(u, 1/det u) L,
      L = [[1, a],[0, 1]], and y = det(u) u, x = a:  tr((M^-1 dM)^2) = G_{1/(n+1)}|_y(dy, dy) + 2 y(dx, dx).
      Exact: n = 2 symbolic; n = 3, 4 at 10 random rational points (sympy).
  E3  scale: fiber length l_F = sqrt|cV/cH| needed for |Lambda_F| = 9/(2 l_F^2) (Z2, n = 4) to equal the observed
      cosmological constant (input: Lambda_obs = 1.1056e-52 m^-2, Planck 2018 VI); floating point.
"""
import itertools
import math
import random
import sys
import time
from fractions import Fraction as F
from F_jets import (Ctx, jzero, finv, fmatmul, base_from_coeffs, base_geometry, total_metric, frame_blocks,
                    ricci_at_point, trace, sym_index)

T0 = time.time()
random.seed(20210401)
ctx = Ctx('lam,cH,cV')
lam, cH, cV = ctx.gens


def P(*a):
    print(*a)
    sys.stdout.flush()


def zero_cub(n):
    return {}


def poly_metric(n, g0, quad):
    lin = [[[ctx.zero] * n for _ in range(n)] for _ in range(n)]
    return base_from_coeffs(g0, lin, quad, {}, n, ctx)


P("F_zorro_einstein.py (exact unless stated)\n")
# ------------------------------------------------------------------ E1
P("E1  n = 1: Met(R) with Z1 metric cH y dx^2 + cV (1-lambda) (dy/y)^2  (hyperbolic plane, K = -1/(4 cV (1-lambda)))")
n = 1
g0 = [[ctx.one]]
gj = poly_metric(1, g0, {})
base = base_geometry(gj, 1, ctx)
G, m, idx, Psi0 = total_metric(gj, base, 1, g0, lam, cH, cV, 'Z1', ctx)
Ric, G0, G0i = ricci_at_point(G, m, ctx)
s = trace(Ric, G0i, ctx)
K = -1 / (4 * cV * (1 - lam))
sigma1 = -ctx.c(F(1, 4)) * ((1 + lam) / (1 - lam) + 1)            # Theorem F3 closed form at n = 1 (s_F = 0)
P(f"  s_Y = {s};  2K = {2 * K};  equal: {s == 2 * K};  closed form sigma_1(n=1)/cV = {sigma1 / cV}: {s == sigma1 / cV}")
P(f"  Ric_HH = {Ric[0][0]} == cH/(4 cV (lambda-1)) g: {Ric[0][0] == cH / (4 * cV * (lam - 1))};  "
  f"Ric_VV = {Ric[1][1]} == -1/4 G_0 = -1/4: {Ric[1][1] == ctx.c(F(-1, 4))}\n")

# ------------------------------------------------------------------ E2
P("E2  Einstein condition along s(X), n = 4, base signature (1,3)")
n = 4
eta = [1, -1, -1, -1]
E = [[ctx.c(eta[i]) if i == j else ctx.zero for j in range(n)] for i in range(n)]


def kn(h, k):
    return {(a, b, c, d): h[a][c] * k[b][d] + h[b][d] * k[a][c] - h[a][d] * k[b][c] - h[b][c] * k[a][d]
            for a, b, c, d in itertools.product(range(n), repeat=4)}


def rsym(s=3):
    M = [[None] * n for _ in range(n)]
    for i in range(n):
        for j in range(i, n):
            M[i][j] = M[j][i] = ctx.c(F(random.randint(-s, s), random.randint(1, 3)))
    return M


R = {key: ctx.zero for key in itertools.product(range(n), repeat=4)}
for _ in range(3):
    T = kn(rsym(), rsym())
    for key in R:
        R[key] += T[key]
Ei = finv(E, ctx)
Ric0 = [[sum((Ei[a][c] * R[(a, b, c, d)] for a in range(n) for c in range(n)), ctx.zero) for d in range(n)] for b in range(n)]
s0 = sum((Ei[b][d] * Ric0[b][d] for b in range(n) for d in range(n)), ctx.zero)
Psch = [[(Ric0[a][b] - s0 / (2 * (n - 1)) * E[a][b]) / (n - 2) for b in range(n)] for a in range(n)]
KNP = kn(Psch, E)
W = {key: R[key] - KNP[key] for key in R}
Wric = [[sum((Ei[a][c] * W[(a, b, c, d)] for a in range(n) for c in range(n)), ctx.zero) for d in range(n)] for b in range(n)]
P(f"  random exact Weyl tensor W: all traces vanish: {all(Wric[b][d] == 0 for b in range(n) for d in range(n))}; "
  f"W != 0: {any(v != 0 for v in W.values())}")
# g_{mu nu}(x) = eta_{mu nu} - (1/3) W_{mu a nu b} x^a x^b
quad = {}
for a in range(n):
    for b in range(a, n):
        M = [[(-(W[(mu, a, nu, b)] + W[(mu, b, nu, a)]) / 3 if a != b else -W[(mu, a, nu, a)] / 3)
              for nu in range(n)] for mu in range(n)]
        quad[(a, b)] = M
g_weyl = poly_metric(n, E, quad)
b_weyl = base_geometry(g_weyl, n, ctx)
P(f"  base metric eta - W x x/3: Ric(g)(0) = 0: {all(b_weyl['Ric'][a][b] == 0 for a in range(n) for b in range(n))}; "
  f"Riem(g)(0) != 0: {any(v != 0 for v in b_weyl['R'].values())}")


def einstein_test(gj, base, y0, lamv, variant):
    G, m, idx, Psi0 = total_metric(gj, base, n, y0, lamv, cH, cV, variant, ctx)
    Ric, G0, G0i = ricci_at_point(G, m, ctx)
    s = trace(Ric, G0i, ctx)
    rho = s / m
    resid = [[Ric[a][b] - rho * G0[a][b] for b in range(m)] for a in range(m)]
    nz = sum(1 for a in range(m) for b in range(m) if resid[a][b] != 0)
    return nz == 0, rho, nz, Ric, G0, idx, Psi0


kappa = -ctx.c(F(n + 1, 4)) / cV
for label, gj, base, lamv, variant in [
        ("(i)   Weyl metric, Z1, lambda = 1/5", g_weyl, b_weyl, ctx.c(F(1, 5)), 'Z1'),
        ("(ii)  Weyl metric, Z1, lambda = 1/2", g_weyl, b_weyl, ctx.c(F(1, 2)), 'Z1'),
        ("(ii)  Weyl metric, Z1, lambda = 1  ", g_weyl, b_weyl, ctx.one, 'Z1'),
        ("(iv)  Weyl metric, Z2, lambda = 1/5", g_weyl, b_weyl, ctx.c(F(1, 5)), 'Z2'),
        ("(iv)  Weyl metric, Z2, lambda = 1/2", g_weyl, b_weyl, ctx.c(F(1, 2)), 'Z2')]:
    t = time.time()
    ok, rho, nz, *_ = einstein_test(gj, base, E, lamv, variant)
    P(f"  {label}: Ric(g_aleph) = rho g_aleph on s(X): {ok}  (rho = s/14 = {rho}; kappa = {kappa}; "
      f"#nonzero residual components {nz})  [{time.time() - t:.1f}s]")
# (iii) non-Ricci-flat random metric
quad2 = {(a, b): rsym() for a in range(n) for b in range(a, n)}
g_rand = poly_metric(n, E, quad2)
b_rand = base_geometry(g_rand, n, ctx)
ok, rho, nz, Ric, G0, idx, Psi0 = einstein_test(g_rand, b_rand, E, ctx.c(F(1, 5)), 'Z1')
RF = frame_blocks(Ric, n, idx, Psi0, ctx)
hh_res = [[RF[a][b] - kappa * cH * E[a][b] for b in range(n)] for a in range(n)]
P(f"  (iii) random metric (Ric(g)(0) != 0: {any(b_rand['Ric'][a][b] != 0 for a in range(n) for b in range(n))}), Z1, "
  f"lambda = 1/5: Einstein on s(X): {ok};  Ric_HH - kappa cH g == Ric(g): "
  f"{all(hh_res[a][b] == b_rand['Ric'][a][b] for a in range(n) for b in range(n))}")
# (v) off the section
while True:
    y0 = [[E[i][j] + (ctx.c(F(random.randint(-2, 2), 5)) if i <= j else 0) for j in range(n)] for i in range(n)]
    for i in range(n):
        for j in range(i):
            y0[i][j] = y0[j][i]
    try:
        finv(y0, ctx)
        break
    except StopIteration:
        pass
t = time.time()
ok, rho, nz, Ric, G0, idx, Psi0 = einstein_test(g_weyl, b_weyl, y0, ctx.c(F(1, 5)), 'Z1')
P(f"  (v)   Weyl metric, Z1, lambda = 1/5, at y0 != g(x0): Ric(g_aleph) proportional to g_aleph: {ok} "
  f"(#nonzero residual components {nz} of {len(G0) ** 2})  [{time.time() - t:.1f}s]")
gflat = poly_metric(n, E, {})
bflat = base_geometry(gflat, n, ctx)
ok, rho, nz, *_ = einstein_test(gflat, bflat, y0, ctx.c(F(1, 5)), 'Z1')
P(f"  control: flat g, Z1, lambda = 1/5, at the same y0: Einstein: {ok}, rho = {rho}\n")

# ------------------------------------------------------------------ E4
import sympy as sp
P("E4  Kaluza-Klein isometry: Z1(flat aleph, lambda = 1/(n+1), cH = 2 cV)  ~  SL(n+1,R)/SO(p',q')")
tt = sp.Symbol('t')
for nn in (2, 3, 4):
    reps = 1 if nn == 2 else 10
    ok = True
    for r in range(reps):
        if nn == 2:
            us = sp.symbols('u0:3'); vs = sp.symbols('v0:3'); a0 = sp.symbols('a0:2'); a1 = sp.symbols('b0:2')
        else:
            rnd = lambda: sp.Rational(random.randint(-5, 5), random.randint(1, 3))
            us = [rnd() for _ in range(nn * (nn + 1) // 2)]; vs = [rnd() for _ in range(nn * (nn + 1) // 2)]
            a0 = [rnd() for _ in range(nn)]; a1 = [rnd() for _ in range(nn)]
        U = sp.zeros(nn, nn); V = sp.zeros(nn, nn)
        for A, (i, j) in enumerate(sym_index(nn)):
            U[i, j] = U[j, i] = us[A]
            V[i, j] = V[j, i] = vs[A]
        if U.det() == 0:
            continue
        u_t = U + tt * V
        a_t = sp.Matrix(a0) + tt * sp.Matrix(a1)
        L = sp.eye(nn + 1)
        for i in range(nn):
            L[i, nn] = a_t[i]
        Dm = sp.zeros(nn + 1, nn + 1)
        Dm[:nn, :nn] = u_t
        Dm[nn, nn] = 1 / u_t.det()
        M = L.T * Dm * L
        dM = M.diff(tt).subs(tt, 0)
        M0 = M.subs(tt, 0)
        X = M0.inv() * dM
        lhs = (X * X).trace()
        y_t = u_t.det() * u_t
        y0 = y_t.subs(tt, 0)
        dy = y_t.diff(tt).subs(tt, 0)
        yi = y0.inv()
        lamv = sp.Rational(1, nn + 1)
        rhs = (yi * dy * yi * dy).trace() - lamv * (yi * dy).trace() ** 2 + 2 * (sp.Matrix(a1).T * y0 * sp.Matrix(a1))[0, 0]
        ok &= sp.cancel(sp.together(lhs - rhs)) == 0
        ok &= sp.simplify(M0.det() - 1) == 0
    P(f"  n={nn}: tr((M^-1 dM)^2) == G_(1/(n+1))(dy,dy) + 2 y(dx,dx) with det M = 1: {ok}")
P("  dim SL(n+1)/SO = n(n+3)/2 = dim Y; Einstein constant of tr((M^-1 dM)^2) on SL(n+1)/SO is -(n+1)/4 (Prop. F1 at n+1)")
P("")

# ------------------------------------------------------------------ E3
P("E3  fiber scale for the observed cosmological constant (floating point)")
Lobs = 1.1056e-52      # m^-2
lF = math.sqrt(9 / (2 * Lobs))
P(f"  Z2, n = 4: |Lambda_F| = 9/(2 l_F^2) = Lambda_obs  =>  l_F = {lF:.3e} m = {lF / 9.4607e15 / 1e9:.2f} Gly")
P(f"  (Planck length 1.616e-35 m would give |Lambda_F| = {9 / (2 * 1.616e-35 ** 2):.3e} m^-2, "
  f"ratio to Lambda_obs {9 / (2 * 1.616e-35 ** 2) / Lobs:.2e})")
P(f"\nTotal time {time.time() - T0:.1f}s")
