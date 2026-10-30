"""GU round 2 (F): geometry of the fiber of Y = Met(X^n) with the DeWitt metric c_V G_lambda.
Write-up: formalization/geometric-unity/F-zorro-reduction.md (Prop. F1). Gaps GU-G05, GU-G06.

  G_lambda|_y(h,k) = tr(y^-1 h y^-1 k) - lambda tr(y^-1 h) tr(y^-1 k)    on Sym^2(R^n), y of signature (p,q).

Checks (all exact; lambda, c_V symbolic in Q(lambda, cH, cV)):
  FC1  Ricci at y0:  Ric = -(n/4) G_{1/n}  (traceless projection), independent of lambda and c_V;
       scalar s_F = -n(n-1)(n+2) / (8 c_V);  for n = 2,3,4, every signature p >= q, plus non-diagonal y0.
  FC2  full curvature (n = 2, 3):  c_V G_lambda(R(h_c,h_d)h_b, h_a) = -(c_V/4) tr([[H_c,H_d],H_b] H_a),
       H = y^-1 h  (convention R_abcd = g(R(d_c,d_d)d_b, d_a)); lambda-independent.
  FC3  geodesic symmetry S_y0(z) = y0 z^-1 y0 is an isometry of G_lambda (symbolic for n = 2; exact at random
       rational points for n = 3, 4; lambda symbolic), fixes y0, dS_y0 = -Id
       (=> (Met(R^n), G_lambda) is a pseudo-Riemannian symmetric space for every lambda != 1/n).
  FC4  product splitting y = e^s u, |det u| = 1:  G_lambda = n(1 - n lambda) ds^2 + tr((u^-1 du)^2)  (n = 2, 3, 4).
"""
import itertools
import random
import time
from fractions import Fraction
import sympy as sp
from F_jets import Ctx, jconst, jinv, dewitt_matrix, ricci_at_point, riemann_at_point, trace, sym_index, supp, \
    fmatmul, finv

random.seed(20210401)
ctx = Ctx('lam,cH,cV')
lam, cH, cV = ctx.gens
T0 = time.time()


def fiber_metric_jets(n, y0, cVv):
    idx = sym_index(n)
    pos = {}
    for A, (i, j) in enumerate(idx):
        pos[(i, j)] = A
        pos[(j, i)] = A
    Y = [[None] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            jt = jconst(y0[i][j], 2)
            jt[1][(pos[(i, j)],)] = ctx.one
            Y[i][j] = jt
    W = jinv(Y, 2, ctx)
    Gm = dewitt_matrix(W, n, lam, 2, ctx)
    Gm = [[[{k: v * cVv for k, v in d.items()} for d in Gm[A][B]] for B in range(len(idx))] for A in range(len(idx))]
    return Gm, idx


def basis_mats(n, idx):
    mats = []
    for (i, j) in idx:
        M = [[ctx.zero] * n for _ in range(n)]
        for (a, b) in supp((i, j)):
            M[a][b] = ctx.one
        mats.append(M)
    return mats


def tr(M):
    return sum((M[i][i] for i in range(len(M))), ctx.zero)


def dewitt_value(y0i, h, k, lv):
    return tr(fmatmul(fmatmul(y0i, h, ctx), fmatmul(y0i, k, ctx), ctx)) - lv * tr(fmatmul(y0i, h, ctx)) * tr(fmatmul(y0i, k, ctx))


def diag(eps):
    n = len(eps)
    return [[ctx.c(eps[i]) if i == j else ctx.zero for j in range(n)] for i in range(n)]


def random_nondiag(eps):
    n = len(eps)
    while True:
        P = [[ctx.c(Fraction(random.randint(-3, 3), random.randint(1, 2))) for _ in range(n)] for _ in range(n)]
        try:
            finv(P, ctx)
            break
        except StopIteration:
            continue
    PT = [list(r) for r in zip(*P)]
    return fmatmul(fmatmul(PT, diag(eps), ctx), P, ctx)


print("F_fiber_curvature.py: fiber (Met(R^n), c_V G_lambda); exact over Q(lambda, cH, cV)\n")
print("FC1  Ricci and scalar curvature at y0")
allok = True
cases = []
for n in (2, 3, 4):
    for q in range(0, n // 2 + 1):
        p = n - q
        cases.append((n, (p, q), [1] * p + [-1] * q, 'diag'))
cases += [(3, (2, 1), [1, 1, -1], 'nondiag'), (4, (1, 3), [1, -1, -1, -1], 'nondiag')]
for n, sig, eps, kind in cases:
    y0 = diag(eps) if kind == 'diag' else random_nondiag(eps)
    Gm, idx = fiber_metric_jets(n, y0, cV)
    N = len(idx)
    Ric, g0, g0i = ricci_at_point(Gm, N, ctx)
    s = trace(Ric, g0i, ctx)
    y0i = finv(y0, ctx)
    E = basis_mats(n, idx)
    pred = [[-Fraction(n, 4) * dewitt_value(y0i, E[A], E[B], ctx.c(Fraction(1, n))) for B in range(N)] for A in range(N)]
    ok_ric = all(Ric[A][B] == pred[A][B] for A in range(N) for B in range(N))
    s_pred = ctx.c(Fraction(-n * (n - 1) * (n + 2), 8)) / cV
    ok_s = (s == s_pred)
    allok &= ok_ric and ok_s
    print(f"  n={n} sig={sig} y0={kind:7s}: Ric == -(n/4) G_(1/n): {ok_ric};  s_F = {s}  (pred {s_pred}): {ok_s}")
print(f"  FC1 all cases pass: {allok}   [{time.time() - T0:.1f}s]\n")

print("FC2  full curvature tensor vs  -(c_V/4) tr([[H_c,H_d],H_b] H_a)")
for n, eps, kind in [(2, [1, 1], 'diag'), (2, [1, -1], 'diag'), (3, [1, 1, 1], 'diag'), (3, [1, 1, -1], 'nondiag')]:
    y0 = diag(eps) if kind == 'diag' else random_nondiag(eps)
    Gm, idx = fiber_metric_jets(n, y0, cV)
    N = len(idx)
    R, g0, g0i = riemann_at_point(Gm, N, ctx)
    y0i = finv(y0, ctx)
    H = [fmatmul(y0i, M, ctx) for M in basis_mats(n, idx)]

    def comm(a, b):
        ab = fmatmul(a, b, ctx)
        ba = fmatmul(b, a, ctx)
        return [[ab[i][j] - ba[i][j] for j in range(n)] for i in range(n)]
    ok = True
    for a, b, c, d in itertools.product(range(N), repeat=4):
        pred = -cV / 4 * tr(fmatmul(comm(comm(H[c], H[d]), H[b]), H[a], ctx))
        if R[(a, b, c, d)] != pred:
            ok = False
    print(f"  n={n} eps={eps} y0={kind}: all {N**4} components match the lambda-free formula: {ok}")
print(f"  [{time.time() - T0:.1f}s]\n")

print("FC3  geodesic symmetry S(z) = y0 z^-1 y0: isometry of G_lambda, fixes y0, dS_y0 = -Id")
print("     (n = 2: symbolic in z, k, lambda; n = 3, 4: exact at 20 random rational (z, k) each, lambda symbolic)")
L = sp.Symbol('lam')


def G_sp(y, h, k):
    yi = y.inv()
    return (yi * h * yi * k).trace() - L * (yi * h).trace() * (yi * k).trace()


def sym_from(vals, n):
    M = sp.zeros(n, n)
    for A, (i, j) in enumerate(sym_index(n)):
        M[i, j] = M[j, i] = vals[A]
    return M


for n, eps in [(2, [1, -1]), (3, [1, 1, -1]), (4, [1, -1, -1, -1])]:
    Y0 = sp.diag(*eps)
    if n == 2:
        Z = sym_from(sp.symbols('z0:3'), n)
        Kmat = sym_from(sp.symbols('k0:3'), n)
        samples = [(Z, Kmat)]
    else:
        samples = []
        while len(samples) < 20:
            Z = sym_from([sp.Rational(random.randint(-5, 5), random.randint(1, 4)) for _ in range(n * (n + 1) // 2)], n)
            if Z.det() == 0:
                continue
            Kmat = sym_from([sp.Rational(random.randint(-5, 5), random.randint(1, 4)) for _ in range(n * (n + 1) // 2)], n)
            samples.append((Z, Kmat))
    iso = True
    for Z, Kmat in samples:
        Zi = Z.inv()
        Sz = Y0 * Zi * Y0
        dS = -Y0 * Zi * Kmat * Zi * Y0          # derivative of z -> y0 z^-1 y0 in direction k
        iso &= sp.cancel(sp.together(G_sp(Sz, dS, dS) - G_sp(Z, Kmat, Kmat))) == 0
    Kmat = samples[0][1]
    fix = (Y0 * Y0.inv() * Y0 - Y0) == sp.zeros(n, n)
    dS0 = (-Y0 * Y0.inv() * Kmat * Y0.inv() * Y0 + Kmat) == sp.zeros(n, n)
    print(f"  n={n} y0=diag{tuple(eps)}: S^*G_lambda = G_lambda: {iso};  S(y0) = y0: {fix};  dS_y0 = -Id: {dS0}")
print()

print("FC4  splitting y = e^s u (|det u| = 1): y^-1 dy = ds I + u^-1 du with tr(u^-1 du) = 0, so")
print("     G_lambda = n(1 - n lambda) ds^2 + tr((u^-1 du)^2)   (n = 2 symbolic; n = 3, 4 at 20 random rational points)")
ds_ = sp.Symbol('ds')
for n in (2, 3, 4):
    ok = True
    reps = 1 if n == 2 else 20
    for r in range(reps):
        if n == 2:
            U = sym_from(sp.symbols('m0:3'), n)
            dU = sym_from(sp.symbols('dm0:3'), n)
            dms = sp.symbols('dm0:3')
        else:
            while True:
                U = sym_from([sp.Rational(random.randint(-5, 5), random.randint(1, 4)) for _ in range(n * (n + 1) // 2)], n)
                if U.det() != 0:
                    break
            dms = sp.symbols(f'dm0:{n * (n + 1) // 2}')
            dU = sym_from(list(dms[:-1]) + [dms[-1]], n)
            vals = {dms[i]: sp.Rational(random.randint(-5, 5), random.randint(1, 4)) for i in range(len(dms) - 1)}
            dU = dU.subs(vals)
        cons = sp.solve(sp.Eq((U.inv() * dU).trace(), 0), dms[-1])
        if not cons:
            continue
        dU = dU.subs(dms[-1], cons[0])
        X = ds_ * sp.eye(n) + U.inv() * dU       # = y^-1 dy for y = e^s U
        Gl = (X * X).trace() - L * X.trace() ** 2
        pred = n * (1 - n * L) * ds_ ** 2 + (U.inv() * dU * U.inv() * dU).trace()
        ok &= sp.cancel(sp.together(sp.expand(Gl - pred))) == 0
    print(f"  n={n}: identity holds: {ok}")
print(f"\nTotal time {time.time() - T0:.1f}s")
