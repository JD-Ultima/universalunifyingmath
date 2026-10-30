"""Exact jet calculus for curvature at a point (library for F_fiber_curvature.py, F_zorro_reduction.py).

A jet of order D at the origin of R^m is the Taylor polynomial truncated at total degree D, stored as a
list of dicts by degree: jet[k] = {monomial: coeff}, monomial = sorted tuple of variable indices of length k.
Coefficients live in an exact field (sympy rational-function field, e.g. Q(lambda, cH, cV)).

The curvature of a metric at the origin needs only its 2-jet; formulas used (any coordinates):
  Gamma_{e,db} = 1/2 (d_d g_eb + d_b g_ed - d_e g_db),  Gamma^a_{db} = g^{ae} Gamma_{e,db}
  R^a_{bcd} = d_c Gamma^a_{db} - d_d Gamma^a_{cb} + Gamma^a_{ce} Gamma^e_{db} - Gamma^a_{de} Gamma^e_{cb}
  (R(d_c, d_d) d_b = R^a_{bcd} d_a,  R_{abcd} = g_ae R^e_{bcd} = g(R(d_c,d_d)d_b, d_a),  Ric_bd = R^a_{bad}).
With this convention the unit sphere has R_{1212} > 0 and Ric = (n-1) g (control test in the scripts).
No floating point anywhere.
"""
import itertools
from fractions import Fraction
from sympy import QQ, Rational
from sympy.polys.fields import field


class Ctx:
    """Coefficient field. names: generator names, e.g. 'lam,cH,cV'."""
    def __init__(self, names):
        out = field(names, QQ)
        self.K = out[0]
        self.gens = out[1:]
        self.zero = self.K.zero
        self.one = self.K.one

    def c(self, x):
        if isinstance(x, Fraction):
            return self.K(Rational(x.numerator, x.denominator))
        if isinstance(x, int):
            return self.K(x)
        return x  # already a field element


# ------------------------------------------------------------------ jets
def jzero(D):
    return [dict() for _ in range(D + 1)]


def jconst(c, D):
    j = jzero(D)
    j[0][()] = c
    return j


def jvar(i, D, ctx):
    j = jzero(D)
    if D >= 1:
        j[1][(i,)] = ctx.one
    return j


def jadd(a, b):
    D = min(len(a), len(b)) - 1
    out = jzero(D)
    for k in range(D + 1):
        o = dict(a[k])
        for m, c in b[k].items():
            o[m] = o[m] + c if m in o else c
        out[k] = o
    return out


def jscale(a, s):
    return [{m: c * s for m, c in d.items()} for d in a]


def jsub(a, b):
    return jadd(a, jscale(b, -1))


def jmul(a, b, D=None):
    if D is None:
        D = min(len(a), len(b)) - 1
    out = jzero(D)
    for i in range(min(D, len(a) - 1) + 1):
        if not a[i]:
            continue
        for j in range(min(D - i, len(b) - 1) + 1):
            if not b[j]:
                continue
            o = out[i + j]
            for ma, ca in a[i].items():
                for mb, cb in b[j].items():
                    m = tuple(sorted(ma + mb)) if (ma and mb) else (ma or mb)
                    v = ca * cb
                    o[m] = o[m] + v if m in o else v
    return out


def jtrunc(a, D):
    return [dict(a[k]) if k < len(a) else dict() for k in range(D + 1)]


def jderiv(a, i):
    """d/dz_i ; result has order len(a)-2."""
    D = len(a) - 2
    out = jzero(max(D, 0))
    for k in range(1, len(a)):
        for m, c in a[k].items():
            cnt = m.count(i)
            if cnt:
                l = list(m)
                l.remove(i)
                mm = tuple(l)
                v = c * cnt
                out[k - 1][mm] = out[k - 1][mm] + v if mm in out[k - 1] else v
    return out


def jval(a, ctx):
    return a[0].get((), ctx.zero)


def jgrad(a, i, ctx):
    return a[1].get((i,), ctx.zero) if len(a) > 1 else ctx.zero


def jhess(a, i, j, ctx):
    if len(a) < 3:
        return ctx.zero
    m = (i, j) if i <= j else (j, i)
    c = a[2].get(m, ctx.zero)
    return c * 2 if i == j else c


def jiszero(a):
    return all(c == 0 for d in a for c in d.values())


# ------------------------------------------------------------------ matrices of jets / of field elements
def mmul(A, B, D):
    n, k, p = len(A), len(B), len(B[0])
    out = [[None] * p for _ in range(n)]
    for i in range(n):
        for j in range(p):
            acc = jzero(D)
            for l in range(k):
                acc = jadd(acc, jmul(A[i][l], B[l][j], D))
            out[i][j] = acc
    return out


def finv(M, ctx):
    """Exact inverse of a square matrix over the field (Gauss-Jordan)."""
    n = len(M)
    A = [list(row) + [ctx.one if i == j else ctx.zero for j in range(n)] for i, row in enumerate(M)]
    for col in range(n):
        piv = next(r for r in range(col, n) if A[r][col] != 0)
        A[col], A[piv] = A[piv], A[col]
        p = A[col][col]
        A[col] = [x / p for x in A[col]]
        for r in range(n):
            if r != col and A[r][col] != 0:
                f = A[r][col]
                A[r] = [x - f * y for x, y in zip(A[r], A[col])]
    return [row[n:] for row in A]


def fmatmul(A, B, ctx):
    return [[sum((A[i][l] * B[l][j] for l in range(len(B))), ctx.zero) for j in range(len(B[0]))]
            for i in range(len(A))]


def jinv(M, D, ctx):
    """Inverse of a matrix of jets with invertible constant part: Neumann series to order D."""
    n = len(M)
    M0 = [[jval(M[i][j], ctx) for j in range(n)] for i in range(n)]
    M0i = finv(M0, ctx)
    M0iJ = [[jconst(M0i[i][j], D) for j in range(n)] for i in range(n)]
    E = [[jsub(M[i][j], jconst(M0[i][j], D)) for j in range(n)] for i in range(n)]
    X = mmul(M0iJ, E, D)                       # M0^{-1} E
    X = [[jscale(x, -1) for x in row] for row in X]
    term = M0iJ
    acc = M0iJ
    for _ in range(D):
        term = mmul(X, term, D)
        acc = [[jadd(acc[i][j], term[i][j]) for j in range(n)] for i in range(n)]
    return acc


# ------------------------------------------------------------------ curvature of a metric at the origin
def metric_point_data(G, m, ctx):
    """G: m x m matrix of jets (order >= 2). Returns (g0, dG[c][a][b], ddG[c][d][a][b])."""
    g0 = [[jval(G[a][b], ctx) for b in range(m)] for a in range(m)]
    dG = [[[jgrad(G[a][b], c, ctx) for b in range(m)] for a in range(m)] for c in range(m)]
    ddG = [[[[jhess(G[a][b], c, d, ctx) for b in range(m)] for a in range(m)] for d in range(m)]
           for c in range(m)]
    return g0, dG, ddG


def christoffel(g0i, dG, m, ctx):
    G1 = [[[(dG[d][e][b] + dG[b][e][d] - dG[e][d][b]) / 2 for b in range(m)] for d in range(m)]
          for e in range(m)]
    Gam = [[[sum((g0i[a][e] * G1[e][d][b] for e in range(m) if g0i[a][e] != 0), ctx.zero)
             for b in range(m)] for d in range(m)] for a in range(m)]
    return G1, Gam


def ricci_at_point(G, m, ctx):
    """Ricci tensor (coordinate components) and inverse metric at the origin, from the 2-jet of G."""
    g0, dG, ddG = metric_point_data(G, m, ctx)
    g0i = finv(g0, ctx)
    G1, Gam = christoffel(g0i, dG, m, ctx)
    # M_c = g^{-1} d_c g
    Mc = [fmatmul(g0i, dG[c], ctx) for c in range(m)]
    # sum_a d_a Gamma^a_{db} = sum_a [ d_a g^{ae} Gamma_{e,db} + g^{ae} d_a Gamma_{e,db} ]
    # d_a g^{ae} = -(g^{-1} d_a g g^{-1})^{ae}
    dginv_diag = [[-sum((Mc[a][a][p] * g0i[p][e] for p in range(m) if Mc[a][a][p] != 0), ctx.zero)
                   for e in range(m)] for a in range(m)]          # [a][e] = d_a g^{ae}
    Ric = [[ctx.zero] * m for _ in range(m)]
    for b in range(m):
        for d in range(b, m):
            t1 = ctx.zero
            for a in range(m):
                for e in range(m):
                    if dginv_diag[a][e] != 0:
                        t1 += dginv_diag[a][e] * G1[e][d][b]
                    if g0i[a][e] != 0:
                        dG1 = (ddG[a][d][e][b] + ddG[a][b][e][d] - ddG[a][e][d][b]) / 2
                        t1 += g0i[a][e] * dG1
            # d_d Gamma^a_{ab} = 1/2 d_d (g^{ae} d_b g_ae) = 1/2 ( -tr(M_d M_b) + g^{ae} dd_{db} g_ae )
            t2 = -sum((Mc[d][p][q] * Mc[b][q][p] for p in range(m) for q in range(m)), ctx.zero)
            t2 += sum((g0i[a][e] * ddG[d][b][a][e] for a in range(m) for e in range(m) if g0i[a][e] != 0),
                      ctx.zero)
            t2 = t2 / 2
            t3 = ctx.zero
            for a in range(m):
                for e in range(m):
                    t3 += Gam[a][a][e] * Gam[e][d][b] - Gam[a][d][e] * Gam[e][a][b]
            Ric[b][d] = t1 - t2 + t3
            Ric[d][b] = Ric[b][d]
    return Ric, g0, g0i


def riemann_at_point(G, m, ctx):
    """All components R_{abcd} = g(R(d_c,d_d)d_b, d_a) at the origin (use for small m)."""
    g0, dG, ddG = metric_point_data(G, m, ctx)
    g0i = finv(g0, ctx)
    G1, Gam = christoffel(g0i, dG, m, ctx)
    # d_c Gamma^a_{db}
    dginv = [[[-sum((g0i[a][p] * dG[c][p][q] * g0i[q][e] for p in range(m) for q in range(m)), ctx.zero)
               for e in range(m)] for a in range(m)] for c in range(m)]

    def dGam(c, a, d, b):
        s = ctx.zero
        for e in range(m):
            s += dginv[c][a][e] * G1[e][d][b]
            if g0i[a][e] != 0:
                s += g0i[a][e] * (ddG[c][d][e][b] + ddG[c][b][e][d] - ddG[c][e][d][b]) / 2
        return s

    Rup = {}
    for a, b, c, d in itertools.product(range(m), repeat=4):
        if c >= d:
            continue
        v = dGam(c, a, d, b) - dGam(d, a, c, b)
        for e in range(m):
            v += Gam[a][c][e] * Gam[e][d][b] - Gam[a][d][e] * Gam[e][c][b]
        Rup[(a, b, c, d)] = v
        Rup[(a, b, d, c)] = -v
    for a, b, c in itertools.product(range(m), repeat=3):
        Rup[(a, b, c, c)] = ctx.zero
    R = {}
    for a, b, c, d in itertools.product(range(m), repeat=4):
        R[(a, b, c, d)] = sum((g0[a][e] * Rup[(e, b, c, d)] for e in range(m) if g0[a][e] != 0), ctx.zero)
    return R, g0, g0i


def ricci_from_riemann(R, g0i, m, ctx):
    return [[sum((g0i[a][c] * R[(c, b, a, d)] for a in range(m) for c in range(m) if g0i[a][c] != 0),
                 ctx.zero) for d in range(m)] for b in range(m)]


def trace(M, g0i, ctx):
    m = len(M)
    return sum((g0i[a][b] * M[a][b] for a in range(m) for b in range(m) if g0i[a][b] != 0), ctx.zero)


# ------------------------------------------------------------------ Sym^2 coordinates and DeWitt metric
def sym_index(n):
    """Coordinates y_A, A = (i,j) with i <= j. Basis tangent vectors E^A = E_ii or E_ij + E_ji."""
    return [(i, j) for i in range(n) for j in range(i, n)]


def supp(A):
    i, j = A
    return [(i, j)] if i == j else [(i, j), (j, i)]


def dewitt_matrix(W, n, lam, D, ctx):
    """Jets of G_lambda(E^A, E^B) at y, given W = jets of y^{-1}:
       tr(W E^A W E^B) - lam tr(W E^A) tr(W E^B)."""
    idx = sym_index(n)
    trWE = [None] * len(idx)
    for A, (i, j) in enumerate(idx):
        acc = jzero(D)
        for (a, b) in supp((i, j)):
            acc = jadd(acc, W[b][a])
        trWE[A] = acc
    Gm = [[None] * len(idx) for _ in idx]
    for A, IA in enumerate(idx):
        for B, IB in enumerate(idx):
            if B < A:
                Gm[A][B] = Gm[B][A]
                continue
            acc = jzero(D)
            for (a, b) in supp(IA):
                for (c, p) in supp(IB):
                    acc = jadd(acc, jmul(W[p][a], W[b][c], D))
            acc = jsub(acc, jscale(jmul(trWE[A], trWE[B], D), lam))
            Gm[A][B] = acc
    return Gm


# ------------------------------------------------------------------ base metric (polynomial) and its curvature
def base_from_coeffs(g0, lin, quad, cub, n, ctx, D=3):
    """g(x) = g0 + sum lin[mu] x^mu + sum quad[(mu,nu)] x^mu x^nu + sum cub[(mu,nu,rho)] x x x
       (coefficient matrices symmetric n x n with field entries). Returns n x n matrix of order-D jets."""
    G = [[jzero(D) for _ in range(n)] for _ in range(n)]
    for i in range(n):
        for j in range(n):
            G[i][j][0][()] = g0[i][j]
            for mu in range(n):
                if lin[mu][i][j] != 0 and D >= 1:
                    G[i][j][1][(mu,)] = lin[mu][i][j]
            for key, M in quad.items():
                if M[i][j] != 0 and D >= 2:
                    G[i][j][2][key] = G[i][j][2].get(key, ctx.zero) + M[i][j]
            for key, M in cub.items():
                if M[i][j] != 0 and D >= 3:
                    G[i][j][3][key] = G[i][j][3].get(key, ctx.zero) + M[i][j]
    return G


def base_geometry(g, n, ctx):
    """From a 3-jet of the base metric: Christoffel jets (order 2), and at the origin:
       Riemann R_{abcd}, Ric, scalar R, nabla Ric, d R, g0, g0^{-1}."""
    D = 3
    ginv = jinv(g, D, ctx)
    dg = [[[jderiv(g[a][b], c) for b in range(n)] for a in range(n)] for c in range(n)]  # order 2
    half = ctx.c(Fraction(1, 2))
    G1 = [[[jscale(jsub(jadd(dg[d][e][b], dg[b][e][d]), dg[e][d][b]), half)
            for b in range(n)] for d in range(n)] for e in range(n)]
    Gam = [[[None] * n for _ in range(n)] for _ in range(n)]
    for a in range(n):
        for d in range(n):
            for b in range(n):
                acc = jzero(2)
                for e in range(n):
                    acc = jadd(acc, jmul(jtrunc(ginv[a][e], 2), G1[e][d][b], 2))
                Gam[a][d][b] = acc
    # Riemann as order-1 jets: R^a_{bcd}
    Rup = {}
    for a, b, c, d in itertools.product(range(n), repeat=4):
        v = jsub(jderiv(Gam[a][d][b], c), jderiv(Gam[a][c][b], d))
        for e in range(n):
            v = jadd(v, jsub(jmul(Gam[a][c][e], Gam[e][d][b], 1), jmul(Gam[a][d][e], Gam[e][c][b], 1)))
        Rup[(a, b, c, d)] = jtrunc(v, 1)
    Ricj = [[None] * n for _ in range(n)]
    for b in range(n):
        for d in range(n):
            acc = jzero(1)
            for a in range(n):
                acc = jadd(acc, Rup[(a, b, a, d)])
            Ricj[b][d] = acc
    Rsj = jzero(1)
    for a in range(n):
        for b in range(n):
            Rsj = jadd(Rsj, jmul(jtrunc(ginv[a][b], 1), Ricj[a][b], 1))
    g0 = [[jval(g[a][b], ctx) for b in range(n)] for a in range(n)]
    g0i = [[jval(ginv[a][b], ctx) for b in range(n)] for a in range(n)]
    Gam0 = [[[jval(Gam[a][d][b], ctx) for b in range(n)] for d in range(n)] for a in range(n)]
    R = {}
    for a, b, c, d in itertools.product(range(n), repeat=4):
        R[(a, b, c, d)] = sum((g0[a][e] * jval(Rup[(e, b, c, d)], ctx) for e in range(n)), ctx.zero)
    Ric = [[jval(Ricj[b][d], ctx) for d in range(n)] for b in range(n)]
    Rs = jval(Rsj, ctx)
    dRic = [[[jgrad(Ricj[b][d], r, ctx) for d in range(n)] for b in range(n)] for r in range(n)]
    nablaRic = [[[dRic[r][b][d]
                  - sum((Gam0[e][r][b] * Ric[e][d] + Gam0[e][r][d] * Ric[b][e] for e in range(n)), ctx.zero)
                  for d in range(n)] for b in range(n)] for r in range(n)]
    dR = [jgrad(Rsj, r, ctx) for r in range(n)]
    return dict(Gam=Gam, R=R, Ric=Ric, Rs=Rs, nablaRic=nablaRic, dR=dR, g0=g0, g0i=g0i, Gam0=Gam0)


# ------------------------------------------------------------------ total space Y = Met(X): the Zorro metrics
def total_metric(g, base, n, y0, lam, cH, cV, variant, ctx):
    """2-jet at z0 = (x=0, y=y0) of the metric on Y in coordinates z = (x^mu, y_A):
         g_aleph = cH * h_H(dx,dx) + cV * G_lambda|_y(Dy, Dy),
         Dy_ij = dy_ij - (Gamma^k_{mu i} y_kj + Gamma^k_{mu j} y_ik) dx^mu   (Levi-Civita of g, induced on Sym^2 T*)
       variant 'Z1': h_H = y (tautological: the point of Y itself; draft (3.7))
       variant 'Z2': h_H = g(x) (base metric pulled back by pi; Riemannian submersion)
       Returns (G, m, idx, Psi0) with Psi0[A][mu] the value at z0 of the connection coefficients."""
    D = 2
    idx = sym_index(n)
    N = len(idx)
    m = n + N
    pos = {}
    for A, (i, j) in enumerate(idx):
        pos[(i, j)] = A
        pos[(j, i)] = A
    Y = [[None] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            jt = jconst(y0[i][j], D)
            jt[1][(n + pos[(i, j)],)] = ctx.one
            Y[i][j] = jt
    W = jinv(Y, D, ctx)
    Gm = dewitt_matrix(W, n, lam, D, ctx)
    Gam = [[[jtrunc(base['Gam'][a][d][b], D) for b in range(n)] for d in range(n)] for a in range(n)]
    Psi = [[None] * n for _ in range(N)]
    for A, (i, j) in enumerate(idx):
        for mu in range(n):
            acc = jzero(D)
            for k in range(n):
                acc = jadd(acc, jmul(Gam[k][mu][i], Y[k][j], D))
                acc = jadd(acc, jmul(Gam[k][mu][j], Y[i][k], D))
            Psi[A][mu] = acc
    Q = [[None] * n for _ in range(N)]          # Q[B][mu] = sum_A G_AB Psi_A mu
    for B in range(N):
        for mu in range(n):
            acc = jzero(D)
            for A in range(N):
                acc = jadd(acc, jmul(Gm[A][B], Psi[A][mu], D))
            Q[B][mu] = acc
    G = [[None] * m for _ in range(m)]
    for mu in range(n):
        for nu in range(mu, n):
            h = Y[mu][nu] if variant == 'Z1' else jtrunc(g[mu][nu], D)
            acc = jscale(h, cH)
            s = jzero(D)
            for B in range(N):
                s = jadd(s, jmul(Q[B][mu], Psi[B][nu], D))
            acc = jadd(acc, jscale(s, cV))
            G[mu][nu] = acc
            G[nu][mu] = acc
    for mu in range(n):
        for B in range(N):
            v = jscale(Q[B][mu], -cV)
            G[mu][n + B] = v
            G[n + B][mu] = v
    for A in range(N):
        for B in range(N):
            G[n + A][n + B] = jscale(Gm[A][B], cV)
    Psi0 = [[jval(Psi[A][mu], ctx) for mu in range(n)] for A in range(N)]
    return G, m, idx, Psi0


def frame_blocks(T, n, idx, Psi0, ctx):
    """Components of a symmetric 2-tensor T (coordinate components on Y at z0) in the adapted frame
       e_mu = d_mu + Psi0[A][mu] d_A (horizontal lifts), e_A = d_A (vertical)."""
    N = len(idx)
    m = n + N
    E = [[ctx.zero] * m for _ in range(m)]   # E[frame][coord]
    for mu in range(n):
        E[mu][mu] = ctx.one
        for A in range(N):
            E[mu][n + A] = Psi0[A][mu]
    for A in range(N):
        E[n + A][n + A] = ctx.one
    out = [[sum((E[p][a] * T[a][b] * E[q][b] for a in range(m) for b in range(m)
                 if E[p][a] != 0 and E[q][b] != 0), ctx.zero) for q in range(m)] for p in range(m)]
    return out


def solve_fit(columns, target, ctx):
    """Exact solution x of sum_k x_k columns[k] = target (vectors over the field), with consistency check.
       Returns (x, residual_is_zero)."""
    K = len(columns)
    L = len(target)
    A = [[columns[k][r] for k in range(K)] + [target[r]] for r in range(L)]
    piv_cols = []
    row = 0
    for col in range(K):
        pr = next((r for r in range(row, L) if A[r][col] != 0), None)
        if pr is None:
            continue
        A[row], A[pr] = A[pr], A[row]
        p = A[row][col]
        A[row] = [x / p for x in A[row]]
        for r in range(L):
            if r != row and A[r][col] != 0:
                f = A[r][col]
                A[r] = [x - f * y for x, y in zip(A[r], A[row])]
        piv_cols.append(col)
        row += 1
    consistent = all(A[r][K] == 0 for r in range(row, L))
    x = [ctx.zero] * K
    for r, col in enumerate(piv_cols):
        x[col] = A[r][K]
    return x, consistent, len(piv_cols)
