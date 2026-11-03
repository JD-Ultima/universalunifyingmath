"""Exact controls for AU vertical Ricci geometry and critical values.

Uses the accepted F_jets metric library, without running any older suite.
The universal geometric derivation is in AU-opposite-ricci-spectrum.md.
"""
from pathlib import Path
from time import perf_counter
import sys
import sympy as sp

sys.path.insert(0, str(Path(__file__).resolve().parent))
import F_jets as F


CHECKS = []


def check(label, condition):
    good = bool(condition)
    CHECKS.append(good)
    print(f'[{len(CHECKS):02d}] {label}: {good}', flush=True)
    assert good, label


def metric_control(n, lam):
    started = perf_counter()
    ctx = F.Ctx('z')
    z = ctx.zero

    def matrix():
        return [[z for _ in range(n)] for _ in range(n)]

    g0 = matrix()
    for i in range(n):
        g0[i][i] = ctx.one
    lin = [matrix() for _ in range(n)]
    quad = {}
    for k in range(n):
        q = matrix()
        for i in range(n):
            q[i][i] = ctx.c((i+2)*(k+1) % 7-3)
        quad[(k, k)] = q
    q = matrix()
    q[0][n-1] = q[n-1][0] = ctx.c(2)
    quad[(0, 1)] = q
    g = F.base_from_coeffs(g0, lin, quad, {}, n, ctx)
    base = F.base_geometry(g, n, ctx)
    y = sp.diag(*range(2, n+2))
    y[0, 1] = y[1, 0] = 1
    y0 = [[ctx.c(int(y[i, j])) for j in range(n)] for i in range(n)]
    cH, cV = sp.Integer(2), sp.Integer(3)
    total, m, idx, _ = F.total_metric(
        g, base, n, y0, ctx.c(lam), ctx.c(int(cH)), ctx.c(int(cV)), 'Z1', ctx)
    ric, _, _ = F.ricci_at_point(total, m, ctx)
    yi = y.inv()
    hi = yi/cH
    basis = []
    for i, j in idx:
        e = sp.zeros(n)
        e[i, j] = e[j, i] = 1
        basis.append(e)

    def rational(value):
        return sp.Rational(str(value))

    def vertical(u, v):
        return cV*(sp.trace(yi*u*yi*v)-lam*sp.trace(yi*u)*sp.trace(yi*v))

    As = {}
    for a in range(n):
        for b in range(n):
            r = sp.Matrix(n, n, lambda i, j: rational(base['R'][i, j, a, b]))
            As[a, b] = (r.T*y+y*r)/2
    pair = [(a, b) for a in range(n) for b in range(n)]
    av = sp.Matrix(len(pair), len(basis),
                   lambda row, col: vertical(As[pair[row]], basis[col]))
    hin = sp.Matrix(len(pair), len(pair),
                    lambda i, j: hi[pair[i][0], pair[j][0]]*hi[pair[i][1], pair[j][1]])
    correction = av.T*hin*av
    flat = sp.Matrix(
        len(basis), len(basis),
        lambda i, j: -(n+1)*sp.trace(yi*basis[i]*yi*basis[j])/4
        +sp.trace(yi*basis[i])*sp.trace(yi*basis[j])/4)
    actual = sp.Matrix(len(basis), len(basis), lambda i, j: rational(ric[n+i][n+j]))
    residual = actual-flat-correction
    prefix = f'n={n}, lambda={lam}'
    check(prefix+': positive off-section fibre and block parameters',
          all(y[:k, :k].det() > 0 for k in range(1, n+1))
          and y != sp.eye(n) and cH > 0 and cV > 0 and lam < sp.Rational(1, n))
    check(prefix+f': all {len(basis)**2} vertical Ricci entries match',
          residual == sp.zeros(len(basis)))
    check(prefix+': the curvature-action correction is nonzero',
          any(v != 0 for value in As.values() for v in value))
    check(prefix+': every curvature action is orthogonal to Euler',
          all(sp.trace(yi*value) == 0 for value in As.values()))
    check(prefix+': exact correction rank respects the skew-horizontal bound',
          correction.rank() <= n*(n-1)//2)
    print(f'  correction rank={correction.rank()}; metric control {perf_counter()-started:.3f}s',
          flush=True)


def critical_curve_controls():
    a, m = sp.symbols('a m', real=True)
    t, z = sp.symbols('t z', positive=True)
    delta = -sp.Rational(8, 15)*(208*a+1)*(880*a+3)
    U = (208*a+1)*(21840*a+113)/480
    V = -sp.Rational(64, 285)*(402896*a+1719)
    W = sp.Rational(1024, 285)*(53-676*a)
    D = V-sp.Rational(11, 2)*W
    qm = 24+2*m
    A = qm*(delta*sp.diff(U, a)-2*W*U)-8*D*U
    B = 48*D+12*qm*(W-48)
    H = ((4629031680*m+48048349440)*a*a
         +(48524736*m+404087008)*a+128681*m+822893)
    N = (43264*m+3266432)*a-2537*m+20620
    check('critical-curve left coefficient has the exact universal factorization',
          sp.expand(A+sp.Rational(16, 4275)*(208*a+1)*H) == 0)
    check('critical-curve right coefficient has the exact universal factorization',
          sp.expand(B+sp.Rational(128, 95)*N) == 0)
    polys = [sp.Poly((t*(208*a+1)*H-360*N).subs(m, rank), a) for rank in range(4, 15)]
    check('all eleven retained rank equations are genuine cubics at positive scale',
          all(p.degree() == 3 and p.LC().is_positive for p in polys))
    left = -sp.Rational(1, 160)-z
    positive_left = True
    for rank in range(4, 15):
        hp = sp.Poly(sp.expand(H.subs({a: left, m: rank})), z)
        np = sp.Poly(sp.expand(-N.subs({a: left, m: rank})), z)
        positive_left &= all(v > 0 for v in hp.all_coeffs()+np.all_coeffs())
    check('scalar critical-rank curves still admit positive left-arc scale candidates',
          positive_left)


if __name__ == '__main__':
    started = perf_counter()
    for n, lam in [(2, sp.Rational(1, 3)), (3, sp.Rational(1, 4)),
                   (4, sp.Rational(1, 5)), (4, sp.Rational(1, 7))]:
        metric_control(n, lam)
    critical_curve_controls()
    print('Finite metric controls support the hand proof; scalar candidates do not construct full fields.',
          flush=True)
    print(f'{len(CHECKS)} exact controls passed; runtime {perf_counter()-started:.3f} s.', flush=True)
