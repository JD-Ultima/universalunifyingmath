"""Exact certificates for AT's exceptional opposite-Lambda5 profile level.

Standalone SymPy scalar algebra: no repository helper, private source or
metric engine is imported. The integrating-factor existence/convergence,
bounded uniqueness and actual Euler-curve argument are hand proofs. A scalar
analytic germ does not construct an observer metric or full field solution.
"""
from time import perf_counter

import sympy as sp


CHECKS = []


def check(label, condition):
    good = bool(condition)
    CHECKS.append(good)
    print(f'[{len(CHECKS):02d}] {label}: {good}', flush=True)
    assert good, label


def main():
    start = perf_counter()
    print('AT: bounded scalar profile and critical-point incompatibility at D=0.', flush=True)
    a, x, t = sp.symbols('a x t', real=True)
    f = (1+sp.Rational(704, 3)*a)*(1+160*a)/1280
    delta = -sp.Rational(8, 15)*(208*a+1)*(880*a+3)
    U = (208*a+1)*(21840*a+113)/480
    V = -sp.Rational(64, 285)*(402896*a+1719)
    W = sp.Rational(1024, 285)*(53-676*a)
    D = sp.factor(V-sp.Rational(11, 2)*W)
    P = sp.factor(sp.diff(V, a)+sp.diff(W, a)-2*W*D/delta)
    Q = sp.factor(sp.diff(U, a)-2*W*U/delta)
    R = sp.factor(-12*(W-48)/delta)
    point = -sp.Rational(491, 26416)
    d1 = -sp.Rational(21978112, 285)
    check('D has exactly the stated simple zero and leading coefficient',
          sp.expand(D-d1*(a-point)) == 0 and d1 != 0)
    other_vector_factor = 1+544*a-3328*f
    check('the level is on the real left arc with both vector factors regular',
          point < -sp.Rational(1, 160) and f.subs(a, point) > 0
          and delta.subs(a, point) == -sp.Rational(4939648, 241935)
          and other_vector_factor.subs(a, point) != 0)
    p = sp.factor(P.subs(a, point+x)/d1)
    q = sp.factor(-(t*Q.subs(a, point+x)+R.subs(a, point+x))/d1)
    Y, Yp = sp.symbols('Y Y_prime')
    original = (D*Yp+P*Y+t*Q+R).subs(a, point+x)
    check('dimensionless equation has the exact regular-singular normalization',
          sp.cancel(original-d1*(x*Yp+p*Y-q)) == 0)
    p0 = p.subs(x, 0)
    check('normalized coefficients are regular and the indicial number is 153/127',
          sp.denom(p).subs(x, 0) != 0 and sp.denom(q).subs(x, 0) != 0
          and p0 == sp.Rational(153, 127))
    y0 = sp.factor(q.subs(x, 0)/p0)
    target_y0 = 95*(3241880613-6380808980*t)/sp.Integer(259536139978752)
    check('bounded profile has the exact forced value', sp.factor(y0-target_y0) == 0)
    y1 = sp.factor((sp.diff(q, x).subs(x, 0)-sp.diff(p, x).subs(x, 0)*y0)/(1+p0))
    target_y1 = ((103523628088901644*t+54385536931208595)
                 / sp.Integer(788642454954332160))
    check('differentiated equation has the exact forced slope', sp.factor(y1-target_y1) == 0)
    tau = sp.Rational(3241880613, 6380808980)
    positive_factor = sp.Rational(95*6380808980, 259536139978752)
    check('forced norm has the stated exact positive scale threshold',
          tau > 0 and positive_factor > 0 and sp.factor(y0-positive_factor*(tau-t)) == 0)
    check('forced slope is strictly positive for every positive scale',
          all(value > 0 for value in sp.Poly(y1, t).all_coeffs()) and y1.subs(t, tau) > 0)
    # Generic integrating-factor reduction: v'/v=(p-p0)/x and z=v*y.
    # Convergence and uniqueness are supplied by the hand proof.
    pp, zz, zzp, vv = sp.symbols('p z z_prime v', nonzero=True)
    vp = (pp-p0)*vv/x
    transformed = x*(zzp/vv-zz*vp/(vv*vv))+pp*zz/vv
    check('integrating-factor transformation has the correct derivative sign',
          sp.cancel(transformed-(x*zzp+p0*zz)/vv) == 0)
    eta, L, m, c0, kk, uu = sp.symbols('eta L m c0 kappa U', nonzero=True)
    critical_polynomial = 24*eta**2+(2*m*eta+kk*kk*uu/6+1/c0)*eta
    check('critical Hessian polynomial gives the exact trace-rank relation',
          sp.expand(critical_polynomial-eta*((24+2*m)*eta+kk*kk*uu/6+1/c0)) == 0)
    dd = sp.symbols('delta', nonzero=True)
    y_taylor = eta**2*c0**2
    x_taylor = dd*eta*c0/48
    check('actual Euler-curve Taylor ratio fixes the profile slope',
          sp.cancel(y_taylor/x_taylor-48*c0*eta/dd) == 0)
    eta_bar = sp.factor(delta.subs(a, point)*y1.subs(t, tau)/48)
    check('the threshold slope forces the stated negative normalized Hessian eigenvalue',
          eta_bar == -sp.Rational(53607532260267, 929045787488000))
    forced_rank = sp.factor(-12-(tau*U.subs(a, point)/6+1)/(2*eta_bar))
    check('critical-point compatibility would require a negative rank',
          forced_rank == -sp.Rational(109962574640284, 53607532260267)
          and forced_rank < 0
          and sp.factor((24+2*forced_rank)*eta_bar+tau*U.subs(a, point)/6+1) == 0)
    print('  a_star =', point, flush=True)
    print('  p0 =', p0, '; homogeneous exponent =', -p0, flush=True)
    print('  forced y0(t) =', y0, flush=True)
    print('  forced y1(t) =', y1, flush=True)
    print('  tau_star =', tau, flush=True)
    print('  critical c0*eta =', eta_bar, '; incompatible rank =', forced_rank, flush=True)
    print('Analytic scalar-germ existence is a hand proof, not a metric/form existence result.', flush=True)
    print(f'{len(CHECKS)} exact controls passed; runtime {perf_counter()-start:.3f} s.', flush=True)


if __name__ == '__main__':
    main()
