"""Exact invariant traces on the four real WA exterior conics.

Uses T's compact-real-form helpers through the independently reconstructed
Clifford engine xcheck_round14; no J/R evaluator is imported. The proof of
the differential reduction and its limits is AA-conic-traces.md. The symbolic
checks concern both exterior pairs, all four conics and unrestricted trace
coefficients, rather than selected points on a curve.
"""
import pathlib
import sys
import time

import sympy as sp

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import T_singular_residual as T

X, ALG = T.X, T.ALG
MEMBER = (1, 0, 0, 0)
DATA = (
    ('Lambda3-c/wGamma', 3, 0, 2, sp.Rational(112, 3), 400, sp.Rational(352, 3), 32),
    ('Lambda3-w/cGamma', 3, 1, 3, 240, sp.Rational(272, 3), 32, sp.Rational(352, 3)),
    ('Lambda5-c/wGamma', 5, 0, 2, 96, sp.Rational(832, 3), 48, sp.Rational(80, 3)),
    ('Lambda5-w/cGamma', 5, 1, 3, sp.Rational(704, 3), 160, sp.Rational(80, 3), 48),
)
CHECKS = []


def check(label, value):
    ok = bool(value)
    CHECKS.append(ok)
    print(f'[{len(CHECKS):02d}] {label}: {ok}', flush=True)
    assert ok, label


def zero(field):
    return not any(X.normalized(field))


def on_conic(expr, w, f):
    return sp.cancel(sp.rem(sp.expand(expr), w*w-f, w)) == 0


def main():
    start = time.time()
    print('AA. Exact invariant traces on the four WA exterior conics', flush=True)
    print(f'sympy {sp.__version__}; compact (14,0), normalized -Tr pairing', flush=True)
    phi, partner = T.scale(ALG.phi, sp.I), T.scale(ALG.phig, sp.I)
    p2, p14 = ALG.forms(1, (0,))[1::2]
    c, b = sp.symbols('c b', real=True)
    vector = X.fplus(T.scale(p2, c), T.scale(p14, b))
    qv = T.quadratic(vector, MEMBER)
    check('vector phases, norms and both quadratic traces',
          T.antihermitian(vector) and T.pair(vector, vector) == 13*c*c+b*b
          and T.pair(phi, qv) == 1456*c*c
          and T.pair(partner, qv) == -208*b*c)
    check('vector derivative trace is 312 div(C); P14 trace is zero',
          T.pair(phi, T.derivative(p2, 0, MEMBER)) == 312
          and T.pair(phi, T.derivative(p2, 1, MEMBER)) == 0
          and all(T.pair(phi, T.derivative(p14, j, MEMBER)) == 0 for j in range(X.N)))
    check('Gamma source/current vanishes in all 14 directions',
          all(zero(T.derivative(partner, j, MEMBER)) for j in range(X.N)))

    a, w, xx, yy = sp.symbols('a w x y', real=True)
    kappa = sp.symbols('kappa', real=True, nonzero=True)
    trace_g, lap, grad = sp.symbols('traceG laplacian_a grad_a_squared', real=True)
    print('  A=1+d*a, B=1+e*a, P=p*q; f=A*B/P; delta=1+224*a-3328*f.', flush=True)
    for name, m, i, j, d, e, p, q in DATA:
        copies = ALG.forms(m, tuple(range(m)))
        one, two = T.real_even(copies[i]), T.real_even(copies[j])
        n1, n2 = T.pair(one, one), T.pair(two, two)
        extra = X.fplus(T.scale(one, xx), T.scale(two, yy))
        qk = T.quadratic(extra, MEMBER)
        check(name+' compact pair and unrestricted quadratic traces',
              T.antihermitian(extra) and p*n1 == q*n2
              and T.pair(phi, qk) == (d*n1*xx*xx+e*n2*yy*yy)/2
              and T.pair(partner, qk) == -p*n1*xx*yy)
        P, aa, bb = p*q, 1+d*a, 1+e*a
        h = aa+bb
        f = sp.expand(aa*bb/P)
        delta = sp.expand(1+224*a-3328*f)
        H = sp.cancel(P*sp.diff(f, a)/(2*h))
        kernel = X.fplus(T.scale(one, p*w), T.scale(two, aa))
        to = T.invariant(a, w)
        lk = X.fplus(kernel, T.mix(kernel, to, MEMBER))
        norm = T.pair(kernel, kernel)
        check(name+' symbolic kernel and reduced trace formulas',
              all(on_conic(value, w, f) for item in lk for value in item.values())
              and on_conic(2*h*T.pair(phi, T.quadratic(kernel, MEMBER))
                           -(d+e+2*d*e*a)*norm, w, f)
              and on_conic(h*T.pair(partner, T.quadratic(kernel, MEMBER))+P*w*norm, w, f))
        cross = T.mix(vector, extra, MEMBER)
        check(name+' both cross traces and all 28 non-vector currents vanish',
              T.pair(phi, cross) == T.pair(partner, cross) == 0
              and all(T.pair(phi, T.derivative(copy, direction, MEMBER)) == 0
                      and T.pair(partner, T.derivative(copy, direction, MEMBER)) == 0
                      for copy in (one, two) for direction in range(X.N)))

        N = h/P*(14*kappa*kappa*(1+208*a)-24920064*grad/delta**2)
        F1 = a+312*a*a+104*f
        phi_trace = (trace_g+14*kappa*kappa*F1+7488*lap/delta
                     +7488*(112-sp.diff(delta, a))*grad/delta**2+H*N)
        potential = sp.expand(14*F1+7*sp.diff(f, a)*(1+208*a))
        reduced = trace_g+kappa*kappa*potential+7488/delta*(lap-sp.diff(delta, a)*grad/(2*delta))
        residual = (a*trace_g+14*kappa*kappa*(1+312*a)*(a*a+f)+7488*a*lap/delta
                    +7488*(-112*a+3328*(a*sp.diff(f, a)-f))*grad/delta**2)
        check(name+' Gamma norm elimination and T2 residual identity',
              sp.cancel(phi_trace-reduced) == 0
              and sp.cancel(residual-N/2-a*reduced) == 0)
        left, right = sorted((-sp.Rational(1, d), -sp.Rational(1, e)))
        check(name+' exact left/right branch sign premises',
              left < -sp.Rational(1, 312) and 1+224*left < 0 and 1+208*left < 0
              and sp.diff(f, a).subs(a, left) < 0 and h.subs(a, left) < 0
              and right > -sp.Rational(1, 208)
              and sp.diff(f, a).subs(a, right) > 0 and h.subs(a, right) > 0)
        print(f'  {name}: d={d}, e={e}, P={P}; left<={left}, right>={right}', flush=True)
        print(f'    f={sp.factor(f)}; delta={sp.factor(delta)}; R={sp.factor(potential)}', flush=True)

    print('Gamma equation for w!=0: P*N/h +24920064*|da|^2/delta^2 =14*kappa^2*(1+208*a).', flush=True)
    print('Reduced trace: Delta(a)-delta_prime*|da|^2/(2*delta) =-delta*(traceG+kappa^2*R)/7488.', flush=True)
    print('On every left branch, delta<0 and R>0; traceG>=0 forbids an interior local maximum of a.', flush=True)
    print('The undivided Phi1 trace also excludes reaching a left endpoint u=0. No other odd tensor equation is solved.', flush=True)
    print(f'\n{len(CHECKS)} checks, all True; runtime {time.time()-start:.1f} s', flush=True)


if __name__ == '__main__':
    main()
