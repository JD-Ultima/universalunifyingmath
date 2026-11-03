"""Exact controls for AQ's opposite moving-Lambda5 conservation reduction.

Only the w/cGamma pair is treated. Imports accepted T/AH operations, never
their test mains. Complete symbolic grade-one calibrations and invariant
contraction arguments establish the stress formula. The hand covariant and
density proofs establish local isoparametricity; no PDE solution is produced.
"""
from pathlib import Path
import sys
from time import perf_counter

import sympy as sp

sys.path.insert(0, str(Path(__file__).resolve().parent))
import T_singular_residual as T
from AH_lambda5_gram import kernel, stress, grade_one, form

X, ALG = T.X, T.ALG
MEMBER = (1, 0, 0, 0)
CHECKS = []


def check(label, condition):
    good = bool(condition)
    CHECKS.append(good)
    print(f'[{len(CHECKS):02d}] {label}: {good}', flush=True)
    assert good, label


def vector(C, B):
    result = [{} for _ in range(X.N)]
    for a in range(X.N):
        copies = ALG.forms(1, (a,))
        result = X.fplus(result, X.fplus(T.scale(copies[1], C[a]),
                                        T.scale(copies[3], B[a])))
    return X.normalized(result)


def stress_controls():
    r, s, c, b, t, u, kappa = sp.symbols('r s c b t u kappa', real=True)
    identity = sp.eye(X.N)
    C = sp.Matrix([c] + [0] * 13)
    value = vector(list(C), [b] + [0] * 13)
    check('complete symbolic aligned-vector grade-one quadratic',
          grade_one(T.quadratic(value, MEMBER)) == 224*C*C.T + 88*c*c*identity)
    # End(V) -> End(V) has the three invariant linear patterns H,H^T,tr(H)I.
    for direction, component in ((0, 0), (0, 1), (1, 0)):
        matrix = sp.zeros(X.N)
        matrix[direction, component] = 1
        for copy in (1, 3):
            jet = ALG.forms(1, (component,))[copy]
            expected = (-24*(matrix.T-sp.trace(matrix)*identity)
                        if copy == 1 else sp.zeros(X.N))
            check(f'complete vector derivative orbit ({direction},{component}), copy{copy}',
                  grade_one(T.derivative(jet, direction, MEMBER)) == expected)
    simple = form((0, 1, 2, 3, 4))
    full = X.plus(X.plus(X.plus(simple, form((0, 1, 5, 6, 7), -2)),
                         form((0, 8, 9, 10, 11), 3)), form((1, 3, 6, 12, 13), -1))
    q = sp.Rational(112, 3)*r*r + 40*s*s
    ell = sp.Rational(320, 3)*r*r - 32*s*s
    for label, omega in (('unit-simple symbolic calibration', simple),
                         ('nondecomposable full-support control', full)):
        extra = kernel(omega, 1, r, s)
        norm, gram = stress(omega)
        check(label + ': all196 stress coefficients',
              grade_one(T.quadratic(extra, MEMBER)) == q*norm*identity + ell*gram)
        check(label + ': entire mixed grade-one output vanishes',
              grade_one(T.mix(value, extra, MEMBER)) == sp.zeros(X.N))
    check('complete invariant grade-one baseline',
          grade_one(T.baseline(t, u, kappa, [0]*X.N))
          == (kappa*t+312*t*t+104*u*u)*identity)
    extra = kernel(full, 1, r, s)
    check('all14 five-form derivative directions have zero grade-one output',
          all(grade_one(T.derivative(extra, j, MEMBER)) == sp.zeros(X.N)
              for j in range(X.N)))


def scalar_controls():
    a, w, c, kappa = sp.symbols('a w c kappa', real=True)
    A, B = 1+sp.Rational(704, 3)*a, 1+160*a
    f = sp.expand(A*B/1280)
    delta = sp.expand(1+224*a-3328*f)
    r, s = sp.Rational(80, 3)*w, A
    q = sp.Rational(112, 3)*r*r+40*s*s
    ell = sp.Rational(320, 3)*r*r-32*s*s
    on_conic = lambda expression: sp.factor(sp.cancel(expression).subs(w*w, f))
    norm = w/(r*s)*(14*kappa*kappa*(1+208*a)-43264*c)/240
    beta = sp.Rational(624, 19)*w/(r*s)
    alpha = (5*norm+beta*c)/14
    F = a+312*a*a+104*f
    scalar = on_conic(kappa*kappa*F+88*c+q*norm+ell*alpha)
    U = (208*a+1)*(21840*a+113)/480
    V = -sp.Rational(64, 285)*(402896*a+1719)
    W = sp.Rational(1024, 285)*(53-676*a)
    Ac = sp.Rational(32, 285)*(816608*a+5155)
    check('exact grade-one scalar and rank-one elimination',
          sp.factor(scalar-kappa*kappa*U-c*V) == 0
          and on_conic(224-ell*beta-W) == 0)
    check('symmetric-Gamma scalar equation retained',
          on_conic(-kappa*kappa*w*(1+208*a)+3328*w*c
                   +sp.Rational(160, 3)*r*s*norm-sp.Rational(304, 3)*r*s*alpha) == 0)
    check('symmetric-Gamma rank-one equation retained',
          on_conic(-3328*w+sp.Rational(304, 3)*r*s*beta) == 0)
    check('contraction-Gram trace retained', on_conic(14*alpha-beta*c-5*norm) == 0)
    ac_zero = -sp.Rational(5155, 816608)
    check('transverse Bianchi coefficient and its single zero',
          sp.factor(288-W/2-V-Ac) == 0
          and Ac.subs(a, ac_zero) == 0 and sp.diff(Ac, a) != 0)
    check('the transverse zero is on the real left arc, so cannot be divided away',
          f.subs(a, ac_zero) == sp.Rational(36853, 9768290415)
          and ac_zero < -sp.Rational(1, 160))
    h, P, z, cc = sp.symbols('h P W c')
    trace_G = -312*h-14*P-z*cc
    ric_scalar = sp.expand(-24*h-P-trace_G/12)
    check('dimension14 trace reversal',
          sp.expand(ric_scalar-(2*h+P/6+z*cc/12)) == 0)
    da, dc, dP, Cslot = sp.symbols('da dc dP Cslot')
    direct = (24*(12*dc+(2*h+P/6-sp.Rational(11, 12)*z*cc)*Cslot)
              -dP-z*(h*Cslot+dc/2)-sp.diff(W, a)*cc*da)
    collected = ((288-z/2)*dc+((48-z)*h+4*P-22*z*cc)*Cslot
                 -dP-sp.diff(W, a)*cc*da)
    check('complete divergence collection', sp.expand(direct-collected) == 0)
    # Substitute dP=kappa^2 U' da+V' c da+V dc and da=delta C/24.
    actual = -collected.subs({z: W, cc: c, P: kappa*kappa*U+V*c,
                              dP: (kappa*kappa*sp.diff(U, a)+sp.diff(V, a)*c)*da+V*dc})
    actual = actual.subs(da, delta*Cslot/24)
    expected = ((V+W/2-288)*dc
                +((W-48)*h+kappa*kappa*(delta*sp.diff(U, a)/24-4*U)
                  +c*(delta*(sp.diff(V, a)+sp.diff(W, a))/24+22*W-4*V))*Cslot)
    check('full covector Bianchi identity with every coefficient derivative',
          sp.expand(actual-expected) == 0)
    w_zero = sp.Rational(2537, 43264)
    check('the mean-curvature coefficient has exactly one exceptional level',
          sp.factor(W-48-sp.Rational(16, 285)*(2537-43264*a)) == 0
          and W.subs(a, w_zero) == 48 and sp.diff(W, a) != 0)
    print('  U =', sp.factor(U), flush=True)
    print('  V =', sp.factor(V), flush=True)
    print('  W =', sp.factor(W), flush=True)
    print('  A_c =', sp.factor(Ac), flush=True)
    print('  Exceptional a levels:', ac_zero, 'and', w_zero, flush=True)


def main():
    start = perf_counter()
    print('AQ: opposite moving Lambda5 stress and Bianchi compatibility.', flush=True)
    print('Real compact WA; all derivatives retained; no existence claim.', flush=True)
    stress_controls()
    scalar_controls()
    print('Local isoparametricity uses the hand covariant/density proof, not finite sampling.', flush=True)
    print(f'{len(CHECKS)} exact controls passed; runtime {perf_counter()-start:.3f} s.', flush=True)


if __name__ == '__main__':
    main()
