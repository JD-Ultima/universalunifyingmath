"""Exact full-system reduction on the opposite moving Lambda5 conic.

Uses accepted T and AL helpers, without
running their suites. AS-opposite-lambda5-system.md supplies the full basis proof.
"""
from pathlib import Path
import sys
import time
import sympy as sp

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import AL_opposite_lambda5_system as L
T, X, ALG = L.T, L.X, L.ALG
CHECKS = []


def check(label, condition):
    good = bool(condition)
    CHECKS.append(good)
    print(f'[{len(CHECKS):02d}] {label}: {good}', flush=True)
    assert good, label


def kernel(omega, r, s):
    return X.fplus(T.scale(L.lift(omega, 'w'), r), T.scale(L.lift(omega, 'cG'), s))


def same(left, right):
    return not any(X.normalized(X.fplus(left, right, -1)))


def tensor_bases(omega, axis):
    wedge = L.grade(ALG.mul({1 << axis: 1}, omega), 6)
    ctr = L.grade(ALG.mul({1 << axis: 1}, omega), 4)
    return [[omega if a == axis else {} for a in range(14)],
            [L.grade(ALG.mul({1 << a: 1}, wedge), 5) for a in range(14)],
            [L.grade(ALG.mul({1 << a: 1}, ctr), 5) for a in range(14)]]


def combine(bases, coefficients):
    out = [{} for _ in range(14)]
    for basis, coefficient in zip(bases, coefficients):
        out = X.fplus(out, T.scale(basis, coefficient))
    return X.normalized(out)


def mixed_coefficients(r, s, c, b, degree):
    if degree == 5:
        return (-sp.Rational(16, 3)*b*s+304*c*r,
                sp.Rational(16, 3)*b*s-sp.Rational(208, 3)*c*r,
                sp.Rational(16, 3)*b*s+32*c*r)
    return (-sp.Rational(32, 3)*b*r-160*c*s,
            sp.Rational(32, 3)*b*r+sp.Rational(80, 3)*c*s,
            sp.Rational(32, 3)*b*r)


def predicted_Q(omega, r, s, degree):
    square = ALG.mul(omega, omega)
    c1, c2, c3, c4 = ((52*r*r-sp.Rational(44, 3)*s*s, 48*r*r+sp.Rational(8, 3)*s*s,
                       56*r*s, 40*r*s)
                      if degree == 5 else
                      (sp.Rational(152, 3)*r*r-sp.Rational(40, 3)*s*s,
                       sp.Rational(16, 3)*r*r-sp.Rational(56, 3)*s*s,
                       sp.Rational(160, 3)*r*s, sp.Rational(56, 3)*r*s))
    other = 14-degree
    pure, cross = [], []
    for a in range(14):
        sandwich = ALG.mul(ALG.mul(omega, {1 << a: 1}), omega)
        wedge = ALG.mul({1 << a: 1}, square)
        pure.append(X.plus(X.scaled(L.grade(sandwich, degree), c1),
                           X.scaled(L.grade(wedge, degree), c2)))
        cross.append(X.plus(X.scaled(L.grade(sandwich, other), c3),
                            X.scaled(L.grade(wedge, other), c4)))
    return X.normalized(X.fplus(pure, ALG.gamma_right(cross)))


def actual_derivative(jet_r, jet_s):
    out = [{} for _ in range(14)]
    for a in range(14):
        out = X.fplus(out, T.derivative(kernel(jet_r[a], 1, 0), a, L.MEMBER))
        out = X.fplus(out, T.derivative(kernel(jet_s[a], 0, 1), a, L.MEMBER))
    return X.normalized(out)


def main():
    start = time.time()
    print('AS: opposite moving Lambda5 w/cGamma, full tensor and jet system.', flush=True)
    A, B, r, s, c, b = sp.symbols('A B r s c b', real=True)
    for overlap in range(6):
        omega = X.plus(L.form(range(5), A),
                       L.form(tuple(range(overlap))+tuple(range(5, 10-overlap)), B))
        field = kernel(omega, r, s)
        Q = T.quadratic(field, L.MEMBER)
        check(f'quadratic pair intersection {overlap}: correct phase and full support',
              T.antihermitian(field)
              and {mask.bit_count() for item in X.normalized(Q) for mask in item} <= {1, 5, 9, 13})
        norm, gram = L.stress(omega)
        converted = ALG.gamma_right(L.real_grade(Q, 13))
        full_gamma = sp.Matrix([[sp.expand(value.get(1 << column, 0))
                                for column in range(X.N)] for value in converted])
        target_gamma = (sp.Rational(160, 3)*r*s*norm*sp.eye(X.N)
                        -sp.Rational(304, 3)*r*s*gram).applyfunc(sp.expand)
        check(f'intersection {overlap}: entire grade-thirteen tensor, without symmetrization',
              full_gamma == target_gamma)
        for degree in (5, 9):
            actual = L.real_grade(Q, degree)
            check(f'intersection {overlap}: full grade {degree} identity including rs terms',
                  same(actual, predicted_Q(omega, r, s, degree)))
            square = ALG.mul(omega, omega)
            factor, cross = ((272*r*r+sp.Rational(256, 3)*s*s, 16*r*s)
                             if degree == 5 else
                             (-sp.Rational(512, 3)*r*r-sp.Rational(176, 3)*s*s,
                              -sp.Rational(80, 3)*r*s))
            target = X.plus(X.scaled(L.grade(square, degree-1), factor),
                            X.scaled(ALG.mul(L.grade(square, 8 if degree == 5 else 4), {X.FULL: 1}), cross))
            check(f'intersection {overlap}: full grade {degree} gamma contraction',
                  L.contraction(actual) == {mask: sp.expand(value) for mask, value in target.items()
                                            if sp.expand(value) != 0})
    omega = L.form(range(5), A)
    for axis in (0, 5):
        p2, p14 = ALG.forms(1, (axis,))[1::2]
        v = X.fplus(T.scale(p2, c), T.scale(p14, b))
        mixed = T.mix(v, kernel(omega, r, s), L.MEMBER)
        check(f'mixed vector/blade orbit {axis}: full support only 5/9',
              {mask.bit_count() for item in X.normalized(mixed) for mask in item} <= {5, 9})
        for degree in (5, 9):
            actual = L.real_grade(mixed, degree)
            if degree == 9:
                actual = ALG.gamma_right(actual)
            check(f'mixed orbit {axis}: exact converted grade {degree} invariant coefficients',
                  same(actual, combine(tensor_bases(omega, axis), mixed_coefficients(r, s, c, b, degree))))
        jet = [omega if a == axis else {} for a in range(14)]
        check(f'derivative orbit {axis}: full weighted derivative and right-Gamma sign',
              same(T.derivative(kernel(omega, r, s), axis, L.MEMBER),
                   X.fplus(T.scale(L.K(jet), -16*sp.I*r),
                           T.scale(ALG.gamma_right(L.K(jet)), -8*sp.I*s))))

    a, w, delta = sp.symbols('a w delta', real=True)
    da, ds, dr = sp.symbols('da ds dr', real=True)
    A0, B0 = 1+sp.Rational(704, 3)*a, 1+160*a
    f = sp.expand(A0*B0/1280)
    # Chart-rescaling-invariant Wronskian identity for a smooth nonzero-w chart.
    W0 = sp.Rational(80, 3)*(w*sp.diff(A0, a)-A0*sp.diff(f, a)/(2*w))
    check('weight Wronskian identity retains all variable-weight derivatives',
          sp.cancel(sp.cancel(W0-sp.Rational(7, 240)*(sp.Rational(80, 3)*w*A0)/w**2).subs(w*w, f)) == 0)
    effective = sp.Rational(3328, 3)*w*(s*s+4*r*r)/(r*s)-sp.Rational(7, 360)*delta/w**2
    combined = [sp.expand(s*x+2*r*y) for x, y in zip(mixed_coefficients(r, s, 1, 208*w, 5),
                                                   mixed_coefficients(r, s, 1, 208*w, 9))]
    weight = sp.Rational(7, 360)*delta/w**2*r*s
    adjusted = [combined[0]+weight, combined[1]-weight, combined[2]-weight]
    target = [r*s*(-16-effective), r*s*(-16+effective), r*s*(32+effective)]
    check('elimination combines both mixed terms and the Wronskian exactly',
          all(sp.cancel(x-y) == 0 for x, y in zip(adjusted, target)))
    check('eliminated quadratic tensor has the exact four invariant coefficients',
          all(sp.expand(x-y) == 0 for x, y in zip(
              [s*(52*r*r-sp.Rational(44, 3)*s*s)-2*r*(sp.Rational(160, 3)*r*s),
               s*(48*r*r+sp.Rational(8, 3)*s*s)-2*r*(sp.Rational(56, 3)*r*s),
               56*r*s*s+2*r*(sp.Rational(152, 3)*r*r-sp.Rational(40, 3)*s*s),
               40*r*s*s+2*r*(sp.Rational(16, 3)*r*r-sp.Rational(56, 3)*s*s)],
              [-sp.Rational(4, 3)*s*(41*r*r+11*s*s),
               sp.Rational(8, 3)*s*(4*r*r+s*s),
               sp.Rational(8, 3)*r*(38*r*r+11*s*s),
               sp.Rational(8, 3)*r*(4*r*r+s*s)])))
    check('eliminated gamma traces retain separate grade-four and grade-six channels',
          sp.expand(s*(272*r*r+sp.Rational(256, 3)*s*s)+2*r*(sp.Rational(80, 3)*r*s)
                    -sp.Rational(16, 3)*s*(61*r*r+16*s*s)) == 0
          and sp.expand(16*r*s*s+2*r*(-sp.Rational(512, 3)*r*r-sp.Rational(176, 3)*s*s)
                        +sp.Rational(16, 3)*r*(64*r*r+19*s*s)) == 0)
    check('combined gradient traces give the exact contraction/exterior coefficients',
          sp.expand(target[0]+10*target[2]-r*s*(304+9*effective)) == 0
          and sp.expand(target[0]+6*target[1]-r*s*(-112+5*effective)) == 0)

    # Full nonlinear control with independent nonzero covariant and weight jets.
    omega = X.plus(X.plus(L.form(range(5)), L.form((1, 2, 5, 6, 7), 2)),
                   L.form((0, 8, 10, 12, 13), -3))
    jet = [{} for _ in range(14)]
    jet[1] = L.form((0, 3, 6, 9, 12), sp.Rational(2, 7))
    jet[9] = X.plus(L.form((2, 3, 5, 7, 11)), L.form((0, 1, 4, 8, 13), -2))
    r0, s0 = sp.Rational(2, 3), -sp.Rational(3, 5)
    dr0 = [sp.Rational(1, 11) if a == 0 else 0 for a in range(14)]
    ds0 = [-sp.Rational(2, 13) if a == 8 else 0 for a in range(14)]
    jr = [X.plus(X.scaled(jet[a], r0), X.scaled(omega, dr0[a])) for a in range(14)]
    js = [X.plus(X.scaled(jet[a], s0), X.scaled(omega, ds0[a])) for a in range(14)]
    actual = actual_derivative(jr, js)
    check('general full-support jet with independent dr/ds: full derivative identity',
          same(actual, X.fplus(T.scale(L.K(jr), -16*sp.I),
                               T.scale(ALG.gamma_right(L.K(js)), -8*sp.I))))
    d5 = L.real_grade(actual, 5)
    d9 = ALG.gamma_right(L.real_grade(actual, 9))
    eliminated = X.fplus(T.scale(d5, s0), T.scale(d9, 2*r0))
    wronskian_jet = [X.scaled(omega, r0*ds0[a]-s0*dr0[a]) for a in range(14)]
    check('general full-support jet: exact common-jet cancellation with nonzero Wronskian',
          same(eliminated, T.scale(L.K(wronskian_jet), 16)))
    check('the full K inverse retains all jet components in the moving system',
          same(L.K(L.K(jet), inverse=True), jet))
    print('Exact higher system: -16 K(nabla(r omega))+Z5=0; +8 K(nabla(s omega))+Z9=0.', flush=True)
    print('Algebraic compatibility: s Z5+2r Z9+16 K((r ds-s dr) tensor omega)=0.', flush=True)
    print('Full equivalence also retains AQ1.1, the actual form Gram tensor and the entire AH2 equation.', flush=True)
    print('No simplicity, observer-metric realization or full-solution existence is inferred.', flush=True)
    print(f'{len(CHECKS)} exact controls passed; runtime {time.time()-start:.3f} s.', flush=True)


if __name__ == '__main__':
    main()
