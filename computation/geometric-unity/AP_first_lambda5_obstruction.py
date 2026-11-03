"""Exact first Lambda5 c/wGamma conic obstruction, including moving odd fields.

Uses accepted T's Clifford engine; the complete finite basis-orbit
coverage is proved in AP-first-lambda5-obstruction.md. All slots use real symbolic
weights, not sampled values. No metric or PDE solution is constructed.
"""
from pathlib import Path
import sys
import time
import sympy as sp

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import T_singular_residual as T

X, ALG = T.X, T.ALG
MEMBER = (1, 0, 0, 0)
CHECKS = []


def check(label, condition):
    ok = bool(condition)
    CHECKS.append(ok)
    print(f'[{len(CHECKS):02d}] {label}: {ok}', flush=True)
    assert ok, label


def grade(value, degree):
    return {mask: sp.expand(coefficient) for mask, coefficient in value.items()
            if mask.bit_count() == degree and sp.expand(coefficient) != 0}


def form(indices, coefficient=1):
    return {sum(1 << index for index in indices): coefficient}


def lift(omega, degree):
    return [grade(ALG.mul({1 << a: 1}, omega), degree) for a in range(X.N)]


def kernel(omega, r, s):
    return T.scale(X.fplus(T.scale(lift(omega, 4), r),
                          T.scale(ALG.gamma_right(lift(omega, 6)), s)), sp.I)


def outgrade(field, degree, factor=1):
    return X.normalized(T.scale([grade(value, degree) for value in field], factor))


def predicted(omega, r, s, degree):
    square = ALG.mul(omega, omega)
    sandwich, wedge = ((sp.Rational(28, 3)*r*r-36*s*s,
                        -sp.Rational(8, 3)*r*r+sp.Rational(112, 3)*s*s)
                       if degree == 5 else
                       (sp.Rational(40, 3)*r*r-40*s*s,
                        -24*r*r-sp.Rational(16, 3)*s*s))
    pure = [
        X.plus(X.scaled(grade(ALG.mul(ALG.mul(omega, {1 << a: 1}), omega), degree), sandwich),
               X.scaled(grade(ALG.mul({1 << a: 1}, square), degree), wedge))
        for a in range(X.N)]
    other = 14-degree
    cross_sandwich, cross_wedge = ((sp.Rational(104, 3), -sp.Rational(136, 3))
                                   if degree == 5 else (sp.Rational(128, 3), -sp.Integer(24)))
    cross = ALG.gamma_right([
        X.plus(X.scaled(grade(ALG.mul(ALG.mul(omega, {1 << a: 1}), omega), other), cross_sandwich),
               X.scaled(grade(ALG.mul({1 << a: 1}, square), other), cross_wedge))
        for a in range(X.N)])
    return X.normalized(X.fplus(pure, T.scale(cross, r*s)))


def contract(field):
    result = {}
    for a, value in enumerate(field):
        result = X.plus(result, ALG.mul({1 << a: 1}, value))
    return {mask: sp.expand(value) for mask, value in result.items() if sp.expand(value) != 0}


def support(field):
    return {mask.bit_count() for value in X.normalized(field) for mask in value}


def main():
    start = time.time()
    print('AP: first Lambda5 c/wGamma conic, exact all-slot orbit certificates.', flush=True)
    A, B, r, s, c, b = sp.symbols('A B r s c b', real=True)
    for intersection in range(6):
        I = tuple(range(5))
        J = tuple(range(intersection)) + tuple(range(5, 10-intersection))
        omega = X.plus(form(I, A), form(J, B))
        field = kernel(omega, r, s)
        output = T.quadratic(field, MEMBER)
        check(f'omega pair intersection {intersection}: compact phase and complete odd support',
              T.antihermitian(field) and support(output) <= {1, 5, 9, 13})
        for degree, factor in ((5, -64*r*r+sp.Rational(1552, 3)*s*s),
                               (9, -sp.Rational(592, 3)*r*r+128*s*s)):
            actual = outgrade(output, degree, 1/sp.I)
            check(f'intersection {intersection}: complete grade {degree} invariant identity',
                  actual == predicted(omega, r, s, degree))
            square = ALG.mul(omega, omega)
            square_component = X.scaled(grade(square, degree-1), factor)
            other_square = grade(square, 8 if degree == 5 else 4)
            converted_square = ALG.mul(other_square, {X.FULL: 1})
            square_component = X.plus(square_component,
                                      X.scaled(converted_square, -sp.Rational(1232, 3)*r*s))
            check(f'intersection {intersection}: gamma trace at grade {degree}',
                  contract(actual) == {mask: sp.expand(value)
                                       for mask, value in square_component.items()
                                       if sp.expand(value) != 0})

    omega = form(tuple(range(5)), A)
    field = kernel(omega, r, s)
    for direction in (0, 5):
        p2, p14 = ALG.forms(1, (direction,))[1::2]
        vector = X.fplus(T.scale(p2, c), T.scale(p14, b))
        cross = T.mix(vector, field, MEMBER)
        check(f'C axis {direction}: arbitrary vector weights and kernel mixed support',
              T.antihermitian(vector) and support(cross) <= {3, 7, 11})
        check(f'jet axis {direction}: all weighted-kernel derivatives retain only grades 3/7/11',
              support(T.derivative(field, direction, MEMBER)) <= {3, 7, 11})
    for direction in (0, 1):
        p2, p14 = ALG.forms(1, (0,))[1::2]
        vector = X.fplus(T.scale(p2, c), T.scale(p14, b))
        check(f'vector jet axis {direction}: no grade 5/9',
              not (support(T.derivative(vector, direction, MEMBER)) & {5, 9}))
    p2, p14 = ALG.forms(1, (0,))[1::2]
    other_p2, other_p14 = ALG.forms(1, (1,))[1::2]
    vector = X.fplus(X.fplus(T.scale(p2, c), T.scale(p14, b)),
                    X.fplus(T.scale(other_p2, A), T.scale(other_p14, B)))
    check('quadratic vector field with independent covectors has no grade 5/9',
          not (support(T.quadratic(vector, MEMBER)) & {5, 9}))

    matrix = sp.Matrix([[sp.Rational(28, 3)*r*r-36*s*s, sp.Rational(104, 3)*r*s],
                        [-sp.Rational(128, 3)*r*s, sp.Rational(40, 3)*r*r-40*s*s]])
    determinant = sp.factor(matrix.det())
    positive = sp.Rational(32, 9)*(35*r**4+176*r*r*s*s+405*s**4)
    check('scalar-square reduction has a strictly positive two-by-two determinant',
          sp.expand(determinant-positive) == 0)
    print('Residual determinant =', determinant, flush=True)
    check('cross traces extract complementary square components with nonzero coefficient',
          -4*sp.Rational(104, 3)+6*(-sp.Rational(136, 3)) == -sp.Rational(1232, 3)
          and -4*sp.Rational(128, 3)+10*(-24) == -sp.Rational(1232, 3))
    check('at either zero-weight endpoint both ordinary trace factors remain nonzero',
          all(value != 0 for value in (-64, sp.Rational(1552, 3),
                                       -sp.Rational(592, 3), 128)))
    print('No derivative was set to zero before the grade 5/9 separation.', flush=True)
    print('Universality uses the complete SO14 equivariance/polarization proof in AP-first-lambda5-obstruction.md.', flush=True)
    print(f'{len(CHECKS)} exact controls passed; runtime {time.time()-start:.3f} s.', flush=True)


if __name__ == '__main__':
    main()
