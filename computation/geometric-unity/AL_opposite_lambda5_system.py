"""Exact controls for AL's opposite-Lambda5 constant endpoint system.

Uses the accepted T Clifford engine and AH's positive form Gram calculation.
Two exhaustive coordinate-jet orbits certify the linear derivative identity;
six symbolic basis-pair orbits certify the quadratic identities. Necessary
point jets and first prolongations are not observer metrics or PDE solutions.
AL4's conditional uniqueness follows from the pathwise smooth ODE proof;
these finite algebraic controls are not an existence/well-posedness test.
"""
from pathlib import Path
import sys
import time

import sympy as sp

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import T_singular_residual as T
from AH_lambda5_gram import stress
X, ALG = T.X, T.ALG
MEMBER = (1, 0, 0, 0)
CHECKS = []
# endpoint, sandwich5, wedge5, sandwich9, wedge9, trace4, trace8
DATA = {'w': (-sp.Rational(3, 704), 52, 48, sp.Rational(152, 3), sp.Rational(16, 3), 272, -sp.Rational(512, 3)),
        'cG': (-sp.Rational(1, 160), -sp.Rational(44, 3), sp.Rational(8, 3), -sp.Rational(40, 3), -sp.Rational(56, 3), sp.Rational(256, 3), -sp.Rational(176, 3))}


def check(label, condition):
    good = bool(condition)
    CHECKS.append(good)
    print(f'[{len(CHECKS):02d}] {label}: {good}', flush=True)
    assert good, label


def grade(poly, degree):
    return {m: sp.expand(c) for m, c in poly.items() if m.bit_count() == degree and sp.expand(c) != 0}


def form(indices, value=1):
    return {sum(1 << j for j in indices): value}


def lift(omega, mode):
    field = [grade(ALG.mul({1 << a: 1}, omega), 6 if mode == 'w' else 4) for a in range(X.N)]
    return field if mode == 'w' else ALG.gamma_right(field)


def real_grade(value, degree):
    return [{m: sp.expand(c/sp.I) for m, c in grade(item, degree).items()} for item in value]


def traces(field):
    alt, ctr = {}, {}
    for a, item in enumerate(field):
        product = ALG.mul({1 << a: 1}, item)
        alt, ctr = X.plus(alt, grade(product, 6)), X.plus(ctr, grade(product, 4))
    return alt, ctr


def K(field, inverse=False):
    alt, ctr = traces(field)
    c1, c2 = (-sp.Rational(1, 5), -sp.Rational(1, 9)) if inverse else (-1, -1)
    return [X.plus(X.plus(item, X.scaled(grade(ALG.mul({1 << a: 1}, alt), 5), c1)),
                   X.scaled(grade(ALG.mul({1 << a: 1}, ctr), 5), c2)) for a, item in enumerate(field)]


def predicted_Q(omega, mode, degree):
    square = ALG.mul(omega, omega)
    b, a = DATA[mode][1:3] if degree == 5 else DATA[mode][3:5]
    return [X.plus(X.scaled(grade(ALG.mul(ALG.mul(omega, {1 << j: 1}), omega), degree), b),
                   X.scaled(grade(ALG.mul({1 << j: 1}, square), degree), a)) for j in range(X.N)]


def contraction(field):
    result = {}
    for a, item in enumerate(field):
        result = X.plus(result, ALG.mul({1 << a: 1}, item))
    return {m: sp.expand(c) for m, c in result.items() if sp.expand(c) != 0}


def forced_jet(omega, mode):
    Q = T.quadratic(lift(omega, mode), MEMBER)
    source = real_grade(Q, 5) if mode == 'w' else ALG.gamma_right(real_grade(Q, 9))
    return T.scale(K(source, inverse=True), sp.Rational(1, 16) if mode == 'w' else -sp.Rational(1, 8))


def norm_gradient(omega, jet):
    return [sp.expand(2*sum(omega.get(mask, 0)*value for mask, value in item.items())) for item in jet]


def norm_current(omega, mode):
    cube = ALG.mul(ALG.mul(omega, omega), omega)
    if mode == 'w':
        value, factor = grade(cube, 1), sp.Rational(157, 18)
    else:
        value = grade(ALG.mul(grade(cube, 13), {(1 << X.N)-1: 1}), 1)
        factor = sp.Rational(8, 5)
    return [sp.expand(factor*value.get(1 << a, 0)) for a in range(X.N)]


def prolongation(omega, mode, jet):
    degree = 9 if mode == 'w' else 5
    field = lift(omega, mode)
    return [real_grade(T.mix(field, lift(item, mode), MEMBER), degree) for item in jet]


def higher_residual(omega, mode, jet):
    Q = T.quadratic(lift(omega, mode), MEMBER)
    derivative = [dict() for _ in range(X.N)]
    for a, item in enumerate(jet):
        derivative = X.fplus(derivative, T.derivative(lift(item, mode), a, MEMBER))
    residual = X.fplus(derivative, Q)
    return [grade(item, 5) for item in residual], [grade(item, 9) for item in residual]


def main():
    start = time.time()
    print('AL: opposite Lambda5 constant endpoint systems; real compact WA, all jets retained.', flush=True)
    simple = form(tuple(range(5)))
    alpha, beta = sp.symbols('A B', real=True)
    kappa = sp.symbols('kappa', real=True, nonzero=True)
    for mode in DATA:
        for direction in (0, 5):
            jet = [simple if a == direction else {} for a in range(X.N)]
            actual = T.derivative(lift(simple, mode), direction, MEMBER)
            expected = T.scale(K(jet), -16*sp.I) if mode == 'w' else T.scale(ALG.gamma_right(K(jet)), -8*sp.I)
            check(f'{mode} jet orbit {direction}: complete derivative in all slots and grades',
                  X.normalized(actual) == X.normalized(expected)
                  and all(mask.bit_count() == (5 if mode == 'w' else 9) for item in actual for mask in item))
            check(f'{mode} jet orbit {direction}: exact full jet inverse', X.normalized(K(K(jet), inverse=True)) == X.normalized(jet))
        for overlap in range(6):
            I, J = tuple(range(5)), tuple(range(overlap))+tuple(range(5, 10-overlap))
            omega = X.plus(form(I, alpha), form(J, beta))
            field = lift(omega, mode)
            Q = T.quadratic(field, MEMBER)
            mass = X.fplus(T.scale(field, kappa), T.mix(field, T.invariant(kappa*DATA[mode][0], 0), MEMBER))
            check(f'{mode} pair overlap{overlap}: compact exact endpoint kernel',
                  T.antihermitian(field) and not any(X.normalized(mass))
                  and all(mask.bit_count() in (1, 5, 9) for item in Q for mask in item))
            check(f'{mode} pair overlap{overlap}: both universal quadratic all-slot formulas',
                  all(X.normalized(real_grade(Q, degree)) == X.normalized(predicted_Q(omega, mode, degree)) for degree in (5, 9)))
            square = ALG.mul(omega, omega)
            check(f'{mode} pair overlap{overlap}: both complete gamma contractions',
                  all(contraction(real_grade(Q, degree)) == X.scaled(grade(square, degree-1), DATA[mode][5+index])
                      for index, degree in enumerate((5, 9))))
        print(f'  {mode} coefficients:', DATA[mode], flush=True)
    cartan = [dict() for _ in range(X.N)]
    cartan[0], cartan[1] = form((1, 2, 3, 4, 5)), form((0, 2, 3, 4, 5))
    check('Cartan jet has zero exterior traces but survives K unchanged',
          not any(traces(cartan)) and X.normalized(K(cartan)) == X.normalized(cartan) and any(cartan))
    pairs = []
    for mode, overlaps in (('w', (2, 3)), ('cG', (0, 1))):
        for overlap in overlaps:
            omega = X.plus(simple, form(tuple(range(overlap))+tuple(range(5, 10-overlap))))
            jet = forced_jet(omega, mode)
            residual = higher_residual(omega, mode, jet)
            check(f'{mode} compatible pair{overlap}: complete higher residual and first constraint prolongation',
                  not any(any(field) for field in residual)
                  and not any(any(field) for field in prolongation(omega, mode, jet)))
            pairs.append(omega)
    ordinary = X.plus(X.plus(form((0, 1, 2, 3, 4)), form((1, 2, 5, 6, 7))), form((3, 4, 5, 6, 7)))
    axial = X.plus(X.plus(form((0, 1, 2, 3, 4)), form((0, 5, 6, 7, 8))), form((0, 9, 10, 11, 12)))
    full = X.plus(axial, form((1, 5, 8, 11, 13), -2))
    for label, omega in [('ordinary triple', ordinary), ('axial triple', axial), ('general full support', full)]:
        for mode in DATA:
            check(f'{mode} {label}: complete cubic norm-current formula',
                  norm_gradient(omega, forced_jet(omega, mode)) == norm_current(omega, mode))
    for mode, omega, direction, value in (('w', ordinary, 0, -sp.Rational(157, 9)), ('cG', axial, 13, sp.Rational(48, 5))):
        jet = forced_jet(omega, mode)
        residual = higher_residual(omega, mode, jet)
        gradient = norm_gradient(omega, jet)
        check(f'{mode} nonconstant-norm real-form jet meets both higher equations and first algebraic prolongation',
              not any(any(field) for field in residual)
              and not any(any(field) for field in prolongation(omega, mode, jet))
              and gradient == [value if a == direction else 0 for a in range(X.N)])
        endpoint = DATA[mode][0]
        b = endpoint*(1+312*endpoint)  # kappa=1 for this algebraic control
        norm, gram = stress(omega)
        q, ell = (sp.Rational(112, 3), sp.Rational(320, 3)) if mode == 'w' else (40, -32)
        synthetic_G = -(b+q*norm)*sp.eye(X.N)-ell*gram
        Q = T.quadratic(lift(omega, mode), MEMBER)
        derivative = [dict() for _ in range(X.N)]
        for a, item in enumerate(jet):
            derivative = X.fplus(derivative, T.derivative(lift(item, mode), a, MEMBER))
        curvature = [{1 << column: sp.I*synthetic_G[row, column]
                      for column in range(X.N) if synthetic_G[row, column] != 0} for row in range(X.N)]
        to = T.invariant(endpoint, 0)
        baseline = X.fplus(curvature, X.fplus(to, T.quadratic(to, MEMBER)))
        check(f'{mode} complete odd residual agrees with the channel reduction for synthetic algebraic G',
              not any(X.normalized(X.fplus(baseline, X.fplus(derivative, Q)))))
        print(f'  {mode} countercontrol: grad N={value}*e{direction}; no observer metric or solution constructed.', flush=True)
    print(f'{len(CHECKS)} exact controls passed; runtime {time.time()-start:.3f} s.', flush=True)
    print('Exact channel reduction and conditional uniqueness; curvature integrability and actual positive-Z1 realization remain open.', flush=True)


if __name__ == '__main__':
    main()
