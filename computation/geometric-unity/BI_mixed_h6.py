"""Exact obstruction for a rational mixed harmonic-degree-three/five H6 family.

The polynomial family is fixed by the alternating coordinate-permutation
subgroup A5 of SO5 and has zero infinitesimal stabilizer. Exact quarter
rotations reconstruct its tensors from BF's two zonal seeds. Complete
diagonal moment polynomials, including mixed terms, have a positive
difference. This helper does not classify the full mixed module.

All inputs and imports are adjacent public files. No imported main is run;
running writes stdout only.
"""
from fractions import Fraction as F
from functools import reduce
from itertools import combinations
from math import prod
from pathlib import Path
from time import perf_counter
import json
import sys
import sympy as s

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import BC_affine_hooks as BC
import BF_h6_casimir as BF
import BC_quadratic_xcheck as X

N = 14
VARIABLES = s.symbols('x0:5')
PAIRS = tuple(combinations(range(5), 2))


def poly_action(value, i, j):
    """Differential of F(x) -> F(k^-1 x), for J_ij=E_ij-E_ji."""
    return s.expand(-VARIABLES[j] * s.diff(value, VARIABLES[i])
                    + VARIABLES[i] * s.diff(value, VARIABLES[j]))


def tangent_action(value, i):
    return {key: 2 * coefficient
            for key, coefficient in BC.action(6 + i, value).items()}


def positive_casimir(value):
    result = BF.casimir(value)
    for which in range(6):
        for key, coefficient in BC.action(which, BC.action(which, value)).items():
            BC.plus(result, key, -coefficient)
    return result


def quarter_turn(seed, polynomial, degree, axis):
    """Exact pi/2 rotation on the odd zonal cyclic representation."""
    tensors = [seed]
    polynomials = [polynomial]
    for _ in range(degree + 1):
        tensors.append(tangent_action(tensors[-1], axis))
        polynomials.append(poly_action(polynomials[-1], axis, 4))
    if degree == 3:
        assert not BF.combine((1, tensors[4]), (10, tensors[2]), (9, seed))
        result = BF.combine((F(7, 6), tensors[1]), (F(1, 6), tensors[3]))
        actual_poly = s.expand((7 * polynomials[1] + polynomials[3]) / 6)
    else:
        assert degree == 5
        assert not BF.combine((1, tensors[6]), (35, tensors[4]),
                              (259, tensors[2]), (225, seed))
        result = BF.combine((F(149, 120), tensors[1]),
                            (F(1, 4), tensors[3]), (F(1, 120), tensors[5]))
        actual_poly = s.expand(
            (149 * polynomials[1] + 30 * polynomials[3] + polynomials[5]) / 120)
    expected = polynomial.xreplace({VARIABLES[4]: VARIABLES[axis],
                                    VARIABLES[axis]: -VARIABLES[4]})
    assert s.expand(actual_poly - expected) == 0
    return result


def main():
    started = perf_counter()
    radius = sum(v * v for v in VARIABLES)
    t = VARIABLES[4]
    zonal3 = t**3 - s.Rational(3, 7) * t * radius
    zonal5 = (t**5 - s.Rational(10, 11) * t**3 * radius
              + s.Rational(5, 33) * t * radius**2)
    family3 = sum(v**3 - s.Rational(3, 7) * v * radius for v in VARIABLES)
    family5 = sum(v**5 - s.Rational(10, 11) * v**3 * radius
                  + s.Rational(5, 33) * v * radius**2 for v in VARIABLES)
    assert all(s.expand(sum(s.diff(value, v, 2) for v in VARIABLES)) == 0
               for value in (zonal3, zonal5, family3, family5))
    ranks = []
    for value in (family3, family5):
        actions = [s.Poly(poly_action(value, i, j), *VARIABLES) for i, j in PAIRS]
        monomials = sorted(set().union(*(set(p.monoms()) for p in actions)))
        matrix = s.Matrix([[p.coeff_monomial(m) for p in actions] for m in monomials])
        rank = matrix.rank()
        assert rank == 10
        ranks.append(rank)
        cycle = value.xreplace(dict(zip(VARIABLES, VARIABLES[1:] + VARIABLES[:1])))
        transposition = value.xreplace(
            {VARIABLES[0]: VARIABLES[1], VARIABLES[1]: VARIABLES[0]})
        assert s.expand(value - cycle) == s.expand(value - transposition) == 0

    data = json.loads(
        (HERE / 'BC_h6_closed_generator.json').read_text(encoding='utf-8'))
    assert data['degree'] == 6
    assert tuple(map(F, data['metric_gram'])) == tuple(BC.METRIC)
    original = {(a, mask): F(value) for a, mask, value in data['entries']}
    assert len(original) == len(data['entries']) == 1358
    first = BF.casimir(original)
    seed3 = BF.combine((F(40, 22), original), (-F(1, 22), first))
    seed5 = BF.combine((F(1, 22), first), (-F(18, 22), original))
    mapped = []
    for degree, polynomial, seed in ((3, zonal3, seed3), (5, zonal5, seed5)):
        rotations = [seed] + [quarter_turn(seed, polynomial, degree, i)
                              for i in range(4)]
        value = BF.combine(*((1, item) for item in rotations))
        assert BC.hook_traces(value) == ({}, {})
        assert not BF.combine((1, positive_casimir(value)),
                              (-degree * (degree + 3), value))
        assert BF.inner(value, value) > 0
        mapped.append(value)
    assert BF.inner(*mapped) == 0
    print('Flat positive Z1: lambda=1/5, cH=2, cV=1; real mixed H6 point values.')
    print('Declared family: B6=A*iota3(H3(sum xi^3))+B*iota5(H5(sum xi^5)).')
    print('Both harmonic polynomials are permutation invariant; A5 acts inside SO5.')
    print('Polynomial Lie-action ranks:', tuple(ranks), '; both infinitesimal stabilizers are zero.')
    print('Every nonzero family value is outside all conjugate-SO4-fixed spaces.')
    print('Quarter turns use (7D+D^3)/6 and (149D+30D^3+D^5)/120; full tensor annihilators checked.')
    print('Mapped supports:', tuple(map(len, mapped)), '; both hooks and full Casimir eigenvalues 18,40 checked.')

    gram = tuple(BC.METRIC)
    determinant = prod(gram)
    full = (1 << N) - 1
    algebra = BC.Clifford(gram)
    assert algebra.mul({full: F(1)}, {full: F(1)}) == {0: -1 / determinant}

    def converted(value):
        return [algebra.mul({mask: c for (b, mask), c in value.items() if a == b},
                            {full: F(1)}) for a in range(N)]

    def moments(value):
        actual = converted(value)
        assert all(mask.bit_count() == 8 for item in actual for mask in item)
        return tuple(determinant * BC.quadratic_coefficient(actual, gram, a, 1 << a)
                     / gram[a] for a in range(N))

    first_q, second_q = map(moments, mapped)
    sum_q = moments(BF.combine((1, mapped[0]), (1, mapped[1])))
    A, B = s.symbols('A B', real=True)
    values = tuple(s.expand(s.Rational(a) * A * A + s.Rational(c - a - b) * A * B
                            + s.Rational(b) * B * B)
                   for a, b, c in zip(first_q, second_q, sum_q))
    norm = (s.Rational(BF.inner(mapped[0], mapped[0])) * A * A
            + s.Rational(BF.inner(mapped[1], mapped[1])) * B * B)
    assert s.expand(sum(values) + s.Rational(52, 3) * norm) == 0
    assert not any(value.atoms(s.Float) for value in values)
    print('Actual compact field i*B6*Gamma; coordinate-volume quadratic multiplied by det Gram.')
    print('Positive norm polynomial:', norm)
    for a, value in enumerate(values):
        print('GRADE1 DIAGONAL', a, s.factor(value))
    print('Exact full diagonal trace=-(52/3)*norm; every A*B coefficient retained.')

    independent = X.BladeAlgebra(gram)
    control = converted(BF.combine((2, mapped[0]), (-3, mapped[1])))
    for axis in (0, 4, 13):
        observed = (determinant * independent.quadratic_coefficient(control, axis, 1 << axis)
                    / gram[axis])
        assert observed == values[axis].subs({A: 2, B: -3})
    print('Three independent nested Clifford targets agree at amplitudes (2,-3).')

    differences = [s.expand(value - values[0]) for value in values[1:]]
    infinity = [value.subs({A: 0, B: 1}) for value in differences]
    z = s.Symbol('z')
    polynomials = [s.Poly(value.subs({A: 1, B: z}), z, domain=s.QQ)
                   for value in differences if value != 0]
    common = reduce(s.gcd, polynomials)
    positive = s.expand(values[13] - values[0])
    polynomial = 892928 * A * A + 163520 * A * B + 15975 * B * B
    assert s.expand(positive - s.Rational(375, 495616) * polynomial) == 0
    positive_matrix = s.Matrix([[892928, 81760], [81760, 15975]])
    assert positive_matrix[0, 0] > 0 and positive_matrix.det() == 7579827200 > 0
    squares = (892928 * (A + s.Rational(2555, 27904) * B)**2
               + s.Rational(7402175, 872) * B * B)
    assert s.expand(polynomial - squares) == 0
    assert common.degree() == 2 and common.discriminant() < 0 and any(infinity)
    intervals = s.polys.polytools.intervals(common)
    assert not intervals
    assert F(3, 104) + 312 * F(3, 104)**2 == F(15, 52)
    print('Projective chart A=0,B!=0 differences:', infinity)
    print('For A!=0, z=B/A, exact common diagonal gcd:', common.as_expr())
    print('Exact real-root intervals:', intervals)
    print('Diagonal13-diagonal0=(375/495616)(892928A^2+163520AB+15975B^2).')
    print('Positive coefficient-matrix determinant:', positive_matrix.det())
    print('Verified squares: 892928(A+2555B/27904)^2+(7402175/872)B^2.')
    print('Every nonzero real amplitude pair fails scalar grade-one moments.')
    print('Zero fails the positive pure-branch baseline; no field is constructed.')
    print('No off-diagonal or higher-grade check is needed for this family obstruction.')
    print('This helper tests only the declared mixed family, not complete modules or general H6 fields.')
    print('No imported main, invariant solve, full gradient or source-pin change was run.')
    print(f'Runtime {perf_counter() - started:.3f} seconds.')


if __name__ == '__main__':
    main()
