"""Exact ordered cubic/quintic grade-one cross calibration.

Cross means Q(B3+B5)-Q(B3)-Q(B5), with no factor one-half. The hand pairing
argument supplies universality. Adjacent libraries are imported without
their mains; running writes stdout only.
"""
from fractions import Fraction as F
from itertools import product
from math import factorial, prod
from pathlib import Path
from time import perf_counter
import json
import sys
import sympy as s

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import BC_affine_hooks as BC
import BF_h6_casimir as BF
import BG_h6_curvature as BG
import BC_quadratic_xcheck as X

x = s.symbols('x0:5')
G = tuple(BC.METRIC)
DET = prod(G)
ALG = BC.Clifford(G)
FULL = (1 << 14)-1
CROSS_COEFFICIENTS = (-s.Rational(2131875, 2816), s.Rational(2874375, 11264),
                      s.Rational(2874375, 11264))
NORM_RATIOS = (s.Rational(2197125, 352), s.Rational(1366875, 32768))


def polynomial_tensor(expression, degree):
    polynomial = s.Poly(s.expand(expression), *x)
    tensor = {}
    for indices in product(range(5), repeat=degree):
        powers = tuple(indices.count(i) for i in range(5))
        value = polynomial.coeff_monomial(powers)*prod(factorial(p) for p in powers)/factorial(degree)
        if value:
            tensor[indices] = F(value)
    assert all(sum(tensor.get((i, i)+tail, 0) for i in range(5)) == 0
               for tail in product(range(5), repeat=degree-2))
    return tensor


def derivative_polynomial(value):
    return s.expand(-x[4]*s.diff(value, x[0])+x[0]*s.diff(value, x[4]))


def derivative_tensor(value):
    return {key: 2*coefficient for key, coefficient in BC.action(6, value).items()}


def full_casimir(value):
    result = BF.casimir(value)
    for j in range(6):
        for key, coefficient in BC.action(j, BC.action(j, value)).items():
            BC.plus(result, key, -coefficient)
    return result


def cross_data(T, U):
    by_pair = {}
    for (i, j, b, c, d), value in U.items():
        by_pair.setdefault((i, j), []).append((b, c, d, value))
    Z = {}
    for (i, j, a), value in T.items():
        for b, c, d, other in by_pair.get((i, j), ()):
            BC.plus(Z, (a, b, c, d), value*other)
    C = s.Matrix(5, 5, lambda c, d: sum(s.Rational(Z.get((a, a, c, d), 0)) for a in range(5)))
    assert C == C.T and C.trace() == 0
    return C, Z


def patterns(data, S, R):
    C, Z = data
    second = sum(s.Rational(value)*S[a, b]*R[c, d] for (a, b, c, d), value in Z.items())
    third = sum(s.Rational(value)*R[a, b]*S[c, d] for (a, b, c, d), value in Z.items())
    return s.Matrix([s.trace(C*S*R), second, third])


def converted(value):
    result = [ALG.mul({mask: coefficient for (b, mask), coefficient in value.items() if b == a}, {FULL: F(1)})
              for a in range(14)]
    assert all(mask.bit_count() == 8 for item in result for mask in item)
    return result


def main():
    start = perf_counter()
    data = json.loads((HERE/'BC_h6_closed_generator.json').read_text(encoding='utf-8'))
    assert data['degree'] == 6 and tuple(map(F, data['metric_gram'])) == G
    original = {(a, mask): F(value) for a, mask, value in data['entries']}
    assert len(original) == len(data['entries']) == 1358
    casimir = BF.casimir(original)
    seeds = (BF.combine((F(40, 22), original), (-F(1, 22), casimir)),
             BF.combine((F(1, 22), casimir), (-F(18, 22), original)))
    radius = sum(v*v for v in x)
    polynomials = (x[4]**3-s.Rational(3, 7)*x[4]*radius,
                   x[4]**5-s.Rational(10, 11)*x[4]**3*radius+s.Rational(5, 33)*x[4]*radius**2)
    tensors = ((seeds[0], derivative_tensor(seeds[0])), (seeds[1], derivative_tensor(seeds[1])))
    symbolic = tuple(tuple(polynomial_tensor(p, degree) for p in (polynomial, derivative_polynomial(polynomial)))
                     for polynomial, degree in zip(polynomials, (3, 5)))
    for group, forms, eigenvalue, ratio in zip(tensors, symbolic, (18, 40), NORM_RATIOS):
        for value, tensor in zip(group, forms):
            assert BC.hook_traces(value) == ({}, {})
            assert not BF.combine((1, full_casimir(value)), (-eigenvalue, value))
            assert s.Rational(BF.inner(value, value))/sum(s.Rational(v*v) for v in tensor.values()) == ratio
    phi = [s.Matrix(mat)/5 for mat in BG.integer_phi()]
    pairs = ((0, 0), (0, 1), (1, 0), (1, 1))
    mixed = {(i, j): BF.combine((1, tensors[0][i]), (1, tensors[1][j])) for i, j in pairs}
    inputs = tuple(tuple(converted(value) for value in group) for group in tensors)
    combined = {pair: converted(value) for pair, value in mixed.items()}
    contractions = {(i, j): cross_data(symbolic[0][i], symbolic[1][j]) for i, j in pairs}
    cache = {}

    def actual(pair, a, b):
        key = pair, a, b
        if key not in cache:
            i, j = pair
            raw = (BC.quadratic_coefficient(combined[pair], G, a, 1 << b)
                   -BC.quadratic_coefficient(inputs[0][i], G, a, 1 << b)
                   -BC.quadratic_coefficient(inputs[1][j], G, a, 1 << b))
            cache[key] = s.Rational(DET*raw)
        return cache[key]

    rows = (((0, 0), 0, 0), ((0, 1), 0, 13), ((0, 1), 13, 0))
    calibration = s.Matrix.vstack(*(patterns(contractions[pair], phi[a], phi[b]).T for pair, a, b in rows))
    assert calibration.det() == -s.Rational(6144, 33275)
    rhs = s.Matrix([actual(pair, a, b) for pair, a, b in rows])
    coefficients = calibration.inv()*rhs
    assert coefficients == s.Matrix(CROSS_COEFFICIENTS)
    assert not calibration.atoms(s.Float) and not coefficients.atoms(s.Float)
    print('BL cross polarization: Q(B3+B5)-Q(B3)-Q(B5), no half factor.')
    print('BI/BK zonal normalizations, weighted hooks, full Casimirs and norm ratios verified.')
    print('Three-pattern calibration rows:', rows)
    print('Exact calibration matrix:', calibration.tolist())
    print('Exact calibration determinant:', calibration.det())
    print('Coefficients v1,v2,v3:', tuple(coefficients))
    tests = (((0, 0), 4, 4), ((0, 0), 13, 13),
             ((0, 1), 4, 1), ((0, 1), 1, 4),
             ((1, 0), 0, 13), ((1, 0), 13, 0),
             ((1, 0), 4, 1), ((1, 0), 1, 4),
             ((1, 1), 0, 0), ((1, 1), 1, 1), ((1, 1), 4, 4), ((1, 1), 13, 13))
    for pair, a, b in tests:
        predicted = (coefficients.T*patterns(contractions[pair], phi[a], phi[b]))[0]
        assert predicted == actual(pair, a, b)
        assert not predicted.atoms(s.Float)
    print('Twelve additional ordered coefficient controls match.')
    independent = X.BladeAlgebra(G)
    for pair, a, b in (((0, 1), 0, 13), ((0, 1), 13, 0), ((1, 0), 4, 1), ((1, 0), 1, 4)):
        i, j = pair
        raw = (independent.quadratic_coefficient(combined[pair], a, 1 << b)
               -independent.quadratic_coefficient(inputs[0][i], a, 1 << b)
               -independent.quadratic_coefficient(inputs[1][j], a, 1 << b))
        assert s.Rational(DET*raw) == actual(pair, a, b)
    print('Four ordered cross entries agree with the independent nested-gradient implementation.')
    assert coefficients[1] == coefficients[2]
    assert actual((0, 1), 0, 13) == actual((0, 1), 13, 0)
    print('Antisymmetric coefficient v2-v3:', s.factor((coefficients[1]-coefficients[2])/2))
    print('Observed ordered difference q_0,13-q_13,0 at(z3,D0z5):', actual((0, 1), 0, 13)-actual((0, 1), 13, 0))
    print('The hand pairing count supplies universality; scalarity and field equations require separate checks.')
    print(f'Runtime {perf_counter()-start:.3f} seconds.')


if __name__ == '__main__':
    main()
