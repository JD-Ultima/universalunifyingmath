"""Exact four-contraction Gram reduction for the pure quintic H6 module.

The hand proof supplies the universal pairing count and symmetry. This
helper calibrates the four coefficients, verifies the Gram channels and
checks physical Clifford conventions. Adjacent libraries are imported
without their mains; running writes stdout only. The contraction target
has dimension 30, so no cubic-style rank contradiction is claimed.
"""
from fractions import Fraction as F
from itertools import combinations_with_replacement, product
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
METRIC = tuple(BC.METRIC)
DET = prod(METRIC)
ALG = BC.Clifford(METRIC)
FULL = (1 << 14)-1


def multinomial(indices):
    return factorial(len(indices))//prod(factorial(indices.count(i)) for i in range(5))


def tensor_blocks(polynomial):
    """T_ijabc in 35 symmetric triple blocks, retaining ordered-index weights."""
    p = s.Poly(s.expand(polynomial), *x)
    blocks = []
    for triple in combinations_with_replacement(range(5), 3):
        matrix = s.zeros(5)
        for i, j in product(range(5), repeat=2):
            indices = (i, j)+triple
            powers = tuple(indices.count(k) for k in range(5))
            matrix[i, j] = p.coeff_monomial(powers)/multinomial(indices)
        assert matrix == matrix.T and matrix.trace() == 0
        assert not matrix.atoms(s.Float)
        if matrix:
            blocks.append((s.Integer(multinomial(triple)), matrix))
    A = sum((weight*matrix*matrix for weight, matrix in blocks), s.zeros(5))
    norm = s.trace(A)
    assert norm > 0 and not norm.atoms(s.Float)
    return blocks, A, norm


def patterns(data, S):
    blocks, A, norm = data
    return s.Matrix([norm*s.trace(S*S), s.trace(A*S*S),
                     sum(weight*s.trace(matrix*S)**2 for weight, matrix in blocks),
                     sum(weight*s.trace(S*matrix*S*matrix) for weight, matrix in blocks)])


def converted(tensor):
    result = [ALG.mul({mask: value for (b, mask), value in tensor.items() if b == a}, {FULL: F(1)})
              for a in range(14)]
    assert all(mask.bit_count() == 8 for component in result for mask in component)
    return result


def coefficient(field, a, b):
    return s.Rational(DET*BC.quadratic_coefficient(field, METRIC, a, 1 << b))


def full_casimir(tensor):
    result = BF.casimir(tensor)
    for j in range(6):
        for key, value in BC.action(j, BC.action(j, tensor)).items():
            BC.plus(result, key, -value)
    return result


def check_module(tensor):
    assert BC.hook_traces(tensor) == ({}, {})
    assert not BF.combine((1, full_casimir(tensor)), (-40, tensor))


def main():
    start = perf_counter()
    receipt = json.loads((HERE/'BC_h6_closed_generator.json').read_text(encoding='utf-8'))
    assert receipt['degree'] == 6 and tuple(map(F, receipt['metric_gram'])) == METRIC
    original = {(a, mask): F(v) for a, mask, v in receipt['entries']}
    assert len(original) == len(receipt['entries']) == 1358
    seed = BF.combine((F(1, 22), BF.casimir(original)), (-F(18, 22), original))
    radius = sum(v*v for v in x)
    seed_poly = x[4]**5-s.Rational(10, 11)*x[4]**3*radius+s.Rational(5, 33)*x[4]*radius**2
    product_tensor = seed
    product_poly = seed_poly
    for i in range(4):
        product_tensor = {key: 2*value for key, value in BC.action(6+i, product_tensor).items()}
        product_poly = s.expand(-x[4]*s.diff(product_poly, x[i])+x[i]*s.diff(product_poly, x[4]))
    assert product_poly == 120*prod(x)
    product_tensor = {key: value/120 for key, value in product_tensor.items()}
    # Three rotation derivatives provide the independent fourth calibration
    # row, without constructing a basis of the 91-dimensional module.
    extra_tensor = seed
    extra_poly = seed_poly
    for i in range(3):
        extra_tensor = {key: 2*value for key, value in BC.action(6+i, extra_tensor).items()}
        extra_poly = s.expand(-x[4]*s.diff(extra_poly, x[i])+x[i]*s.diff(extra_poly, x[4]))
    assert s.expand(extra_poly-s.Rational(60, 11)*x[0]*x[1]*x[2]*(11*x[4]**2-radius)) == 0
    extra_poly *= s.Rational(11, 60)
    extra_tensor = {key: value*F(11, 60) for key, value in extra_tensor.items()}
    polynomials = (seed_poly, prod(x), extra_poly)
    values = (seed, product_tensor, extra_tensor)
    for value in values:
        check_module(value)
    data = [tensor_blocks(polynomial) for polynomial in polynomials]
    phi = [s.Matrix(matrix)/5 for matrix in BG.integer_phi()]
    assert all(s.trace(phi[a]*phi[b]) == (METRIC[a] if a == b else 0)
               for a in range(14) for b in range(14))
    actual = [converted(value) for value in values]
    base = ((0, 0), (0, 4), (0, 13))
    product_matrix = s.Matrix.vstack(*(patterns(data[k], phi[a]).T for k, a in base+((1, 4),)))
    assert product_matrix.det() == 0
    selected = base+((2, 0),)
    calibration = s.Matrix.vstack(*(patterns(data[k], phi[a]).T for k, a in selected))
    rhs = s.Matrix([coefficient(actual[k], a, a) for k, a in selected])
    assert calibration.det() == s.Rational(131072, 49413375)
    assert not calibration.atoms(s.Float) and not rhs.atoms(s.Float)
    coefficients = calibration.inv()*rhs
    assert coefficients == s.Matrix([-s.Rational(2311875, 32768), s.Rational(2244375, 32768),
                                     s.Rational(1873125, 32768), s.Rational(860625, 16384)])
    assert not coefficients.atoms(s.Float)
    alpha, beta, gamma, delta = coefficients
    channels = (14*alpha+s.Rational(14, 5)*beta+gamma+s.Rational(3, 10)*delta,
                s.Rational(27, 20)*beta+gamma+delta/10, gamma+delta, gamma-delta/2)
    assert channels == (-s.Rational(5923125, 8192), s.Rational(20300625, 131072),
                        s.Rational(3594375, 32768), s.Rational(253125, 8192))
    assert all(value != 0 for value in channels)
    print('BK: pure harmonic-quintic W40; normalized flat positive Z1, lambda=1/5,cH=2,cV=1.')
    print('Full weighted hooks and all ten Casimir terms checked on all three declared values.')
    print('Polynomial bridge D3 D2 D1 D0 z5=120*x0*x1*x2*x3*x4 checked exactly.')
    print('Zonal/product calibration minor is zero; xyz*(11*x4^2-|x|^2) supplies an independent row.')
    print('Calibration input/axis rows:', selected)
    print('Calibration determinant:', calibration.det())
    print('Four coefficients alpha,beta,gamma,delta:', tuple(coefficients))
    for name, channel in zip(('Scalar', 'H2', 'H4', 'Young22'), channels):
        print(name+' multiplier:', channel)
    ratios = []
    for k in range(3):
        blocks, A, norm = data[k]
        ratios.append(s.Rational(BF.inner(values[k], values[k]))/norm)
        for a in range(14):
            assert (coefficients.T*patterns(data[k], phi[a]))[0] == coefficient(actual[k], a, a)
        gram = s.Matrix(14, 14, lambda a, b: sum(
            weight*s.trace(matrix*phi[a])*s.trace(matrix*phi[b]) for weight, matrix in blocks))
        operator = s.diag(*(1/s.Rational(v) for v in METRIC))*gram
        assert s.trace(operator) == norm
        assert gram.rank() == (14, 10, 10)[k]
        assert not gram.atoms(s.Float)
        spectrum = operator.eigenvals()
        if k == 0:
            assert spectrum == {s.Rational(28, 1089): 4, s.Rational(4, 3267): 9,
                                s.Rational(140, 1089): 1}
            assert all(gram[a, b] == 0 for a in range(14) for b in range(14) if a != b)
            assert all(gram[a, a] > 0 for a in range(14))
        print('Calibration Gram operator', k, 'rank:', gram.rank(), '; spectrum:', spectrum)
    assert ratios[0] == ratios[1] == ratios[2] == s.Rational(1366875, 32768)
    assert channels[0] == -s.Rational(52, 3)*ratios[0]
    print('Common positive intertwiner norm ratio:', ratios[0])
    print('All 42 diagonal coefficient controls match the universal four-pattern formula.')
    derivative = {key: 2*value for key, value in BC.action(6, seed).items()}
    mix = BF.combine((1, seed), (1, derivative))
    check_module(mix)
    mixed_poly = s.expand(seed_poly-x[4]*s.diff(seed_poly, x[0])+x[0]*s.diff(seed_poly, x[4]))
    mixed_data = tensor_blocks(mixed_poly)
    mixed_field = converted(mix)
    a, b = 0, 13
    cross = (patterns(mixed_data, phi[a]+phi[b])-patterns(mixed_data, phi[a])-patterns(mixed_data, phi[b]))/2
    expected = (coefficients.T*cross)[0]
    forward = coefficient(mixed_field, a, b)
    backward = coefficient(mixed_field, b, a)
    assert forward == backward == expected == -s.Rational(568125, 11264)
    assert not forward.atoms(s.Float) and not backward.atoms(s.Float)
    independent = X.BladeAlgebra(METRIC)
    for i, j, target in ((0, 13, forward), (13, 0, backward)):
        assert s.Rational(DET*independent.quadratic_coefficient(mixed_field, i, 1 << j)) == target
    print('Ordered off-diagonal entries q_ab=q_ba:', forward, '; both independent nested targets agree.')
    print('Exact remaining condition: K_T=(||T||^2/14)*Id on Sym0(R5).')
    print('L_T maps 14 inputs to H3 of dimension 30: no cubic-style rank contradiction follows.')
    print('No scalar-Gram candidate, closed-field lift or universal quintic exclusion is inferred.')
    print('No earlier main, recorded output, invariant solve or full quadratic gradient was run/imported.')
    print(f'Runtime {perf_counter()-start:.3f} seconds.')


if __name__ == '__main__':
    main()
