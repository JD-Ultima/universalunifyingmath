"""Exact four-contraction certificate for the harmonic-cubic H6 stress.

The hand invariant-contraction proof supplies universality. This helper
calibrates all four coefficients, checks the trace/reshuffle channels and
validates additional diagonal and mixed stress entries. Adjacent public
libraries are imported without their mains; running writes stdout only.
"""
from fractions import Fraction as F
from itertools import product
from math import prod, factorial
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


def tensors(polynomial):
    p = s.Poly(s.expand(polynomial), *x)
    result = []
    for a in range(5):
        mat = s.zeros(5)
        for b, c in product(range(5), repeat=2):
            powers = tuple((a, b, c).count(i) for i in range(5))
            multiplicity = 6//prod(factorial(k) for k in powers)
            mat[b, c] = p.coeff_monomial(powers)/multiplicity
        result.append(mat)
    assert all(mat == mat.T and mat.trace() == 0 for mat in result)
    return result


def patterns(T, S):
    A = sum((mat*mat for mat in T), s.zeros(5))
    norm = s.trace(A)
    return s.Matrix([norm*s.trace(S*S), s.trace(A*S*S),
                     sum(s.trace(mat*S)**2 for mat in T),
                     sum(s.trace(S*mat*S*mat) for mat in T)])


def converted(tensor):
    result = [ALG.mul({mask: value for (b, mask), value in tensor.items() if b == a}, {FULL: F(1)})
              for a in range(14)]
    assert all(mask.bit_count() == 8 for item in result for mask in item)
    return result


def positive_casimir(tensor):
    """All ten SO5 terms, including the six stabilizer actions."""
    result = BF.casimir(tensor)
    for j in range(6):
        for key, value in BC.action(j, BC.action(j, tensor)).items():
            BC.plus(result, key, -value)
    return result


def check_hook_module(tensor):
    assert BC.hook_traces(tensor) == ({}, {})
    assert not BF.combine((1, positive_casimir(tensor)), (-18, tensor))


def covariant_coefficient(field, a, b):
    return s.Rational(DET*BC.quadratic_coefficient(field, METRIC, a, 1 << b))


def linear_channel_controls(T):
    """Exact component checks of the trace and symmetrization decomposition."""
    indices = tuple(product(range(5), repeat=4))
    eye = lambda a, b: int(a == b)
    P = {(a, b, c, d): s.Rational(eye(a, c)*eye(b, d)+eye(a, d)*eye(b, c), 2)
         -s.Rational(eye(a, b)*eye(c, d), 5) for a, b, c, d in indices}

    def cross(K):
        return s.Matrix(5, 5, lambda a, c: sum(K[a, b, c, b] for b in range(5)))

    def embedding(H):
        return {(a, b, c, d):
                (eye(a, c)*H[b, d]+eye(a, d)*H[b, c]+eye(b, c)*H[a, d]+eye(b, d)*H[a, c])/4
                -(eye(a, b)*H[c, d]+H[a, b]*eye(c, d))/5 for a, b, c, d in indices}

    def reshuffle(K):
        A = cross(K)
        norm = s.trace(A)
        return {(a, b, c, d): (K[a, c, b, d]+K[a, d, b, c])/2
                -(eye(a, b)*A[c, d]+A[a, b]*eye(c, d))/5
                +s.Rational(eye(a, b)*eye(c, d), 25)*norm for a, b, c, d in indices}

    H = s.diag(1, -1, 0, 0, 0)
    KH = embedding(H)
    assert cross(P) == s.Rational(14, 5)*s.eye(5)
    assert cross(KH) == s.Rational(27, 20)*H
    for K, factor in ((P, s.Rational(3, 10)), (KH, s.Rational(1, 10))):
        J = reshuffle(K)
        assert all(s.expand(J[key]-factor*K[key]) == 0 for key in indices)
    K = {(a, b, c, d): sum(mat[a, b]*mat[c, d] for mat in T) for a, b, c, d in indices}
    A = cross(K)
    norm = s.trace(A)
    H = s.Rational(20, 27)*(A-norm*s.eye(5)/5)
    KH = embedding(H)
    remainder = {key: K[key]-norm*P[key]/14-KH[key] for key in indices}
    H4 = {(a, b, c, d): (remainder[a, b, c, d]+remainder[a, c, b, d]
                          +remainder[a, d, b, c])/3 for a, b, c, d in indices}
    Y22 = {key: remainder[key]-H4[key] for key in indices}
    assert any(H4.values()) and any(Y22.values())
    for channel, factor in ((H4, s.Integer(1)), (Y22, -s.Rational(1, 2))):
        assert cross(channel) == s.zeros(5)
        assert all(sum(channel[a, b, c, c] for c in range(5)) == 0 for a, b in product(range(5), repeat=2))
        J = reshuffle(channel)
        assert all(s.expand(J[key]-factor*channel[key]) == 0 for key in indices)
    assert all(not value.atoms(s.Float) for tensor in (P, KH, H4, Y22) for value in tensor.values())
    print('625-component reshuffle identities: scalar 3/10, H2 1/10, H4 1, Young22 -1/2.')


def main():
    start = perf_counter()
    data = json.loads((HERE/'BC_h6_closed_generator.json').read_text(encoding='utf-8'))
    assert data['degree'] == 6 and tuple(map(F, data['metric_gram'])) == METRIC
    original = {(a, mask): F(v) for a, mask, v in data['entries']}
    assert len(original) == len(data['entries']) == 1358
    seed = BF.combine((F(40, 22), original), (-F(1, 22), BF.casimir(original)))
    # iota(t^3-3t|x|^2/7)=v18 and D_i=2*A_i.
    seed_poly = x[4]**3-s.Rational(3, 7)*x[4]*sum(v*v for v in x)
    xyz_tensor = seed
    xyz_poly = seed_poly
    for i in (0, 1, 2):
        xyz_tensor = {key: 2*value for key, value in BC.action(6+i, xyz_tensor).items()}
        xyz_poly = s.expand(-x[4]*s.diff(xyz_poly, x[i])+x[i]*s.diff(xyz_poly, x[4]))
    assert xyz_poly == 6*x[0]*x[1]*x[2]
    xyz_tensor = {key: value/6 for key, value in xyz_tensor.items()}
    cubics = (seed_poly, x[0]*x[1]*x[2])
    values = (seed, xyz_tensor)
    for value in values:
        check_hook_module(value)
    T = [tensors(p) for p in cubics]
    assert all(not entry.atoms(s.Float) for value in T for mat in value for entry in mat)
    phi = [s.Matrix(mat)/5 for mat in BG.integer_phi()]
    assert all(s.trace(phi[a]*phi[b]) == (METRIC[a] if a == b else 0)
               for a in range(14) for b in range(14))
    actual = [converted(value) for value in values]
    selected = ((0, 0), (0, 4), (0, 13), (1, 6))
    matrix = s.Matrix.vstack(*(patterns(T[k], phi[a]).T for k, a in selected))
    rhs = s.Matrix([covariant_coefficient(actual[k], a, a) for k, a in selected])
    assert matrix.det() == -s.Rational(4096, 60025)
    assert not matrix.atoms(s.Float) and not rhs.atoms(s.Float)
    coefficients = matrix.inv()*rhs
    assert coefficients == s.Matrix([-s.Rational(9124875, 968), s.Rational(4095875, 484),
                                     -s.Rational(268625, 484), s.Rational(1031625, 484)])
    alpha, beta, gamma, delta = coefficients
    assert not coefficients.atoms(s.Float)
    print('BI: harmonic-cubic H6 stress; normalized flat positive Z1, lambda=1/5,cH=2,cV=1.')
    print('Calibration values satisfy both weighted hooks and the full ten-term Casimir eigenvalue 18.')
    print('Four-contraction calibration determinant:', matrix.det())
    print('Coefficients alpha,beta,gamma,delta:', tuple(coefficients))
    channels = (14*alpha+s.Rational(14, 5)*beta+gamma+s.Rational(3, 10)*delta,
                s.Rational(27, 20)*beta+gamma+delta/10, gamma+delta, gamma-delta/2)
    assert channels == (-s.Rational(9520875, 88), s.Rational(21455875, 1936),
                        s.Rational(190750, 121), -s.Rational(142625, 88))
    assert all(value != 0 for value in channels)
    for name, value in zip(('Scalar', 'H2', 'H4', 'Young22'), channels):
        print(name+' channel:', value)
    linear_channel_controls(T[1])
    ratios = []
    for k in range(2):
        norm_T = sum(s.trace(mat*mat) for mat in T[k])
        ratios.append(s.Rational(BF.inner(values[k], values[k]))/norm_T)
        for a in range(14):
            expected = (coefficients.T*patterns(T[k], phi[a]))[0]
            assert expected == covariant_coefficient(actual[k], a, a)
    assert ratios[0] == ratios[1] == s.Rational(2197125, 352)
    assert channels[0] == -s.Rational(52, 3)*ratios[0]
    print('Both cubic intertwiners have identical positive norm ratio:', ratios[0])
    print('All 28 diagonal coefficient controls match the four-pattern formula.')
    # A nontrivial off-diagonal stress control on a new non-zonal mixture.
    derivative = {key: 2*value for key, value in BC.action(6, seed).items()}
    mix = BF.combine((1, seed), (1, derivative))
    check_hook_module(mix)
    mixed_poly = s.expand(seed_poly-x[4]*s.diff(seed_poly, x[0])+x[0]*s.diff(seed_poly, x[4]))
    mixed_T = tensors(mixed_poly)
    mixed_field = converted(mix)
    a, b = 0, 13
    cross_patterns = (patterns(mixed_T, phi[a]+phi[b])-patterns(mixed_T, phi[a])-patterns(mixed_T, phi[b]))/2
    expected = (coefficients.T*cross_patterns)[0]
    observed = (covariant_coefficient(mixed_field, a, b)+covariant_coefficient(mixed_field, b, a))/2
    assert expected == observed == -s.Rational(59625, 11)
    assert not expected.atoms(s.Float) and not observed.atoms(s.Float)
    independent = X.BladeAlgebra(METRIC)
    for i, j in ((0, 13), (13, 0)):
        assert s.Rational(DET*independent.quadratic_coefficient(mixed_field, i, 1 << j)) == covariant_coefficient(mixed_field, i, j)
    print('Nonzero symmetric off-diagonal mixture control:', observed)
    print('Both nested-gradient off-diagonal controls agree.')
    print('Universal metric-contraction count and Gram-rank obstruction are supplied by the hand proof.')
    print('No earlier main, invariant solve, full quadratic gradient or curvature rank was run.')
    print(f'Runtime {perf_counter()-start:.3f} seconds.')


if __name__ == '__main__':
    main()
