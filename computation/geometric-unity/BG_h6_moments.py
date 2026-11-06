"""Exact grade-one obstruction on SO(4)-fixed H6 curvature values.

The adjacent curvature receipt supplies three pointwise tensors. This helper
checks all six quadratic amplitude coefficients, the compact phase, norm,
trace and positive-square obstruction. Completeness of the curvature kernel
is certified separately. No field reconstruction or full-module claim follows.

Imports adjacent public libraries without their mains; running writes stdout
only and uses no private source, saved output or external workspace path.
"""
from fractions import Fraction as F
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
FULL = (1 << N) - 1
GRAM = tuple(BC.METRIC)
DET = prod(GRAM)
ALG = BC.Clifford(GRAM)
BASIS_ORDER = (
    'H0', 'H1', 'H2', 'H3', 'O01', 'O02', 'O03',
    'O12', 'O13', 'O23', 'D1', 'D2', 'D3', 'I',
)


def converted(tensor):
    """Right-multiply B6 by the coordinate volume, giving real grade eight."""
    result = []
    for a in range(N):
        value = {mask: coefficient for (b, mask), coefficient in tensor.items()
                 if a == b}
        value = ALG.mul(value, {FULL: F(1)})
        assert all(mask.bit_count() == 8 for mask in value)
        result.append(value)
    return result


def moments(tensor):
    """Metric-relative Q1(i B6 Gamma)/i on H0, O01 and the scalar axis."""
    value = converted(tensor)
    return tuple(DET * BC.quadratic_coefficient(value, GRAM, a, 1 << a) / GRAM[a]
                 for a in (0, 4, 13))


def main():
    started = perf_counter()
    receipt = json.loads(
        (HERE / 'BG_h6_curvature_kernel.json').read_text(encoding='utf-8'))
    assert receipt['degree'] == 6
    assert tuple(map(F, receipt['metric_gram'])) == GRAM
    assert tuple(receipt['basis_order']) == BASIS_ORDER
    raw = []
    for entry in receipt['kernel_tensors']:
        value = {(a, mask): F(coefficient)
                 for a, mask, coefficient in entry['entries']}
        assert len(value) == len(entry['entries'])
        assert all(0 <= a < N and mask.bit_count() == 6 for a, mask in value)
        raw.append(value)
    assert len(raw) == 3

    data = json.loads(
        (HERE / 'BC_h6_closed_generator.json').read_text(encoding='utf-8'))
    assert data['degree'] == 6
    assert tuple(map(F, data['metric_gram'])) == GRAM
    original = {(a, mask): F(value) for a, mask, value in data['entries']}
    assert len(original) == len(data['entries']) == 1358
    first = BF.casimir(original)
    assert BF.combine((F(1, 48), raw[0]), (-F(5, 24), raw[1])) == original
    assert BF.combine((-F(17, 18), raw[0]), (-F(35, 9), raw[1])) == first
    v18 = BF.combine((F(40, 22), original), (-F(1, 22), first))
    v40 = BF.combine((F(1, 22), first), (-F(18, 22), original))
    basis = (v18, v40, raw[2])
    for eigenvalue, value in zip((18, 40, 10), basis):
        assert BC.hook_traces(value) == ({}, {})
        assert all(not BC.action(i, value) for i in range(6))
        assert not BF.combine((1, BF.casimir(value)), (-eigenvalue, value))
    gram = s.Matrix(3, 3, lambda i, j: s.Rational(BF.inner(basis[i], basis[j])))
    assert all(gram[i, i] > 0 for i in range(3))
    assert all(gram[i, j] == 0 for i in range(3) for j in range(3) if i != j)
    assert gram[2, 2] == s.Rational(18225, 16)
    assert ALG.mul({FULL: F(1)}, {FULL: F(1)}) == {0: -1 / DET}

    diagonal = [moments(value) for value in basis]
    matrices = [s.zeros(3) for _ in range(3)]
    for i in range(3):
        for j in range(3):
            matrices[j][i, i] = s.Rational(diagonal[i][j])
    for i, j in combinations(range(3), 2):
        mixed = moments(BF.combine((1, basis[i]), (1, basis[j])))
        for block in range(3):
            cross = mixed[block] - diagonal[i][block] - diagonal[j][block]
            matrices[block][i, j] = matrices[block][j, i] = s.Rational(cross) / 2

    A, B, Z = s.symbols('A B Z', real=True)
    u = s.Matrix([A, B, Z])
    polynomials = tuple(s.expand((u.T * matrix * u)[0]) for matrix in matrices)
    norm = s.expand((u.T * gram * u)[0])
    trace = 4 * polynomials[0] + 9 * polynomials[1] + polynomials[2]
    assert s.expand(trace) == -s.Rational(52, 3) * norm
    known = tuple(v.subs({A: 1, B: 1, Z: 0}) for v in polynomials)
    assert known == (F(-15287625, 4096), F(-10178375, 2048), F(-4793625, 2048))
    assert matrices[2][:2, :2] - matrices[0][:2, :2] == (
        s.Rational(1125, 45056) * s.Matrix([[54272, 480], [480, 505]]))

    print('Flat positive Z1: lambda=1/5, cH=2, cV=1; all tensor slots covariant.')
    print('SO(4)-fixed point values: B6=A*v18+B*v40+Z*v10; Casimir eigenvalues 18,40,10.')
    print('Actual field i*B6*Gamma. Gamma=sqrt(det Gram)*E; Q(i*B6*E)/i is multiplied by det Gram.')
    print('Receipt tensor supports:', tuple(map(len, raw)))
    print('Orthogonal positive norm Gram:', gram.tolist())
    print('Quadratic coefficient order: A^2, A*B, A*Z, B^2, B*Z, Z^2.')
    monomials = (A * A, A * B, A * Z, B * B, B * Z, Z * Z)
    for name, matrix, polynomial in zip(('H4', 'Sym0_9', 'I1'), matrices, polynomials):
        coefficients = tuple(s.Poly(polynomial, A, B, Z).coeff_monomial(m)
                             for m in monomials)
        print('MATRIX', name, matrix.tolist())
        print('COEFFICIENTS', name, coefficients)
        print('EIGENVALUE', name, s.factor(polynomial))
    print('B*=v18+v40 recovers BC moments:', known)
    print('Exact weighted trace=-(52/3)*N; N=', norm)

    independent = X.BladeAlgebra(GRAM)
    for weights in ((1, 1, 0), (1, 0, 1), (0, 1, 1), (2, -3, 1)):
        actual = converted(BF.combine(*zip(weights, basis)))
        for block, axis in enumerate((0, 4, 13)):
            observed = (DET * independent.quadratic_coefficient(actual, axis, 1 << axis)
                        / GRAM[axis])
            assert observed == polynomials[block].subs(dict(zip((A, B, Z), weights)))
        print('PASS independent nested Clifford moments at', weights)
    print('Independent comparison targets: 12.')

    differences = (matrices[1] - matrices[0], matrices[2] - matrices[0])
    assert tuple(matrix[2, 2] for matrix in matrices) == (
        -s.Rational(48375, 32), -s.Rational(11625, 8), -s.Rational(2475, 4))
    assert all(matrix[i, 2] == matrix[2, i] == 0
               for matrix in matrices for i in range(2))
    for name, matrix in zip(('S-H', 'I-H'), differences):
        principal = tuple(matrix.extract(indices, indices).det()
                          for size in range(1, 4)
                          for indices in combinations(range(3), size))
        print('DIFFERENCE', name, matrix.tolist())
        print('Principal minors, ordered by sizes 1,2,3:', principal)
        print('Difference polynomial:', s.factor((u.T * matrix * u)[0]))
    positive = s.expand((u.T * differences[1] * u)[0])
    squares = (
        s.Rational(1125, 45056) * (
            54272 * (A + s.Rational(15, 1696) * B)**2
            + s.Rational(26540, 53) * B * B)
        + s.Rational(28575, 32) * Z * Z)
    assert s.expand(positive - squares) == 0
    assert all(differences[1][:size, :size].det() > 0 for size in (1, 2, 3))
    assert not any(v.atoms(s.Float) for v in polynomials)
    assert F(3, 104) + 312 * F(3, 104)**2 == F(15, 52)
    print('All I-H leading principal minors are positive; coefficients contain no Float.')
    print('Verified I-H=(1125/45056)[54272(A+15B/1696)^2+(26540/53)B^2]+(28575/32)Z^2.')
    print('Scalar grade-one moments force A=B=Z=0 over the reals.')
    print('Zero fails the pure-branch baseline L=15/(2cV)+(15/52)kappa^2>0.')
    print('Pointwise obstruction on this fixed three-space only; no full-module or field-construction claim.')
    print('No imported main, invariant solve, curvature rank or full quadratic gradient was run.')
    print(f'Runtime {perf_counter() - started:.3f} seconds.')


if __name__ == '__main__':
    main()
