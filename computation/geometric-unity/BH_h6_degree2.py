"""Exact grade-one obstruction on the degree-two H6 point-value module.

The complete diagonal section has four amplitudes.  All ten off-diagonal
stress polynomials and their cross terms are computed exactly.  SO5 spectral
diagonalization and the hand lemma in BH-h6-degree2.md supply the universal
pointwise result.  No closed-field lift or full curvature-kernel claim is
made.  Imports use adjacent public libraries; no earlier main is run.
"""
from fractions import Fraction as F
from itertools import combinations
from math import prod
from pathlib import Path
from time import perf_counter
import hashlib
import json
import sys
import sympy as s

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import BC_affine_hooks as BC
import BF_h6_casimir as BF
import BC_quadratic_xcheck as X


def positive_casimir(tensor):
    """All ten SO5 terms: six stabilizer and four doubled Nomizu actions."""
    result = BF.casimir(tensor)
    for j in range(6):
        for key, value in BC.action(j, BC.action(j, tensor)).items():
            BC.plus(result, key, -value)
    return result


def main():
    start = perf_counter()
    path = HERE/'BG_h6_curvature_kernel.json'
    receipt = json.loads(path.read_text(encoding='utf-8'))
    assert receipt['degree'] == 6
    assert tuple(map(F, receipt['metric_gram'])) == tuple(BC.METRIC)
    assert receipt['rank'] == 20 and receipt['kernel_dimension'] == 3
    entries = receipt['kernel_tensors'][2]['entries']
    seed = {(a, mask): F(c) for a, mask, c in entries}
    assert len(seed) == len(entries) == 878
    assert all(0 <= a < BC.N and mask.bit_count() == 6 for a, mask in seed)
    assert BC.hook_traces(seed) == ({}, {})
    assert all(not BC.action(j, seed) for j in range(6))
    assert not BF.combine((1, positive_casimir(seed)), (-10, seed))

    # A_i=rho(J_i4)/2.  The zonal degree-two seed satisfies A_i^3 v=-A_i v;
    # hence its actual pi/2 rotation is v+(rho(J_i4)^2 v)/2=v+2A_i^2 v.
    zonal = []
    for i in range(4):
        first = BC.action(6+i, seed)
        second = BC.action(6+i, first)
        third = BC.action(6+i, second)
        assert not BF.combine((1, third), (1, first))
        zonal.append(BF.combine((1, seed), (2, second)))
    zonal.append(seed)
    assert not BF.combine(*((1, value) for value in zonal))
    for value in zonal:
        assert BC.hook_traces(value) == ({}, {})
        assert not BF.combine((1, positive_casimir(value)), (-10, value))
    basis = zonal[:4]
    gram = s.Matrix(4, 4, lambda i, j: s.Rational(BF.inner(basis[i], basis[j])))
    assert gram.det() > 0
    assert all(gram[i, j] == (s.Rational(18225, 16) if i == j else -s.Rational(18225, 64))
               for i in range(4) for j in range(4))

    metric = tuple(BC.METRIC)
    determinant = prod(metric)
    algebra = BC.Clifford(metric)
    volume = {(1 << BC.N)-1: F(1)}

    def converted(tensor):
        result = [algebra.mul({mask: c for (b, mask), c in tensor.items() if b == a}, volume)
                  for a in range(BC.N)]
        assert all(mask.bit_count() == 8 for item in result for mask in item)
        return result

    def moments(tensor):
        value = converted(tensor)
        return [determinant*BC.quadratic_coefficient(value, metric, a, 1 << a)/metric[a]
                for a in range(10)]

    # Four pure amplitudes and all six mixed pairs determine every coefficient
    # of these exact homogeneous quadratics; no numerical polynomial fit.
    diagonal = [moments(tensor) for tensor in basis]
    matrices = [s.zeros(4) for _ in range(10)]
    for i in range(4):
        for axis in range(10):
            matrices[axis][i, i] = s.Rational(diagonal[i][axis])
    for i, j in combinations(range(4), 2):
        mixed = moments(BF.combine((1, basis[i]), (1, basis[j])))
        for axis in range(10):
            matrices[axis][i, j] = matrices[axis][j, i] = s.Rational(
                mixed[axis]-diagonal[i][axis]-diagonal[j][axis])/2
    z = s.symbols('z0:4', real=True)
    vector = s.Matrix(z)
    stresses = [s.expand((vector.T*m*vector)[0]) for m in matrices]
    mean = sum(z)/5
    eigenvalues = [value-mean for value in z]+[-mean]
    norm = sum(value*value for value in eigenvalues)
    aa, bb, cc = s.symbols('a b c')
    pairs = [(i, 4) for i in range(4)]+list(combinations(range(4), 2))
    equations = []
    for (i, j), stress in zip(pairs, stresses):
        expression = s.Poly(stress-aa*(eigenvalues[i]**2+eigenvalues[j]**2)
                            -bb*eigenvalues[i]*eigenvalues[j]-cc*norm, *z)
        equations.extend(expression.coeffs())
    solutions = s.solve(equations, (aa, bb, cc), dict=True)
    assert len(solutions) == 1
    coefficients = solutions[0]
    assert coefficients == {aa: -s.Rational(2375, 32), bb: s.Rational(1125, 16),
                            cc: -s.Rational(3625, 2)}
    assert all(expression.subs(coefficients) == 0 for expression in equations)
    assert coefficients[bb] > 0 and coefficients[bb]-3*coefficients[aa] > 0
    assert not any(value.atoms(s.Float) for value in stresses)

    # Exact algebraic steps of the hand spectrum lemma.  The proof's case
    # split and SO5 diagonalization are not replaced by finite test spectra.
    si, sj, sk, r = s.symbols('si sj sk r', real=True)

    def pair_polynomial(x, y):
        return aa*(x*x+y*y)+bb*x*y

    assert s.expand(pair_polynomial(si, sj)-pair_polynomial(si, sk)
                    -(sj-sk)*(aa*(sj+sk)+bb*si)) == 0
    exceptional = s.expand(pair_polynomial(-4*r, r)-pair_polynomial(r, r))
    assert s.expand(exceptional+5*(bb-3*aa)*r*r) == 0
    assert exceptional.subs(coefficients) == -s.Rational(46875, 32)*r*r

    # Different sorted-generator multiplication and nested-gradient code
    # checks all ten targets on a generic exact mixture.
    weights = (1, -2, 3, -4)
    generic = converted(BF.combine(*zip(weights, basis)))
    independent = X.BladeAlgebra(metric)
    for axis in range(10):
        observed = determinant*independent.quadratic_coefficient(generic, axis, 1 << axis)/metric[axis]
        assert observed == stresses[axis].subs(dict(zip(z, weights)))

    print('BH: flat positive Z1, lambda=1/5,cH=2,cV=1; real degree-two H6 point-value module.')
    print('Actual field i*B6*Gamma; B6*E quadratic multiplied by det Gram, with no extra sign.')
    print('Seed: all six SO4 generators and all ten positive SO5 Casimir terms give eigenvalue10.')
    print('Four exact quarter rotations checked; five zonals sum to zero and remain hook/Casimir10.')
    print('Independent diagonal section dimension4; Gram=', gram.tolist())
    print('All four amplitudes and all six cross terms retained on each of ten off-diagonal targets.')
    for pair, stress in zip(pairs, stresses):
        print('PAIR', pair, 'POLYNOMIAL', s.factor(stress))
    print('For trace-zero S=diag(s0,...,s4): qij=a(si^2+sj^2)+b si sj+c tr(S^2).')
    print('Exact coefficients:', coefficients)
    print('Nonzero lemma coefficients b,b-3a:', coefficients[bb], coefficients[bb]-3*coefficients[aa])
    print('Pair-difference factorization and exceptional spectrum(-4r,r,r,r,r) identity verified.')
    print('All ten independent nested-Clifford targets agree at (1,-2,3,-4).')
    print('Hand lemma and SO5 diagonalization: scalar Q1 forces the entire module value to zero.')
    print('Zero fails pure H6Gamma baseline L=15/(2cV)+(15/52)kappa^2>0, either kappa sign.')
    print('Positive scales preserve scalarity: horizontal dilation and common cV^(-7) factor.')
    print('No closed lift, full curvature-kernel classification or mixed-module no-go inferred.')
    print('Input receipt SHA-256:', hashlib.sha256(path.read_bytes()).hexdigest())
    print('No previous helper main or recorded output was run/imported.')
    print(f'Runtime {perf_counter()-start:.3f} seconds.')


if __name__ == '__main__':
    main()
