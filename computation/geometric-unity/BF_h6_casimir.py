"""Exact pointwise SO5 Casimir certificate for a contained H6 curvature module.

The hand argument in BF-h6-curvature-module.md converts the two spectral
components into a 121-dimensional SO5 orbit span in the curvature kernel.
This is not a full kernel classification or a field-equation solution.
Running writes no files and does not run BC's main.
"""
from fractions import Fraction as F
from math import prod
from pathlib import Path
from time import perf_counter
import json
import sys
import sympy as s

sys.path.insert(0, str(Path(__file__).resolve().parent))
import BC_affine_hooks as BC


def inner(a, b):
    return sum(value * b.get(key, 0) / (BC.METRIC[key[0]] * prod(BC.METRIC[j] for j in BC.bit_indices(key[1])))
               for key, value in a.items())


def casimir(tensor):
    result = {}
    for i in range(4):
        for key, value in BC.action(6+i, BC.action(6+i, tensor)).items():
            BC.plus(result, key, -4*value)
    return result


def calibration():
    sym = []
    from itertools import combinations
    for i, j in combinations(range(4), 2):
        mat = s.zeros(4); mat[i, j] = mat[j, i] = 1; sym.append(mat)
    for k in range(1, 4):
        sym.append(s.diag(*([1]*k + [-k] + [0]*(3-k))))
    sym.append(s.eye(4))
    phi = []
    for a in range(14):
        value = s.zeros(5)
        if a < 4:
            value[a, 4] = value[4, a] = -1
        else:
            mat = sym[a-4]
            value[:4, :4] = -2*mat + 2*mat.trace()*s.eye(4)/5
            value[4, 4] = 2*mat.trace()/5
        phi.append(value)
    assert all((phi[a]*phi[b]).trace() == (BC.METRIC[a] if a == b else 0) for a in range(14) for b in range(14))
    for i in range(4):
        J = s.zeros(5); J[i, 4] = 1; J[4, i] = -1
        M = s.Matrix(14, 14, lambda a, b: -dict(BC.COV[6+i][a]).get(b, 0))
        for a in range(14):
            lifted = sum((M[d, a]*phi[d] for d in range(14)), s.zeros(5))
            assert 2*lifted == J*phi[a] - phi[a]*J
    for which, (i, j) in enumerate(combinations(range(4), 2)):
        J = s.zeros(5); J[i, j] = 1; J[j, i] = -1
        M = s.Matrix(14, 14, lambda a, b: -dict(BC.COV[which][a]).get(b, 0))
        for a in range(14):
            lifted = sum((M[d, a]*phi[d] for d in range(14)), s.zeros(5))
            assert lifted == J*phi[a] - phi[a]*J
    print('Phi is metric-isometric; all140 tangent columns calibrate the full SO5 action.')
    print('On SO4-fixed tensors, C=-4 sum(horizontal Nomizu tensor actions squared).')



def combine(*terms):
    result = {}
    for coefficient, tensor in terms:
        for key, value in tensor.items():
            BC.plus(result, key, coefficient*value)
    return result


def main():
    start = perf_counter()
    calibration()
    path = Path(__file__).with_name('BC_h6_closed_generator.json')
    data = json.loads(path.read_text(encoding='utf-8'))
    assert data['degree'] == 6
    tensor = {(a, mask): F(value) for a, mask, value in data['entries']}
    assert len(tensor) == len(data['entries']) == 1358
    assert tuple(map(F, data['metric_gram'])) == tuple(BC.METRIC)
    assert inner(tensor, tensor) > 0
    assert BC.hook_traces(tensor) == ({}, {})
    assert all(not BC.action(i, tensor) for i in range(6))
    assert not BC.exterior_covariant(tensor)
    first = casimir(tensor)
    second = casimir(first)
    for value in (first, second):
        assert BC.hook_traces(value) == ({}, {})
        assert all(not BC.action(i, value) for i in range(6))
    assert inner(first, tensor) == 4*sum(inner(BC.action(6+i, tensor), BC.action(6+i, tensor)) for i in range(4))
    vectors = (tensor, first, second)
    gram = s.Matrix(3, 3, lambda i, j: s.Rational(inner(vectors[i], vectors[j])))
    assert gram[:2, :2].det() > 0 and gram.rank() == 2
    assert not combine((1, second), (-58, first), (720, tensor))
    v18 = combine((F(40, 22), tensor), (-F(1, 22), first))
    v40 = combine((F(1, 22), first), (-F(18, 22), tensor))
    assert inner(v18, v18) > 0 and inner(v40, v40) > 0
    assert inner(v18, v40) == 0
    assert not combine((1, casimir(v18)), (-18, v18))
    assert not combine((1, casimir(v40)), (-40, v40))
    assert not combine((1, tensor), (-1, v18), (-1, v40))
    from math import comb
    assert 3*(3+3) == 18 and 5*(5+3) == 40
    dimensions = (comb(7, 4)-comb(5, 4), comb(9, 4)-comb(7, 4))
    assert dimensions == (30, 91) and sum(dimensions) == 121
    print('Explicit nonzero H6 generator: both hooks, SO4 invariance and entire covariant closure verified.')
    print('Positive tensor norm:', inner(tensor, tensor))
    print('Positive 2x2 Gram determinant:', gram[:2, :2].det())
    print('Gram ranks on (B,CB) and (B,CB,C^2B): 2,2.')
    print('Every component verifies C^2B-58CB+720B=0; minimal polynomial=(c-18)(c-40).')
    print('Nonzero orthogonal projection norms:', inner(v18, v18), inner(v40, v40))
    print('Every projection component verifies Cv18=18v18 and Cv40=40v40.')
    print('Sphere degrees3,5 have harmonic dimensions30,91; the hand intertwining proof gives orbit span121.')
    print('Contained curvature-kernel module only; no whole-kernel rank or nonlinear solution is asserted.')
    print('The pointwise polynomial is not asserted as a global differential-Casimir field identity.')
    print(f'Runtime {perf_counter()-start:.3f} seconds.')


if __name__ == '__main__':
    main()
