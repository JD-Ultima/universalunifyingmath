"""Reproduce the exact real partial root and a separating grade-nine coefficient."""

import argparse
import gzip
import hashlib
import json
from fractions import Fraction as F
from itertools import combinations, combinations_with_replacement, permutations
from math import lcm, prod
from pathlib import Path
from time import perf_counter

import sympy as s

import BN_h6_invariant_grade9 as BN
import BM_mixed_grade9 as BM
import BQ_h6_invariant_forms as BQ


GR = BN.GRAM
N = len(GR)
HERE = Path(__file__).resolve().parent
READ = ((0, 0), (0, 1), (0, 7), (0, 10))
WEDGE = 15358
POWERS = {mask: prod(GR[j] for j in BN.BC.bit_indices(mask))
          for mask in range(BN.FULL + 1)}


def add(table, key, value):
    BN.BC.plus(table, key, value)


def read_gzip(name):
    payload = gzip.decompress((HERE / name).read_bytes())
    return json.loads(payload), hashlib.sha256(payload).hexdigest()


def dict_tensor(series):
    # Integer entries with a common divisor are the on-disk representation.
    divisor = int(series['denominator'])
    data = {(a, m): F(c, divisor) for a, m, c in series['entries']}
    assert divisor > 0 and len(data) == len(series['entries'])
    assert all(0 <= a < N and m.bit_count() == 6 and c for (a, m), c in data.items())
    return data


def verify_detection(scope):
    def character(p):
        f = sum(i == p[i] for i in range(5))
        f2 = sum(i == p[p[i]] for i in range(5))
        return (f * f + f2) // 2 - 1

    gg = tuple(permutations(range(5)))
    plus = minus = 0
    for p in gg:
        chi, square = character(p), character(tuple(p[p[i]] for i in range(5)))
        plus += chi**2 + square
        minus += chi**2 - square
    assert (plus, minus) == (1200, 240)

    gram = s.diag(*GR)
    symmetric = [s.Matrix([[s.Rational(v) for v in row] for row in m])
                 for m in scope['symmetric_basis']]
    alternating = s.Matrix([[s.Rational(v) for v in row] for row in scope['alternating_basis']])
    assert len(symmetric) == 4
    assert all(m == m.T and s.trace(gram.inv() * m) == 0 for m in symmetric)
    assert alternating == -alternating.T
    assert s.Matrix([[m[a, b] for m in symmetric] for a, b in READ]).det() == 1
    assert WEDGE == BN.FULL ^ 1 ^ (1 << 10)
    assert -alternating[0, 10] / (GR[0] * GR[10]) == -F(5, 16)
    for i in range(4):
        action = s.zeros(N)
        for a, row in enumerate(BM.swap_rows(i, i + 1)):
            for b, value in row:
                action[a, b] = value
        assert action.det() == 1
        assert all(action.T * m * action == m for m in (*symmetric, alternating))


def contractions(tensor):
    by_index = [{} for _ in range(N)]
    inverse = [[{} for _ in range(N)] for _ in range(N)]
    for (a, m), value in tensor.items():
        by_index[a][m] = value
        for k in BN.BC.bit_indices(m):
            removed = m ^ (1 << k)
            sign = (-1) ** (removed & ((1 << k) - 1)).bit_count()
            inverse[a][k][removed] = sign * value
    return by_index, inverse


def metric_dot(a, b):
    if len(a) > len(b):
        a, b = b, a
    return sum((v * b[m] / POWERS[m] for m, v in a.items() if m in b), F(0))


def exterior_dot(a, b):
    total = F(0)
    for c, (left, right) in enumerate(zip(a, b)):
        for m, value in left.items():
            if m & ~WEDGE:
                continue
            n = WEDGE ^ m
            if n in right:
                total += BN.BC.sign(m, n) * value * right[n] / GR[c]
    return total


def lower_rows(a, b, gram):
    ax, ai = a
    bx, bi = b
    output = []
    for i, j in READ:
        ext = metric_dot(ax[i], bx[j]) + metric_dot(bx[i], ax[j])
        inside = sum((metric_dot(ai[c][i], bi[c][j]) + metric_dot(bi[c][i], ai[c][j]))
                     / (2 * GR[c]) for c in range(N))
        output.append(ext + inside - (F(4, 7) * GR[i] * gram if i == j else 0))
    return (*output, exterior_dot(ax, bx))


def sphere_certificate(poly, fixed):
    assert len(poly) == 8 and all(p == p.T and p.shape == (9, 9) for p in poly)
    radius = F(fixed['box_radius'])
    inverse = s.Matrix([[s.Rational(v) for v in row] for row in fixed['inverse']])
    assert radius == F(1, 10**8) and inverse.det()
    equations = [sum((inverse[i, j] * poly[j] for j in range(8)), s.zeros(9))
                 for i in range(8)]
    contractions = []
    images = []
    for i, q in enumerate(equations):
        defect = sum(abs(F(j == i) - 2 * q[0, j + 1]) for j in range(8))
        lipschitz = 2 * sum(abs(q[j, k]) for j in range(1, 9) for k in range(1, 9))
        derivative = defect + radius * lipschitz
        contractions.append(derivative)
        images.append(abs(q[0, 0]) + radius * derivative)
    assert max(contractions) < F(7, 10**7) and max(images) < F(105, 10**11) < radius
    print('Rational contraction < 7e-7; image radius < 1.05e-9 inside box 1e-8.')
    return radius


def physical_field(tensor, algebra):
    output = [{} for _ in GR]
    for (a, mask), value in tensor.items():
        target, factor = algebra.multiply(mask, BN.FULL)
        output[a][target] = value * factor
    return output


def check_separating(tensors, data, radius):
    a, mask = data['target']
    assert mask.bit_count() == 9
    algebra = BN.X.BladeAlgebra(GR)
    converted = [physical_field(t, algebra) for t in tensors]
    pairs = tuple(combinations_with_replacement(range(9), 2))
    assert tuple(map(tuple, data['pairs'])) == pairs
    diagonal = [BN.PhysicalQuadratic(field).coefficient(a, mask) for field in converted]
    assert diagonal[0] == BN.DETERMINANT * algebra.quadratic_coefficient(converted[0], a, mask)
    assert diagonal[0] > 100
    coefficients = []
    for i, j in pairs:
        if i == j:
            coefficients.append(diagonal[i])
        else:
            combined = [BN.BC.add(x, y) for x, y in zip(converted[i], converted[j])]
            mixed = BN.PhysicalQuadratic(combined).coefficient(a, mask)
            coefficients.append(mixed - diagonal[i] - diagonal[j])
    assert list(map(F, data['coefficients'])) == coefficients
    variation = sum(abs(c) * radius**(int(i != 0) + int(j != 0))
                    for c, (i, j) in zip(coefficients[1:], pairs[1:]))
    assert variation == F(data['exact_variation_bound'])
    assert variation < F(584, 10**8) and diagonal[0] > variation
    print('All 45 coefficients verified for grade-nine receiver', (a, mask))
    print('Separating coefficient at the center > 100; variation < 5.84e-6 on the entire box.')


def check_moments(tensors, poly, full):
    """Use exact Casimir components and independently calibrated scalar ratios."""
    pieces, _ = read_gzip('partial_isotypes.json.gz')
    components = []
    components.append({tuple(block['weight']): dict_tensor(block)
                       for block in pieces['anchor_isotypic_components']})
    assert BN.BF.combine(*((1, tensor) for tensor in components[0].values())) == tensors[0]
    assert len(pieces['direction_weights']) == 8
    components.extend({tuple(weight): tensor}
                      for weight, tensor in zip(pieces['direction_weights'], tensors[1:]))
    source = Path(BQ.__file__).with_name('BQ_h6_invariant_forms.json')
    published = json.loads(source.read_text(encoding='utf-8'))['blocks']
    scalar = {tuple(r['weight']): tuple(F(r['forms'][str(q)][0][0]) / F(r['gram'][0][0])
                                       for q in (5, 9)) for r in published if r['copies'] == 1}
    for component in components:
        for weight, tensor in component.items():
            l, m = weight
            assert weight in scalar
            if full:
                assert BM.full_casimir(tensor) == BN.BF.combine((l * (l + 3) + m * (m + 1), tensor))
    for i, j in combinations_with_replacement(range(9), 2):
        shared = components[i].keys() & components[j].keys()
        pairings = {w: F(BN.BF.inner(components[i][w], components[j][w])) for w in shared}
        assert sum(pairings.values()) == BN.BF.inner(tensors[i], tensors[j])
        for k in range(2):
            assert sum(v * scalar[w][k] for w, v in pairings.items()) == F(poly[5 + k][i, j])
    print('Both invariant rows agree with exact single-copy Casimir calibration.'
          if full else 'Both invariant row matrices agree with the packaged Casimir components.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--full', action='store_true',
                        help='also check all affine generators, all Casimir actions and all 45 higher coefficients')
    options = parser.parse_args()
    start = perf_counter()
    data = json.loads((HERE / 'partial_polynomials.json').read_text(encoding='utf-8'))
    packaged, digest = read_gzip('partial_fields.json.gz')
    assert digest == data['raw_sha256']
    assert tuple(map(F, packaged['metric_gram'])) == GR
    tensors = [dict_tensor(f) for f in packaged['fields']]
    assert len(tensors) == 9
    verify_detection(data['receiver_scope'])
    if options.full:
        curvature = BN.BG.curvature_model()
        for idx, t in enumerate(tensors):
            assert BN.BC.hook_traces(t) == ({}, {})
            assert not BN.BG.full_sparse_curvature(curvature, t)
            assert all(BM.swap(t, i, i + 1) == BN.BF.combine((-1, t)) for i in range(4))
            print('Rational affine field', idx, 'passes hooks, curvature and permutation sign.', flush=True)
    poly = [s.Matrix([[s.Rational(v) for v in row] for row in matrix])
            for matrix in data['exact_quadratic_row_matrices']]
    contracted = list(map(contractions, tensors))
    selected = tuple(combinations_with_replacement(range(9), 2)) if options.full else ((0, 0), (1, 3), (6, 8))
    for i, j in selected:
        gram = F(BN.BF.inner(tensors[i], tensors[j]))
        rows = lower_rows(contracted[i], contracted[j], gram)
        assert tuple(F(p[i, j]) for p in poly[:5]) == rows, (i, j, rows)
        assert F(poly[-1][i, j]) == gram - int(i == j == 0)
    check_moments(tensors, poly, options.full)
    radius = sphere_certificate(poly, data['fixed_point'])
    if options.full:
        check_separating(tensors, data['separating_receiver'], radius)
    else:
        receiver = data['separating_receiver']
        coefficients = list(map(F, receiver['coefficients']))
        pairs = tuple(combinations_with_replacement(range(9), 2))
        assert tuple(map(tuple, receiver['pairs'])) == pairs
        bound = sum(abs(c) * radius**(int(i != 0) + int(j != 0))
                    for c, (i, j) in zip(coefficients[1:], pairs[1:]))
        assert coefficients[0] > 100 and bound == F(receiver['exact_variation_bound']) < F(584, 10**8)
        print('Packaged grade-nine polynomial stays positive throughout the box.')
    print(f'Runtime {perf_counter() - start:.3f} seconds.')


if __name__ == '__main__':
    main()
