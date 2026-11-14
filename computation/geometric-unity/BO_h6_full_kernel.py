"""Exact upper and lower dimension bounds for the full H6 curvature kernel."""

import json
from collections import defaultdict
from fractions import Fraction as F
from itertools import combinations
from math import comb, isqrt, prod
from pathlib import Path
from time import perf_counter

import numpy as np
import sympy as s


N = 14
PRIME = 65537
PAIRS = tuple(combinations(range(N), 2))


def commutator(a, b):
    return a * b - b * a


def matrix_units(i, j):
    result = s.zeros(5)
    result[i, j] = 1
    return result


def build_model():
    metric = s.zeros(5)
    metric[0, 1] = metric[1, 0] = metric[2, 3] = metric[3, 2] = 2
    metric[4, 4] = 1
    axis_weights = [(1, 0), (-1, 0), (0, 1), (0, -1), (0, 0)]
    basis, weights = [], []
    for i in range(5):
        for j in range(i, 5):
            if (i, j) in ((0, 1), (2, 3), (4, 4)):
                continue
            matrix = matrix_units(i, j)
            basis.append((matrix if i == j else matrix + matrix.T) * metric)
            weights.append(tuple(a + b for a, b in zip(axis_weights[i], axis_weights[j])))
    for i, j in ((0, 1), (2, 3)):
        basis.append((matrix_units(i, j) + matrix_units(j, i)
                      - 4 * matrix_units(4, 4)) * metric)
        weights.append((0, 0))

    assert len(basis) == N
    assert all(matrix.trace() == 0 and matrix.T * metric == metric * matrix
               for matrix in basis)
    flat = s.Matrix.hstack(*(s.Matrix(list(matrix)) for matrix in basis))
    gram = flat.T * s.Matrix.hstack(*(s.Matrix(list(matrix.T)) for matrix in basis))
    inverse = gram.inv()

    def adjoint_action(matrix):
        def column(a):
            value = commutator(matrix, a)
            return s.Matrix(list(value.T))
        return inverse * flat.T * s.Matrix.hstack(*(column(a) for a in basis))

    covectors = {}
    for i, j in PAIRS:
        curvature = -adjoint_action(commutator(basis[i], basis[j])) / 4
        assert gram * curvature + curvature.T * gram == s.zeros(N)
        for a in range(N):
            covectors[a, i, j] = tuple((b, -F(curvature[a, b]))
                                        for b in range(N) if curvature[a, b])

    h1, h2 = s.diag(1, -1, 0, 0, 0), s.diag(0, 0, 1, -1, 0)
    for axis, h in enumerate((h1, h2)):
        assert all(commutator(h, a) == weight[axis] * a
                   for weight, a in zip(weights, basis))
    f1 = matrix_units(2, 0) - matrix_units(1, 3)
    f2 = matrix_units(3, 4) - 2 * matrix_units(4, 2)
    f12 = commutator(f1, f2)
    f122 = commutator(f12, f2)
    negative = (f1, f2, f12, f122)
    swap = s.zeros(5)
    for i, j in enumerate((1, 0, 3, 2, 4)):
        swap[i, j] = 1
    generators = (h1, h2, *negative, *(swap * f * swap for f in negative))
    assert all(q.T * metric + metric * q == s.zeros(5) for q in generators)
    assert s.Matrix.hstack(*(s.Matrix(list(q)) for q in generators)).rank() == 10
    for matrix, weight in zip(negative, ((-1, 1), (0, -1), (-1, 0), (-1, -1))):
        assert all(commutator(h, matrix) == w * matrix for h, w in zip((h1, h2), weight))
    root_actions = []
    for root in (f1, f2):
        vector_action = adjoint_action(root)
        root_actions.append(tuple(tuple((b, -F(vector_action[a, b]))
                                         for b in range(N) if vector_action[a, b])
                                  for a in range(N)))
    return weights, tuple(tuple(F(c) for c in row) for row in inverse.tolist()), covectors, root_actions


def axes(mask):
    return [i for i in range(N) if mask >> i & 1]


def wedge_sign(index, mask):
    return (-1) ** (mask & ((1 << index) - 1)).bit_count()


def add(table, key, value):
    combined = table.get(key, 0) + value
    if combined:
        table[key] = combined
    elif key in table:
        del table[key]


def constraint_entries(a, mask, inverse, curvature):
    """Yield one column of the alternating, contraction and curvature stack."""
    if not mask >> a & 1:
        yield (0, 0, mask | (1 << a)), wedge_sign(a, mask)
    for k in axes(mask):
        if inverse[a][k]:
            yield (1, 0, mask ^ (1 << k)), wedge_sign(k, mask) * inverse[a][k]
    missing = [i for i in range(N) if not mask >> i & 1]
    for i, j in combinations(missing, 2):
        target = mask | (1 << i) | (1 << j)
        sign = wedge_sign(i, mask) * wedge_sign(j, mask)
        for b, c in curvature[a, i, j]:
            yield (2, b, target), sign * c


def root_action(tensor, covectors):
    out = {}
    for (a, mask), value in tensor.items():
        for b, c in covectors[a]:
            add(out, (b, mask), value * c)
        for i in axes(mask):
            reduced = mask ^ (1 << i)
            remove = wedge_sign(i, mask)
            for j, c in covectors[i]:
                if not reduced >> j & 1:
                    add(out, (a, reduced | (1 << j)), value * c * remove * wedge_sign(j, reduced))
    return out


def inventory(weights):
    forms = {d: [sum(1 << i for i in js) for js in combinations(range(N), d)]
             for d in (5, 6, 7, 8)}
    weight = {mask: tuple(sum(weights[i][k] for i in axes(mask)) for k in range(2))
              for candidates in forms.values() for mask in candidates}

    def slot(a, mask):
        return tuple(weights[a][k] + weight[mask][k] for k in range(2))

    columns = defaultdict(list)
    rows = defaultdict(list)
    for a in range(N):
        for mask in forms[6]:
            columns[slot(a, mask)].append((a, mask))
        for mask in forms[8]:
            rows[slot(a, mask)].append((2, a, mask))
    for kind, degree in ((0, 7), (1, 5)):
        for mask in forms[degree]:
            rows[weight[mask]].append((kind, 0, mask))
    assert len(columns) == 185 and sum(map(len, columns.values())) == 42042
    assert sum(map(len, rows.values())) == 47476
    return columns, rows, slot


def modular_rank(matrix):
    rank = 0
    rows, cols = matrix.shape
    for col in range(cols):
        choices = np.flatnonzero(matrix[rank:, col])
        if not len(choices):
            continue
        pivot_row = rank + int(choices[0])
        if pivot_row != rank:
            matrix[[rank, pivot_row]] = matrix[[pivot_row, rank]]
        active = np.flatnonzero(matrix[rank + 1:, col]) + rank + 1
        inverse = pow(int(matrix[rank, col]), -1, PRIME)
        normalized = matrix[rank, col:] * inverse % PRIME
        matrix[active, col:] = (matrix[active, col:]
                                - matrix[active, col, None] * normalized) % PRIME
        rank += 1
        if rank == min(rows, cols):
            break
    return rank


def rank_lower_bound(columns, rows, inverse, curvature):
    assert all(PRIME % d for d in range(2, isqrt(PRIME) + 1))

    def reduce_number(number):
        return number.numerator * pow(number.denominator, -1, PRIME) % PRIME

    mod_inverse = tuple(tuple(reduce_number(v) for v in row) for row in inverse)
    mod_curvature = {key: tuple((b, reduce_number(c)) for b, c in row)
                     for key, row in curvature.items()}
    total = 0
    for weight, block in columns.items():
        labels = {label: i for i, label in enumerate(rows[weight])}
        matrix = np.zeros((len(labels), len(block)), dtype=np.int64)
        for j, (a, mask) in enumerate(block):
            for label, value in constraint_entries(a, mask, mod_inverse, mod_curvature):
                matrix[labels[label], j] += value
        matrix %= PRIME
        total += modular_rank(matrix)
    return total


def kernel_lower_bound(entries, slot, inverse, curvature, root_actions):
    total, copies = 0, 0
    seen = set()
    for item in entries:
        l, m = item['weight']
        assert l >= m >= 0 and (l, m) not in seen and item['vectors']
        seen.add((l, m))
        leading = {}
        for vector in item['vectors']:
            tensor = {}
            for a, mask, value in vector:
                assert 0 <= a < N and 0 < mask < 1 << N and mask.bit_count() == 6
                assert value and (a, mask) not in tensor and slot(a, mask) == (l, m)
                tensor[a, mask] = F(value)
            constraints = {}
            for (a, mask), value in tensor.items():
                for label, coefficient in constraint_entries(a, mask, inverse, curvature):
                    add(constraints, label, coefficient * value)
            assert not constraints and all(not root_action(tensor, root) for root in root_actions)
            remainder = tensor.copy()
            for key, previous in leading.items():
                factor = remainder.get(key, 0)
                if factor:
                    for label, value in previous.items():
                        add(remainder, label, -factor * value)
            assert remainder
            key = min(remainder)
            factor = remainder[key]
            leading[key] = {label: value / factor for label, value in remainder.items()}
        dimension = (l - m + 1) * (l + m + 2) * (2 * l + 3) * (2 * m + 1) // 6
        copies += len(leading)
        total += len(leading) * dimension
        print('Weight', (l, m), 'copies', len(leading), 'dimension each', dimension)
    assert len(seen) == 15 and copies == 18
    return total


def check_negative_grade_nine(entries):
    import BN_h6_invariant_grade9 as BN
    import BG_h6_curvature as BG

    tensor = {(a, mask): F(value) for a, mask, value in entries}
    assert len(tensor) == len(entries) == 96
    assert BN.BC.hook_traces(tensor) == ({}, {})
    assert not BG.full_sparse_curvature(BG.curvature_model(), tensor)
    norm = BN.BF.inner(tensor, tensor)
    assert norm == F(1, 2048)
    receiver = BN.invariant_receiver()
    algebra = BN.X.BladeAlgebra(BN.GRAM)
    field = [{} for _ in range(N)]
    for (a, mask), value in tensor.items():
        target, factor = algebra.multiply(mask, BN.FULL)
        assert mask.bit_count() == 6 and target.bit_count() == 8
        field[a][target] = value * factor
    quadratic = BN.PhysicalQuadratic(field)
    value = sum(c * quadratic.coefficient(a, mask)
                / (BN.GRAM[a] * prod(BN.GRAM[j] for j in axes(mask)))
                for (a, mask), c in receiver.items())
    assert value / norm == -F(1, 384)
    print('Real negative grade-nine witness: support', len(tensor),
          'norm squared', norm, 'projection / norm squared', value / norm)


def main():
    start = perf_counter()
    path = Path(__file__).with_name('BO_h6_full_kernel_witnesses.json')
    certificate = json.loads(path.read_text(encoding='utf-8'))
    weights, inverse, curvature, root_actions = build_model()
    columns, rows, slot = inventory(weights)
    lower = kernel_lower_bound(certificate['extremal_vectors'], slot, inverse,
                               curvature, root_actions)
    print('Exact extremal-vector lower bound:', lower, flush=True)
    rank = rank_lower_bound(columns, rows, inverse, curvature)
    upper = N * comb(N, 6) - rank
    assert lower == upper == 6834
    print('Rank modulo', PRIME, 'on 185 weight blocks:', rank)
    print('H6 dimension 36608; curvature rank on H6 29774; kernel dimension 6834.')
    check_negative_grade_nine(certificate['negative_grade_nine_witness'])
    print(f'Runtime {perf_counter() - start:.3f} seconds.')


if __name__ == '__main__':
    main()
