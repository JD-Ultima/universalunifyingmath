"""Reconstruct the two invariant quadratic forms on every H6 kernel type."""

import argparse
import json
from fractions import Fraction as F
from itertools import combinations
from math import gcd, lcm, prod
from pathlib import Path
from time import perf_counter

import sympy as s

import BN_h6_invariant_grade9 as BN
import BO_h6_full_kernel as BO
import BP_h6_cartan_grade9 as BP


ALLOWED = ((2, 0), (3, 0), (5, 0))
LARGE = ((2, 0), (3, 2), (5, 0), (5, 2), (6, 1),
         (6, 3), (7, 1), (7, 2), (7, 3))
BOUND = F(753, 1285)


def primitive(tensor):
    denominator = lcm(*(v.denominator for v in tensor.values()))
    divisor = gcd(*(int(v * denominator) for v in tensor.values()))
    return {key: F(v * denominator, divisor) for key, v in tensor.items()}


def receiver(degree):
    """Construct by dividing the trace of the ordered exterior power by its gcd."""
    forms = {0: s.eye(5)}
    generators = [s.Matrix(p) for p in BN.BG.integer_phi()]
    for k in range(1, degree + 1):
        next_forms = {}
        for indices in combinations(range(14), k):
            mask = sum(1 << j for j in indices)
            coefficient = s.zeros(5)
            for j, i in enumerate(indices):
                tail = forms.get(mask ^ (1 << i))
                if tail is not None:
                    coefficient += (-1)**j * (generators[i] * tail)
            if any(coefficient):
                next_forms[mask] = coefficient
        forms = next_forms
    raw = {}
    for mask, q in forms.items():
        for a, generator in enumerate(generators):
            value = int(s.trace(generator * q))
            if value:
                raw[a, mask] = F(value)
    return primitive(raw)


def weighted_receivers():
    data = {5: receiver(5), 9: BN.invariant_receiver()}
    expected = {5: F(14175, 512), 9: F(2835, 2097152)}
    for degree, coordinates in data.items():
        assert len(coordinates) == 1143
        assert BN.BF.inner(coordinates, coordinates) == expected[degree]
        assert all(not BN.BC.action(i, coordinates) for i in range(10))
    return {
        d: [(a, m, value / (BN.GRAM[a] * prod(BN.GRAM[i] for i in BN.BC.bit_indices(m))))
            for (a, m), value in r.items()]
        for d, r in data.items()
    }


def quadratic_forms(tensor, receivers, algebra):
    field = [{} for _ in range(14)]
    for (a, mask), value in tensor.items():
        target, factor = algebra.multiply(mask, BN.FULL)
        field[a][target] = value * factor
    engine = BN.PhysicalQuadratic(field)
    output = []
    for grade in (5, 9):
        total = F(0)
        for index, (a, mask, weight) in enumerate(receivers[grade]):
            coefficient = engine.coefficient(a, mask)
            if index in (0, 1, 1141, 1142):
                assert coefficient == BN.DETERMINANT * algebra.quadratic_coefficient(field, a, mask)
            total += weight * coefficient
        output.append(total)
    algebra.multiply.cache_clear()
    return output


def exact_block(vectors, change, receivers, algebra):
    seeds = [primitive(BP.transport(v, change)[0]) for v in vectors]
    products = [[BN.BF.inner(a, b) for b in seeds] for a in seeds]
    diagonal = [quadratic_forms(seed, receivers, algebra) for seed in seeds]
    matrices = [s.diag(*(row[i] for row in diagonal)) for i in range(2)]
    if len(seeds) == 2:
        a, b = seeds
        combined = quadratic_forms(BN.BF.combine((1, a), (1, b)), receivers, algebra)
        for i, matrix in enumerate(matrices):
            matrix[0, 1] = matrix[1, 0] = (combined[i] - diagonal[0][i] - diagonal[1][i]) / 2
    return s.Matrix(products), *matrices


def positive(matrix, strict=False):
    minors = [matrix[i, i] for i in range(matrix.rows)] + [matrix.det()]
    return min(minors) > 0 if strict else min(minors) >= 0


def support_checks(blocks, pencils):
    large_dimension = 0
    for weight, (dim, gram, a, b) in blocks.items():
        if weight in LARGE:
            assert positive(2 * a - b, strict=True)
            large_dimension += dim
        p = int(weight in ALLOWED)
        assert positive((BOUND - p) * gram - F(72, 6425) * a + F(51, 257) * b)

    assert large_dimension == 5290
    supports = []
    for item in pencils:
        weights = set(map(tuple, item['excluded_weights']))
        alpha, beta = F(item['alpha_L5']), F(item['beta_384L9'])
        assert item['dimension'] == sum(blocks[w][0] for w in weights)
        assert all(positive(alpha * blocks[w][2] + beta * blocks[w][3], strict=True)
                   for w in weights)
        supports.append(weights)
    assert all(any(w in support for support in supports) for w in blocks)
    count = sum(any({a, b} <= chosen for chosen in supports)
                for a, b in combinations(blocks, 2))
    assert count == 94

    # A point on the relaxation boundary, with scalar blocks for all three types.
    weights = ((3, 0), (5, 0), (6, 3))
    ratios = [(blocks[w][2][0, 0] / blocks[w][1][0, 0],
               blocks[w][3][0, 0] / blocks[w][1][0, 0]) for w in weights]
    constraints = s.Matrix([[1, 1, 1], [a for a, b in ratios], [b for a, b in ratios]])
    energies = constraints.inv() * s.Matrix([1, 0, 0])
    assert all(v > 0 for v in energies) and sum(energies[:2]) == BOUND
    print('Nine exact pencils exclude all 15 single types and 94 of 105 two-type supports.')
    print('Largest displayed positive pencil: 2 L5 - 384 L9, dimension 5290.')
    print('Maximum old-module energy in the invariant relaxation:', BOUND)
    print('Attaining energies for weights (3,0), (5,0), (6,3):', *energies)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--full', action='store_true',
                        help='recalibrate all 15 types instead of three representative blocks')
    options = parser.parse_args()
    start = perf_counter()
    here = Path(__file__).resolve().parent
    record = json.loads((here / 'BQ_h6_invariant_forms.json').read_text(encoding='utf-8'))
    data_path = Path(BO.__file__).with_name('BO_h6_full_kernel_witnesses.json')
    data = json.loads(data_path.read_text(encoding='utf-8'))
    by_weight = {tuple(item['weight']): item['vectors'] for item in data['extremal_vectors']}
    change = BP.split_to_physical_covectors()
    algebra = BN.X.BladeAlgebra(BN.GRAM)
    receivers = weighted_receivers()
    blocks = {}
    for item in record['blocks']:
        weight = tuple(item['weight'])
        gram = s.Matrix(item['gram'])
        a = s.Matrix(item['forms']['5'])
        b = s.Matrix(item['forms']['9'])
        assert positive(gram, strict=True)
        if options.full or weight in ((7, 3), (6, 1), (3, 2)):
            assert exact_block(by_weight[weight], change, receivers, algebra) == (gram, a, b)
            print('Weight', weight, 'copies', gram.rows, 'Gram and both forms agree', flush=True)
        assert gram.rows == item['copies']
        l, m = weight
        dimension = (l - m + 1) * (l + m + 2) * (2 * l + 3) * (2 * m + 1) // 6
        assert dimension == item['dimension'] and weight not in blocks
        blocks[weight] = (dimension * gram.rows, gram, a, 384 * b)
    assert len(blocks) == 15 and sum(block[0] for block in blocks.values()) == 6834
    support_checks(blocks, record['positive_pencils'])
    print(f'Runtime {perf_counter() - start:.3f} seconds.')


if __name__ == '__main__':
    main()
