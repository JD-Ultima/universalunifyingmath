"""Reconstruct the grade-nine receiver and calibrate it on all three BG modules."""

import argparse
import json
from fractions import Fraction as F
from itertools import combinations
from math import gcd, lcm, prod
from pathlib import Path
from time import perf_counter

import BC_affine_hooks as BC
import BC_quadratic_xcheck as X
import BF_h6_casimir as BF
import BG_h6_curvature as BG


GRAM = tuple(BC.METRIC)
FULL = (1 << BC.N) - 1
DETERMINANT = prod(GRAM)


def invariant_receiver():
    """Coefficients of tr(Q_a (sum Q_b e^b)^9), divided by their gcd."""
    generators = BG.integer_phi()
    sparse = [
        [(i, j, value) for i, row in enumerate(matrix)
         for j, value in enumerate(row) if value]
        for matrix in generators
    ]
    powers = {0: tuple(int(i == j) for i in range(5) for j in range(5))}
    for degree in range(1, 10):
        next_power = {}
        for indices in combinations(range(BC.N), degree):
            mask = sum(1 << i for i in indices)
            matrix = [0] * 25
            for position, generator in enumerate(indices):
                previous = powers.get(mask ^ (1 << generator))
                if previous is None:
                    continue
                sign = (-1) ** (degree - 1 - position)
                for i, j, value in sparse[generator]:
                    for row in range(5):
                        matrix[5 * row + j] += sign * previous[5 * row + i] * value
            if any(matrix):
                next_power[mask] = matrix
        powers = next_power

    raw = {}
    for mask, matrix in powers.items():
        for a, entries in enumerate(sparse):
            value = sum(c * matrix[5 * j + i] for i, j, c in entries)
            if value:
                raw[a, mask] = value
    divisor = gcd(*raw.values())
    assert len(raw) == 1143 and divisor == 10**10
    return {key: F(value, divisor) for key, value in raw.items()}


def seeds():
    """Use independent public inputs for the cubic/quintic projections and H2."""
    source = Path(__file__).with_name('BC_h6_closed_generator.json')
    closed = {(a, mask): F(value)
              for a, mask, value in json.loads(source.read_text(encoding='utf-8'))['entries']}
    casimir = BF.casimir(closed)
    source = Path(__file__).with_name('BG_h6_curvature_kernel.json')
    quadratic = json.loads(source.read_text(encoding='utf-8'))['kernel_tensors'][2]
    return {
        2: {(a, mask): F(value) for a, mask, value in quadratic['entries']},
        3: BF.combine((F(40, 22), closed), (-F(1, 22), casimir)),
        5: BF.combine((F(1, 22), casimir), (-F(18, 22), closed)),
    }


def lower_prefix(mask):
    """Bit j records the parity of generators strictly below j."""
    for shift in (1, 2, 4, 8):
        mask ^= mask << shift
    return (mask << 1) & FULL


class PhysicalQuadratic:
    """BC's exact single-coefficient formula with integer inner loops."""

    def __init__(self, field):
        self.denominator = lcm(*(v.denominator for part in field for v in part.values()))
        self.field = [
            {mask: int(value * self.denominator) for mask, value in part.items()}
            for part in field
        ]
        self.prefix = [lower_prefix(m) for m in range(FULL + 1)]
        # det(g) divided by the inverse-Clifford repetition denominator.
        # Multiplying by 5 clears the only nonintegral Gram factor.
        self.weights = [
            int(5 * prod(GRAM[j] for j in range(BC.N) if not mask >> j & 1))
            for mask in range(FULL + 1)
        ]

    def coefficient(self, a, target):
        total = 0
        for b, c in combinations(range(BC.N), 2):
            if a in (b, c):
                continue
            triple = (1 << a) | (1 << b) | (1 << c)
            ordering = -1 if b < a < c else 1
            for x, vx in self.field[b].items():
                y = triple ^ x ^ target
                vy = self.field[c].get(y)
                if vy is None:
                    continue
                tx = (triple & x).bit_count() % 2
                ty = (triple & y).bit_count() % 2
                xy = (x & y).bit_count() % 2
                parity = (1 - (-1)**xy + (-1)**tx + (-1)**(tx + ty)
                          - (-1)**(xy + ty) - (-1)**(tx + ty + xy))
                if parity:
                    permutation = (triple & self.prefix[x]) ^ ((triple ^ x) & self.prefix[y])
                    sign = -1 if permutation.bit_count() % 2 else 1
                    repeated = (triple & x) | (triple & y) | (x & y)
                    total += ordering * parity * sign * self.weights[repeated] * vx * vy
        return -F(2, 15) * GRAM[a] * F(total, self.denominator**2)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--full-nested', action='store_true',
                        help='compare every receiver coefficient with the separate nested Clifford engine')
    args = parser.parse_args()
    started = perf_counter()
    receiver = invariant_receiver()
    assert all(not BC.action(j, receiver) for j in range(10))
    receiver_norm = BF.inner(receiver, receiver)
    assert receiver_norm == F(2835, 2097152)
    weights = {
        (a, mask): value / (GRAM[a] * prod(GRAM[j] for j in BC.bit_indices(mask)))
        for (a, mask), value in receiver.items()
    }
    print('Receiver: 1143 entries, gcd 10000000000, norm squared', receiver_norm)
    print('All ten infinitesimal SO(5) actions vanish.')

    algebra = X.BladeAlgebra(GRAM)
    assert algebra.multiply(FULL, FULL) == (0, -1 / DETERMINANT)
    expected = {
        2: (F(18225, 16), F(26325, 2048), F(13, 1152)),
        3: (F(313875, 88), F(25875, 1408), F(23, 4464)),
        5: (F(455625, 45056), F(354375, 5767168), F(7, 1152)),
    }
    guard = set(list(receiver)[:2] + list(receiver)[-2:])
    for degree, seed in seeds().items():
        assert BC.hook_traces(seed) == ({}, {})
        assert all(not BC.action(j, seed) for j in range(6))
        assert BF.casimir(seed) == BF.combine((degree * (degree + 3), seed))
        field = [{} for _ in range(BC.N)]
        for (a, mask), value in seed.items():
            dual, factor = algebra.multiply(mask, FULL)
            assert mask.bit_count() == 6 and dual.bit_count() == 8
            field[a][dual] = value * factor
        engine = PhysicalQuadratic(field)
        norm = BF.inner(seed, seed)
        projection = F(0)
        for (a, mask), weight in weights.items():
            coefficient = engine.coefficient(a, mask)
            if args.full_nested or (a, mask) in guard:
                separate = DETERMINANT * algebra.quadratic_coefficient(field, a, mask)
                assert coefficient == separate, (degree, a, mask, coefficient, separate)
            projection += weight * coefficient
        ratio = projection / norm
        assert (norm, projection, ratio) == expected[degree]
        print('Harmonic degree', degree, 'norm squared', norm,
              'projection', projection, 'ratio', ratio, flush=True)
        algebra.multiply.cache_clear()

    print('Ordered nested Clifford comparisons:', 3429 if args.full_nested else 12)
    print('35712 L9 = 403 ||B2||^2 + 184 ||B3||^2 + 217 ||B5||^2.')
    print(f'Runtime {perf_counter() - started:.3f} seconds.')


if __name__ == '__main__':
    main()
