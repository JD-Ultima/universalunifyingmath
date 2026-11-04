"""Independent Clifford moments on the explicit affine H6 closure generator.

Sorted-generator multiplication and nested cubic-gradient lookup provide a
separate implementation from BC_affine_hooks.py. The input certificate is
independently reconstructed and compared there. Only three diagonal targets
are computed here; running writes no files and imports no earlier engine.
"""
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations
from math import prod
from pathlib import Path
from time import perf_counter
import hashlib
import json

N = 14
FULL = (1 << N) - 1


class BladeAlgebra:
    def __init__(self, gram):
        self.gram = tuple(map(F, gram))

    @lru_cache(None)
    def multiply(self, left, right):
        """Append each right generator to the sorted left blade."""
        current, factor = left, F(1)
        for j in range(N):
            if right & (1 << j):
                if (current >> (j + 1)).bit_count() % 2:
                    factor = -factor
                if current & (1 << j):
                    factor /= self.gram[j]
                current ^= 1 << j
        assert current == left ^ right
        return current, factor

    def commutator_coefficient(self, left, right, target):
        result = F(0)
        for x, vx in left.items():
            y = x ^ target
            if y in right:
                _, xy = self.multiply(x, y)
                _, yx = self.multiply(y, x)
                result += vx * right[y] * (xy - yx)
        return result

    def monomial_anticommutator(self, mask, coefficient, field):
        result = {}
        for x, value in field.items():
            target, ax = self.multiply(mask, x)
            _, xa = self.multiply(x, mask)
            term = coefficient * value * (ax + xa)
            if term:
                result[target] = term
        return result

    def quadratic_coefficient(self, field, a, target):
        total = F(0)
        for b, c in combinations(range(N), 2):
            if a in (b, c):
                continue
            m, first = self.multiply(1 << a, 1 << b)
            m, second = self.multiply(m, 1 << c)
            triple = first * second
            needed = m ^ target
            bc = self.commutator_coefficient(field[b], field[c], needed)
            _, am = self.multiply(m, needed)
            _, ma = self.multiply(needed, m)
            total += triple * (am + ma) * bc
            ab = self.monomial_anticommutator(m, triple, field[b])
            ac = self.monomial_anticommutator(m, triple, field[c])
            total += self.commutator_coefficient(ab, field[c], target)
            total -= self.commutator_coefficient(ac, field[b], target)
        return -self.gram[a] * total / 3


def main():
    started = perf_counter()
    source = Path(__file__).with_name('BC_h6_closed_generator.json')
    data = json.loads(source.read_text(encoding='utf-8'))
    assert data['degree'] == 6 and len(data['entries']) == 1358
    gram = tuple(map(F, data['metric_gram']))
    assert gram == (F(2),) * 4 + (F(8),) * 7 + (F(24), F(48), F(16, 5))
    algebra = BladeAlgebra(gram)
    field = [{} for _ in range(N)]
    for a, mask, value in data['entries']:
        assert 0 <= a < N and mask.bit_count() == 6
        assert mask not in field[a]
        field[a][mask] = F(value)

    alt, ctr = {}, {}
    norm = F(0)
    converted = [{} for _ in range(N)]
    for a, item in enumerate(field):
        for mask, value in item.items():
            target, factor = algebra.multiply(1 << a, mask)
            trace = ctr if mask & (1 << a) else alt
            trace[target] = trace.get(target, F(0)) + value * factor
            norm += value * value / (gram[a] * prod(gram[j] for j in range(N)
                                                   if mask & (1 << j)))
            target, factor = algebra.multiply(mask, FULL)
            assert target.bit_count() == 8
            converted[a][target] = value * factor
    assert not any(alt.values()) and not any(ctr.values())
    determinant = prod(gram)
    assert algebra.multiply(FULL, FULL) == (0, -1 / determinant)

    # Independent hand normalization: each admissible pair gives 4e_a.
    unit = BladeAlgebra([1] * N)
    phi = [{1 << a: F(1)} for a in range(N)]
    assert unit.quadratic_coefficient(phi, 0, 1) == 312

    representatives = (0, 4, 13)
    values = tuple(determinant * algebra.quadratic_coefficient(converted, a, 1 << a)
                   / gram[a] for a in representatives)
    expected = (F(-15287625, 4096), F(-10178375, 2048), F(-4793625, 2048))
    assert norm == F(14650875, 4096)
    assert values == expected
    assert len(set(values)) == 3
    # The weighted trace uses SO4 invariance of the supplied tensor; the nine
    # Sym0 directions are one irreducible representation, not nine new targets.
    trace = 4 * values[0] + 9 * values[1] + values[2]
    assert trace == -F(52, 3) * norm
    print('1358 real all-covariant entries: degree6, unique slots, Alt=0, Ctr=0.')
    print('Independent sorted-generator products; nested cubic gradient; only H0,O01,I targets.')
    print('Volume square and hand 312Phi normalization verified.')
    print('Primitive norm^2:', norm)
    print('Q1 / metric, H0 O01 I:', *values)
    print('Norm-normalized:', *(value / norm for value in values))
    print('SO4-weighted trace:', trace, '= -(52/3) norm^2.')
    print('All three moments differ: no nonzero real amplitude has isotropic Q1.')
    print('Input SHA-256:', hashlib.sha256(source.read_bytes()).hexdigest())
    print('No closure rank, full Clifford output or field solution was recomputed.')
    print(f'Runtime {perf_counter() - started:.3f} seconds.')


if __name__ == '__main__':
    main()
