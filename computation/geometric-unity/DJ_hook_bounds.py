"""Verify the exact Cartan incidence bounds and Clifford calibrations in DJ."""

from collections import Counter
from fractions import Fraction
from itertools import combinations, combinations_with_replacement
import json
from math import comb
from pathlib import Path

import CD_clifford_hook_identities as clifford
from CN_h4_s5_cubic import weight_multiplicities


HERE = Path(__file__).resolve().parent


def cartan_bounds():
    defining_weights = (1, 0, 0, 0, -1)
    tangent = [a + b for a, b in combinations_with_replacement(defining_weights, 2)]
    tangent.remove(0)
    assert Counter(tangent) == {2: 1, 1: 3, 0: 6, -1: 3, -2: 1}
    dimensions = Counter(map(abs, tangent))
    result = {}
    for degree, sign, scalar, certificates in (
        (4, -2, Fraction(1, 7), [
            (6, (2, 4, 3), 7, Fraction(6, 7)),
            (7, (1, 2, 1), 5, Fraction(4, 7)),
        ]),
        (6, 2, Fraction(4, 7), [
            (6, (0, 2, 3), 13, Fraction(72, 91)),
            (7, (0, 0, 1), 3, Fraction(8, 21)),
        ]),
    ):
        checked = 0
        for form in combinations(tangent, degree):
            weight = sum(form)
            counts = Counter(map(abs, form))
            for external in tangent:
                checked += 1
                incidence = [counts[j] + sign * (abs(external) == j) for j in range(3)]
                for threshold, dual, denominator, _ in certificates:
                    lhs = sum(coefficient * value for coefficient, value in zip(dual, incidence))
                    assert lhs >= denominator * (abs(weight + external) >= threshold)
        assert checked == 14 * comb(14, degree)
        for _, dual, denominator, bound in certificates:
            assert scalar * sum(dimensions[j] * value for j, value in enumerate(dual)) == denominator * bound
        result[degree] = {threshold: str(bound) for threshold, _, _, bound in certificates}

    for degree, cap in ((1, 4), (5, 7), (9, 7), (13, 4)):
        assert max(tangent) + sum(sorted(tangent, reverse=True)[:degree]) == cap
    for document, key, excluded, dimension, high_weight in (
        ("BX_data/h4_character.json", "components", None, 10388, 480),
        ("BQ_h6_invariant_forms.json", "blocks", [2, 0], 6820, 573),
    ):
        blocks = json.loads((HERE / document).read_text(encoding="utf-8"))[key]
        character = Counter()
        for block in blocks:
            if block["weight"] != excluded:
                for weight, multiplicity in weight_multiplicities(*block["weight"]).items():
                    character[weight] += block["copies"] * multiplicity
        assert sum(character.values()) == dimension
        assert sum(count for (weight, _), count in character.items() if weight >= 4) == high_weight
    return result


def commutator(left, right):
    return clifford.add(clifford.mul(left, right), clifford.scale(clifford.mul(right, left), -1))


def norm_squared(element):
    return sum(value * value for value in element.values())


def clifford_calibrations():
    pairs = list(combinations(range(14), 2))
    triples = list(combinations(range(14), 3))
    gamma = clifford.GAMMA
    spin = {(a, b): clifford.mul(gamma[a], gamma[b])
            for a, b in [*pairs, *((b, a) for a, b in pairs)]}
    records = []
    for degree, equation_cap, expected_baseline, expected_equation in (
        (4, Fraction(2), Fraction(5833, 264992), Fraction(191, 8281)),
        (6, Fraction(78, 43), Fraction(2501, 113568), Fraction(3257, 138411)),
    ):
        k = degree if degree == 4 else 8
        ratio = 4 * (k - 1) * (13 - k)
        baseline = Fraction(2, 91) + Fraction(ratio, 6 * 91**2 * 64)
        sharper = Fraction(2, 91) + Fraction(ratio, 6 * 91**2) / equation_cap
        assert (baseline, sharper) == (expected_baseline, expected_equation)
        for seed in (18004 + degree, 18014 + degree):
            hooks = clifford.hook_sample(degree, seed, terms=24)
            values = hooks if degree == 4 else [clifford.mul(b, {clifford.FULL: 1}) for b in hooks]
            size = sum(map(norm_squared, values))
            assert size
            assert not clifford.add(*(clifford.mul(gamma[a], values[a]) for a in range(14)))
            assert not clifford.add(*(clifford.mul(values[a], gamma[a]) for a in range(14)))
            brackets = {(a, b): commutator(values[a], values[b])
                        for a, b in [*pairs, *((b, a) for a, b in pairs)]}
            scalar_product = sum(
                sum(value * spin[a, b].get(mask, 0)
                    for mask, value in brackets[a, b].items())
                for a, b in pairs
            )
            assert scalar_product == -size
            external = [[sum(value * values[c].get(mask, 0)
                             for mask, value in values[a].items())
                         for c in range(14)] for a in range(14)]
            for a in range(14):
                for c in range(14):
                    ricci = sum(
                        sum(value * brackets[a, b].get(mask, 0)
                            for mask, value in spin[c, b].items())
                        for b in range(14) if b != a and b != c
                    )
                    assert ricci == -2 * external[a][c]
            jacobi_norm = sum(
                norm_squared(clifford.add(
                    commutator(values[a], spin[b, c]),
                    commutator(values[b], spin[c, a]),
                    commutator(values[c], spin[a, b]),
                )) for a, b, c in triples
            )
            assert jacobi_norm == ratio * size
            variance = sum((external[a][c] - Fraction(size, 14) * (a == c))**2
                           for a in range(14) for c in range(14))
            energy = sum(map(norm_squared, brackets.values()))
            assert energy >= Fraction(2, 91) * size**2 + Fraction(2, 3) * variance
        records.append((k, ratio, str(baseline), str(sharper)))
    return records


if __name__ == "__main__":
    print("Exact Cartan incidence bounds:", cartan_bounds())
    print("Exact Clifford calibrations:", clifford_calibrations())
