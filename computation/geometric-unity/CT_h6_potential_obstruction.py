"""Verify the absence of type (5,5) in curvature-compatible potential jets."""

from fractions import Fraction
from itertools import combinations
import json
from pathlib import Path

import BO_h6_full_kernel as BO


def curvature_and_raising_matrix():
    weights, inverse, curvature, roots = BO.build_model()
    five_blades = [sum(1 << i for i in indices) for indices in combinations(range(14), 5)]
    form_weights = {
        mask: tuple(sum(weights[i][axis] for i in BO.axes(mask)) for axis in range(2))
        for mask in five_blades
    }
    columns = [
        (i, a, mask)
        for i in range(14)
        for a in range(14)
        for mask in five_blades
        if tuple(weights[i][j] + weights[a][j] + form_weights[mask][j]
                 for j in range(2)) == (5, 5)
    ]
    rows = {}
    for column, (i, a, mask) in enumerate(columns):
        for (kind, b, blade), value in BO.constraint_entries(a, mask, inverse, curvature):
            if kind == 2:
                BO.add(rows.setdefault((0, kind, i, b, blade), {}), column, value)
        for root, action in enumerate(roots):
            for b, value in action[i]:
                BO.add(rows.setdefault((1, root, b, a, mask), {}), column, value)
            for (b, blade), value in BO.root_action({(a, mask): Fraction(1)}, action).items():
                BO.add(rows.setdefault((1, root, i, b, blade), {}), column, value)
    return columns, rows


def determinant_modular(matrix, prime):
    determinant = 1
    for j in range(len(matrix)):
        pivot = next(i for i in range(j, len(matrix)) if matrix[i][j])
        if pivot != j:
            matrix[j], matrix[pivot] = matrix[pivot], matrix[j]
            determinant = -determinant
        value = matrix[j][j]
        determinant = determinant * value % prime
        inverse = pow(value, -1, prime)
        for i in range(j + 1, len(matrix)):
            scale = matrix[i][j] * inverse % prime
            if scale:
                for k in range(j, len(matrix)):
                    matrix[i][k] = (matrix[i][k] - scale * matrix[j][k]) % prime
    return determinant


def main():
    known = json.loads(Path(__file__).with_suffix(".json").read_text(encoding="utf-8"))
    prime = known["prime"]
    assert BO.s.isprime(prime)
    columns, rows = curvature_and_raising_matrix()
    assert len(columns) == known["raw_weight_columns"] == 88
    labels = [tuple(row) for row in known["selected_equation_labels"]]
    assert len(set(labels)) == len(columns)
    matrix = []
    for label in labels:
        values = [Fraction(rows[label].get(column, 0)) for column in range(len(columns))]
        matrix.append([value.numerator * pow(value.denominator, -1, prime) % prime
                       for value in values])
    determinant = determinant_modular(matrix, prime)
    assert determinant == known["minor_mod_prime"] != 0
    print(f"All 88 type-(5,5) potential columns are independent. Minor {determinant} mod {prime}.")


if __name__ == "__main__":
    main()
