"""Check the exact tensor and Lie algebra of the H6 analytic local lift."""

from fractions import Fraction
from itertools import combinations
import json
from pathlib import Path

import sympy as sp

import BO_h6_full_kernel as BO
from CF_h6_fourth_jet import verify_real_part


def commutator(a, b):
    return a * b - b * a


def tangent_matrices(metric):
    matrices = []
    for i in range(5):
        for j in range(i, 5):
            if (i, j) in ((0, 1), (2, 3), (4, 4)):
                continue
            unit = BO.matrix_units(i, j)
            matrices.append((unit if i == j else unit + unit.T) * metric)
    for i, j in [(0, 1), (2, 3)]:
        matrices.append((BO.matrix_units(i, j) + BO.matrix_units(j, i)
                         - 4 * BO.matrix_units(4, 4)) * metric)
    return matrices


def main():
    weights, inverse, curvature, roots = BO.build_model()
    metric = sp.zeros(5)
    metric[0, 1] = metric[1, 0] = metric[2, 3] = metric[3, 2] = 2
    metric[4, 4] = 1
    tangent = tangent_matrices(metric)
    compact = [
        (BO.matrix_units(i, j) - BO.matrix_units(j, i)) * metric
        for i, j in combinations(range(5), 2)
    ]
    plane = (1, 3)
    outside = (0, 2, 4)
    off_plane = [(i, j) for i in outside for j in plane]

    def preserves_plane(matrix):
        return all(matrix[i, j] == 0 for i, j in off_plane)

    def plane_trace(matrix):
        return sum(matrix[i, i] for i in plane)

    leaves = sp.Matrix([[matrix[i, j] for matrix in tangent] for i, j in off_plane])
    assert leaves.rank() == 5
    leaf_axes = [i for i, matrix in enumerate(tangent) if preserves_plane(matrix)]
    transverse_axes = [i for i in range(14) if i not in leaf_axes]
    leaf = [tangent[i] for i in leaf_axes]
    conditions = sp.Matrix([[matrix[i, j] for matrix in compact] for i, j in off_plane]
                           + [[plane_trace(matrix) for matrix in compact]])
    stabilizer = [
        sum((x * matrix for x, matrix in zip(column, compact)), sp.zeros(5))
        for column in conditions.nullspace()
    ]
    coefficient = (BO.matrix_units(1, 2) - BO.matrix_units(3, 0)) / 4
    assert all(commutator(matrix, coefficient) == sp.zeros(5) for matrix in stabilizer)
    for a, b in (*combinations(leaf, 2), *combinations(stabilizer, 2)):
        bracket = commutator(a, b)
        assert bracket.T * metric == -metric * bracket
        assert preserves_plane(bracket) and plane_trace(bracket) == 0
    for a in stabilizer:
        for b in leaf:
            bracket = commutator(a, b)
            assert bracket.T * metric == metric * bracket
            assert bracket.trace() == 0 and preserves_plane(bracket)

    volume = sum(1 << axis for axis in transverse_axes)
    tensor = {}
    for i in leaf_axes:
        component = commutator(tangent[i], coefficient) / 2
        for a, matrix in enumerate(tangent):
            value = Fraction(sp.trace(matrix * component)) * BO.wedge_sign(i, volume)
            BO.add(tensor, (a, volume | (1 << i)), value)
    residual = {}
    for (a, blade), value in tensor.items():
        assert tuple(weights[a][j] + sum(weights[i][j] for i in BO.axes(blade))
                     for j in range(2)) == (5, 5)
        for label, scale in BO.constraint_entries(a, blade, inverse, curvature):
            BO.add(residual, label, value * scale)
    assert tensor and not residual
    assert all(not BO.root_action(tensor, root) for root in roots)
    norm = verify_real_part(tensor)
    known = json.loads(Path(__file__).with_suffix(".json").read_text(encoding="utf-8"))
    assert leaf_axes == known["D_indices"] and transverse_axes == known["transverse_indices"]
    assert len(stabilizer) == known["kc_dimension"]
    assert [[a, blade, str(value)] for (a, blade), value in sorted(tensor.items())] == known["B_covariant_entries"]
    assert norm == Fraction(known["compact_real_value_norm_squared"]) == Fraction(3, 131072)
    print(f"Exact leaf bracket and hook checks passed. Real base squared norm {norm}.")


if __name__ == "__main__":
    main()
