"""Verify two exact H4 grade-five exclusions using cubic multiples."""

import argparse
from fractions import Fraction
from functools import lru_cache
import gzip
from itertools import combinations_with_replacement, permutations, product
import json
import os
from pathlib import Path

os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"

import numpy as np
import sympy as sp

import BC_affine_hooks as BC
import BG_h6_curvature as BG
import CK_modular_clifford as CK

HERE = Path(__file__).resolve().parent
PRIME = CK.P
PANEL_WIDTH = 64
assert PANEL_WIDTH * (PRIME - 1) ** 2 + PRIME < 2**53


def panel_pivots(panel):
    reduced = panel.astype(np.int64)
    row_ids = np.arange(len(reduced))
    columns = []
    for column in range(reduced.shape[1]):
        rank = len(columns)
        hits = np.flatnonzero(reduced[rank:, column])
        if not len(hits):
            continue
        pivot = rank + int(hits[0])
        reduced[[rank, pivot]] = reduced[[pivot, rank]]
        row_ids[[rank, pivot]] = row_ids[[pivot, rank]]
        row = reduced[rank, column:]
        row[:] = row * pow(int(row[0]), -1, PRIME) % PRIME
        trailing = reduced[rank + 1 :, column:]
        trailing[:] = (trailing - trailing[:, :1] * row) % PRIME
        columns.append(column)
        if len(columns) == len(reduced):
            break
    return row_ids[: len(columns)], columns


def block_inverse(block):
    dimension = len(block)
    identity = np.eye(dimension, dtype=np.int64)
    reduced = np.c_[block.astype(np.int64), identity]
    determinant = 1
    for column in range(dimension):
        pivot = column + int(np.flatnonzero(reduced[column:, column])[0])
        if pivot != column:
            reduced[[column, pivot]] = reduced[[pivot, column]]
            determinant = -determinant
        value = int(reduced[column, column])
        determinant = determinant * value % PRIME
        reduced[column] = reduced[column] * pow(value, -1, PRIME) % PRIME
        others = np.r_[np.arange(column), np.arange(column + 1, dimension)]
        reduced[others] = (
            reduced[others] - reduced[others, column, None] * reduced[column]
        ) % PRIME
    inverse = reduced[:, dimension:]
    assert np.array_equal(block.astype(np.int64) @ inverse % PRIME, identity)
    return inverse.astype(float), determinant


def modular_rank(matrix):
    """Find a nonzero minor by exact Schur elimination modulo PRIME."""
    assert np.issubdtype(matrix.dtype, np.integer)
    assert np.min(matrix) >= 0 and np.max(matrix) < PRIME
    remaining = matrix.astype(float)
    row_ids = np.arange(len(remaining))
    offset = 0
    panels = []
    determinant = 1
    rank = 0
    while all(remaining.shape):
        width = min(PANEL_WIDTH, remaining.shape[1])
        rows, columns = panel_pivots(remaining[:, :width])
        keep = np.ones(len(remaining), dtype=bool)
        keep[rows] = False
        if len(rows):
            inverse, block_det = block_inverse(remaining[np.ix_(rows, columns)])
            multipliers = np.remainder(remaining[keep][:, columns] @ inverse, PRIME)
            eliminated = multipliers @ remaining[rows, width:]
            remaining = np.remainder(remaining[keep, width:] - eliminated, PRIME)
            panels.append({
                "rows": row_ids[rows].tolist(),
                "columns": [offset + column for column in columns],
                "determinant": block_det,
            })
            determinant = (determinant * block_det) % PRIME
        else:
            remaining = remaining[:, width:]
        # Every nonpivot column in this panel vanishes in the residual rows.
        rank += len(rows)
        offset += width
        row_ids = row_ids[keep]
    assert determinant
    return {"rank": rank, "panels": panels, "determinant_product_mod_prime": determinant}


def adjacent_permutation(index):
    phi = [sp.Matrix(matrix) / 5 for matrix in BG.integer_phi()]
    permutation = sp.eye(5)
    permutation.row_swap(index, index + 1)
    action = sp.Matrix(
        14, 14,
        lambda a, b: sp.trace(phi[a] * permutation * phi[b] * permutation.T)
        / sp.Rational(BC.METRIC[a]),
    )
    metric = sp.diag(*map(sp.Rational, BC.METRIC))
    assert action**2 == sp.eye(14) and action.T * metric * action == metric
    rows = [
        [(column, Fraction(value)) for column, value in enumerate(row) if value]
        for row in action.tolist()
    ]

    @lru_cache(None)
    def exterior(mask):
        output = {0: Fraction(1)}
        for axis in BC.bit_indices(mask):
            extended = {}
            for part, coefficient in output.items():
                for target, value in rows[axis]:
                    if not (part >> target) & 1:
                        sign = (-1) ** (part >> (target + 1)).bit_count()
                        BC.plus(extended, part | (1 << target), sign * coefficient * value)
            output = extended
        return tuple(output.items())

    def apply(tensor):
        output = {}
        for (axis, mask), coefficient in tensor.items():
            for target, factor in rows[axis]:
                for blade, value in exterior(mask):
                    BC.plus(output, (target, blade), coefficient * factor * value)
        return output

    return apply


def compact_casimir(tensor):
    output = {}
    for generator in range(10):
        scale = -1 if generator < 6 else -4
        for slot, value in BC.action(generator, BC.action(generator, tensor)).items():
            BC.plus(output, slot, scale * value)
    return output


def eigenvalue(weight):
    first, second = weight
    return first * (first + 3) + second * (second + 1)


@lru_cache(None)
def kostant_partition(a, b):
    # Counts nonnegative sums of (1,-1), (0,1), (1,0), (1,1).
    return sum(
        b + negative - positive >= 0
        for negative in range(a + 1)
        for positive in range(a - negative + 1)
    )


def weight_multiplicities(first, second):
    highest = (2 * first + 3, 2 * second + 1)
    weyl = [
        (signs[0] * highest[perm[0]], signs[1] * highest[perm[1]],
         signs[0] * signs[1] * (1 if perm == (0, 1) else -1))
        for perm in permutations((0, 1))
        for signs in product((-1, 1), repeat=2)
    ]
    weights = {}
    for a, b in product(range(-first, first + 1), repeat=2):
        multiplicity = sum(
            sign * kostant_partition((u - 2 * a - 3) // 2, (v - 2 * b - 1) // 2)
            for u, v, sign in weyl
        )
        assert multiplicity >= 0
        if multiplicity:
            weights[a, b] = multiplicity
    dimension = (2 * first + 3) * (2 * second + 1) * (first - second + 1) * (first + second + 2) // 6
    assert sum(weights.values()) == dimension
    return weights


@lru_cache(None)
def restriction(weight):
    # (Root of unity order, two rotation exponents, class size, sign).
    classes = [(1, 0, 0, 1, 1), (2, 1, 1, 10, -1), (2, 1, 0, 15, 1),
               (3, 1, 0, 20, 1), (6, 3, 1, 20, -1), (4, 2, 1, 30, -1),
               (5, 1, 2, 24, 1)]
    weights = weight_multiplicities(*weight)
    variable = sp.Symbol("z")
    values = []
    for order, a, b, size, sign in classes:
        polynomial = sum(
            count * variable ** ((a * i + b * j) % order)
            for (i, j), count in weights.items()
        )
        value = sp.rem(polynomial, sp.cyclotomic_poly(order, variable), variable)
        assert value.is_Integer
        values.append(int(value))
    totals = {
        character: sum(value * size * (sign if character == -1 else 1)
                       for value, (_, _, _, size, sign) in zip(values, classes))
        for character in (1, -1)
    }
    assert all(total % 120 == 0 and total >= 0 for total in totals.values())
    return values, {character: total // 120 for character, total in totals.items()}


def original_coefficient(tensor, axis, mask):
    field = [{} for _ in range(14)]
    for (a, blade), value in tensor.items():
        field[a][blade] = value
    return BC.quadratic_coefficient(field, BC.METRIC, axis, mask)


def verify_domain(domain, kernel, curvature, permutations):
    character = domain["character"]
    tensors = []
    groups = []
    claimed_weights = set()
    for group in domain["groups"]:
        weight = tuple(group["weight"])
        value = eigenvalue(weight)
        assert weight not in claimed_weights
        claimed_weights.add(weight)
        values, multiplicities = restriction(weight)
        copies = kernel[weight]
        dimension = copies * multiplicities[character]
        assert copies == group["copies"]
        assert dimension == group["dimension"] == len(group["tensors"])
        assert dimension == sum(
            count * restriction(other)[1][character]
            for other, count in kernel.items() if eigenvalue(other) == value
        )
        for entries in group["tensors"]:
            tensor = {(axis, mask): Fraction(coefficient) for axis, mask, coefficient in entries}
            assert len(tensor) == len(entries)
            assert all(0 <= axis < 14 and 0 < mask < 2**14 and mask.bit_count() == 4
                       for axis, mask in tensor)
            assert BC.hook_traces(tensor) == ({}, {})
            assert not BG.full_sparse_curvature(curvature, BC.hodge_scaled(tensor))
            assert compact_casimir(tensor) == {slot: value * entry for slot, entry in tensor.items()}
            signed = {slot: character * entry for slot, entry in tensor.items()}
            assert all(apply(tensor) == signed for apply in permutations)
            tensors.append(tensor)
        groups.append({"weight": weight, "copies": copies, "dimension": dimension,
                       "S5_character_values": values})

    dimension = len(tensors)
    assert dimension == domain["dimension"]
    slots = sorted(set().union(*(tensor.keys() for tensor in tensors)))
    basis = np.array([[CK.scalar(tensor.get(slot, 0)) for slot in slots]
                      for tensor in tensors], dtype=np.int64)
    basis_rank = modular_rank(basis)
    assert basis_rank["rank"] == dimension

    receivers = [(axis, mask, 0) for axis, mask in domain["receivers"]]
    assert len(set(receivers)) == len(receivers)
    assert all(0 <= axis < 14 and 0 < mask < 2**14 and mask.bit_count() == 5
               for axis, mask, _ in receivers)
    guards = 0

    def coefficients(tensor, audit):
        nonlocal guards
        result = CK.coefficients(tensor, receivers)
        if audit:
            for axis, mask, _ in receivers[:2] + receivers[-2:]:
                assert result[axis, mask] == CK.scalar(original_coefficient(tensor, axis, mask))
                guards += 1
        return np.array([result[axis, mask] for axis, mask, _ in receivers], dtype=np.int64)

    diagonals = [coefficients(tensor, i in (0, dimension - 1)) for i, tensor in enumerate(tensors)]
    pairs = list(combinations_with_replacement(range(dimension), 2))
    columns = []
    for i, j in pairs:
        if i == j:
            columns.append(diagonals[i])
        else:
            total = tensors[i].copy()
            for slot, value in tensors[j].items():
                BC.plus(total, slot, value)
            mixed = coefficients(total, (i, j) == (0, 1)) - diagonals[i] - diagonals[j]
            columns.append(mixed % PRIME)

    quadratics = np.array(columns).T
    quadratic_rank = modular_rank(quadratics)
    assert quadratic_rank["rank"] == len(receivers) == domain["quadratic_rank"]
    # The rank is exact also over Q, since all rows are already independent mod PRIME.
    cubics = list(combinations_with_replacement(range(dimension), 3))
    cubic_slots = {monomial: column for column, monomial in enumerate(cubics)}
    multiples = np.zeros((dimension * len(receivers), len(cubics)), dtype=np.int64)
    for k in range(dimension):
        indices = [cubic_slots[tuple(sorted((i, j, k)))] for i, j in pairs]
        assert len(set(indices)) == len(indices)
        multiples[k::dimension, indices] = quadratics
    cubic_rank = modular_rank(multiples)
    assert cubic_rank["rank"] == len(cubics)
    return {
        "character": character, "dimension": dimension, "types": groups,
        "basis_rank": basis_rank, "quadratic_rank": quadratic_rank, "cubic_rank": cubic_rank,
        "actual_grade_five_coefficients": quadratics.size, "rational_controls": guards,
        "quadratic_nullity": len(pairs) - len(receivers), "cubic_matrix_shape": multiples.shape,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--wide", action="store_true", help="Verify the 28- and 31-dimensional domains.")
    parser.add_argument("--sign-only", action="store_true")
    parser.add_argument("--receipt", type=Path)
    args = parser.parse_args()
    character = json.loads((HERE / "BX_data" / "h4_character.json").read_text(encoding="utf-8"))
    kernel = {tuple(component["weight"]): component["copies"] for component in character["components"]}
    filename = "CN_h4_s5_wide_witness.json.gz" if args.wide else "CN_h4_s5_witness.json.gz"
    with gzip.open(HERE / filename, "rt", encoding="utf-8") as source:
        data = json.load(source)
    assert data["kernel_character"] == [
        {"weight": component["weight"], "copies": component["copies"]}
        for component in character["components"]
    ]
    expected_dimensions = {-1: 28, 1: 31} if args.wide else {-1: 16, 1: 21}
    assert {domain["character"]: domain["dimension"] for domain in data["domains"]} == expected_dimensions
    assert {
        sign: sum(copies * restriction(weight)[1][sign] for weight, copies in kernel.items())
        for sign in (-1, 1)
    } == {-1: 110, 1: 68}

    curvature = BG.curvature_model()
    actions = [adjacent_permutation(index) for index in range(4)]
    records = []
    for domain in data["domains"]:
        if args.sign_only and domain["character"] == 1:
            continue
        record = verify_domain(domain, kernel, curvature, actions)
        records.append(record)
        print(json.dumps({key: value for key, value in record.items() if not key.endswith("_rank")}),
              flush=True)
    if args.receipt:
        args.receipt.write_text(json.dumps({"prime": PRIME, "records": records}, indent=2), encoding="utf-8")
    print("CN exact cubic certificates passed.", flush=True)


if __name__ == "__main__":
    main()
