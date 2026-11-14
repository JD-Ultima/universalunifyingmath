"""Verify the exact H6 curvature source and positive invariant pencils."""

from fractions import Fraction as F
from itertools import combinations_with_replacement
import json
from pathlib import Path

import sympy as s

import BC_affine_hooks as BC
import BG_h6_curvature as BG
import BO_h6_full_kernel as BO
import BP_h6_cartan_grade9 as BP
import BQ_h6_invariant_forms as BQ
import BS_h6_partial_value as BS
import CH_h4_second_order_stress as CH


HERE = Path(__file__).resolve().parent
METRIC = tuple(BC.METRIC)
DETERMINANT = CH.gram((1 << 14) - 1)


def inner(a, b):
    return sum((a[key] * b.get(key, 0) / (METRIC[key[0]] * CH.gram(key[1]))
                for key in a), F(0))


def source_form(basis, covariant_operations):
    scaled_duals = [BC.hodge_scaled(t) for t in basis]
    operators = [CH.polarization(t, covariant_operations) for t in scaled_duals]
    entries = [[CH.evaluate(op, t) for t in scaled_duals] for op in operators]
    n = len(basis)
    source = s.Matrix(n, n, lambda i, j: DETERMINANT * (entries[i][j] + entries[j][i]) / 2)
    norm = s.Matrix(n, n, lambda i, j: inner(basis[i], basis[j]))
    assert norm == DETERMINANT * s.Matrix(n, n, lambda i, j: inner(scaled_duals[i], scaled_duals[j]))
    return source, norm


def verify_types(covariant_operations, records):
    source = json.loads((HERE / "BO_h6_full_kernel_witnesses.json").read_text(encoding="utf-8"))
    moments = json.loads((HERE / "BQ_h6_invariant_forms.json").read_text(encoding="utf-8"))
    exact_source = {tuple(row["weight"]): row for row in records}
    bq = {tuple(row["weight"]): row for row in moments["blocks"]}
    transform = BP.split_to_physical_covectors()
    dimension = 0
    for entry in source["extremal_vectors"]:
        weight = tuple(entry["weight"])
        record = exact_source.pop(weight)
        seeds = [BQ.primitive(BP.transport(v, transform)[0]) for v in entry["vectors"]]
        computed, norm = source_form(seeds, covariant_operations)
        assert norm == s.Matrix(record["gram"]) == s.Matrix(bq[weight]["gram"])
        assert computed == s.Matrix(record["moment"])
        assert (computed - F(19, 100) * norm).is_positive_definite
        l5 = s.Matrix(bq[weight]["forms"]["5"])
        l9 = 384 * s.Matrix(bq[weight]["forms"]["9"])
        pencil = computed + F(2303, 25000) * l5 - F(8647, 50000) * l9
        assert (pencil - F(5193, 25000) * norm).is_positive_definite
        dimension += bq[weight]["dimension"] * len(seeds)
        print(f"{weight}: {len(seeds)} copies, both exact shifted forms positive.", flush=True)
    assert not exact_source and dimension == 6834


def verify_box(covariant_operations, record):
    package, _ = BS.read_gzip("partial_fields.json.gz")
    basis = [BS.dict_tensor(t) for t in package["fields"]]
    H, _ = source_form(basis, covariant_operations)
    assert H == s.Matrix(record["matrix"]) and H[0, 0] == F(record["center"])
    r = F(record["radius"])
    assert r == F(1, 10**8)
    error = sum(abs(F(H[i, j])) * (1 if i == j else 2) * r ** ((i != 0) + (j != 0))
                for i, j in list(combinations_with_replacement(range(len(basis)), 2))[1:])
    assert error == F(record["variation_bound"])
    assert H[0, 0] - error > F(167097264, 10**8)
    assert (H[0, 0] - error) / 3 > F(55699088, 10**8)
    print("BS box: exact curvature source > 1.67097264. Gradient square must exceed 0.55699088.",
          flush=True)


def main():
    data = json.loads((HERE / "CJ_h6_stress_certificate.json").read_text(encoding="utf-8"))
    curvature = CH.covariant_operations(BG.curvature_model())
    verify_box(curvature, data["bs_box"])
    verify_types(curvature, data["blocks"])
    print("All 15 H6 invariant stress blocks and the entire BS box verified.")


if __name__ == "__main__":
    main()
