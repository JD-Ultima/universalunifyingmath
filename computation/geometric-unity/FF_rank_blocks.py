"""Constant rank certificates for FF's two fixed electric/mixed directions.

The main supplies six fresh coefficient matrices in total-degree-two node order.
The prior hand degree bound turns exact node identities into polynomial ones.
All block certificates and J7 completion are constant rational identities; no
sample ranks, parameter searches or off-repository caches are used.
"""
import json
import pathlib

import sympy as sp
import Z_cubic_rank_one as Z

NODES = ((0, 0), (0, 1), (1, 0), (0, 2), (1, 1), (2, 0))
LAYOUTS = ((2, 3, 4, 5), (2, 3, 5, 7), (2, 4, 5, 9))
PAIRS = tuple((i, j) for i in range(4) for j in range(i, 4))
RAW = tuple(1 << i for i in range(4)) + tuple((1 << i) ^ (1 << j) for i, j in PAIRS)


def character(value):
    return min(value, value ^ 15)


CHAR = tuple(character(value) for value in RAW)
COLCHAR = CHAR * 3
ROWCHAR = tuple(character(RAW[a] ^ RAW[b] ^ RAW[d]) for a, b, d in Z.META)
CERTIFICATES = pathlib.Path(__file__).with_name("FF_rank_certificates.json")


def sparse(rows, cols, entries):
    M = sp.zeros(rows, cols)
    for i, j, value in entries:
        assert 0 <= i < rows and 0 <= j < cols and M[i, j] == 0
        M[i, j] = sp.Rational(value)
    return M


def verify_rank_case(rho, layout, matrices, emit=print):
    layout = tuple(layout)
    assert rho in (1, -1) and layout in LAYOUTS
    assert len(matrices) == 6 and all(len(A) == 1274 and all(len(row) == 42 for row in A) for A in matrices)
    data = json.loads(CERTIFICATES.read_text(encoding="utf-8"))
    assert data["schema_version"] == 1 and tuple(map(tuple, data["nodes"])) == NODES
    cases = [case for case in data["cases"] if case["rho"] == rho and tuple(case["layout"]) == layout]
    assert len(cases) == 1
    case, count = cases[0], 0

    def check(label, condition):
        nonlocal count
        count += 1
        assert bool(condition), f"FF rank rho={rho}, y={layout}: {label}"

    K = sp.zeros(42, 4)
    for i, di in enumerate(layout):
        K[14 + i, i] = sp.Rational(24 * di, 5)
        K[28 + i, i] = 1
    check("even-sign character support and universal four-vector kernel",
          all(not A[r][c] for A in matrices for r in Z.INVARIANT for c in range(42)
              if ROWCHAR[r] != COLCHAR[c])
          and all(sp.Matrix([A[r] for r in Z.INVARIANT]) * K == sp.zeros(722, 4) for A in matrices))

    assert sorted(item["character"] for item in case["blocks"]) == [0, 1, 2, 3, 4, 5, 6, 7]
    ranks = {0: 12, 1: 2, 2: 2, 3: 6, 4: 2, 5: 6, 6: 6, 7: 2}
    for item in case["blocks"]:
        char = item["character"]
        cols = [i for i in range(42) if COLCHAR[i] == char]
        rows = [r for r in Z.INVARIANT if ROWCHAR[r] == char and any(A[r][c] for A in matrices for c in cols)]
        assert item["rows"] == rows and item["columns"] == cols and item["constant_rank"] == ranks[char]
        node_blocks = [sp.Matrix([[A[r][c] for c in cols] for r in rows]) for A in matrices]
        if "left_inverse" in item:
            T = sparse(len(cols), len(rows), item["left_inverse"])
            check(f"character {char}: constant left inverse", all(T * M == sp.eye(len(cols)) for M in node_blocks))
        else:
            assert char in (1, 2, 4, 7) and len(cols) == 3
            L = sparse(2, len(rows), item["constant_combinations"])
            C = sparse(2, 3, item["constant_output"])
            check(f"horizontal character {char}: constant rank-two output",
                  C.rank() == 2 and all(L * M == C for M in node_blocks))

    completion = case["completion"]
    rows = [Z.META.index((i, 7, 7)) for i in range(4)]
    assert completion["vertical"] == 7 and completion["rows"] == rows
    diagonal = list(map(sp.Rational, completion["diagonal"]))
    check("J7 has a constant invertible restriction on K",
          len(diagonal) == 4 and all(diagonal)
          and all(sp.Matrix([A[r] for r in rows]) * K == sp.diag(*diagonal) for A in matrices))
    emit(f"  rank certificates rho={rho}, y={layout}: {count} exact checks; invariant rank 38, J7-completed rank 42 at every parameter.")
    return dict(checks=count, invariant_rank=38, completed_rank=42, kernel_dimension=4)
