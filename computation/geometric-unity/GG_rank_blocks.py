"""New Ricci-zero rank certificates, separate from all nonzero-Ricci data.

For alpha=+/-1, three t nodes certify degree-two polynomial identities. The
pure-mixed signs are separate exact cases. No observer-engine call or cache.
"""
import json
import pathlib

import sympy as sp
from sympy.polys.matrices import DomainMatrix
import Z_cubic_rank_one as Z

LAYOUT = (2, 3, 5, 7)
NODES = (-1, 0, 1)
PAIRS = tuple((i, j) for i in range(4) for j in range(i, 4))
RAW = tuple(1 << i for i in range(4)) + tuple((1 << i) ^ (1 << j) for i, j in PAIRS)
CHAR = tuple(min(c, c ^ 15) for c in RAW)
COLCHAR = CHAR * 3
ROWCHAR = tuple(min(RAW[i] ^ RAW[j] ^ RAW[k], (RAW[i] ^ RAW[j] ^ RAW[k]) ^ 15) for i, j, k in Z.META)
CERTIFICATES = pathlib.Path(__file__).with_name("GG_rank_certificates.json")


def rational_matrix(rows):
    return DomainMatrix.from_list_sympy(len(rows), len(rows[0]), rows).convert_to(sp.QQ)


def kernel():
    K = sp.zeros(42, 4)
    for i, di in enumerate(LAYOUT):
        K[14 + i, i], K[28 + i, i] = sp.Rational(24 * di, 5), 1
    return K


def sparse(rows, cols, entries):
    M = sp.zeros(rows, cols)
    for i, j, value in entries:
        assert 0 <= i < rows and 0 <= j < cols and M[i, j] == 0
        M[i, j] = sp.Rational(value)
    return M


def verify_normalized_case(sign, matrices, emit=print):
    assert sign in (1, -1) and len(matrices) == 3
    assert all(len(A) == 1274 and all(len(row) == 42 for row in A) for A in matrices)
    data = json.loads(CERTIFICATES.read_text(encoding="utf-8"))
    assert data["schema_version"] == 1 and tuple(data["layout"]) == LAYOUT and tuple(data["nodes"]) == NODES
    cases = [case for case in data["cases"] if case["alpha"] == sign]
    assert len(cases) == 1
    case, count, K = cases[0], 0, kernel()

    def check(label, condition):
        nonlocal count
        count += 1
        assert bool(condition), f"GG alpha={sign}: {label}"

    check("degree-bounded even-sign support and universal K",
          all(not A[r][c] for A in matrices for r in Z.INVARIANT for c in range(42) if ROWCHAR[r] != COLCHAR[c])
          and all(sp.Matrix([A[r] for r in Z.INVARIANT]) * K == sp.zeros(722, 4) for A in matrices))
    assert sorted(item["character"] for item in case["blocks"]) == list(range(8))
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
            assert char in (1, 2, 4, 7)
            L = sparse(2, len(rows), item["constant_combinations"])
            C = sparse(2, 3, item["constant_output"])
            check(f"horizontal character {char}: rank-two constant output", C.rank() == 2 and all(L * M == C for M in node_blocks))

    completion = data["completion"]
    rows = [Z.META.index((i, 5, 5)) for i in range(4)]
    assert completion["vertical"] == 5 and completion["rows"] == rows
    c, q = sp.Rational(completion["positive_constant"]), sp.Rational(completion["quadratic_coefficient"])
    factors = list(map(sp.Rational, completion["prefactors"]))
    check("J5 is strictly invertible for all real t", c > 0 and q > 0 and len(factors) == 4 and all(x > 0 for x in factors)
          and all(sp.Matrix([A[r] for r in rows]) * K == -sp.diag(*factors) * (c + q * t * t) for t, A in zip(NODES, matrices)))
    j7 = [Z.META.index((i, 7, 7)) for i in range(4)]
    check("the old J7 completion is zero at Ricci zero", all(sp.Matrix([A[r] for r in j7]) * K == sp.zeros(4) for A in matrices))
    emit(f"  normalized alpha={sign}: {count} exact rank checks; invariant rank 38, J5-completed rank 42 for every real t.")
    return count


def verify_pure_case(sign, A, b, emit=print):
    assert sign in (1, -1)
    M, K = rational_matrix([A[r] for r in Z.INVARIANT]), kernel()
    assert M.rank() == 38 and M.to_Matrix() * K == sp.zeros(722, 4)
    rows = [Z.META.index((i, 5, 5)) for i in range(4)]
    assert (sp.Matrix([A[r] for r in rows]) * K).det() != 0
    assert rational_matrix(A).rank() == 42 and rational_matrix([row + [rhs] for row, rhs in zip(A, b)]).rank() == 43
    emit(f"  pure-mixed beta={sign}: 3 exact checks; invariant/full/augmented ranks 38/42/43, J5 invertible.")
    return 3
