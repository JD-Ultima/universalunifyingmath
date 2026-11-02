"""Exact rank certificates for CC's diagonal Weyl plane.

The main driver supplies the six freshly computed full coefficient matrices
at DEG2_NODES. This module never calls the curvature engine or uses caches.
CC's prior jet argument bounds every entry's total parameter degree by two,
so exact checks at these unisolvent nodes certify the stated polynomial
identities. Constant left inverses, two-row Bezout identities and the three
kernel restrictions prove rank uniformly, rather than inferring it from
the ranks of sampled matrices. The data encode rational row combinations
and exact polynomial identities.
"""
import json
import pathlib

import sympy as sp

import Z_cubic_rank_one as Z


DEG2_NODES = ((0, 0), (0, 1), (1, 0), (0, 2), (1, 1), (2, 0))
LAYOUTS = ((2, 3, 4, 5), (2, 3, 5, 7), (2, 4, 5, 9))
PAIRS = tuple((i, j) for i in range(4) for j in range(i, 4))
CHAR = tuple(1 << i for i in range(4)) + tuple((1 << i) ^ (1 << j) for i, j in PAIRS)
ROWCHAR = tuple(CHAR[i] ^ CHAR[j] ^ CHAR[k] for i, j, k in Z.META)
COLCHAR = CHAR*3
HORIZONTAL = (1, 2, 4, 8)
OFF_DIAGONAL = (3, 5, 6, 9, 10, 12)
A, B = sp.symbols('a b')
MONOMIALS = (sp.S.One, B, A, B*B, A*B, A*A)
VANDERMONDE = sp.Matrix([[mon.subs({A: a, B: b}) for mon in MONOMIALS] for a, b in DEG2_NODES])
INVERSE = VANDERMONDE.inv()
CERTIFICATE_PATH = pathlib.Path(__file__).with_name('CC_rank_certificates.json')


def polynomial(values):
    """Reconstruct within the hand-proved total-degree-two bound."""
    coefficients = INVERSE*sp.Matrix([sp.Rational(value) for value in values])
    return sp.expand(sum(coefficient*monomial for coefficient, monomial in zip(coefficients, MONOMIALS)))


def sparse_matrix(rows, cols, entries):
    result = sp.zeros(rows, cols)
    for i, j, value in entries:
        assert 0 <= i < rows and 0 <= j < cols and result[i, j] == 0
        result[i, j] = sp.Rational(value)
    return result


def block_rows(node_A, character, cols):
    """Identically nonzero rows, using the established degree bound."""
    return [row for row in Z.INVARIANT if ROWCHAR[row] == character
            and any(matrix[row][col] for matrix in node_A for col in cols)]


def submatrices(node_A, rows, cols):
    return [sp.Matrix([[matrix[row][col] for col in cols] for row in rows]) for matrix in node_A]


def kernel(layout):
    result = sp.zeros(42, 4)
    for i, value in enumerate(layout):
        result[14+i, i] = sp.Rational(24, 5)*value
        result[28+i, i] = 1
    return result


def load_case(rho, layout):
    payload = json.loads(CERTIFICATE_PATH.read_text(encoding='utf-8'))
    assert payload['schema_version'] == 1
    assert tuple(map(tuple, payload['degree2_nodes'])) == DEG2_NODES
    cases = [case for case in payload['cases'] if case['rho'] == rho and tuple(case['layout']) == layout]
    assert len(cases) == 1
    return cases[0]


def verify_rank_case(rho, layout, node_A, emit=print):
    """Verify rank38+completion at one signed layout from six full matrices.

    Returns proof-summary counts. Input matrices are in DEG2_NODES order,
    with the same 1274 rows and 42 unknown columns as U's relaxed system.
    The main driver remains responsible for the prior polynomial degree bound.
    """
    layout = tuple(layout)
    assert rho in (1, -1) and layout in LAYOUTS
    assert len(node_A) == 6 and all(len(matrix) == 1274 and all(len(row) == 42 for row in matrix)
                                   for matrix in node_A)
    record = load_case(rho, layout)
    checks = 0

    def check(label, condition):
        nonlocal checks
        checks += 1
        assert bool(condition), f'CC rank rho={rho}, layout={layout}: {label}'

    K = kernel(layout)
    check('reflection-character support and universal four-vector kernel',
          all(not matrix[row][col] for matrix in node_A for row in Z.INVARIANT for col in range(42)
              if ROWCHAR[row] != COLCHAR[col])
          and all(sp.Matrix([matrix[row] for row in Z.INVARIANT])*K == sp.zeros(722, 4) for matrix in node_A))

    trivial = record['trivial']
    cols = [col for col in range(42) if COLCHAR[col] == 0]
    rows = block_rows(node_A, 0, cols)
    assert trivial['rows'] == rows and trivial['columns'] == cols
    positions = {row: i for i, row in enumerate(rows)}
    left_inverse = sp.zeros(12, len(rows))
    for i, entries in enumerate(trivial['left_inverse']):
        for row, value in entries:
            left_inverse[i, positions[row]] = sp.Rational(value)
    check('trivial block has a constant left inverse',
          all(left_inverse*M == sp.eye(12) for M in submatrices(node_A, rows, cols)))

    assert sorted(item['character'] for item in record['blocks']) == sorted(HORIZONTAL+OFF_DIAGONAL)
    for item in record['blocks']:
        character = item['character']
        full_cols = [col for col in range(42) if COLCHAR[col] == character]
        cols = full_cols[:2] if character in HORIZONTAL else full_cols
        rows = block_rows(node_A, character, full_cols)
        assert item['cols'] == cols and item['rows'] == rows
        matrices = submatrices(node_A, rows, cols)
        if character in HORIZONTAL:
            left_inverse = sparse_matrix(2, len(rows), item['left_inverse'])
            check(f'horizontal character{character} has a constant reduced left inverse',
                  all(left_inverse*M == sp.eye(2) for M in matrices))
            continue

        L = sparse_matrix(2, len(rows), item['constant_combinations'])
        C0 = sparse_matrix(2, 3, item['constant_output'])
        ko = sp.Matrix([sp.Rational(value) for value in item['constant_kernel']])
        constant_ok = C0.rank() == 2 and ko[0] == 0 and ko[2] == 1 and C0*ko == sp.zeros(2, 1)
        constant_ok &= all(L*M == C0 for M in matrices)
        variable = A if item['variable'] == 'a' else B
        assert item['variable'] in ('a', 'b')
        u, v = 3*variable+rho, 6*variable-rho
        residual = []
        for row in item['residual_rows']:
            assert row in rows
            residual.append(polynomial([sum(sp.Rational(matrix[row][col])*ko[j] for j, col in enumerate(cols))
                                        for matrix in node_A]))
        cu, cv = map(sp.Rational, item['residual_scales'])
        assert cu and cv
        fu, fv = (20*rho-48*variable)/(27*rho**3*cu), (12*variable+7*rho)/(27*rho**3*cv)
        check(f'off-diagonal character{character} has a unit residual ideal',
              constant_ok and sp.expand(residual[0]-cu*u*u) == 0
              and sp.expand(residual[1]-cv*v*v) == 0
              and sp.expand(fu*residual[0]+fv*residual[1]-1) == 0)

    scalars = {5: rho+3*A, 6: rho+3*B, 7: rho-3*A-3*B}
    assert sorted(item['vertical'] for item in record['completion']) == [5, 6, 7]
    for item in record['completion']:
        vertical = item['vertical']
        rows = [Z.META.index((i, vertical, vertical)) for i in range(4)]
        assert item['rows'] == rows
        diagonal = list(map(sp.Rational, item['diagonal']))
        assert len(diagonal) == 4 and all(diagonal)
        expected = sp.diag(*diagonal)
        check(f'four-row completion at vertical{vertical}',
              all(sp.Matrix([matrix[row] for row in rows])*K
                  == expected*scalars[vertical].subs({A: a, B: b})**2
                  for (a, b), matrix in zip(DEG2_NODES, node_A)))
    check('the completion scalar factors have no simultaneous zero',
          sp.expand(sum(scalars.values())) == 3*rho and rho != 0)
    emit(f'  rank certificates rho={rho}, y={layout}: {checks} exact checks; invariant rank 38, completed rank 42, no exceptional parameter.')
    return dict(checks=checks, invariant_rank=38, completed_rank=42, kernel_dimension=4)
