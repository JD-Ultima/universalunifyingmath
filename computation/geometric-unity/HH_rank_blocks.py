"""Small exact rank identities for all Ricci-zero chiral eigenvalues, p=+/-1.

Ten fresh full coefficient matrices are supplied by HH_ricci_zero.py. The
coefficient degree-two bound precedes interpolation. This helper checks only
constant row identities, a positive two-column residual and J5 completion;
there are no engine calls, stored node matrices or parameter rank samples.
"""
import json
import pathlib

import sympy as sp
import Z_cubic_rank_one as Z

NODES = ((0, 0, 0), (0, 0, 1), (0, 1, 0), (1, 0, 0),
         (0, 0, 2), (0, 1, 1), (0, 2, 0), (1, 0, 1), (1, 1, 0), (2, 0, 0))
LAYOUT = (2, 3, 5, 7)
PAIRS = tuple((i, j) for i in range(4) for j in range(i, 4))
RAW = tuple(1 << i for i in range(4)) + tuple((1 << i) ^ (1 << j) for i, j in PAIRS)
CHAR = tuple(min(x, x ^ 15) for x in RAW)
COLCHAR = CHAR * 3
ROWCHAR = tuple(min(RAW[a] ^ RAW[b] ^ RAW[d], (RAW[a] ^ RAW[b] ^ RAW[d]) ^ 15) for a, b, d in Z.META)
DATA = pathlib.Path(__file__).with_name('HH_rank_certificates.json')


def sparse(rows, cols, entries):
    M = sp.zeros(rows, cols)
    for i, j, value in entries:
        assert 0 <= i < rows and 0 <= j < cols and M[i, j] == 0
        M[i, j] = sp.Rational(value)
    return M


def verify_rank_case(sign, matrices, emit=print):
    """Verify rank38 + positive J5 completion for every real q,r,s at fixed p."""
    assert sign in (1, -1) and len(matrices) == 10
    assert all(len(A) == 1274 and all(len(row) == 42 for row in A) for A in matrices)
    data = json.loads(DATA.read_text(encoding='utf-8'))
    assert data['schema_version'] == 1 and tuple(data['layout']) == LAYOUT
    assert tuple(map(tuple, data['nodes'])) == NODES
    cases = [case for case in data['cases'] if case['p'] == sign]
    assert len(cases) == 1
    case, count = cases[0], 0
    q, r, s = sp.symbols('q r s', real=True)
    variables = (q, r, s)

    def check(label, condition):
        nonlocal count
        count += 1
        assert bool(condition), f'HH rank p={sign}: {label}'

    monomials = (sp.S.One, s, r, q, s*s, r*s, r*r, q*s, q*r, q*q)
    V = sp.Matrix([[m.subs(dict(zip(variables, node))) for m in monomials] for node in NODES])
    check('ten degree-two nodes are unisolvent', abs(V.det()) == 8)
    K = sp.zeros(42, 4)
    for i, di in enumerate(LAYOUT):
        K[14 + i, i] = sp.Rational(24*di, 5)
        K[28 + i, i] = 1
    check('universal K and baseline even-sign support',
          all(not A[row][col] for A in matrices for row in Z.INVARIANT for col in range(42)
              if ROWCHAR[row] != COLCHAR[col])
          and all(sp.Matrix([A[row] for row in Z.INVARIANT])*K == sp.zeros(722, 4) for A in matrices))
    assert sorted(block['character'] for block in case['blocks']) == list(range(8))
    for item in case['blocks']:
        char = item['character']
        cols = [i for i in range(42) if COLCHAR[i] == char]
        rows = [row for row in Z.INVARIANT if ROWCHAR[row] == char
                and any(A[row][c] for A in matrices for c in cols)]
        assert item['rows'] == rows and item['columns'] == cols
        node_blocks = [sp.Matrix([[A[row][c] for c in cols] for row in rows]) for A in matrices]
        if 'left_inverse' in item:
            assert char in (0, 5, 6) and len(cols) == (12 if char == 0 else 6)
            T = sparse(len(cols), len(rows), item['left_inverse'])
            check(f'character {char}: constant rational left inverse', all(T*M == sp.eye(len(cols)) for M in node_blocks))
        elif char in (1, 2, 4, 7):
            L = sparse(2, len(rows), item['constant_combinations'])
            C = sparse(2, 3, item['constant_output'])
            check(f'horizontal character {char}: constant rank-two output', C.rank() == 2 and all(L*M == C for M in node_blocks))
        else:
            assert char == 3 and item['constant_rank'] == 4 and tuple(item['kernel_shape']) == (6, 2)
            L = sparse(4, len(rows), item['constant_combinations'])
            C = sparse(4, 6, item['constant_output'])
            Q = sparse(6, 2, item['kernel'])
            check('character 3: constant rank-four output and exact two-vector kernel',
                  C.rank() == 4 and Q.rank() == 2 and C*Q == sp.zeros(4, 2)
                  and all(L*M == C for M in node_blocks))
            residual = item['positive_residual_certificate']
            selected = (776, 791, 958, 1043)
            divisors = tuple(map(sp.Rational, residual['divisors']))
            assert tuple(residual['rows']) == selected
            assert divisors == (sp.Rational(-27, 14000), sp.Rational(-3, 784),
                                sp.Rational(81, 245000), sp.Rational(9, 13720))
            h = 31*q*q + 22*q*s + 31*s*s + 29*(q+sign)**2 + 2*(q+sign)*(r+s) + 29*(r+s)**2
            sos = sp.Rational(35, 3) + 48*(q+sp.Rational(7*sign, 12))**2 + 20*s*s + 28*(r+s)**2 + 11*(q+s)**2 + (q+r+s+sign)**2
            check('the two-column residual scalar has a strict real SOS bound', sp.expand(h-sos) == 0)
            identities = []
            for node, M in zip(NODES, node_blocks):
                MQ = M*Q
                pieces = [MQ[rows.index(row), :]/divisor for row, divisor in zip(selected, divisors)]
                combined = sp.Matrix.vstack(pieces[0]+pieces[1], pieces[2]+pieces[3])
                identities.append(combined == sp.eye(2)*h.subs(dict(zip(variables, node))))
            check('four residual rows combine to h I2, hence the O6 block has full rank', all(identities))

    completion_rows = [Z.META.index((i, 5, 5)) for i in range(4)]
    diagonal = sp.diag(sp.Rational(1, 700), sp.Rational(3, 1400), sp.Rational(1, 280), sp.Rational(1, 200))
    scalar = 41*r*r + 58*sign*r + 41
    check('J5K is its exact degree-two scalar times a nonzero diagonal',
          all(sp.Matrix([A[row] for row in completion_rows])*K == -diagonal*scalar.subs(r, node[1])
              for node, A in zip(NODES, matrices)))
    check('the J5 scalar is strictly positive for all real r',
          sp.expand(scalar-(41*(r+sp.Rational(29*sign, 41))**2+sp.Rational(840, 41))) == 0)
    emit(f'  normalized p={sign}: {count} exact rank checks; invariant rank 38 and J5-completed rank 42 for every real q,r,s.')
    return dict(checks=count, invariant_rank=38, completed_rank=42, kernel_dimension=4)
