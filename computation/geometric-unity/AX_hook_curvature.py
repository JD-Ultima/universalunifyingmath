"""Exact AX curvature-integrability kernel on SL(5,R)/SO(5).

Fresh rational construction; no repository or supporting helper imports.
The curvature is full ambient curvature on p=Sym^2_0(R^5), not vertical
fibre curvature. Overall nonzero curvature scale does not change this kernel.
"""
from collections import defaultdict
from fractions import Fraction as F
from itertools import combinations
from math import gcd, lcm
from time import perf_counter
import sympy as s
from sympy.polys.matrices import DomainMatrix

N = 14
PAIRS = list(combinations(range(N), 2))
TRIPLES = list(combinations(range(N), 3))
QUADS = list(combinations(range(N), 4))


def mm(a, b):
    return tuple(tuple(sum(a[i][k] * b[k][j] for k in range(5))
                       for j in range(5)) for i in range(5))


def comm(a, b):
    ab, ba = mm(a, b), mm(b, a)
    return tuple(tuple(ab[i][j] - ba[i][j] for j in range(5)) for i in range(5))


def dot(a, b):
    return sum(a[i][j] * b[j][i] for i in range(5) for j in range(5))


def plus(row, col, val):
    if val:
        row[col] = row.get(col, 0) + val
        if not row[col]:
            del row[col]


def primitive(row):
    den = lcm(*(v.denominator if isinstance(v, F) else 1 for v in row.values()))
    out = {j: int(v * den) for j, v in row.items() if v}
    div = gcd(*out.values())
    if div:
        out = {j: v // div for j, v in out.items()}
    if out and out[min(out)] < 0:
        out = {j: -v for j, v in out.items()}
    return out


def coordinates(matrix, basis, metric):
    return tuple(F(dot(matrix, u), g) for u, g in zip(basis, metric))


def build():
    basis, chars, names = [], [], []
    for i, j in combinations(range(5), 2):
        mat = [[0] * 5 for _ in range(5)]
        mat[i][j] = mat[j][i] = 1
        basis.append(tuple(map(tuple, mat)))
        chars.append((1 << i) ^ (1 << j))
        names.append(f'O{i}{j}')
    for k in range(1, 5):
        diag = [1] * k + [-k] + [0] * (4 - k)
        basis.append(tuple(tuple(diag[i] if i == j else 0 for j in range(5)) for i in range(5)))
        chars.append(0)
        names.append(f'D{k}')
    metric = tuple(dot(u, u) for u in basis)
    assert metric == (2,) * 10 + (2, 6, 12, 20)
    assert all(dot(u, v) == 0 for i, u in enumerate(basis) for v in basis[i+1:])

    # R_ij maps p to p by -[[e_i,e_j],.]. Column a is R(e_i,e_j)e_a.
    curvature = {}
    for i, j in PAIRS:
        bracket = comm(basis[i], basis[j])
        matrix = [[F(0) for _ in range(N)] for _ in range(N)]
        for a in range(N):
            raw = comm(bracket, basis[a])
            coeffs = coordinates(raw, basis, metric)
            for d, value in enumerate(coeffs):
                matrix[d][a] = -value
        assert all(metric[d] * matrix[d][a] + metric[a] * matrix[a][d] == 0
                   for a in range(N) for d in range(N))
        curvature[i, j] = matrix

    def R(i, j, d, a):
        return curvature[i, j][d][a] if i < j else -curvature[j, i][d][a] if j < i else F(0)

    assert all(R(i, j, d, k) + R(j, k, d, i) + R(k, i, d, j) == 0
               for i, j, k in TRIPLES for d in range(N))
    ric = [[sum(R(i, j, i, k) for i in range(N)) for k in range(N)] for j in range(N)]
    assert all(ric[j][k] == (-5 * metric[j] if j == k else 0)
               for j in range(N) for k in range(N))

    variables = [(a, b, c) for a in range(N) for b, c in PAIRS]
    pos = {v: j for j, v in enumerate(variables)}
    columns = defaultdict(list)
    for index, (a, b, c) in enumerate(variables):
        columns[chars[a] ^ chars[b] ^ chars[c]].append(index)

    hook_rows, curv_rows = defaultdict(list), defaultdict(list)
    constant_rows = defaultdict(list)

    def insert(groups, row):
        if row:
            row = primitive(row)
            rowchars = {chars[variables[j][0]] ^ chars[variables[j][1]] ^ chars[variables[j][2]] for j in row}
            assert len(rowchars) == 1
            groups[rowchars.pop()].append(row)

    # Alt B: B_a,bc-B_b,ac+B_c,ab=0 for a<b<c.
    for a, b, c in TRIPLES:
        insert(hook_rows, {pos[a, b, c]: 1, pos[b, a, c]: -1, pos[c, a, b]: 1})
    # Ctr B: sum_a g^{aa} (i_a B_a)_j=0.
    for j in range(N):
        row = {}
        for a in range(N):
            if a < j:
                row[pos[a, a, j]] = F(1, metric[a])
            elif a > j:
                row[pos[a, j, a]] = F(-1, metric[a])
        insert(hook_rows, row)

    # Curvature on the COEFFICIENT covector index is F_ij=-R_ij^T.
    # For i<j<k<l, (F wedge B)_a,ijkl has the six wedge-shuffle terms.
    def add_curvature_rows(curv, groups):
        for i, j, k, ell in QUADS:
            shuffles = ((i, j, k, ell, 1), (i, k, j, ell, -1),
                        (i, ell, j, k, 1), (j, k, i, ell, 1),
                        (j, ell, i, k, -1), (k, ell, i, j, 1))
            for a in range(N):
                row = {}
                for u, v, b, c, sign in shuffles:
                    for d in range(N):
                        plus(row, pos[d, b, c], -sign * curv(u, v, d, a))
                insert(groups, row)

    add_curvature_rows(R, curv_rows)

    def constant_R(i, j, d, a):
        # R(X,Y)Z=<Y,Z>X-<X,Z>Y, in the same non-unit orthogonal basis.
        return metric[j] * int(d == i and a == j) - metric[i] * int(d == j and a == i)

    add_curvature_rows(constant_R, constant_rows)
    return basis, metric, names, chars, variables, columns, hook_rows, curv_rows, constant_rows


def rational_kernel(rows, cols):
    index = {j: k for k, j in enumerate(cols)}
    data = {i: {index[j]: v for j, v in row.items()} for i, row in enumerate(rows)}
    dm = DomainMatrix.from_dict_sympy(len(rows), len(cols), data).convert_to(s.QQ)
    reduced, pivots = dm.rref()
    rd = reduced.to_sdm()
    free = [j for j in range(len(cols)) if j not in pivots]
    vectors = []
    for f in free:
        vector = {f: F(1)}
        for i, p in enumerate(pivots):
            value = rd.get(i, {}).get(f, 0)
            if value:
                value = s.Rational(value)
                vector[p] = F(-int(value.p), int(value.q))
        assert all(sum(F(v) * vector.get(index[j], 0) for j, v in row.items()) == 0 for row in rows)
        vectors.append({cols[j]: v for j, v in vector.items()})
    return len(pivots), vectors


def square_witness(rows, cols):
    index = {j: k for k, j in enumerate(cols)}
    data = {i: {index[j]: v for j, v in row.items()} for i, row in enumerate(rows)}
    dm = DomainMatrix.from_dict_sympy(len(rows), len(cols), data).convert_to(s.QQ)
    _, independent_rows = dm.transpose().rref()
    assert len(independent_rows) == len(cols)
    witness = dm.extract(independent_rows, range(len(cols)))
    determinant = witness.det()
    assert determinant != 0
    return tuple(independent_rows), determinant


def main():
    started = perf_counter()
    _, metric, names, chars, variables, columns, hook_rows, curv_rows, constant_rows = build()
    print('Full ambient p=Sym0^2(R5), rational orthogonal basis:', names)
    print('Metric diagonal:', metric)
    print('R(X,Y)Z=-[[X,Y],Z]; covector curvature=-R^T.')
    print('Bianchi, metric skewness and Ric=-5*g verified exactly.')
    print('Raw tensor variables:', len(variables), '; Alt/Ctr rows:', sum(map(len, hook_rows.values())))
    print('character | rawcols | hookrank | H2dim | combinedrank | curvature rank on H2 | kernel')
    total_hook = total_rank = total_kernel = 0
    full_kernel = []
    for char in sorted(columns):
        cols = columns[char]
        hr, hook_kernel = rational_kernel(hook_rows[char], cols)
        assert all(sum(F(v) * vector.get(j, 0) for j, v in row.items()) == 0
                   for row in constant_rows[char] for vector in hook_kernel)
        cr, vectors = rational_kernel(hook_rows[char] + curv_rows[char], cols)
        selected, determinant = square_witness(hook_rows[char] + curv_rows[char], cols)
        total_hook += hr
        total_rank += cr
        total_kernel += len(vectors)
        full_kernel.extend(vectors)
        print(f'{char:05b}', len(cols), hr, len(cols)-hr, cr, cr-hr, len(vectors))
        print('  exact square row witness (hook rows, then curvature rows):', selected)
        print('  exact square determinant:', determinant)
    assert total_hook == 378 and len(variables)-total_hook == 896
    assert len(full_kernel) == total_kernel == len(variables)-total_rank == 0
    assert total_rank == 1274
    print('EXACT H2 dimension:', 896)
    print('EXACT curvature rank on H2:', total_rank-total_hook)
    print('EXACT curvature kernel dimension:', total_kernel)
    print('Constant-sectional-curvature control annihilates the entire 896-dimensional rational H2 basis.')
    for i, vector in enumerate(full_kernel):
        print('KERNEL', i, repr([(variables[j], str(v)) for j, v in sorted(vector.items())]))
    print('This computes full linear integrability, not a quadratic solution or physical test.')
    print(f'Runtime {perf_counter()-started:.3f} seconds.')


if __name__ == '__main__':
    main()
