"""Independent full-column-rank certificate for first-hook integrability.

Uses a nonorthogonal matrix basis, direct coordinate extraction and finite-field
Gaussian elimination. Full modular column rank proves full rational/real rank:
an integer maximal minor nonzero modulo p cannot be zero over the integers.
This implication would not determine a deficient rational rank.
"""
from itertools import combinations
from math import isqrt
from time import perf_counter

N, PRIME = 14, 104729
PAIRS = list(combinations(range(N), 2))


def multiply(A, B):
    return [[sum(A[i][k] * B[k][j] for k in range(5))
             for j in range(5)] for i in range(5)]


def bracket(A, B):
    AB, BA = multiply(A, B), multiply(B, A)
    return [[AB[i][j] - BA[i][j] for j in range(5)] for i in range(5)]


def inner(A, B):
    return sum(A[i][j] * B[j][i] for i in range(5) for j in range(5))


def put(row, col, value):
    row[col] = row.get(col, 0) + value
    if not row[col]:
        del row[col]


def main():
    started = perf_counter()
    assert all(PRIME % d for d in range(2, isqrt(PRIME) + 1))
    basis, chars = [], []
    for i, j in combinations(range(5), 2):
        A = [[0] * 5 for _ in range(5)]
        A[i][j] = A[j][i] = 1
        basis.append(A)
        chars.append((1 << i) ^ (1 << j))
    # Nonorthogonal diagonals differ from the rational certificate basis.
    for i in range(4):
        A = [[0] * 5 for _ in range(5)]
        A[i][i], A[4][4] = 1, -1
        basis.append(A)
        chars.append(0)
    G = [[inner(A, B) for B in basis] for A in basis]
    inverse10 = [[5 * int(i == j) if i < 10 and j < 10 else
                  10 * int(i == j) - 2 if i >= 10 and j >= 10 else 0
                  for j in range(N)] for i in range(N)]
    assert all(sum(G[i][k] * inverse10[k][j] for k in range(N))
               == 10 * int(i == j) for i in range(N) for j in range(N))
    R = {}
    for i, j in PAIRS:
        K = bracket(basis[i], basis[j])
        columns = []
        for A in basis:
            C = bracket(K, A)
            assert all(C[a][b] == C[b][a] for a in range(5) for b in range(5))
            assert sum(C[a][a] for a in range(5)) == 0
            columns.append([-C[a][b] for a, b in combinations(range(5), 2)]
                           + [-C[a][a] for a in range(4)])
        R[i, j] = columns

    def curvature(i, j, output, argument):
        if i == j:
            return 0
        return R[i, j][argument][output] if i < j else -R[j, i][argument][output]

    assert all(sum(G[d][b] * curvature(i, j, d, a)
                   + G[a][d] * curvature(i, j, d, b) for d in range(N)) == 0
               for i, j in PAIRS for a in range(N) for b in range(N))
    assert all(curvature(i, j, a, k) + curvature(j, k, a, i)
               + curvature(k, i, a, j) == 0
               for i, j, k in combinations(range(N), 3) for a in range(N))
    assert all(sum(curvature(i, j, i, k) for i in range(N)) == -5 * G[j][k]
               for j in range(N) for k in range(N))
    variables = [(a, b, c) for a in range(N) for b, c in PAIRS]
    position = {v: i for i, v in enumerate(variables)}
    signatures = [chars[a] ^ chars[b] ^ chars[c] for a, b, c in variables]
    sizes = {char: signatures.count(char) for char in set(signatures)}
    pivots = {char: {} for char in sizes}

    def consume(row):
        if not row:
            return
        signatures_used = {signatures[col] for col in row}
        assert len(signatures_used) == 1
        char = signatures_used.pop()
        if len(pivots[char]) == sizes[char]:
            return
        work = {col: value % PRIME for col, value in row.items() if value % PRIME}
        while work:
            pivot = min(work)
            if pivot not in pivots[char]:
                inv = pow(work[pivot], -1, PRIME)
                pivots[char][pivot] = {col: value * inv % PRIME for col, value in work.items()}
                return
            factor = work[pivot]
            for col, value in pivots[char][pivot].items():
                remainder = (work.get(col, 0) - factor * value) % PRIME
                if remainder:
                    work[col] = remainder
                else:
                    work.pop(col, None)

    def component(row, a, b, c, value):
        if b < c:
            put(row, position[a, b, c], value)
        elif c < b:
            put(row, position[a, c, b], -value)

    for a, b, c in combinations(range(N), 3):
        consume({position[a, b, c]: 1, position[b, a, c]: -1, position[c, a, b]: 1})
    for j in range(N):
        row = {}
        for a in range(N):
            for b in range(N):
                if inverse10[a][b]:
                    component(row, a, b, j, inverse10[a][b])
        consume(row)
    assert sum(len(p) for p in pivots.values()) == 378
    for indices in combinations(range(N), 4):
        for a in range(N):
            row = {}
            for first, second in combinations(range(4), 2):
                rest = [j for j in range(4) if j not in (first, second)]
                i, j, k, ell = [indices[q] for q in (first, second, *rest)]
                sign = (-1) ** (first + second - 1)
                for d in range(N):
                    value = -sign * curvature(i, j, d, a)
                    if value:
                        component(row, d, k, ell, value)
            consume(row)
    print('Nonorthogonal Sym0(R5) basis; direct matrix-entry coordinates.')
    print('Exact integer curvature identities: metric skewness, Bianchi, Ric=-5*g.')
    print('Prime verified by trial division:', PRIME)
    print('character | raw columns | full modular rank')
    for char in sorted(sizes):
        assert len(pivots[char]) == sizes[char]
        print(f'{char:05b}', sizes[char], len(pivots[char]))
    assert sum(sizes.values()) == 1274
    print('Combined Alt/Ctr/curvature matrix has full column rank 1274 modulo p.')
    print('Therefore rational/real kernel is zero; curvature is injective on H2.')
    print('This is a maximal-minor certificate, not a deficient-rank inference.')
    print(f'Runtime {perf_counter() - started:.3f} seconds.')


if __name__ == '__main__':
    main()
