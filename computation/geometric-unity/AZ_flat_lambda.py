"""Exact positive-parameter first-hook curvature certificate.

Standalone affine/Nomizu geometry and full polynomial parameter cover.
The hand model and field argument are in AZ-first-hook-classification.md.
No earlier helper or recorded result is imported; running writes no files.
"""
from collections import defaultdict
from itertools import combinations
from math import gcd, lcm
from time import perf_counter
import sympy as s
from sympy.polys.matrices import DomainMatrix
mu = s.Symbol('mu', real=True)
N = 14
PAIRS = tuple(combinations(range(N), 2))

def model():
    sym, chars, names = ([], [1 << i for i in range(4)], [f'H{i}' for i in range(4)])
    for i, j in combinations(range(4), 2):
        u = s.zeros(4)
        u[i, j] = u[j, i] = 1
        sym.append(u)
        chars.append(1 << i ^ 1 << j)
        names.append(f'O{i}{j}')
    for k in range(1, 4):
        sym.append(s.diag(*[1] * k + [-k] + [0] * (3 - k)))
        chars.append(0)
        names.append(f'D{k}')
    sym.append(s.eye(4))
    chars.append(0)
    names.append('I')
    norms = [(u * u).trace() for u in sym]
    metric = [s.Integer(2)] * 4 + [4 * n for n in norms[:-1]] + [16 / mu]

    def coords(u):
        return s.Matrix([(u * v).trace() / n for v, n in zip(sym, norms)])
    L = []
    for i in range(4):
        mat = s.zeros(N)
        for j in range(4):
            v = s.zeros(4)
            v[i, j] += s.Rational(1, 4)
            v[j, i] += s.Rational(1, 4)
            if i == j:
                v += (mu - 1) * s.eye(4) / 8
            mat[4:, j] = coords(v)
        for a, u in enumerate(sym):
            mat[:4, 4 + a] = -u[:, i]
        L.append(mat)
    G = s.diag(*metric)
    assert all((s.cancel(v) == 0 for mat in L for v in G * mat + mat.T * G))
    # Torsion is checked against the reductive p-bracket in this moving frame.
    for i, j in PAIRS:
        torsion = (L[i][:, j] if i < 4 else s.zeros(N, 1)) - (L[j][:, i] if j < 4 else s.zeros(N, 1))
        expected_bracket = s.zeros(N, 1)
        if i < 4 <= j:
            expected_bracket[:4, 0] = -sym[j - 4][:, i]
        assert torsion == expected_bracket
    R = {}
    for i, j in PAIRS:
        if j < 4:
            value = L[i] * L[j] - L[j] * L[i]
        elif i < 4:
            value = sum((sym[j - 4][b, i] * L[b] for b in range(4)), s.zeros(N))
        else:
            K = sym[i - 4] * sym[j - 4] - sym[j - 4] * sym[i - 4]
            value = s.zeros(N)
            value[:4, :4] = -K
            for a, u in enumerate(sym):
                value[4:, 4 + a] = -coords(K * u - u * K)
        R[i, j] = value.applyfunc(s.expand)
        assert all((s.degree(v, mu) <= 1 for v in value if v))
        assert all((s.cancel(v) == 0 for v in G * value + value.T * G))

    def rv(i, j):
        return R[i, j] if i < j else -R[j, i] if i > j else s.zeros(N)
    assert all((rv(i, j)[:, k] + rv(j, k)[:, i] + rv(k, i)[:, j] == s.zeros(N, 1) for i, j, k in combinations(range(N), 3)))
    ric = s.Matrix(N, N, lambda j, k: sum((rv(i, j)[i, k] for i in range(N)))).applyfunc(s.expand)
    expected_ric = s.zeros(N)
    expected_ric[:4, :4] = -mu * s.eye(4) / 2
    for a, u in enumerate(sym):
        for b, v in enumerate(sym):
            expected_ric[4 + a, 4 + b] = -5 * (u * v).trace() + u.trace() * v.trace()
    assert ric == expected_ric
    assert ric.subs(mu, 5) == -s.Rational(5, 4) * G.subs(mu, 5)
    Phi = []
    for a in range(N):
        value = s.zeros(5)
        if a < 4:
            value[a, 4] = value[4, a] = -1
        else:
            u = sym[a - 4]
            value[:4, :4] = -2 * u + 2 * u.trace() * s.eye(4) / 5
            value[4, 4] = 2 * u.trace() / 5
        Phi.append(value)
    assert all(((Phi[a] * Phi[b]).trace() == G.subs(mu, 5)[a, b] for a in range(N) for b in range(N)))
    for i, j in PAIRS:
        K = Phi[i] * Phi[j] - Phi[j] * Phi[i]
        for a in range(N):
            lhs = sum((R[i, j][d, a].subs(mu, 5) * Phi[d] for d in range(N)), s.zeros(5))
            assert lhs == -(K * Phi[a] - Phi[a] * K) / 4
    return (metric, chars, names, R, ric)

def system(metric, chars, R):
    variables = [(a, b, c) for a in range(N) for b, c in PAIRS]
    positions = {v: j for j, v in enumerate(variables)}
    cols, hook, curv = (defaultdict(list), defaultdict(list), defaultdict(list))
    for j, (a, b, c) in enumerate(variables):
        cols[chars[a] ^ chars[b] ^ chars[c]].append(j)

    def insert(group, row):
        row = {j: s.expand(v) for j, v in row.items() if v != 0}
        if not row:
            return
        ch = {chars[variables[j][0]] ^ chars[variables[j][1]] ^ chars[variables[j][2]] for j in row}
        assert len(ch) == 1
        assert all((s.degree(v, mu) <= 1 for v in row.values()))
        # Clear only CONSTANT rational factors; parameter roots are retained.
        den = lcm(*(int(c.q) for v in row.values() for c in s.Poly(v, mu).all_coeffs()))
        row = {j: s.expand(den * v) for j, v in row.items()}
        div = gcd(*(int(c) for v in row.values() for c in s.Poly(v, mu).all_coeffs()))
        row = {j: v / div for j, v in row.items()}
        group[ch.pop()].append(row)
    for a, b, c in combinations(range(N), 3):
        insert(hook, {positions[a, b, c]: 1, positions[b, a, c]: -1, positions[c, a, b]: 1})
    for j in range(N):
        row = {}
        for a in range(N):
            if a != j:
                row[positions[a, min(a, j), max(a, j)]] = (1 if a < j else -1) / metric[a]
        insert(hook, row)
    for i, j, k, ell in combinations(range(N), 4):
        shuffles = ((i, j, k, ell, 1), (i, k, j, ell, -1), (i, ell, j, k, 1), (j, k, i, ell, 1), (j, ell, i, k, -1), (k, ell, i, j, 1))
        for a in range(N):
            row = {}
            for u, v, b, c, sign in shuffles:
                for d in range(N):
                    value = -sign * R[u, v][d, a]
                    if value:
                        idx = positions[d, b, c]
                        row[idx] = row.get(idx, 0) + value
            insert(curv, row)
    return (cols, hook, curv)

def dm(rows, cols, parameter=None):
    loc = {j: k for k, j in enumerate(cols)}
    data = {i: {loc[j]: v if parameter is None else v.subs(mu, parameter) for j, v in row.items() if parameter is None or v.subs(mu, parameter) != 0} for i, row in enumerate(rows)}
    return DomainMatrix.from_dict_sympy(len(rows), len(cols), data)

def main():
    started = perf_counter()
    metric, chars, names, R, ric = model()
    print('Nomizu model: mu=1/(1-4lambda)>0, sym(vw)=(vw^T+wv^T)/2.')
    print('Basis:', names)
    print('Metric diagonal:', metric)
    print('Exact metric compatibility, zero torsion, first Bianchi and affine curvature degree verified.')
    print('mu=5: exact Ric=-5*g/4 and all 1274 curvature columns match the SPD5 Phi identity.')
    print('General Ricci diagonal:', list(ric.diagonal()))
    cols, hook, curv = system(metric, chars, R)
    print('Combined all-covariant Alt/Ctr/(-R^T wedge B) matrix is affine in mu.')
    print('character | rawcolumns | hookrank at5 | combinedrank at5')
    for ch in range(16):
        hr = len(dm(hook[ch], cols[ch], 5).convert_to(s.QQ).rref()[1])
        rank = len(dm(hook[ch] + curv[ch], cols[ch], 5).convert_to(s.QQ).rref()[1])
        assert rank == len(cols[ch])
        print(f'{ch:04b}', len(cols[ch]), hr, rank)
    print('Symbolic square determinants, one representative per sign-weight orbit:')
    # The two numerical points select rows; determinants are computed in QQ[mu].
    for ch in (0, 1, 3, 7, 15):
        rows = hook[ch] + curv[ch]
        at5 = dm(rows, cols[ch], 5).convert_to(s.QQ)
        selected = at5.transpose().rref()[1]
        assert len(selected) == len(cols[ch])
        pencil = dm(rows, cols[ch]).convert_to(s.QQ.poly_ring(mu))
        det = s.factor(pencil.extract(selected, range(len(cols[ch]))).det().as_expr())
        assert det.subs(mu, 5) != 0
        print('CHAR', f'{ch:04b}', 'SIZE', len(cols[ch]), 'DET', det)
        print('ROWS', tuple(selected))
        at2 = dm(rows, cols[ch], 2).convert_to(s.QQ)
        second = at2.transpose().rref()[1]
        assert len(second) == len(cols[ch])
        det2 = s.factor(pencil.extract(second, range(len(cols[ch]))).det().as_expr())
        assert det2.subs(mu, 2) != 0
        common = s.factor(s.gcd(det, det2))
        common_poly = s.Poly(common, mu, domain=s.QQ).monic()
        allowed = s.Poly(mu * (mu - 1) * (mu - 3) * (mu - 7) * (mu + 1), mu, domain=s.QQ)
        assert (allowed ** common_poly.degree()).rem(common_poly).is_zero
        print('SECOND DET', det2)
        print('SECOND ROWS', tuple(second))
        print('COMMON FACTOR', common)
    print('Exact rational exceptional-parameter ranks (no inference from samples):')
    # These are ALL remaining positive common roots, not a sample grid.
    for value in (1, 3, 7):
        total_kernel = 0
        for ch, multiplicity in ((0, 1), (1, 4), (3, 6), (7, 4), (15, 1)):
            rows = hook[ch] + curv[ch]
            mat = dm(rows, cols[ch], value).convert_to(s.QQ)
            _, pivots = mat.rref()
            kernel_dimension = len(cols[ch]) - len(pivots)
            assert kernel_dimension == 0
            total_kernel += multiplicity * kernel_dimension
            print('mu', value, 'CHAR', f'{ch:04b}', 'RAW RANK', len(pivots), 'KERNEL', kernel_dimension, 'MULTIPLICITY', multiplicity)
        print('mu', value, 'FULL H2 KERNEL DIMENSION', total_kernel)
    print('Only displayed exact polynomial certificates support parameter conclusions; no sampled-rank classification.')
    print(f'Runtime {perf_counter() - started:.3f} seconds.')
if __name__ == '__main__':
    main()
