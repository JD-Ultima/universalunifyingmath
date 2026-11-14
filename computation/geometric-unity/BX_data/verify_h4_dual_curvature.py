"""Verify the complete H4 pointwise hook and dual-curvature dimension."""

import argparse
from collections import defaultdict
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations, product
import json
from math import lcm
from time import perf_counter

from verify_h4_an import HERE, V, before, parity, modular_rank, np
import sympy as s
from sympy.polys.matrices import DomainMatrix


def geometry():
    K = s.zeros(5)
    K[0, 1] = K[1, 0] = K[2, 3] = K[3, 2] = 2
    K[4, 4] = 1
    E = [s.eye(5)[:, i] for i in range(5)]
    weight = [(1, 0), (-1, 0), (0, 1), (0, -1), (0, 0)]
    P, weights = [], []
    for i in range(5):
        for j in range(i, 5):
            if (i, j) in ((0, 1), (2, 3), (4, 4)):
                continue
            C = E[i] * E[j].T
            P.append((C + C.T if i != j else C) * K)
            weights.append(tuple(a + b for a, b in zip(weight[i], weight[j])))
    for i, j in ((0, 1), (2, 3)):
        P.append((E[i] * E[j].T + E[j] * E[i].T - 4 * E[4] * E[4].T) * K)
        weights.append((0, 0))
    gram = s.Matrix(14, 14, lambda i, j: s.trace(P[i] * P[j]))
    inv = gram.inv()
    cov = {}
    for i, j in combinations(V, 2):
        C = P[i] * P[j] - P[j] * P[i]
        for a, p in enumerate(P):
            A = (C * p - p * C) / 4
            coeff = inv * s.Matrix([s.trace(p * A) for p in P])
            for b, x in enumerate(coeff):
                if x:
                    assert all(weights[b][k] == weights[i][k] + weights[j][k] + weights[a][k] for k in range(2))
                    cov.setdefault((b, i, j), []).append((a, F(x)))
    return weights, [[F(c) for c in row] for row in inv.tolist()], cov


def inventory(weights):
    columns, rows = defaultdict(list), defaultdict(list)
    for degree, tag in ((3, 1), (4, -1), (5, 0), (12, 2)):
        for indices in combinations(V, degree):
            mask = sum(1 << i for i in indices)
            w = tuple(sum(weights[i][j] for i in indices) for j in (0, 1))
            if tag == -1 or tag == 2:
                for a in V:
                    v = tuple(weights[a][j] + w[j] for j in (0, 1))
                    if tag == -1:
                        columns[v].append((a, mask))
                    else:
                        rows[v].append((tag, a, mask))
            else:
                rows[w].append((tag, 0, mask))
    return columns, rows


def hodge(inverse):
    @lru_cache(None)
    def dual(mask):
        blades = tuple(a for a in V if mask >> a & 1)
        choices = [[b for b in V if inverse[b][a]] for a in blades]
        rows = sorted({tuple(sorted(p)) for p in product(*choices) if len(set(p)) == 4})
        output = []
        for inds in rows:
            det = F(s.Matrix([[inverse[i][j] for j in blades] for i in inds]).det())
            if det:
                rest = ((1 << 14) - 1) ^ sum(1 << k for k in inds)
                output.append((rest, det * parity(sum(before(rest, k) for k in inds))))
        return output

    return dual


def block(columns, rows, inv, cov, dual):
    labels = {key: i for i, key in enumerate(rows)}
    entries = defaultdict(F)
    for j, (a, m) in enumerate(columns):
        def add(key, value):
            entries[labels[key], j] += value

        if not m >> a & 1:
            add((0, 0, m | (1 << a)), parity(before(m, a)))
        for e in V:
            if m >> e & 1 and inv[a][e]:
                add((1, 0, m ^ (1 << e)), parity(before(m, e)) * inv[a][e])
        for h, value in dual(m):
            for i, k in combinations((c for c in V if not h >> c & 1), 2):
                out = h | (1 << i) | (1 << k)
                sign = parity(before(h, i) + before(h, k))
                for b, curvature in cov.get((a, i, k), ()):
                    add((2, b, out), sign * value * curvature)
    d = lcm(*(v.denominator for v in entries.values()))
    out = np.zeros((len(rows), len(columns)), dtype=np.int64)
    for (i, j), q in entries.items():
        out[i, j] = int(d * q)
    return out


def pivots(A):
    p = 65537
    a = A.copy() % p
    vals = np.arange(a.shape[0])
    rows, columns = [], []
    n = 0
    for j in range(a.shape[1]):
        nz = np.flatnonzero(a[n:, j])
        if not len(nz):
            continue
        i = n + int(nz[0])
        a[[n, i]] = a[[i, n]]
        vals[[n, i]] = vals[[i, n]]
        a[n, j:] = a[n, j:] * pow(int(a[n, j]), -1, p) % p
        idx = np.flatnonzero(a[n + 1:, j]) + n + 1
        a[idx, j:] = (a[idx, j:] - a[idx, j, None] * a[n, None, j:]) % p
        rows.append(int(vals[n]))
        columns.append(j)
        n += 1
    return rows, columns


def check_upper(A, rank, rs, cs):
    if rank == 0:
        return not np.any(A)
    if rank == min(A.shape):
        return True
    B = A[np.asarray(rs)[:, None], cs]
    C = A[rs]
    minor = DomainMatrix.from_list_sympy(*B.shape, B.tolist()).convert_to(s.ZZ)
    rhs = DomainMatrix.from_list_sympy(*C.shape, C.tolist()).convert_to(s.ZZ)
    x, denom = minor.solve_den(rhs)
    a = np.asarray(x.to_Matrix().tolist(), dtype=object)
    return np.array_equal(A[:, cs].astype(object) @ a, int(denom) * A.astype(object))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--exhaustive", action="store_true")
    opt = p.parse_args()
    t = perf_counter()
    rec = json.loads((HERE / "h4_dual_curvature_manifest.json").read_text())
    known = {tuple(v["weight"]): v for v in rec["blocks"]}
    weights, inv, cov = geometry()
    cols, rows = inventory(weights)
    dual = hodge(inv)
    dims = []
    for w, b in sorted(cols.items()):
        if w < (0, 0):
            continue
        a = block(b, rows[w], inv, cov, dual)
        assert a.shape == (known[w]["rows"], known[w]["columns"])
        rs, cs = pivots(a)
        rank = len(rs)
        assert rank == known[w]["exact_rank"] and modular_rank(a, 65521) == rank
        if opt.exhaustive or w in ((0, 0), (1, 1), (3, 0)):
            assert check_upper(a, rank, rs, cs)
        dims.append((1 if w == (0, 0) else 2) * (a.shape[1] - rank))
    assert sum(dims) == rec["exact_dimension"] == 10388
    print(f"{'verified exact dimension' if opt.exhaustive else 'modular dimension upper bound'}: {sum(dims)}; checked {len(dims)} conjugate weight blocks; seconds {perf_counter() - t:.2f}")


if __name__ == "__main__":
    main()
