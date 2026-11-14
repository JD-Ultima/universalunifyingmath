"""Verify CA: every pointwise H4 curvature-kernel value has a closed lift.

The BX geometry builders are imported from the same directory.
"""

import argparse
from collections import defaultdict
from fractions import Fraction as F
from functools import lru_cache
import gzip
import hashlib
from itertools import combinations, product
import json
from pathlib import Path
from time import perf_counter

import numpy as np
import sympy as s


HERE = Path(__file__).resolve().parent
from verify_h4_an import METRIC as G, V, before, parity, phi, connection, operator
from verify_h4_dual_curvature import geometry, hodge, inventory, block
from verify_h4_character import dimension


def split_frame():
    K = s.zeros(5)
    K[0, 1] = K[1, 0] = K[2, 3] = K[3, 2] = 2
    K[4, 4] = 1
    units = [s.eye(5)[:, i] for i in range(5)]
    e = units
    frame = []
    for i in range(5):
        for j in range(i, 5):
            if (i, j) not in ((0, 1), (2, 3), (4, 4)):
                z = e[i] * e[j].T
                frame.append((z + z.T if i != j else z) * K)
    for i, j in ((0, 1), (2, 3)):
        frame.append((e[i] * e[j].T + e[j] * e[i].T - 4 * e[4] * e[4].T) * K)
    assert all(x.T * K == K * x and s.trace(x) == 0 for x in frame)
    return K, frame


def root_actions(K, P, weights):
    gram = s.Matrix(14, 14, lambda i, j: s.trace(P[i] * P[j]))
    inv = gram.inv()
    actions = []
    table = [((2, 0, 1), (1, 3, -1)), ((3, 4, 1), (4, 2, -2))]
    for data, root in zip(table, ((1, -1), (0, 1))):
        r = s.zeros(5)
        for i, j, v in data:
            r[i, j] = v
        assert r.T * K + K * r == s.zeros(5)
        t = [inv * s.Matrix([s.trace(p * (r * q - q * r)) for p in P]) for q in P]
        slots = []
        for a in V:
            row = [(b, -F(t[b][a])) for b in V if t[b][a]]
            assert all(weights[b] == tuple(weights[a][i] + root[i] for i in (0, 1)) for b, _ in row)
            slots.append(row)
        actions.append(slots)
    return actions


def action(values, slot):
    out = defaultdict(F)
    for (a, m), value in values.items():
        for b, c in slot[a]:
            out[b, m] += c * value
        for i in V:
            if not m >> i & 1:
                continue
            reduced = m ^ (1 << i)
            for j, c in slot[i]:
                if not reduced >> j & 1:
                    out[a, reduced | (1 << j)] += c * value * parity(before(reduced, i) + before(reduced, j))
    return {key: c for key, c in out.items() if c}


def dual_frame(frame):
    C = s.Matrix([[1, 1, 0, 0, 0], [s.I, -s.I, 0, 0, 0], [0, 0, 1, 1, 0],
                  [0, 0, s.I, -s.I, 0], [0, 0, 0, 0, 1]])
    v = [s.Matrix(x.tolist()) / 5 for x in phi()]
    Z = [C * x * C.inv() for x in frame]
    S = s.Matrix(14, 14, lambda a, b: s.trace(v[a] * Z[b]) / G[a])
    target = s.Matrix(14, 14, lambda i, j: s.trace(frame[i] * frame[j]))
    assert (S.T * s.diag(*G) * S - target).applyfunc(s.simplify) == s.zeros(14)
    return S.inv()


def mod(x, p):
    x = F(x)
    return x.numerator * pow(x.denominator, -1, p) % p


def minor_pairs(dual, vectors, values, p, imag):
    T = np.array([[(mod(s.re(x), p) - imag * mod(s.im(x), p)) % p for x in row] for row in dual.tolist()], dtype=np.int64)
    active = [np.flatnonzero(row).tolist() for row in T]

    @lru_cache(None)
    def blade_power(m):
        src = [i for i in V if m >> i & 1]
        candidates = {tuple(sorted(v)) for v in product(*(active[i] for i in src)) if len(set(v)) == 4}
        answer = []
        for c in sorted(candidates):
            minor = int(s.Matrix(T[np.asarray(src)[:, None], c].tolist()).det()) % p
            if minor:
                answer.append((sum(1 << x for x in c), minor))
        return answer

    output = np.empty((len(vectors), len(values)), dtype=np.int64)
    divided = [{(a, m): mod(F(c, G[a] * np.prod([G[j] for j in V if m >> j & 1])), p)
                for (a, m), c in v.items()} for v in values]
    for i, z in enumerate(vectors):
        transformed = defaultdict(int)
        for (a, m), c in z.items():
            for b in active[a]:
                for n, coeff in blade_power(m):
                    transformed[b, n] += c * int(T[a, b]) * coeff
        for j, q in enumerate(divided):
            output[i, j] = sum((x % p) * q.get(key, 0) for key, x in transformed.items()) % p
    return output.tolist()


def verify_AN(witnesses):
    lookup = {}
    tables = {}
    geom = connection()
    for ob in witnesses:
        c, i = label = tuple(ob["identifier"])
        if c not in tables:
            A, cols, _, _ = operator(geom, c)
            tables[c] = A, {item: j for j, item in enumerate(cols)}
        A, pos = tables[c]
        entries = {(a, b): v for a, b, v in ob["entries"]}
        assert len(entries) == len(ob["entries"]) and entries and label not in lookup
        vec = [0] * A.shape[1]
        for k, value in entries.items():
            vec[pos[k]] = value
        # Python integer arithmetic avoids overflow at arbitrary coefficients.
        assert all(sum(int(row[j]) * vec[j] for j in np.flatnonzero(row)) == 0 for row in A)
        lookup[label] = entries
    return lookup


def run(second):
    t = perf_counter()
    raw = gzip.decompress((HERE / "h4_lift_witness.json.gz").read_bytes())
    assert hashlib.sha256(raw).hexdigest() == "bc3862b6781c5285389356336ebff7befd30b6f4e2d8e4fe66f53a07dd0195fa"
    data = json.loads(raw)
    AN = verify_AN(data["AN"])
    assert len(AN) == 28 and set(k[0] for k in AN) == {0, 5, 17, 20, 29}
    K, P = split_frame()
    weights, inv, cov = geometry()
    roots = root_actions(K, P, weights)
    columns, rows = inventory(weights)
    star = hodge(inv)
    transport = dual_frame(P)
    result = []
    p, imaginary = data["prime"], data["gaussian_unit"]
    assert s.isprime(p) and imaginary ** 2 % p == p - 1
    bx = json.loads((HERE / "h4_character.json").read_text(encoding="utf-8"))
    known = {tuple(x["weight"]): x["copies"] for x in bx["components"]}
    exact_highests = 0
    total = 0
    for group in data["groups"]:
        weight = tuple(group["type"])
        c = group["multiplicity"]
        assert c == known[weight] == len(group["vectors"]) and dimension(*weight) == group["per_copy"]
        M = block(columns[weight], rows.get(weight, []), inv, cov, star)
        z = []
        for item in group["vectors"]:
            v = {(a, m): n for a, m, n in item}
            assert len(v) == len(item) and set(v) <= set(columns[weight])
            q = [v.get(key, 0) for key in columns[weight]]
            assert all(sum(int(row[j]) * q[j] for j in np.flatnonzero(row)) == 0 for row in M)
            assert all(not action(v, slot) for slot in roots)
            z.append(v)
        picks = [AN[tuple(label)] for label in group["pivot_sources"]]
        Q = minor_pairs(transport, z, picks, p, imaginary)
        assert Q == group["minor"]
        determinant = int(s.Matrix(Q).det()) % p
        assert determinant == group["determinant"] and determinant != 0
        if second:
            other = 65521
            i = next(t for t in range(other) if t * t % other == other - 1)
            alternate = minor_pairs(transport, z, picks, other, i)
            det_alternate = int(s.Matrix(alternate).det()) % other
            assert det_alternate != 0
        else:
            det_alternate = None
        total += c * dimension(*weight)
        exact_highests += c
        result.append(weight)
        print(f"{weight}: multiplicity {c}, nonzero determinant {determinant}, alternate {det_alternate}")
    assert len(set(result)) == len(known) == 22 and exact_highests == 57 and total == bx["dimension"] == 10388
    print(f"Verified full analytic H4 lift span: {total}; {exact_highests} exact highests and {len(AN)} exact closed values; seconds {perf_counter() - t:.2f}")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--exhaustive", action="store_true", help="repeat the pairing reductions at a second Gaussian prime")
    run(p.parse_args().exhaustive)
