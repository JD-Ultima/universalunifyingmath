"""Rebuild and verify the flat-observer H4 AN closure certificates.

Python 3.10+ and NumPy are required. Default runs build all exact row matrices,
check their hashes and all triangular independence patterns, and test three
exact null vectors per block. --exhaustive checks all 6,643 null vectors and
both independent prime-rank reductions. Run from any working directory.
"""

import argparse
from collections import defaultdict
from fractions import Fraction
from itertools import combinations
from pathlib import Path
import gzip
import hashlib
import json
from math import lcm
from time import perf_counter

import numpy as np


F = Fraction
V = range(14)
METRIC = [F(2)] * 4 + [F(8)] * 7 + [F(24), F(48), F(16, 5)]
INPUT_MASKS = [sum(1 << i for i in p) for p in combinations(V, 4)]
CHARS = [(1 << a) | 16 for a in range(4)] + [(1 << a) | (1 << b) for a, b in combinations(range(4), 2)] + [0] * 4
HERE = Path(__file__).resolve().parent


def before(mask, a):
    return (mask & ((1 << a) - 1)).bit_count()


def parity(i):
    return 1 - 2 * (i % 2)


def character(mask):
    x = 0
    for i, a in enumerate(CHARS):
        if mask >> i & 1:
            x ^= a
    return x


def phi():
    """Return the 14 integer matrices Q=5Φ in the public BC order."""
    def sym(a, b, q):
        m = np.zeros((5, 5), dtype=np.int64)
        m[a, b] = m[b, a] = q
        return m

    return [sym(i, 4, -5) for i in range(4)] + [sym(i, j, -10) for i, j in combinations(range(4), 2)] + [
        -10 * np.diag(([1] * j + [-j] + [0] * (4 - j))) for j in range(1, 4)
    ] + [np.diag([-2, -2, -2, -2, 8])]


def connection():
    q = phi()
    assert all(F(int(np.trace(A @ B)), 25) == (METRIC[i] if i == j else 0)
               for i, A in enumerate(q) for j, B in enumerate(q))
    result = []
    for x in q:
        k = np.triu(x, 1) - np.tril(x, -1)
        action = [k @ a - a @ k for a in q]
        values = [[-F(int(np.trace(action[d] @ q[s])), 250 * METRIC[s]) for s in V] for d in V]
        result.append([(d, s, v) for d in V for s, v in enumerate(values[d]) if v])
    return result


def operator(table, ch):
    columns = [(a, m) for a in V for m in INPUT_MASKS if CHARS[a] ^ character(m) == ch]
    denom = lcm(*(c.denominator for row in table for _, _, c in row), *(F(1, g).denominator for g in METRIC),
                *((c / g).denominator for row in table for _, _, c in row for g in METRIC))
    entries = defaultdict(int)
    rows = {}
    for j, (a, m) in enumerate(columns):
        def add(key, value):
            i = rows.setdefault(key, len(rows))
            assert value.denominator == 1
            entries[i, j] += int(value)

        if m >> a & 1:
            add((1, 0, m ^ (1 << a)), F(denom * parity(before(m, a)), METRIC[a]))
        else:
            add((0, 0, m | (1 << a)), F(denom * parity(before(m, a))))
        for i, row in enumerate(table):
            factor = -F(denom, METRIC[i])
            for d, s, c in row:
                if s == a and m >> i & 1:
                    add((2, d, m ^ (1 << i)), factor * c * parity(before(m, i)))
                if not m >> s & 1:
                    continue
                tail = m ^ (1 << s)
                if tail >> d & 1:
                    continue
                target = tail | (1 << d)
                if target >> i & 1:
                    add((2, a, target ^ (1 << i)), factor * c * parity(before(m, s) + before(tail, d) + before(target, i)))
    A = np.zeros((len(rows), len(columns)), dtype=np.int64)
    for (i, j), value in entries.items():
        A[i, j] = value
    return A, columns, list(rows), denom


def modular_rank(a, p):
    a = a.copy() % p
    rank = 0
    for j in range(a.shape[1]):
        hits = np.flatnonzero(a[rank:, j])
        if not len(hits):
            continue
        k = rank + int(hits[0])
        a[[rank, k]] = a[[k, rank]]
        a[rank, j:] = a[rank, j:] * pow(int(a[rank, j]), -1, p) % p
        nz = np.flatnonzero(a[rank + 1:, j]) + rank + 1
        a[nz, j:] = (a[nz, j:] - a[nz, j, None] * a[rank, None, j:]) % p
        rank += 1
        if rank == min(a.shape):
            break
    return rank


def verify(vectors, matrix, expected, exhaustive):
    free = expected["free_columns"]
    assert len(free) == len(vectors) == expected["exact_dimension"]
    free_set = set(free)
    assert len(free) == len(free_set)
    A = [[(int(c), int(row[c])) for c in np.flatnonzero(row)] for row in matrix]
    test = set(range(len(vectors))) if exhaustive else {0, len(vectors) // 2, len(vectors) - 1}
    for j, rows in enumerate(vectors):
        d = dict(rows)
        assert len(d) == len(rows) and all(0 <= u < matrix.shape[1] and v for u, v in d.items())
        assert set(d) & free_set == {free[j]}
        if j in test:
            assert all(sum(a * d.get(c, 0) for c, a in row) == 0 for row in A)
    return len(test)


def main():
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--exhaustive", action="store_true")
    opt = cli.parse_args()
    t = perf_counter()
    spec = json.loads((HERE / "h4_an_manifest.json").read_text(encoding="utf-8"))
    table = connection()
    dim = tested = 0
    for b in spec["blocks"]:
        A, cols, rows, denominator = operator(table, b["character"])
        digest = hashlib.sha256(A.astype("<i8", copy=False).tobytes()).hexdigest()
        assert A.shape == tuple(b["shape"]) and digest == b["operator_sha256"] and denominator == b["denominator"]
        raw = gzip.decompress((HERE / b["witness"]).read_bytes())
        assert hashlib.sha256(raw).hexdigest() == b["witness_sha256"]
        data = json.loads(raw)
        assert data["character"] == b["character"]
        tested += verify(data["columns"], A, data, opt.exhaustive)
        if opt.exhaustive:
            assert modular_rank(A, 65537) == modular_rank(A, 65521) == A.shape[1] - data["exact_dimension"]
        dim += data["exact_dimension"]
        print(f"character {b['character']:02d}: {A.shape[1]} inputs, {data['exact_dimension']} independent candidate vectors")
    assert len(spec["blocks"]) == 16 and dim == 6643
    mode = "verified exact dimension" if opt.exhaustive else "recorded dimension, sample verification only"
    print(f"{mode}: {dim}; exact null vectors checked: {tested}; seconds: {perf_counter() - t:.2f}")


if __name__ == "__main__":
    main()
