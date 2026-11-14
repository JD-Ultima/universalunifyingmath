"""Exact grade-five exclusion on the two lowest H4 curvature isotypes."""

import argparse
from collections import defaultdict
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations, combinations_with_replacement
import gzip
import hashlib
import json
from math import prod
from pathlib import Path
from time import perf_counter

import numpy as np


HERE = Path(__file__).resolve().parent
I = range(14)
G = [F(2)] * 4 + [F(8)] * 7 + [F(24), F(48), F(16, 5)]
PAIRS = list(combinations_with_replacement(range(24), 2))
LOCATION = {p: j for j, p in enumerate(PAIRS)}
BLADES = (1 << 14) - 1


def bits(m):
    return tuple(i for i in I if m >> i & 1)


def before(m, i):
    return (m & ((1 << i) - 1)).bit_count()


def sign(v):
    return 1 - 2 * (v % 2)


def check_phi():
    def offdiag(i, j, value):
        a = np.zeros((5, 5), dtype=np.int64)
        a[i, j] = a[j, i] = value
        return a

    q = [offdiag(i, 4, -5) for i in range(4)]
    q.extend(offdiag(i, j, -10) for i, j in combinations(range(4), 2))
    q.extend(-10 * np.diag([1] * i + [-i] + [0] * (4 - i)) for i in range(1, 4))
    q.append(np.diag([-2, -2, -2, -2, 8]))
    assert all(int(np.trace(u @ v)) == (25 * G[a] if a == b else 0) for a, u in enumerate(q) for b, v in enumerate(q))
    return q


def maps():
    q = check_phi()
    R, L = [], []
    for i, j in combinations(I, 2):
        m = q[i] @ q[j] - q[j] @ q[i]
        y = [m @ p - p @ m for p in q]
        R.append((i, j, [[(d, F(int(np.trace(y[d] @ q[a])), 2500 * G[a])) for d in I if np.trace(y[d] @ q[a])]
                          for a in I]))
    for i, j in combinations(range(5), 2):
        a = np.zeros((5, 5), dtype=np.int64)
        a[i, j], a[j, i] = 1, -1
        y = [a @ p - p @ a for p in q]
        L.append([[(d, -F(int(np.trace(y[d] @ q[b])), 25 * G[b])) for d in I if np.trace(y[d] @ q[b])]
                  for b in I])
    return R, L


@lru_cache(None)
def metric(m):
    return prod(G[i] for i in bits(m))


def verify_linear(t, typ, curvature, action):
    alt, inside, curv = defaultdict(F), defaultdict(F), defaultdict(F)
    for (d, m), c in t.items():
        if m >> d & 1:
            inside[m ^ (1 << d)] += F(sign(before(m, d)) * c, G[d])
        else:
            alt[m | (1 << d)] += sign(before(m, d)) * c
        rest = BLADES ^ m
        v = F(c * sign(sum(before(rest, i) for i in bits(m))), metric(m))
        for a, b, rows in curvature:
            include = (1 << a) | (1 << b)
            if not rest & include:
                for j, coeff in rows[d]:
                    curv[j, rest | include] += v * coeff * sign(before(rest, a) + before(rest, b))
    assert all(not any(q.values()) for q in (alt, inside, curv))

    def apply(t, slot):
        v = defaultdict(F)
        for (a, m), c in t.items():
            for d, k in slot[a]:
                v[d, m] += c * k
            for b in bits(m):
                tail = m ^ (1 << b)
                for d, k in slot[b]:
                    if not (tail >> d & 1):
                        v[a, tail | (1 << d)] += c * k * sign(before(tail, b) + before(tail, d))
        return {k: v for k, v in v.items() if v}

    out = defaultdict(F)
    for slot in action:
        for j, v in apply(apply(t, slot), slot).items():
            out[j] -= v
    a, b = typ
    scalar = a * (a + 3) + b * (b + 1)
    assert {k: v for k, v in out.items() if v} == {k: scalar * v for k, v in t.items()}


def independent(tensors, p):
    pivots = {}
    for a in tensors:
        w = {k: int(v) % p for k, v in a.items() if int(v) % p}
        while w:
            key = min(w)
            if key not in pivots:
                scale = pow(w[key], -1, p)
                pivots[key] = {k: (scale * v) % p for k, v in w.items()}
                break
            s = w[key]
            for k, v in pivots[key].items():
                z = (w.get(k, 0) - s * v) % p
                if z:
                    w[k] = z
                elif k in w:
                    del w[k]
    return len(pivots)


@lru_cache(None)
def clifford_sign(a, b):
    return sign(sum(before(b, i) for i in bits(a)))


def coefficient(a, target, U, p):
    m = np.zeros(300, dtype=np.int64)
    den = lambda v: int(v.numerator) * pow(int(v.denominator), -1, p) % p
    common = den(-F(2, 3) * G[a])
    for b, c in combinations(I, 2):
        if a in (b, c):
            continue
        triple = (1 << a) | (1 << b) | (1 << c)
        parity = clifford_sign(1 << a, 1 << b) * clifford_sign((1 << a) ^ (1 << b), 1 << c)
        for left, entries in U[b].items():
            right = triple ^ left ^ target
            if right not in U[c]:
                continue
            total = sum(sg * clifford_sign(x, y) * clifford_sign(x ^ y, z)
                        for sg, (x, y, z) in ((1, (triple, left, right)), (-1, (triple, right, left)),
                                               (1, (left, triple, right)), (1, (left, right, triple)),
                                               (-1, (right, triple, left)), (-1, (right, left, triple))))
            if total:
                factor = common * parity * total * den(F(1, metric((triple & left) | (triple & right) | (left & right)))) % p
                for i, s in entries:
                    for j, t in U[c][right]:
                        location = LOCATION[min(i, j), max(i, j)]
                        m[location] = (int(m[location]) + factor * (s % p) * (t % p)) % p
    return m


@lru_cache(None)
def multiply(a, b):
    return F(clifford_sign(a, b), metric(a & b))


def direct(a, receiver, tensor):
    """Native Fraction ordered products, without shared contraction prefactors."""
    field = [{m: v for (i, m), v in tensor.items() if a1 == i} for a1 in I]
    result = F(0)
    for b, c in combinations(I, 2):
        if a in (b, c):
            continue
        g, h, i = 1 << a, 1 << b, 1 << c
        t = g ^ h ^ i
        factor = -F(2, 3) * G[a] * multiply(g, h) * multiply(g ^ h, i)
        for m, v in field[b].items():
            n = receiver ^ t ^ m
            if n not in field[c]:
                continue
            terms = ((t, m, n), (m, t, n), (m, n, t), (t, n, m), (n, t, m), (n, m, t))
            total = sum((-1 if j >= 3 else 1) * multiply(x, y) * multiply(x ^ y, z) for j, (x, y, z) in enumerate(terms))
            result += factor * total * v * field[c][n]
    return result


def eliminate(M, p):
    A = np.array(M, dtype=np.int64) % p
    rank = 0
    pivots = []
    for j in range(A.shape[1]):
        hits = np.flatnonzero(A[rank:, j])
        if not len(hits):
            continue
        i = rank + int(hits[0])
        A[[rank, i]] = A[[i, rank]]
        A[rank, j:] = A[rank, j:] * pow(int(A[rank, j]), -1, p) % p
        t = np.flatnonzero(A[:, j])
        t = t[t != rank]
        A[t, j:] = (A[t, j:] - A[t, j, None] * A[rank, None, j:]) % p
        pivots.append(j)
        rank += 1
        if rank == A.shape[0]:
            break
    return rank, pivots, A


def quadratic_product(a, b):
    return np.array([a[i] * b[j] + (a[j] * b[i] if i != j else 0) for i, j in combinations_with_replacement(range(3), 2)])


def rank_one(M, p, minors):
    rank, pivots, A = eliminate(M, p)
    assert rank == 297
    rest = sorted(set(range(300)) - set(pivots))
    K = np.zeros((300, 3), dtype=np.int64)
    K[rest] = np.eye(3, dtype=np.int64)
    K[pivots] = -A[:, rest] % p
    assert not np.any(M @ K % p)
    C = []
    for i, j in minors:
        off = K[LOCATION[i, j]]
        C.append((quadratic_product(K[LOCATION[i, i]], K[LOCATION[j, j]]) - quadratic_product(off, off)) % p)
    assert eliminate(C, p)[0] == 6
    return rest, C


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--exhaustive", action="store_true")
    args = p.parse_args()
    t = perf_counter()
    raw = gzip.decompress((HERE / "h4_low_witness.json.gz").read_bytes())
    assert hashlib.sha256(raw).hexdigest() == "1d4669c61c3e79e2e95eba16a6dcd930cd962593772d6324a4fec8b3fa39eecc"
    data = json.loads(raw)
    types = [(1, 0)] * 10 + [(2, 0)] * 14
    assert [tuple(v["type"]) for v in data["tensors"]] == types
    T = [{(a, m): v for a, m, v in c["entries"]} for c in data["tensors"]]
    C, generators = maps()
    for v, lab in zip(T, types):
        verify_linear(v, lab, C, generators)
    assert independent(T, 65537) == 24
    U = [{} for _ in I]
    for j, field in enumerate(T):
        for (a, m), v in field.items():
            U[a].setdefault(m, []).append((j, v))
    receivers = data["receivers"]
    assert len(receivers) == len({tuple(r) for r in receivers}) == 297 and all(m.bit_count() == 5 for a, m in receivers)
    primes = (65537, 65521) if args.exhaustive else (65537,)
    for prime in primes:
        M = np.array([coefficient(a, m, U, prime) for a, m in receivers])
        for j, (a, m) in enumerate(receivers[:5]):
            v = direct(a, m, T[j])
            assert int(v.numerator) * pow(int(v.denominator), -1, prime) % prime == M[j, LOCATION[j, j]]
        for j in (0, 3, 11):
            a, m = receivers[j]
            b, c = j, j + 1
            z = defaultdict(int, T[b])
            for i, v in T[c].items():
                z[i] += v
            mixed = direct(a, m, z) - direct(a, m, T[b]) - direct(a, m, T[c])
            assert mixed.numerator * pow(mixed.denominator, -1, prime) % prime == M[j, LOCATION[b, c]]
        free, minors = rank_one(M, prime, data["independent_minor_indices"])
        print(f"prime {prime}: grade-five rank 297, free columns {free}, quadratic minor rank 6")
    print(f"24 independent exact hook-curvature vectors and grade-five exclusion verified in {perf_counter() - t:.2f}s")


if __name__ == "__main__":
    main()
