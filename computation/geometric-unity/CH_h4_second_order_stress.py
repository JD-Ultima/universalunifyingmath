"""Verify CH's exact curvature functional on the CB box and the complete BX kernel."""

import argparse
from collections import defaultdict
from fractions import Fraction as F
from functools import lru_cache
import gzip
import hashlib
from itertools import combinations, combinations_with_replacement, product
import json
from math import prod
from pathlib import Path
import sys
from time import perf_counter

import sympy as s

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE / "BX_data"))

import BC_affine_hooks as BC
import BG_h6_curvature as BG
import verify_h4_lift as CA

G = tuple(BC.METRIC)
V = range(14)


def read_packed(path, digest):
    raw = gzip.decompress(path.read_bytes())
    assert hashlib.sha256(raw).hexdigest() == digest
    return json.loads(raw)


@lru_cache(None)
def gram(m):
    return prod(G[i] for i in BC.bit_indices(m))


def interior(m, i):
    if m >> i & 1:
        n = m ^ (1 << i)
        return n, BC.sign(1 << i, n)


def covariant_operations(R):
    out = {}
    for (i, j), A in R.items():
        out[i, j] = tuple(tuple((b, -c) for b, c in enumerate(row) if c) for row in A)
    return out


def act(t, maps):
    out = {}
    for (a, m), c in t.items():
        for b, v in maps[a]:
            BC.plus(out, (b, m), c * v)
        for i in BC.bit_indices(m):
            n, e = interior(m, i)
            for b, v in maps[i]:
                if not n >> b & 1:
                    BC.plus(out, (a, n | (1 << b)), c * v * e * BC.sign(1 << b, n))
    return out


def polarization(t, R):
    """A coefficient functional u -> 2E_R(u,t)-I_R(u,t)."""
    result = defaultdict(F)
    for (i, j), operation in R.items():
        factor = 1 / (G[i] * G[j])
        for (a, m), coefficient in act(t, operation).items():
            for p, q, ordering in ((i, j, -1), (j, i, 1)):
                c = coefficient * factor * ordering
                if a == q:
                    result[p, m] += 2 * c / gram(m)
                if (hit := interior(m, q)) is None:
                    continue
                n, e = hit
                if not n >> p & 1:
                    new = n | (1 << p)
                    result[a, new] -= c * e * BC.sign(1 << p, n) / (G[a] * gram(n))
    return result


def evaluate(linear, t):
    return sum((c * linear.get(key, 0) for key, c in t.items()), F(0))


def explicit_calibration(R, cov):
    # A two-entry coordinate hook checks the ordered signs and exterior Gram.
    z = {(0, 30): F(1), (1, 29): F(1)}
    assert BC.hook_traces(z) == ({}, {})
    H = evaluate(polarization(z, cov), z)
    # Use the full ordered curvature sum and explicit interior pairs.
    total1 = total2 = F(0)
    for (i, j), M in R.items():
        for I, J, sgn in ((i, j, -1), (j, i, 1)):
            D = {(a, m): sgn * c for (a, m), c in act(z, cov[i, j]).items()}
            w = 1 / (G[I] * G[J])
            for (b, m), c in z.items():
                if b == I:
                    total1 += w * c * D.get((J, m), 0) / gram(m)
                r = interior(m, I)
                if r is None or r[0] >> J & 1:
                    continue
                n, e = r
                total2 += w * c * D.get((b, n | (1 << J)), 0) * e * BC.sign(1 << J, n) / (G[b] * gram(n))
    assert H == 2 * total1 - total2
    assert (total1, total2) == (F(-11, 512), F(-33, 512))
    print("Exact two-entry rational hook curvature sign control:", H)


def check_box(data, expected, cov):
    z = [{(a, m): F(c, entry["denominator"]) for a, m, c in entry["entries"]} for entry in data["tensors"]]
    assert len(z) == 35
    pol = list(map(F, expected["box_coefficients"]))
    assert len(pol) == 630
    operators = [polarization(t, cov) for t in z]
    pairs = tuple(combinations_with_replacement(range(len(z)), 2))
    exact = [evaluate(operators[i], z[j]) + (evaluate(operators[j], z[i]) if i != j else 0) for i, j in pairs]
    assert exact == pol
    r = F(data["certificate"]["box_radius"])
    error = sum(abs(c) * r**((i != 0) + (j != 0)) for c, (i, j) in zip(exact[1:], pairs[1:]))
    assert error == F(expected["variation"]) and exact[0] == F(expected["anchor"])
    assert exact[0] - error > F(146150093, 10**8)
    assert (exact[0] - error) / 3 > F(48716697, 10**8)
    print("All 630 curvature polynomial coefficients rebuilt directly from the ambient 91 curvature matrices.")
    print("CB box has curvature source > 1.46150093 and requires squared gradient > 0.48716697.")


def real_transport():
    # Independent exterior transport by minors of the complex covector map.
    _, frame = CA.split_frame()
    dual = CA.dual_frame(frame)
    table = dual.tolist()

    @lru_cache(None)
    def wedge(m):
        source = BC.bit_indices(m)
        choices = [[j for j in V if table[i][j]] for i in source]
        target = set(tuple(sorted(z)) for z in product(*choices) if len(set(z)) == 4)
        return [(sum(1 << j for j in inds), s.Matrix([[table[i][j] for j in inds] for i in source]).det())
                for inds in target]

    def convert(z):
        out = defaultdict(F)
        for a, m, c in z:
            for b, v in enumerate(table[a]):
                if v:
                    for n, w in wedge(m):
                        real = F(s.re(v * w))
                        if real:
                            out[b, n] += real * c
        return {key: c for key, c in out.items() if c}

    return convert


def inertia(M):
    """Inertia from exact rational LDL congruences."""
    M = s.Matrix(M)
    assert M == M.T
    pos = neg = 0
    n = M.rows
    while M.rows:
        piv = next((i for i in range(M.rows) if M[i, i]), None)
        if piv is None:
            ij = next(((i, j) for i, j in combinations(range(M.rows), 2) if M[i, j]), None)
            if ij is None:
                return [pos, neg, n - pos - neg]
            piv, j = ij
            T = s.eye(M.rows)
            T[j, piv] = 1
            M = T.T * M * T
        z = M[piv, piv]
        pos += int(bool(z > 0))
        neg += int(bool(z < 0))
        keep = [i for i in range(M.rows) if i != piv]
        col = M.extract(keep, [piv])
        M = M.extract(keep, keep) - col * col.T / z
    return [pos, neg, 0]


def check_types(lift, expected, cov):
    kinds = {tuple(x["weight"]): x for x in expected["blocks"]}
    grab = real_transport()
    signature = [0, 0, 0]
    negative_types = 0
    for x in lift["groups"]:
        key = tuple(x["type"])
        saved = kinds.pop(key)
        tensors = [grab(z) for z in x["vectors"]]
        d = len(tensors)
        assert d == saved["copies"] and x["per_copy"] == saved["per_copy"]
        mats = [polarization(z, cov) for z in tensors]
        gram_mat = [[F(0)] * d for _ in range(d)]
        form = [[F(0)] * d for _ in range(d)]
        for i, j in combinations_with_replacement(range(d), 2):
            u, v = tensors[i], tensors[j]
            gram_mat[i][j] = gram_mat[j][i] = sum((c * v.get((a, m), 0) / (G[a] * gram(m))
                                                 for (a, m), c in u.items()), F(0))
            form[i][j] = form[j][i] = (evaluate(mats[i], v) + evaluate(mats[j], u)) / 2
        # The original CA compact transform can choose the opposite global
        # sign of an extremal. Same real parts are fixed for this certificate.
        assert gram_mat == [list(map(F, z)) for z in saved["gram"]]
        assert form == [list(map(F, z)) for z in saved["form"]]
        assert inertia(gram_mat) == [d, 0, 0]
        by_sign = inertia(form)
        assert by_sign == [saved["inertia"][k] for k in "+-0"]
        negative_types += int(by_sign[1] > 0)
        if key in ((1, 0), (2, 0), (4, 0)):
            eigenvalues = {(1, 0): (F(-5, 8), F(17, 40)), (2, 0): (F(3, 4),), (4, 0): (F(-81, 140),)}[key]
            M = s.Matrix(gram_mat).inv() * s.Matrix(form)
            polynomial = M.charpoly()
            assert s.expand(polynomial.as_expr() - prod(polynomial.gen - c for c in eigenvalues)) == 0
        for j, count in enumerate(by_sign):
            signature[j] += count * x["per_copy"]
        print(f"{key}: {d} copies, inertia {by_sign}")
    assert not kinds and negative_types == 12 and signature == [8822, 1566, 0]
    print("Complete BX curvature-source signature (positive, negative, zero):", signature)


def run(full_ca):
    t = perf_counter()
    expected = json.loads((HERE / "CH_h4_stress_certificate.json").read_text(encoding="utf-8"))
    cb = read_packed(HERE / "CB_h4_partial_witness.json.gz", expected["CB_data_sha256"])
    lift = read_packed(HERE / "BX_data/h4_lift_witness.json.gz", expected["CA_data_sha256"])
    R = BG.curvature_model()
    cov = covariant_operations(R)
    explicit_calibration(R, cov)
    check_box(cb, expected, cov)
    if full_ca:
        CA.run(False)
    check_types(lift, expected, cov)
    print(f"Runtime {perf_counter() - t:.3f} seconds.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check-ca", action="store_true", help="also verify all CA split extremals and overlap witnesses")
    run(parser.parse_args().check_ca)
