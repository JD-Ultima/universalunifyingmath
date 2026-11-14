"""Verify the DS and DP positive receivers on exact H6 axial value spaces."""

import argparse
from fractions import Fraction as F
import gzip
from itertools import product
import json
from math import isqrt, prod
from pathlib import Path
from time import perf_counter

import numpy as np
import sympy as s

HERE = Path(__file__).resolve().parent
import CV_h6_axial as CV

# Curvature dimension, image dimension, joint rank, rank with D, chain counts.
BLOCKS = {
    0: (389, 64, 628, None, (2, 23, 13, 5)),
    3: (423, 83, 653, None, (3, 27, 19, 6)),
    5: (423, 83, 653, None, (3, 27, 19, 6)),
    6: (423, 83, 653, None, (3, 27, 19, 6)),
    9: (423, 83, 653, None, (3, 27, 19, 6)),
    15: (443, 105, 628, 372, (4, 35, 26, 6)),
    17: (423, 83, 653, None, (3, 27, 19, 6)),
    27: (443, 110, 634, None, (4, 31, 29, 7)),
}
THREE_SIGNS = (0, 15, 27)
LOCAL_DIMS = {0: 107, 3: 138, 5: 138, 6: 138, 9: 138, 15: 187, 17: 138, 27: 181}
BITS = 46
RECEIVER_DEN = 1 << 20


def archive(name):
    return json.loads(gzip.decompress((HERE / name).read_bytes()))


def build_pencil(ch):
    coords, _ = CV.BR.stack(ch)
    static, axial = {}, {}
    for j, (a, m) in enumerate(coords):
        if w := CV.BR.wedge(a, m):
            static.setdefault((0, 0, w[0]), {})[j] = CV.BR.den * w[1]
        else:
            reduced = m ^ (1 << a)
            f = F(CV.BR.den, CV.G[a])
            assert f.denominator == 1
            static.setdefault((1, 0, reduced), {})[j] = f.numerator * CV.BC.sign(1 << a, reduced)
        for (b, k), c in CV.BR.closure_column(a, m).items():
            static.setdefault((2, b, k), {})[j] = c
        if w := CV.BR.wedge(13, m):
            axial.setdefault((2, a, w[0]), {})[j] = CV.BR.den * w[1]
    labels = sorted(static.keys() | axial.keys())
    return coords, [static.get(k, {}) for k in labels], [axial.get(k, {}) for k in labels]


def row_dot(row, col):
    return sum(c * col.get(j, 0) for j, c in row.items())


def check_chain(columns, A, D):
    degree = len(columns) - 1
    # Chain archives are primitive integer coefficients; the exact matrix
    # product avoids a Fraction allocation on each absent entry.
    weights = [p for v in columns for p in v.values()]
    if all(p.denominator == 1 for p in weights):
        columns = [{i: int(v) for i, v in c.items()} for c in columns]
    for a, d in zip(A, D):
        assert row_dot(a, columns[0]) == 0
        for j in range(1, degree + 1):
            assert row_dot(a, columns[j]) + row_dot(d, columns[j-1]) == 0
        assert row_dot(d, columns[degree]) == 0


def read_point_space(ch, coords, A, D):
    filename = "CZ_h6_axial_local_raw_chains.json.gz" if ch == 0 else f"DS_chains_{ch:02d}.json.gz"
    raw = archive(filename)
    assert raw["ch"] == ch
    loc = {x: i for i, x in enumerate(coords)}
    R, image_source = [], []
    count = BLOCKS[ch][4]
    assert [d for d, _ in raw["chains"]] == [d for d, n in enumerate(count) for _ in range(n)]
    for degree, entries in raw["chains"]:
        columns = [{} for _ in range(degree + 1)]
        for k, a, m, c in entries:
            assert 0 <= k <= degree and loc[a, m] not in columns[k]
            columns[k][loc[a, m]] = F(c)
        check_chain(columns, A, D)
        image_source.extend(range(len(R), len(R) + degree))
        R.extend(columns)
    E = []
    if ch == 15:
        regular = archive("DS_regular_15.json.gz")
        assert regular["character"] == 15
        assert [(F(z["rate"]), len(z["witnesses"])) for z in regular["regular"]] == [
            (0, 2), (-4, 7), (F(-32, 5), 2)]
        for section in regular["regular"]:
            z = F(section["rate"])
            for entries in section["witnesses"]:
                b = {loc[a, m]: F(c) for a, m, c in entries}
                assert len(b) == len(entries)
                assert all(row_dot(a, b) + z * row_dot(d, b) == 0 for a, d in zip(A, D))
                E.append(b)
    assert len(R) + len(E) == LOCAL_DIMS[ch]
    return R, image_source, E


def exact_dimensions():
    path = Path(CV.BO.__file__).with_name("BO_h6_full_kernel_witnesses.json")
    highest = json.loads(path.read_text(encoding="utf-8"))["extremal_vectors"]

    def partition(a, b):
        return sum(b + x - y >= 0 for x in range(a+1) for y in range(a-x+1)) if a >= 0 else 0

    def weight(l, m, a, b):
        result = 0
        for swap in (False, True):
            L, M = (2*m+1, 2*l+3) if swap else (2*l+3, 2*m+1)
            for i, j in product((-1, 1), repeat=2):
                result += (1 - 2*swap)*i*j*partition((i*L - (2*a+3))//2, (j*M - (2*b+1))//2)
        return result

    total = 0
    output = dict.fromkeys(BLOCKS, 0)
    for item in highest:
        l, m = item["weight"]
        table = {(a, b): weight(l, m, a, b) for a, b in product(range(-l, l+1), repeat=2)}
        dim = sum(table.values())
        assert min(table.values()) >= 0 and dim == (l-m+1)*(l+m+2)*(2*l+3)*(2*m+1)//6
        traces = {
            0: dim,
            2: sum((-1)**(a % 2)*v for (a, b), v in table.items()),
            4: sum((-1)**((a+b) % 2)*v for (a, b), v in table.items()),
        }
        copies = len(item["vectors"])
        total += dim * copies
        for ch in output:
            num = sum((-1)**((ch & g).bit_count() % 2)*traces[g.bit_count()]
                      for g in range(32) if g.bit_count() % 2 == 0)
            assert num % 16 == 0
            output[ch] += copies * (num // 16)
    assert total == 6834 and output == {ch: x[0] for ch, x in BLOCKS.items()}
    return output


def quotient_proof(ch, coords, A, D, R, source, E, dim):
    p = CV.P
    curvature = CV.BG.curvature_model()
    entries = [[(i, j, b, CV.finite(-v)) for (i, j), matrix in curvature.items()
                for b, v in enumerate(matrix[a]) if v] for a in CV.V]
    rows = {}
    for col, (a, mask) in enumerate(coords):
        if w := CV.BR.wedge(a, mask):
            rows.setdefault((0, 0, w[0]), {})[col] = w[1] % p
        else:
            w = mask ^ (1 << a)
            rows.setdefault((1, 0, w), {})[col] = CV.finite(F(CV.BC.sign(1 << a, w), CV.G[a]))
        for i, j, b, v in entries[a]:
            w = (1 << i) | (1 << j)
            if not w & mask and v:
                rows.setdefault((2, b, w | mask), {})[col] = v * CV.BC.sign(w, mask) % p
    K = CV.modular_sparse_kernel(list(rows.values()), len(coords))
    assert K.shape == (len(coords), dim)

    def dense(source_rows):
        X = np.zeros((len(source_rows), len(coords)), dtype=np.int64)
        for i, row in enumerate(source_rows):
            for j, v in row.items():
                X[i, j] = CV.finite(v)
        return X

    Ddense = dense(D)
    X = dense(R + E).T
    Aimage, Dimage = dense(A) @ K % p, Ddense @ K % p
    T = Ddense @ X[:, source] % p
    _, t_dim, joint, just_d, _ = BLOCKS[ch]
    rank = CV.BO.modular_rank
    assert rank(X.copy()) == len(R) + len(E)
    assert rank(T.copy()) == len(source) == t_dim
    assert rank(np.column_stack((T, Aimage, Dimage))) == joint
    d = dim - len(R)
    if just_d is None:
        assert joint - t_dim == 2*d and not E
    else:
        assert rank(np.column_stack((T, Dimage))) == just_d == t_dim + d
        assert 2*d - (joint - t_dim) == len(E) == 11
    print(f"Character {ch}: exact local dimension {len(R)+len(E)}, image {t_dim}, "
          f"quotient joint rank lower bound {joint-t_dim}.", flush=True)


def unit_proxy(ch, coords, raw, proof, prefix):
    change = archive(f"{prefix}_coordinates_{ch:02d}.json.gz")
    n = len(raw)
    assert change["character"] == ch and change["raw_basis_size"] == n
    rows = change["coordinate_matrix"]
    assert len(rows) == n and all(len(x) == n for x in rows)
    grid, scale = 1 << 95, 1 << BITS
    metric = [CV.physical(*x) for x in coords]
    roots = [F(isqrt(g.denominator * grid * grid // g.numerator), grid) for g in metric]
    assert all(h*h <= 1/g < (h + F(1, grid))**2 for h, g in zip(roots, metric))
    W = np.zeros((len(coords), n), dtype=np.int64)
    err = F(0)
    for j, row in enumerate(rows):
        x = {}
        for col, w in zip(raw, row):
            if v := F(w):
                for i, c in col.items():
                    x[i] = x.get(i, 0) + v*c
        for i, c in x.items():
            h = c * roots[i]
            g = round(scale*h)
            assert abs(g) < 1 << 60
            W[i, j] = g
            err += (abs(F(g, scale) - h) + abs(c)/grid)**2
    assert err == F(proof["unit_proxy_sum_error2_exact"]) < F(1, 10**20)
    assert sum(int(x)**2 for x in W.ravel()) < (n+1)*scale**2
    return W


def physical_certificate(ch, coords, raw, Y, full, prefix):
    proof = archive(f"{prefix}_positive_{ch:02d}.json.gz")
    W = unit_proxy(ch, coords, raw, proof, prefix)
    n = W.shape[1]
    assert proof["character"] == ch and proof["C"] == BITS and proof["receiver_gram"] == RECEIVER_DEN
    H = np.asarray(proof["polarization_integer"], dtype=object)
    assert H.shape == (n, n) and np.array_equal(H.T, H)
    factor = np.asarray(proof["factor"], dtype=object)
    assert factor.shape == (n, n)
    unit = proof["denominator"]
    scale = proof["factor_denominator"]
    assert unit == 21*2**(2*BITS+20) and scale == 1 << 36
    trace = sum(val for a, mask, val in Y if mask == 1 << a)
    C = 26*trace*(W.astype(object).T @ W.astype(object)) - H
    assert all(v % 7 == 0 for v in C.ravel())
    C //= 7
    y = sum(abs(v) for _, _, v in Y)
    l1 = max(sum(abs(int(v)) for v in col) for col in W.T)
    cap = 120 * l1*l1*y
    primes = proof["primes"]
    assert len(set(primes)) == len(primes) and all(s.isprime(p) and p % 4 == 1 for p in primes)
    assert cap == int(proof["integer_bound"]) and prod(primes) > 2*cap
    assert max(abs(v) for v in C.ravel()) <= cap
    for p in primes if full else primes[:1]:
        M = CV.BV.Gradient(ch, coords, Y, p)
        right = np.column_stack([M.core(M.matrix(v)) for v in W.T])
        assert np.array_equal(C % p, (W % p).T @ right % p)
        print(f"Character {ch}: exact {n} by {n} Gram modulo {p}.", flush=True)
    shift = F(proof["shift"])
    Z = factor @ factor.T
    defect = max(sum(abs(F(int(H[i, j]), unit) - F(int(Z[i, j]), scale*scale)
                         - (shift if i == j else 0)) for j in range(n)) for i in range(n))
    assert defect == F(proof["exact_defect"])
    tolerance = F(1, 10**10)
    loss = F(16*y, RECEIVER_DEN) * (2*(n+1)*tolerance + tolerance*tolerance)
    assert loss == F(proof["exact_loss"])
    assert F(n+2) > n+1 + 2*(n+1)*tolerance + tolerance*tolerance
    margin = (shift - defect - loss) / (n+2)
    assert margin == F(proof["uniform_receiver_margin"]) > F(1, 1250)
    print(f"Character {ch}: rigorous physical receiver margin {float(margin):.9g}.", flush=True)


def main(full, selected, geometry_only, prefix="DS"):
    t = perf_counter()
    common = archive(f"{prefix}_receiver.json.gz")
    assert common["denominator"] == RECEIVER_DEN
    Y = common["rows"]
    valid = set(CV.BV.RECEIVERS)
    assert len({(a, m) for a, m, v in Y}) == len(Y)
    assert all((a, m) in valid and isinstance(v, int) for a, m, v in Y)
    dim = exact_dimensions()
    for ch in selected:
        coords, A, D = build_pencil(ch)
        R, source, E = read_point_space(ch, coords, A, D)
        quotient_proof(ch, coords, A, D, R, source, E, dim[ch])
        if not geometry_only:
            physical_certificate(ch, coords, R + E, Y, full, prefix)
    if not geometry_only:
        print(f"One receiver is > 1/1250 times squared norm on the direct sum of {sum(LOCAL_DIMS[ch] for ch in selected)} dimensions.", flush=True)
    print(f"Runtime {perf_counter()-t:.3f} seconds.", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--full", action="store_true", help="verify all CRT primes")
    parser.add_argument("--character", type=int, choices=BLOCKS, action="append", help="verify specified sign masks")
    parser.add_argument("--geometry-only", action="store_true", help="skip the quadratic receiver blocks")
    parser.add_argument("--eight-signs", action="store_true", help="verify the full DP eight-character domain")
    args = parser.parse_args()
    available = tuple(BLOCKS) if args.eight_signs else THREE_SIGNS
    selected = args.character or available
    if not set(selected) <= set(available):
        parser.error("new character certificates require --eight-signs")
    main(args.full, selected, args.geometry_only, "DP" if args.eight_signs else "DS")
