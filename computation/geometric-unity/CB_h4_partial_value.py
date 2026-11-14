"""Verify CB's real unit H4 partial zero and its separating grade-five coefficient."""

import argparse
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
import gzip
import hashlib
from itertools import combinations, combinations_with_replacement, permutations
import json
from math import gcd, lcm, prod
from pathlib import Path
from time import perf_counter

import sympy as s

HERE = Path(__file__).resolve().parent

import BC_affine_hooks as BC
import BC_quadratic_xcheck as X
import BG_h6_curvature as BG
import BM_mixed_grade9 as BM
import BN_h6_invariant_grade9 as BN

G = tuple(BC.METRIC)
V = range(14)
FULL = (1 << 14) - 1
SELECTED_STRESS = ((0, 0), (0, 1), (0, 7), (0, 10))
BLADES = [prod(G[i] for i in BC.bit_indices(m)) for m in range(FULL + 1)]


def decode(entry):
    denominator = entry["denominator"]
    contents = entry["entries"]
    tensor = {(a, m): F(c, denominator) for a, m, c in contents}
    assert len(tensor) == len(contents) and denominator > 0
    assert all(a in V and m.bit_count() == 4 and not m & ~FULL and c for (a, m), c in tensor.items())
    return tensor


def perm_power(p, n):
    answer = tuple(range(5))
    for _ in range(n):
        answer = tuple(p[i] for i in answer)
    return answer


def tangent_character(p):
    f = sum(i == p[i] for i in range(5))
    f2 = sum(i == j for i, j in enumerate(perm_power(p, 2)))
    return (f * f + f2) // 2 - 1


def exterior_character(p, n):
    tr = [0] + [tangent_character(perm_power(p, k)) for k in range(1, n + 1)]
    e = [1]
    for k in range(1, n + 1):
        total = sum((-1)**(j - 1) * e[k - j] * tr[j] for j in range(1, k + 1))
        assert total % k == 0
        e.append(total // k)
    return e[-1]


def verify_detectors(scope):
    group = list(permutations(range(5)))
    assert sum(tangent_character(p)**2 + tangent_character(perm_power(p, 2)) for p in group) == 1200
    assert sum(exterior_character(p, 8) for p in group) == 120 * 28
    # Four trace-free invariant stress matrices and 28 independent projected blades.
    basis = [s.Matrix([[s.Rational(v) for v in row] for row in b]) for b in scope["symmetric_basis"]]
    gram = s.diag(*G)
    assert len(basis) == 4 and all(b == b.T and s.trace(gram.inv() * b) == 0 for b in basis)
    assert s.Matrix([[b[a, j] for b in basis] for a, j in SELECTED_STRESS]).det() == 1
    for j in range(4):
        action = s.zeros(14)
        for a, row in enumerate(BM.swap_rows(j, j + 1)):
            for b, value in row:
                action[b, a] = value
        assert action.det() == 1
        assert all(action * b * action.T == b for b in basis)

    def move(z, j):
        u = {}
        for m, c in z.items():
            for n, d in BM.swap_form(j, j + 1, m):
                BC.plus(u, n, c * d)
        return u

    def average(seed):
        z = {seed: F(1)}
        for j in range(4):
            acc = dict(z)
            for k in range(j, -1, -1):
                z = move(z, k)
                for m, c in z.items():
                    BC.plus(acc, m, c)
            z = acc
        assert all(move(z, j) == z for j in range(4))
        return z

    selectors = scope["wedge_detectors"]
    seeds = scope["wedge_seeds"]
    assert len(set(selectors)) == len(set(seeds)) == 28 and all(m.bit_count() == 8 for m in (*selectors, *seeds))
    columns = list(map(average, seeds))
    matrix = s.Matrix([[s.Rational(z.get(m, 0)) for z in columns] for m in selectors])
    det = matrix.det()
    assert det != 0
    print("S5 receiver dimensions: 4 trace-free symmetric, 28 invariant eight-forms. Both detector maps injective.")
    print("Wedge detector determinant (unnormalized group sums):", det)
    return selectors


def invariant_multiplicities():
    weights = [(x, y) for x, y in ((2, 0), (-2, 0), (0, 2), (0, -2), (1, 1), (1, -1),
                                  (-1, 1), (-1, -1), (1, 0), (-1, 0), (0, 1), (0, -1), (0, 0), (0, 0))]
    denominator = Counter({(0, 0): 1})
    for a, b in ((1, -1), (1, 1), (1, 0), (0, 1)):
        for (x, y), c in list(denominator.items()):
            denominator[x - a, y - b] -= c

    def constant(char):
        return sum(c * char[-x, -y] for (x, y), c in denominator.items())

    def wedge(n):
        return Counter((sum(x for x, _ in ws), sum(y for _, y in ws)) for ws in combinations(weights, n))

    assert constant(wedge(8)) == 0
    for n in (1, 5, 9):
        character = Counter()
        for (x, y), c in wedge(n).items():
            for a, b in weights:
                character[x + a, y + b] += c
        assert constant(character) == 1
    print("SO(5) invariant multiplicities: one each in grades 1, 5, 9; zero in eight-forms.")


def reconstruct_l5():
    # Antisymmetrized trace of six symmetric matrices, with the external
    # matrix kept outside the five ordered wedge variables.
    matrices = BG.integer_phi()
    sparse = [[(i, j, c) for i, row in enumerate(A) for j, c in enumerate(row) if c] for A in matrices]
    powers = {0: [int(i == j) for i in range(5) for j in range(5)]}
    for n in range(1, 6):
        accum = {}
        for inds in combinations(V, n):
            m = sum(1 << j for j in inds)
            z = [0] * 25
            for i, j in enumerate(inds):
                prev = powers.get(m ^ (1 << j))
                if prev is None:
                    continue
                sign = (-1)**(n - 1 - i)
                for l, r, c in sparse[j]:
                    for k in range(5):
                        z[5*k + r] += sign * c * prev[5*k + l]
            if any(z):
                accum[m] = z
        powers = accum
    all_entries = {(a, m): c for m, p in powers.items() for a, A in enumerate(sparse)
                   if (c := sum(v * p[5*j + i] for i, j, v in A))}
    g = gcd(*all_entries.values())
    assert g == 50000 and len(all_entries) == 1143
    receiver = {key: v // g for key, v in all_entries.items()}
    assert all(not BC.action(j, receiver) for j in range(10))
    return receiver


def fields(tensor):
    by_index = [dict() for _ in V]
    inside = [[dict() for _ in V] for _ in V]
    for (a, m), c in tensor.items():
        by_index[a][m] = c
        for i in BC.bit_indices(m):
            k = m ^ (1 << i)
            inside[a][i][k] = BC.sign(1 << i, k) * c
    return by_index, inside


def dot(left, right):
    if len(right) < len(left):
        left, right = right, left
    return sum((v * right[m] / BLADES[m] for m, v in left.items() if m in right), F(0))


def lower(a, b, selectors):
    ax, ai = a
    bx, bi = b
    n = sum(dot(ax[i], bx[i]) / G[i] for i in V)
    stress = []
    for i, j in SELECTED_STRESS:
        E = (dot(ax[i], bx[j]) + dot(bx[i], ax[j])) / 2
        I = sum((dot(ai[k][i], bi[k][j]) + dot(bi[k][i], ai[k][j])) / (2 * G[k]) for k in V)
        stress.append(2 * E - I + (n * G[i] / 7 if i == j else 0))
    wedge = []
    for m in selectors:
        w = F(0)
        for k in V:
            z = sum(BC.sign(h, m ^ h) * c * bx[k][m ^ h]
                    for h, c in ax[k].items() if h & m == h and m ^ h in bx[k])
            w += z / G[k]
        wedge.append(w)
    return [*stress, *wedge, n]


def check_formula_patterns():
    algebra = X.BladeAlgebra([F(1)] * 14)
    b = [{30: F(1)}, {29: F(1)}] + [{} for _ in range(12)]
    observed = [algebra.quadratic_coefficient(b, i, 1 << i) for i in V]
    assert observed == [F(0)] * 2 + [F(-4)] * 3 + [F(-4, 3)] * 9
    assert sum(observed) == -24
    # Add a disjoint four-form at index zero with its alternating-trace mate.
    b[0][480] = F(1)
    b[5][449] = F(1)
    alt, ctr = {}, {}
    square = {}
    for a, field in enumerate(b):
        for m, c in field.items():
            k, sgn = algebra.multiply(1 << a, m)
            BC.plus(ctr if m >> a & 1 else alt, k, sgn * c)
            for v, d in field.items():
                if not v & m:
                    BC.plus(square, v | m, BC.sign(m, v) * c * d)
    assert not alt and not ctr and square == {510: F(2)}
    for a in (0, 1, 4, 5, 8, 9, 13):
        mask = 510 | (1 << a)
        if mask.bit_count() == 9:
            value = F(-2, 3) * BC.sign(1 << a, 510) * square[510]
            assert value == algebra.quadratic_coefficient(b, a, mask)
    # The unique grade-nine SO(5) receiver pulls back to zero through wedge.
    pullback = {}
    for (a, m), c in BN.invariant_receiver().items():
        if m >> a & 1:
            k = m ^ (1 << a)
            BC.plus(pullback, k, BC.sign(1 << a, k) * c / G[a])
    assert not pullback
    print("Independent exact Clifford controls pass for the grade-one normalization and nonzero grade-nine samples.")


def prefix(m):
    for k in (1, 2, 4, 8):
        m ^= m << k
    return (m << 1) & FULL


def clifford_polynomials(tensors, targets):
    """Exact coordinate polynomials by direct six-word trace expansion."""
    d = len(tensors)
    pairs = tuple(combinations_with_replacement(range(d), 2))
    indexes = [[0 for _ in range(d)] for _ in range(d)]
    for k, (i, j) in enumerate(pairs):
        indexes[i][j] = indexes[j][i] = k
    denominator = [lcm(*(w.denominator for w in z.values())) for z in tensors]
    integer = [dict() for _ in V]
    for j, z in enumerate(tensors):
        for (a, m), c in z.items():
            integer[a].setdefault(m, []).append((j, int(c * denominator[j])))
    repeat_den = lcm(*(v.numerator for v in BLADES))
    repeat_weights = [int(F(repeat_den, v)) for v in BLADES]
    assert all(F(v) * g == repeat_den for v, g in zip(repeat_weights, BLADES))
    signs = [prefix(m) for m in range(FULL + 1)]
    parity = [1 - (-1)**xy + (-1)**tx + (-1)**(tx + ty) - (-1)**(ty + xy) - (-1)**(tx + ty + xy)
              for tx in range(2) for ty in range(2) for xy in range(2)]
    routes = [[(b, c, (1 << a) | (1 << b) | (1 << c), -1 if b < a < c else 1)
               for b, c in combinations(V, 2) if a not in (b, c)] for a in V]
    numerators = [0] * len(pairs)
    weights = {key: F(value) * G[key[0]] for key, value in targets.items()}
    common = lcm(*(c.denominator for c in weights.values()))
    for (a, m), factor in weights.items():
        fact = int(factor * common)
        for b, c, z, ordering in routes[a]:
            for x, vx in integer[b].items():
                y = x ^ z ^ m
                if y not in integer[c]:
                    continue
                tx, ty, xy = (z & x).bit_count() & 1, (z & y).bit_count() & 1, (x & y).bit_count() & 1
                e = parity[4 * tx + 2 * ty + xy]
                if not e:
                    continue
                e *= ordering * fact * repeat_weights[(z & x) | (z & y) | (x & y)]
                if ((z & signs[x]) ^ ((z ^ x) & signs[y])).bit_count() & 1:
                    e = -e
                for i, u in vx:
                    for j, v in integer[c][y]:
                        numerators[indexes[i][j]] += e * u * v
    scale = -F(2, 3 * common * repeat_den)
    return {key: scale * F(v, denominator[key[0]] * denominator[key[1]]) for key, v in zip(pairs, numerators)}


def check_geometry(tensors):
    R = BG.curvature_model()
    for z in tensors:
        assert BC.hook_traces(z) == ({}, {})
        assert not BG.full_sparse_curvature(R, BC.hodge_scaled(z))
        assert all(BM.swap(z, j, j + 1) == {key: -v for key, v in z.items()} for j in range(4))
    print(f"All {len(tensors)} rational tensors pass both hooks, full dual curvature, and all four signed permutation generators.")


def certify(poly, cert):
    n = len(poly)
    assert all(len(q) == n + 1 and all(len(row) == n + 1 for row in q) for q in poly)
    assert all(q[i][j] == q[j][i] for q in poly for i in range(n + 1) for j in range(i))
    inverse = [list(map(F, row)) for row in cert["inverse"]]
    assert s.Matrix(inverse).det()
    radius = F(cert["box_radius"])
    raw, defect, local_lipschitz = [], [], []
    for i, p in enumerate(inverse):
        def entry(j, k):
            return sum(v * q[j][k] for v, q in zip(p, poly))

        raw.append(abs(entry(0, 0)))
        defect.append(sum(abs((i == j) - 2 * entry(0, j + 1)) for j in range(n)))
        local_lipschitz.append(2 * sum(abs(entry(j, k)) for j in range(1, n + 1) for k in range(1, n + 1)))
    lip = max(defect) + radius * max(local_lipschitz)
    image = max(raw) + radius * lip
    assert lip == F(cert["contraction_bound"]) < 1
    assert image == F(cert["image_radius_bound"]) < radius
    print(f"Rational sup-norm box radius {radius}. Contraction <= {float(lip):.10g}; image radius <= {float(image):.10g}.")
    return radius


def run(full):
    start = perf_counter()
    packed = gzip.decompress((HERE / "CB_h4_partial_witness.json.gz").read_bytes())
    manifest = json.loads((HERE / "CB_h4_partial_manifest.json").read_text(encoding="utf-8"))
    assert hashlib.sha256(packed).hexdigest() == manifest["raw_sha256"]
    data = json.loads(packed)
    z = list(map(decode, data["tensors"]))
    assert len(z) == 35
    n = len(z) - 1
    polys = [[list(map(F, row)) for row in M] for M in data["polynomials"]]
    assert len(polys) == n
    detectors = verify_detectors(data["scope"])
    invariant_multiplicities()
    l5 = reconstruct_l5()
    check_formula_patterns()
    check_geometry(z)
    by_index = list(map(fields, z))
    for i, j in combinations_with_replacement(range(n + 1), 2):
        expected = lower(by_index[i], by_index[j], detectors)
        expected[-1] -= int(i == j == 0)
        stored = [p[i][j] for p in polys[:-1]]
        assert expected == stored, (i, j, [(k, u, v) for k, (u, v) in enumerate(zip(expected, stored)) if u != v])
    scale = {key: F(v, G[key[0]] * BLADES[key[1]]) for key, v in l5.items()}
    pairs = tuple(combinations_with_replacement(range(n + 1), 2))
    selected = tuple(range(n + 1)) if full else (0, 1, 8, 34)
    computed = clifford_polynomials([z[i] for i in selected], scale)
    for (i, j), actual in computed.items():
        assert (1 if i == j else 2) * polys[-1][selected[i]][selected[j]] == actual
    print("All lower polynomial entries and " + ("all 630" if full else "10 selected") + " L5 entries independently reconstructed.")
    radius = certify(polys, data["certificate"])
    w = data["separating_receiver"]
    a, m = w["target"]
    assert m.bit_count() == 5
    separate = clifford_polynomials(z, {(a, m): 1})
    assert [separate[k] for k in pairs] == list(map(F, w["coefficients"]))
    assert separate[0, 0] == X.BladeAlgebra(G).quadratic_coefficient(by_index[0][0], a, m)
    variation = sum(abs(separate[i, j]) * radius**((i != 0) + (j != 0)) for i, j in pairs[1:])
    assert variation == F(w["exact_variation_bound"]) and separate[0, 0] + variation < 0
    print(f"Separating grade-five coordinate {(a, m)} < {float(separate[0, 0] + variation):.12g} everywhere on the box.")
    print(f"Runtime {perf_counter() - start:.3f} seconds.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--full", action="store_true", help="reconstruct all L5 coefficients")
    run(parser.parse_args().full)
