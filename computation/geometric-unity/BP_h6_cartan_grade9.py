"""Check four nonzero Cartan-square grade-nine responses in the physical frame."""

import json
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path
from time import perf_counter

import sympy as s

import BN_h6_invariant_grade9 as BN
import BO_h6_full_kernel as BO


CASES = {
    (3, 2): (10, 14709, F(-399015), F(0)),
    (3, 0): (4, 10367, F(0), -F(12623, 2)),
    (2, 1): (13, 10997, -F(7449791, 30), F(0)),
    (2, 0): (9, 8927, F(0), -F(431, 5)),
}


def split_to_physical_covectors():
    axes = s.Matrix([[1, 1, 0, 0, 0], [s.I, -s.I, 0, 0, 0],
                     [0, 0, 1, 1, 0], [0, 0, s.I, -s.I, 0], [0, 0, 0, 0, 1]])
    compact = [s.Matrix(q) / 5 for q in BN.BG.integer_phi()]
    split = []
    for i in range(5):
        for j in range(i, 5):
            if (i, j) in ((0, 1), (2, 3), (4, 4)):
                continue
            c = BO.matrix_units(i, j)
            split.append(axes * (c if i == j else c + c.T) * axes.T)
    for i, j in ((0, 1), (2, 3)):
        split.append(axes * (BO.matrix_units(i, j) + BO.matrix_units(j, i)
                             - 4 * BO.matrix_units(4, 4)) * axes.T)
    change = s.Matrix(14, 14, lambda a, b: s.trace(compact[a] * split[b]) / BN.GRAM[a])
    inverse = change.inv()
    assert change * inverse == s.eye(14)
    # The inverse vector matrix transforms all covariant components.
    return tuple(tuple((j, F(s.re(v)), F(s.im(v))) for j, v in enumerate(row) if v)
                 for row in inverse.tolist())


def complex_multiply(a, b, c, d):
    return a * c - b * d, a * d + b * c


def transport(entries, transform):
    forms = {}
    for _, mask, _ in entries:
        if mask in forms:
            continue
        wedge = {0: (F(1), F(0))}
        for axis in BO.axes(mask):
            next_wedge = {}
            for partial, (a, b) in wedge.items():
                for j, c, d in transform[axis]:
                    if partial >> j & 1:
                        continue
                    sign = (-1) ** (partial >> (j + 1)).bit_count()
                    target = partial | (1 << j)
                    v, w = complex_multiply(a, b, c, d)
                    re, im = next_wedge.get(target, (0, 0))
                    next_wedge[target] = re + sign * v, im + sign * w
            wedge = next_wedge
        forms[mask] = wedge

    real, imag = {}, {}
    for a, mask, value in entries:
        for b, re, im in transform[a]:
            for target, (v, w) in forms[mask].items():
                c, d = complex_multiply(re, im, v, w)
                BO.add(real, (b, target), value * c)
                BO.add(imag, (b, target), value * d)
    return real, imag


def coefficient(algebra, tensor, a, target):
    field = [{} for _ in range(14)]
    for (b, mask), value in tensor.items():
        blade, factor = algebra.multiply(mask, BN.FULL)
        field[b][blade] = value * factor
    return BN.DETERMINANT * algebra.quadratic_coefficient(field, a, target)


def main():
    start = perf_counter()
    input_path = Path(__file__).with_name('BO_h6_full_kernel_witnesses.json')
    data = json.loads(input_path.read_text(encoding='utf-8'))['extremal_vectors']
    vectors = {tuple(row['weight']): row['vectors'][0] for row in data}
    transform = split_to_physical_covectors()
    curvature = BN.BG.curvature_model()
    algebra = BN.X.BladeAlgebra(BN.GRAM)
    for w, (a, mask, expectation_re, expectation_im) in CASES.items():
        tensor_re, tensor_im = transport(vectors[w], transform)
        assert all(BN.BC.hook_traces(tensor) == ({}, {})
                   and not BN.BG.full_sparse_curvature(curvature, tensor)
                   for tensor in (tensor_re, tensor_im))
        combined = BN.BF.combine((1, tensor_re), (1, tensor_im))
        assert mask.bit_count() == 9
        real_square = coefficient(algebra, tensor_re, a, mask)
        imag_square = coefficient(algebra, tensor_im, a, mask)
        result = (real_square - imag_square,
                  coefficient(algebra, combined, a, mask) - real_square - imag_square)
        assert result == (expectation_re, expectation_im) and any(result)
        algebra.multiply.cache_clear()
        print('Weight', w, 'receiver', (a, mask), 'real and imaginary parts', *result, flush=True)

    for i, w in enumerate(CASES):
        remaining = list(CASES)[i:]
        for j, u in enumerate(remaining):
            for v in remaining[j:]:
                if u == v == w:
                    continue
                # (2*w) cannot occur in any remaining tensor product except w*w.
                assert u[0] + v[0] < 2 * w[0] or sum(u) + sum(v) < 2 * sum(w)
    print('All four descending highest-weight exclusions verified.')
    print(f'Runtime {perf_counter() - start:.3f} seconds.')


if __name__ == '__main__':
    main()
