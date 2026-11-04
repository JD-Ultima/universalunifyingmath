"""Exact nonlinear obstruction for the affine-invariant pure H4 branch.

Reconstructs BC's invariant/closure bases without running its main. Eighteen
genuine higher-grade coefficients and exact rank-one elimination suffice.
No private imports, saved basis receipts, numerical fits or output writes.
"""
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations
from pathlib import Path
import sys
from time import perf_counter
import sympy as sp

sys.path.insert(0, str(Path(__file__).resolve().parent))
import BC_affine_hooks as BC

N = 14
MONOMIALS = tuple((i, j) for i in range(6) for j in range(i, 6))
OUTPUTS = (
    (5, 0, (0, 1, 2, 4, 5)),
    (5, 0, (0, 1, 3, 4, 7)),
    (5, 0, (2, 4, 5, 7, 8)),
    (5, 0, (0, 5, 6, 7, 8)),
    (5, 0, (0, 4, 5, 6, 10)),
    (5, 0, (0, 2, 3, 4, 13)),
    (5, 0, (0, 4, 5, 6, 13)),
    (5, 0, (3, 5, 7, 8, 13)),
    (5, 4, (2, 3, 4, 5, 6)),
    (5, 4, (0, 1, 4, 5, 7)),
    (5, 4, (0, 3, 5, 6, 7)),
    (5, 4, (4, 5, 6, 7, 8)),
    (5, 4, (0, 2, 4, 5, 13)),
    (5, 4, (2, 3, 6, 7, 13)),
    (5, 4, (6, 7, 9, 10, 13)),
    (5, 13, (1, 2, 4, 5, 13)),
    (9, 0, (0, 1, 2, 4, 6, 9, 10, 11, 12)),
    (9, 4, (0, 1, 2, 3, 4, 5, 8, 11, 12)),
)
ALG = BC.Clifford(BC.METRIC)


def as_field(tensor):
    return [{mask: value for (b, mask), value in tensor.items() if a == b}
            for a in range(N)]


def reconstruct():
    invariant = BC.invariant_basis(4)
    hooks = BC.hook_basis(invariant)
    assert (len(invariant), len(hooks)) == (20, 13)
    dual = [BC.hodge_scaled(tensor) for tensor in hooks]
    assert all(BC.hook_traces(v) == ({}, {}) for v in dual)
    images = [BC.exterior_covariant(v) for v in dual]
    rank, _, _, determinant, kernel = BC.rank_witness(images)
    assert rank == 7 and len(kernel) == 6 and determinant != 0
    raw = [BC.combine(vector, hooks) for vector in kernel]
    difference = dict(raw[5])
    for key, value in raw[4].items():
        BC.plus(difference, key, -value)
    raw[5] = difference
    basis = [BC.primitive(tensor) for tensor in raw]
    assert all(BC.hook_traces(v) == ({}, {}) for v in basis)
    assert all(not BC.action(j, v) for v in basis for j in range(6))
    assert all(not BC.exterior_covariant(BC.hodge_scaled(v)) for v in basis)
    print('Reconstructed20 SO4 invariants,13 actual H4 hooks, closure rank7/kernel6.')
    print('Closure minor determinant:', determinant)
    print('Raw closure combinations:', [tuple((j, str(c)) for (j, _), c in sorted(v.items())) for v in kernel])
    print('Final direction replaced by K5-K4; each field then normalized primitively.')
    print('Primitive basis support sizes:', tuple(len(v) for v in basis))
    return hooks, [as_field(tensor) for tensor in basis]


@lru_cache(None)
def triple_factor(a, b, c):
    return ALG.factor(a, b)*ALG.factor(a ^ b, c)


@lru_cache(None)
def ordered_factor(a, b, c, left, right):
    triple = (1 << a) | (1 << b) | (1 << c)
    sign = ALG.factor(1 << a, 1 << b)*ALG.factor((1 << a) ^ (1 << b), 1 << c)
    value = (triple_factor(triple, left, right)-triple_factor(triple, right, left)
             +triple_factor(left, right, triple)-triple_factor(right, left, triple)
             +triple_factor(left, triple, right)-triple_factor(right, triple, left))
    return -F(2, 3)*BC.METRIC[a]*sign*value


def coefficient_matrix(a, target, fields):
    """Full covariant coefficient of Q(iC)/i, not a trace or projection ansatz."""
    indexed = []
    for b in range(N):
        masks = set().union(*(field[b] for field in fields))
        indexed.append({mask: [(j, field[b][mask]) for j, field in enumerate(fields) if mask in field[b]]
                        for mask in masks})
    size = len(fields)
    result = [[F(0)]*size for _ in range(size)]
    for b, c in combinations(range(N), 2):
        if a in (b, c):
            continue
        flip = (1 << a) ^ (1 << b) ^ (1 << c) ^ target
        for left, values in indexed[b].items():
            right = left ^ flip
            if right not in indexed[c]:
                continue
            factor = ordered_factor(a, b, c, left, right)
            if factor:
                for j, value in values:
                    for k, extra in indexed[c][right]:
                        result[j][k] += factor*value*extra
    return sp.Matrix(size, size, lambda i, j:
                     sp.Rational((result[i][j]+result[j][i]).numerator,
                                 2*(result[i][j]+result[j][i]).denominator))


def validate(hooks):
    controls = []
    for terms in (
        ((1, (2, 3, 6, 7), F(3, 7)), (3, (4, 5, 6, 7), F(-5, 11))),
        ((1, (2, 3, 4, 6), F(2, 5)), (3, (9, 10, 11, 12), F(-4, 9)))):
        field = [{} for _ in range(N)]
        for a, indices, value in terms:
            field[a][sum(1 << i for i in indices)] = value
        controls.append(field)
    controls.append(as_field({key: F(7, 13)*value for key, value in hooks[0].items()}))
    for j, field in enumerate(controls):
        complete = ALG.quadratic(field)
        for grade, a, indices in OUTPUTS:
            target = sum(1 << i for i in indices)
            actual = sp.Rational(complete[a].get(target, 0))
            assert coefficient_matrix(a, target, [field])[0, 0] == actual
        if j in (0, 1):
            _, a, indices = OUTPUTS[0 if j == 0 else 16]
            assert complete[a].get(sum(1 << i for i in indices), 0) != 0
        assert all(mask.bit_count() in (1, 5, 9) for item in complete for mask in item)
        print('PASS complete cubic-gradient comparison on18 outputs, control', j)
    weights = sp.Matrix([2, -3])
    mixed = [BC.add(BC.scale(controls[0][a], 2), BC.scale(controls[1][a], -3))
             for a in range(N)]
    complete = ALG.quadratic(mixed)
    for _, a, indices in OUTPUTS:
        target = sum(1 << i for i in indices)
        matrix = coefficient_matrix(a, target, controls[:2])
        assert (weights.T*matrix*weights)[0] == sp.Rational(complete[a].get(target, 0))
    print('PASS mixed-amplitude polarization (2,-3) on18 outputs;72 comparison targets in total.')


def primitive_poly(expression, variable):
    polynomial = sp.Poly(expression, variable, domain=sp.QQ)
    _, polynomial = polynomial.clear_denoms(convert=True)
    _, polynomial = polynomial.primitive()
    if polynomial.LC() < 0:
        polynomial = -polynomial
    return sp.Poly(polynomial.as_expr(), variable, domain=sp.QQ)


def eliminate(matrix):
    _, pivots = matrix.rref()
    kernel = matrix.nullspace()
    assert len(pivots) == 18 and len(kernel) == 3
    assert all(matrix*v == sp.zeros(18, 1) for v in kernel)
    minor = matrix.extract(range(18), pivots).det()
    assert minor != 0
    assert MONOMIALS[18:] == ((4, 4), (4, 5), (5, 5))
    assert sp.Matrix(3, 3, lambda i, j: kernel[j][18+i]) == sp.eye(3)
    print('18 necessary quadratic rows have rank18 on21 raw monomials.')
    print('Monomial order:', MONOMIALS)
    print('Nonzero18-by18 determinant:', minor)
    print('Monomial kernel, free y44,y45,y55:', tuple(tuple(v) for v in kernel))
    for j in range(6):
        if j != 3:
            column = MONOMIALS.index(tuple(sorted((j, 3))))
            assert all(v[column] == 0 for v in kernel)
    coefficient = kernel[2][MONOMIALS.index((3, 3))]
    assert coefficient != 0
    print('Chart x3!=0 excluded by all x3*xj=0 and zero free monomials.')
    print('With x3=0, chart x4=0 excluded by nonzero y55 coefficient in y33:', coefficient)
    z = sp.Symbol('z')

    def relation(pair):
        column = MONOMIALS.index(pair)
        return sp.expand(kernel[0][column]+kernel[1][column]*z+kernel[2][column]*z*z)

    f = primitive_poly(relation((3, 3)), z)
    g = primitive_poly(relation((0, 0))*relation((1, 1))-relation((0, 1))**2, z)
    assert (f.degree(), g.degree()) == (2, 4)
    assert sp.gcd(f, g).degree() == 0
    u, v, h = sp.gcdex(f, g)
    assert h == 1 and u*f+v*g == 1
    denominator = sp.ilcm(*(c.q for p in (u, v) for c in p.all_coeffs()))
    U, V = sp.expand(denominator*u.as_expr()), sp.expand(denominator*v.as_expr())
    assert denominator != 0
    assert sp.expand(U*f.as_expr()+V*g.as_expr()) == denominator
    print('Chart x4!=0, z=x5/x4; actual rank-one minor y00*y11-y01^2 retained.')
    print('F(z)=', f.as_expr())
    print('G(z)=', g.as_expr())
    print('Exact QQ gcd(F,G)=1.')
    print('Integer Bezout U(z)=', U)
    print('Integer Bezout V(z)=', V)
    print('Verified U*F+V*G=', denominator, '!=0.')
    return True


def main():
    start = perf_counter()
    print('Affine-invariant actual H4, flat Z1 lambda1/5,cH2,cV1,t=kappa/24,u0.')
    hooks, fields = reconstruct()
    validate(hooks)
    rows = []
    for grade, a, indices in OUTPUTS:
        target = sum(1 << i for i in indices)
        assert target.bit_count() == grade
        quadratic = coefficient_matrix(a, target, fields)
        row = [quadratic[i, j]*(1 if i == j else 2) for i, j in MONOMIALS]
        rows.append(row)
        print('COEFFICIENT', (grade, a, indices), tuple(row))
    assert eliminate(sp.Matrix(rows))
    assert F(1, 24)+312*F(1, 24)**2 == F(7, 12)
    print('All selected higher equations force zero amplitudes; positive baseline15/(2cV)+7kappa^2/12 excludes zero.')
    print('Only the declared affine-invariant flat-base H4 branch is excluded; no arbitrary-field claim.')
    print('No earlier helper main was run; no private data or output write is used.')
    print(f'Runtime {perf_counter()-start:.3f} seconds.')


if __name__ == '__main__':
    main()
