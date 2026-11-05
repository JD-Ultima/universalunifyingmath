"""Exact affine-invariant actual-H2 obstruction for every positive Z1 parameter.

Uses public BC construction as a library, no saved/private basis. The real
field is B2, with no i or Gamma. Only this new main runs; T/xcheck is used
as a library for small independently implemented convention comparisons.
"""
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations
from math import prod
from pathlib import Path
import sys
from time import perf_counter
import sympy as sp

sys.path.insert(0, str(Path(__file__).resolve().parent))
import BC_affine_hooks as BC
import T_singular_residual as T

N = 14
MU = sp.Symbol('mu', positive=True)
GRAM = tuple(sp.sympify(v) for v in BC.METRIC[:-1])+(16/MU,)


def normalize(value):
    return {key: coefficient for key, raw in value.items()
            if (coefficient := sp.cancel(raw)) != 0}


class ParameterClifford(BC.Clifford):
    def __init__(self, gram):
        self.gram = tuple(map(sp.sympify, gram))

    @lru_cache(None)
    def factor(self, left, right):
        parity = sum((right & ((1 << j)-1)).bit_count()
                     for j in range(N) if left >> j & 1)
        return sp.Integer((-1)**parity)/prod(self.gram[j] for j in range(N)
                                            if left & right & (1 << j))


ALG = ParameterClifford(GRAM)
LAMBDA = []
for r in range(4):
    matrix = sp.Matrix(N, N, lambda a, b: -dict(BC.COV[6+r][a]).get(b, 0))
    matrix[13, r] = MU/8
    assert all(sp.cancel(v) == 0 for v in sp.diag(*GRAM)*matrix+matrix.T*sp.diag(*GRAM))
    LAMBDA.append(matrix)


def field(tensor):
    return [{mask: value for (b, mask), value in tensor.items() if b == a}
            for a in range(N)]


def derivative_from_jets(jets, gram, algebra=None):
    """Full D_WA(B)/i: no projection and no extra inverse b,c Gram factors."""
    algebra = BC.Clifford(gram) if algebra is None else algebra
    exterior = {(b, c): BC.add(jets[b][c], BC.scale(jets[c][b], -1))
                for b, c in combinations(range(N), 2)}
    output = []
    for a in range(N):
        total = {}
        for (b, c), value in exterior.items():
            if a in (b, c) or not value:
                continue
            triple = algebra.mul(algebra.mul({1 << a: 1}, {1 << b: 1}), {1 << c: 1})
            total = BC.add(total, algebra.anti(triple, value))
        output.append(normalize(BC.scale(total, gram[a])))
    return output


def parameter_action(matrix, tensor):
    """Own immutable-parameter action; never modifies BC globals or caches."""
    output = {}
    for (a, mask), value in tensor.items():
        for b in range(N):
            if matrix[a, b]:
                BC.plus(output, (b, mask), -matrix[a, b]*value)
        for i in BC.bit_indices(mask):
            rest = mask ^ (1 << i)
            removed = (-1)**((mask & ((1 << i)-1)).bit_count())
            for j in range(N):
                if matrix[i, j] and not rest >> j & 1:
                    inserted = (-1)**((rest & ((1 << j)-1)).bit_count())
                    BC.plus(output, (a, rest | (1 << j)), -matrix[i, j]*value*removed*inserted)
    return normalize(output)


def homogeneous_derivative(tensor):
    jets = [field(parameter_action(LAMBDA[r], tensor)) if r < 4 else [{} for _ in range(N)]
            for r in range(N)]
    return derivative_from_jets(jets, GRAM, ALG)


def quadratic(value):
    # BC returns Q(iC)/i. Real actual B has the opposite quadratic phase.
    return [normalize(BC.scale(v, -1)) for v in ALG.quadratic(value)]


def to_orthonormal(value, scales, extra=1):
    return [{mask: sp.Rational(coefficient)/(extra*scales[a]*prod(scales[j] for j in BC.bit_indices(mask)))
             for mask, coefficient in item.items()} for a, item in enumerate(value)]


def to_raw_divide_i(value, scales):
    return [{mask: sp.expand(coefficient*scales[a]*prod(scales[j] for j in BC.bit_indices(mask))/sp.I)
             for mask, coefficient in item.items() if sp.expand(coefficient) != 0}
            for a, item in enumerate(value)]


def compare(left, right):
    for a in range(N):
        for mask in set(left[a]) | set(right[a]):
            assert sp.expand(left[a].get(mask, 0)-right[a].get(mask, 0)) == 0


def conventions():
    scales = tuple(1+j % 3 for j in range(N))
    gram = tuple(F(s*s) for s in scales)
    real_hook = [{} for _ in range(N)]
    for a, indices, value in ((0, (1, 2), F(1)), (1, (0, 2), F(1)),
                               (0, (3, 4), F(2, 3)), (3, (0, 4), F(2, 3))):
        real_hook[a][sum(1 << j for j in indices)] = value
    single = [{} for _ in range(N)]
    single[0][(1 << 1) | (1 << 2)] = F(3, 7)
    for j, (value, direction) in enumerate(((real_hook, 5), (single, 3))):
        jets = [[{} for _ in range(N)] for _ in range(N)]
        jets[direction] = value
        actual = derivative_from_jets(jets, gram)
        independent = T.derivative(to_orthonormal(value, scales, scales[direction]),
                                   direction, (1, 0, 0, 0))
        compare(actual, to_raw_divide_i(independent, scales))
        grade5 = any(mask.bit_count() == 5 for item in actual for mask in item)
        assert grade5 == (j == 1)
        print('PASS independent full derivative/non-unit-weight comparison, control', j)
    orthonormal = to_orthonormal(real_hook, scales)
    assert T.antihermitian(orthonormal)
    independent = T.quadratic(orthonormal, (1, 0, 0, 0))
    actual = [BC.scale(item, -1) for item in BC.Clifford(gram).quadratic(real_hook)]
    compare(actual, to_raw_divide_i(independent, scales))
    assert any(mask.bit_count() == 5 for item in actual for mask in item)
    print('PASS independent real-B quadratic phase/complete-output comparison, including nonzero grade5.')


def main():
    started = perf_counter()
    conventions()
    invariant = BC.invariant_basis(2)
    basis = BC.hook_basis(invariant)
    assert len(invariant) == 3 and len(basis) == 2
    # Completeness for ALL mu: only inverse weight13 changes, and it never
    # occurs on the coefficient slot of ANY of the three raw invariants.
    assert all(a != 13 for tensor in invariant for a, _ in tensor)
    print('All three raw SO4 invariants have no external13 component; hook constraints are unchanged for every mu.')
    for j, tensor in enumerate(basis):
        assert BC.hook_traces(tensor) == ({}, {})
        assert all(not BC.action(g, tensor) for g in range(6))
        alt, ctr = {}, {}
        for (a, mask), value in tensor.items():
            sign = (-1)**((mask & ((1 << a)-1)).bit_count())
            if mask >> a & 1:
                BC.plus(ctr, mask ^ (1 << a), value*sign/GRAM[a])
            else:
                BC.plus(alt, mask | (1 << a), value*sign)
        assert not normalize(alt) and not normalize(ctr)
        print('BASIS', j, [(a, BC.bit_indices(mask), str(v)) for (a, mask), v in sorted(tensor.items())])
    values = [field(v) for v in basis]
    derivatives = [homogeneous_derivative(v) for v in basis]
    assert all(mask.bit_count() == 1 for output in derivatives for item in output for mask in item)
    print('Full all-slot derivative computed first; grade5 cancels on both hook basis directions.')
    squares = [quadratic(v) for v in values]
    total_square = quadratic([BC.add(a, b) for a, b in zip(*values)])
    cross = [normalize(BC.add(total_square[a], BC.scale(squares[0][a], -1), BC.scale(squares[1][a], -1)))
             for a in range(N)]
    assert all(mask.bit_count() in (1, 5) for output in squares+[cross] for item in output for mask in item)
    x, y = sp.symbols('x y', real=True)
    d, q = [], []
    for a in (0, 4, 13):
        gram = GRAM[a]
        d.append(sp.expand((x*derivatives[0][a].get(1 << a, 0)
                           +y*derivatives[1][a].get(1 << a, 0))/gram))
        q.append(sp.expand((x*x*squares[0][a].get(1 << a, 0)
                           +y*y*squares[1][a].get(1 << a, 0)
                           +x*y*cross[a].get(1 << a, 0))/gram))
    assert d == [(-9*x+27*MU*y)/16, x/4-MU*y, 9*MU*y/4]
    assert q == [-(3*x*x+51*MU*y*y)/16, -(11*x*x+109*MU*y*y)/24,
                 -(3*x*x+117*MU*y*y)/8]
    assert [v.subs(MU, 5) for v in d] == [(-9*x+135*y)/16, x/4-5*y, 45*y/4]
    assert [v.subs(MU, 5) for v in q] == [-(3*x*x+255*y*y)/16,
                                         -(11*x*x+545*y*y)/24, -(3*x*x+585*y*y)/8]
    # Check every grade-one entry, including off-diagonals, on all fourteen slots.
    groups = (0,)*4+(1,)*9+(2,)
    for a in range(N):
        masks = set(derivatives[0][a]) | set(derivatives[1][a]) | set(squares[0][a]) | set(squares[1][a]) | set(cross[a])
        for mask in masks | {1 << a}:
            if mask.bit_count() != 1:
                continue
            actual_d = sp.expand(x*derivatives[0][a].get(mask, 0)+y*derivatives[1][a].get(mask, 0))
            actual_q = sp.expand(x*x*squares[0][a].get(mask, 0)
                                 +y*y*squares[1][a].get(mask, 0)+x*y*cross[a].get(mask, 0))
            assert sp.cancel(actual_d-(GRAM[a]*d[groups[a]] if mask == 1 << a else 0)) == 0
            assert sp.cancel(actual_q-(GRAM[a]*q[groups[a]] if mask == 1 << a else 0)) == 0
    target5 = sum(1 << j for j in (0, 1, 2, 4, 5))
    selected5 = sp.expand(x*x*squares[0][0].get(target5, 0)+y*y*squares[1][0].get(target5, 0)
                          +x*y*cross[0].get(target5, 0))
    assert selected5 == -sp.Rational(4, 3)*x*x
    norm = sp.Rational(9, 16)*(x*x+13*MU*y*y)
    actual_norm = 0
    for a in range(N):
        for mask in set(values[0][a]) | set(values[1][a]):
            coefficient = (x*sp.Rational(values[0][a].get(mask, 0))
                           +y*sp.Rational(values[1][a].get(mask, 0)))
            denominator = GRAM[a]*prod(GRAM[j] for j in BC.bit_indices(mask))
            actual_norm += coefficient*coefficient/denominator
    assert sp.expand(actual_norm) == norm
    assert sp.expand(4*d[0]+9*d[1]+d[2]) == 0
    assert sp.expand(4*q[0]+9*q[1]+q[2]) == -sp.Rational(28, 3)*norm
    print('Metric-relative derivative H,Sym0,I:', tuple(d))
    print('Metric-relative quadratic H,Sym0,I:', tuple(q))
    print('All fourteen grade1 diagonals and every off-diagonal checked.')
    print('Selected genuine grade5 coefficient (external0;0,1,2,4,5):', selected5)
    print('Trace checks: trD=0, trQ=-(28/3)N, N=(9/16)(x^2+13mu*y^2).')
    rho = sp.Symbol('rho', positive=True)
    ric = (-MU/4, -sp.Rational(5, 4), -MU/4)
    scalar = 5*ric[0]+9*ric[1]
    einstein = tuple(sp.expand(v-scalar/2) for v in ric)
    assert einstein == (3*(MU+15)/8, 5*(MU+7)/8, 3*(MU+15)/8)
    residual = tuple(sp.expand(einstein[j]+sp.Rational(93, 98)*rho+d[j]+q[j]) for j in range(3))
    hi = sp.factor((residual[0]-residual[2]).subs(x, 0))
    assert sp.expand(hi-3*MU*y*(61*y-3)/16) == 0
    candidate0 = sp.factor(residual[0].subs({x: 0, y: 0}))
    candidate1 = sp.factor(residual[0].subs({x: 0, y: sp.Rational(3, 61)}))
    positive_part = sp.simplify(candidate1-candidate0)
    assert positive_part == sp.Rational(2241, 29768)*MU
    assert candidate0.is_positive is True and candidate1.is_positive is True
    assert F(3, 56)+312*F(3, 56)**2 == F(93, 98)
    assert not any(v.atoms(sp.Float) for data in derivatives+squares+[cross]
                   for item in data for v in item.values())
    print('Normalized Einstein H,Sym0,I:', einstein, '; rho=cV*kappa^2>0.')
    print('After grade5 forces x=0, H-I=', hi, '; hence y=0 or3/61.')
    print('Horizontal residual at y=0:', candidate0, '>0.')
    print('Horizontal residual at y=3/61:', candidate1, '>0; extra contribution=', positive_part)
    print('Every exact coefficient is free of Float; mu=5 calibration agrees with the original tables.')
    print('All positive observer scales use B/cV and rho=cV*kappa^2; no parameter sampling or solve is used.')
    print('Only affine-invariant actual H2 on flat positive Z1 is excluded, not arbitrary fields or curved bases.')
    print('No old helper main was run; basis and equations use no private data or saved receipts.')
    print(f'Runtime {perf_counter()-started:.3f} seconds.')


if __name__ == '__main__':
    main()
