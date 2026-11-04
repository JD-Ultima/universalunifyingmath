"""Exact affine-invariant hook reduction on the flat positive Z1 model.

Real SO4 tensors are constructed by finite-group reduction and six Lie
algebra equations, with inverse-Gram hooks and all-slot Nomizu derivatives.
The H6 closure line is excluded by exact nonlinear grade-one moments.
BC-affine-hooks.md gives the proof and the affine-invariance boundary.
Running writes no files; no earlier helper or output is imported.
"""
from collections import defaultdict
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations, permutations
from math import gcd, lcm, prod
from pathlib import Path
import json
from time import perf_counter
import sympy as sp
from sympy.polys.matrices import DomainMatrix

N = 14


def fraction(value):
    if isinstance(value, F):
        return value
    value = sp.Rational(value)
    return F(int(value.p), int(value.q))


def mm(a, b):
    return [[sum(a[i][k]*b[k][j] for k in range(4)) for j in range(4)] for i in range(4)]


def bracket(a, b):
    ab, ba = mm(a, b), mm(b, a)
    return [[ab[i][j]-ba[i][j] for j in range(4)] for i in range(4)]


def dot(a, b):
    return sum(a[i][j]*b[j][i] for i in range(4) for j in range(4))


def plus(out, key, value):
    if value:
        out[key] = out.get(key, 0)+value
        if not out[key]:
            del out[key]


def setup():
    sym, chars = [], [1 << i for i in range(4)]
    for i, j in combinations(range(4), 2):
        s = [[0]*4 for _ in range(4)]
        s[i][j] = s[j][i] = 1
        sym.append(s)
        chars.append((1 << i) ^ (1 << j))
    for k in range(1, 4):
        diag = [1]*k+[-k]+[0]*(3-k)
        sym.append([[diag[i] if i == j else 0 for j in range(4)] for i in range(4)])
        chars.append(0)
    sym.append([[int(i == j) for j in range(4)] for i in range(4)])
    chars.append(0)
    norms = [dot(s, s) for s in sym]
    metric = [F(2)]*4+[F(4*n) for n in norms[:-1]]+[F(16, 5)]

    def coords(s):
        return [F(dot(s, u), norm) for u, norm in zip(sym, norms)]

    matrices = []
    for i, j in combinations(range(4), 2):
        J = [[0]*4 for _ in range(4)]
        J[i][j], J[j][i] = 1, -1
        M = [[F(0)]*N for _ in range(N)]
        for a in range(4):
            for b in range(4):
                M[a][b] = F(J[a][b])
        for b, S in enumerate(sym):
            for a, value in enumerate(coords(bracket(J, S))):
                M[4+a][4+b] = value
        matrices.append(M)
    for i in range(4):
        M = [[F(0)]*N for _ in range(N)]
        for j in range(4):
            S = [[F(0)]*4 for _ in range(4)]
            S[i][j] += F(1, 4)
            S[j][i] += F(1, 4)
            if i == j:
                for a in range(4):
                    S[a][a] += F(1, 2)  # (mu-1)/8 at mu=5
            for a, value in enumerate(coords(S)):
                M[4+a][j] = value
        for a, S in enumerate(sym):
            for b in range(4):
                M[b][4+a] = -F(S[b][i])
        matrices.append(M)
    assert all(metric[a]*M[a][b]+metric[b]*M[b][a] == 0
               for M in matrices for a in range(N) for b in range(N))
    # On covectors, an input e^a is sent to -sum_b M[a,b]e^b.
    cov = [tuple(tuple((b, -M[a][b]) for b in range(N) if M[a][b]) for a in range(N))
           for M in matrices]
    return metric, chars, cov


METRIC, CHARS, COV = setup()


def bit_indices(mask):
    return tuple(i for i in range(N) if mask >> i & 1)


def char_of(a, mask):
    value = CHARS[a]
    for i in bit_indices(mask):
        value ^= CHARS[i]
    return value


@lru_cache(None)
def form_action(which, mask):
    out = {}
    for i in bit_indices(mask):
        reduced = mask ^ (1 << i)
        remove_sign = (-1)**((mask & ((1 << i)-1)).bit_count())
        for j, coefficient in COV[which][i]:
            if not reduced >> j & 1:
                insert_sign = (-1)**((reduced & ((1 << j)-1)).bit_count())
                plus(out, reduced | (1 << j), coefficient*remove_sign*insert_sign)
    return tuple(out.items())


def action_basis(which, a, mask):
    out = {}
    for b, coefficient in COV[which][a]:
        plus(out, (b, mask), coefficient)
    for other, coefficient in form_action(which, mask):
        plus(out, (a, other), coefficient)
    return out


def action(which, tensor):
    out = {}
    for (a, mask), value in tensor.items():
        for key, coefficient in action_basis(which, a, mask).items():
            plus(out, key, value*coefficient)
    return out


def primitive(tensor):
    den = lcm(*(F(v).denominator for v in tensor.values()))
    out = {key: int(v*den) for key, v in tensor.items() if v}
    common = gcd(*out.values())
    out = {key: F(v//common) for key, v in out.items()}
    if out[min(out)] < 0:
        out = {key: -v for key, v in out.items()}
    return out


def nullspace(rows, columns):
    loc = {key: j for j, key in enumerate(columns)}
    data = {i: {loc[key]: sp.Rational(value.numerator, value.denominator)
                for key, value in row.items()} for i, row in enumerate(rows)}
    dm = DomainMatrix.from_dict_sympy(len(rows), len(columns), data).convert_to(sp.QQ)
    reduced, pivots = dm.rref()
    reduced = reduced.to_sdm()
    basis = []
    for free in range(len(columns)):
        if free in pivots:
            continue
        vector = {columns[free]: F(1)}
        for i, pivot in enumerate(pivots):
            value = reduced.get(i, {}).get(free, 0)
            if value:
                vector[columns[pivot]] = -fraction(value)
        basis.append(primitive(vector))
    return len(pivots), basis


def finite_actions():
    pairs = list(combinations(range(4), 2))
    index = {pair: j for j, pair in enumerate(pairs)}
    actions = []
    for perm in permutations(range(4)):
        odd = sum(perm[i] > perm[j] for i, j in combinations(range(4), 2)) % 2
        signs = [-1 if odd and perm[i] == 0 else 1 for i in range(4)]
        images = [(perm[i], signs[i]) for i in range(4)]
        images += [(4+index[tuple(sorted((perm[i], perm[j])))], signs[i]*signs[j]) for i, j in pairs]
        images += [(10+perm[i], 1) for i in range(4)]
        actions.append(tuple(images))
    return actions


FINITE = finite_actions()


def finite_image(images, a, mask):
    new_a, sign = images[a]
    ordered = []
    for i in bit_indices(mask):
        j, factor = images[i]
        sign *= factor
        ordered.append(j)
    sign *= (-1)**sum(ordered[i] > ordered[j] for i, j in combinations(range(len(ordered)), 2))
    return (new_a, sum(1 << i for i in ordered)), sign


def finite_basis(k):
    """SO4 signed-permutation invariants in monomial H/O/Eii coordinates."""
    masks = [sum(1 << i for i in indices) for indices in combinations(range(N), k)]
    columns = [(a, mask) for a in range(N) for mask in masks if char_of(a, mask) in (0, 15)]
    seen, basis = set(), []
    for key in columns:
        if key in seen:
            continue
        average = {}
        for images in FINITE:
            dest, sign = finite_image(images, *key)
            seen.add(dest)
            plus(average, dest, sign)
        if average:
            tensor = primitive(average)
            for images in FINITE:
                transformed = {}
                for state, value in tensor.items():
                    dest, sign = finite_image(images, *state)
                    plus(transformed, dest, value*sign)
                assert transformed == tensor
            basis.append(tensor)
    assert seen == set(columns)
    print('Even-sign candidate columns', len(columns), 'finite subgroup invariant dimension', len(basis), flush=True)
    return basis


# Old Eii diagonal covectors expanded in the orthogonal D1,D2,D3,I basis.
DIAGONAL_CHANGE = ((1, 1, 1, 1), (-1, 1, 1, 1), (0, -2, 1, 1), (0, 0, -3, 1))


def coform_change(i):
    return ((i, 1),) if i < 10 else tuple((10+j, v) for j, v in enumerate(DIAGONAL_CHANGE[i-10]) if v)


@lru_cache(None)
def change_form(mask):
    out = {0: 1}
    for i in bit_indices(mask):
        after = {}
        for previous, value in out.items():
            for j, coefficient in coform_change(i):
                if not previous >> j & 1:
                    sign = (-1)**((previous >> (j+1)).bit_count())
                    plus(after, previous | (1 << j), sign*coefficient*value)
        out = after
    return tuple(out.items())


def change_tensor(tensor):
    out = {}
    for (a, mask), value in tensor.items():
        for b, coefficient in coform_change(a):
            for other, extra in change_form(mask):
                plus(out, (b, other), value*coefficient*extra)
    return out


def invariant_basis(k):
    finite = [change_tensor(tensor) for tensor in finite_basis(k)]
    equations = defaultdict(dict)
    for j, tensor in enumerate(finite):
        for J in range(6):
            for (a, mask), value in action(J, tensor).items():
                plus(equations[J, a, mask], (j, 0), value)
    rank, coefficients = nullspace(list(equations.values()), [(j, 0) for j in range(len(finite))])
    invariant = [primitive(combine(v, finite)) for v in coefficients]
    assert all(not action(J, v) for v in invariant for J in range(6))
    print('Exact six-generator invariance rank', rank, 'invariant dimension', len(invariant), flush=True)
    return invariant


def hook_traces(tensor):
    alt, ctr = {}, {}
    for (a, mask), value in tensor.items():
        sign = (-1)**((mask & ((1 << a)-1)).bit_count())
        if mask >> a & 1:
            plus(ctr, mask ^ (1 << a), value*sign/METRIC[a])
        else:
            plus(alt, mask | (1 << a), value*sign)
    return alt, ctr


def combine(coefficients, basis):
    out = {}
    for (j, _), coefficient in coefficients.items():
        for key, value in basis[j].items():
            plus(out, key, coefficient*value)
    return out


def hook_basis(invariant):
    equations = defaultdict(dict)
    for j, tensor in enumerate(invariant):
        alt, ctr = hook_traces(tensor)
        for kind, field in (('Alt', alt), ('Ctr', ctr)):
            for mask, coefficient in field.items():
                plus(equations[kind, mask], (j, 0), coefficient)
    rank, coefficients = nullspace(list(equations.values()), [(j, 0) for j in range(len(invariant))])
    basis = [primitive(combine(vector, invariant)) for vector in coefficients]
    assert all(hook_traces(v) == ({}, {}) for v in basis)
    assert all(not action(J, v) for v in basis for J in range(6))
    print('Invariant hook constraints rank', rank, 'hook dimension', len(basis), flush=True)
    return basis


def hodge_scaled(tensor):
    """Hodge star divided by sqrt(det G), a nonzero constant factor."""
    out = {}
    full = (1 << N)-1
    for (a, mask), value in tensor.items():
        indices = bit_indices(mask)
        parity = sum(i-j for j, i in enumerate(indices))
        factor = F((-1)**parity)
        for i in indices:
            factor /= METRIC[i]
        out[a, full ^ mask] = value*factor
    return out


def exterior_covariant(tensor):
    out = {}
    for r in range(4):  # Lambda_S=0 for the ten vertical directions.
        for (a, mask), value in action(6+r, tensor).items():
            if not mask >> r & 1:
                sign = (-1)**((mask & ((1 << r)-1)).bit_count())
                plus(out, (a, mask | (1 << r)), sign*value)
    return out


def codifferential(tensor):
    out = {}
    for r in range(4):
        for (a, mask), value in action(6+r, tensor).items():
            if mask >> r & 1:
                sign = (-1)**((mask & ((1 << r)-1)).bit_count())
                plus(out, (a, mask ^ (1 << r)), -sign*value/METRIC[r])
    return out


def rank_witness(images):
    keys = sorted(set().union(*(v.keys() for v in images)))
    data = {i: {j: sp.Rational(images[j][key].numerator, images[j][key].denominator)
                for j in range(len(images)) if key in images[j]} for i, key in enumerate(keys)}
    dm = DomainMatrix.from_dict_sympy(len(keys), len(images), data).convert_to(sp.QQ)
    rank = len(dm.rref()[1])
    selected = dm.transpose().rref()[1]
    assert len(selected) == rank
    pivots = dm.extract(selected, range(len(images))).rref()[1]
    witness = dm.extract(selected, pivots)
    determinant = witness.det()
    assert determinant != 0
    kernel_rank, kernel = nullspace([{(j, 0): fraction(v) for j, v in data[i].items()} for i in data],
                                  [(j, 0) for j in range(len(images))])
    assert kernel_rank == rank
    return rank, tuple(keys[i] for i in selected), tuple(pivots), determinant, kernel



def add(*values):
    result = {}
    for value in values:
        for mask, coefficient in value.items():
            result[mask] = result.get(mask, 0) + coefficient
    return {mask: value for mask, value in result.items() if value}


def scale(value, factor):
    return {mask: coefficient * factor for mask, coefficient in value.items()
            if coefficient * factor}


class Clifford:
    def __init__(self, gram):
        self.gram = tuple(map(F, gram))

    @lru_cache(None)
    def factor(self, a, b):
        parity = sum((b & ((1 << j) - 1)).bit_count()
                     for j in range(N) if a >> j & 1)
        return F((-1) ** parity, prod(self.gram[j] for j in range(N)
                                     if a & b & (1 << j)))

    def mul(self, left, right):
        result = {}
        for a, av in left.items():
            for b, bv in right.items():
                key = a ^ b
                result[key] = result.get(key, 0) + av * bv * self.factor(a, b)
        return {key: value for key, value in result.items() if value}

    def comm(self, a, b):
        return add(self.mul(a, b), scale(self.mul(b, a), -1))

    def anti(self, a, b):
        return add(self.mul(a, b), self.mul(b, a))

    def quadratic(self, field):
        """Return all-covariant coefficients of Q(iC)/i for real C."""
        output = []
        for a in range(N):
            total = {}
            for b, c in combinations(range(N), 2):
                if a in (b, c) or (not field[b] and not field[c]):
                    continue
                triple = self.mul(self.mul({1 << a: 1}, {1 << b: 1}),
                                  {1 << c: 1})
                term = add(self.anti(triple, self.comm(field[b], field[c])),
                           self.comm(self.anti(triple, field[b]), field[c]),
                           scale(self.comm(self.anti(triple, field[c]), field[b]), -1))
                total = add(total, term)
            output.append(scale(total, -self.gram[a] / 3))
        return output



@lru_cache(None)
def sign(a, b):
    return (-1) ** sum((b & ((1 << j) - 1)).bit_count()
                      for j in range(N) if a >> j & 1)


def quadratic_coefficient(field, gram, a, target):
    """One exact coefficient of the complete cubic gradient, without projection assumptions."""
    total = F(0)
    for b, c in combinations(range(N), 2):
        if a in (b, c):
            continue
        m = (1 << a) | (1 << b) | (1 << c)
        triple_sign = sign(1 << a, 1 << b) * sign((1 << a) | (1 << b), 1 << c)
        for x, vx in field[b].items():
            y = m ^ x ^ target
            vy = field[c].get(y, 0)
            if not vy:
                continue
            parity_sum = (sign(m, x) * sign(m ^ x, y)
                          - sign(m, y) * sign(m ^ y, x)
                          + sign(x, m) * sign(x ^ m, y)
                          + sign(x, y) * sign(x ^ y, m)
                          - sign(y, m) * sign(y ^ m, x)
                          - sign(y, x) * sign(y ^ x, m))
            if parity_sum:
                repeated = (m & x) | (m & y) | (x & y)
                denominator = prod(gram[j] for j in range(N) if repeated >> j & 1)
                total += F(2 * triple_sign * parity_sum * vx * vy, denominator)
    return -gram[a] * total / 3



def main():
    started = perf_counter()
    print('Flat positive Z1: lambda=1/5, cH=2, cV=1; all tensor slots covariant.')
    print('Gram diagonal:', tuple(map(str, METRIC)))
    closed_h6 = None
    for degree, expected in ((4, (20, 13, 7, 6)), (6, (30, 23, 22, 1))):
        print('DEGREE', degree, flush=True)
        invariant = invariant_basis(degree)
        basis = hook_basis(invariant)
        assert (len(invariant), len(basis)) == expected[:2]
        if degree == 4:
            dual = [hodge_scaled(v) for v in basis]
            assert all(hook_traces(v) == ({}, {}) for v in dual)
            assert all(not action(j, v) for v in dual for j in range(6))
            for v in basis:
                assert hodge_scaled(hodge_scaled(v)) == {key: c / prod(METRIC) for key, c in v.items()}
                assert all(action(6+r, hodge_scaled(v)) == hodge_scaled(action(6+r, v)) for r in range(4))
                assert exterior_covariant(hodge_scaled(v)) == hodge_scaled(codifferential(v))
            images = [exterior_covariant(v) for v in dual]
        else:
            images = [exterior_covariant(v) for v in basis]
        rank, rows, columns, determinant, kernel = rank_witness(images)
        assert (rank, len(kernel)) == expected[2:]
        closed = [primitive(combine(v, basis)) for v in kernel]
        assert all(not (codifferential(v) if degree == 4 else exterior_covariant(v)) for v in closed)
        print('Closure rank:', rank, '; kernel dimension:', len(kernel))
        print('Nonzero minor rows (external index, form bitmask):', rows)
        print('Nonzero minor columns:', columns)
        print('Nonzero minor determinant:', determinant)
        print('Kernel amplitude vectors:', [tuple((j, str(c)) for (j, _), c in sorted(v.items())) for v in kernel])
        print('All six isotropy generators, both weighted hook traces and full derivative slots verified.')
        if degree == 4:
            print('H4 uses dual H10 closure: Hodge square, connection commutation and delta equivalence verified.')
        else:
            closed_h6 = closed[0]
    certificate = json.loads(Path(__file__).with_name('BC_h6_closed_generator.json').read_text(encoding='utf-8'))
    assert certificate['degree'] == 6
    assert tuple(map(F, certificate['metric_gram'])) == tuple(METRIC)
    explicit = {(a, mask): F(c) for a, mask, c in certificate['entries']}
    assert len(explicit) == len(certificate['entries']) == 1358
    assert explicit == closed_h6
    print('Independently reconstructed H6 generator equals the explicit JSON certificate.')
    field = [{} for _ in range(N)]
    for (a, mask), c in closed_h6.items():
        field[a][mask] = c
    norm = sum(c*c / (METRIC[a] * prod(METRIC[j] for j in bit_indices(mask)))
               for (a, mask), c in closed_h6.items())
    algebra = Clifford(METRIC)
    volume = {(1 << N)-1: 1}
    determinant = prod(METRIC)
    assert algebra.mul(volume, volume) == {0: -1 / determinant}
    converted = [algebra.mul(v, volume) for v in field]
    moments = [determinant * quadratic_coefficient(converted, METRIC, a, 1 << a) / METRIC[a]
               for a in range(N)]
    assert len(set(moments[:4])) == len(set(moments[4:13])) == 1
    assert sum(moments) == -F(52, 3) * norm
    ratios = tuple(c / norm for c in (moments[0], moments[4], moments[13]))
    assert ratios == (-F(13589, 13023), -F(162854, 117207), -F(8522, 13023))
    assert moments[0] != moments[4]
    print('H6 closed generator: nonzero coefficients', len(closed_h6), '; norm squared', norm)
    print('Q1 metric-relative diagonal on H4, Sym0(9), I:', tuple(map(str, (moments[0], moments[4], moments[13]))))
    print('Norm-normalized diagonal:', tuple(map(str, ratios)))
    print('All14 diagonal coefficients and trace Q1=-(52/3)*norm squared verified.')
    control = [{} for _ in range(N)]
    control[0][sum(1 << j for j in range(8))] = F(2, 3)
    control[1][sum(1 << j for j in range(1, 9))] = -F(3, 5)
    complete = algebra.quadratic(control)
    for a in range(N):
        for target in (1 << a, 1 << ((a+3) % N), 31):
            assert quadratic_coefficient(control, METRIC, a, target) == complete[a].get(target, 0)
    print('42 selected coefficients agree with the complete cubic-gradient evaluator.')
    print('Affine-invariant pure H6Gamma: nonzero fields have anisotropic Q1; zero field has positive baseline.')
    print('H4 remains a six-dimensional linear candidate family; nonlinear sufficiency is not asserted.')
    print('No exclusion of arbitrary nonhomogeneous fields, other parameters or coincident kernels.')
    print(f'Runtime {perf_counter()-started:.3f} seconds.')


if __name__ == '__main__':
    main()
