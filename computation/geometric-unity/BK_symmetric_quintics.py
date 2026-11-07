"""Exact scalar-Gram obstructions for two complete symmetric quintic families.

The coordinate-S5 family in R5 has four harmonic amplitudes. The standard-S6
family on the sum-zero subspace of R6 has one. All contractions use
T=(1/120)*d^5 F and the Euclidean metric. This helper verifies these finite
families, not the full 91-dimensional harmonic-quintic space or a field PDE.
It uses only SymPy and the standard library, and writes only to stdout.
"""

from fractions import Fraction
from functools import reduce
from itertools import combinations_with_replacement, permutations, product
from math import factorial, gcd, lcm, prod
from time import perf_counter

import sympy as s


PARTITIONS = (
    (5, 0, 0, 0, 0),
    (4, 1, 0, 0, 0),
    (3, 2, 0, 0, 0),
    (3, 1, 1, 0, 0),
    (2, 2, 1, 0, 0),
    (2, 1, 1, 1, 0),
    (1, 1, 1, 1, 1),
)
HARMONIC_COLUMNS = (
    (4, 5, -10, 0, 0, 0, 0),
    (0, 1, 0, 0, -2, 0, 0),
    (0, 0, 0, 2, 0, -9, 0),
    (0, 0, 0, 0, 0, 0, 1),
)
GRAM_ROWS = (
    (46800, 5280, 0, 0, 384, 0, 0, -2070, 0, -1),
    (1200, 240, 240, 0, 16, 72, 0, 54, -3, 0),
    (600, 480, -240, 0, 56, 216, 0, 264, -9, 0),
    (0, 0, 660, 0, 0, 84, -2, 111, 0, 0),
)


def primitive(vector):
    """Normalize a nonzero rational row by a constant rational factor."""
    denominator = lcm(*(value.q for value in vector))
    integers = [int(value * denominator) for value in vector]
    common = reduce(gcd, integers)
    integers = [value // common for value in integers]
    if next(value for value in integers if value) < 0:
        integers = [-value for value in integers]
    return tuple(integers)


def rational_identity(equations, variables, degree):
    """Construct h_i with sum h_i*f_i=1 and verify by exact expansion."""
    if not variables:
        index = next(i for i, value in enumerate(equations) if value != 0)
        multipliers = [s.Integer(0)] * len(equations)
        multipliers[index] = 1 / equations[index]
    else:
        powers = []
        for total in range(degree + 1):
            for indices in combinations_with_replacement(range(len(variables)), total):
                powers.append(tuple(indices.count(i) for i in range(len(variables))))
        monomials = [prod(v**p for v, p in zip(variables, power)) for power in powers]
        columns = [
            s.Poly(equation * monomial, *variables, domain=s.QQ).as_dict()
            for equation in equations
            for monomial in monomials
        ]
        zero = (0,) * len(variables)
        output_powers = sorted({zero}.union(*(set(column) for column in columns)))
        matrix = s.Matrix(
            [[column.get(power, 0) for column in columns] for power in output_powers]
        )
        target = s.Matrix([int(power == zero) for power in output_powers])
        solution, parameters = matrix.gauss_jordan_solve(target)
        solution = solution.subs({parameter: 0 for parameter in parameters})
        assert matrix * solution == target
        multipliers = [
            s.expand(
                sum(
                    solution[i * len(monomials) + j] * monomial
                    for j, monomial in enumerate(monomials)
                )
            )
            for i in range(len(equations))
        ]
    assert all(not value.atoms(s.Float) for value in multipliers)
    assert s.expand(sum(h * f for h, f in zip(multipliers, equations))) == 1
    return multipliers


def coordinate_s5():
    x = s.symbols("x0:5")
    orbits = [
        sum(
            prod(v**power for v, power in zip(x, powers))
            for powers in sorted(set(permutations(partition)))
        )
        for partition in PARTITIONS
    ]
    laps = [s.Poly(sum(s.diff(value, v, 2) for v in x), *x) for value in orbits]
    cubic_monomials = sorted(set().union(*(set(value.monoms()) for value in laps)))
    trace_matrix = s.Matrix(
        [[value.coeff_monomial(monomial) for value in laps] for monomial in cubic_monomials]
    )
    kernel = tuple(primitive(value) for value in trace_matrix.nullspace())
    assert trace_matrix.rank() == 3 and kernel == HARMONIC_COLUMNS
    basis = [
        s.Poly(s.expand(sum(c * p for c, p in zip(vector, orbits))), *x, domain=s.QQ)
        for vector in kernel
    ]
    for polynomial in basis:
        expression = polynomial.as_expr()
        assert s.expand(sum(s.diff(expression, v, 2) for v in x)) == 0
        assert s.expand(expression - expression.xreplace(dict(zip(x, x[1:] + x[:1])))) == 0
        assert s.expand(expression - expression.xreplace({x[0]: x[1], x[1]: x[0]})) == 0
        assert not expression.atoms(s.Float)

    # The seven degree-five partitions give all coordinate-permutation orbit sums.
    print("Coordinate S5 on R5: 7 invariant quintics, trace rank 3, harmonic dimension 4.")
    print("Monomial-symmetric partition order:", PARTITIONS)
    print("Harmonic basis columns:", HARMONIC_COLUMNS)

    indices = tuple(product(range(5), repeat=5))
    tensors = []
    for polynomial in basis:
        coefficients = {
            powers: int(value) * prod(factorial(k) for k in powers)
            for powers, value in polynomial.terms()
        }
        tensor = {
            index: coefficients.get(tuple(index.count(i) for i in range(5)), 0)
            for index in indices
        }
        assert all(
            sum(tensor[(i, i) + rest] for i in range(5)) == 0
            for rest in product(range(5), repeat=3)
        )
        tensors.append(tensor)
    norm_gram = s.Matrix(
        4, 4, lambda i, j: sum(tensors[i][key] * tensors[j][key] for key in indices)
    )
    assert norm_gram == s.Matrix(
        [[5472000, 288000, 0, 0], [288000, 115200, 0, 0],
         [0, 0, 475200, 0], [0, 0, 0, 120]]
    )
    assert all(norm_gram[:i, :i].det() > 0 for i in range(1, 5))
    amplitudes = s.symbols("a b c d")
    norm = s.expand((s.Matrix(amplitudes).T * norm_gram * s.Matrix(amplitudes))[0] / 120**2)
    assert norm == (
        380 * amplitudes[0]**2 + 40 * amplitudes[0] * amplitudes[1]
        + 8 * amplitudes[1]**2 + 33 * amplitudes[2]**2 + amplitudes[3]**2 / 120
    )
    print("Norm Gram of 120*T:", norm_gram.tolist())
    print("Norm squared of T:", norm)

    triples = tuple(product(range(5), repeat=3))
    vectors = [
        {(a, b): tuple(tensor[(a, b) + tail] for tail in triples)
         for a, b in product(range(5), repeat=2)}
        for tensor in tensors
    ]
    monomials = tuple(combinations_with_replacement(range(4), 2))
    norm_coefficients = [int(norm_gram[i, j]) * (1 if i == j else 2) for i, j in monomials]
    rows = {}
    trace_coefficients = [0] * len(monomials)
    projector_trace = 0
    for a, b, c, d in product(range(5), repeat=4):
        projector10 = (
            5 * (int(a == c) * int(b == d) + int(a == d) * int(b == c))
            - 2 * int(a == b) * int(c == d)
        )
        row = []
        for index, (norm_coefficient, (i, j)) in enumerate(zip(norm_coefficients, monomials)):
            coefficient = sum(v * w for v, w in zip(vectors[i][a, b], vectors[j][c, d]))
            if i != j:
                coefficient += sum(v * w for v, w in zip(vectors[j][a, b], vectors[i][c, d]))
            row.append(140 * coefficient - norm_coefficient * projector10)
            if (c, d) == (a, b):
                trace_coefficients[index] += coefficient
        if (c, d) == (a, b):
            projector_trace += projector10
        if any(row):
            rows.setdefault(primitive([s.Rational(value) for value in row]), (a, b, c, d))
    assert tuple(rows) == GRAM_ROWS and s.Matrix(GRAM_ROWS).rank() == 4
    assert trace_coefficients == norm_coefficients and projector_trace == 140
    quadrics = [
        s.expand(sum(c * amplitudes[i] * amplitudes[j] for c, (i, j) in zip(row, monomials)))
        for row in GRAM_ROWS
    ]
    print("All 625 Gram-deviation slots are equivalent to four independent quadrics.")
    print("The Sym0 projector has trace 14; every amplitude coefficient of tr(K)=||T||^2 agrees.")
    print("Amplitude monomial order:", monomials)
    print("Quadratic coefficient rows:", GRAM_ROWS)
    for index, quadric in enumerate(quadrics):
        print(f"f{index} = {quadric}")

    # The highest nonzero amplitude is set to one. These four charts cover
    # every nonzero vector; no division by an unknown polynomial is used.
    for chart, degree in ((3, 3), (2, 1), (1, 0), (0, 0)):
        substitutions = {amplitudes[chart]: 1}
        substitutions.update({amplitudes[j]: 0 for j in range(chart + 1, 4)})
        equations = [s.expand(value.subs(substitutions)) for value in quadrics]
        multipliers = rational_identity(equations, amplitudes[:chart], degree)
        term_counts = tuple(
            0 if value == 0 else len(s.Poly(value, *amplitudes[:chart]).terms()) if chart else 1
            for value in multipliers
        )
        print(f"Projective chart {amplitudes[chart]}=1, later amplitudes zero: multiplier degree <= {degree}.")
        print("Multiplier term counts:", term_counts, "; expanded identity sum(h_i*f_i)=1 verified.")
    print("No nonzero coordinate-S5-invariant real harmonic quintic satisfies the scalar Gram criterion.")


def standard_s6():
    x = s.symbols("x0:6")
    y = tuple(value - sum(x) / 6 for value in x)
    radius = sum(value**2 for value in y)
    p3 = sum(value**3 for value in y)
    p5 = sum(value**5 for value in y)
    reductions = (
        (1, 0), (-1, 0), (-1, 1), (1, -s.Rational(1, 2)),
        (1, -1), (-1, s.Rational(5, 6)), (s.Rational(1, 5), -s.Rational(1, 6)),
    )
    for partition, (coefficient5, coefficient23) in zip(PARTITIONS, reductions):
        orbit = sum(
            prod(value**power for value, power in zip(y, powers))
            for powers in sorted(set(permutations(partition + (0,))))
        )
        assert s.expand(orbit - coefficient5 * p5 - coefficient23 * radius * p3) == 0
    assert s.expand(sum(s.diff(p5, value, 2) for value in x) - s.Rational(50, 3) * p3) == 0
    assert s.expand(sum(s.diff(radius * p3, value, 2) for value in x) - 22 * p3) == 0
    expression = s.expand(p5 - s.Rational(25, 33) * radius * p3)
    polynomial = s.Poly(expression, *x, domain=s.QQ)
    assert expression != 0 and polynomial.total_degree() == 5
    assert s.expand(sum(s.diff(expression, value) for value in x)) == 0
    assert s.expand(sum(s.diff(expression, value, 2) for value in x)) == 0
    assert s.expand(expression - expression.xreplace(dict(zip(x, x[1:] + x[:1])))) == 0
    assert s.expand(expression - expression.xreplace({x[0]: x[1], x[1]: x[0]})) == 0
    assert not expression.atoms(s.Float)
    print("Standard S6 on E={sum(y_i)=0}: invariant degree-five span p5, p2*p3; harmonic dimension 1.")
    print("All seven monomial-symmetric orbit sums reduce exactly to the stated two-dimensional span.")
    print("Delta(p5)=(50/3)*p3 and Delta(p2*p3)=22*p3; generator F=p5-(25/33)*p2*p3.")

    components = {
        powers: Fraction(value) * prod(factorial(k) for k in powers) / factorial(5)
        for powers, value in polynomial.terms()
    }
    denominator = lcm(*(value.denominator for value in components.values()))
    tensor = {}
    for index in product(range(6), repeat=5):
        value = components.get(tuple(index.count(i) for i in range(6)), Fraction(0)) * denominator
        assert value.denominator == 1
        tensor[index] = value.numerator
    assert len(tensor) == 7776 and denominator == 198
    assert all(
        sum(tensor[(i,) + rest] for i in range(6)) == 0
        for rest in product(range(6), repeat=4)
    )
    assert all(
        sum(tensor[(i, i) + rest] for i in range(6)) == 0
        for rest in product(range(6), repeat=3)
    )
    norm_integer = sum(value**2 for value in tensor.values())
    norm = Fraction(norm_integer, denominator**2)
    assert norm == Fraction(10, 33)
    projector = tuple(
        tuple(Fraction(int(i == j)) - Fraction(1, 6) for j in range(6)) for i in range(6)
    )
    assert all(sum(row) == 0 for row in projector)
    assert all(
        sum(projector[i][k] * projector[k][j] for k in range(6)) == projector[i][j]
        for i, j in product(range(6), repeat=2)
    )
    triples = tuple(product(range(6), repeat=3))
    vectors = {
        (a, b): tuple(tensor[(a, b) + tail] for tail in triples)
        for a, b in product(range(6), repeat=2)
    }
    P, K = {}, {}
    for a, b, c, d in product(range(6), repeat=4):
        P[a, b, c, d] = (
            (projector[a][c] * projector[b][d] + projector[a][d] * projector[b][c]) / 2
            - projector[a][b] * projector[c][d] / 5
        )
        K[a, b, c, d] = sum(v * w for v, w in zip(vectors[a, b], vectors[c, d]))
    assert len(P) == len(K) == 1296
    assert sum(P[a, b, a, b] for a, b in product(range(6), repeat=2)) == 14
    assert sum(K[a, b, a, b] for a, b in product(range(6), repeat=2)) == norm_integer
    scale = norm / 14
    differences = {
        key: Fraction(value, denominator**2) - scale * P[key] for key, value in K.items()
    }
    assert all(differences.values())
    key = (0, 0, 0, 0)
    actual = Fraction(K[key], denominator**2)
    expected = scale * P[key]
    assert actual == Fraction(35, 3267)
    assert expected == Fraction(25, 2079) and differences[key] == Fraction(-10, 7623)
    print("Symmetric tensor: 7776 ordered components; 198*T is integral.")
    print("All-ones annihilation, trace-free contractions, translation and permutation invariance verified.")
    print("Norm squared:", norm, "; required scalar ||T||^2/14:", scale)
    print("The projected Sym0(E) identity has trace 14; tr(K)=||T||^2 verified.")
    print("All 1296 comparison slots are nonzero for the generator.")
    print("Slot (0,0,0,0): K =", actual, "; scalar*projector =", expected, "; difference =", differences[key])
    print("No nonzero standard-S6-invariant real harmonic quintic satisfies the scalar Gram criterion.")


def main():
    started = perf_counter()
    coordinate_s5()
    standard_s6()
    print("Scope: these two complete invariant families; no classification of all harmonic quintics or field solutions.")
    print(f"Runtime {perf_counter() - started:.3f} seconds.")


if __name__ == "__main__":
    main()
