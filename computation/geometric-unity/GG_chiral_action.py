"""Small exact chiral-frame controls for GG's pointwise geometric lemma.

This has no observer-metric engine. The orbit proof additionally uses the
spectral theorem and connected-group argument in the formalization.
"""
from itertools import combinations

import sympy as sp

PAIRS = tuple(combinations(range(4), 2))


def verify_chiral_action(electric_tensor, mixed_tensor, emit=print):
    count = 0

    def check(label, condition):
        nonlocal count
        count += 1
        assert bool(condition), label

    def wedge(u, v):
        return sp.Matrix([u[i] * v[j] - u[j] * v[i] for i, j in PAIRS])

    def e(index):
        return sp.eye(4)[:, index]

    def J(a, b):
        result = sp.zeros(4)
        result[b, a], result[a, b] = 1, -1
        return result

    def induced(generator):
        return sp.Matrix.hstack(*[wedge(generator * e(i), e(j)) + wedge(e(i), generator * e(j)) for i, j in PAIRS])

    theta = sp.Matrix.hstack(wedge(e(0), e(1)), wedge(e(0), e(2)), wedge(e(0), e(3)))
    tau = sp.Matrix.hstack(wedge(e(2), e(3)), wedge(e(3), e(1)), wedge(e(1), e(2)))
    P = sp.Matrix.hstack(theta + tau, theta - tau) / sp.sqrt(2)
    star = sp.Matrix([[sp.LeviCivita(i, j, k, l) for k, l in PAIRS] for i, j in PAIRS])
    check("orthonormal chiral bases and Hodge signs", P.T * P == sp.eye(6) and star * P == P * sp.diag(1, 1, 1, -1, -1, -1))

    generators, paired = [], []
    for i, (j, k) in enumerate(((2, 3), (3, 1), (1, 2)), start=1):
        cross = sp.Matrix([[sp.LeviCivita(i - 1, b, a) for b in range(3)] for a in range(3)])
        plus, minus = (J(0, i) + J(j, k)) / 2, (-J(0, i) + J(j, k)) / 2
        generators.extend((plus, minus))
        paired.extend((J(0, i), J(j, k)))
        check(f"independent chiral Lie actions for index {i}",
              sp.simplify(P.T * induced(plus) * P) == sp.diag(cross, sp.zeros(3))
              and sp.simplify(P.T * induced(minus) * P) == sp.diag(sp.zeros(3), cross))
    coefficients = []
    for generator in paired:
        value = sp.simplify(P.T * induced(generator) * P)
        coefficients.append([value[3 * sign + (i + 2) % 3, 3 * sign + (i + 1) % 3] for sign in (0, 1) for i in range(3)])
    action = sp.Matrix.hstack(*[sp.Matrix(column) for column in coefficients])
    check("the chiral differential is onto", action.det() == -8
          and sp.Matrix.hstack(*[generator.reshape(16, 1) for generator in generators]).rank() == 6)

    D = sp.diag(1, -1, 0)
    E0 = sp.Matrix([[electric_tensor[p + q] for q in PAIRS] for p in PAIRS])
    H = sp.Matrix([[mixed_tensor[p + q] for q in PAIRS] for p in PAIRS])
    check("the actual curvature tensors have the stated chiral profiles",
          sp.simplify(P.T * E0 * P) == sp.diag(D, D)
          and sp.simplify(P.T * H * P) == sp.diag(D, -D))
    emit(f"  chiral geometric controls: {count} exact checks; independent action determinant -8.")
    return count
