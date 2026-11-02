"""Fresh exact reproduction for nonzero rank-one Ricci and zero electric Weyl.

The default assembles 60 exact node systems once in memory and verifies compact
constant row-combination certificates against their first six coefficient
matrices per case. Total degree A<=2 and RHS<=3 is established by the prior hand
jet argument before interpolation. Completion nonvanishing is proved for REAL
parameters using positive constants plus positive semidefinite quadratic forms.
No off-repository paths, stored node matrices or scratch caches are used.
"""
import itertools
import os
import pathlib
import sys
import time
from fractions import Fraction as Fr

for _key in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_key, "1")

import sympy as sp

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import U_rank_one_frontier as U
import EE_rank_blocks as Rank

A, B = sp.symbols("a b")
ETA = (1, 1, 1, 1)
ROWS = (168, 337, 492)
LAYOUTS = ((2, 3, 4, 5), (2, 3, 5, 7), (2, 4, 5, 9))
NODES = tuple((i, total - i) for total in range(4) for i in range(total + 1))
MON3 = tuple(A**i * B**(total - i) for total in range(4) for i in range(total + 1))
MON2 = MON3[:6]
V2 = sp.Matrix([[mon.subs({A: i, B: j}) for mon in MON2] for i, j in NODES[:6]])
V3 = sp.Matrix([[mon.subs({A: i, B: j}) for mon in MON3] for i, j in NODES])
INV2, INV3 = V2.inv(), V3.inv()
CHECKS = []


def check(label, condition):
    ok = bool(condition)
    CHECKS.append(ok)
    print(f"[{len(CHECKS):02d}] {label}: {ok}", flush=True)
    assert ok, label


def theta(i, j):
    return [[Fr(int((a, b) == (i, j)) - int((a, b) == (j, i))) for b in range(4)] for a in range(4)]


def tensors():
    g = U.X.gmat(ETA)
    P = [[(Fr(int(i == j == 0)) - g[i][j] / 6) / 2 for j in range(4)] for i in range(4)]
    R0 = U.X.kn(P, g)

    def cross(p, q, key):
        x, y = theta(*p), theta(*q)
        return x[key[0]][key[1]] * y[key[2]][key[3]] + y[key[0]][key[1]] * x[key[2]][key[3]]

    HA = {key: cross((0, 1), (2, 3), key) - cross((0, 3), (1, 2), key) for key in R0}
    HB = {key: -cross((0, 2), (1, 3), key) - cross((0, 3), (1, 2), key) for key in R0}
    return R0, HA, HB


R0, HA, HB = tensors()


def normal_form_checks():
    check("both pure-mixed basis tensors satisfy Bianchi, zero Ricci and zero electric Weyl",
          (HA[(0, 1, 2, 3)], HA[(0, 2, 3, 1)], HA[(0, 3, 1, 2)]) == (1, 0, -1)
          and (HB[(0, 1, 2, 3)], HB[(0, 2, 3, 1)], HB[(0, 3, 1, 2)]) == (0, 1, -1)
          and all(W[(i, j, k, l)] + W[(i, k, l, j)] + W[(i, l, j, k)] == 0
                  for W in (HA, HB) for i, j, k, l in itertools.product(range(4), repeat=4))
          and all(sum(W[(i, j, i, k)] for i in range(4)) == 0
                  for W in (HA, HB) for j, k in itertools.product(range(4), repeat=2))
          and all(W[(0, i, 0, j)] == 0 for W in (HA, HB) for i, j in itertools.product(range(4), repeat=2)))

    dual = (theta(2, 3), theta(3, 1), theta(1, 2))
    electric = (theta(0, 1), theta(0, 2), theta(0, 3))
    variables = sp.symbols("b11 b12 b13 b21 b22 b23 b31 b32 b33")
    magnetic = sp.Matrix(3, 3, variables)
    W = {key: sp.expand(sum(magnetic[i, j] * (electric[i][key[0]][key[1]] * dual[j][key[2]][key[3]]
                                             + dual[j][key[0]][key[1]] * electric[i][key[2]][key[3]])
                                for i, j in itertools.product(range(3), repeat=2))) for key in R0}
    Ric = sp.Matrix(4, 4, lambda i, j: sum(W[(k, i, k, j)] for k in range(4)))
    expected = sp.zeros(4)
    for j, value in enumerate((magnetic[1, 2] - magnetic[2, 1], magnetic[2, 0] - magnetic[0, 2],
                               magnetic[0, 1] - magnetic[1, 0]), 1):
        expected[0, j] = expected[j, 0] = value
    recovered = sp.Matrix(3, 3, lambda i, j: sum(dual[j][k][l] * W[(0, i + 1, k, l)]
                                                for k, l in itertools.product(range(1, 4), repeat=2)) / 2)
    check("magnetic extraction is exact; Weyl Ricci gives symmetry and Bianchi gives trace zero",
          Ric == expected and recovered == magnetic
          and W[(0, 1, 2, 3)] + W[(0, 2, 3, 1)] + W[(0, 3, 1, 2)] == sp.trace(magnetic))

    q = sp.symbols("q0:6")
    spatial = (theta(1, 2), theta(1, 3), theta(2, 3))
    C, position = sp.zeros(3), 0
    for i in range(3):
        for j in range(i, 3):
            C[i, j] = C[j, i] = q[position]
            position += 1
    S = {key: sp.expand(sum(C[i, j] * spatial[i][key[0]][key[1]] * spatial[j][key[2]][key[3]]
                            for i, j in itertools.product(range(3), repeat=2)))
         for key in itertools.product(range(1, 4), repeat=4)}
    Ric3 = sp.Matrix(3, 3, lambda i, j: sum(S[(k, i + 1, k, j + 1)] for k in range(1, 4)))
    outputs = [Ric3[i, j] for i in range(3) for j in range(i, 3)]
    contraction = sp.Matrix([[sp.diff(value, qi) for qi in q] for value in outputs])
    check("the spatial three-dimensional curvature-to-Ricci map is invertible", contraction.det() == -2)
    check("the entire pure-mixed plane has precisely the baseline even-sign symmetry",
          all(sp.prod(signs[i] for i in key) * W[key] == sp.prod(signs) * W[key]
              for W in (HA, HB) for signs in itertools.product((-1, 1), repeat=4) for key in W))


def system(rho, node, layout):
    av, bv = node
    curvature = {key: Fr(rho) * R0[key] + Fr(av) * HA[key] + Fr(bv) * HB[key] for key in R0}
    base = U.X.weyl_base(curvature, ETA)
    D = U.X.point_data(base, U.diagonal(layout), Fr(1, 5), Fr(2), Fr(3), riem=True)
    return U.relaxed(D)


def interpolate(values, degree):
    if all(value == values[0] for value in values):
        return sp.Rational(values[0])
    monomials, inverse = (MON2, INV2) if degree == 2 else (MON3, INV3)
    coefficients = inverse * sp.Matrix([sp.Rational(value) for value in values])
    return sp.factor(sum(coefficient * mon for coefficient, mon in zip(coefficients, monomials)))


def analyse(rho, layout, systems):
    prefix = f"rho={rho}, y={layout}"
    coeff = [[interpolate([M[row][col] for M, _ in systems[:6]], 2) for col in range(42)] for row in ROWS]
    rhs = [interpolate([vector[row] for _, vector in systems], 3) for row in ROWS]
    allowed = {27, 32, 36, 39, 41}
    check(prefix + ": selected rows have only xi_13 and diagonal n coefficients",
          all(poly == 0 for row in coeff for col, poly in enumerate(row) if col not in allowed))
    check(prefix + ": n coefficients are common and independent of curvature",
          all(row[col] == coeff[0][col] and not row[col].has(A, B)
              for row in coeff for col in (32, 36, 39, 41)))
    check(prefix + ": four additional nodes check every selected coefficient polynomial",
          all(poly.subs({A: i, B: j}) == M[row][col]
              for (i, j), (M, _) in zip(NODES[6:], systems[6:])
              for row, polynomials in zip(ROWS, coeff) for col, poly in enumerate(polynomials)))
    check(prefix + ": the exact RHS interpolation has zero cubic coefficients",
          all(sp.Poly(poly, A, B).total_degree() <= 2 for poly in rhs))
    determinant = sp.factor(sp.det(sp.Matrix([[coeff[i][27], 1, rhs[i]] for i in range(3)])))
    numerator = sp.Poly(sp.together(determinant).as_numer_denom()[0], A, B)
    _, numerator = numerator.primitive()
    print("  xi coefficients:", [row[27] for row in coeff], flush=True)
    print("  RHS:", rhs, flush=True)
    print("  primitive compatibility numerator:", numerator.as_expr(), flush=True)
    return numerator.as_expr()


def main():
    start = time.time()
    print("EE. Nonzero rank-one Ricci tensor and zero electric Weyl part", flush=True)
    print("Fresh portable path: all 60 node systems once in memory; no external caches.", flush=True)
    normal_form_checks()
    check("total-degree node sets are unisolvent and agree with the rank helper",
          abs(V2.det()) == 4 and abs(V3.det()) == 576 and tuple(Rank.NODES) == NODES[:6])
    assembled = 0
    for rho in (1, -1):
        polynomials = []
        for layout in LAYOUTS:
            systems = [system(rho, node, layout) for node in NODES]
            assembled += len(systems)
            polynomials.append(analyse(rho, layout, systems))
            summary = Rank.verify_rank_case(rho, layout, [M for M, _ in systems[:6]],
                                            emit=lambda text: print(text, flush=True))
            print("  rank certificate summary:", summary, flush=True)
        G = sp.groebner(polynomials, A, B, order="grevlex", domain=sp.QQ)
        check(f"rho={rho}: the three compatibility numerators generate the unit ideal", list(G) == [sp.S.One])
        print("  exact compatibility Groebner basis:", list(G), flush=True)
    check("all 60 exact node systems were assembled once", assembled == 60)
    print("The three finite compatibility fibres have full coefficient rank at the same configurations for all REAL parameters.", flush=True)
    print("Real completion positivity is separate from the complex compatibility unit ideals. Higher-jet scope uses the fixed-curvature-independent DD3 proof.", flush=True)
    print("The intrinsic result covers zero electric Weyl relative to the Ricci eigenline; simultaneous nonzero electric and mixed Weyl and Ricci-flat nonflat points remain outside it.", flush=True)
    print(f"\n{len(CHECKS)} main checks, all True; rank-helper checks also passed; runtime {time.time() - start:.1f} s", flush=True)


if __name__ == "__main__":
    main()
