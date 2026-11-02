"""Exact diagonal-Weyl compatibility and finite same-witness rank certificates.

Default reproduction freshly assembles all 60 exact node systems once in memory.
No off-repository path, cache or fitted determinant degree is used. BB1's hand
jet argument supplies total degrees A<=2 and RHS<=3 before interpolation.
CC_rank_blocks verifies compact rational row-combination certificates against
the same freshly computed matrices. Higher-jet scope rests on BB3's accepted
curvature-independent support/parity/gradient proof, not node sampling.
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
import CC_rank_blocks as Rank

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


def tensors():
    g = U.X.gmat(ETA)
    P = [[(Fr(int(i == j == 0)) - g[i][j] / 6) / 2 for j in range(4)] for i in range(4)]
    R0 = U.X.kn(P, g)

    def weyl(weights):
        forms = {pair: [[Fr(int((i, j) == pair) - int((j, i) == pair)) for j in range(4)]
                        for i in range(4)] for pair in weights}
        return {key: sum(Fr(weights[pair]) * forms[pair][key[0]][key[1]]
                         * forms[pair][key[2]][key[3]] for pair in forms) for key in R0}

    WA = weyl({(0, 1): 1, (2, 3): 1, (0, 3): -1, (1, 2): -1})
    WB = weyl({(0, 2): 1, (1, 3): 1, (0, 3): -1, (1, 2): -1})
    return R0, WA, WB


R0, WA, WB = tensors()


def curvature_checks():
    check("both Weyl basis tensors have zero Ricci contraction and satisfy Bianchi",
          all(sum(W[(i, j, i, k)] for i in range(4)) == 0
              for W in (WA, WB) for j, k in itertools.product(range(4), repeat=2))
          and all(W[(i, j, k, l)] + W[(i, k, l, j)] + W[(i, l, j, k)] == 0
                  for W in (WA, WB) for i, j, k, l in itertools.product(range(4), repeat=4)))
    permutation, signs = (0, 1, 3, 2), (-1, 1, 1, 1)
    S = sp.zeros(4)
    for i in range(4):
        S[permutation[i], i] = signs[i]

    def transformed(T, key):
        return sp.prod(signs[i] for i in key) * T[tuple(permutation[i] for i in key)]

    check("the signed perpendicular swap is in SO(4) and acts as (a,b)->(a,-a-b)",
          S.det() == 1 and S.T * S == sp.eye(4)
          and all(transformed(R0, key) == R0[key]
                  and transformed(WA, key) == WA[key] - WB[key]
                  and transformed(WB, key) == -WB[key] for key in R0))


def system(rho, node, layout):
    av, bv = node
    curvature = {key: Fr(rho) * R0[key] + Fr(av) * WA[key] + Fr(bv) * WB[key] for key in R0}
    base = U.X.weyl_base(curvature, ETA)
    data = U.X.point_data(base, U.diagonal(layout), Fr(1, 5), Fr(2), Fr(3), riem=True)
    return U.relaxed(data)


def interpolate(values, degree):
    monomials, inverse = (MON2, INV2) if degree == 2 else (MON3, INV3)
    coefficients = inverse * sp.Matrix([sp.Rational(value) for value in values])
    return sp.factor(sum(coefficient * mon for coefficient, mon in zip(coefficients, monomials)))


def compatibility_case(rho, layout, systems):
    prefix = f"rho={rho}, y={layout}"
    allowed = {27, 32, 36, 39, 41}
    check(prefix + ": selected rows have only xi_13 and diagonal n coefficients",
          all(matrix[row][col] == 0 for matrix, _ in systems[:6] for row in ROWS
              for col in range(42) if col not in allowed))
    check(prefix + ": the n combination is common to the three rows and curvature independent",
          all(matrix[row][col] == systems[0][0][ROWS[0]][col]
              for matrix, _ in systems for row in ROWS for col in (32, 36, 39, 41)))
    xi = [interpolate([matrix[row][27] for matrix, _ in systems[:6]], 2) for row in ROWS]
    rhs = [interpolate([vector[row] for _, vector in systems], 3) for row in ROWS]
    check(prefix + ": four additional nodes check every coefficient polynomial",
          all(interpolate([matrix[row][col] for matrix, _ in systems[:6]], 2).subs({A: i, B: j})
              == matrix[row][col]
              for (i, j), (matrix, _) in zip(NODES[6:], systems[6:])
              for row in ROWS for col in range(42)))
    scales = {
        LAYOUTS[0]: (Fr(1, 500), Fr(1, 4500), Fr(1, 12000)),
        LAYOUTS[1]: (Fr(5, 4116), Fr(1, 6174), Fr(1, 20580)),
        LAYOUTS[2]: (Fr(7, 8748), Fr(5, 69984), Fr(1, 21870)),
    }[layout]
    expected = [sp.Rational(scales[0]) * (3 * A + 3 * B - rho) ** 2,
                sp.Rational(scales[1]) * (6 * B - rho) ** 2,
                sp.Rational(scales[2]) * (6 * A - rho) ** 2]
    check(prefix + ": exact RHS polynomials have the stated quadratic square formulas",
          all(sp.expand(x - y) == 0 for x, y in zip(rhs, expected)))
    determinant = sp.factor(sp.det(sp.Matrix([[xi[i], 1, rhs[i]] for i in range(3)])))
    numerator = sp.Poly(sp.together(determinant).as_numer_denom()[0], A, B)
    _, primitive = numerator.primitive()
    print("  xi coefficients:", xi, flush=True)
    print("  RHS:", rhs, flush=True)
    print("  primitive compatibility numerator:", primitive.as_expr(), flush=True)
    return primitive.as_expr()


def ideal_checks(rho, polynomials):
    ideal = sp.groebner(polynomials, A, B, order="grevlex", domain=sp.QQ)
    expected = sp.groebner([(6 * A - rho) ** 2, (6 * A - rho) * (6 * B - rho),
                           (6 * B - rho) ** 2], A, B, order="grevlex", domain=sp.QQ)
    check(f"rho={rho}: the three compatibility polynomials generate the squared point ideal", ideal == expected)
    print("  exact compatibility Groebner basis:", list(ideal), flush=True)
    rotated = sp.expand(polynomials[0].subs(B, -A - B))
    point = {A: sp.Rational(rho, 6), B: sp.Rational(rho, 6)}
    value = rotated.subs(point)
    unit = sp.groebner(polynomials + [rotated], A, B, order="grevlex", domain=sp.QQ)
    check(f"rho={rho}: adding the rotated first witness gives the unit ideal", value != 0 and list(unit) == [sp.S.One])
    print("  rotated primitive numerator at the sole original exception:", value, flush=True)


def main():
    start = time.time()
    print("CC. Diagonal Weyl curvature with arbitrary higher metric jets", flush=True)
    print("Fresh exact reproduction: all 60 node systems are assembled once in memory; no external caches.", flush=True)
    curvature_checks()
    check("total-degree node sets are unisolvent and agree with the rank verifier",
          abs(V2.det()) == 4 and abs(V3.det()) == 576 and tuple(Rank.DEG2_NODES) == NODES[:6])
    print("  degree-three nodes:", NODES, flush=True)
    assembled = 0
    for rho in (1, -1):
        polynomials = []
        for layout in LAYOUTS:
            case_start = time.time()
            systems = [system(rho, node, layout) for node in NODES]
            assembled += len(systems)
            polynomials.append(compatibility_case(rho, layout, systems))
            summary = Rank.verify_rank_case(rho, layout, [matrix for matrix, _ in systems[:6]],
                                            emit=lambda text: print(text, flush=True))
            print("  rank certificate summary:", summary, flush=True)
            print(f"  fresh case runtime {time.time() - case_start:.1f} s", flush=True)
        ideal_checks(rho, polynomials)
    check("all 60 exact node systems were freshly assembled once", assembled == 60)
    print("The three numerical fibres and the signed-permutation fourth configuration provide finite same-witness rank and incompatibility.", flush=True)
    print("Higher-jet scope uses the proved cubic-invariant subsystem, unchanged diagonal-vertical certificate rows and completion restrictions; no reflection symmetry of arbitrary higher jets is assumed.", flush=True)
    print("This covers the declared simultaneous Ricci-axis/diagonal-Weyl frame orbit, not arbitrary Weyl tensors or Ricci-flat nonflat points.", flush=True)
    print(f"\n{len(CHECKS)} main checks, all True; rank-helper assertions also passed; runtime {time.time() - start:.1f} s", flush=True)


if __name__ == "__main__":
    main()
