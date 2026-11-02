"""Fresh exact reproduction for the TWO fixed joint Weyl directions FF.

The default assembles 60 exact node systems once in memory. Prior hand jet
bounds give total degree A<=2 and RHS<=3 before interpolation. Constant rational
block certificates and J7 completion are checked against these fresh matrices.
No scratch paths, caches, parameter searches or sampled-rank inference are used.
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
import FF_rank_blocks as Rank

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
    weights = {(0, 1): 1, (2, 3): 1, (0, 2): -1, (1, 3): -1}
    forms = {pair: theta(*pair) for pair in weights}
    E0 = {key: sum(weights[p] * forms[p][key[0]][key[1]] * forms[p][key[2]][key[3]] for p in forms) for key in R0}
    H = {key: sum(forms[p][key[0]][key[1]] * forms[q][key[2]][key[3]]
                  + forms[q][key[0]][key[1]] * forms[p][key[2]][key[3]]
                  for p, q in (((0, 1), (2, 3)), ((0, 2), (1, 3)))) for key in R0}
    return R0, E0, H


R0, E0, H = tensors()


def curvature_checks():
    check("the fixed electric and mixed tensors satisfy curvature symmetries, Bianchi and zero Ricci",
          all(W[(i, j, k, l)] == -W[(j, i, k, l)] == -W[(i, j, l, k)]
              and W[(i, j, k, l)] == W[(k, l, i, j)]
              and W[(i, j, k, l)] + W[(i, k, l, j)] + W[(i, l, j, k)] == 0
              for W in (E0, H) for i, j, k, l in itertools.product(range(4), repeat=4))
          and all(sum(W[(i, j, i, k)] for i in range(4)) == 0
                  for W in (E0, H) for j, k in itertools.product(range(4), repeat=2)))
    dual = (theta(2, 3), theta(3, 1), theta(1, 2))

    def electric(W):
        return sp.Matrix(3, 3, lambda i, j: W[(0, i + 1, 0, j + 1)])

    def magnetic(W):
        return sp.Matrix(3, 3, lambda i, j: sum(dual[j][k][l] * W[(0, i + 1, k, l)]
                                                for k, l in itertools.product(range(1, 4), repeat=2)) / 2)

    profile = sp.diag(1, -1, 0)
    e = A * electric(E0) + B * electric(H)
    m = A * magnetic(E0) + B * magnetic(H)
    kernel = sp.Matrix([0, 0, 1])
    check("electric/magnetic profiles have the exact simultaneous normal form and common kernel",
          electric(E0) == magnetic(H) == profile and electric(H) == magnetic(E0) == sp.zeros(3)
          and e * m == m * e and e * kernel == m * kernel == sp.zeros(3, 1))
    check("even sign flips fix the plane pointwise, while odd flips reverse only its mixed basis",
          all(sp.prod(signs[i] for i in key) * E0[key] == E0[key]
              and sp.prod(signs[i] for i in key) * H[key] == sp.prod(signs) * H[key]
              for signs in itertools.product((-1, 1), repeat=4) for key in R0))


def system(rho, node, layout):
    av, bv = node
    curvature = {key: Fr(rho) * R0[key] + Fr(av) * E0[key] + Fr(bv) * H[key] for key in R0}
    D = U.X.point_data(U.X.weyl_base(curvature, ETA), U.diagonal(layout), Fr(1, 5), Fr(2), Fr(3), riem=True)
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
    check(prefix + ": four additional nodes check the selected coefficient polynomials",
          all(poly.subs({A: i, B: j}) == M[row][col]
              for (i, j), (M, _) in zip(NODES[6:], systems[6:])
              for row, polynomials in zip(ROWS, coeff) for col, poly in enumerate(polynomials)))
    check(prefix + ": exact RHS interpolation has zero cubic coefficients",
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
    print("FF. A two-parameter joint electric/mixed Weyl plane", flush=True)
    print("Fresh portable path: all 60 exact node systems once in memory; no external caches.", flush=True)
    curvature_checks()
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
    print("Each finite compatibility configuration has constant coefficient rank42 at the SAME fibre for every parameter pair.", flush=True)
    print("Higher-jet transfer uses the fixed-curvature-independent DD3 argument; baseline even-sign symmetry is not imposed on actual higher jets.", flush=True)
    print("Intrinsic scope requires commuting electric/magnetic Weyl matrices AND a common nonzero kernel vector. General commuting pairs and general Weyl tensors remain outside it.", flush=True)
    print(f"\n{len(CHECKS)} main checks, all True; rank-helper checks also passed; runtime {time.time() - start:.1f} s", flush=True)


if __name__ == "__main__":
    main()
