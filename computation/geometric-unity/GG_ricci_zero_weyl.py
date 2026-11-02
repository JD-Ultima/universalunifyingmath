"""Fresh Ricci-zero certificates with the flat origin kept separate.

The required union contains 15 node systems at one positive fibre: ten bivariate
compatibility nodes, alpha=+/-1 rank nodes and the second pure-mixed sign. Every
system is assembled once in memory. No nonzero-Ricci rank data or scratch caches
are reused. Prior jet bounds are coefficient total degree<=2 and RHS<=3.
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
import GG_rank_blocks as Rank
from GG_chiral_action import verify_chiral_action

A, B = sp.symbols("alpha beta", real=True)
LAYOUT = (2, 3, 5, 7)
ROWS = (168, 337, 492)
COMPAT_NODES = tuple((i, total - i) for total in range(4) for i in range(total + 1))
ALL_NODES = tuple(dict.fromkeys(COMPAT_NODES + ((1, -1), (-1, -1), (-1, 0), (-1, 1), (0, -1))))
MON3 = tuple(A**i * B**(total - i) for total in range(4) for i in range(total + 1))
MON2 = MON3[:6]
V2 = sp.Matrix([[mon.subs({A: i, B: j}) for mon in MON2] for i, j in COMPAT_NODES[:6]])
V3 = sp.Matrix([[mon.subs({A: i, B: j}) for mon in MON3] for i, j in COMPAT_NODES])
INV2, INV3 = V2.inv(), V3.inv()
CHECKS = []


def check(label, condition):
    ok = bool(condition)
    CHECKS.append(ok)
    print(f"[{len(CHECKS):02d}] {label}: {ok}", flush=True)
    assert ok, label


def tensors():
    def theta(i, j):
        return [[Fr(int((a, b) == (i, j)) - int((a, b) == (j, i))) for b in range(4)] for a in range(4)]

    weights = {(0, 1): 1, (2, 3): 1, (0, 2): -1, (1, 3): -1}
    forms = {pair: theta(*pair) for pair in weights}
    keys = tuple(itertools.product(range(4), repeat=4))
    E0 = {key: sum(weights[p] * forms[p][key[0]][key[1]] * forms[p][key[2]][key[3]] for p in forms) for key in keys}
    H = {key: sum(forms[p][key[0]][key[1]] * forms[q][key[2]][key[3]]
                  + forms[q][key[0]][key[1]] * forms[p][key[2]][key[3]]
                  for p, q in (((0, 1), (2, 3)), ((0, 2), (1, 3)))) for key in keys}
    return E0, H


E0, H = tensors()


def system(alpha, beta):
    R = {key: Fr(alpha) * E0[key] + Fr(beta) * H[key] for key in E0}
    base = U.X.weyl_base(R, (1, 1, 1, 1))
    data = U.X.point_data(base, U.diagonal(LAYOUT), Fr(1, 5), Fr(2), Fr(3), riem=True)
    return U.relaxed(data)


def interpolate(values, degree):
    if all(value == values[0] for value in values):
        return sp.Rational(values[0])
    monomials, inverse = (MON2, INV2) if degree == 2 else (MON3, INV3)
    coefficients = inverse * sp.Matrix([sp.Rational(value) for value in values])
    return sp.factor(sum(coefficient * mon for coefficient, mon in zip(coefficients, monomials)))


def compatibility(systems):
    nodes = [systems[node] for node in COMPAT_NODES]
    coeff = [[interpolate([M[row][col] for M, _ in nodes[:6]], 2) for col in range(42)] for row in ROWS]
    rhs = [interpolate([b[row] for _, b in nodes], 3) for row in ROWS]
    allowed = {27, 32, 36, 39, 41}
    check("the selected rows have only xi_13 and a common curvature-independent n term",
          all(poly == 0 for row in coeff for col, poly in enumerate(row) if col not in allowed)
          and all(row[col] == coeff[0][col] and not row[col].has(A, B) for row in coeff for col in (32, 36, 39, 41)))
    check("all unused coefficient nodes check the reconstructed selected rows",
          all(poly.subs({A: i, B: j}) == M[row][col] for (i, j), (M, _) in systems.items()
              if (i, j) not in COMPAT_NODES[:6] for row, polynomials in zip(ROWS, coeff) for col, poly in enumerate(polynomials)))
    check("the five extra normalized nodes independently check the RHS interpolation",
          all(poly.subs({A: i, B: j}) == b[row] for (i, j), (_, b) in systems.items()
              if (i, j) not in COMPAT_NODES for row, poly in zip(ROWS, rhs)))
    expected_q = [-(371 * A**2 + 300 * B**2 + 1750) / 4200,
                  -(565 * A**2 + 546 * B**2 + 4900) / 11760,
                  -(4113 * A**2 + 3570 * B**2 + 24500) / 58800]
    expected_b = [8 * B**2 / 245, (16 * A**2 + 21 * B**2) / 2744,
                  3 * (8 * A**2 + 91 * B**2) / 13720]
    check("the exact selected coefficient and RHS formulas match",
          all(sp.expand(row[27] - expected) == 0 for row, expected in zip(coeff, expected_q))
          and all(sp.expand(poly - expected) == 0 for poly, expected in zip(rhs, expected_b)))
    determinant = sp.factor(sp.det(sp.Matrix([[coeff[i][27], 1, rhs[i]] for i in range(3)])))
    sos = (529 * A**4 + 1052 * A**2 * B**2 + 735 * B**4) / 14406000
    check("the compatibility determinant is the exact positive SOS polynomial", sp.expand(determinant - sos) == 0)
    check("the SOS is strictly positive away from the REAL origin and vanishes there",
          all(coefficient > 0 for coefficient in (529, 1052, 735, 14406000))
          and determinant.subs({A: 0, B: 0}) == 0)
    print("  xi coefficients:", [row[27] for row in coeff], flush=True)
    print("  RHS:", rhs, flush=True)
    print("  exact compatibility determinant:", determinant, flush=True)


def main():
    start = time.time()
    print("GG. Ricci zero with nonzero Weyl curvature and singular chiral blocks", flush=True)
    print("Fresh portable path: the required union of 15 node systems once in memory; no external caches.", flush=True)
    check("both fixed tensors have zero Ricci and satisfy Bianchi",
          all(sum(W[(i, j, i, k)] for i in range(4)) == 0 for W in (E0, H) for j, k in itertools.product(range(4), repeat=2))
          and all(W[(i, j, k, l)] + W[(i, k, l, j)] + W[(i, l, j, k)] == 0
                  for W in (E0, H) for i, j, k, l in itertools.product(range(4), repeat=4)))
    chiral_checks = verify_chiral_action(E0, H, emit=lambda text: print(text, flush=True))
    check("compatibility and normalized-rank node sets agree with their proved interpolation bounds",
          abs(V2.det()) == 4 and abs(V3.det()) == 576 and tuple(Rank.NODES) == (-1, 0, 1)
          and Rank.LAYOUT == LAYOUT and len(ALL_NODES) == 15)
    systems = {node: system(*node) for node in ALL_NODES}
    compatibility(systems)
    rank_checks = 0
    for sign in (1, -1):
        matrices = [systems[(sign, t)][0] for t in Rank.NODES]
        rank_checks += Rank.verify_normalized_case(sign, matrices, emit=lambda text: print(text, flush=True))
    for sign in (1, -1):
        M, b = systems[(0, sign)]
        rank_checks += Rank.verify_pure_case(sign, M, b, emit=lambda text: print(text, flush=True))
    flat, flat_b = systems[(0, 0)]
    check("the flat origin retains coefficient/augmented ranks 28/28",
          Rank.rational_matrix(flat).rank() == 28 and not any(flat_b)
          and Rank.rational_matrix([row + [rhs] for row, rhs in zip(flat, flat_b)]).rank() == 28)
    check("all 15 exact node systems were assembled once", len(systems) == len(ALL_NODES) == 15)
    print("The Ricci-independent covariance proof transports each normalized rank42 and incompatible witness TOGETHER to one positive fibre of the original metric.", flush=True)
    print("No claim that the original fixed numeric fibre works at every amplitude is used. The flat origin is excluded from the nonzero-Weyl theorem.", flush=True)
    print("Higher-jet and local scope use the fixed-curvature-independent support argument and nonzero Riemann curvature, without assuming Ricci remains zero nearby.", flush=True)
    print(f"\n{len(CHECKS)} main + {rank_checks} rank + {chiral_checks} chiral checks, all True; runtime {time.time() - start:.1f} s", flush=True)


if __name__ == "__main__":
    main()
