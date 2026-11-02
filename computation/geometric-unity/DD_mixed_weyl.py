"""Fresh exact reproduction for one genuinely mixed-Weyl pencil.

All 16 node systems are assembled once in memory. Prior hand jet bounds give
degree A<=2 and RHS<=3 before interpolation. Constant block left inverses and
J7 completion are checked against these same fresh coefficient matrices.
No scratch paths, caches or amplitude searches are used.
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
import DD_rank_blocks as Rank

W = sp.Symbol("w")
ETA = (1, 1, 1, 1)
NODES = (-1, 0, 1, 2)
LAYOUTS = ((2, 3, 4, 5), (2, 3, 5, 7))
ROWS = (168, 337, 492)
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
    pairs = ((0, 1), (2, 3), (0, 2), (1, 3))
    forms = {pair: [[Fr(int((i, j) == pair) - int((j, i) == pair)) for j in range(4)]
                    for i in range(4)] for pair in pairs}
    H = {key: sum(forms[p][key[0]][key[1]] * forms[q][key[2]][key[3]]
                  + forms[q][key[0]][key[1]] * forms[p][key[2]][key[3]]
                  for p, q in (((0, 1), (2, 3)), ((0, 2), (1, 3)))) for key in R0}
    return R0, H


R0, H = tensors()


def curvature_checks():
    check("H has the stated mixed entries, curvature symmetries, Bianchi and zero Ricci/electric part",
          (H[(0, 1, 2, 3)], H[(0, 2, 3, 1)], H[(0, 3, 1, 2)]) == (1, -1, 0)
          and all(H[(a, b, c, d)] == -H[(b, a, c, d)] == -H[(a, b, d, c)]
                  and H[(a, b, c, d)] == H[(c, d, a, b)]
                  and H[(a, b, c, d)] + H[(a, c, d, b)] + H[(a, d, b, c)] == 0
                  for a, b, c, d in itertools.product(range(4), repeat=4))
          and all(sum(H[(a, b, a, d)] for a in range(4)) == 0 for b, d in itertools.product(range(4), repeat=2))
          and all(H[(0, i, 0, j)] == 0 for i, j in itertools.product(range(4), repeat=2)))
    check("even sign flips fix H and odd sign flips reverse it",
          all(sp.prod(signs[i] for i in key) * H[key] == sp.prod(signs) * H[key]
              for signs in itertools.product((-1, 1), repeat=4) for key in H))


def system(rho, node, layout):
    R = {key: Fr(rho) * R0[key] + Fr(node) * H[key] for key in R0}
    base = U.X.weyl_base(R, ETA)
    D = U.X.point_data(base, U.diagonal(layout), Fr(1, 5), Fr(2), Fr(3), riem=True)
    return U.relaxed(D)


def polynomial(minus, zero, plus):
    return sp.factor(sp.Rational(zero) + (sp.Rational(plus) - minus) * W / 2
                     + ((sp.Rational(plus) + minus) / 2 - zero) * W**2)


def analyse(rho, layout, values):
    prefix = f"rho={rho}, y={layout}"
    Am, bm = values[0]
    A0, b0 = values[1]
    Ap, bp = values[2]
    A2, b2 = values[3]
    allowed = {27, 32, 36, 39, 41}
    check(prefix + ": selected rows have only xi_13 and the diagonal n combination",
          all(A[row][col] == 0 for A in (Am, A0, Ap) for row in ROWS for col in range(42) if col not in allowed))
    check(prefix + ": the n combination is common and curvature independent",
          all(A[row][col] == A0[ROWS[0]][col] for A, _ in values for row in ROWS for col in (32, 36, 39, 41)))
    check(prefix + ": fourth-node control checks the entire degree-two coefficient matrix",
          all(A2[row][col] == 3 * Ap[row][col] + Am[row][col] - 3 * A0[row][col]
              for row in range(1274) for col in range(42)))
    xi = [polynomial(Am[row][27], A0[row][27], Ap[row][27]) for row in ROWS]
    rhs = [sp.factor(sp.interpolate([(node, sp.Rational(vector[row])) for node, (_, vector) in zip(NODES, values)], W))
           for row in ROWS]
    if layout == LAYOUTS[0]:
        expected_rhs = [(125 * W**2 + 8) / 4000, (225 * W**2 + 8) / 36000, (225 * W**2 + 1) / 12000]
        expected_det = ((432000 * rho - 12170) * W**2 + 867 - 4800 * rho) / 138240000
    else:
        expected_rhs = [(672 * W**2 + 25) / 20580, (189 * W**2 + 4) / 24696, (819 * W**2 + 2) / 41160]
        expected_det = (39690 * W**4 + (2381400 * rho - 60594) * W**2 + 5183 - 22050 * rho) / 777924000
    check(prefix + ": exact RHS interpolation gives the stated quadratic formulas",
          all(sp.expand(x - y) == 0 for x, y in zip(rhs, expected_rhs)))
    determinant = sp.factor(sp.det(sp.Matrix([[xi[i], 1, rhs[i]] for i in range(3)])))
    check(prefix + ": exact compatibility determinant has the stated polynomial",
          sp.expand(determinant - expected_det) == 0)
    print("  xi coefficients:", xi, flush=True)
    print("  RHS:", rhs, flush=True)
    print("  compatibility determinant:", determinant, flush=True)
    return sp.Poly(sp.together(determinant).as_numer_denom()[0], W)


def main():
    start = time.time()
    print("DD. One mixed-Weyl pencil, finite same-witness obstruction", flush=True)
    print("Fresh portable path: 16 exact node systems, once in memory; no external caches.", flush=True)
    curvature_checks()
    check("rank helper uses the same three unisolvent coefficient nodes", tuple(Rank.NODES) == NODES[:3])
    assembled = 0
    for rho in (1, -1):
        polynomials = []
        for layout in LAYOUTS:
            values = [system(rho, node, layout) for node in NODES]
            assembled += len(values)
            polynomials.append(analyse(rho, layout, values))
            summary = Rank.verify_rank_case(rho, layout, [A for A, _ in values[:3]],
                                            emit=lambda text: print(text, flush=True))
            print("  rank certificate summary:", summary, flush=True)
        left, right, unit = sp.gcdex(polynomials[0], polynomials[1])
        check(f"rho={rho}: the two compatibility numerators have an exact unit Bezout identity",
              unit.as_expr() == 1
              and sp.expand(left.as_expr() * polynomials[0].as_expr()
                            + right.as_expr() * polynomials[1].as_expr() - 1) == 0)
        print("  Bezout coefficients:", left.as_expr(), right.as_expr(), flush=True)
    check("all 16 exact node systems were assembled once", assembled == 16)
    print("At least one of the two finite fibres is incompatible for every real w, with coefficient rank 42 at the same fibre.", flush=True)
    print("Higher-jet transfer uses the curvature-independent BB3 support/parity/gradient proof, not even-sign symmetry of arbitrary higher jets.", flush=True)
    print("The result concerns one declared mixed-Weyl frame orbit and S.0's slaved-vector sector; no general-Weyl theorem is asserted.", flush=True)
    print(f"\n{len(CHECKS)} main checks, all True; rank-helper checks also passed; runtime {time.time() - start:.1f} s", flush=True)


if __name__ == "__main__":
    main()
