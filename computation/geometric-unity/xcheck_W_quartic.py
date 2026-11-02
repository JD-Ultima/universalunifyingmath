"""Independent mixed-quartic inputs for the W internal review.

Uses the accepted U relaxed-system assembly and xcheck_round10 curvature
engine. This is not a second curvature implementation. The hand gradient
argument, rather than these finite controls, proves the quartic statement.
Run from the repository root with Python -X utf8 -I.
"""
import pathlib
import sys
from fractions import Fraction as F

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import U_rank_one_frontier as U

ROWS = (168, 337, 492)
WEIGHTS = {1: (-445, 41538, -41093), -1: (445, 44862, -45307)}


def main():
    print("W internal review: independent mixed-quartic inputs; same accepted curvature engine")
    for sign in (1, -1):
        A0, b0 = U.relaxed(U.point(sign, F(2), F(3)))
        base = U.base_metric(sign)
        # Q4 = x0*x1*x2*x3 M1 + x0^2*x3^2 M2 + x1^3*x2 M3.
        # Each matrix is symmetric, with nonzero diagonal and off-diagonal entries.
        base[(1, 1, 1, 1)] = [[F((i + 1) * (j + 1), 7) for j in range(4)]
                               for i in range(4)]
        base[(2, 0, 0, 2)] = [[F((i + j + 1) * (-1)**(i + j), 11) for j in range(4)]
                               for i in range(4)]
        base[(0, 3, 1, 0)] = [[F(2 * i + 3, 5) if i == j else F(i + j + 1, 13)
                               for j in range(4)] for i in range(4)]
        data = U.X.point_data(base, U.diagonal(U.YSTD), F(1, 5), F(2), F(3), riem=True)
        A, b = U.relaxed(data)
        assert A == A0
        assert all(b[row] == b0[row] for row in ROWS)
        changed = sum(x != y for x, y in zip(b, b0))
        assert changed > 0
        defect = sum(F(weight) * b[row] for row, weight in zip(ROWS, WEIGHTS[sign]))
        assert defect == (F(3933, 800) if sign == 1 else F(5667, 800))
        print(f"rho={sign}: A and three certificate RHS rows fixed; "
              f"{changed} other RHS rows changed; defect {defect}")
    print("8 assertions passed; no arbitrary-cubic or arbitrary-Weyl conclusion follows.")


if __name__ == "__main__":
    main()
