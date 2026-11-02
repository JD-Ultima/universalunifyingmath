"""Exact fibre obstruction for a quartic metric with zero curvature at 0.

Base four-jet: g=(1+beta*(p.x)^4)*I, |p|=1, beta!=0.
Only an origin-fibre obstruction is claimed. Coefficient rank is 28, so this
does not supply W's full-rank base-neighbourhood persistence argument.
Uses the reviewed W point helper, U assembly and xcheck_round10 engine.
Run from the repository root with Python -X utf8 -I.
"""
import pathlib
import sys
import time
from fractions import Fraction as F

import sympy as sp

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import W_quartic_rank_one as W

ROWS = (5, 20, 191)
WEIGHTS = (4, -3, 2)
MINOR_ROWS = (0, 1, 4, 5, 6, 7, 9, 10, 14, 20, 21, 28, 35, 42, 46,
              47, 48, 49, 56, 57, 70, 84, 98, 112, 126, 140, 154, 168, 191)
MINOR_COLS = tuple(range(28)) + (42,)
DETERMINANT = F(30517578125, 13249474533898474022240256)
CHECKS = []


def check(label, condition):
    ok = bool(condition)
    CHECKS.append(ok)
    print(f"[{len(CHECKS):02d}] {label}: {ok}", flush=True)
    assert ok, label


def minor(A, b):
    augmented = [row + [val] for row, val in zip(A, b)]
    selected = [[augmented[row][col] for col in MINOR_COLS] for row in MINOR_ROWS]
    return sp.QQ.to_sympy(W.U.qmatrix(selected).det())


def main():
    start = time.time()
    print("Y. A quartic fibre obstruction at a zero-curvature base point", flush=True)
    print(f"sympy {sp.__version__}; rho=0, lambda=1/5, cH=2, cV=3, y=diag{W.U.YSTD}", flush=True)
    check("rows encode (0,1,5), (0,2,6), (1,2,9)",
          ROWS == (W.row_index(0, 1, 5), W.row_index(0, 2, 6), W.row_index(1, 2, 9)))
    A, b0 = W.U.relaxed(W.point(0, 0))
    A1, b1 = W.U.relaxed(W.point(0, 1))
    Am2, bm2 = W.U.relaxed(W.point(0, -2))
    slope = W.predicted_slope()
    check("coefficient matrix fixed at beta=0,1,-2", A == A1 == Am2)
    check("flat control has zero RHS", not any(b0))
    check("beta=1 RHS is the proved six-row slope", b1 == slope)
    check("beta=-2 control equals -2 times the slope", bm2 == [-2 * val for val in slope])
    check("coefficient rank is 28", W.U.qmatrix(A).rank() == 28)
    check("augmented rank is 29 at beta=1",
          W.U.qmatrix([row + [val] for row, val in zip(A, b1)]).rank() == 29)
    sparse = ({4: F(1, 4), 8: F(-1, 6)},
              {4: F(1, 3), 11: F(-1, 6)},
              {8: F(1, 3), 11: F(-1, 4)})
    check("all 42 coefficients match the three sparse equations",
          all(A[row] == [values.get(col, F(0)) for col in range(42)]
              for row, values in zip(ROWS, sparse)))
    relation = [sum(F(weight) * A[row][col] for row, weight in zip(ROWS, WEIGHTS))
                for col in range(42)]
    defect = sum(F(weight) * b1[row] for row, weight in zip(ROWS, WEIGHTS))
    check("weights (4,-3,2) annihilate every coefficient column", not any(relation))
    check("the RHS combination is -3 at beta=1", defect == -3)
    determinant = minor(A, b1)
    check("recorded augmented 29-minor is nonzero", determinant == sp.Rational(DETERMINANT))
    check("the beta=-2 minor scales by -2", minor(A, bm2) == -2 * determinant)
    print(f"  rows {ROWS}; weights {WEIGHTS}; RHS defect = -3*beta", flush=True)
    print(f"  29-minor rows: {MINOR_ROWS}", flush=True)
    print(f"  29-minor columns: {MINOR_COLS}", flush=True)
    print(f"  29-minor determinant at beta=1: {determinant}", flush=True)
    print("The jet-order proof gives A(beta)=A(0), b(beta)=beta*b(1); every nonzero beta is inconsistent at this witness.", flush=True)
    print("GL(4) covariance gives coefficient rank 28 throughout the positive fibre; the 29-minor excludes a fibre-open consistent patch.", flush=True)
    print("No base-neighbourhood persistence or local Upsilon exclusion is inferred from ranks (28,29).", flush=True)
    print(f"\n{len(CHECKS)} checks, all True; runtime {time.time() - start:.1f} s", flush=True)


if __name__ == "__main__":
    main()
