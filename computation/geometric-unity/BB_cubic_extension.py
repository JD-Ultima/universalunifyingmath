"""Exact higher-jet support for the declared BB Weyl direction.

The all-amplitude compatibility and coefficient-minor polynomials remain in
BB_weyl_direction.py. This supplement proves polynomial kernel/restriction
identities within BB1's hand degree bound and resolves invariant rank at its
twelve finite exceptions. Mixed cubic/quartic controls supplement the universal
Z support, parity and gradient argument; they do not infer higher-jet scope.
"""
import pathlib
import sys
import time
from fractions import Fraction as Fr

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import BB_weyl_direction as B
import Z_cubic_rank_one as Z

sp, U = B.sp, B.U
W = B.W
LAYOUTS = ((2, 3, 4, 5), (2, 3, 5, 7))
EXCEPTIONS = {
    1: (Fr(-13, 21), Fr(-1, 3), Fr(-7, 39), Fr(-1, 6), Fr(1, 6), Fr(1, 3)),
    -1: (Fr(-1, 3), Fr(-1, 6), Fr(1, 21), Fr(1, 6), Fr(1, 3), Fr(19, 39)),
}
CHECKS = []


def check(label, condition):
    ok = bool(condition)
    CHECKS.append(ok)
    print(f"[{len(CHECKS):02d}] {label}: {ok}", flush=True)
    assert ok, label


def kernel_product(A, rows, layout):
    """Multiply by the four sparse K columns without dense SymPy conversion."""
    return [[A[row][14 + i] * Fr(24 * layout[i], 5) + A[row][28 + i]
             for i in range(4)] for row in rows]


def polynomial_matrix(minus, zero, plus):
    return [[sp.factor(B.polynomial_entry(a, b, c)) for a, b, c in zip(rm, rz, rp)]
            for rm, rz, rp in zip(minus, zero, plus)]


def expected_restriction(rho, layout, vertical):
    if vertical == 5:
        diagonal = [Fr(d, 90) for d in layout]
        scalar = -(3 * W + rho) ** 2
    else:
        diagonal = [Fr(9 * d, 250 if layout[2] == 5 else 360) for d in layout]
        scalar = -(3 * W - rho) ** 2
    return [[scalar * sp.Rational(diagonal[i]) if i == j else sp.S.Zero
             for j in range(4)] for i in range(4)]


def kernel_checks(rho, layout):
    values = {node: B.system(rho, node, layout) for node in (-1, 0, 1)}
    check(f"rho={rho}, y={layout}: all 722 rows annihilate K at three exact nodes",
          all(x == 0 for node in values
              for row in kernel_product(values[node][0], Z.INVARIANT, layout) for x in row))
    print("  K sparse columns:", [[(14 + i, str(Fr(24 * layout[i], 5))), (28 + i, "1")]
                                  for i in range(4)], flush=True)
    for vertical, rows in ((5, B.EXTRA5), (6, B.EXTRA6)):
        restrictions = [kernel_product(values[node][0], rows, layout) for node in (-1, 0, 1)]
        polynomial = polynomial_matrix(*restrictions)
        expected = expected_restriction(rho, layout, vertical)
        check(f"rho={rho}, y={layout}: J{vertical}K has its exact quadratic formula",
              all(sp.expand(a - b) == 0 for ra, rb in zip(polynomial, expected) for a, b in zip(ra, rb)))
        print(f"  J{vertical}K diagonal:", [polynomial[i][i] for i in range(4)], flush=True)
    return dict(rho=rho, ydiag=layout, values=values)


def exceptional_checks(rho, case):
    for value in EXCEPTIONS[rho]:
        A = B.rational_coefficient_matrix(case, value)
        rank = U.qmatrix([A[r] for r in Z.INVARIANT]).rank()
        check(f"rho={rho}, w={value}: the full invariant subsystem has rank 38", rank == 38)
    check(f"rho={rho}: J5 and J6 restriction factors have no common root",
          sp.gcd(sp.Poly((3 * W + rho) ** 2, W), sp.Poly((3 * W - rho) ** 2, W)).degree() == 0)


def mixed_base(rho, w, cubic_sign=1, quartic=False):
    curvature = {key: Fr(rho) * B.R0[key] + Fr(w) * B.WEYL[key] for key in B.R0}
    base = U.X.weyl_base(curvature, B.ETA)
    for k, exponent in enumerate(((1, 1, 1, 0), (0, 0, 2, 1), (2, 0, 0, 1))):
        base[exponent] = [[cubic_sign * Fr((i + 1) * (j + 1) + (k + 1) * (i + j + 2), 5 + 2 * k)
                           for j in range(4)] for i in range(4)]
    if quartic:
        for k, exponent in enumerate(((1, 1, 1, 1), (2, 0, 0, 2), (0, 3, 1, 0))):
            base[exponent] = [[Fr(2 * (i + 1) * (j + 1) - (k + 2) * (i + j + 2), 7 + 2 * k)
                               for j in range(4)] for i in range(4)]
    return base


def mixed_system(rho, layout, sign=1, quartic=False):
    D = U.X.point_data(mixed_base(rho, 1, sign, quartic), U.diagonal(layout),
                       Fr(1, 5), Fr(2), Fr(3), riem=True)
    return U.relaxed(D)


def mixed_checks(case):
    rho, layout = case["rho"], case["ydiag"]
    A0, b0 = case["values"][1]
    Ap, bp = mixed_system(rho, layout)
    Am, bm = mixed_system(rho, layout, -1)
    Aq, bq = mixed_system(rho, layout, quartic=True)
    prefix = f"rho={rho}, w=1, y={layout}"
    check(prefix + ": mixed cubic changes other coefficients and RHS rows",
          Ap != A0 and bp != b0)
    check(prefix + ": all 722 invariant rows and all three certificate rows are fixed",
          all(Ap[r] == A0[r] for r in Z.INVARIANT)
          and all(bp[r] == b0[r] for r in B.ROWS))
    check(prefix + ": both extra-row variations occur only in xi_5 or xi_6",
          all(Ap[r][c] == A0[r][c] for vertical, rows in ((5, B.EXTRA5), (6, B.EXTRA6))
              for r in rows for c in range(42) if c != 14 + vertical))
    check(prefix + ": both extra-row changes annihilate K",
          all(a == b for rows in (B.EXTRA5, B.EXTRA6)
              for ra, rb in zip(kernel_product(Ap, rows, layout), kernel_product(A0, rows, layout))
              for a, b in zip(ra, rb)))
    selected = B.PIVOT38 + B.EXTRA5
    det0 = U.qmatrix([A0[r] for r in selected]).det()
    detp = U.qmatrix([Ap[r] for r in selected]).det()
    check(prefix + ": the existing generic 42-minor survives the cubic control", det0 != 0 and detp == det0)
    row_sign = [(-1) ** sum(index < 4 for index in triple) for triple in Z.META]
    col_sign = [(-1 if c % 14 < 4 else 1) for c in range(42)]
    check(prefix + ": full-system reflection parity and cubic affinity hold",
          all(Am[r][c] == row_sign[r] * col_sign[c] * Ap[r][c]
              and Ap[r][c] + Am[r][c] == 2 * A0[r][c]
              for r in range(1274) for c in range(42))
          and all(bm[r] == row_sign[r] * bp[r] and bp[r] + bm[r] == 2 * b0[r] for r in range(1274)))
    check(prefix + ": added mixed quartics preserve all coefficients and selected RHS",
          Aq == Ap and all(bq[r] == b0[r] for r in B.ROWS))


def main():
    start = time.time()
    print("BB. Arbitrary cubic and quartic jets in the declared Weyl direction", flush=True)
    check("the inherited invariant subsystem contains 722 rows", len(Z.INVARIANT) == 722)
    print("  inherited 38-row list:", B.PIVOT38, flush=True)
    print("  J5:", B.EXTRA5, "; J6:", B.EXTRA6, flush=True)
    cases = {(rho, layout): kernel_checks(rho, layout) for rho in (1, -1) for layout in LAYOUTS}
    for rho in (1, -1):
        exceptional_checks(rho, cases[(rho, LAYOUTS[0])])
    mixed_checks(cases[(1, LAYOUTS[0])])
    mixed_checks(cases[(-1, LAYOUTS[1])])
    print("Universal higher-jet scope follows from the hand support/parity/gradient proof and rank-nullity, not mixed-input sampling.", flush=True)
    print("Existing BB polynomial minors handle generic amplitudes; the twelve invariant-rank checks and J5/J6 restrictions handle every remaining amplitude.", flush=True)
    print(f"\n{len(CHECKS)} checks, all True; runtime {time.time() - start:.1f} s", flush=True)


if __name__ == "__main__":
    main()
