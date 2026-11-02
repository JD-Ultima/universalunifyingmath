"""Exact cubic-jet extension of the rank-one, zero-Weyl obstruction.

The quadratic normal-coordinate jet is U's rho*R_p. All 200 cubic and 350
quartic metric coefficients are free. The hand support and parity proof is in
Z-cubic-rank-one.md; the calculation supplies the rank certificates and checks
the generic cubic diagonal-vertical cancellation symbolically.

Only U helpers and the earlier independent referee xcheck_round10 engine are
imported. Mixed cubic, reflection and cubic-plus-quartic controls supplement the
proof, not a sampled degree inference. No randomness; exact rational arithmetic.
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

DIAGONAL_V = {4, 8, 11, 13}
PIVOT38 = [4, 18, 32, 42, 46, 50, 60, 74, 88, 98, 140, 168, 186, 200, 211, 214,
           267, 309, 337, 354, 366, 368, 422, 464, 508, 633, 648, 662, 676, 704,
           718, 746, 775, 776, 789, 817, 831, 859]
EXTRA4 = [61, 229, 383, 523]
CERT_ROWS = (168, 337, 492)
CERT_WEIGHTS = {1: (-445, 41538, -41093), -1: (445, 44862, -45307)}
CERT_DEFECT = {1: Fr(3933, 800), -1: Fr(5667, 800)}
CHECKS = []


def check(label, condition):
    ok = bool(condition)
    CHECKS.append(ok)
    print(f"[{len(CHECKS):02d}] {label}: {ok}", flush=True)
    assert ok, label


def row_metadata():
    return [(a, b, d) for a in range(14) for b in range(a + 1, 14) for d in range(14)]


META = row_metadata()


def row_class(a, b, d):
    if a >= 4 and b >= 4 and d >= 4:
        return "VVV"
    if a < 4 and b in DIAGONAL_V and d < 4:
        return "H-D-H"
    if a < 4 and b >= 4 and d in DIAGONAL_V:
        return "H-V-D"
    if a < 4 and b < 4 and d in DIAGONAL_V:
        return "H-H-D"
    if a in DIAGONAL_V and b in DIAGONAL_V and d < 4:
        return "D-D-H"
    return None


INVARIANT = [r for r, triple in enumerate(META) if row_class(*triple)]


def symbolic_cubic_gradient():
    xs = sp.symbols("X0:4")
    monomials = [sp.prod(x**power for x, power in zip(xs, e))
                 for e in itertools.product(range(4), repeat=4) if sum(e) == 3]
    coefficients = sp.symbols("q0:200")
    Q = [[None] * 4 for _ in range(4)]
    pos = 0
    for i in range(4):
        for j in range(i, 4):
            Q[i][j] = Q[j][i] = sum(coefficients[pos + k] * mon for k, mon in enumerate(monomials))
            pos += len(monomials)
    cv, lam = sp.symbols("cV lambda", nonzero=True)
    ds = sp.symbols("d0:4", nonzero=True)

    def gamma(a, j, b):
        return (sp.diff(Q[a][b], xs[j]) + sp.diff(Q[a][j], xs[b]) - sp.diff(Q[j][b], xs[a])) / 2

    gradient_ok = pos == 200
    curvature_ok = True
    for vertical in range(4):
        F = -cv / ds[vertical] * (Q[vertical][vertical] - lam * sum(Q[a][a] for a in range(4)))
        Hv = [-2 * cv / ds[vertical] * (gamma(vertical, j, vertical)
               - lam * sum(gamma(a, j, a) for a in range(4))) for j in range(4)]
        gradient_ok &= all(sp.expand(Hv[j] - sp.diff(F, xs[j])) == 0 for j in range(4))
        for i, j, k in itertools.product(range(4), repeat=3):
            curvature_ok &= sp.expand(sp.diff(Hv[i], xs[j], xs[k]) - sp.diff(Hv[j], xs[i], xs[k])) == 0
    check("all 200 cubic coefficients: every diagonal-vertical HV variation is a gradient", gradient_ok)
    check("all 200 cubic coefficients: its three-horizontal Riemann variation vanishes", curvature_ok)


def kernel_basis():
    K = sp.zeros(42, 4)
    c = Fr(12, 5)
    for i, di in enumerate(U.YSTD):
        K[14 + i, i] = sp.Rational(c * 2 * di)
        K[28 + i, i] = 1
    return K


def mixed_base(sign, cubic_sign=1, quartic=False):
    base = U.base_metric(sign)
    cubic_exponents = ((1, 1, 1, 0), (0, 0, 2, 1), (2, 0, 0, 1))
    for k, exponent in enumerate(cubic_exponents):
        base[exponent] = [[cubic_sign * Fr((i + 1) * (j + 1) + (k + 1) * (i + j + 2), 5 + 2 * k)
                           for j in range(4)] for i in range(4)]
    if quartic:
        for k, exponent in enumerate(((1, 1, 1, 1), (2, 0, 0, 2), (0, 3, 1, 0))):
            base[exponent] = [[Fr(2 * (i + 1) * (j + 1) - (k + 2) * (i + j + 2), 7 + 2 * k)
                               for j in range(4)] for i in range(4)]
    return base


def data_for(base):
    return U.X.point_data(base, U.diagonal(U.YSTD), Fr(1, 5), Fr(2), Fr(3), riem=True)


def determinant42(A):
    return sp.QQ.to_sympy(U.qmatrix([A[i] for i in PIVOT38 + EXTRA4]).det())


def signed_checks(sign):
    data0 = U.point(sign, 2, 3)
    A0, b0 = U.relaxed(data0)
    M = U.qmatrix([A0[r] for r in INVARIANT])
    K = kernel_basis()
    check(f"rho={sign}: invariant subsystem has exact rank 38", M.rank() == 38)
    check(f"rho={sign}: the recorded 38 rows span that subsystem", U.qmatrix([A0[r] for r in PIVOT38]).rank() == 38)
    check(f"rho={sign}: the sparse four-vector basis is its kernel", M.to_Matrix() * K == sp.zeros(722, 4) and K.rank() == 4)
    restriction = U.qmatrix([A0[r] for r in EXTRA4]).to_Matrix() * K
    expected = sp.diag(-sp.Rational(1, 45), -sp.Rational(1, 30), -sp.Rational(2, 45), -sp.Rational(1, 18))
    check(f"rho={sign}: four extra rows restrict to the nonzero diagonal matrix", restriction == expected)
    Rh = U.X.mmul(data0["gi"], data0["Ric"])
    check(f"rho={sign}: vertical Einstein defect S^5_5 is 1/432", Rh[5][5] + Fr(5, 12) == Fr(1, 432))
    det0 = determinant42(A0)
    check(f"rho={sign}: the recorded 42-by-42 coefficient minor is nonzero", det0 != 0)
    print(f"  rho={sign}: coefficient determinant {det0}", flush=True)
    left = [sum(Fr(w) * A0[r][col] for r, w in zip(CERT_ROWS, CERT_WEIGHTS[sign])) for col in range(42)]
    defect = sum(Fr(w) * b0[r] for r, w in zip(CERT_ROWS, CERT_WEIGHTS[sign]))
    check(f"rho={sign}: three-row certificate has defect {CERT_DEFECT[sign]}",
          all(x == 0 for x in left) and defect == CERT_DEFECT[sign])

    Ap, bp = U.relaxed(data_for(mixed_base(sign, 1)))
    Am, bm = U.relaxed(data_for(mixed_base(sign, -1)))
    changed_A = sum(x != y for ra, rb in zip(A0, Ap) for x, y in zip(ra, rb))
    changed_b = sum(x != y for x, y in zip(b0, bp))
    check(f"rho={sign}: mixed cubic control changes other entries ({changed_A} coefficients, {changed_b} RHS rows)",
          changed_A > 0 and changed_b > 0)
    check(f"rho={sign}: all 722 invariant coefficient rows are fixed", all(Ap[r] == A0[r] for r in INVARIANT))
    check(f"rho={sign}: extra-row variations occur only in xi_5 (column 19)",
          all(Ap[r][c] == A0[r][c] for r in EXTRA4 for c in range(42) if c != 19))
    check(f"rho={sign}: the fixed 42-minor is unchanged by the mixed cubic tensor", determinant42(Ap) == det0)
    check(f"rho={sign}: all selected coefficient rows and RHS are fixed",
          all(Ap[r] == A0[r] and bp[r] == b0[r] for r in CERT_ROWS))
    row_sign = [(-1) ** sum(int(index < 4) for index in triple) for triple in META]
    col_sign = [(-1 if c % 14 < 4 else 1) for c in range(42)]
    check(f"rho={sign}: full-system base-reflection parity holds",
          all(Am[r][c] == row_sign[r] * col_sign[c] * Ap[r][c] for r in range(1274) for c in range(42))
          and all(bm[r] == row_sign[r] * bp[r] for r in range(1274)))
    check(f"rho={sign}: cubic dependence is affine in the independent mixed control",
          all(Ap[r][c] + Am[r][c] == 2 * A0[r][c] for r in range(1274) for c in range(42))
          and all(bp[r] + bm[r] == 2 * b0[r] for r in range(1274)))
    Aq, bq = U.relaxed(data_for(mixed_base(sign, 1, quartic=True)))
    check(f"rho={sign}: additional mixed quartics preserve the cubic coefficient matrix", Aq == Ap)
    check(f"rho={sign}: cubic-plus-quartic selected RHS stays fixed", all(bq[r] == b0[r] for r in CERT_ROWS))


def main():
    start = time.time()
    print("Z. Arbitrary cubic and quartic jets with rank-one Ricci and zero Weyl", flush=True)
    counts = {kind: sum(row_class(*triple) == kind for triple in META)
              for kind in ("VVV", "H-D-H", "H-V-D", "H-H-D", "D-D-H")}
    check("the invariant row classes have counts 450,64,160,24,24", list(counts.values()) == [450, 64, 160, 24, 24])
    print(f"  diagonal vertical coordinates {sorted(DIAGONAL_V)}; row counts {counts}", flush=True)
    print(f"  38 independent invariant rows: {PIVOT38}", flush=True)
    print(f"  extra rows: {EXTRA4}; triples {[META[r] for r in EXTRA4]}", flush=True)
    K = kernel_basis()
    print("  sparse kernel basis:", [[(i, str(K[i, j])) for i in range(42) if K[i, j]] for j in range(4)], flush=True)
    symbolic_cubic_gradient()
    signed_checks(1)
    signed_checks(-1)
    print("The full-rank extension rests on the proved invariant-row support and sparse kernel argument, not on sampling cubic coefficients.", flush=True)
    print("Arbitrary cubic and quartic jets are covered with the prescribed nonzero rank-one Ricci/zero-Weyl quadratic jet; nonzero Weyl and Ricci-flat nonflat points remain outside this result.", flush=True)
    print(f"\n{len(CHECKS)} checks, all True; runtime {time.time() - start:.1f} s", flush=True)


if __name__ == "__main__":
    main()
