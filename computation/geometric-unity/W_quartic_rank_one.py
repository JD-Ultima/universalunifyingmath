"""Exact three-row obstruction with arbitrary symmetric quartic base terms.

Base four-jet: I-rho R0(p)_iajb x^a x^b/3 + Q4(x), |p|=1, rho!=0,
where Q4 is an arbitrary homogeneous symmetric quartic matrix and cubic terms
are zero. Higher remainders with zero derivatives through order four are allowed.
The hand jet-order proof makes the relaxed coefficient matrix independent of
Q4. For each diagonal vertical coordinate, its cubic HV variation is a gradient;
the principal third-jet terms therefore cancel in the three certificate RHS rows.
The script checks this with a symbolic 350-coefficient Q4 and mixed-monomial
referee-engine controls. The one-parameter beta*(p.x)^4*I calculation is retained
as a control with a nonzero predicted six-row RHS variation.

Only accepted U helpers and the earlier referee xcheck_round10 engine are used.
The MD shift in that engine retains the required fourth base derivatives as
coefficients of epsilon*x^3; the nonzero predicted slope checks this explicitly.
No random inputs. Run from the repository root with Python -X utf8 -I.
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

ROWS = (168, 337, 492)
WEIGHTS = {1: (-445, 41538, -41093), -1: (445, 44862, -45307)}
DEFECTS = {1: Fr(3933, 800), -1: Fr(5667, 800)}
CHECKS = []


def check(label, condition):
    ok = bool(condition)
    CHECKS.append(ok)
    print(f"[{len(CHECKS):02d}] {label}: {ok}", flush=True)
    assert ok, label


def base_metric(rho, beta):
    base = U.base_metric(Fr(rho))
    base[(4, 0, 0, 0)] = U.X.gmat((1, 1, 1, 1), Fr(beta))
    return base


def point(rho, beta, cH=Fr(2), cV=Fr(3), yscale=Fr(1)):
    return U.X.point_data(base_metric(rho, beta), U.diagonal([yscale * val for val in U.YSTD]),
                          Fr(1, 5), Fr(cH), Fr(cV), riem=True)


def row_index(a, b, d):
    assert 0 <= a < b < 14 and 0 <= d < 14
    return 14 * (13 * a - a * (a - 1) // 2 + b - a - 1) + d


def predicted_slope(cH=Fr(2), cV=Fr(3), ydiag=U.YSTD):
    """The principal-jet formula proved in W2, not finite-difference fitting."""
    slope = [Fr(0)] * 1274
    d0 = Fr(ydiag[0])
    for i in (1, 2, 3):
        di = Fr(ydiag[i])
        vertical = 4 + U.X.AOF[(0, i)]
        slope[row_index(0, i, vertical)] = 6 * (di - d0) / (cH * d0)
        slope[row_index(0, vertical, i)] = 12 * cV * (di - d0) / (cH**2 * d0**2 * di**2)
    return slope


def signed_certificate(sign):
    systems = []
    for beta in (Fr(0), Fr(1), Fr(-2)):
        data = point(sign, beta)
        systems.append(U.relaxed(data))
    A, b0 = systems[0]
    A1, b1 = systems[1]
    Am2, bm2 = systems[2]
    slope = predicted_slope()
    check(f"rho={sign}: coefficient matrix unchanged at beta=0,1,-2", A == A1 == Am2)
    check(f"rho={sign}: entire beta=1 RHS variation equals the proved six-row slope",
          [y - x for x, y in zip(b0, b1)] == slope)
    check(f"rho={sign}: independent beta=-2 control agrees with affine formula",
          bm2 == [x - 2 * y for x, y in zip(b0, slope)])
    check(f"rho={sign}: all three certificate RHS rows have zero beta slope",
          all(slope[row] == 0 for row in ROWS))
    weights = WEIGHTS[sign]
    left = [sum(Fr(weight) * A[row][col] for row, weight in zip(ROWS, weights)) for col in range(42)]
    defect = sum(Fr(weight) * b0[row] for row, weight in zip(ROWS, weights))
    check(f"rho={sign}: integer combination of the three coefficient rows is zero", all(val == 0 for val in left))
    check(f"rho={sign}: the same RHS combination is {DEFECTS[sign]} != 0", defect == DEFECTS[sign] and defect != 0)
    check(f"rho={sign}: coefficient rank remains 42", U.qmatrix(A).rank() == 42)
    print(f"  rows {ROWS}; weights {weights}; RHS defect {defect}", flush=True)
    print("  the hand affinity and zero-slope arguments extend this exact contradiction to every real beta", flush=True)
    return A, b0


def general_gradient_identity():
    """Exact symbolic cancellation with all 350 homogeneous quartic coefficients."""
    xs = sp.symbols("X0:4")
    exponents = [e for e in itertools.product(range(5), repeat=4) if sum(e) == 4]
    monomials = [sp.prod(x**k for x, k in zip(xs, e)) for e in exponents]
    coefficients = sp.symbols("q0:350")
    Q = [[None] * 4 for _ in range(4)]
    pos = 0
    for a in range(4):
        for b in range(a, 4):
            Q[a][b] = Q[b][a] = sum(coefficients[pos + k] * mono for k, mono in enumerate(monomials))
            pos += len(monomials)
    cv, d3, lam = sp.symbols("cV d3 lambda", nonzero=True)
    # Cubic connection variation at a Euclidean normal base point.
    def gamma(a, j, b):
        return (sp.diff(Q[a][b], xs[j]) + sp.diff(Q[a][j], xs[b]) - sp.diff(Q[j][b], xs[a])) / 2
    F = -cv / d3 * (Q[3][3] - lam * sum(Q[a][a] for a in range(4)))
    Hv = [-2 * cv / d3 * (gamma(3, j, 3) - lam * sum(gamma(a, j, a) for a in range(4))) for j in range(4)]
    check("all 350 quartic coefficients: the diagonal-vertical HV variation is a gradient",
          pos == 350 and all(sp.expand(Hv[j] - sp.diff(F, xs[j])) == 0 for j in range(4)))
    weights = sp.symbols("inverse_H0:4")
    cancellation = []
    for i in range(3):
        term = sum(weights[c] * (sp.diff(Hv[c], xs[c], xs[i], xs[i])
                                 - sp.diff(Hv[i], xs[c], xs[c], xs[i])) / 2 for c in range(4))
        cancellation.append(sp.expand(term))
    check("all 350 quartic coefficients: the three principal mixed Ricci derivatives vanish",
          all(val == 0 for val in cancellation))


def mixed_quartic_control(sign, baseline):
    """Full symmetric matrices and mixed monomials, outside the beta control line."""
    A0, b0 = baseline
    base = U.base_metric(sign)
    exponents = ((1, 1, 1, 1), (2, 0, 0, 2), (0, 3, 1, 0))
    matrices = (
        [[Fr((i + 1) * (j + 1) - 3 * (i + j + 1), 7) for j in range(4)] for i in range(4)],
        [[Fr(2 * (i + j + 2) + int(i == j), 5) for j in range(4)] for i in range(4)],
        [[Fr((i + 1) * (j + 1) + 3 * int(i == j), 11) for j in range(4)] for i in range(4)],
    )
    for exponent, matrix in zip(exponents, matrices):
        base[exponent] = matrix
    data = U.X.point_data(base, U.diagonal(U.YSTD), Fr(1, 5), Fr(2), Fr(3), riem=True)
    A, b = U.relaxed(data)
    changed = sum(left != right for left, right in zip(b, b0))
    check(f"rho={sign}: mixed quartics leave the full coefficient matrix fixed", A == A0)
    check(f"rho={sign}: mixed quartics leave all three certificate RHS rows fixed",
          all(b[row] == b0[row] for row in ROWS))
    check(f"rho={sign}: mixed quartics change other RHS rows (count {changed})", changed > 0)
    defect = sum(Fr(weight) * b[row] for row, weight in zip(ROWS, WEIGHTS[sign]))
    check(f"rho={sign}: mixed-quartic defect remains {DEFECTS[sign]}", defect == DEFECTS[sign])


def scaling_control(rho, beta, cH, cV, s):
    rho, beta, cH, cV, s = map(Fr, (rho, beta, cH, cV, s))
    signed = 1 if rho > 0 else -1
    canonical_beta = beta / rho**2
    alpha = Fr(3) / cV
    q = Fr(2) * cV / (Fr(3) * cH)
    yscale = q / s**2
    canonical = point(signed, canonical_beta)
    target = point(rho, beta, cH, cV, yscale)
    jac = [s] * 4 + [yscale] * 10
    check(f"rho={rho}, beta={beta}: quartic invariant is beta/rho^2={canonical_beta}",
          rho * s**2 == signed and beta * s**4 == canonical_beta)
    check(f"rho={rho}: pulled-back metric and covariant Ricci agree",
          all(canonical["g"][i][j] == alpha * jac[i] * jac[j] * target["g"][i][j]
              and canonical["Ric"][i][j] == jac[i] * jac[j] * target["Ric"][i][j]
              for i in range(14) for j in range(14)))
    keys = set(canonical["R"]) | set(target["R"])
    check(f"rho={rho}: pulled-back mixed Riemann agrees",
          all(canonical["R"].get(key, Fr(0)) == jac[key[1]] * jac[key[2]] * jac[key[3]] / jac[key[0]]
              * target["R"].get(key, Fr(0)) for key in keys))
    check(f"rho={rho}: pulled-back covariant derivative of Ricci agrees",
          all(canonical["nRic"][k][i][j] == jac[k] * jac[i] * jac[j] * target["nRic"][k][i][j]
              for k in range(14) for i in range(14) for j in range(14)))


def main():
    start = time.time()
    print("W. Arbitrary symmetric quartic base terms: an exact three-row obstruction", flush=True)
    print(f"sympy {sp.__version__}; canonical lambda=1/5, cH=2, cV=3, y=diag{U.YSTD}", flush=True)
    check("row indices encode (0,13,0), (1,13,1), (2,13,2)",
          ROWS == tuple(row_index(i, 13, i) for i in range(3)))
    print("  nonzero predicted RHS beta slopes:", [(i, str(val)) for i, val in enumerate(predicted_slope()) if val], flush=True)
    general_gradient_identity()
    canonical = {sign: signed_certificate(sign) for sign in (1, -1)}
    for sign in (1, -1):
        mixed_quartic_control(sign, canonical[sign])
    scaling_control(Fr(4), Fr(48, 5), Fr(5), Fr(7), Fr(1, 2))
    scaling_control(Fr(-9), Fr(-162, 7), Fr(3, 2), Fr(5, 2), Fr(1, 3))
    print("Every homogeneous symmetric quartic Q4 is excluded with the fixed U quadratic jet and zero cubic terms at rho!=0. Arbitrary quadratic/Weyl and cubic jets, and the rho=0 frontier, are not classified.", flush=True)
    print(f"\n{len(CHECKS)} checks, all True; runtime {time.time() - start:.1f} s", flush=True)


if __name__ == "__main__":
    main()
