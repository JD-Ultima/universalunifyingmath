"""Exact algebra for one declared diagonal Weyl direction.

The base metric is the complete normal quadratic polynomial for R=rho*R0+w*W0.
The hand jet-degree bound is A(w) of degree at most two and RHS of degree at most
three, at a fixed fibre point. Matrix entries are reconstructed only within these
proved bounds; coefficient determinants are computed over QQ[w], not fitted.
All results remain limited to this declared curvature direction until reviewed.
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
from sympy.polys.matrices import DomainMatrix

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import U_rank_one_frontier as U

W = sp.Symbol("w")
ETA = (1, 1, 1, 1)
ROWS = (168, 337, 492)
PIVOT38 = [4, 18, 32, 42, 46, 50, 60, 74, 88, 98, 140, 168, 186, 200, 211, 214,
           267, 309, 337, 354, 366, 368, 422, 464, 508, 633, 648, 662, 676, 704,
           718, 746, 775, 776, 789, 817, 831, 859]
EXTRA5 = [61, 229, 383, 523]
EXTRA6 = [76, 244, 398, 538]


def curvature_tensors():
    g = U.X.gmat(ETA)
    P = [[(Fr(int(i == j == 0)) - g[i][j] / 6) / 2 for j in range(4)] for i in range(4)]
    R0 = U.X.kn(P, g)
    pairs = {(0, 1): 1, (2, 3): 1, (0, 2): -1, (1, 3): -1}
    forms = {pair: [[Fr(int((i, j) == pair) - int((j, i) == pair)) for j in range(4)]
                     for i in range(4)] for pair in pairs}
    Weyl = {key: sum(Fr(pairs[pair]) * forms[pair][key[0]][key[1]] * forms[pair][key[2]][key[3]]
                     for pair in pairs) for key in R0}
    assert all(sum(Weyl[(a, b, a, d)] for a in range(4)) == 0 for b in range(4) for d in range(4))
    assert all(Weyl[(a, b, c, d)] + Weyl[(a, c, d, b)] + Weyl[(a, d, b, c)] == 0
               for a, b, c, d in itertools.product(range(4), repeat=4))
    return R0, Weyl


R0, WEYL = curvature_tensors()


def system(rho, w, ydiag):
    R = {key: Fr(rho) * R0[key] + Fr(w) * WEYL[key] for key in R0}
    base = U.X.weyl_base(R, ETA)
    data = U.X.point_data(base, U.diagonal(ydiag), Fr(1, 5), Fr(2), Fr(3), riem=True)
    return U.relaxed(data)


def polynomial_entry(minus, zero, plus):
    return sp.Rational(zero) + W * sp.Rational((plus - minus) / 2) \
        + W**2 * sp.Rational((plus + minus) / 2 - zero)


def analyse(rho, ydiag, coefficient_minor=True):
    start = time.time()
    values = {node: system(rho, node, ydiag) for node in (-1, 0, 1, 2)}
    Am, bm = values[-1]
    A0, b0 = values[0]
    Ap, bp = values[1]
    A2, b2 = values[2]
    allowed = {27, 32, 36, 39, 41}
    # Three nodes decide the unwanted coefficient polynomials (proved degree <=2).
    assert all(A[row][col] == 0 for A in (Am, A0, Ap) for row in ROWS for col in range(42) if col not in allowed)
    a = [sp.expand(polynomial_entry(Am[row][27], A0[row][27], Ap[row][27])) for row in ROWS]
    b = [sp.expand(sp.interpolate([(node, sp.Rational(values[node][1][row])) for node in (-1, 0, 1, 2)], W)) for row in ROWS]
    assert all(sp.degree(poly, W) <= 2 for poly in a)
    assert all(sp.degree(poly, W) <= 3 for poly in b)
    assert all(poly.subs(W, 2) == A2[row][27] for poly, row in zip(a, ROWS))
    assert all(A2[row][col] == 3 * Ap[row][col] + Am[row][col] - 3 * A0[row][col]
               for row in range(1274) for col in range(42))
    assert all(A[row][col] == A[ROWS[0]][col] for A in (Am, A0, Ap) for row in ROWS for col in (32, 36, 39, 41))
    compatibility = sp.factor(sp.det(sp.Matrix([[a[i], 1, b[i]] for i in range(3)])))
    print(f"rho={rho}, y=diag{ydiag}", flush=True)
    print(f"  xi coefficients: {[sp.factor(poly) for poly in a]}", flush=True)
    print(f"  RHS: {[sp.factor(poly) for poly in b]}", flush=True)
    print(f"  compatibility determinant: {compatibility}", flush=True)
    result = dict(rho=rho, ydiag=tuple(ydiag), a=a, b=b, compatibility=compatibility, values=values)
    if coefficient_minor:
        ring = sp.QQ.poly_ring(W)
        for tag, extra in (("V5", EXTRA5), ("V6", EXTRA6)):
            selected = PIVOT38 + extra
            polyrows = [[polynomial_entry(Am[row][col], A0[row][col], Ap[row][col]) for col in range(42)] for row in selected]
            matrix = DomainMatrix.from_list_sympy(42, 42, polyrows).convert_to(ring)
            det = sp.factor(ring.to_sympy(matrix.det()))
            print(f"  coefficient determinant {tag}: {det}", flush=True)
            result[tag] = det
    print(f"  case runtime {time.time()-start:.1f} s", flush=True)
    return result


def rational_coefficient_matrix(case, value):
    """Evaluation of the exactly reconstructed full coefficient polynomial."""
    value = Fr(value)
    Am = case["values"][-1][0]
    A0 = case["values"][0][0]
    Ap = case["values"][1][0]
    return [[zero + value * (plus - minus) / 2 + value**2 * ((plus + minus) / 2 - zero)
             for minus, zero, plus in zip(rm, r0, rp)] for rm, r0, rp in zip(Am, A0, Ap)]


def numerator(expression):
    return sp.Poly(sp.together(expression).as_numer_denom()[0], W).monic()


def combined_cases(rho):
    cases = [analyse(rho, layout) for layout in ((2, 3, 4, 5), (2, 3, 5, 7))]
    compatibility_gcd = sp.gcd(numerator(cases[0]["compatibility"]), numerator(cases[1]["compatibility"]))
    print(f"rho={rho}: compatibility gcd = {compatibility_gcd.as_expr()}", flush=True)
    assert compatibility_gcd.degree() == 0
    products = [numerator(case[tag] * case["compatibility"]) for case in cases for tag in ("V5", "V6")]
    common = products[0]
    for product in products[1:]:
        common = sp.gcd(common, product)
    common = common.sqf_part()
    print(f"rho={rho}: exceptional product gcd = {sp.factor(common.as_expr())}", flush=True)
    factors = sp.factor_list(common.as_expr(), W)[1]
    assert all(sp.degree(factor, W) == 1 for factor, _ in factors)
    roots = sorted([sp.solve(factor, W)[0] for factor, _ in factors])
    for root in roots:
        accepted = False
        for case in cases:
            if case["compatibility"].subs(W, root) == 0:
                continue
            A = rational_coefficient_matrix(case, Fr(root))
            rank = U.qmatrix(A).rank()
            print(f"  w={root}, y=diag{case['ydiag']}: exact coefficient rank {rank}, nonzero compatibility determinant", flush=True)
            if rank == 42:
                accepted = True
                break
        assert accepted, f"unresolved rank at rho={rho}, w={root}"
    print(f"rho={rho}: every exceptional value is resolved; the two witnesses exclude every real w", flush=True)
    return cases


if __name__ == "__main__":
    start = time.time()
    combined_cases(1)
    combined_cases(-1)
    print("All exact assertions passed for both signs and every real w in the declared quadratic family.", flush=True)
    print(f"runtime {time.time()-start:.1f} s", flush=True)
