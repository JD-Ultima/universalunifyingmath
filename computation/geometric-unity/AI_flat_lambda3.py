"""Exact controls for the flat-base Z1 single-Lambda3 endpoint classification.

AF and AG supply constant-odd endpoints. T2 excludes the three endpoints below
-1/312 when traceG>=0. At the remaining endpoint the grade-one equation,
AA7 simplicity and a rank argument force the extra even field to vanish.
The accepted flat Z1 Ricci formula supplies G=15/(2cV)*I; this helper does not
recompute that metric. It verifies the Clifford coefficients, derivatives,
scale and complete odd-baseline converse for arbitrary nonzero real kappa.
"""
import pathlib
import sys
import time

import sympy as sp

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import T_singular_residual as T

X, ALG = T.X, T.ALG
MEMBER = (1, 0, 0, 0)
CHECKS = []


def check(label, condition):
    ok = bool(condition)
    CHECKS.append(ok)
    print(f"[{len(CHECKS):02d}] {label}: {ok}", flush=True)
    assert ok, label


def grade(value, degree):
    return {mask: coefficient for mask, coefficient in value.items() if mask.bit_count() == degree}


def grade_one(field):
    return sp.Matrix([[sp.expand(item.get(1 << b, 0) / sp.I) for b in range(X.N)]
                      for item in field])


def lift(omega):
    return ALG.gamma_right([grade(ALG.mul({1 << a: 1}, omega), 4) for a in range(X.N)])


def stress(omega):
    contraction = [grade(ALG.mul({1 << a: 1}, omega), 2) for a in range(X.N)]
    norm = sp.expand(sum(value * value for value in omega.values()))
    gram = sp.Matrix([[sum(left.get(mask, 0) * right.get(mask, 0)
                           for mask in left.keys() | right.keys())
                       for right in contraction] for left in contraction])
    return norm, gram


def main():
    start = time.time()
    print("AI: flat-base Z1, both single WA Lambda3 kernel branches", flush=True)
    print(f"SymPy {sp.__version__}; real compact (14,0), cV>0, kappa!=0.", flush=True)
    kappa = sp.symbols("kappa", real=True, nonzero=True)
    coefficients = sp.symbols("G0:14", real=True)
    trace = sum(coefficients)
    endpoints = (-sp.Rational(3, 112), -sp.Rational(3, 272), -sp.Rational(1, 240))
    for a in endpoints:
        to = T.invariant(kappa * a, 0)
        residual = sp.cancel(T.pair(to, T.baseline(kappa * a, 0, kappa, coefficients)) / kappa)
        expected = a * trace + 14 * kappa**2 * (1 + 312 * a) * a**2
        check(f"endpoint a={a}: exact residual and strict-negative coefficients",
              sp.expand(residual - expected) == 0 and a < -sp.Rational(1, 312)
              and a < 0 and 14 * (1 + 312 * a) * a**2 < 0)
        print("  pairing / kappa =", sp.factor(residual), flush=True)

    q = sp.symbols("q", real=True)
    coordinate = {(1 << 0) | (1 << 1) | (1 << 2): 1}
    rotated = grade(ALG.mul(ALG.mul({1 << 0: 1, 1 << 4: 1},
                                    {1 << 1: 1, 1 << 5: 1}),
                            {1 << 2: 1, 1 << 6: -1}), 3)
    for label, omega in (("coordinate", coordinate), ("non-coordinate", rotated)):
        norm, gram = stress(omega)
        k = T.scale(lift(omega), q)
        Q = T.quadratic(k, MEMBER)
        expected = q**2 * (180 * norm * sp.eye(X.N) - sp.Rational(320, 3) * gram)
        actual = grade_one([grade(item, 1) for item in Q])
        check(label + " simple form: complete grade-one quadratic tensor",
              T.antihermitian(k) and gram.rank() == 3 and gram * gram == norm * gram
              and (actual - expected).applyfunc(sp.expand) == sp.zeros(X.N))
        check(label + " derivative controls retain only grades nine/thirteen",
              all(all(mask.bit_count() in (9, 13) for item in T.derivative(k, j, MEMBER) for mask in item)
                  for j in (0, 3, 7, 13)))

    a = -sp.Rational(1, 400)
    g = sp.symbols("g", real=True)
    baseline = T.baseline(kappa * a, 0, kappa, [g] * X.N)
    scalar = g - sp.Rational(11, 20000) * kappa**2
    check("remaining endpoint: full baseline is the stated grade-one scalar",
          all(all(mask.bit_count() == 1 for mask in item) for item in baseline)
          and (grade_one(baseline) - scalar * sp.eye(X.N)).applyfunc(sp.expand) == sp.zeros(X.N))
    cV = sp.symbols("cV", positive=True)
    equation = sp.Rational(15, 2) / cV - sp.Rational(11, 20000) * kappa**2
    scale = sp.Rational(150000, 11) / kappa**2
    check("flat Einstein normalization gives the unique positive scale",
          sp.solve(equation, cV) == [scale] and scale.is_positive is True)
    solved = T.baseline(kappa * a, 0, kappa, [sp.Rational(11, 20000) * kappa**2] * X.N)
    check("known zero-even solution has identically zero complete odd residual",
          not any(X.normalized(solved)))
    all_endpoints = endpoints + (a,)
    check("both vector determinant factors are nonzero at every endpoint",
          all(1 + 224 * value != 0 and 1 + 544 * value != 0 for value in all_endpoints))

    variables, blocks = X.templates(1)
    normalized = {key: block.subs(dict(zip(variables, (1, a, 0, 0))))
                  for key, block in blocks.items()}
    check("known solution lies on only the intended non-vector factor",
          normalized["Lam3"].rank() == 3
          and all(block.rank() == block.cols for key, block in normalized.items() if key != "Lam3"))

    print("First branch: only Te=0, t=-kappa/400, u=0, at cV*kappa²=150000/11.", flush=True)
    print("Second branch: no solution at positive scales. Nonflat endpoint equations remain open.", flush=True)
    print(f"\n{len(CHECKS)} checks, all True; runtime {time.time() - start:.1f} s", flush=True)


if __name__ == "__main__":
    main()
