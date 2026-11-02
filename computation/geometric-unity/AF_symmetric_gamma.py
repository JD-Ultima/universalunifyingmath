"""Exact controls for AF's symmetric grade-thirteen equation.

Compact real (14,0), fixed WA, invariant odd field, kappa!=0 and det M!=0.
The proof uses AA7's classification of the Lambda3 c/wGamma hook zero locus.
SO(14) equivariance makes the complete symbolic simple-form and aligned-vector
calibrations universal on their declared orbits. The rank argument is a hand
proof; no range of parameters or metrics is inferred from finite controls.

Right multiplication by Gamma converts grade thirteen to grade one; Gamma²=-1.
All coefficients below use the anti-Hermitian compact phases from T.
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
R, S, C, B = sp.symbols("r s c b", real=True)


def check(label, condition):
    ok = bool(condition)
    CHECKS.append(ok)
    print(f"[{len(CHECKS):02d}] {label}: {ok}", flush=True)
    assert ok, label


def grade(value, degree):
    return {mask: coefficient for mask, coefficient in value.items() if mask.bit_count() == degree}


def gamma_symmetric(field):
    """Symmetric grade13*Gamma coefficient matrix divided by i."""
    converted = ALG.gamma_right([grade(item, 13) for item in field])
    M = sp.Matrix([[sp.expand(item.get(1 << b, 0) / sp.I) for b in range(X.N)]
                   for item in converted])
    return (M + M.T).applyfunc(lambda value: sp.expand(value / 2))


def lift(omega):
    contraction = [grade(ALG.mul({1 << a: 1}, omega), 2) for a in range(X.N)]
    exterior = [grade(ALG.mul({1 << a: 1}, omega), 4) for a in range(X.N)]
    return contraction, ALG.gamma_right(exterior)


def kernel(omega, r=R, s=S):
    contraction, exterior_gamma = lift(omega)
    return X.fplus(T.scale(contraction, r), T.scale(exterior_gamma, s))


def form(indices, coefficient=1):
    return {sum(1 << index for index in indices): coefficient}


def stress(omega):
    contraction, _ = lift(omega)
    norm = sum(value * value for value in omega.values())
    gram = sp.Matrix([[sum(left.get(mask, 0) * right.get(mask, 0)
                           for mask in left.keys() | right.keys())
                       for right in contraction] for left in contraction])
    return norm, gram


def vector(direction, c, b):
    p2, p14 = ALG.forms(1, (direction,))[1::2]
    return X.fplus(T.scale(p2, c), T.scale(p14, b))


def main():
    start = time.time()
    print("AF: symmetric right-Gamma equation on the single Lambda3 c/wGamma conic", flush=True)
    print(f"SymPy {sp.__version__}; real compact phases, no observer-metric samples.", flush=True)
    omega = form((0, 1, 2))
    k = kernel(omega)
    Q = gamma_symmetric(T.quadratic(k, MEMBER))
    expected = sp.diag(*([-sp.Rational(88, 3) * R * S] * 3 + [40 * R * S] * 11))
    check("complete symbolic simple-form quadratic has the stated 3+11 eigenvalues",
          T.antihermitian(k) and Q == expected)

    v = vector(0, C, B)
    expected = sp.diag(*([0] + [16 * B * C] * 13))
    check("complete symbolic aligned-vector quadratic has coefficient 16bc",
          T.antihermitian(v) and gamma_symmetric(T.quadratic(v, MEMBER)) == expected)

    # These two vector-jet types determine the symmetric equivariant map:
    # it has form alpha*Sym(H)+beta*tr(H)*I.
    copies = ALG.forms(1, (0,))[1::2]
    check("both vector derivatives have zero symmetric converted component",
          all(gamma_symmetric(T.derivative(copy, direction, MEMBER)) == sp.zeros(X.N)
              for copy in copies for direction in (0, 1)))

    # A simple form with a non-coordinate three-plane. Pairwise orthogonal
    # factors are e0+e4, e1+2e5 and e2-e6; norm²=2*5*2.
    rotated = grade(ALG.mul(ALG.mul({1 << 0: 1, 1 << 4: 1},
                                    {1 << 1: 1, 1 << 5: 2}),
                            {1 << 2: 1, 1 << 6: -1}), 3)
    norm, gram = stress(rotated)
    expected = 40 * R * S * norm * sp.eye(X.N) - sp.Rational(208, 3) * R * S * gram
    check("non-coordinate simple form obeys the full stress formula",
          norm == 20 and gram.rank() == 3 and gram * gram == norm * gram
          and gamma_symmetric(T.quadratic(kernel(rotated), MEMBER)) == expected)

    vk = X.fplus(vector(0, 2, 6), vector(5, -1, -3))
    check("vector/kernel mixed control has no symmetric converted component",
          gamma_symmetric(T.mix(vk, kernel(rotated), MEMBER)) == sp.zeros(X.N))
    derivative = X.fplus(T.derivative(kernel(omega, 2, 3), 0, MEMBER),
                         T.derivative(kernel(rotated, -1, 4), 7, MEMBER))
    check("kernel derivative control, including independent weight jets, has no such component",
          gamma_symmetric(derivative) == sp.zeros(X.N))

    t, u, kappa = sp.symbols("t u kappa", real=True)
    baseline = T.baseline(t, u, kappa, list(range(1, X.N + 1)))
    check("invariant baseline has scalar -u*(kappa+208t)",
          (gamma_symmetric(baseline) + u * (kappa + 208 * t) * sp.eye(X.N)).applyfunc(sp.expand) == sp.zeros(X.N))

    a, w, c = sp.symbols("a w c", real=True)
    qv = gamma_symmetric(T.quadratic(vector(0, c, 208 * w * c), MEMBER))
    expected = 3328 * w * c * c * sp.diag(*([0] + [1] * 13))
    check("slaved aligned-vector specialization has coefficient 3328w",
          qv == expected)

    f = sp.expand((1 + sp.Rational(112, 3) * a) * (1 + 400 * a) / sp.Rational(11264, 3))
    check("the scalar Gamma candidate lies outside the real conic",
          f.subs(a, -sp.Rational(1, 208)) == -sp.Rational(3, 14872))
    check("the only zero-w endpoints are -3/112 and -1/400",
          sp.solve(f, a) == [-sp.Rational(3, 112), -sp.Rational(1, 400)])

    print("Symmetric equation: [-kappa²w(1+208a)+3328w|C|²+40nu]I ="
          "3328w C*C^T +(208/3)nu P3, nu=rs|omega|².", flush=True)
    print("AA7 plus rank14 versus rank<=4 rules out every w!=0 point on this single branch.", flush=True)
    print("Only two constant-odd endpoints remain; their even-field/PDE equations are not solved.", flush=True)
    print(f"\n{len(CHECKS)} checks, all True; runtime {time.time() - start:.1f} s", flush=True)


if __name__ == "__main__":
    main()
