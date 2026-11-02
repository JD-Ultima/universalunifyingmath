"""Exact controls for the two single-Lambda3 branches on positive Z1.

F4/K, N1 and S4 supply the metric geometry; AA/AF/AG supply the endpoint and
simple-form reductions. This helper checks Clifford calibrations, retained
covariant derivative channels, and the exact radial multiplicity algebra.
It shares T's accepted Clifford engine and does not reconstruct the metric.
"""
from pathlib import Path
import sys

import sympy as sp

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import T_singular_residual as T

X, ALG = T.X, T.ALG
MEMBER = (1, 0, 0, 0)
CHECKS = []
# Label, odd ratio, lift, compact phase, Q1 eigenvalues (multiplicities 3/11).
CASES = (
    ("c", -sp.Rational(3, 112), "c", 1, sp.Integer(4), sp.Integer(4)),
    ("wGamma", -sp.Rational(1, 400), "wG", 1, sp.Rational(220, 3), sp.Integer(180)),
    ("i*w", -sp.Rational(1, 240), "w", sp.I, sp.Integer(220), sp.Integer(60)),
    ("i*cGamma", -sp.Rational(3, 272), "cG", sp.I, sp.Rational(4, 3), sp.Integer(12)),
)


def check(label, condition):
    ok = bool(condition)
    CHECKS.append(ok)
    print(f"[{len(CHECKS):02d}] {label}: {ok}", flush=True)
    assert ok, label


def grade(value, degree):
    return {mask: coefficient for mask, coefficient in value.items()
            if mask.bit_count() == degree}


def lift(form, kind, phase):
    degree = next(iter(form)).bit_count() if form else 3
    target = degree - 1 if kind.startswith("c") else degree + 1
    field = [grade(ALG.mul({1 << a: 1}, form), target) for a in range(X.N)]
    if kind.endswith("G"):
        field = ALG.gamma_right(field)
    return T.scale(field, phase)


def grade_one(field):
    return sp.Matrix([[sp.expand(item.get(1 << b, 0) / sp.I)
                       for b in range(X.N)] for item in field])


def norm_gram(form):
    contraction = lift(form, "c", 1)
    norm = sp.expand(sum(value**2 for value in form.values()))
    gram = sp.Matrix([[sum(left.get(mask, 0)*right.get(mask, 0)
                           for mask in left.keys() | right.keys())
                       for right in contraction] for left in contraction])
    return norm, gram


def spectral_polynomial(H, parallel, perpendicular):
    D = 3*parallel + 11*perpendicular
    tau = sp.trace(H)
    return (H*H - (parallel+perpendicular)*tau*H/D
            + parallel*perpendicular*tau**2*sp.eye(X.N)/D**2).applyfunc(sp.expand)


def main():
    print("AK: complete single-Lambda3 classification on positive Z1")
    print(f"SymPy {sp.__version__}; real compact (14,0); imported metric lemmas.")
    kappa = sp.symbols("kappa", real=True, nonzero=True)
    amplitude = sp.symbols("amplitude", real=True)
    coordinate = {(1 << 0) | (1 << 1) | (1 << 2): 1}
    rotated = grade(ALG.mul(ALG.mul({1 << 0: 1, 1 << 4: 1},
                                    {1 << 1: 1, 1 << 5: 2}),
                            {1 << 2: 1, 1 << 6: -1}), 3)
    for label, a, kind, phase, parallel, perpendicular in CASES:
        for name, form in (("coordinate", coordinate), ("non-coordinate", rotated)):
            field = T.scale(lift(form, kind, phase), amplitude)
            norm, gram = norm_gram(form)
            expected = amplitude**2*(perpendicular*norm*sp.eye(X.N)
                                     + (parallel-perpendicular)*gram)
            mass = X.fplus(T.scale(field, kappa),
                          T.mix(field, T.invariant(kappa*a, 0), MEMBER))
            check(f"{label}, {name}: compact kernel and full grade-one stress",
                  T.antihermitian(field) and not any(X.normalized(mass))
                  and gram*gram == norm*gram and gram.rank() == 3
                  and (grade_one(T.quadratic(field, MEMBER))-expected).applyfunc(sp.expand)
                  == sp.zeros(X.N))
        degrees = {"c": {1, 5}, "wG": {9, 13}, "w": {3, 7}, "cG": {11}}[kind]
        jets = [T.derivative(lift(coordinate, kind, phase), direction, MEMBER)
                for direction in range(X.N)]
        check(f"{label}: all fourteen derivative directions retain the declared grades",
              all(all(mask.bit_count() in degrees for item in field for mask in item)
                  for field in jets)
              and all(grade_one(field)+grade_one(field).T == sp.zeros(X.N)
                      for field in jets))
        b = sp.factor(a*(1+312*a))
        baseline = T.baseline(kappa*a, 0, kappa, [0]*X.N)
        check(f"{label}: complete constant-odd baseline coefficient b/kappa^2={b}",
              (grade_one(baseline)-b*kappa**2*sp.eye(X.N)).applyfunc(sp.expand)
              == sp.zeros(X.N)
              and all(mask.bit_count() == 1 for item in baseline for mask in item))

    z, delta = sp.symbols("z delta", real=True)
    for label, _, _, _, parallel, perpendicular in CASES[1:]:
        H = sp.diag(*([-parallel*z]*3 + [-perpendicular*z]*11))
        D = 3*parallel + 11*perpendicular
        coefficient = -33*(parallel-perpendicular)**2/D**2
        check(f"{label}: spectral elimination includes z=0 and scalar limits",
              spectral_polynomial(H, parallel, perpendicular) == sp.zeros(X.N)
              and spectral_polynomial(delta*sp.eye(X.N), parallel, perpendicular)
              == coefficient*delta**2*sp.eye(X.N) and coefficient < 0)
        check(f"{label}: both 5/9 eigenvalue assignments violate the trace",
              5*parallel+9*perpendicular != D
              and 9*parallel+5*perpendicular != D)

    A, B = sp.symbols("A B", real=True)
    D = 3*A + 11*B
    check("universal multiplicity differences are exactly 2(A-B) and 6(A-B)",
          sp.expand(5*A+9*B-D) == 2*(A-B)
          and sp.expand(9*A+5*B-D) == 6*(A-B))
    # v=1-4lambda>0 parametrizes the full positive vertical-metric range.
    v, cV = sp.symbols("v cV", positive=True)
    c = 4*cV*v
    ric_h, ric_v = -1/c, -sp.Rational(5, 4)/cV
    scalar_flat = 5*ric_h+9*ric_v
    p, q = sp.simplify(ric_h-scalar_flat/2), sp.simplify(ric_v-scalar_flat/2)
    check("general-lambda radial Einstein limits have positive 5/9 eigenvalues",
          sp.simplify(p-(3+45*v)/(8*cV*v)) == 0
          and sp.simplify(q-(5+35*v)/(8*cV*v)) == 0
          and p.is_positive is True and q.is_positive is True)
    check("radial eigenvalues coincide exactly at lambda=1/5",
          sp.solve(p-q, v) == [sp.Rational(1, 5)]
          and (1-sp.Rational(1, 5))/4 == sp.Rational(1, 5)
          and sp.simplify(p.subs(v, sp.Rational(1, 5))-sp.Rational(15, 2)/cV) == 0)
    norm2 = sp.symbols("norm2", nonnegative=True)
    positive_b = sp.Rational(309, 1568)*kappa**2
    check("first-left symmetric equation contradicts the negative Euler Ricci eigenvalue",
          (positive_b+4*norm2+6/c).is_positive is True)
    check("opposite endpoint b values are strictly positive for either kappa sign",
          all((a*(1+312*a)*kappa**2).is_positive is True for _, a, *_ in CASES[2:]))
    scale = sp.Rational(150000, 11)/kappa**2
    solved = T.baseline(-kappa/400, 0, kappa,
                        [sp.Rational(15, 2)/scale]*X.N)
    check("surviving constant odd field has zero complete odd residual at the exact scale",
          scale.is_positive is True and not any(X.normalized(solved)))
    print(f"{len(CHECKS)} exact controls passed.")
    print("First branch: flat base, lambda=1/5, cV*kappa^2=150000/11, t=-kappa/400, u=0, Te=0.")
    print("Opposite branch: no solution. Other kernels, operators and field sectors remain open.")


if __name__ == "__main__":
    main()
