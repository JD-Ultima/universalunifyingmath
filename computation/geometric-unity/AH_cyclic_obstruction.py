"""Exact certificate for AH5's single explicit cyclic five-form orbit.

The first Lambda5 c/wGamma pair has compact phase i. For WA,
Q(i*u)=-i*ALG.mixed('A',u,u)/6. Integer polarizations determine the
quadratic weight polynomials; a direct phased T call checks every relevant
output at AH4's regular point. No all-five-form invariant ansatz is used.
"""
from pathlib import Path
import sys
import time

import sympy as sp

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import T_singular_residual as T

X, ALG = T.X, T.ALG
CHECKS = []


def check(label, condition):
    ok = bool(condition)
    CHECKS.append(ok)
    print(f"[{len(CHECKS):02d}] {label}: {ok}", flush=True)
    assert ok, label


def grade(item, degree):
    return {mask: coefficient for mask, coefficient in item.items() if mask.bit_count() == degree}


def mask(indices):
    return sum(1 << index for index in indices)


def gamma_symmetric(field):
    converted = ALG.gamma_right([grade(item, 13) for item in field])
    matrix = sp.Matrix([[sp.expand(item.get(1 << j, 0)/sp.I) for j in range(14)] for item in converted])
    return (matrix+matrix.T).applyfunc(lambda value: sp.expand(value/2))


def main():
    start = time.time()
    print("AH5. The explicit AH4 cyclic orbit fails the first Lambda5 branch at nonzero w.", flush=True)
    block = {0, 1, 3, 6, 10}
    supports = [frozenset((index+shift) % 14 for index in block) for shift in range(14)]
    omega = {mask(support): 1 for support in supports}
    c = [grade(ALG.mul({1 << a: 1}, omega), 4) for a in range(14)]
    w_lift = [grade(ALG.mul({1 << a: 1}, omega), 6) for a in range(14)]
    d = ALG.gamma_right(w_lift)
    gram = sp.Matrix([[sum(left.get(key, 0)*right.get(key, 0) for key in left.keys() | right.keys())
                       for right in c] for left in c])
    check("the declared cyclic form has N=14 and contraction Gram 5I",
          len(omega) == 14 and sum(value*value for value in omega.values()) == 14 and gram == 5*sp.eye(14))
    check("both first-pair copies have the correct compact phase i",
          T.antihermitian(T.scale(c, sp.I)) and T.antihermitian(T.scale(d, sp.I)))

    # These are raw integer Clifford polarizations, before compact phases.
    cc, cd, dd = (ALG.mixed('A', c, c), ALG.mixed('A', c, d), ALG.mixed('A', d, d))
    r, s = sp.symbols('r s', real=True)
    indices = ((2, 3, 4, 5, 6), (1, 2, 4, 7, 9), (0, 1, 3, 4, 6, 8, 9, 10, 11))
    expected = (-sp.Rational(8, 3)*(7*r*r-27*s*s), -160*r*s,
                -sp.Rational(16, 3)*(4*r*r+17*s*s))
    polynomials = []
    for coordinate, target in zip(indices, expected):
        key = mask(coordinate)
        polynomial = sp.factor(-(cc[0].get(key, 0)*r*r+2*cd[0].get(key, 0)*r*s+dd[0].get(key, 0)*s*s)/6)
        polynomials.append(polynomial)
        check(f"Q/i at external index 0, Clifford indices {coordinate}", sp.expand(polynomial-target) == 0)
        print("  exact polynomial:", polynomial, flush=True)
    ideal = sp.groebner(polynomials, r, s, domain=sp.QQ)
    check("the three displayed coefficients generate (r^2,rs,s^2)", list(ideal) == [r*r, r*s, s*s])
    check("the grade-nine coefficient is strictly negative at every real nonzero chart",
          sp.Poly(-sp.Rational(3, 16)*polynomials[2], r, s) == sp.Poly(4*r*r+17*s*s, r, s))

    rr, ss = sp.sqrt(14), sp.Integer(2)
    control = T.scale(X.fplus(T.scale(c, rr), T.scale(d, ss)), sp.I)
    direct = T.quadratic(control, (1, 0, 0, 0))
    agreement = True
    for a in range(14):
        for key in set(direct[a]) | set(cc[a]) | set(cd[a]) | set(dd[a]):
            if key.bit_count() in (5, 9):
                prediction = -sp.I*(rr*rr*cc[a].get(key, 0)+2*rr*ss*cd[a].get(key, 0)+ss*ss*dd[a].get(key, 0))/6
                agreement &= sp.expand(direct[a].get(key, 0)-prediction) == 0
    check("direct phased T.quadratic agrees in every grade-five and grade-nine control slot", agreement)
    check("the actual AH4 control has the displayed nonzero coefficients",
          [sp.simplify(p.subs({r: rr, s: ss})) for p in polynomials]
          == [sp.Rational(80, 3), -320*sp.sqrt(14), -sp.Rational(1984, 3)])
    check("that same control still satisfies AH4's isotropic quadratic Gamma channel",
          gamma_symmetric(direct) == 480*sp.sqrt(14)*sp.eye(14))

    a, w, kappa = sp.Rational(1, 96), sp.sqrt(14)/48, 96*sp.sqrt(sp.Rational(15, 19))
    check("AH4's first-conic chart and scalar channel balance hold at the regular control",
          w*w == (1+96*a)*(1+sp.Rational(832, 3)*a)/1280
          and (1+96*a)*rr == 48*w*ss and (1+sp.Rational(832, 3)*a)*ss == sp.Rational(80, 3)*w*rr
          and sp.simplify(kappa*kappa*w*(1+208*a)-480*sp.sqrt(14)) == 0)

    for name, copy in (("c", c), ("wGamma", d)):
        derivatives = [T.derivative(T.scale(copy, sp.I), direction, (1, 0, 0, 0)) for direction in (0, 3, 7, 13)]
        check(name+": nonzero derivative controls contain no grade-five or grade-nine term",
              any(any(field) for field in derivatives)
              and all(key.bit_count() in (3, 7, 11) for field in derivatives for item in field for key in item))
    copies = ALG.forms(1, (0,))
    vector_derivatives = [T.derivative(copies[index], direction, (1, 0, 0, 0))
                          for index in (1, 3) for direction in (0, 1)]
    check("slaved-vector derivative controls have no grade-five or grade-nine output",
          all(key.bit_count() not in (5, 9) for field in vector_derivatives for item in field for key in item))

    print("The Gram channel forces C=0 for the cyclic scaling/rotation orbit at nonzero w.", flush=True)
    print("The remaining grade-nine polynomial cannot vanish for real (r,s)!=(0,0), independently of derivative choices.", flush=True)
    print("Scope: this explicit orbit, the first c/wGamma pair, nonzero w only; no all-five-form or endpoint classification.", flush=True)
    print(f"\n{len(CHECKS)} checks, all True; runtime {time.time()-start:.1f} s", flush=True)


if __name__ == '__main__':
    main()
