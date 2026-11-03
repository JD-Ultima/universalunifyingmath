"""Independent Clifford trace controls for AV's partial Gram reduction.

No repository or Clifford helper is imported. The finite forms are convention
checks, not AS-admissible data or a universal basis certificate. The AV note
proves universality by cyclic trace, reversion and the covariant product rule.
"""
from fractions import Fraction as Q
from time import perf_counter
import sympy as sp

N = 14
FULL = (1 << N)-1
CHECKS = []


def check(label, condition):
    assert bool(condition), label
    CHECKS.append(label)
    print(f'[{len(CHECKS):02d}] {label}: True')


def add(*values):
    out = {}
    for value in values:
        for key, val in value.items():
            out[key] = out.get(key, 0)+val
    return {key: val for key, val in out.items() if val}


def scale(value, coefficient):
    return {key: val*coefficient for key, val in value.items() if val*coefficient}


def mul(left, right):
    out = {}
    for a, av in left.items():
        for b, bv in right.items():
            crossings = sum((b & ((1 << i)-1)).bit_count() for i in range(N) if a >> i & 1)
            key = a ^ b
            out[key] = out.get(key, 0)+(-1 if crossings & 1 else 1)*av*bv
    return {key: val for key, val in out.items() if val}


def grade(value, degree):
    return {key: val for key, val in value.items() if key.bit_count() == degree}


def blade(indices, coefficient=1):
    return {sum(1 << i for i in indices): Q(coefficient)}


def ctr(i, value):
    return grade(mul({1 << i: 1}, value), 4)


def dot(left, right):
    return sum(val*right.get(key, 0) for key, val in left.items())


def gram_variation(field, omega, i, j, k):
    return dot(ctr(j, field[i]), ctr(k, omega))+dot(ctr(j, omega), ctr(k, field[i]))


def trace_controls():
    examples = (
        ('simple', blade(range(5))),
        ('disjoint and shared', add(blade(range(5)), blade(range(5, 10), 2), blade((0, 10, 11, 12, 13), -3))),
        ('epsilon mixed', add(blade(range(5)), blade((0, 5, 6, 7, 8), 2), blade((1, 10, 11, 12, 13), 3))),
    )
    for label, omega in examples:
        gammas = [{1 << i: 1} for i in range(N)]
        sandwiches = [mul(mul(omega, gamma), omega) for gamma in gammas]
        y5 = [grade(value, 5) for value in sandwiches]
        y9G = [mul(grade(value, 9), {FULL: 1}) for value in sandwiches]
        ordinary, axial = [], []
        for i in range(N):
            ordinary.append(dot(y5[i], omega))
            axial.append(dot(y9G[i], omega))
        good5 = good9 = symmetry = True
        t3, t3G = {}, {}
        for i in range(N):
            for j in range(N):
                middle = mul(mul(sandwiches[i], gammas[j]), omega)
                for k in range(N):
                    raw = mul(middle, gammas[k])
                    t3[i, j, k] = raw.get(0, 0)
                    t3G[i, j, k] = mul(raw, {FULL: 1}).get(0, 0)
                    good5 &= gram_variation(y5, omega, i, j, k) == (ordinary[i] if j == k else 0)+t3[i, j, k]
                    good9 &= gram_variation(y9G, omega, i, j, k) == (axial[i] if j == k else 0)-t3G[i, j, k]
        for key in t3:
            i, j, k = key
            symmetry &= (t3[key] == t3[j, i, k] == t3[k, j, i]
                         and t3G[key] == t3G[j, i, k] == t3G[k, j, i])
        check(label+' convention check: ordinary Gram identity in all 2744 slots', good5)
        check(label+' convention check: right-Gamma Gram identity in all 2744 slots', good9)
        check(label+' convention check: both cubic tensors are fully symmetric', symmetry)


def coefficient_controls():
    t, theta = sp.symbols('tau theta', real=True)
    x = t*t
    alpha = -3*(304+9*theta)/(16*(61+16*x))
    zeta = 3*t*(-112+5*theta)/(16*(64+19*x))
    Ac = -sp.Rational(4, 3)*t*(41+11*x)
    Bc = sp.Rational(8, 3)*(38+11*x)
    d = sp.factor(sp.Rational(8, 3)*(4+x)*(zeta-t*alpha)-48*t)
    leading = sp.factor(sp.Rational(13, 4)-sp.Rational(11, 12)*x-sp.Rational(7, 2)*t*Ac/Bc)
    check('projected compatibility leaves the positive symmetric cubic coefficient',
          sp.cancel(leading-(741+436*x+55*x*x)/(6*(38+11*x))) == 0)
    a_poly, b_poly = sp.Poly(sp.cancel(-Ac/t), t), sp.Poly(Bc, t)
    check('compatibility coefficient factors are positive in tau squared',
          all(value >= 0 for value in a_poly.all_coeffs()+b_poly.all_coeffs())
          and a_poly.TC() > 0 and b_poly.TC() > 0)
    a = (5*x-9)/(32*(66-25*x))
    w = 7*t/(16*(66-25*x))
    delta = sp.factor(1+224*a-3328*w*w)
    th = sp.factor(sp.Rational(3328, 3)*w*(x+4)/t-7*delta/(360*w*w))
    dd = sp.factor(d.subs(theta, th))
    z = sp.symbols('z', real=True)
    Pd = sp.Poly(1568750*z**5+565335*z**4-50705892*z**3-81149677*z**2
                 +138437232*z+46051632, z)
    expected_d = 16*Pd.as_expr().subs(z, x)/(315*t*(16*x+61)*(19*x+64)*(25*x-66))
    check('the exact hook coefficient has a nonzero finite-level numerator',
          sp.cancel(dd-expected_d) == 0 and Pd.degree() == 5 and Pd.TC() != 0)
    aa = sp.symbols('a', real=True)
    x_of_a = 9*(1+sp.Rational(704, 3)*aa)/(5*(1+160*aa))
    check('the canonical squared ratio has strictly positive a derivative',
          sp.factor(sp.diff(x_of_a, aa)-672/(5*(1+160*aa)**2)) == 0
          and sp.cancel(x_of_a.subs(aa, a)-x) == 0)
    print('  d(tau) =', dd)
    print('  Positive cubic coefficient =', sp.factor(leading))
    print('  Finite exceptional levels are retained by the density proof; no root search is needed.')


if __name__ == '__main__':
    start = perf_counter()
    print('AV: opposite-conic Gram transport, independent exact Clifford convention controls.')
    trace_controls()
    coefficient_controls()
    print('No control form is asserted to satisfy the full AS zero-order constraints.')
    print('Cyclic trace/reversion and the stated hand proof supply universality.')
    print('Zero-order redundancy of the residual symmetric cubic remains unresolved.')
    print(f'{len(CHECKS)} new controls passed; runtime {perf_counter()-start:.3f} s')
