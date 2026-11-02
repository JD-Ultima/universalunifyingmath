"""Exact controls for AH's necessary Gram constraints on both Lambda5 conics.

For real compact WA, each pair has
 Sym(Q(k)_13*Gamma/i)=(160/3)rs*N*I-(304/3)rs*S,
where N=|omega5|^2 and S_ab=<i_a omega5,i_b omega5>.
The contraction proof has only norm*I and Gram patterns: ten input-form
indices plus two output indices lie below dimension14, excluding epsilon;
metric contractions within one alternating form vanish. A complete symbolic
unit-simple calibration fixes both constants, not a sampled parameter range.

The two real phases are i on c/wGamma (grades4/8) and 1 on w/cGamma
(grades6/10). Derivative and mixed terms have input V tensor Lambda5,
containing no Sym2V, so their symmetric converted output is zero, including
variable weights. The remaining vector/invariant terms are AF's accepted
formula. No five-form hook classification or full PDE solution is inferred.

At w!=0, the full symmetric equation requires N>0 and S=alpha*I-b*C*C^T,
b=(624/19)w/(rs). Gram positivity and rank>=5 imply alpha>0 and rank13/14.
On left conic branches b<0, so rank14; on right branches b>0, rank13
requires C to be the unique annihilator direction. This excludes support
dimension<=12, not every five-form. The cyclic isotropic control below
exhibits the remaining algebraic gap without an expensive full Clifford run.
"""
import pathlib
import sys
import time
from itertools import combinations

import sympy as sp

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import T_singular_residual as T

X, ALG = T.X, T.ALG
MEMBER = (1, 0, 0, 0)
R, S = sp.symbols('r s', real=True)
CHECKS = []
DATA = (('c/wGamma', 96, sp.Rational(832, 3), 48, sp.Rational(80, 3)),
        ('w/cGamma', sp.Rational(704, 3), 160, sp.Rational(80, 3), 48))


def check(label, condition):
    ok = bool(condition)
    CHECKS.append(ok)
    print(f'[{len(CHECKS):02d}] {label}: {ok}', flush=True)
    assert ok, label


def grade(value, degree):
    return {mask: coefficient for mask, coefficient in value.items() if mask.bit_count() == degree}


def form(indices, coefficient=1):
    return {sum(1 << index for index in indices): coefficient}


def lifts(omega):
    c = [grade(ALG.mul({1 << a: 1}, omega), 4) for a in range(X.N)]
    w = [grade(ALG.mul({1 << a: 1}, omega), 6) for a in range(X.N)]
    return c, w


def kernel(omega, pair, r=R, s=S):
    c, w = lifts(omega)
    if pair == 0:
        return T.scale(X.fplus(T.scale(c, r), T.scale(ALG.gamma_right(w), s)), sp.I)
    return X.fplus(T.scale(w, r), T.scale(ALG.gamma_right(c), s))


def stress(omega):
    c, _ = lifts(omega)
    norm = sum(coefficient**2 for coefficient in omega.values())
    gram = sp.Matrix([[sum(left.get(mask, 0)*right.get(mask, 0) for mask in left.keys() | right.keys())
                       for right in c] for left in c])
    return norm, gram


def symmetric_gamma(value):
    converted = ALG.gamma_right([grade(item, 13) for item in value])
    matrix = sp.Matrix([[sp.expand(item.get(1 << j, 0)/sp.I) for j in range(X.N)] for item in converted])
    return (matrix+matrix.T).applyfunc(lambda coefficient: sp.expand(coefficient/2))


def main():
    start = time.time()
    print('AH: necessary Lambda5 Gram spectrum, real compact WA', flush=True)
    simple = form((0, 1, 2, 3, 4))
    examples = (
        ('unit simple', simple),
        ('two disjoint five-forms', X.plus(simple, form((5, 6, 7, 8, 9)))),
        ('overlap with off-diagonal Gram', X.plus(X.plus(simple, form((0, 1, 2, 3, 5), 2)), form((4, 6, 7, 8, 9), -1))),
        ('full support', X.plus(X.plus(simple, form((5, 6, 7, 8, 9), 2)), form((0, 10, 11, 12, 13), -1))),
    )
    vector = [dict() for _ in range(X.N)]
    for index, direction in enumerate((0, 4, 8, 12)):
        vector = X.fplus(vector, T.scale(T.real_even(ALG.forms(1, (direction,))[index]), index+1))
    a, w = sp.symbols('a w', real=True)
    for pair, (name, d, e, p, q) in enumerate(DATA):
        for label, omega in examples:
            norm, gram = stress(omega)
            k = kernel(omega, pair)
            actual = symmetric_gamma(T.quadratic(k, MEMBER))
            expected = sp.Rational(160, 3)*R*S*norm*sp.eye(X.N)-sp.Rational(304, 3)*R*S*gram
            check(name+', '+label+': exact symbolic Gram formula and compact phase',
                  T.antihermitian(k) and (actual-expected).applyfunc(sp.expand) == sp.zeros(X.N)
                  and gram.trace() == 5*norm and gram.rank() >= 5)
        check(name+': nondecomposable control has a genuinely off-diagonal Gram entry',
              stress(examples[2][1])[1][4, 5] != 0)
        omega = examples[-1][1]
        check(name+': vector/kernel mixed term has no symmetric converted component',
              symmetric_gamma(T.mix(vector, kernel(omega, pair), MEMBER)) == sp.zeros(X.N))
        dk = X.fplus(T.derivative(kernel(simple, pair, 2, 3), 0, MEMBER),
                     T.derivative(kernel(omega, pair, -1, 4), 7, MEMBER))
        check(name+': derivative and independent weight-jet control has no such component',
              any(dk) and symmetric_gamma(dk) == sp.zeros(X.N))

        f = sp.expand((1+d*a)*(1+e*a)/(p*q))
        for tag, rr, ss in (('first', p*w, 1+d*a), ('alternate', 1+e*a, q*w)):
            chart = kernel(simple, pair, rr, ss)
            lk = X.fplus(chart, T.mix(chart, T.invariant(a, w), MEMBER))
            check(name+', '+tag+': exact real kernel chart modulo w squared minus f',
                  all(sp.expand(sp.rem(sp.expand(value), w*w-f, w)) == 0 for item in lk for value in item.values()))
        left, right = sorted((-1/sp.Rational(d), -1/sp.Rational(e)))
        check(name+': branch signs and negative scalar-candidate value',
              p > 0 and d > 0 and e > 0 and left < -sp.Rational(1, 208) < right
              and f.subs(a, -sp.Rational(1, 208)) < 0)
        print('  '+name+': left endpoint', left, 'right endpoint', right,
              '; f(-1/208)=', f.subs(a, -sp.Rational(1, 208)), flush=True)

    block = {0, 1, 3, 6, 10}
    supports = [frozenset((index+shift) % 14 for index in block) for shift in range(14)]
    check('cyclic limitation: fourteen distinct supports with overlap at most two',
          len(set(supports)) == 14 and max(len(left & right) for left, right in combinations(supports, 2)) == 2)
    omega = {}
    for indices in supports:
        omega = X.plus(omega, form(indices))
    norm, gram = stress(omega)
    check('cyclic limitation: exact full-support N=14 and S=5I', norm == 14 and gram == 5*sp.eye(14))

    # This is algebraic channel data only, not an observer metric or PDE
    # solution. The nonzero-a point also avoids the intersection at a=0.
    aa, ww = sp.Rational(1, 96), sp.sqrt(14)/48
    rr, ss, kk = sp.sqrt(14), sp.Integer(2), 96*sp.sqrt(sp.Rational(15, 19))
    f1 = (1+96*aa)*(1+sp.Rational(832, 3)*aa)/1280
    f2 = (1+sp.Rational(704, 3)*aa)*(1+160*aa)/1280
    check('cyclic limitation: regular first-conic chart, separate from the other Lambda5 factor',
          ww**2 == f1 and f1 != f2 and rr == 48*ww and ss == 1+96*aa
          and 1+224*aa-3328*f1 == -sp.Rational(152, 9)
          and 1+544*aa-3328*f1 == -sp.Rational(122, 9))
    (kap, sigma, theta, delta), blocks = X.templates(1)
    substitutions = {kap: 1, sigma: aa, theta: 0, delta: -ww}
    check('cyclic limitation: all other non-vector factors are regular',
          all(sp.simplify(blocks[name].subs(substitutions).det()) != 0
              for name in ('Lam3', 'Lam7', 'H2', 'H4', 'H6'))
          and sp.simplify(blocks['Lam5'].extract((1, 3), (1, 3)).subs(substitutions).det()) != 0)
    scalar = -kk*kk*ww*(1+208*aa)+sp.Rational(160, 3)*rr*ss*norm
    check('cyclic limitation: the full symmetric Gamma channel balances with C=0',
          sp.simplify(scalar-sp.Rational(304, 3)*rr*ss*5) == 0)
    print('Necessary reduction: S=alpha I-b C*C^T, b=(624/19)w/(rs), trace S=5N.', flush=True)
    print('N>0, rank S=13 or 14. Left branch: rank 14. Right: rank 13 only with C as annihilator.', flush=True)
    print('The isotropic full-support form passes this selected channel; no full Lambda5 exclusion or PDE solution follows.', flush=True)
    print(f'\n{len(CHECKS)} checks, all True; runtime {time.time()-start:.1f} s', flush=True)


if __name__ == '__main__':
    main()
