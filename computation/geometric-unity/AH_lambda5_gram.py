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

AH6 adds all four constant-odd zero-u endpoint copies. Their full grade-one
quadratics have the universal form q*N*I+ell*S. N1's Euler Ricci identity
then gives exact lower norm bounds and requires the Euler direction to be
a Gram eigenvector. Actual real-form Grams attain the selected algebraic
bounds; no observer metric or full field solution is constructed. An isolated
zero-u point on a varying conic is not a constant-odd endpoint domain.

AH7 excludes a constant normalized Gram spectrum on any nonempty relative-open
fibre patch of such a domain. Homogeneous Cayley-Hamilton continues only a
metric identity. The positive self-adjoint zero-trace barrier and connectedness
preserve the actual negative trace sign, contradicting F4's positive radial
trace limit. No field or inequality is analytically continued.
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
# name, pair, pure chart weights, endpoint, q, ell, b/kappa^2, norm-bound denominator
ENDPOINTS = (
    ('i*c5', 0, 1, 0, -sp.Rational(1, 96), sp.Rational(40, 3), sp.Rational(32, 3), sp.Rational(3, 128), 24),
    ('i*w5Gamma', 0, 0, 1, -sp.Rational(3, 832), sp.Integer(112), -sp.Integer(64), sp.Rational(3, 6656), 48),
    ('w5', 1, 1, 0, -sp.Rational(3, 704), sp.Rational(112, 3), sp.Rational(320, 3), sp.Rational(87, 61952), 336),
    ('c5Gamma', 1, 0, 1, -sp.Rational(1, 160), sp.Integer(40), -sp.Integer(32), sp.Rational(19, 3200), 40),
)


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


def grade_one(value):
    return sp.Matrix([[sp.expand(item.get(1 << j, 0)/sp.I) for j in range(X.N)] for item in value])


def first_pair_derivative(omega, direction, contraction_copy):
    """AA4 at m=5; a coordinate jet includes arbitrary absorbed-weight jets."""
    product = ALG.mul({1 << direction: 1}, omega)
    exterior, codifferential = grade(product, 6), X.scaled(grade(product, 4), -1)
    wedge = [grade(ALG.mul({1 << a: 1}, exterior), 7) for a in range(X.N)]
    contraction = [grade(ALG.mul({1 << a: 1}, codifferential), 3) for a in range(X.N)]
    left, right = (-10, 6) if contraction_copy else (-14, 18)
    result = X.fplus(T.scale(wedge, left), T.scale(contraction, right))
    return result if contraction_copy else ALG.gamma_right(result)


def endpoint_controls(examples, cyclic):
    """AH6: exact constants and channel sharpness, not metric attainability."""
    kappa = sp.symbols('kappa', real=True, nonzero=True)
    norm, b, inv_c, mu, q, ell = sp.symbols('N b inv_c mu q ell', real=True)
    ric_e = b/6+(2*q+5*ell)*norm/12-ell*mu
    residual = 2*b+12*inv_c+(2*q+5*ell)*(norm-mu)+(2*q-7*ell)*mu
    check('AH6: trace reversal and the full Euler-component residual agree',
          sp.expand(12*(ric_e+inv_c)-residual) == 0)
    off_diagonal = examples[2][1]
    samples = (('zero', {}),)+examples+(('cyclic isotropic full support', cyclic),)
    cartan_first, cartan_second = form((1, 2, 3, 4, 5)), form((0, 2, 3, 4, 5))
    for name, pair, r, s, endpoint, qq, ll, bcoef, denominator in ENDPOINTS:
        d, e = DATA[pair][1:3]
        A, B = 1+d*endpoint, 1+e*endpoint
        variables, blocks = X.templates(1)
        substitutions = dict(zip(variables, (1, endpoint, 0, 0)))
        check('AH6 '+name+': pure endpoint chart, det M and every other block regular',
              A*r == 0 and B*s == 0 and (r != 0 or s != 0)
              and 1+224*endpoint != 0 and 1+544*endpoint != 0
              and blocks['Lam5'].subs(substitutions).rank() == 3
              and all(matrix.subs(substitutions).rank() == matrix.cols
                      for key, matrix in blocks.items() if key != 'Lam5'))
        baseline = T.baseline(kappa*endpoint, 0, kappa, [0]*X.N)
        check('AH6 '+name+': positive baseline bI for either nonzero real kappa sign',
              endpoint*(1+312*endpoint) == bcoef and bcoef > 0
              and grade_one(baseline) == bcoef*kappa*kappa*sp.eye(X.N)
              and all(mask.bit_count() == 1 for item in baseline for mask in item))
        stress_ok = True
        for _, omega in samples:
            nn, gram = stress(omega)
            field = kernel(omega, pair, r, s)
            mass = X.fplus(T.scale(field, kappa), T.mix(field, T.invariant(kappa*endpoint, 0), MEMBER))
            stress_ok &= (T.antihermitian(field) and not any(X.normalized(mass))
                          and gram.trace() == 5*nn
                          and grade_one(T.quadratic(field, MEMBER)) == qq*nn*sp.eye(X.N)+ll*gram)
        check('AH6 '+name+': full all-slot stress on zero/simple/nondecomposable/off-diagonal/full-support controls', stress_ok)
        allowed = ({3, 7} if r else {7, 11}) if pair == 0 else ({5, 9} if r else {9, 13})
        derivatives = [T.derivative(kernel(off_diagonal, pair, r, s), j, MEMBER) for j in range(X.N)]
        check('AH6 '+name+': all fourteen arbitrary-weight jet directions retain only the stated non-one grades',
              any(any(value) for value in derivatives)
              and all(mask.bit_count() in allowed for value in derivatives for item in value for mask in item)
              and all(grade_one(value) == sp.zeros(X.N) for value in derivatives))
        if pair == 0:
            check('AH6 '+name+': complete first-pair AA4 covariant derivative formula',
                  all(X.normalized(value) == X.normalized(first_pair_derivative(off_diagonal, j, bool(r)))
                      for j, value in enumerate(derivatives)))
        else:
            cartan = X.fplus(T.derivative(kernel(cartan_first, pair, r, s), 0, MEMBER),
                              T.derivative(kernel(cartan_second, pair, r, s), 1, MEMBER))
            check('AH6 '+name+': nonzero Cartan jet remains in the other odd equations', any(cartan))
        c0, c1 = 2*qq+5*ll, 2*qq-7*ll
        required = sp.cancel((b+6*inv_c+c0*norm/2)/(6*ll))
        bound_residual = sp.cancel(required-norm) if ll > 0 else required
        expected = ((b+6*inv_c-denominator*norm)/(6*ll) if ll > 0
                    else (denominator*norm-b-6*inv_c)/(-6*ll))
        check('AH6 '+name+': exact Euler eigenvalue and sharp lower norm bound with mixed coefficients',
              min(c0, c1) < 0 < max(c0, c1)
              and sp.expand(ric_e.subs({q: qq, ell: ll, mu: required})+inv_c) == 0
              and sp.cancel(bound_residual-expected) == 0)
        # c=kappa^2=1. These tensors meet the selected algebraic equations;
        # they are not asserted to be curvatures of reconstructed Z1 metrics.
        L0 = bcoef+6
        sharp_norm = L0/denominator
        sharp_omega = form((0, 1, 2, 3, 4), sp.sqrt(sharp_norm))
        nn, gram = stress(sharp_omega)
        direction = 0 if ll > 0 else 5
        unit = sp.eye(X.N)[:, direction]
        synthetic_G = -(bcoef+qq*nn)*sp.eye(X.N)-ll*gram
        synthetic_Ric = synthetic_G-sp.trace(synthetic_G)*sp.eye(X.N)/12
        check('AH6 '+name+': bound attained by a real-form Gram satisfying the full selected Euler vector equation',
              sharp_norm > 0 and nn == sharp_norm and synthetic_Ric*unit == -unit)
        nn = 2*sharp_norm
        mm = sp.cancel(required.subs({norm: nn, b: bcoef, inv_c: 1}))
        omega = X.plus(form((0, 1, 2, 3, 4), sp.sqrt(mm)),
                       form((5, 6, 7, 8, 9), sp.sqrt(nn-mm)))
        actual_norm, gram = stress(omega)
        synthetic_G = -(bcoef+qq*actual_norm)*sp.eye(X.N)-ll*gram
        synthetic_Ric = synthetic_G-sp.trace(synthetic_G)*sp.eye(X.N)/12
        check('AH6 '+name+': a nondecomposable two-support Gram also meets that vector equation',
              0 < mm < nn and actual_norm == nn and synthetic_Ric[:, 0] == -sp.eye(X.N)[:, 0])
        check('AH6 '+name+': isotropic nonzero Gram is incompatible with positive Euler b,c', qq+5*ll/14 > 0)
        print(f'  AH6 {name}: a={endpoint}, b/kappa^2={bcoef}, q={qq}, ell={ll}; '
              f'N >= (b+6/c)/{denominator}; derivative grades {sorted(allowed)}.', flush=True)
    print('AH6 bounds are sharp only for the selected algebraic Gram/Euler channel; no Z1 metric or full solution is constructed.', flush=True)
    print('AH6 requires a constant-odd endpoint domain; isolated zero-u points on a moving conic may retain Dv.', flush=True)


def homogeneous_characteristic(matrix, trace, coefficient_sum, coefficients):
    """Evaluate sum c_j*(D H)^(14-j)*tau^j without symbolic eigenvalues."""
    identity = sp.eye(matrix.rows)
    result = coefficients[0]*identity
    for degree, coefficient in enumerate(coefficients[1:], 1):
        result = coefficient_sum*matrix*result+coefficient*trace**degree*identity
    return result


def spectrum_controls(examples):
    """AH7's algebraic normalizations; F4/S4/N1 remain geometric inputs."""
    omega = examples[2][1]
    norm, gram = stress(omega)
    scaled_norm, scaled_gram = stress(X.scaled(omega, -3))
    shape = gram/norm
    check('AH7: a nondecomposable normalized Gram shape is unchanged by nonzero scaling',
          scaled_norm == 9*norm and scaled_gram/scaled_norm == shape
          and shape.trace() == 5)
    simple_norm, simple_gram = stress(examples[0][1])
    disjoint_norm, disjoint_gram = stress(examples[1][1])
    variable = sp.symbols('z')
    check('AH7: genuinely different real five-form shapes have different normalized spectra',
          (simple_gram/simple_norm).charpoly(variable).as_expr()
          != (disjoint_gram/disjoint_norm).charpoly(variable).as_expr())
    zero_trace_matrix = sp.diag(1, -1, *([0]*12))
    for name, _, _, _, _, qq, ll, _, _ in ENDPOINTS:
        A = qq*sp.eye(X.N)+ll*shape
        D = 14*qq+5*ll
        H = -norm*A
        trace = -D*norm
        coefficients = A.charpoly(variable).all_coeffs()
        check('AH7 '+name+': full homogeneous characteristic identity with repeated spectral roots',
              min(qq, qq+ll) > 0 and D > 0 and A.trace() == D
              and H.trace() == trace
              and homogeneous_characteristic(H, trace, D, coefficients) == sp.zeros(X.N))
        check('AH7 '+name+': at zero trace only the nonzero self-adjoint fourteenth-power barrier remains',
              homogeneous_characteristic(zero_trace_matrix, 0, D, coefficients)
              == D**14*zero_trace_matrix**14 != sp.zeros(X.N))
    b, c, cV, cH, radial = sp.symbols('b c cV cH radial', positive=True)
    base_scalar = sp.symbols('R_g', real=True)
    ricci_at_zero = -b*sp.eye(X.N)-(-14*b)*sp.eye(X.N)/12
    check('AH7: H=0 would give the positive Ricci scalar b/6, incompatible with N1',
          ricci_at_zero == b*sp.eye(X.N)/6 and b/6+1/c > 0)
    scalar = base_scalar/(radial*cH)-5/c-sp.Rational(45, 4)/cV
    trace = sp.expand(-6*scalar+14*b)
    limit = sp.limit(trace, radial, sp.oo)
    check('AH7: exact F4 radial trace and its strictly positive general-lambda limit',
          sp.expand(trace-(14*b+30/c+sp.Rational(135, 2)/cV-6*base_scalar/(radial*cH))) == 0
          and limit == 14*b+30/c+sp.Rational(135, 2)/cV and limit.is_positive)
    print('AH7: no nonempty fibre-open patch can have constant normalized Gram spectrum on a constant-odd endpoint domain.', flush=True)
    print('Only the metric polynomial continues; connectedness and the zero-trace barrier carry the sign.', flush=True)
    print('Fixed SO14/scaling families are excluded; general varying shapes and the full PDE remain unresolved.', flush=True)


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
    endpoint_controls(examples, omega)
    spectrum_controls(examples)
    print(f'\n{len(CHECKS)} checks, all True; runtime {time.time()-start:.1f} s', flush=True)


if __name__ == '__main__':
    main()
