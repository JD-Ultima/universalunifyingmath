"""Exact controls for AG's other Lambda3 w/cGamma branch, compact phase i.

Use real compact (14,0) fields, the real WA member and Tr(1)=1. Let
 k=i[r*w(omega)+s*c(omega)Gamma], with omega a real 3-form and r,s real.
Z_a=(interior_a omega) wedge omega and F4=sum_a(interior_a omega)^2_wedge.
The full grade identities are
 Q5=(8i/3)(3s^2-59r^2)Z-8i*r^2*w(F4),
 Q9*Gamma=(224i/3)rs*Z.
The conversion is RIGHT multiplication Q9*Gamma, with Gamma^2=-1; it is
not the inverse conversion -Q9*Gamma. The hook projection on V tensor Lambda5
subtracts c(sum e^a wedge Z_a)/6 and w(sum interior_a Z_a)/10.

Universality by invariant contractions. The input consists of two alternating
3-tensors and the output has one free vector slot and five alternating slots.
After simplifying the paired Gamma factors, Clifford normal ordering uses
metric contractions. There are twelve tensor slots in total, below dimension14,
so an epsilon invariant could not contribute either. Metric pairings cannot contract two indices in
one input 3-form, or pair two of the five output alternating slots. If no
input indices contract each other, the unique surviving pattern is Z. If
one cross-input pair contracts, the free output slot must pair with one
alternating output slot, giving w(F4). Two such cross-input contractions
would force a forbidden output-output pairing. Thus the full quadratic map
has only these two patterns. After hook projection only H5(Z) remains.
Right-Gamma conversion has the same tensor types; two Gamma factors leave
no additional orientation-dependent contraction. Calibration on two forms
where Z and w(F4) are independent fixes all symbolic r,s coefficients.

AA4 is not assumed on these copies. The derivative of input grade4 lies
in grades3/7, while input grade12 gives grade11. Thus D(k) has no grade5,
grade9 or grade13 component, including all variable-weight derivatives.
Other derivative components, including Cartan components, are retained. The slaved
baseline terms have input types scalar or V tensor V and curvature grade1,
so no H5. The vector/kernel mixed term has input type V tensor Lambda3,
which is Lambda4 + Lambda2 + H3, again no H5. These statements use parallel
equivariant lifts and remain valid when weights and omega vary. Hence a
full solution requires H5(Z)=0, since the two projected coefficients have
no common zero with real (r,s)!=(0,0). This is only a necessary pointwise
condition. AA7 classifies its zero locus as zero or simple, allowing AF's
symmetric rank argument to be applied to this other conic. No endpoint PDE
sufficiency follows.

Exact checks below use the independent Clifford engine through T. No AA4
derivative formula is used for these other copies. The universal proof and
scope are in formalization/geometric-unity/AG-other-lambda3.md.
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
R, S = sp.symbols('r s', real=True)


def check(label, value):
    ok = bool(value)
    CHECKS.append(ok)
    print(f'[{len(CHECKS):02d}] {label}: {ok}', flush=True)
    assert ok, label


def grade(value, degree):
    return {key: coefficient for key, coefficient in value.items() if key.bit_count() == degree}


def form(indices, coefficient=1):
    return {sum(1 << index for index in indices): coefficient}


def lift(value, degree, exterior):
    target = degree+1 if exterior else degree-1
    return [grade(ALG.mul({1 << a: 1}, value), target) for a in range(X.N)]


def same(left, right):
    return not any(X.normalized(X.fplus(left, right, -1)))


def traces(field):
    wedge, contraction = {}, {}
    for a, value in enumerate(field):
        product = ALG.mul({1 << a: 1}, value)
        wedge = X.plus(wedge, grade(product, 6))
        contraction = X.plus(contraction, grade(product, 4))
    return wedge, contraction


def hook(field):
    wedge, contraction = traces(field)
    return X.normalized(X.fplus(X.fplus(field, T.scale(lift(wedge, 6, False), -sp.Rational(1, 6))),
                                T.scale(lift(contraction, 4, True), -sp.Rational(1, 10))))


def projections(field):
    q5 = [grade(value, 5) for value in field]
    q9_gamma = ALG.gamma_right([grade(value, 9) for value in field])
    return q5, q9_gamma


def candidate(omega):
    contraction = lift(omega, 3, False)
    z = [grade(ALG.mul(value, omega), 5) for value in contraction]
    f4 = {}
    for value in contraction:
        f4 = X.plus(f4, grade(ALG.mul(value, value), 4))
    return z, f4


def field(omega, r=R, s=S):
    return T.scale(X.fplus(T.scale(lift(omega, 3, True), r),
                           T.scale(ALG.gamma_right(lift(omega, 3, False)), s)), sp.I)


def gamma_symmetric(value):
    converted = ALG.gamma_right([grade(item, 13) for item in value])
    matrix = sp.Matrix([[sp.expand(item.get(1 << j, 0)/sp.I) for j in range(X.N)] for item in converted])
    return (matrix+matrix.T).applyfunc(lambda coefficient: sp.expand(coefficient/2))


def derivative_jet(alpha, beta):
    out = [dict() for _ in range(X.N)]
    for direction in range(X.N):
        value = T.scale(X.fplus(lift(alpha[direction], 3, True),
                                ALG.gamma_right(lift(beta[direction], 3, False))), sp.I)
        out = X.fplus(out, T.derivative(value, direction, MEMBER))
    return out


def main():
    start = time.time()
    print('AG: other Lambda3 w/cGamma branch, compact phase i, real WA', flush=True)
    print(f'sympy {sp.__version__}; Q9 conversion is Q9 right-multiplied by Gamma.', flush=True)
    examples = (
        ('decomposable control', form((0, 1, 2)), False),
        ('two disjoint triples', X.plus(form((0, 1, 2)), form((3, 4, 5))), True),
        ('shared-index rank-four two-form factor', X.plus(form((0, 1, 2)), form((0, 3, 4), 2)), True),
        ('three disjoint triples', X.plus(X.plus(form((0, 1, 2)), form((3, 4, 5), 2)), form((6, 7, 8), -1)), True),
        ('overlapping triples', X.plus(X.plus(form((0, 1, 2)), form((0, 3, 4), 2)), form((0, 5, 6), -3)), True),
        ('full support in fourteen dimensions', X.plus(X.plus(X.plus(X.plus(form((0, 1, 2)), form((3, 4, 5), 2)),
             form((6, 7, 8), -1)), form((9, 10, 11), 3)), form((0, 12, 13), -2)), True),
    )
    alpha = sp.Rational(8, 3)*(3*S*S-59*R*R)
    beta = -8*R*R
    gamma, eta = sp.Rational(224, 3)*R*S, 0
    for name, omega, nonzero_hook in examples:
        z, f4 = candidate(omega)
        wz = lift(f4, 4, True)
        hz = hook(z)
        wz_trace, contraction_trace = traces(z)
        check(name+': hook projector and expected zero/nonzero status',
              not wz_trace and contraction_trace == f4 and same(hz, hook(hz))
              and not any(traces(hz))
              and ((T.pair(T.scale(hz, sp.I), T.scale(hz, sp.I)) > 0) if nonzero_hook else not any(hz)))
        k = field(omega)
        q5, q9_gamma = projections(T.quadratic(k, MEMBER))
        expected5 = T.scale(X.fplus(T.scale(z, alpha), T.scale(wz, beta)), sp.I)
        expected9 = T.scale(X.fplus(T.scale(z, gamma), T.scale(wz, eta)), sp.I)
        check(name+': both full grade identities for symbolic r,s',
              T.antihermitian(k) and same(q5, expected5) and same(q9_gamma, expected9))
        check(name+': both projected coefficients for symbolic r,s',
              same(hook(q5), T.scale(hz, sp.I*alpha))
              and same(hook(q9_gamma), T.scale(hz, sp.I*gamma)))

    # A mixed nondecomposable form and an unrestricted real vector-isotype
    # field. The representation proof above, rather than this control,
    # supplies the universal absence of H5 in mixed/baseline terms.
    omega = examples[-1][1]
    vector = [dict() for _ in range(X.N)]
    for index, direction in enumerate((0, 4, 8, 12)):
        copy = ALG.forms(1, (direction,))[index]
        vector = X.fplus(vector, T.scale(T.real_even(copy), index+1))
    check('unrestricted vector/kernel mixed term has zero H5 in both grades',
          all(not any(hook(value)) for value in projections(T.mix(vector, field(omega), MEMBER))))
    # These diagonal curvature coefficients are algebraic convention data,
    # not a constructed observer metric.
    baseline = X.fplus(T.baseline(sp.Rational(2, 7), -sp.Rational(3, 11), 5, list(range(1, 15))),
                       T.quadratic(vector, MEMBER))
    for direction in (1, 5, 9):
        baseline = X.fplus(baseline, T.derivative(vector, direction, MEMBER))
    check('baseline control has no H5 in either grade',
          all(not any(hook(value)) for value in projections(baseline)))

    jet = [{} for _ in range(X.N)]
    jet[0] = examples[2][1]
    jet[7] = examples[1][1]
    weighted = []
    for value, derivative in ((2, {0: 5}), (3, {4: 7})):
        weighted.append([X.plus(X.scaled(v, value), X.scaled(omega, derivative.get(j, 0))) for j, v in enumerate(jet)])
    dk = derivative_jet(*weighted)
    check('covariant derivative retains weight jets but has no grade5/9/13 output',
          any(dk) and all(not any(grade(item, degree) for item in dk) for degree in (5, 9, 13)))
    check('two hook coefficients have no common nonzero real weight pair',
          alpha.subs(R, 0) == 8*S*S
          and alpha.subs(S, 0) == -sp.Rational(472, 3)*R*R and gamma == sp.Rational(224, 3)*R*S)

    simple_q = gamma_symmetric(T.quadratic(field(form((0, 1, 2))), MEMBER))
    expected = sp.diag(*([-sp.Rational(88, 3)*R*S]*3+[40*R*S]*11))
    check('simple-form symmetric Gamma tensor has the same 3+11 eigenvalues as AF', simple_q == expected)
    rotated = grade(ALG.mul(ALG.mul({1 << 0: 1, 1 << 4: 1}, {1 << 1: 1, 1 << 5: 2}),
                            {1 << 2: 1, 1 << 6: -1}), 3)
    contraction = lift(rotated, 3, False)
    norm = sum(coefficient**2 for coefficient in rotated.values())
    gram = sp.Matrix([[sum(left.get(mask, 0)*right.get(mask, 0) for mask in left.keys() | right.keys())
                       for right in contraction] for left in contraction])
    check('non-coordinate simple form has the universal Gamma stress tensor',
          norm == 20 and gram.rank() == 3 and gram*gram == norm*gram
          and gamma_symmetric(T.quadratic(field(rotated), MEMBER))
          == 40*R*S*norm*sp.eye(X.N)-sp.Rational(208, 3)*R*S*gram)
    twistor = [{} for _ in range(X.N)]
    twistor[0], twistor[1] = form((1, 2, 3)), form((0, 2, 3))
    empty = [{} for _ in range(X.N)]
    check('Cartan derivative controls survive on both copies, without the three projected grades',
          all(any(value) and all(not any(grade(item, degree) for item in value) for degree in (5, 9, 13))
              for value in (derivative_jet(twistor, empty), derivative_jet(empty, twistor))))
    check('vector/kernel mixed control has zero symmetric Gamma tensor',
          gamma_symmetric(T.mix(vector, field(omega), MEMBER)) == sp.zeros(X.N))
    a = sp.symbols('a', real=True)
    f = (1+240*a)*(1+sp.Rational(272, 3)*a)/sp.Rational(11264, 3)
    w = sp.symbols('w', real=True)
    for name, rr, ss in (('first', 32*w, 1+240*a),
                         ('alternate', 1+sp.Rational(272, 3)*a, sp.Rational(352, 3)*w)):
        chart = field(form((0, 1, 2)), rr, ss)
        lk = X.fplus(chart, T.mix(chart, T.invariant(a, w), MEMBER))
        check(name+' compact-phase kernel chart is exact modulo w squared minus f',
              all(sp.expand(sp.rem(sp.expand(value), w*w-f, w)) == 0 for item in lk for value in item.values()))
    check('the Gamma scalar candidate is outside this real conic', f.subs(a, -sp.Rational(1, 208)) == -sp.Rational(1, 43264))
    check('the only zero-w endpoints are -3/272 and -1/240', sp.solve(f, a) == [-sp.Rational(3, 272), -sp.Rational(1, 240)])
    print('Q5 = (8i/3)(3s^2-59r^2)Z-8i*r^2*w(F4).', flush=True)
    print('Q9*Gamma = (224i/3)rs*Z; right conversion sign fixed as displayed.', flush=True)
    print('AA7 and the symmetric Gamma rank argument exclude w!=0 on this single branch; endpoint equations remain open.', flush=True)
    print(f'\n{len(CHECKS)} checks, all True; runtime {time.time()-start:.1f} s', flush=True)


if __name__ == '__main__':
    main()
