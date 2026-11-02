"""Exact nonlinear H5 obstruction on the Lambda3 c/wGamma branch (AA6).

Use real compact (14,0) fields, the real WA member and Tr(1)=1. Let
 k=r*c(omega)+s*w(omega)Gamma, with omega a real 3-form and r,s real.
Z_a=(interior_a omega) wedge omega and F4=sum_a(interior_a omega)^2_wedge.
The full grade identities are
 Q5=(8i/3)[(43s^2-3r^2)Z+(4r^2-49s^2)w(F4)],
 Q9*Gamma=(32i/3)rs[5Z-7w(F4)].
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

On the full single-kernel branch, D(k) has no H5 by AA4: its grade5 and
converted grade9 parts are pure exterior lifts of 4-forms. The slaved
baseline terms have input types scalar or V tensor V and curvature grade1,
so no H5. The vector/kernel mixed term has input type V tensor Lambda3,
which is Lambda4 + Lambda2 + H3, again no H5. These statements use parallel
equivariant lifts and remain valid when weights and omega vary. Hence a
full solution requires H5(Z)=0, since the two projected coefficients have
no common zero with real (r,s)!=(0,0). This is only a necessary pointwise
condition; it does not characterize decomposable forms or solve the PDE.

Exact checks below use the independent Clifford engine through T and the
reviewed derivative helper. The proof and limits are in AA6 of
formalization/geometric-unity/AA-conic-traces.md.
"""
import pathlib
import sys
import time

import sympy as sp

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import T_singular_residual as T
import AA_exterior_derivative as D

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
    return X.fplus(T.scale(lift(omega, 3, False), r),
                   T.scale(ALG.gamma_right(lift(omega, 3, True)), s))


def main():
    start = time.time()
    print('AA6: nonlinear Lambda3 c/wGamma hook projection, compact (14,0), real WA', flush=True)
    print(f'sympy {sp.__version__}; Q9 conversion is Q9 right-multiplied by Gamma.', flush=True)
    examples = (
        ('decomposable control', form((0, 1, 2)), False),
        ('two disjoint triples', X.plus(form((0, 1, 2)), form((3, 4, 5))), True),
        ('explicit excluded form e012+e034', X.plus(form((0, 1, 2)), form((0, 3, 4))), True),
        ('shared-index rank-four two-form factor', X.plus(form((0, 1, 2)), form((0, 3, 4), 2)), True),
        ('three disjoint triples', X.plus(X.plus(form((0, 1, 2)), form((3, 4, 5), 2)), form((6, 7, 8), -1)), True),
        ('overlapping triples', X.plus(X.plus(form((0, 1, 2)), form((0, 3, 4), 2)), form((0, 5, 6), -3)), True),
        ('full support in fourteen dimensions', X.plus(X.plus(X.plus(X.plus(form((0, 1, 2)), form((3, 4, 5), 2)),
             form((6, 7, 8), -1)), form((9, 10, 11), 3)), form((0, 12, 13), -2)), True),
    )
    alpha = sp.Rational(8, 3)*(43*S*S-3*R*R)
    beta = sp.Rational(8, 3)*(4*R*R-49*S*S)
    gamma, eta = sp.Rational(160, 3)*R*S, -sp.Rational(224, 3)*R*S
    for name, omega, nonzero_hook in examples:
        z, f4 = candidate(omega)
        wz = lift(f4, 4, True)
        hz = hook(z)
        wz_trace, contraction_trace = traces(z)
        check(name+': normalized hook and expected zero/nonzero status',
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
    explicit = X.plus(form((0, 1, 2)), form((0, 3, 4)))
    explicit_z, explicit_f4 = candidate(explicit)
    check('explicit excluded form has H5(Z)_5=-e12345/5',
          explicit_f4 == form((1, 2, 3, 4), 2)
          and hook(explicit_z)[5] == form((1, 2, 3, 4, 5), -sp.Rational(1, 5)))
    omega = examples[-1][1]
    vector = [dict() for _ in range(X.N)]
    for index, direction in enumerate((0, 4, 8, 12)):
        copy = ALG.forms(1, (direction,))[index]
        vector = X.fplus(vector, T.scale(T.real_even(copy), index+1))
    check('unrestricted vector/kernel mixed term has zero H5 in both grades',
          all(not any(hook(value)) for value in projections(T.mix(vector, field(omega), MEMBER))))
    # The diagonal curvature coefficients are algebraic control data, not an
    # assertion that an observer metric realizing this diagonal has been found.
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
    dk = X.fplus(D.actual(weighted[0], 3, False), D.actual(weighted[1], 3, True))
    check('covariant derivative control retains weight derivatives and has zero H5',
          all(not any(hook(value)) for value in projections(dk)))
    check('two hook coefficients have no common nonzero real weight pair',
          alpha.subs(R, 0) == sp.Rational(344, 3)*S*S
          and alpha.subs(S, 0) == -8*R*R and gamma == sp.Rational(160, 3)*R*S)
    print('Q5 = (8i/3)[(43s^2-3r^2)Z+(4r^2-49s^2)w(F4)].', flush=True)
    print('Q9*Gamma = (32i/3)rs[5Z-7w(F4)]; right conversion sign fixed as displayed.', flush=True)
    print('Full single-kernel solutions therefore require H5(Z)=0 pointwise; no decomposability or PDE sufficiency claim.', flush=True)
    print(f'\n{len(CHECKS)} checks, all True; runtime {time.time()-start:.1f} s', flush=True)


if __name__ == '__main__':
    main()
