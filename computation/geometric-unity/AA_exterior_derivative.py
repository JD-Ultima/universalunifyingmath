"""Derivative of the c/wGamma exterior kernels supporting AA.

For p=3,5, eps_p=1,i respectively, D=*^-1 WA d0, and real p-forms alpha,beta:
 D(eps_p*c(alpha)) = i*eps_p*[2p*w(d alpha)-2(p-2)*c(delta alpha)],
 D(eps_p*w(beta)Gamma) = i*eps_p*[2(12-p)*w(d beta)-2(14-p)*c(delta beta)]Gamma.
Here delta=-sum_j interior(e_j)*nabla_j and Gamma is parallel, Gamma^2=-1.
Maps c and w are the Clifford contraction and exterior lifts from R1. These
formulas concern the c/wGamma conics only. They do not solve Q(k) or the PDE.

Proof. The anticommutator of a grade-3 element with c(omega) has grades
p-2,p+2 only. For w(omega)Gamma, move Gamma right to turn the anticommutator
into a commutator; its pre-Gamma grades are again p-2,p+2. The output modules
therefore contain no copy of the Cartan summand H_p in T* tensor Lambda^p,
for p=3,5. Equivariance makes D factor through d and delta. On the d-only jet
e0 tensor (e1 wedge ... wedge ep), p terms for c and 14-p-2 terms for wGamma
each contribute twice the exterior lift. On the delta-only jet e0 tensor
(e0 wedge ... wedge e_(p-1)), the corresponding counts are p-2 and 14-p,
with a minus sign from delta=-interior*nabla. Other components vanish by the
same Clifford parity rule. This gives the four stated coefficients.

The maps and Gamma are Levi-Civita parallel, so this proof holds in a normal
frame for arbitrary covariant jets. For variable weights use alpha=A*omega,
beta=B*omega; the derivatives of A and B must be included. On a nonzero-u
conic arc they can be written as functions of a. A kernel chart is
(A,B)=(p_coupling*w,1+d*a), with the alternative (1+e*a,q_coupling*w)
where the first pair vanishes. Here w=u/kappa and the coupling constants are
those in AA's table. At u=0 use a smooth conic chart rather than presuming
w is a smooth function of a.

For fixed real weights A,B the principal-symbol norm is
 4(13-p)[p^2 A^2+(12-p)^2 B^2] |xi wedge omega|^2
 +4(p-1)[(p-2)^2 A^2+(14-p)^2 B^2] |interior(xi)omega|^2.
The exterior/contraction lifts have norms 14-r and r on r-forms. Distinct
grades are orthogonal. In the p=5 middle-grade term, right Gamma is an
orthogonal complex structure and hence skew, so its cross term is zero.
Both coefficients are positive when A,B are not both zero, and the two
form norms sum to |xi|^2 |omega|^2. This proves injectivity for real xi!=0;
it does not establish ellipticity or solvability of the coupled full system.

This helper supplies exact full-output checks with the independent
xcheck_round14 product through T's compact helpers. No J/R field-equation
evaluator is imported. The identities are integrated as AA4-AA5 and internally
reviewed in reports/reviews/REVIEW23-GU-AA-DERIV-2026-11-02.md.
"""
import pathlib
import sys
import time

import sympy as sp

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import T_singular_residual as T

X, ALG = T.X, T.ALG
MEMBER = (1, 0, 0, 0)
CHECKS = []


def check(label, condition):
    ok = bool(condition)
    CHECKS.append(ok)
    print(f'[{len(CHECKS):02d}] {label}: {ok}', flush=True)
    assert ok, label


def homogeneous(value, grade):
    return {key: coefficient for key, coefficient in value.items() if key.bit_count() == grade}


def form(indices, coefficient=1):
    return {sum(1 << j for j in indices): coefficient}


def lift(value, grade, kind):
    target_grade = grade-1 if kind == 'c' else grade+1
    return [homogeneous(ALG.mul({1 << j: 1}, value), target_grade) for j in range(X.N)]


def exterior_and_codifferential(jet, p):
    exterior, codifferential = {}, {}
    for direction, value in enumerate(jet):
        product = ALG.mul({1 << direction: 1}, value)
        exterior = X.plus(exterior, homogeneous(product, p+1))
        codifferential = X.plus(codifferential, homogeneous(product, p-1), -1)
    return exterior, codifferential


def actual(jet, p, partner):
    out = [dict() for _ in range(X.N)]
    eps = 1 if p == 3 else sp.I
    for direction, value in enumerate(jet):
        if not value:
            continue
        field = lift(value, p, 'w' if partner else 'c')
        if partner:
            field = ALG.gamma_right(field)
        out = X.fplus(out, T.derivative(T.scale(field, eps), direction, MEMBER))
    return X.normalized(out)


def predicted(jet, p, partner):
    exterior, codifferential = exterior_and_codifferential(jet, p)
    a, b = ((2*(12-p), -2*(14-p)) if partner else (2*p, -2*(p-2)))
    out = X.fplus(T.scale(lift(exterior, p+1, 'w'), a),
                  T.scale(lift(codifferential, p-1, 'c'), b))
    if partner:
        out = ALG.gamma_right(out)
    return T.scale(out, sp.I*(1 if p == 3 else sp.I))


def same(left, right):
    return not any(X.normalized(X.fplus(left, right, -1)))


def form_norm(value):
    return sp.expand(sum(coefficient**2 for coefficient in value.values()))


def main():
    start = time.time()
    print('AA supporting calculation: exterior derivative of the two c/wGamma kernels', flush=True)
    print(f'sympy {sp.__version__}; compact (14,0), real WA, delta=-interior*nabla', flush=True)
    for p in (3, 5):
        cases = 0
        for indices in (tuple(range(p)), tuple(range(1, 2*p, 2))):
            omega = form(indices)
            for direction in range(X.N):
                jet = [{} for _ in range(X.N)]
                jet[direction] = omega
                for partner in (False, True):
                    assert same(actual(jet, p, partner), predicted(jet, p, partner))
                    cases += 1
        check(f'p={p}: all {cases} coordinate-jet full-output identities', cases == 56)

        jet = [{} for _ in range(X.N)]
        jet[0] = X.plus(form(tuple(range(1, p+1)), sp.Rational(2, 3)),
                        form(tuple(range(p)), -sp.Rational(3, 5)))
        jet[4] = form(tuple(range(2, p+2)), sp.Rational(7, 11))
        jet[9] = form(tuple(range(3, p+3)), -sp.Rational(5, 13))
        check(f'p={p}: both identities on a mixed covariant jet',
              all(same(actual(jet, p, partner), predicted(jet, p, partner)) for partner in (False, True)))

        # A symmetric mixed jet with zero exterior derivative and codifferential.
        twistor = [{} for _ in range(X.N)]
        twistor[0] = form((1,)+tuple(range(2, p+1)))
        twistor[1] = form((0,)+tuple(range(2, p+1)))
        ed, co = exterior_and_codifferential(twistor, p)
        check(f'p={p}: a nonzero Cartan jet is annihilated on both declared copies',
              not ed and not co and all(not any(actual(twistor, p, partner)) for partner in (False, True)))
        excluded = []
        for partner in (False, True):
            out = [dict() for _ in range(X.N)]
            for direction, value in enumerate(twistor):
                field = lift(value, p, 'c' if partner else 'w')
                if partner:
                    field = ALG.gamma_right(field)
                out = X.fplus(out, T.derivative(T.real_even(field), direction, MEMBER))
            excluded.append(any(X.normalized(out)))
        check(f'p={p}: the same Cartan jet survives on both excluded w/cGamma copies', all(excluded))

        # The variable weights are arbitrary; no conic equation is needed for
        # the differential product rule. This includes their gradient terms.
        a = sp.symbols('a', real=True)
        weight_c, weight_g = 1+a*a, 2*a-1
        omega = form(tuple(range(p)))
        gradient = {0: sp.Rational(2, 3), 6: -sp.Rational(3, 5)}
        weighted = []
        for weight in (weight_c, weight_g):
            weighted.append([X.plus(X.scaled(value, weight),
                                    X.scaled(omega, sp.diff(weight, a)*gradient.get(j, 0)))
                             for j, value in enumerate(jet)])
        actual_weighted = X.fplus(actual(weighted[0], p, False), actual(weighted[1], p, True))
        predicted_weighted = X.fplus(predicted(weighted[0], p, False), predicted(weighted[1], p, True))
        check(f'p={p}: variable-weight derivative includes both chain-rule terms',
              same(actual_weighted, predicted_weighted) and T.antihermitian(actual_weighted))

        # Principal-symbol norm: the middle-grade right-Gamma action for p=5
        # is orthogonal and skew, so its cross term is zero, not a cancellation
        # assumption about two independent grade-7 outputs.
        cc, gg = sp.symbols('C G', real=True)
        omega = X.plus(form(tuple(range(p)), sp.Rational(2, 3)),
                       form(tuple(range(1, p+1)), sp.Rational(3, 5)))
        symbol_jet = [{} for _ in range(X.N)]
        symbol_jet[0] = omega
        out = X.fplus(T.scale(actual(symbol_jet, p, False), cc),
                      T.scale(actual(symbol_jet, p, True), gg))
        ed, co = exterior_and_codifferential(symbol_jet, p)
        norm = (4*(13-p)*(p*p*cc*cc+(12-p)**2*gg*gg)*form_norm(ed)
                +4*(p-1)*((p-2)**2*cc*cc+(14-p)**2*gg*gg)*form_norm(co))
        check(f'p={p}: positive principal-symbol norm identity', sp.expand(T.pair(out, out)-norm) == 0)
        check(f'p={p}: compact phases and wedge/contraction norm decomposition',
              T.antihermitian(T.scale(lift(omega, p, 'c'), 1 if p == 3 else sp.I))
              and T.antihermitian(T.scale(ALG.gamma_right(lift(omega, p, 'w')), 1 if p == 3 else sp.I))
              and form_norm(ed)+form_norm(co) == form_norm(omega))
        print(f'  p={p}: D(c) coefficients d={2*p}, delta={-2*(p-2)}; '
              f'D(wGamma) coefficients d={2*(12-p)}, delta={-2*(14-p)}.', flush=True)
    print('For alpha=A*omega and beta=B*omega use d(A*omega)=A*domega+dA wedge omega and delta(A*omega)=A*deltaomega-interior(grad A)omega.', flush=True)
    print('The derivative symbol in omega is injective for real nonzero covectors when real A,B are not both zero.', flush=True)
    print('This concerns the two c/wGamma conics only; Q(k) and the untraced full odd PDE remain unsolved.', flush=True)
    print(f'\n{len(CHECKS)} checks, all True; runtime {time.time()-start:.1f} s', flush=True)


if __name__ == '__main__':
    main()
