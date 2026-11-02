"""Exact checks for T's residual/norm identity in the compact Euclidean real form.

Uses only the earlier independent Clifford implementation xcheck_round14.py;
does not import J, R, or the original field-equation computation. All phases
are applied explicitly: a real member has code coefficients -i*(x,y,z,w),
the invariant odd field is i*(t Phi1+u Phi1Gamma), and even monomials receive
the anti-Hermitian phase dictated by their Clifford adjoint.

These checks verify algebraic identities and specified jets. The proof and
the distinction between a residual constraint and a full solution are in
formalization/geometric-unity/T-singular-residual.md.
"""
from fractions import Fraction as F
from itertools import combinations
import sympy as sp

import xcheck_round14 as X


ALG = X.Algebra((1,) * X.N)


def scale(field, factor):
    return [X.scaled(value, factor) for value in field]


def real_even(field):
    """Turn each real Clifford monomial into an anti-Hermitian matrix."""
    return [{key: value * (sp.I if ((key.bit_count() * (key.bit_count()-1)//2) % 2 == 0) else 1)
             for key, value in item.items()} for item in field]


def antihermitian(field):
    return all(sp.expand(sp.conjugate(value) * (-1)**(key.bit_count()*(key.bit_count()-1)//2) + value) == 0
               for item in field for key, value in item.items())


def pair(a, b):
    return sp.expand(-ALG.pair(a, b))


def mix(a, b, member):
    """Actual real-form 2B(a,b), including the -i coefficient of each member."""
    out = [dict() for _ in range(X.N)]
    for index, coefficient in enumerate(member):
        if coefficient:
            factor = -sp.I * X.PHASE[index] * coefficient / 3
            out = X.fplus(out, scale(ALG.mixed(X.OPS[index], a, b), factor))
    return X.normalized(out)


def quadratic(a, member):
    return scale(mix(a, a, member), sp.Rational(1, 2))


def derivative(a, direction, member):
    zeta = {(b,c): X.plus(X.scaled(a[c], int(b == direction)),
                         X.scaled(a[b], -int(c == direction)))
            for b,c in combinations(range(X.N), 2)}
    out = [dict() for _ in range(X.N)]
    for index, coefficient in enumerate(member):
        if coefficient:
            out = X.fplus(out, scale(ALG.dual(X.OPS[index], zeta), -sp.I*X.PHASE[index]*coefficient))
    return X.normalized(out)


def invariant(t, u):
    return scale(X.fplus(scale(ALG.phi, t), scale(ALG.phig, u)), sp.I)


def baseline(t, u, kappa, g_diagonal):
    """WA baseline, dt=0: Y_odd=i G gamma + kappa To + Q(To)."""
    to = invariant(t, u)
    curvature = [{1 << a: sp.I*g_diagonal[a]} for a in range(X.N)]
    return X.normalized(X.fplus(X.fplus(curvature, scale(to, kappa)), quadratic(to, (1,0,0,0))))


def main():
    # An arbitrary invariant To, vector-sector even field, and non-vector field.
    # The kernel is not assumed in this first algebraic check.
    to = invariant(F(2,7), F(-3,11))
    vector_copies = ALG.forms(1, (5,))
    vector = [dict() for _ in range(X.N)]
    for coefficient, value in zip((1,2,-1,3), vector_copies):
        vector = X.fplus(vector, scale(real_even(value), coefficient))
    extra = real_even(ALG.forms(3, (0,4,11))[2])
    member = (1,2,3,5)
    assert all(antihermitian(value) for value in (to, vector, extra))
    cross = mix(vector, extra, member)
    assert pair(to, cross) == pair(extra, mix(to, vector, member)) == 0
    kappa = sp.Rational(7,5)
    lk = X.fplus(scale(extra, kappa), mix(extra, to, member))
    assert pair(extra, lk) == kappa*pair(extra, extra) + 2*pair(to, quadratic(extra, member))
    current_count = 0
    for direction in range(X.N):
        assert pair(to, derivative(extra, direction, member)) == 0
        current_count += 1
    print(f"T1: compact real phases, cubic polarization, cross-term orthogonality, and {current_count} zero-current contractions passed.")

    # R10 control: the previously solved flat Z1 background has zero scalar residual.
    control_g = [sp.Rational(11,20000)] * X.N  # 15/(2*cV), cV=150000/11.
    control_ev = baseline(F(-1,400), 0, 1, control_g)
    assert not any(X.normalized(control_ev))
    assert pair(invariant(F(-1,400), 0), control_ev) == 0
    print("R10 control: t=-1/400, u=0, cV=150000/11 gives <To,E_V>=0.")

    # WA's two Lambda7 singular lines: u can vary, while dt=0 kills its source.
    # G=diag(1,...,14) is an algebraic coefficient check, not an asserted observerse metric.
    g_diagonal = list(range(1, X.N+1))
    trace_g = sum(g_diagonal)
    k,u = sp.symbols('kappa u', real=True)
    for denominator in (176,208):
        a = -sp.Rational(1,denominator)
        scalar_over_kappa = a*sp.Symbol('traceG', nonnegative=True) + 14*(1+312*a)*(a*a*k*k+u*u)
        assert a < 0 and 1+312*a < 0
        print(f"T4 line t=-kappa/{denominator}: <To,E_V>/kappa = {sp.factor(scalar_over_kappa)}; strictly negative for kappa!=0, traceG>=0.")
        for kv in (1,-2):
            tv,uv = F(-kv,denominator),F(2,13)
            determinant = (kv*kv+224*kv*tv-3328*uv*uv)*(kv*kv+544*kv*tv-3328*uv*uv)
            assert determinant > 0
            to = invariant(tv,uv)
            # The c and w Lambda7 copies have roots -kappa/176 and -kappa/208.
            copy = ALG.forms(7,(0,2,4,6,8,10,12))[0 if denominator == 176 else 2]
            extra = real_even(copy)
            lk = X.fplus(scale(extra,kv),mix(extra,to,(1,0,0,0)))
            assert not any(X.normalized(lk))
            assert pair(extra,extra) > 0 and antihermitian(extra)
            assert not any(derivative(scale(ALG.phig,sp.I),1,(1,0,0,0)))
            ev = baseline(tv,uv,kv,g_diagonal)
            expected = tv*trace_g + 14*(kv+312*tv)*(tv*tv+uv*uv)
            assert pair(to,ev) == expected and sp.Rational(expected)/kv < 0
            dq = X.fplus(derivative(extra,3,(1,0,0,0)),quadratic(extra,(1,0,0,0)))
            assert pair(to,dq) == -kv*pair(extra,extra)/2
            full = X.fplus(ev,dq)
            assert pair(to,full) == expected-kv*pair(extra,extra)/2
    print("T2/T4: four real kernel cases (both kappa signs), WA du-source=0, scalar residual and full-minus-baseline pairing verified.")

    # A pure-hook addition cannot change the Gamma-partner scalar equation.
    partner = scale(ALG.phig,sp.I)
    hook_lines = []
    (kappa_symbol,sigma,theta,delta),blocks = X.templates(1)
    for grade in (2,4,6):
        for copy in ALG.hook(grade):
            extra = real_even(copy)
            assert not any(mix(extra,partner,(1,0,0,0)))
            assert pair(partner,quadratic(extra,(1,0,0,0))) == 0
            for direction in range(X.N):
                assert pair(partner,derivative(extra,direction,(1,0,0,0))) == 0
        for diagonal in blocks[f'H{grade}'].diagonal():
            slope = sp.diff(diagonal,sigma)
            ratio = -1/slope
            assert ratio > 0 and 1+208*ratio != 0
            seven = blocks['Lam7'].subs({sigma:ratio*kappa_symbol,theta:0,delta:0})
            assert sp.cancel(seven.det()/kappa_symbol**4) > 0
            hook_lines.append(ratio)
    tv,uv = F(3,136),F(2,13)
    ev = baseline(tv,uv,1,g_diagonal)
    assert pair(partner,ev) == 14*uv*(1+208*tv)
    symbolic_t,symbolic_u = sp.symbols('t u')
    for name in ('Lam3','Lam5'):
        determinant = sp.Poly(sp.expand(blocks[name].det().subs(
            {sigma:symbolic_t,theta:0,delta:-symbolic_u})),symbolic_u)
        assert determinant.degree() == 4 and determinant.LC() > 0
    assert len(set(hook_lines)) == 6
    print("T5: six hook copies have B(k,iPhi13)=0 and zero Gamma current; six distinct positive hook lines; other non-vector degeneracies occur at finitely many u values.")
    print("T4 excludes its two Lambda7 lines under the stated metric hypotheses. T5 forces u constant on the positive hook lines; constant-field cases remain open.")


if __name__ == '__main__':
    main()
