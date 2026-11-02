"""Exact positive-hook intersection catalogue and a grade-five obstruction.

The six lines and four factor families are the rational formulas of R2/R6.
They are checked against the earlier independently verified block templates.
Selected kernel and quadratic checks use only T's compact-real-form helpers
and the independent xcheck_round14 Clifford engine, not R/J's evaluator.
No field-equation solution or general hook classification is inferred.
"""
import pathlib
import sys
import sympy as sp

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
import T_singular_residual as T

X,ALG = T.X,T.ALG
R = sp.Rational
LINES = ((R(3,136),'H2G',896),(R(1,40),'H4G',11648),(R(3,104),'H6G',36608),
         (R(3,88),'H6',36608),(R(1,24),'H4',11648),(R(3,56),'H2',896))
# name, exterior degree, copy indices, diagonal slopes, upper/lower u coefficients.
PAIRS = (('L3c',3,(0,2),R(112,3),400,R(352,3),32),
         ('L3w',3,(1,3),240,R(272,3),32,R(352,3)),
         ('L5c',5,(0,2),96,R(832,3),48,R(80,3)),
         ('L5w',5,(1,3),R(704,3),160,R(80,3),48))


def squared_substitution(expression, variable, square):
    polynomial = sp.Poly(sp.expand(expression),variable)
    assert all(power[0] % 2 == 0 for power in polynomial.monoms())
    return sp.factor(sum(coefficient*square**(power[0]//2)
                         for power,coefficient in polynomial.terms()))


def catalogue():
    (k,sigma,theta,delta),blocks = X.templates(1)
    u = sp.Symbol('u')
    records = []
    print('a=t/kappa | hook kernel | r3c | r3w | r5c | r5w; r=(u/kappa)^2')
    for a,hook,dimension in LINES:
        roots = []
        norms = []
        m = sp.factor(blocks['Lam1'].det().subs({k:1,sigma:a,theta:0,delta:-u}))
        seven = blocks['Lam7'].subs({k:1,sigma:a,theta:0,delta:-u})
        assert seven.det() != 0
        zero_hook = 0
        for name in ('H2','H4','H6'):
            hb = blocks[name].subs({k:1,sigma:a,theta:0,delta:-u})
            zero_hook += len(hb.nullspace())
        assert zero_hook == 1
        for name,degree,copies,d,e,p,q in PAIRS:
            aa,bb = 1+d*a,1+e*a
            root = sp.factor(aa*bb/(p*q))
            assert root > 0
            determinant = blocks[f'Lam{degree}'].det().subs({k:1,sigma:a,theta:0,delta:-u})
            assert squared_substitution(determinant,u,root) == 0
            assert squared_substitution(m,u,root) != 0
            # The singular two-copy block has rank one, not zero: aa is positive.
            assert aa > 0 and bb > 0
            roots.append(root)
            # Required Lambda norm from the Gamma trace, if a full solution exists.
            required = sp.factor(14*(1+208*a)*(aa+bb)/(p*q))
            g1,g2 = (degree,14-degree) if copies == (0,2) else (14-degree,degree)
            assert g1*p == g2*q
            test_norm = g1*p*p*root+g2*aa*aa
            test_dnorm = g1*aa*p*p*root+g2*bb*aa*aa
            assert sp.factor(test_norm*28*root*(1+208*a)/test_dnorm) == required
            norms.append(required)
        assert len(set(roots)) == 4
        assert m.subs(u,0) != 0
        # At each root the other exterior-degree block, and the other pair, stay regular.
        for index,(name,degree,copy,d,e,p,q) in enumerate(PAIRS):
            other_degree = 5 if degree == 3 else 3
            determinant = blocks[f'Lam{other_degree}'].det().subs({k:1,sigma:a,theta:0,delta:-u})
            assert squared_substitution(determinant,u,roots[index]) != 0
        records.append((a,hook,dimension,roots,norms))
        print(f'{a} | {hook} ({dimension}) | '+ ' | '.join(map(str,roots)))
    print('All 24 squared intersections are positive; each line has four distinct roots, no detM collision, and one hook copy.')
    print('At u=0 only that hook copy is singular; at each nonzero point add one Lambda3 (364) or Lambda5 (2002) copy.')
    return records


def kernel_controls(records):
    partner = T.scale(ALG.phig,sp.I)
    z1 = T.scale(ALG.phi,sp.I)
    a = records[0][0]
    for index,(name,degree,indices,d,e,p,q) in enumerate(PAIRS):
        u = sp.sqrt(records[0][3][index])
        pats = [T.real_even(value) for value in ALG.forms(degree,(0,2,4,8,11)[:degree])]
        aa,bb = 1+d*a,1+e*a
        field = X.fplus(T.scale(pats[indices[0]],p*u),T.scale(pats[indices[1]],aa))
        to = T.invariant(a,u)
        image = X.fplus(field,T.mix(field,to,(1,0,0,0)))
        assert not any(X.normalized(image))
        assert T.antihermitian(field)
        norm = T.pair(field,field)
        dnorm = norm + a*T.pair(field,T.mix(field,z1,(1,0,0,0)))
        assert T.pair(partner,T.quadratic(field,(1,0,0,0))) == -dnorm/(2*u)
        required_dnorm = 28*u*u*(1+208*a)
        required_norm = sp.factor(norm*required_dnorm/dnorm)
        assert required_norm == records[0][4][index]
    print('Four independent Clifford controls verify the Lambda kernels and the Gamma-trace norm formula.')


def hook_pattern(p,q,r):
    value = [dict() for _ in range(X.N)]
    value[p] = ALG.mul({1<<q:1},{1<<r:1})
    value[q] = ALG.mul({1<<p:1},{1<<r:1})
    contraction = {}
    for a in range(X.N):
        contraction = X.plus(contraction,ALG.mul({1<<a:1},value[a]))
    assert contraction == {}
    return T.real_even(ALG.gamma_right(value))


def grade_obstruction():
    a = R(3,136)
    field = X.fplus(hook_pattern(0,1,2),hook_pattern(0,3,4))
    assert T.antihermitian(field) and T.pair(field,field) == 4
    to = T.invariant(a,0)
    assert not any(X.normalized(X.fplus(field,T.mix(field,to,(1,0,0,0)))))
    quadratic = T.quadratic(field,(1,0,0,0))
    assert {key.bit_count() for item in quadratic for key in item} == {1,5}
    fifth = [{key:value for key,value in item.items() if key.bit_count()==5} for item in quadratic]
    expected = T.scale(ALG.forms(4,(1,2,3,4))[1],4*sp.I)
    assert X.normalized(fifth) == X.normalized(expected) and any(fifth)
    for direction in range(X.N):
        derivative = T.derivative(field,direction,(1,0,0,0))
        assert all(key.bit_count()==11 for item in derivative for key in item)
    trace_g = sp.factor((T.pair(field,field)/2-14*(1+312*a)*a*a)/a)
    cv = sp.factor(105/trace_g)
    assert trace_g == R(152995,1734) and cv == R(36414,30599)
    baseline = T.baseline(a,0,1,[trace_g/14]*X.N)
    assert T.pair(to,baseline) == T.pair(field,field)/2
    scalar_remainder = X.fplus(baseline,quadratic)
    assert T.pair(T.scale(ALG.phi,sp.I),scalar_remainder) == 0
    assert T.pair(T.scale(ALG.phig,sp.I),scalar_remainder) == 0
    print('H2G witness: k=i(h(0,1;2)+h(0,3;4))Gamma, norm^2=4, Lk=0.')
    print('Its grade-five quadratic term is exactly 4i*w(E1234), nonzero; all 14 tested derivative symbols have only grade11.')
    print(f'Flat Z1 control: lambda=1/5, kappa=1, cV={cv}, traceG={trace_g}; both invariant odd traces vanish, but grade5 obstructs the full equation.')
    print('The derivative d0k remains in the necessary grade11 equation; no assumption of a parallel k is made.')


def main():
    records = catalogue()
    kernel_controls(records)
    grade_obstruction()
    print('Exact catalogue and bounded obstructions passed. Constant hook solutions and the other tensor equations remain open.')


if __name__ == '__main__':
    main()
