"""Symbolic quadratic gradient on the common-axis symmetric hook family.

For a fixed unit r, B_a=sum_b S_ab gamma_b gamma_r on r-perp, B_r=0,
S real symmetric and traceless, k=i B Gamma. Uses the earlier independent
Clifford implementation through T's compact-real-form helpers. No imports
from the original R/J evaluator. Thirteen independent symbolic eigenvalues
give a polynomial identity, followed by the trace-zero restriction.

Equivariance and SO(13) diagonalization extend the identity to every S in
this declared family. They do not cover arbitrary H2Gamma fields.
"""
import pathlib
import sys
import sympy as sp

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
import T_singular_residual as T

X,ALG = T.X,T.ALG
R = sp.Rational
AXIS = 13
WA = (1,0,0,0)


def field_from_matrix(matrix):
    bivectors = [dict() for _ in range(14)]
    for a in range(13):
        for b in range(13):
            if matrix[a,b] != 0:
                bivectors[a] = X.plus(bivectors[a],ALG.mul({1<<b:matrix[a,b]},{1<<AXIS:1}))
    return bivectors,T.real_even(ALG.gamma_right(bivectors))


def zero(field):
    return not any(X.normalized(field))


def main():
    eigenvalues = sp.symbols('lambda0:13',real=True)
    matrix = sp.diag(*eigenvalues)
    bivectors,field = field_from_matrix(matrix)
    total = sum(eigenvalues)
    norm = sum(value*value for value in eigenvalues)
    contraction = {}
    for a in range(14):
        contraction = X.plus(contraction,ALG.mul({1<<a:1},bivectors[a]))
    assert X.normalized([contraction])[0] == {1<<AXIS:total}
    assert T.antihermitian(field)
    assert sp.expand(T.pair(field,field)-norm) == 0

    quadratic = T.quadratic(field,WA)
    expected = []
    for a in range(13):
        coefficient = 2*((total-eigenvalues[a])**2-(norm-eigenvalues[a]**2))
        expected.append({1<<a:sp.I*coefficient})
    expected.append({1<<AXIS:2*sp.I*(total*total-norm)/3})
    assert X.normalized(quadratic) == X.normalized(expected)
    assert all(key.bit_count()==1 for item in quadratic for key in item)
    print('13-symbol identity: Q_a/i=2[(L-lambda_a)^2-(N-lambda_a^2)] on r-perp; Q_r/i=(2/3)(L^2-N).')
    print('All 14 component polynomials match; off-diagonal and grade-five components vanish. L=sum(lambda), N=sum(lambda^2).')

    pair_values = {value:int(index in (0,1)) for index,value in enumerate(eigenvalues)}
    pair_field = X.normalized([{key:value.subs(pair_values) for key,value in item.items()} for item in field])
    pair_q = T.quadratic(pair_field,WA)
    first = T.scale(ALG.dual('A',ALG.wedge(pair_field,pair_field)),sp.I)
    second = X.fplus(T.scale(pair_q,3),first,-1)
    expected_first = [{1<<a:4*sp.I} if a not in (0,1) else {} for a in range(14)]
    expected_second = [{1<<a:8*sp.I} if a not in (0,1,AXIS) else {} for a in range(14)]
    assert X.normalized(first) == X.normalized(expected_first)
    assert X.normalized(second) == X.normalized(expected_second)
    print('Pair coefficient check in 3Q: first term 4i off the two input indices; adjoint term 8i on the remaining W indices and zero at r.')

    restriction = {eigenvalues[-1]:-sum(eigenvalues[:-1])}
    restricted_field = [{key:sp.expand(value.subs(restriction)) for key,value in item.items()} for item in field]
    restricted_norm = sp.expand(norm.subs(restriction))
    restricted_q = [{key:sp.expand(value.subs(restriction)) for key,value in item.items()} for item in quadratic]
    trace_w = sp.expand(sum(restricted_q[a].get(1<<a,0)/sp.I for a in range(13)))
    rr = sp.expand(restricted_q[AXIS].get(1<<AXIS,0)/sp.I)
    assert sp.expand(trace_w+22*restricted_norm) == 0
    assert sp.expand(rr+2*restricted_norm/3) == 0
    assert sp.expand(trace_w-13*rr+40*restricted_norm/3) == 0
    assert zero(X.fplus(restricted_field,T.mix(restricted_field,T.invariant(R(3,136),0),WA)))
    print('Trace-zero family: Q_W/i=4*S^2-2*N*I_13; Q_rr/i=-2*N/3; norm(k)^2=N; Lk=0 at t=3*kappa/136,u=0.')
    print('Anisotropy trace: tr_W(Q/i)-13*(Q_rr/i)=-40*N/3. A scalar quadratic tensor therefore requires N=0.')

    # Rational SO(13) rotation of diag(1,-1,0,...), checked independently.
    rotation = sp.eye(13)
    rotation[0,0],rotation[0,1] = R(3,5),-R(4,5)
    rotation[1,0],rotation[1,1] = R(4,5),R(3,5)
    assert rotation.T*rotation == sp.eye(13) and rotation.det() == 1
    sample = rotation*sp.diag(1,-1,*([0]*11))*rotation.T
    _,rotated_field = field_from_matrix(sample)
    rotated_q = T.quadratic(rotated_field,WA)
    sample_norm = sp.trace(sample*sample)
    coefficients = 4*sample*sample-2*sample_norm*sp.eye(13)
    rotated_expected = [{1<<b:sp.I*coefficients[a,b] for b in range(13) if coefficients[a,b] != 0}
                        for a in range(13)]
    rotated_expected.append({1<<AXIS:-2*sp.I*sample_norm/3})
    assert X.normalized(rotated_q) == X.normalized(rotated_expected)
    assert T.antihermitian(rotated_field) and sample.trace() == 0
    for direction in range(14):
        derivative = T.derivative(rotated_field,direction,WA)
        assert all(key.bit_count()==11 for item in derivative for key in item)
    print('Nondiagonal rational-rotation control passed; all 14 derivative symbols have grade11, not grade1.')

    kappa = sp.Symbol('kappa',real=True,nonzero=True)
    cv = sp.Symbol('cV',positive=True)
    t = 3*kappa/136
    baseline = sp.factor(15/(2*cv)+kappa*t+312*t*t)
    assert sp.expand(baseline-15/(2*cv)-R(201,1156)*kappa*kappa) == 0
    determinant = sp.factor(kappa*kappa*(kappa+224*t)*(kappa+544*t))
    assert sp.expand(determinant-R(1313,17)*kappa**4) == 0
    assert baseline.is_positive and determinant.is_positive
    print('Flat Z1: baseline=15/(2*cV)+201*kappa^2/1156>0; detM=(1313/17)*kappa^4>0.')
    print('The common-axis symmetric hook family cannot solve the grade-one equation, regardless of derivative jets.')
    print('This excludes only the declared family. The general H2Gamma kernel and its grade11 differential equation remain open.')


if __name__ == '__main__':
    main()
