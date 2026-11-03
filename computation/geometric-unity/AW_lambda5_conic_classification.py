"""Exact controls for AW: both single Lambda5 conics on positive Z1.

The forms deliberately need not satisfy the actual Gram or AS constraints.
They test cyclic-trace signs on noncommuting inputs, not admissible data or
full-system existence. Scalar rational certificates support the separate hand proof.
No old helper, private file or metric engine is imported.
"""
from fractions import Fraction as F
from time import perf_counter
import sympy as sp
N = 14
FULL = (1 << N) - 1
COUNT = 0

def add(a, b, factor=1):
    out = dict(a)
    for mask, value in b.items():
        out[mask] = out.get(mask, 0) + factor * value
    return {mask: value for mask, value in out.items() if value != 0}

def scale(a, factor):
    return {mask: value * factor for mask, value in a.items() if value * factor != 0}

def mul(a, b):
    out = {}
    for left, lv in a.items():
        for right, rv in b.items():
            parity = 0
            bits = left
            while bits:
                j = (bits & -bits).bit_length() - 1
                parity += (right & (1 << j) - 1).bit_count()
                bits &= bits - 1
            mask = left ^ right
            out[mask] = out.get(mask, 0) + lv * rv * (-1) ** parity
    return {mask: value for mask, value in out.items() if value != 0}

def grade(a, d):
    return {mask: value for mask, value in a.items() if mask.bit_count() == d}

def trace(a):
    return a.get(0, 0)

def pair(a, b):
    return trace(mul(a, b))

def reverse(a):
    return {mask: value * (-1) ** (mask.bit_count() * (mask.bit_count() - 1) // 2) for mask, value in a.items()}

def form(indices):
    return {sum((1 << i for i in indices)): 1}

def linear(basis, coefs):
    out = {}
    for indices, value in zip(basis, coefs):
        out = add(out, scale(form(indices), value))
    return out

def check(label, value):
    global COUNT
    COUNT += 1
    assert bool(value), label
    print(f'[{COUNT:02}] {label}: True', flush=True)

def clifford_controls():
    e0 = {1: 1}
    J = {FULL ^ 1: 1}
    Gamma = {FULL: 1}
    ub = ((2, 3, 4, 5), (2, 3, 6, 7), (4, 5, 6, 7), (6, 7, 8, 9), (2, 4, 5, 6))
    pb = ((1, 2, 3, 4, 5), (1, 6, 7, 8, 9), (1, 10, 11, 12, 13), (1, 2, 3, 6, 7))
    cases = (((1, 2, -1, 3, 2), (1, -2, 3, 2)), ((2, -1, 3, 1, -2), (-2, 1, 2, 3)), ((-1, 3, 2, -2, 1), (3, 2, -1, 1)))
    for i, (uc, pc) in enumerate(cases, 1):
        U, psi = (linear(ub, uc), linear(pb, pc))
        V = mul(psi, J)
        Omega = add(mul(e0, U), psi)
        e0U = mul(e0, U)
        square = mul(Omega, Omega)
        sandwich = mul(mul(Omega, e0), Omega)
        Y5 = grade(sandwich, 5)
        Y9G = mul(grade(sandwich, 9), Gamma)
        E5 = grade(mul(e0, square), 5)
        E9G = mul(grade(mul(e0, square), 9), Gamma)
        m, ell = (pair(U, U), pair(V, V))
        M30 = trace(mul(mul(U, U), U))
        M12 = trace(mul(U, mul(V, V)))
        M21 = trace(mul(mul(U, U), V))
        M03 = trace(mul(mul(V, V), V))
        check(f'case{i}: positive Hermitian form phases and genuinely noncommuting U,V', all((reverse(v) == v for v in (U, V, psi, Omega, e0U))) and m > 0 and (ell > 0) and bool(add(mul(U, V), mul(V, U), -1)) and (M12 != 0) and (M21 != 0))
        actual = (pair(e0U, Y5), pair(e0U, Y9G), pair(psi, Y5), pair(psi, Y9G))
        target = (M30 - M12, -2 * M21, 2 * M12, M21 - M03)
        check(f'case{i}: all four full normal sandwich pairings, including factor2', actual == target)
        check(f'case{i}: square pairings and the normal E-channel pairings are exact', pair(U, square) == M30 + M12 and pair(V, square) == M21 + M03 and (pair(e0U, E5) == M30 + M12) and (pair(psi, E9G) == M21 + M03) and (pair(psi, E5) == pair(e0U, E9G) == pair(U, V) == 0))
        for tau in (F(1), F(2, 3)):
            rho = F(5, 7)
            Ac = -F(4, 3) * tau * (41 + 11 * tau * tau)
            Bc = F(8, 3) * (38 + 11 * tau * tau)
            ec = F(8, 3) * tau * (4 + tau * tau)
            fc = F(8, 3) * (4 + tau * tau)
            normal = {}
            for value, coeff in ((Y5, Ac), (E5, ec), (Y9G, Bc), (E9G, fc), (e0U, 16 * tau * rho), (psi, -32 * tau * rho)):
                normal = add(normal, scale(value, coeff))
            expected1 = Ac * (M30 - M12) - 2 * Bc * M21 + ec * (M30 + M12) + 16 * tau * rho * m
            expected2 = 2 * Ac * M12 + Bc * (M21 - M03) + fc * (M21 + M03) - 32 * tau * rho * ell
            assert (pair(e0U, normal), pair(psi, normal)) == (expected1, expected2)
        check(f'case{i}: full AS normal component has the two stated signed moment equations', True)
        print('  moments M30,M12,M21,M03 =', (M30, M12, M21, M03), flush=True)


def scalar_controls():
    tau, theta, rho = sp.symbols('tau theta rho', nonzero=True, real=True)
    x = tau*tau
    alpha = -3*(304+9*theta)/(16*(61+16*x))
    zeta = 3*tau*(-112+5*theta)/(16*(64+19*x))
    Ac = -sp.Rational(4,3)*tau*(41+11*x)
    Bc = sp.Rational(8,3)*(38+11*x)
    Ec = sp.Rational(8,3)*tau*(4+x)
    Fc = sp.Rational(8,3)*(4+x)
    p, q = Ec*alpha+16*tau, 32*tau-Fc*zeta
    R1, R2 = sp.factor(Ac*alpha+p), sp.factor(Bc*zeta+q)

    M, L, u30, u12, u21, u03 = sp.symbols('M L u30 u12 u21 u03', real=True)
    normal4 = Ac*(u30-u12)-2*Bc*u21+rho*p*M
    normal8 = Bc*(u21-u03)+2*Ac*u12-rho*q*L
    squareU = u30+u12-alpha*rho*M
    squareV = u21+u03-zeta*rho*L
    check('normal plus square trace residuals cancel every third moment',
          sp.factor(normal4+normal8-Ac*squareU+Bc*squareV-rho*(R1*M-R2*L)) == 0)

    a, c, kk, t = sp.symbols('a c kappa_squared t', real=True)
    w_over_tau = (1+160*a)/48
    xs = 9*(1+sp.Rational(704,3)*a)/(5*(1+160*a))
    w_squared = (1+sp.Rational(704,3)*a)*(1+160*a)/1280
    delta = -sp.Rational(8,15)*(208*a+1)*(880*a+3)
    thet = sp.factor(sp.Rational(3328,3)*w_over_tau*(xs+4)-7*delta/(360*w_squared))
    check('delta equals the exact conic expression', sp.factor(delta-(1+224*a-3328*w_squared)) == 0)
    ratio1 = sp.factor((R1/tau).subs(theta,thet).subs(tau*tau,xs))
    ratio2 = sp.factor((R2/tau).subs(theta,thet).subs(tau*tau,xs))
    check('normal moment coefficient ratios have no remaining tau', not (ratio1.has(tau) or ratio2.has(tau)))
    N = sp.factor(w_over_tau*(sp.Rational(7,120)*kk*(1+208*a)-sp.Rational(2704,15)*c))
    beta = 13*(1+160*a)/19
    M_actual, L_actual = (5*N-13*beta*c)/14, (9*N+13*beta*c)/14
    check('normal and tangent norm pieces sum to actual norm', sp.factor(M_actual+L_actual-N) == 0)
    check('actual Gram normal eigenvalue is M', sp.factor((5*N+beta*c)/14-beta*c-M_actual) == 0)
    relation = sp.factor(ratio1*M_actual-ratio2*L_actual)
    num, den = sp.fraction(relation)
    check('moment relation is exactly linear in c and kappa squared', sp.Poly(num,c,kk).total_degree() == 1)
    F = sp.factor(sp.solve(num,c)[0]/kk)
    P5 = (7456962482284789760*a**5+205388793727942656*a**4
          +2208655001075712*a**3+11621448908288*a**2+29961444000*a+30290493)
    Q5 = (151851622538936320*a**5+11761525204713472*a**4
          +206238452785152*a**3+1516403340544*a**2+5097590656*a+6480165)
    targetF = -19*(208*a+1)*P5/(173056*Q5)
    check('actual Gram and full normal compatibility force stated rational norm', sp.factor(F-targetF) == 0)
    forced_polynomial = 173056*Q5*c+19*kk*(208*a+1)*P5
    relation_factor = sp.factor(relation/forced_polynomial)
    check('uncancelled norm identity has only a scalar rational prefactor',
          not relation_factor.has(c,kk) and sp.factor(relation-relation_factor*forced_polynomial)==0)
    check('norm denominator is a nonzero polynomial of degree five', sp.degree(Q5,a)==5 and sp.LC(sp.Poly(Q5,a)) != 0)
    check('rational norm numerator and denominator are coprime', sp.gcd(sp.Poly(sp.fraction(F)[0],a),sp.Poly(sp.fraction(F)[1],a)).degree()==0)

    U = (208*a+1)*(21840*a+113)/480
    V = -sp.Rational(64,285)*(402896*a+1719)
    W = sp.Rational(1024,285)*(53-676*a)
    D = sp.factor(V-sp.Rational(11,2)*W)
    P = sp.factor(sp.diff(V,a)+sp.diff(W,a)-2*W*D/delta)
    Q = sp.factor(sp.diff(U,a)-2*W*U/delta)
    R = sp.factor(-12*(W-48)/delta)
    G = sp.factor(D*sp.diff(F,a)+P*F+Q)
    flead = sp.limit(F/a,a,sp.oo)
    glead = sp.limit(G/a,a,sp.oo)
    check('forced norm has exact nonzero linear rational coefficient', flead == -sp.Rational(35575543,31723328))
    check('AR reduced profile has exact nonzero linear rational coefficient', glead == sp.Rational(3138842736144,15324925))
    check('the remaining AR source divided by a has zero rational limit', sp.limit(R/a,a,sp.oo) == 0)
    residual, residual_den = sp.fraction(sp.factor(t*G+R))
    respoly = sp.Poly(residual,a)
    check('cleared AR identity has exact degree thirteen', respoly.degree()==13)
    check('highest cleared coefficient is nonzero times t', respoly.LC() == 3202906730881248700955630955569647392741811814400*t)
    check('same contradiction follows directly from rational leading orders', sp.limit((t*G+R)/a,a,sp.oo) == t*glead)


    # Independently recover the actual normal Gram norm from the full Gamma trace.
    gamma_norm = sp.symbols('gamma_norm')
    wt = (1+160*a)/48
    trace_equation = (-14*kk*wt*(1+208*a)+14*3328*wt*c
                      +sp.Rational(2240,3)*gamma_norm-3328*wt*c
                      -sp.Rational(1520,3)*gamma_norm)
    check('full converted Gamma trace recovers the actual form norm',
          sp.factor(sp.solve(trace_equation,gamma_norm)[0]-N)==0)
    Mc = sp.factor(sp.diff(relation,c))
    check('the norm coefficient is nonzero as a rational function',
          Mc.subs(a,0)==-sp.Rational(4672630976,37698489))
    # Independent elimination of div C from AQ's tensor Bianchi and fixed Ricci alignment.
    cp,h,c0 = sp.symbols('cp h c0',real=True)
    alignment = (delta*cp/2+2*h+kk*U/6+V*c/6-sp.Rational(11,12)*W*c+1/c0)
    h0 = sp.solve(alignment,h)[0]
    bianchi = ((V+W/2-288)*delta*cp/24+(W-48)*h
               +kk*(delta*sp.diff(U,a)/24-4*U)
               +c*(delta*(sp.diff(V,a)+sp.diff(W,a))/24+22*W-4*V))
    profile = D*cp+P*c+kk*Q+R/c0
    check('direct AQ elimination recovers AR5b with its sign and normalization',
          sp.factor(24*bianchi.subs(h,h0)/delta-profile)==0)
    G0=sp.factor(G.subs(a,0))
    R0=sp.factor(R.subs(a,0))
    check('the rational profile coefficient has the exact finite value at zero',
          G0==sp.Rational(35018717715292549351,51860784957622875))
    check('the rational source has the exact finite value at zero',
          R0==sp.Rational(20296,19))
    check('finite-value certificate also contradicts every positive constant scale',G0>0 and R0>0)
    print('  Moment norm coefficient at a=0 =',Mc.subs(a,0),flush=True)
    print('  Leading F/a =',flead,flush=True)
    print('  Leading (D F_prime+P F+Q)/a =',glead,flush=True)
    print('  Leading source R/a = 0',flush=True)
    print('  Finite coefficient H(0) =',G0,'; source R(0) =',R0,flush=True)
    print('  Cleared residual degree =',respoly.degree(),'; leading coefficient =',respoly.LC(),flush=True)


if __name__ == '__main__':
    started = perf_counter()
    clifford_controls()
    scalar_controls()
    print('Universality, finite-exception avoidance and rational-identity contradiction use the stated hand proof.', flush=True)
    print('No observer field or inequality is continued to zero or infinity.', flush=True)
    print(f'{COUNT} exact controls passed; runtime {perf_counter()-started:.3f} s.', flush=True)
