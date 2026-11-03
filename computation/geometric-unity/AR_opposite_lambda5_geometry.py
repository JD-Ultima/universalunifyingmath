"""Exact scalar and frame controls for AR's opposite-Lambda5 geometry.

No Clifford/metric engine, existing check main or private source is imported.
The Euler controls use algebraic tensors satisfying N1. The ordered-frame
controls verify the h/A contraction signs, not curvature of a constructed
observer metric. Reconstruction, ODE uniqueness and rational continuation
remain hand proofs with the geometric dependencies stated in the AR note.
"""
from time import perf_counter

import sympy as sp


CHECKS = []


def check(label, condition):
    good = bool(condition)
    CHECKS.append(good)
    print(f'[{len(CHECKS):02d}] {label}: {good}', flush=True)
    assert good, label


def euler_controls():
    a, x, b, L, W, zeta = sp.symbols('a x b L W zeta', real=True)
    c0 = sp.symbols('c0', positive=True)
    U = (208*a+1)*(21840*a+113)/480
    left, right = -sp.Rational(1, 160), -sp.Rational(3, 704)
    roots = (-sp.Rational(113, 21840), -sp.Rational(1, 208))
    check('U roots lie strictly inside the forbidden conic gap',
          left < roots[0] < roots[1] < right
          and all(U.subs(a, root) == 0 for root in roots))
    left_poly = sp.Poly(sp.expand(U.subs(a, left-x)), x)
    right_poly = sp.Poly(sp.expand(U.subs(a, right+x)), x)
    check('U is strictly positive on the entire real left arc',
          all(value > 0 for value in left_poly.all_coeffs()))
    check('U is strictly positive on the entire real right arc',
          all(value > 0 for value in right_poly.all_coeffs()))
    print('  U(left-x) =', left_poly.as_expr(), flush=True)
    print('  U(right+x) =', right_poly.as_expr(), flush=True)
    check('dimension14 constant negative Einstein tensor has positive Ricci',
          -b-(-14*b)/12 == b/6)
    # This is an algebraic frame with horizontal dimension4, vertical
    # dimension10 and the prescribed Euler direction. It is not a metric jet.
    identity = sp.eye(14)
    PH = sp.diag(*([1]*4+[0]*10))
    C = sp.Matrix(sp.symbols('C0:14', real=True))
    E = sp.zeros(14, 1)
    E[13] = sp.sqrt(c0)
    rho = sp.symbols('rho0:13', real=True)
    Ric = sp.diag(*rho, -1/c0)
    R04, R17, R912 = sp.symbols('R04 R17 R912', real=True)
    for j, k, value in ((0, 4, R04), (1, 7, R17), (9, 12, R912)):
        Ric[j, k] = Ric[k, j] = value
    check('algebraic Euler data have the N1 length and Ricci eigenvalue',
          (E.T*E)[0] == c0 and PH*E == sp.zeros(14, 1)
          and Ric*E == -E/c0 and Ric == Ric.T)
    c = (C.T*C)[0]
    p = (C.T*E)[0]
    CH = PH*C
    Hess = (Ric-L*identity+W*C*C.T)/24
    direct = Hess*E+CH/2
    expected = -(1/c0+L)*E/24+W*p*C/24+CH/2
    check('all14 components of the covector Euler transport',
          (direct-expected).applyfunc(sp.expand) == sp.zeros(14, 1))
    Ep = (E.T*direct)[0]
    check('Euler derivative of p retains its full coefficient',
          sp.expand(24*Ep-(W*p*p-1-c0*L)) == 0)
    Ec = 2*(C.T*Hess*E)[0]
    check('Euler norm transport from the full Hessian',
          sp.expand(Ec-p*(-1/c0-L+W*c)/12) == 0)
    Zp = (C.T*direct)[0]
    check('gradient transport of p includes the horizontal norm term',
          sp.expand(Zp-Ec/2-(CH.T*CH)[0]/2) == 0)
    check('self-adjoint Euler contraction gives the Ricci dichotomy',
          sp.expand((E.T*(Ric*C-zeta*C))[0]+(zeta+1/c0)*p) == 0)
    check('Euler Cauchy defect is the full perpendicular norm',
          sp.expand(c0*c-p*p-c0*sum(C[j]**2 for j in range(13))) == 0)
    c_a, H, h, T = sp.symbols('c_a H h T', real=True)
    eigenvalue = 12*c_a/H+2*h+T/6-sp.Rational(11, 12)*W*c
    metric_L = 2*h+T/6+W*c/12
    check('AQ Ricci-gradient eigenvalue agrees with Euler transport',
          sp.expand(eigenvalue-(12*c_a/H+metric_L-W*c)) == 0)


def reparametrization_controls():
    a, k, B, L = sp.symbols('a k B L', real=True)
    mu, fp = sp.symbols('mu fp', positive=True)
    W = sp.Rational(1024, 285)*(53-676*a)
    f = (1+sp.Rational(704, 3)*a)*(1+160*a)/1280
    delta = sp.expand(1+224*a-3328*f)
    check('opposite-conic delta factorization',
          sp.factor(delta+sp.Rational(8, 15)*(208*a+1)*(880*a+3)) == 0)
    # F'_chi=fp=24mu; this F'' cancels the full C tensor C term.
    fpp = -W*fp/24
    check('full quadratic-gradient Hessian term cancels',
          sp.expand(fp*W/24+fpp) == 0)
    check('Ricci and free scalar coefficients are mu and minus mu*L',
          (fp/24).subs(fp, 24*mu) == mu
          and (-fp*L/24).subs(fp, 24*mu) == -mu*L)
    dR_da = 576*mu/delta
    dmu_da = -W*mu/delta
    kk = sp.Rational(16, 2565)*(676*a-53)
    check('mu derivative in R includes both chain-rule factors',
          sp.factor(dmu_da/dR_da-kk) == 0 and sp.factor(kk+W/576) == 0)
    ak = (2565*k+848)/10816
    check('affine k has a genuine rational inverse',
          sp.expand(kk.subs(a, ak)-k) == 0)
    check('da equals delta times (dR/mu)/576',
          sp.factor(1/dR_da-delta/(576*mu)) == 0)
    Q = sp.diff(kk, a)*delta/576
    check('mu times its second R derivative',
          sp.factor(Q-sp.Rational(169, 23085)*delta) == 0)
    Qk = -(57*k+20)*(141075*k+48668)/sp.Integer(200070)
    check('Q is the exact rational quadratic in k', sp.factor(Q.subs(a, ak)-Qk) == 0)
    check('Q has nonzero value at the rational-square pole',
          Qk.subs(k, 0) == -sp.Rational(97336, 20007))
    X = 1-1/(B*k)
    rho = 1+4*k*k*X/Qk
    check('Q7 D gives precisely the claimed rational rho',
          sp.cancel(Qk*(1-rho)+4*k*k*X) == 0)
    denominator = sp.Integer(200070)
    Y = sp.expand(denominator*Qk)
    A = B*Y+4*denominator*B*k*k-4*denominator*k
    polynomial = k*A*A-B*(B*k-1)*Y*Y
    check('denominator-cleared square relation is exact',
          sp.cancel(polynomial-B*B*k*Y*Y*(rho*rho-X)) == 0)
    check('cleared identity has nonzero constant coefficient for B!=0',
          sp.expand(polynomial.subs(k, 0)-B*973360**2) == 0)
    print('  k(a) =', sp.factor(kk), flush=True)
    print('  a(k) =', sp.factor(ak), flush=True)
    print('  Q(k) =', sp.factor(Qk), flush=True)
    print('  Cleared constant coefficient =', sp.factor(polynomial.subs(k, 0)), flush=True)


def horizontal_trace_controls():
    # Four horizontal and ten vertical orthonormal directions. b is any
    # symmetric vertical-valued horizontal tensor; A is any skew one.
    # All 160 independent components remain symbolic in these exact sums.
    b = [[sp.zeros(10, 1) for _ in range(4)] for _ in range(4)]
    A = [[sp.zeros(10, 1) for _ in range(4)] for _ in range(4)]
    for i in range(4):
        for j in range(i, 4):
            value = sp.Matrix(sp.symbols(f'b{i}{j}_0:10', real=True))
            b[i][j] = b[j][i] = value
        for j in range(i+1, 4):
            value = sp.Matrix(sp.symbols(f'A{i}{j}_0:10', real=True))
            A[i][j], A[j][i] = value, -value
    B = [[b[i][j]+A[i][j] for j in range(4)] for i in range(4)]
    dot = lambda u, v: (u.T*v)[0]
    mean = sum((b[i][i] for i in range(4)), sp.zeros(10, 1))
    bnorm = sum(dot(b[i][j], b[i][j]) for i in range(4) for j in range(4))
    Anorm = sum(dot(A[i][j], A[i][j]) for i in range(4) for j in range(4))
    mnorm = dot(mean, mean)
    cross = sum(2*dot(b[i][j], A[i][j]) for i in range(4) for j in range(4))
    check('ordered horizontal symmetric/skew cross term cancels identically',
          sp.expand(cross) == 0)
    raw_HH = sum(-dot(B[j][j], B[i][i])+dot(B[i][j], B[j][i])
                 +2*dot(A[i][j], B[j][i])
                 for i in range(4) for j in range(4))
    check('complete ordered HH contraction has the minus3 integrability norm',
          sp.expand(raw_HH-(-mnorm+bnorm-3*Anorm)) == 0)
    divmean, base_trace, c0 = sp.symbols('divmean base_trace c0', real=True)
    walczak = divmean-bnorm+mnorm+Anorm
    actual = (base_trace+raw_HH+walczak).subs(divmean, -4/c0)
    expected = base_trace-4/c0-2*Anorm
    check('Walczak addition gives minus4/c0 minus2 integrability norm',
          sp.expand(actual-expected) == 0)
    print('  Frame signs checked with all160 independent b/A components.', flush=True)


def profile_controls():
    a, c, cp, kappa = sp.symbols('a c c_prime kappa', real=True)
    c0 = sp.symbols('c0', positive=True)
    U = (208*a+1)*(21840*a+113)/480
    V = -sp.Rational(64, 285)*(402896*a+1719)
    W = sp.Rational(1024, 285)*(53-676*a)
    f = (1+sp.Rational(704, 3)*a)*(1+160*a)/1280
    delta = sp.expand(1+224*a-3328*f)
    D = -sp.Rational(64, 285)*(343408*a+6383)
    check('implicit profile coefficient D equals V minus11W/2',
          sp.expand(D-V+sp.Rational(11, 2)*W) == 0)
    # This h is obtained from Ric(C)=-C/c0 on a regular patch, without
    # dividing by D. The original Bianchi covector is proportional to C.
    h = -delta*cp/4-kappa*kappa*U/12-(2*V-11*W)*c/24-1/(2*c0)
    bianchi = ((V+W/2-288)*delta*cp/24+(W-48)*h
               +kappa*kappa*(delta*sp.diff(U, a)/24-4*U)
               +c*(delta*(sp.diff(V, a)+sp.diff(W, a))/24+22*W-4*V))
    profile = (D*cp+(sp.diff(V, a)+sp.diff(W, a)-2*W*D/delta)*c
               +kappa*kappa*(sp.diff(U, a)-2*W*U/delta)
               -12*(W-48)/(c0*delta))
    check('fixed Ricci-gradient equation reduces Bianchi to the exact linear profile',
          sp.cancel(bianchi-delta*profile/24) == 0)
    zero = -sp.Rational(6383, 343408)
    other_factor = 1+544*a-3328*f
    check('D zero is a real regular left-arc level, so no division by D is allowed',
          D.subs(a, zero) == 0 and zero < -sp.Rational(1, 160)
          and f.subs(a, zero) > 0 and delta.subs(a, zero) != 0
          and other_factor.subs(a, zero) != 0)
    print('  D(a) =', sp.factor(D), '; zero =', zero, flush=True)


def moving_endpoint_controls():
    a = sp.symbols('a', real=True)
    f = (1+sp.Rational(704, 3)*a)*(1+160*a)/1280
    delta = sp.expand(1+224*a-3328*f)
    U = (208*a+1)*(21840*a+113)/480
    # a, delta, f', delta*f', Hess(chi)/(dw tensor dw), U.
    rows = ((-sp.Rational(1, 160), -sp.Rational(2, 5), -sp.Rational(7, 120),
             sp.Rational(7, 300), sp.Rational(14400, 7), sp.Rational(47, 3200)),
            (-sp.Rational(3, 704), sp.Rational(1, 22), sp.Rational(7, 120),
             sp.Rational(7, 2640), sp.Rational(126720, 7), sp.Rational(877, 185856)))
    check('both opposite endpoints have the exact nonzero delta and conic slope',
          all(f.subs(a, endpoint) == 0 and delta.subs(a, endpoint) == dd
              and sp.diff(f, a).subs(a, endpoint) == slope and dd*slope == product
              for endpoint, dd, slope, product, factor, uu in rows))
    check('twice-differentiated conic gives a positive semidefinite Hessian of chi',
          all(24/dd*2/slope == factor and factor > 0 and product > 0
              for endpoint, dd, slope, product, factor, uu in rows))
    check('both endpoint scalar Ricci contributions are strictly positive',
          all(U.subs(a, endpoint) == uu and uu > 0
              for endpoint, dd, slope, product, factor, uu in rows))
    c0 = sp.symbols('c0', positive=True)
    euler_dw_squared, norm_dw, kappa_squared = sp.symbols(
        'euler_dw_squared norm_dw kappa_squared', nonnegative=True)
    # At the endpoint C=0, so Ric=24 Hess(chi)+(2h+kappa^2 U/6)g.
    # Subtracting N1's Ric(E,E)=-1 yields the polynomial below. Its
    # coefficients are strictly positive, proving the contradiction even
    # at a higher-order contact dw=0. No divided norm formula is extended.
    check('endpoint Ricci and the Euler identity have a strictly positive gap',
          all(all(coefficient > 0 for coefficient in sp.Poly(
              24*factor*euler_dw_squared+(2*factor*norm_dw+kappa_squared*uu/6)*c0+1,
              euler_dw_squared, norm_dw, kappa_squared).coeffs())
              for endpoint, dd, slope, product, factor, uu in rows))
    print('  Moving endpoint Hessian factors =', tuple(row[4] for row in rows), flush=True)


def main():
    start = perf_counter()
    print('AR: opposite moving Lambda5 Euler and Ricci-Hessian geometry.', flush=True)
    print('Exact scalar/frame controls; no observer metric or PDE solution is constructed.', flush=True)
    euler_controls()
    reparametrization_controls()
    horizontal_trace_controls()
    profile_controls()
    moving_endpoint_controls()
    print('Density, ODE uniqueness, global reconstruction and continuation use the stated hand proofs.', flush=True)
    print(f'{len(CHECKS)} exact controls passed; runtime {perf_counter()-start:.3f} s.', flush=True)


if __name__ == '__main__':
    main()
