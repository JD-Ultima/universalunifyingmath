"""Independent exact review of R: full even blocks and the bounded R10 continuation.

No repository modules are imported. The Clifford product, Shiab maps, cubic
gradient and its polarization are constructed here from D.1/D5 and I.0.
Sparse monomials use integers/Fractions; i is introduced only as an operator
phase. The block templates are transcribed from R2 and checked against the
independently differentiated cubic, not used to compute its action.

With Tr(1)=1, Psi(alpha,zeta)=sum_{a,b<c}Tr(alpha_a M_abc(zeta_bc)), where
M_A=-{gamma^a gamma^b gamma^c,.}, M_A'=-{gamma^a gamma^b gamma^c Gamma,.},
M_CG=i[gamma^a gamma^b gamma^c,Ad_Gamma(.)] and similarly for its Gamma partner.
The coefficient of 3*Q is the derivative of Psi(T,T wedge T). The coefficient
of 3*2B(X,T) is the cross term of this quadratic gradient.

Checks: every R2 constituent/copy, both T_o maps, all four operators, both
signatures; one independently chosen test vector per constituent. R10 checks
the cubic symmetry and the zero current at the R9 kernel root in u(128).
The representation-theoretic and differential arguments remain in the note;
finite tests do not substitute for them. No claim about nonconstant T_o.
"""
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations
import sympy as sp

N = 14
FULL = (1 << N) - 1
OPS = ('A', 'Ap', 'CG', 'CpG')
PHASE = (-1, -1, sp.I, sp.I)


def plus(a, b, scale=1):
    out = dict(a)
    for key, value in b.items():
        out[key] = out.get(key, 0) + scale * value
    return {key: value for key, value in out.items() if value}


def scaled(a, scale):
    return {key: value * scale for key, value in a.items() if value * scale}


def fplus(a, b, scale=1):
    return [plus(x, y, scale) for x, y in zip(a, b)]


class Algebra:
    def __init__(self, eta):
        self.eta = eta
        self.negative = sum(1 << a for a in range(N) if eta[a] < 0)
        self.phi = [{1 << a: 1} for a in range(N)]
        self.phig = [self.mul(value, {FULL: 1}) for value in self.phi]

    @lru_cache(maxsize=None)
    def sign(self, a, b):
        inversions = sum((a >> (j + 1)).bit_count() for j in range(N) if b & (1 << j))
        return (-1) ** (inversions + (a & b & self.negative).bit_count())

    def mul(self, a, b):
        out = {}
        for i, ai in a.items():
            for j, bj in b.items():
                key = i ^ j
                out[key] = out.get(key, 0) + ai * bj * self.sign(i, j)
        return {key: value for key, value in out.items() if value}

    def tr(self, a, b):
        return sum(value * b.get(key, 0) * self.sign(key, key) for key, value in a.items())

    def pair(self, a, b):
        return sum(self.eta[j] * self.tr(a[j], b[j]) for j in range(N))

    @lru_cache(maxsize=None)
    def triple(self, a, b, c, partner):
        value = self.mul(self.mul({1 << a: self.eta[a]}, {1 << b: self.eta[b]}),
                         {1 << c: self.eta[c]})
        return self.mul(value, {FULL: 1}) if partner else value

    def op(self, op, a, b, c, value, transpose=False):
        gamma = self.triple(a, b, c, op in ('Ap', 'CpG'))
        if op in ('A', 'Ap'):
            return plus(self.mul(gamma, value), self.mul(value, gamma))
        if transpose:
            result = plus(self.mul(value, gamma), self.mul(gamma, value), -1)
            return {key: (-v if key.bit_count() % 2 else v) for key, v in result.items()}
        hatted = {key: (-v if key.bit_count() % 2 else v) for key, v in value.items()}
        return plus(self.mul(gamma, hatted), self.mul(hatted, gamma), -1)

    def wedge(self, x, y):
        return {(b, c): plus(self.mul(x[b], y[c]), self.mul(x[c], y[b]), -1)
                for b, c in combinations(range(N), 2)}

    def dual(self, op, zeta):
        out = [dict() for _ in range(N)]
        for (b, c), value in zeta.items():
            if value:
                for a in range(N):
                    if a not in (b, c):
                        out[a] = plus(out[a], self.op(op, a, b, c, value))
        return out

    def mixed(self, op, x, y):
        """Y coefficient of 3*2B(x,y), without the common operator phase."""
        xy, yx = self.wedge(x, y), self.wedge(y, x)
        out = self.dual(op, {key: plus(xy[key], yx[key]) for key in xy})
        for b, c in combinations(range(N), 2):
            for a in range(N):
                if a in (b, c):
                    continue
                for left, right in ((x, y), (y, x)):
                    if not left[a]:
                        continue
                    transpose = self.op(op, a, b, c, left[a], transpose=True)
                    if not transpose:
                        continue
                    out[c] = plus(out[c], plus(self.mul(transpose, right[b]), self.mul(right[b], transpose), -1))
                    out[b] = plus(out[b], plus(self.mul(transpose, right[c]), self.mul(right[c], transpose), -1), -1)
        return [scaled(value, self.eta[a]) for a, value in enumerate(out)]

    def gamma_right(self, x):
        return [self.mul(value, {FULL: 1}) for value in x]

    def forms(self, m, indices):
        omega = {sum(1 << a for a in indices): 1}
        raw = [self.mul(value, omega) for value in self.phi]
        c = [{key: value for key, value in item.items() if key.bit_count() == m - 1} for item in raw]
        w = [{key: value for key, value in item.items() if key.bit_count() == m + 1} for item in raw]
        if m == 7:
            # c(omega Gamma)=w(omega) Gamma; w(omega Gamma)=c(omega) Gamma.
            return [c, self.gamma_right(w), w, self.gamma_right(c)]
        return [c, w, self.gamma_right(w), self.gamma_right(c)]

    def hook(self, k):
        p, q = 1, 9
        mask = sum(1 << a for a in (0, 3, 5, 8, 12)[:k - 1])
        x = [dict() for _ in range(N)]
        x[p] = self.mul({1 << q: 1}, {mask: 1})
        candidate = self.mul({1 << p: 1}, {mask: 1})
        left = self.mul({1 << p: self.eta[p]}, x[p])
        right = self.mul({1 << q: self.eta[q]}, candidate)
        key = next(iter(left))
        x[q] = scaled(candidate, F(-left[key], right[key]))
        assert plus(left, self.mul({1 << q: self.eta[q]}, x[q])) == {}
        return [x, self.gamma_right(x)]


def templates(epsilon):
    k, sigma, theta, delta = sp.symbols('kappa sigma theta delta')
    r, i, e = sp.Rational, sp.I, epsilon
    blocks = {
        'Lam1': sp.Matrix([[k, 0, 208*e*delta, 0], [0, k+224*sigma, 384*i*e*theta, 16*e*delta],
                           [16*delta, -384*i*theta, k+544*sigma, 0], [0, 208*delta, 0, k]]),
        'Lam3': sp.Matrix([[k+r(112,3)*sigma, 0, r(352,3)*e*delta, 64*i*e*theta],
                           [0,k+240*sigma,320*i*e*theta,32*e*delta],
                           [32*delta,-320*i*theta,k+400*sigma,0],
                           [-64*i*theta,r(352,3)*delta,0,k+r(272,3)*sigma]]),
        'Lam5': sp.Matrix([[k+96*sigma,0,48*e*delta,128*i*e*theta],
                           [0,k+r(704,3)*sigma,256*i*e*theta,r(80,3)*e*delta],
                           [r(80,3)*delta,-256*i*theta,k+r(832,3)*sigma,0],
                           [-128*i*theta,48*delta,0,k+160*sigma]]),
        'Lam7': sp.Matrix([[k+176*sigma,0,0,192*i*e*theta],[0,k+176*sigma,-192*i*theta,0],
                           [0,192*i*e*theta,k+208*sigma,0],[-192*i*theta,0,0,k+208*sigma]])}
    for key, a, b in (('H2',r(-56,3),r(-136,3)),('H4',-24,-40),('H6',r(-88,3),r(-104,3))):
        blocks[key] = sp.Matrix([[k+a*sigma,-32*i*e*theta],[32*i*theta,k+b*sigma]])
    return (k,sigma,theta,delta), blocks


def combine(copies, coefficients):
    out = [dict() for _ in range(N)]
    for copy, coefficient in zip(copies, coefficients):
        out = fplus(out, copy, coefficient)
    return out


def normalized(value):
    return [{key: sp.expand(v) for key,v in item.items() if sp.expand(v) != 0} for item in value]


def positive_quadratic(value, variable):
    poly = sp.Poly(value, variable)
    assert poly.degree() == 2 and poly.LC() > 0
    assert sp.discriminant(poly.as_expr(), variable) < 0


def positivity_certificate():
    """R7: rational interval certificate, without the source's root-isolation route."""
    (k,sigma,theta,delta), blocks = templates(-1)
    m,r = sp.symbols('m r')
    m0 = sp.Rational(25,144)
    for name in ('Lam3','Lam5'):
        polynomial = sp.Poly(sp.expand(9*blocks[name].det().subs(
            {k:1,theta:sp.sqrt(m)*sigma,delta:sp.sqrt(r)})),r)
        a,b,c = polynomial.all_coeffs()
        dm = sp.cancel(sp.diff(polynomial.as_expr(),m)/sigma**2).expand()
        terms = sp.Poly(dm,m,r,sigma)
        alpha = terms.coeff_monomial(m*sigma**2)
        beta = terms.coeff_monomial(r)
        quad = dm.subs({m:0,r:0})
        assert sp.expand(dm-alpha*m*sigma**2-beta*r-quad) == 0
        assert alpha > 0 and beta > 0
        positive_quadratic(quad,sigma)
        b0,c0 = sp.expand(b.subs(m,m0)),sp.expand(c.subs(m,m0))
        assert a > 0
        if name == 'Lam5':
            f = 163840*sigma**2+2304*sigma+9
            g = 229376*sigma**2+1536*sigma+3
            assert sp.expand(c0-f*g/3) == 0
            for value in (b0,f,g):
                positive_quadratic(value,sigma)
        else:
            f = 1024000*sigma**2+5760*sigma+9
            assert sp.expand(c0-(64*sigma+1)**2*f) == 0
            positive_quadratic(f,sigma)
            # B<0 can occur only between two roots, both inside [lo,hi].
            lo,hi = -sp.Rational(149,10000),-sp.Rational(3,1000)
            vertex = sp.solve(sp.diff(b0,sigma),sigma)[0]
            assert lo < vertex < hi and b0.subs(sigma,vertex) < 0
            assert b0.subs(sigma,lo) > 0 and b0.subs(sigma,hi) > 0
            dquad = 43264*sigma**2-1152*sigma-27
            assert sp.expand(4*a*c0-b0*b0+sp.Rational(12992276070400,9)*sigma**2*dquad) == 0
            # Convex quadratic below zero at both endpoints is below zero throughout.
            assert sp.Poly(dquad,sigma).LC() > 0
            assert dquad.subs(sigma,lo) < 0 and dquad.subs(sigma,hi) < 0
        assert b.subs(sigma,0) > 0 and c.subs(sigma,0) > 0
    print("R7 positivity: Lam3/Lam5 exact monotonicity and threshold certificates passed; Lam3 B<0 confined to [-149/10000,-3/1000].")


def main():
    checks = 0
    for eta in ((1,)*14, (1,)*7+(-1,)*7):
        algebra = Algebra(eta)
        epsilon = eta[7]
        variables, blocks = templates(epsilon)
        k, sigma, theta, delta = variables
        copies = {f'Lam{m}': algebra.forms(m, (0,2,4,6,8,10,12)[:m]) for m in (1,3,5,7)}
        copies.update({f'H{j}': algebra.hook(j) for j in (2,4,6)})
        count = 0
        for name, pats in copies.items():
            for op_index, op in enumerate(OPS):
                member = [int(j == op_index) for j in range(4)]
                x,y,z,w = member
                for lab, to in enumerate((algebra.phi, algebra.phig)):
                    t,u = int(lab == 0), int(lab == 1)
                    matrix = blocks[name].subs({k:0, sigma:t*x+epsilon*u*y, theta:t*w-u*z, delta:t*y-u*x})
                    for column, value in enumerate(pats):
                        actual = [scaled(item, sp.sympify(PHASE[op_index])/3) for item in algebra.mixed(op,value,to)]
                        expected = combine(pats, matrix[:,column])
                        assert normalized(actual) == normalized(expected), (epsilon,name,op,lab,column)
                        count += 1
        checks += count
        print(f"signature ({sum(a>0 for a in eta)},{sum(a<0 for a in eta)}): {count} independent block-column checks passed")

    positivity_certificate()

    # R9: real u(128) field k=w(omega)Gamma has real coefficients in grade 10.
    # The primed real member sharp has code coefficients (-i,0,-i,0), while To=i*t*Phi1.
    algebra = Algebra((1,)*14)
    pats = algebra.forms(3,(0,4,11))
    kernel = pats[2]
    t = F(-1,400)
    q = [dict() for _ in range(N)]
    even = [dict() for _ in range(N)]
    for op_index in (0,2):
        op = OPS[op_index]
        # Mixed term: (-i operator)*(i To) = 1; Q(k): (-i operator)*real k^2.
        mixed = [scaled(item,PHASE[op_index]*sp.Rational(t.numerator,t.denominator)/3)
                 for item in algebra.mixed(op,kernel,algebra.phi)]
        even = fplus(even,mixed)
        quadratic = [scaled(item,-sp.I*PHASE[op_index]/6)
                     for item in algebra.mixed(op,kernel,kernel)]
        q = fplus(q,quadratic)
    even = fplus(even,kernel)
    assert not any(normalized(even))
    norm = -algebra.pair(kernel,kernel)
    pairing = -algebra.pair([scaled(item,sp.I*sp.Rational(t.numerator,t.denominator)) for item in algebra.phi],q)
    assert norm > 0 and sp.expand(norm+2*pairing) == 0
    print(f"R9 kernel: kappa=1, t=-1/400; Lk=0; -Tr norm(k)^2={norm}; <To,Q(k)>={pairing}")
    assert any(normalized(q))
    current_checks = 0
    for direction in range(N):
        zeta = {(b,c): plus(scaled(kernel[c],int(b==direction)),scaled(kernel[b],-int(c==direction)))
                for b,c in combinations(range(N),2)}
        for op_index in (0,2):
            out = algebra.dual(OPS[op_index],zeta)
            assert sum(algebra.tr(algebra.phi[a],out[a]) for a in range(N)) == 0
            current_checks += 1
    print(f"R10 current: {current_checks} zero contractions (all 14 derivative directions, both sharp terms)")
    print("R10 polarization: kappa||k||^2 + 2<To,Q(k)> = 0 exactly; Q(k) is nonzero.")
    print("Conclusion uses the note's proof: fixed constant real To background, detM!=0, positive u(128) pairing; no boundary/integration assumptions.")
    print(f"Total: {checks} block-column checks, kernel/current/polarization checks passed; no repository imports.")


if __name__ == '__main__':
    main()
