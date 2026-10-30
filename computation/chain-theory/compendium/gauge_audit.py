#!/usr/bin/env python3
"""Audit of compendium Theorems 64.2 ("The gauge group as the face group of the tether") and
65.1 ("The gauge group as a corollary of one construction"); abstract: "We obtain the gauge
group SU(3) x SU(2) x U(1) ... as theorems".

Usage: python3 -I gauge_audit.py
Exact (sympy) where possible. Statements transcribed from
compendium/tex/parts/05-gauge-structure.tex lines 25-83 and the genre matrix,
compendium/tex/parts/01-framework.tex lines 52-56.
"""
import itertools, math
import sympy as sp

print(f"sympy {sp.__version__}")
GENRES = [(2, 2), (3, 2), (4, 2), (3, 3), (4, 3), (4, 4)]

def herm_basis(d):
    """basis of u(d) as i * Hermitian matrices (real span)"""
    B = []
    for j in range(d):
        M = sp.zeros(d); M[j, j] = sp.I; B.append(M)
    for j in range(d):
        for k in range(j+1, d):
            M = sp.zeros(d); M[j, k] = 1; M[k, j] = -1; B.append(M)
            M = sp.zeros(d); M[j, k] = sp.I; M[k, j] = sp.I; B.append(M)
    return B

def realvec(M):
    return [sp.re(x) for x in M] + [sp.im(x) for x in M]

def rank_real(mats):
    return sp.Matrix([realvec(M) for M in mats]).rank()

def kron(A, B):
    return sp.kronecker_product(A, B)

print("\n[1] Image of u(d1)+u(d2) acting on C^d1 (x) C^d2 by X(x)1 + 1(x)Y (the 'face group of the tether')")
print("| genre | dim image | d1^2+d2^2-1 | dim derived | rank (d1+d2-1) | central scalar of s(u+u) generator | kernel S(U x U)->U(d1 d2) | |ker SU x SU x U(1) -> S(U x U)| | readable d1*d2<=6 | acting centre |")
print("|---|---|---|---|---|---|---|---|---|---|")
results = {}
for d1, d2 in GENRES:
    I1, I2 = sp.eye(d1), sp.eye(d2)
    gens = [kron(X, I2) for X in herm_basis(d1)] + [kron(I1, Y) for Y in herm_basis(d2)]
    dim = rank_real(gens)
    # derived algebra: span of commutators of the su parts
    su1 = [kron(X, I2) for X in herm_basis(d1) if X.trace() == 0 or d1 == 1]
    su2 = [kron(I1, Y) for Y in herm_basis(d2) if Y.trace() == 0]
    comms = [A*B - B*A for A, B in itertools.product(gens, gens)]
    comms = [C for C in comms if C != sp.zeros(d1*d2)]
    ddim = rank_real(comms)
    rank = (d1 - 1) + (d2 - 1) + 1
    # s(u(d1)+u(d2)) centre generator (i d2 I_d1, -i d1 I_d2): action on tensor product
    Z = kron(sp.I*d2*I1, I2) + kron(I1, -sp.I*d1*I2)
    scal = sp.simplify(Z[0, 0]/sp.I)
    is_scalar = Z == Z[0, 0]*sp.eye(d1*d2)
    # kernel of S(U(d1)xU(d2)) -> U(d1 d2): (lam I, mu I), lam*mu = 1, lam^d1 mu^d2 = 1 -> lam^(d1-d2) = 1
    ker_tensor = "U(1)" if d1 == d2 else f"Z_{abs(d1-d2)}" if abs(d1-d2) > 1 else "trivial"
    ker_cover = math.lcm(d1, d2)
    readable = d1*d2 <= 6
    acting = d1 != d2
    results[(d1, d2)] = dict(dim=dim, ddim=ddim, rank=rank, scal=scal, readable=readable, acting=acting, kc=ker_cover)
    print(f"| ({d1},{d2}) | {dim} | {d1*d1+d2*d2-1} | {ddim} | {rank} | {scal} (scalar: {is_scalar}) | {ker_tensor} | {ker_cover} | {readable} | {acting} |")

PRINTED = {(2, 2): (7, 0), (3, 2): (12, 1), (4, 2): (19, 2), (3, 3): (17, 0), (4, 3): (24, 1), (4, 4): (31, 0)}
print("  genre matrix (01-framework.tex:52-56) prints (dimension, central scalar):", PRINTED)
print("  reproduced as (dim image, |central scalar|):",
      all(PRINTED[g] == (r['dim'], abs(r['scal'])) for g, r in results.items()))
print("\n  dim 12 and rank 4 single out: ",[g for g, r in results.items() if r['dim'] == 12 and r['rank'] == 4])
print("  readable AND acting centre:    ", [g for g, r in results.items() if r['readable'] and r['acting']])
print("  acting centre alone:           ", [g for g, r in results.items() if r['acting']])
print("  readable alone:                ", [g for g, r in results.items() if r['readable']])
for cap in (4, 6, 8, 9, 12):
    sel = [g for g in GENRES if g[0]*g[1] <= cap and g[0] != g[1]]
    print(f"  sensitivity: read cap d1*d2 <= {cap:>2} with acting centre selects {sel}")
for alph in ([2, 3], [2, 3, 4], [2, 3, 4, 5], list(range(2, 9))):
    gs = [(a, b) for a in alph for b in alph if a >= b]
    sel = [g for g in gs if g[0]*g[1] <= 6 and g[0] != g[1]]
    print(f"  sensitivity: alphabet {alph} (cap 6, acting centre) -> {len(gs)} genres, selects {sel}")

print("\n[2] Z_6 covering for (3,2): SU(2) x SU(3) x U(1) -> S(U(2) x U(3)), (A,B,z) -> (z^3 A, z^-2 B)")
w = sp.exp(2*sp.pi*sp.I/6)
ker = []
for k in range(6):
    z = w**k
    A = sp.simplify(z**-3) * sp.eye(2); B = sp.simplify(z**2) * sp.eye(3)
    ok = sp.simplify(A.det()) == 1 and sp.simplify(B.det()) == 1 and \
         sp.simplify(z**3*A - sp.eye(2)) == sp.zeros(2) and sp.simplify(z**-2*B - sp.eye(3)) == sp.zeros(3)
    ker.append(ok)
print(f"  z = exp(2 pi i k/6), k=0..5, with A = z^-3 I2 in SU(2), B = z^2 I3 in SU(3), maps to identity: {ker}")
# no other z: z^6 must be 1
zz = sp.symbols("zz")
print(f"  necessary condition from det: z^6 = 1 -> {len(sp.solve(zz**6 - 1, zz))} solutions; kernel = Z_6")
# map is a homomorphism into S(U(2)xU(3)): det(z^3 A) det(z^-2 B) = z^6 z^-6 = 1
print("  det(z^3 A) det(z^-2 B) = z^6 * z^-6 = 1 for A in SU(2), B in SU(3): image lies in S(U(2) x U(3))")

print("\n[3] The u(1) 'slice representative' (3,3,-2,-2,-2) printed in 65.1")
Yh = sp.diag(3, 3, -2, -2, -2)
blocks = [sp.diag(X, sp.zeros(3)) for X in herm_basis(2) if X.trace() == 0] + \
         [sp.diag(sp.zeros(2), Y) for Y in herm_basis(3) if Y.trace() == 0]
print(f"  traceless: {Yh.trace() == 0}; commutes with su(2)+su(3) block generators: {all(Yh*M == M*Yh for M in blocks)}")
print(f"  equals 6 x Y on the 5 of SU(5) = (1,2)_(1/2) + (3,1)_(-1/3), Y = diag(1/2,1/2,-1/3,-1/3,-1/3): "
      f"{sp.diag(3,3,-2,-2,-2) == 6*sp.diag(sp.Rational(1,2), sp.Rational(1,2), sp.Rational(-1,3), sp.Rational(-1,3), sp.Rational(-1,3))}")

print("\n[4] 'fourth-root monomial groups at the orders four, thirty-two, three hundred eighty-four and six thousand one hundred forty-four'")
def monomial_group_order(n):
    roots = [1, sp.I, -1, -sp.I]
    count = 0
    for perm in itertools.permutations(range(n)):
        for phases in itertools.product(range(4), repeat=n):
            count += 1
    return count
orders = [monomial_group_order(n) for n in range(1, 5)]
print(f"  |mu_4 wr S_n| by enumeration, n=1..4: {orders}; formula 4^n n!: {[4**n*math.factorial(n) for n in range(1,5)]}; "
      f"printed (4, 32, 384, 6144): {orders == [4, 32, 384, 6144]}")

print("\n[5] Z_6 compatibility of the Standard Model spectrum (Q = T3 + Y): phase exp(2 pi i (t/3 + d/2 + Y)) must be 1")
SM = [("Q_L", 1, 1, sp.Rational(1, 6)), ("u_R", 1, 0, sp.Rational(2, 3)), ("d_R", 1, 0, sp.Rational(-1, 3)),
      ("L_L", 0, 1, sp.Rational(-1, 2)), ("e_R", 0, 0, -1), ("nu_R?", 0, 0, 0), ("H", 0, 1, sp.Rational(1, 2)),
      ("gluon", 0, 0, 0), ("W", 0, 0, 0)]
for name, t, d, Y in SM:
    ph = sp.Rational(t, 3) + sp.Rational(d, 2) + Y
    print(f"  {name:6s} triality {t} duality {d} Y {str(Y):>5}: t/3+d/2+Y = {ph} -> trivial: {ph.is_integer}")
print("  Every SM field is Z_6-neutral, so (SU3 x SU2 x U1)/Z_6 acts faithfully; the quotients by Z_1, Z_2, Z_3 are")
print("  equally consistent with the perturbative spectrum (Tong, JHEP 07 (2017) 104): the global form is not yet measured.")

print("\n[6] Discrete choices on which the identification depends (count):")
choices = ["dimension alphabet {2,3,4} (inherited from the 64 listed webs; 9.1)",
           "read cap value d1*d2 <= 6 (9.1, CONDITIONAL)",
           "'acting centre' criterion d1 != d2 (9.2)",
           "assignment of SU(2) to hub-leg qubit copies (65.1)",
           "assignment of SU(3) to qutrit ends open at the hub leg (65.1)",
           "assignment of U(1) to the dock qubit end slice (65.1)"]
for c in choices: print("   -", c)
print(f"  {len(choices)} discrete choices; the Lie-algebra output (dims 8,3,1; centre Z_6) is the standard structure of S(U(2)xU(3)).")
