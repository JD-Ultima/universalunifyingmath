"""Independent coordinator cross-check of Theorem D1 (dim S = 200 for g = u(128)), by a method
different from D_rep_counts.py.

As so(14,C)-modules, g (x) C = End(S_Dirac) = Cl(14) = Lambda^* V, and Lambda^13 V = Lambda^1 V. Hence
    dim_C Hom_{so(14)}(Lambda^2 (x) g, Lambda^13 (x) g) = mult of the trivial rep in
    Lambda^2 V (x) V (x) Lambda^* V (x) Lambda^* V      (Lambda^* is self-dual).
The multiplicity of the trivial representation is computed by the Brauer/Racah-Speiser formula
    mult_0(rho) = sum_{w in W(D7)} eps(w) m_rho(w(delta) - delta),  delta = (6,5,4,3,2,1,0),
with weight multiplicities m_rho from the torus character. Lambda^* V has character
prod_i (1+x_i)(1+1/x_i), so the Lambda^*(x)Lambda^* factor is prod_i (x_i^-1 + 2 + x_i)^2.
Sanity checks: mult_0(V(x)V)=1, mult_0(Lambda^2(x)Lambda^2)=1, dim End_{so}(Lambda^*)=30.
Run: python3 xcheck_shiab_dim.py
"""
import itertools
from collections import defaultdict
n = 7
basis = []
for i in range(n):
    for s in (1, -1):
        v = [0] * n; v[i] = s; basis.append(tuple(v))
def wedge(k):
    d = defaultdict(int)
    for c in itertools.combinations(range(2 * n), k):
        w = [0] * n
        for j in c:
            for t in range(n): w[t] += basis[j][t]
        d[tuple(w)] += 1
    return d
def tensor(a, b):
    d = defaultdict(int)
    for x, m in a.items():
        for y, l in b.items():
            d[tuple(p + q for p, q in zip(x, y))] += m * l
    return d
# per-coordinate factor polynomials
f1 = {-1: 1, 0: 2, 1: 1}                        # Lambda^* restricted to one coordinate pair
f2 = {-2: 1, -1: 4, 0: 6, 1: 4, 2: 1}           # (Lambda^* (x) Lambda^*) per coordinate
def mult_fn(finite, per_coord):
    def m(mu):
        tot = 0
        for nu, c in finite.items():
            p = c
            for t in range(n):
                p *= per_coord.get(mu[t] - nu[t], 0)
                if p == 0: break
            tot += p
        return tot
    return m
delta = (6, 5, 4, 3, 2, 1, 0)
def perm_sign(p):
    s, seen = 1, [False] * len(p)
    for i in range(len(p)):
        if not seen[i]:
            j, L = i, 0
            while not seen[j]:
                seen[j] = True; j = p[j]; L += 1
            if L % 2 == 0: s = -s
    return s
W = []
for p in itertools.permutations(range(n)):
    sp = perm_sign(p)
    for signs in itertools.product((1, -1), repeat=n):
        if signs.count(-1) % 2: continue
        wd = tuple(signs[i] * delta[p[i]] for i in range(n))
        W.append((sp, tuple(a - b for a, b in zip(wd, delta))))
assert len(W) == 322560
def mult0(m, bound):
    return sum(e * m(mu) for e, mu in W if max(abs(x) for x in mu) <= bound)
V = wedge(1); L2 = wedge(2)
unit = {(0,) * n: 1}
print("sanity  mult_0(V x V)        =", mult0(mult_fn(tensor(V, V), {0: 1}), 2))
print("sanity  mult_0(L2 x L2)      =", mult0(mult_fn(tensor(L2, L2), {0: 1}), 4))
print("sanity  dim End_so(Lambda^*) =", mult0(mult_fn(unit, f2), 2), "(expected 30: 2*(L0..L6) + L7+ + L7-)")
dimS = mult0(mult_fn(tensor(L2, V), f2), 5)
print("dim_C Hom_so(14)(Lambda^2 (x) gl(128), Lambda^13 (x) gl(128)) =", dimS, "(Theorem D1: 200)")

# so(64,64) case. so(S, C) for the symmetric invariant form C_- of Cl(7,7) consists of the Clifford grades k with
# C gamma^(k) antisymmetric: sign(k) = (-1)^k (-1)^{k(k-1)/2} (C_- symmetric)  ->  k in {1,2,5,6,9,10,13,14}.
grades = [k for k in range(15) if (-1) ** k * (-1) ** (k * (k - 1) // 2) == -1]
from math import comb
assert sum(comb(14, k) for k in grades) == 128 * 127 // 2, grades
g = defaultdict(int)
for k in grades:
    for wgt, c in wedge(k).items():
        g[wgt] += c
gg = tensor(g, g)
L2V = tensor(L2, V)
def m_so(mu):
    return sum(c * gg.get(tuple(a - b for a, b in zip(mu, nu)), 0) for nu, c in L2V.items())
dim_so = sum(e * m_so(mu) for e, mu in W if max(abs(x) for x in mu) <= 7)
print("so(64,64): Clifford grades", grades, "; dim_C Hom_so(14)(Lambda^2 (x) g, Lambda^13 (x) g) =", dim_so, "(Theorem D1: 52)")

# Self-adjoint part (D2): circledast = Psi ^ with Psi in (Lambda^3 (x) S^2 g)^{so(14)}  (assuming Psi -> Psi^ injective).
# chi_{S^2 g}(x) = (chi_g(x)^2 + chi_g(x^2)) / 2.
L3 = wedge(3)
gsq_u = {-2: 1, 0: 2, 2: 1}                      # Lambda^* at x -> x^2, per coordinate
m_a = mult_fn(L3, f2); m_b = mult_fn(L3, gsq_u)
sa_u = sum(e * (m_a(mu) + m_b(mu)) for e, mu in W if max(abs(x) for x in mu) <= 5) // 2
print("u(64,64): mult_0(Lambda^3 (x) S^2 g) =", sa_u, "(Theorem D1: dim S_sa = 48)")
g2 = defaultdict(int)
for wgt, c in g.items():
    g2[tuple(2 * x for x in wgt)] += c
L3 = dict(L3)
def m_so_sa(mu):
    a = sum(c * gg.get(tuple(p - q for p, q in zip(mu, nu)), 0) for nu, c in L3.items())
    b = sum(c * g2.get(tuple(p - q for p, q in zip(mu, nu)), 0) for nu, c in L3.items())
    return a + b
sa_so = sum(e * m_so_sa(mu) for e, mu in W if max(abs(x) for x in mu) <= 7) // 2
print("so(64,64): mult_0(Lambda^3 (x) S^2 g) =", sa_so, "(Theorem D1: dim S_sa = 12)")
