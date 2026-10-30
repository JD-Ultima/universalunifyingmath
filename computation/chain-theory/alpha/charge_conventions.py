"""Convention dependence of the 'acting central charge' that selects the (2,3) link.

The Letter (Eq. 3, Lemma S-center) takes the identification algebra of a
(d1,d2) link to be s(u(d1) + u(d2)), traceless on the DIRECT SUM C^d1 + C^d2,
acting on the TENSOR PRODUCT C^d1 (x) C^d2 as xi1 (x) 1 + 1 (x) xi2.
Its centre (i x 1_d1, i y 1_d2), with d1 x + d2 y = 0, then acts on the
link as i(x+y) = i(d2-d1) (normalised x=d2): nonzero iff d1 != d2.

We evaluate the same question under four readings:
  A  (Letter)   trace on direct sum, far end acts by xi2
  A* (conj.)    trace on direct sum, far end acts by conj(xi2) (the far end
                transforms by U-bar under a matching, Eq. transpose identity)
  B             'no global phase' imposed on the matched directions: x + y = 0
  C             centre of u(d1)+u(d2) modulo the kernel of its action on the link
and, for each, which readable types (d1 d2 <= 6) carry an acting centre.
We also verify the stabiliser claim of Theorem S-stabilizer (dim 4 at (2,3)).

Run: python3 charge_conventions.py
"""
import itertools
import sympy as sp

x, y = sp.symbols('x y')

def act(d1, d2, conj):
    """scalar by which centre element (ix, iy) acts on C^d1 (x) C^d2 (divided by i)."""
    return x - y if conj else x + y

def acting(d1, d2, reading):
    if reading == 'A':
        sol = {y: -sp.Rational(d1, d2) * x}
        return sp.simplify(act(d1, d2, False).subs(sol)) != 0
    if reading == 'A*':
        # trace on direct sum of the *transformation* of each block.
        # Under conj on the far block the Lie-algebra element there is conj(xi2);
        # imposing the trace on (xi1, xi2) and acting by conj(xi2):
        sol = {y: -sp.Rational(d1, d2) * x}
        return sp.simplify(act(d1, d2, True).subs(sol)) != 0
    if reading == 'B':
        return False   # x + y = 0 identically kills the action i(x+y)
    if reading == 'C':
        return True    # (1,0) acts as i on every link
    raise ValueError

types = [(a, b) for a in range(2, 6) for b in range(a, 6)]
readable = [(a, b) for (a, b) in types if a * b <= 6]
print("types (d1<=d2<=5); readable = d1 d2 <= 6:", readable)
for R in ['A', 'A*', 'B', 'C']:
    S = [t for t in readable if acting(*t, R)]
    print(f"  reading {R:2s}: readable types with an acting central charge = {S}")

# hypercharge of the link under A and A*: Q = i diag(3,3;-2,-2,-2) = 6 i Y
for R, conj in [('A', False), ('A*', True)]:
    val = act(2, 3, conj).subs({x: 3, y: -2})
    print(f"  reading {R}: Q acts on the (2,3) link as {val}*i, i.e. 6Y = {val}, Y = {sp.Rational(val, 6)}")

# ---------- stabiliser of the pair ket at (2,3), both conventions ----------
def basis_u(d):
    B = []
    for j in range(d):
        for k in range(j + 1, d):
            m = sp.zeros(d); m[j, k] = 1; m[k, j] = -1; B.append(m)
            m = sp.zeros(d); m[j, k] = sp.I; m[k, j] = sp.I; B.append(m)
    for j in range(d):
        m = sp.zeros(d); m[j, j] = sp.I; B.append(m)
    return B

def stabiliser_dim(d1, d2, conj, traceless_sum=True):
    m = min(d1, d2)
    omega = sp.zeros(d1 * d2, 1)
    for i in range(m):
        omega[i * d2 + i] = 1
    gens = [(A, sp.zeros(d2)) for A in basis_u(d1)] + [(sp.zeros(d1), B) for B in basis_u(d2)]
    cs = sp.symbols(f'c0:{len(gens)}', real=True)
    xi1 = sum((c * g[0] for c, g in zip(cs, gens)), sp.zeros(d1))
    xi2 = sum((c * g[1] for c, g in zip(cs, gens)), sp.zeros(d2))
    far = xi2.conjugate() if conj else xi2
    L = sp.kronecker_product(xi1, sp.eye(d2)) + sp.kronecker_product(sp.eye(d1), far)
    v = L * omega
    eqs = []
    for e in v:
        e = sp.expand(e)
        eqs += [sp.re(e), sp.im(e)]
    if traceless_sum:
        t = sp.expand(xi1.trace() + xi2.trace())
        eqs.append(sp.im(t))
    Mx = sp.Matrix([[sp.diff(e, c) for c in cs] for e in eqs])
    alg_dim = len(cs) - (1 if traceless_sum else 0)
    return alg_dim, len(cs) - Mx.rank()

for conj in (False, True):
    a, s = stabiliser_dim(2, 3, conj)
    print(f"  (2,3), far end by {'conj(xi2)' if conj else 'xi2'}: algebra dim {a}, stabiliser of pair ket dim {s}")
for conj in (False, True):
    a, s = stabiliser_dim(2, 2, conj)
    print(f"  (2,2), far end by {'conj(xi2)' if conj else 'xi2'}: algebra dim {a}, stabiliser of pair ket dim {s}")
