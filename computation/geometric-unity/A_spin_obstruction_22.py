"""GU, Prop. A7 for a base of signature (2,2) (second-round review 2026-10-30, N5; reports/reviews/REVIEW2-2026-10-30.md).
Theorem A4 / note Theorem 1(1): a base of signature (2,2) gives a chimeric signature (7,7) for every lambda > 1/4 and
either block sign. A_spin_obstruction.py treats only the Lorentzian case. Here pi^*TX = P (+) N, with P and N the
y-positive and y-negative 2-planes of the tautological metric y.

S1  Exact Gram check at y = diag(1,1,-1,-1), c_V = 1, c_H = +1 and c_H = -1, several lambda > 1/4:
    V = Sym^2_0 P* (+) Sym^2_0 N* (+) (trace plane) (+) P* (x) N*, with signatures (2,0), (2,0), (1,1), (0,4);
    H* = P* (+) N*, with signatures (2,0), (0,2) for c_H = +1 and (0,2), (2,0) for c_H = -1. Hence
      c_H > 0:  C+ = Sym^2_0 P* + Sym^2_0 N* + R y^perp + P*,   C- = R y + P*(x)N* + N*;
      c_H < 0:  C+ = Sym^2_0 P* + Sym^2_0 N* + R y^perp + N*,   C- = R y + P*(x)N* + P*.
S2  Splitting principle, mod 2: w1(C+-) and o = w2(C+) + w2(C-) (the Karoubi obstruction, A_spin_obstruction.py K1).
    Expected: o = w2(X) + w1(N)^2 for c_H > 0 and o = w2(X) + w1(P)^2 for c_H < 0. The overall sign flip
    (c_V, c_H) -> (-c_V, -c_H) swaps C+ and C-, and o is symmetric.
S3  A family RP^2 -> (one fiber of Met_(2,2)): y_u = (2 u u^T - I_3) (+) (1), u a unit vector of R^3. Then P = R u + R e4
    and N = u^perp in R^3, so P = gamma (+) R and N = gamma^perp, with w1(P) = w1(N) = w1(gamma) =: a and a^2 != 0 in
    H^2(RP^2; Z2). pi^* w2(X) restricts to 0 on a fiber, so o != 0 on every fiber for both block signs:
    the chimeric bundle has no Spin(7,7) structure over Met_(2,2)(X), for any X.
Exact arithmetic throughout (sympy rationals, integer polynomials mod 2).
"""
import itertools
import sympy as sp


def mod2(expr, gens):
    p = sp.Poly(sp.expand(expr), *gens)
    return sp.Poly.from_dict({m: c % 2 for m, c in p.terms() if c % 2}, *gens) if p.terms() else sp.Poly(0, *gens)


def degree_part(expr, gens, d):
    p = mod2(expr, gens)
    return sum((c * sp.prod([v ** k for v, k in zip(gens, m)]) for m, c in p.terms() if sum(m) == d), sp.Integer(0))


def E(i, j):
    M = sp.zeros(4, 4)
    M[i, j] += 1
    M[j, i] += 1 if i != j else 0
    return M


def s1():
    print("=== S1. Exact Gram check of C+- at y = diag(1,1,-1,-1), c_V = 1 ===")
    y = sp.diag(1, 1, -1, -1)
    yi = y.inv()
    S0P = [E(0, 0) - E(1, 1), E(0, 1)]
    S0N = [E(2, 2) - E(3, 3), E(2, 3)]
    PN = [E(0, 2), E(0, 3), E(1, 2), E(1, 3)]
    yP, yN = E(0, 0) + E(1, 1), -(E(2, 2) + E(3, 3))       # trace plane; y = yP + yN
    ok = True
    for lam in [sp.Rational(3, 10), sp.Rational(1, 2), sp.Integer(1), sp.Integer(3)]:
        G = lambda h, k: (yi * h * yi * k).trace() - lam * (yi * h).trace() * (yi * k).trace()
        gram = lambda B: sp.Matrix([[G(a, b) for b in B] for a in B])
        g = yP + yN
        perp = G(g, g) * yP - G(g, yP) * g                    # G-orthogonal complement of y in the trace plane
        blocks = {'Sym^2_0 P*': S0P, 'Sym^2_0 N*': S0N, 'P*(x)N*': PN, 'R y': [g], 'R y^perp': [perp]}
        sig = {}
        for name, B in blocks.items():
            ev = gram(B).eigenvals()
            sig[name] = (sum(m for v, m in ev.items() if v > 0), sum(m for v, m in ev.items() if v < 0))
        cross = all(G(a, b) == 0 for (n1, B1), (n2, B2) in itertools.combinations(blocks.items(), 2)
                    for a in B1 for b in B2)
        expV = {'Sym^2_0 P*': (2, 0), 'Sym^2_0 N*': (2, 0), 'P*(x)N*': (0, 4), 'R y': (0, 1), 'R y^perp': (1, 0)}
        for cH in (1, -1):
            Psig = (2, 0) if cH * yi[0, 0] > 0 else (0, 2)
            Nsig = (2, 0) if cH * yi[2, 2] > 0 else (0, 2)
            tot = (sum(v[0] for v in sig.values()) + Psig[0] + Nsig[0], sum(v[1] for v in sig.values()) + Psig[1] + Nsig[1])
            expH = ((2, 0), (0, 2)) if cH == 1 else ((0, 2), (2, 0))
            good = sig == expV and cross and (Psig, Nsig) == expH and tot == (7, 7)
            ok &= good
            print(f"  lambda = {str(lam):5s} c_H = {cH:+d}: V-blocks {sig}; G-orthogonal: {cross}; "
                  f"H*: P* {Psig}, N* {Nsig}; total {tot}; as stated in the docstring: {good}")
    return ok


def s2():
    print("\n=== S2. Stiefel-Whitney classes of C+ and C- (splitting principle, mod 2) ===")
    p1, p2, n1, n2 = sp.symbols('p1 p2 n1 n2')      # w1 of the line summands of P and N (duals have equal classes)
    gens = (p1, p2, n1, n2)
    w_sym2P = 1 + p1 + p2                         # roots 2p1, 2p2, p1+p2; Sym^2_0 P (+) R = Sym^2 P
    w_sym2N = 1 + n1 + n2
    w_PN = sp.prod([1 + a + b for a in (p1, p2) for b in (n1, n2)])
    w_P, w_N = (1 + p1) * (1 + p2), (1 + n1) * (1 + n2)
    w_X = w_P * w_N
    w1P, w1N = degree_part(w_P, gens, 1), degree_part(w_N, gens, 1)
    w2X = degree_part(w_X, gens, 2)
    ok = True
    for cH, (hp, hm), t, tname in [(+1, (w_P, w_N), w1N, 'w1(N)'), (-1, (w_N, w_P), w1P, 'w1(P)')]:
        w_Cp = w_sym2P * w_sym2N * hp                 # Sym^2_0 P + Sym^2_0 N + R y^perp + (P* or N*)
        w_Cm = w_PN * hm                              # R y + P*(x)N* + (N* or P*)
        w1p, w1m = degree_part(w_Cp, gens, 1), degree_part(w_Cm, gens, 1)
        o = degree_part(w_Cp, gens, 2) + degree_part(w_Cm, gens, 2)
        okw1 = mod2(w1p - t, gens).is_zero and mod2(w1m - t, gens).is_zero
        oko = mod2(o - (w2X + t ** 2), gens).is_zero
        w2C = o + w1p * w1m
        okC = mod2(w2C - w2X, gens).is_zero
        ok &= okw1 and oko and okC
        print(f"  c_H = {cH:+d}: w1(C+) = {w1p}, w1(C-) = {w1m}; w1(C+) = w1(C-) = {tname}: {okw1}")
        print(f"           o = w2(C+) + w2(C-) == w2(X) + {tname}^2 (mod 2): {oko};  w2(C) == w2(X): {okC}")
    print("  On a fiber P (+) N is trivial, so w1(P) = w1(N) there and pi^* w2(X) = 0: o restricts to w1(N)^2 in both cases.")
    return ok


def s3():
    print("\n=== S3. An RP^2 inside one fiber of Met_(2,2) ===")
    print("  y_u = (2 u u^T - I_3) (+) (1); y_u = y_{-u}, so u -> y_u factors through RP^2.")
    ok = True
    e4 = sp.Matrix([0, 0, 0, 1])
    for u in [(1, 0, 0), (sp.Rational(3, 5), sp.Rational(4, 5), 0), (sp.Rational(2, 3), sp.Rational(2, 3), sp.Rational(1, 3)),
              (0, sp.Rational(5, 13), sp.Rational(12, 13))]:
        v3 = sp.Matrix(u)
        Y = sp.diag(2 * v3 * v3.T - sp.eye(3), sp.Matrix([[1]]))
        v = sp.Matrix(list(u) + [0])
        ev = Y.eigenvals()
        sig22 = ev == {1: 2, -1: 2}
        Pok = (Y * v == v) and (Y * e4 == e4)
        perp = [w for w in (sp.Matrix([1, 0, 0, 0]), sp.Matrix([0, 1, 0, 0]), sp.Matrix([0, 0, 1, 0]))]
        Nbasis = [w - (w.dot(v)) * v for w in perp]       # spans u^perp in R^3 (rank 2)
        Nok = all(Y * w == -w for w in Nbasis) and sp.Matrix.hstack(*Nbasis).rank() == 2
        same = sp.diag(2 * (-v3) * (-v3).T - sp.eye(3), sp.Matrix([[1]])) == Y
        ok &= sig22 and Pok and Nok and same
        print(f"  u = {tuple(u)}: eigenvalues {ev} (signature (2,2): {sig22}); P = Ru + Re4: {Pok}; N = u^perp in R^3: {Nok}")
    print("  So on this RP^2, P = gamma (+) R and N = gamma^perp: w1(P) = w1(N) = a, and a^2 != 0 in H^2(RP^2; Z2).")
    print("  Hence o = pi^* w2(X) + w1(N)^2 (or + w1(P)^2) != 0 on every fiber: no Spin(7,7) structure on C over")
    print("  Met_(2,2)(X), for any X and either block sign.")
    return ok


if __name__ == '__main__':
    r = [s1(), s2(), s3()]
    print(f"\nall checks pass: {all(r)}")
