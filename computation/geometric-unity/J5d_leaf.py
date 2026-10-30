"""GU round 6, follow-up: the indefinite case of Theorem J5(d) (the leaf condition of the radial Z2 ansatz).

Formalization: formalization/geometric-unity/J5d-indefinite-leaf.md.  Builds on J5 (J-slaved-sector.md), F3 and L-F5
(F-zorro-reduction.md), A4/A5 (A-chimeric-signature.md).  Imports no repository code.

Setting.  Base (R^n, g), g = diag(eta), signature (p, q); e_a the coordinate (orthonormal) frame, theta^a its dual.
Curvature convention of F: R_abcd = g(R(e_c, e_d) e_b, e_a), Ric_bd = g^ac R_abcd.  An endomorphism B acts on symmetric
bilinear forms by  B.y = -(B^T y + y B)  (the induced action on Sym^2 T*X, as in L-F5 and J_slaved_sector.py section J5).
With A_X Y = -1/2 (R(X,Y).y)^vert (L-F5) and the Z2 fibre metric c_V G_lambda, the O'Neill sum of J5(c)-(d) at (x0, y) is
   (A_X, A_Y) = c_V/(4 c_H) * Phi_y(X, Y),   Phi_y(X, Y) = sum_ij g^ij tr( y^-1 (R(X,e_i).y) y^-1 (R(Y,e_j).y) ),
and its order-eps^2 part at y = g + eps h is the quadratic form of J5(d),
   F(X, Y; h) = sum_ij g^ij G_0(R(X,e_i).h, R(Y,e_j).h),   G_0(h, k) = tr(g^-1 h g^-1 k).
The leaf condition gives F(X, Y; h, h) = 0 for all X, Y, h (Lemma L0).  Its polarisation in h only gives the part of
the bilinear form F(X, Y; h, h') that is symmetric under h <-> h' (equivalently X <-> Y); F itself is not symmetric.
Contractions (theta^p.theta^q = (theta^p (x) theta^q + theta^q (x) theta^p)/2, h_aq = theta^a.theta^q):
   C1     = sum_{a,p,q} eta_a eta_p eta_q F(e_a, e_a; theta^p.theta^q, theta^p.theta^q)     (both traces; quadratic values)
   C2sym  = sum_{a,b,q} eta_q (1/2)[Q_ab(h_aq + h_bq) - Q_ab(h_aq) - Q_ab(h_bq)],  Q_ab(h) = F(e_a, e_b; h, h)
   C2, C2' = sum_{a,b,q} eta_q F(e_a, e_b; h_aq, h_bq), resp. F(e_a, e_b; h_bq, h_aq)   (bilinear; C2sym = (C2 + C2')/2)
Only C1 and C2sym are forced to vanish by the leaf condition.  (Corrected after review 2026-10-30, R7-1: the first
version used C2, which is not; see reports/CORRECTIONS.md.)
Sections (labels as the results of the note)
  L0  Lemma L0, the reduction of J5(d), re-derived: R.y is y-traceless (the lambda-term drops at every y); the closed
      form Phi_y = 2 sum g^ij [tr(B_i C_j) - tr(M^-1 B_i M C_j)], M = g^-1 y; Phi_{g+eps h} = eps^2 F(X,Y;h,h) + O(eps^3);
      F(h + c g) = F(h) and F(X,Y;h,h) symmetric in X, Y; the constant-curvature polynomials of J_slaved_sector.out
      section J5 reproduced in the same parametrisation.
  L1  Theorem L1, for generic algebraic curvature tensors (n = 3, 4, 5; nine signatures): C1 = (n+2)|R|^2 and
      C2sym = (n+8)/8 |R|^2 + (n+4)/4 |Ric|^2, both from values of the quadratic form; hence
      8(n+2) C2sym - (n+8) C1 = 2(n+2)(n+4)|Ric|^2, i.e. 4 C2sym - C1 = 8|Ric|^2 at n = 4.  Side identities (not
      consequences of the leaf condition): C2 = 3/4 |R|^2 + (n+3)/2 |Ric|^2, C2' = (n+5)/4 |R|^2 + 1/2 |Ric|^2,
      C2 != C2' (expected), and the antisymmetric part F(X,Y;h,h') - F(X,Y;h',h) = -sum_i eta_i tr([B_i, C_i][H', H]).
  L2  Corollary L2: Einstein R = W + (Lambda/3)(g_ac g_bd - g_ad g_bc), n = 4, four signatures: C1 = 6|W|^2 +
      16 Lambda^2, C2sym = 3/2 |W|^2 + 12 Lambda^2, 4 C2sym - C1 = 32 Lambda^2; with J5(b), Lambda_g = -c_H/c_V.
  L3  Proposition L3, positivity for a Lorentzian (or definite) base: F(T, T; T^b (x) T^b) = 2 g(T,T) sum_{i,j}
      eta_i eta_j R(T,e_i,T,e_j)^2 (i, j over an orthonormal frame of T^perp), generic R; R(T,X,T,Y) determines R.
  L4  Proposition L4, from the coefficients of the quadratic form F(X,X;h,h) (n = 4, Einstein): their span has rank 47
      of 66 and contains Lambda^2 in (1,3), (2,2), (4,0) ((3,1), (0,4) follow by g -> -g, which sends F -> -F).  Base
      (2,2): the curvature operator commutes with *, the 81 products of the blocks R+ and R- lie in the span; on R+ = 0
      (resp. R- = 0) the span equals the span of the entries of (R-)^2 (resp. (R+)^2); the null example
      R = n (x) n (n an anti-self-dual null 2-form) has F = 0 and satisfies the full condition det(y)^2 Phi_y = 0
      identically in y.
Exact arithmetic (sympy rationals).  Fixed seed random.Random(20261030).  One core.  Every check is asserted.
"""
import itertools
import os
import random
import sys
import time

for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import sympy as sp  # noqa: E402

T0 = time.time()
RNG = random.Random(20261030)
RES = []
Lam = sp.Symbol('Lambda')
SIG4 = {'(4,0)': (1, 1, 1, 1), '(1,3)': (1, -1, -1, -1), '(3,1)': (1, 1, 1, -1), '(2,2)': (1, 1, -1, -1)}


def say(*a):
    print(*a)
    sys.stdout.flush()


def check(label, ok):
    ok = bool(ok)
    RES.append((label, ok))
    say(f"  [{ok}] {label}")
    assert ok, label


def rat(lo=-5, hi=5):
    return sp.Rational(RNG.randint(lo, hi), RNG.randint(1, 4))


# ------------------------------------------------------------------------------------------ algebraic curvature tensors
def curvature_family(n, sig, einstein=False):
    """generic algebraic curvature tensor on (R^n, diag(sig)): pair-symmetric components, first Bianchi solved over Q;
    with einstein=True also Ric = Lambda g (Lambda the symbol Lam).  Returns (R dict, free parameters)."""
    g = sp.diag(*sig)
    gi = g.inv()
    pairs = [(a, b) for a in range(n) for b in range(a + 1, n)]
    pidx = {P: i for i, P in enumerate(pairs)}
    comp, syms = {}, []
    for i in range(len(pairs)):
        for j in range(i, len(pairs)):
            s = sp.Symbol(f"r{i}_{j}")
            syms.append(s)
            comp[(i, j)] = comp[(j, i)] = s

    def c(a, b, cc, d):
        if a == b or cc == d:
            return sp.S(0)
        sg = 1
        if a > b:
            a, b, sg = b, a, -sg
        if cc > d:
            cc, d, sg = d, cc, -sg
        return sg * comp[(pidx[(a, b)], pidx[(cc, d)])]
    R = {(a, b, cc, d): c(a, b, cc, d) for a, b, cc, d in itertools.product(range(n), repeat=4)}
    eqs = [R[a, b, cc, d] + R[a, cc, d, b] + R[a, d, b, cc] for a, b, cc, d in itertools.product(range(n), repeat=4)]
    eqs = [e for e in eqs if e != 0]
    sol = sp.solve(eqs, syms, dict=True)[0] if eqs else {}
    R = {key: sp.expand(v.subs(sol)) for key, v in R.items()}
    free = sorted(set().union(*[v.free_symbols for v in R.values()]), key=str)
    if einstein:
        Ric = ricci(R, n, gi)
        es = [sp.expand(Ric[b, d] - Lam * g[b, d]) for b in range(n) for d in range(b, n)]
        sol = sp.solve(es, free, dict=True)[0]
        R = {key: sp.expand(v.subs(sol)) for key, v in R.items()}
        free = sorted(set().union(*[v.free_symbols for v in R.values()]) - {Lam}, key=str)
    return R, free


def ricci(R, n, gi):
    return sp.Matrix(n, n, lambda b, d: sum(gi[a, c] * R[a, b, c, d] for a in range(n) for c in range(n)))


def norm2(R, n, sig):
    return sp.expand(sum(R[a, b, c, d] ** 2 * sig[a] * sig[b] * sig[c] * sig[d]
                         for a, b, c, d in itertools.product(range(n), repeat=4)))


def norm2_ric(R, n, sig):
    Ric = ricci(R, n, sp.diag(*sig).inv())
    return sp.expand(sum(Ric[a, b] ** 2 * sig[a] * sig[b] for a in range(n) for b in range(n)))


def const_curv(n, sig):
    """R_1: R_1(X,Y)Z = g(Y,Z)X - g(X,Z)Y, i.e. (R_1)_abcd = g_ac g_bd - g_ad g_bc; Ric = (n-1) g"""
    g = sp.diag(*sig)
    return {(a, b, c, d): g[a, c] * g[b, d] - g[a, d] * g[b, c] for a, b, c, d in itertools.product(range(n), repeat=4)}


def rend(R, n, gi):
    """the endomorphisms R(e_c, e_d): (R(e_c,e_d))^p_q = g^pe R_eqcd"""
    return {(c, d): sp.Matrix(n, n, lambda p, q: sum(gi[p, e] * R[e, q, c, d] for e in range(n)))
            for c in range(n) for d in range(n)}


def act(B, y):
    return -(B.T * y + y * B)


def sym_basis(n, p, q):
    """theta^p . theta^q as a symmetric matrix"""
    return sp.Matrix(n, n, lambda i, j: sp.Rational(int((i, j) == (p, q)) + int((i, j) == (q, p)), 2))


def sym_generic(n, name):
    s = sp.symbols(f"{name}0:{n * (n + 1) // 2}")
    h = sp.zeros(n, n)
    it = iter(s)
    for i in range(n):
        for j in range(i, n):
            h[i, j] = h[j, i] = next(it)
    return h, s


def F_value(Rm, n, gi, X, Y, h, hp):
    BX = [sum((X[a] * Rm[(a, i)] for a in range(n)), sp.zeros(n, n)) for i in range(n)]
    BY = [sum((Y[a] * Rm[(a, i)] for a in range(n)), sp.zeros(n, n)) for i in range(n)]
    return sp.expand(sum(gi[i, j] * (gi * act(BX[i], h) * gi * act(BY[j], hp)).trace()
                         for i in range(n) for j in range(n) if gi[i, j] != 0))


def Phi_value(Rm, n, gi, X, Y, y):
    BX = [sum((X[a] * Rm[(a, i)] for a in range(n)), sp.zeros(n, n)) for i in range(n)]
    BY = [sum((Y[a] * Rm[(a, i)] for a in range(n)), sp.zeros(n, n)) for i in range(n)]
    yi = y.inv()
    return sp.expand(sum(gi[i, j] * (yi * act(BX[i], y) * yi * act(BY[j], y)).trace()
                         for i in range(n) for j in range(n) if gi[i, j] != 0)), BX, BY


def contractions(R, n, sig):
    """C1 and C2sym from values of the quadratic form Q_ab(h) = F(e_a, e_b; h, h) only; C2, C2' from the bilinear form"""
    g = sp.diag(*sig)
    gi = g.inv()
    Rm = rend(R, n, gi)
    Sb = {(p, q): sym_basis(n, p, q) for p in range(n) for q in range(n)}
    A = {(a, i, p, q): gi * act(Rm[(a, i)], Sb[(p, q)]) for a in range(n) for i in range(n) for p in range(n)
         for q in range(n)}

    def Q(a, b, hs):
        """Q_ab(h) = F(e_a, e_b; h, h) for h = sum of the basis forms theta^p.theta^q, (p, q) in hs (A is linear in h)"""
        Aa = [sum((A[(a, i) + pq] for pq in hs), sp.zeros(n, n)) for i in range(n)]
        Ab = [sum((A[(b, i) + pq] for pq in hs), sp.zeros(n, n)) for i in range(n)]
        return sum(gi[i, i] * (Aa[i] * Ab[i]).trace() for i in range(n))

    def Fb(a, b, pq, rs):   # bilinear F(e_a, e_b; theta^p.theta^q, theta^r.theta^s), g diagonal
        return sum(gi[i, i] * (A[(a, i) + pq] * A[(b, i) + rs]).trace() for i in range(n))
    trip = [(a, b, q) for a in range(n) for b in range(n) for q in range(n)]
    C1 = sp.expand(sum(gi[a, a] * g[p, p] * g[q, q] * Q(a, a, [(p, q)])
                       for a in range(n) for p in range(n) for q in range(n)))
    C2s = sp.expand(sum(g[q, q] * sp.Rational(1, 2) * (Q(a, b, [(a, q), (b, q)]) - Q(a, b, [(a, q)]) - Q(a, b, [(b, q)]))
                        for a, b, q in trip))
    C2 = sp.expand(sum(g[q, q] * Fb(a, b, (a, q), (b, q)) for a, b, q in trip))
    C2p = sp.expand(sum(g[q, q] * Fb(a, b, (b, q), (a, q)) for a, b, q in trip))
    return C1, C2s, C2, C2p


def random_point(R, free):
    return {s: rat() for s in free}


def random_sym_invertible(n):
    while True:
        y = sp.Matrix(n, n, lambda i, j: 0)
        for i in range(n):
            for j in range(i, n):
                y[i, j] = y[j, i] = rat()
        if y.det() != 0:
            return y


# ================================================================================================ L0
def L0():
    say("\n=== L0. Lemma L0: the reduction of J5(d), re-derived (n = 4) ===")
    n = 4
    for lab, sig in SIG4.items():
        g = sp.diag(*sig)
        gi = g.inv()
        say(f"  base {lab}:")
        # (a) lambda-term: tr(y^-1 (B.y)) = -2 tr B, = 0 for B in so(g)
        ok_a = True
        for _ in range(3):
            y = random_sym_invertible(n)
            Bgen = sp.Matrix(n, n, lambda i, j: rat())
            K = sp.Matrix(n, n, lambda i, j: 0)
            for i in range(n):
                for j in range(i + 1, n):
                    K[i, j] = rat()
                    K[j, i] = -K[i, j]
            Bso = gi * K
            ok_a &= (y.inv() * act(Bgen, y)).trace() == -2 * Bgen.trace()
            ok_a &= (Bso.T * g + g * Bso) == sp.zeros(n, n) and (y.inv() * act(Bso, y)).trace() == 0
        check(f"{lab} tr(y^-1 (B.y)) = -2 tr B, = 0 for B in so(g), at 3 random rational y: the lambda-term of G_lambda"
              f" drops at every y, and A_X Y is tangent to the unimodular slice", ok_a)
        # (b) closed form of Phi_y with M = g^-1 y, random curvature tensor (not Einstein), random y, X, Y
        Rg, free = curvature_family(n, sig)
        pt = random_point(Rg, free)
        R = {k_: v.subs(pt) for k_, v in Rg.items()}
        Rm = rend(R, n, gi)
        ok_b = True
        for _ in range(3):
            y = random_sym_invertible(n)
            X = [rat() for _ in range(n)]
            Y = [rat() for _ in range(n)]
            Phi, BX, BY = Phi_value(Rm, n, gi, X, Y, y)
            M = gi * y
            closed = 2 * sum(gi[i, i] * ((BX[i] * BY[i]).trace() - (M.inv() * BX[i] * M * BY[i]).trace()) for i in range(n))
            Phi2, _, _ = Phi_value(Rm, n, gi, X, Y, sp.Rational(7, 3) * y)
            ok_b &= sp.simplify(Phi - closed) == 0 and Phi2 == Phi
        check(f"{lab} Phi_y(X,Y) = 2 sum_ij g^ij [tr(B_i C_j) - tr(M^-1 B_i M C_j)], M = g^-1 y, B_i = R(X,e_i),"
              f" C_j = R(Y,e_j); Phi_{{cy}} = Phi_y (3 random rational y, X, Y; random curvature tensor)", ok_b)
        # (c) Taylor expansion at y = g + eps h, h g-traceless
        eps = sp.Symbol('epsilon')
        h = sp.Matrix(n, n, lambda i, j: 0)
        for i in range(n):
            for j in range(i, n):
                h[i, j] = h[j, i] = rat()
        h[0, 0] = h[0, 0] - (gi * h).trace() * g[0, 0]
        assert (gi * h).trace() == 0
        X = [rat() for _ in range(n)]
        Y = [rat() for _ in range(n)]
        y = g + eps * h
        adj = y.adjugate()
        BX = [sum((X[a] * Rm[(a, i)] for a in range(n)), sp.zeros(n, n)) for i in range(n)]
        BY = [sum((Y[a] * Rm[(a, i)] for a in range(n)), sp.zeros(n, n)) for i in range(n)]
        P = sp.Poly(sp.expand(sum(gi[i, i] * (adj * act(BX[i], y) * adj * act(BY[i], y)).trace() for i in range(n))), eps)
        Fv = F_value(Rm, n, gi, X, Y, h, h)
        cs = [P.coeff_monomial(eps ** m) for m in range(3)]
        check(f"{lab} det(y)^2 Phi_y at y = g + eps h (h g-traceless) = det(g)^2 F(X,Y;h) eps^2 + O(eps^3): eps^0, eps^1"
              f" coefficients 0, eps^2 coefficient = F(X,Y;h,h) (random rational data; value {Fv})",
              cs[0] == 0 and cs[1] == 0 and cs[2] == g.det() ** 2 * Fv)
        # (d) F(h + c g) = F(h); F symmetric in X, Y
        cc = sp.Symbol('c')
        hh, _ = sym_generic(n, 'h')
        Xs = sp.symbols('X0:4')
        Ys = sp.symbols('Y0:4')
        Fa = F_value(Rm, n, gi, Xs, Ys, hh + cc * g, hh + cc * g)
        Fb = F_value(Rm, n, gi, Xs, Ys, hh, hh)
        Fc_ = F_value(Rm, n, gi, Ys, Xs, hh, hh)
        check(f"{lab} R.g = 0, so F(X,Y; h + c g, h + c g) = F(X,Y; h, h) (the traceless restriction is immaterial), and"
              f" F(X,Y;h,h) = F(Y,X;h,h) (symbolic h, X, Y, c)", sp.expand(Fa - Fb) == 0 and sp.expand(Fb - Fc_) == 0)
    # (e) reproduction of J_slaved_sector.out section J5 (constant curvature +1, h0 eliminated by the trace)
    printed = {
        (1, 1, 1, 1): "4*(3*h1**2 + 3*h2**2 + 3*h3**2 + 3*h4**2 + 5*h4*h7 + 5*h4*h9 + h5**2 + h6**2 + 3*h7**2 + 5*h7*h9"
                      " + h8**2 + 3*h9**2)",
        (1, -1, -1, -1): "-4*(3*h1**2 + 3*h2**2 + 3*h3**2 - 3*h4**2 - 5*h4*h7 - 5*h4*h9 - h5**2 - h6**2 - 3*h7**2"
                         " - 5*h7*h9 - h8**2 - 3*h9**2)",
        (1, 1, -1, -1): "4*(3*h1**2 - 3*h2**2 - 3*h3**2 + 3*h4**2 - 5*h4*h7 - 5*h4*h9 - h5**2 - h6**2 + 3*h7**2"
                        " + 5*h7*h9 + h8**2 + 3*h9**2)"}
    ok_e = True
    for sig, txt in printed.items():
        g = sp.diag(*sig)
        gi = g.inv()
        Rm = rend(const_curv(n, sig), n, gi)
        h, hs = sym_generic(n, 'h')
        h = h.subs(hs[0], sp.solve((gi * h).trace(), hs[0])[0])
        e0 = [1, 0, 0, 0]
        ok_e &= sp.expand(F_value(Rm, n, gi, e0, e0, h, h) - sp.sympify(txt, locals={f"h{i}": hs[i] for i in range(10)})) == 0
    check("same conventions as J_slaved_sector.py section J5: F(e0,e0;h) for constant curvature +1 equals the three"
          " polynomials printed in J_slaved_sector.out (bases (4,0), (1,3), (2,2))", ok_e)


# ================================================================================================ L1
def L1():
    say("\n=== L1. Theorem L1: traces of F for generic algebraic curvature tensors ===")
    say("  C1 and C2sym are built from values Q_ab(h) = F(e_a, e_b; h, h) of the quadratic form only, which the leaf")
    say("  condition forces to vanish; C2 and C2' use the bilinear form at h != h' and are side identities.")
    cases = ((3, ((1, 1, 1), (1, -1, -1))),
             (4, tuple(SIG4.values())),
             (5, ((1, 1, 1, 1, 1), (1, -1, -1, -1, -1), (1, 1, -1, -1, -1))))
    for n, sigs in cases:
        for sig in sigs:
            R, free = curvature_family(n, sig)
            C1, C2s, C2, C2p = contractions(R, n, sig)
            R2, Ric2 = norm2(R, n, sig), norm2_ric(R, n, sig)
            ok1 = sp.expand(C1 - (n + 2) * R2) == 0
            ok2 = sp.expand(C2s - sp.Rational(n + 8, 8) * R2 - sp.Rational(n + 4, 4) * Ric2) == 0
            ok3 = sp.expand(8 * (n + 2) * C2s - (n + 8) * C1 - 2 * (n + 2) * (n + 4) * Ric2) == 0
            check(f"n = {n}, signature {sig.count(1)},{sig.count(-1)} ({len(free)} free components): C1 = (n+2)|R|^2,"
                  f" C2sym = (n+8)/8 |R|^2 + (n+4)/4 |Ric|^2, hence 8(n+2) C2sym - (n+8) C1 = 2(n+2)(n+4)|Ric|^2"
                  f"  [{time.time() - T0:.0f} s]", ok1 and ok2 and ok3)
            side = (sp.expand(C2 - sp.Rational(3, 4) * R2 - sp.Rational(n + 3, 2) * Ric2) == 0
                    and sp.expand(C2p - sp.Rational(n + 5, 4) * R2 - sp.Rational(1, 2) * Ric2) == 0
                    and sp.expand(C2s - (C2 + C2p) / 2) == 0
                    and sp.expand(C2 - C2p - sp.Rational(n + 2, 2) * Ric2 + sp.Rational(n + 2, 4) * R2) == 0)
            check(f"   side identities (bilinear form, NOT consequences of the leaf condition): C2 = 3/4 |R|^2 +"
                  f" (n+3)/2 |Ric|^2, C2' = (n+5)/4 |R|^2 + 1/2 |Ric|^2, C2sym = (C2 + C2')/2, and the expected asymmetry"
                  f" C2 - C2' = (n+2)/2 |Ric|^2 - (n+2)/4 |R|^2, a nonzero polynomial: {sp.expand(C2 - C2p) != 0}",
                  side and sp.expand(C2 - C2p) != 0)
    say("  n = 4: 4 C2sym - C1 = 8 |Ric|^2 for every algebraic curvature tensor, every signature.")
    # the antisymmetric part of the bilinear form (Jacobi), at random rational data
    n = 4
    for lab, sig in SIG4.items():
        g = sp.diag(*sig)
        gi = g.inv()
        Rg, free = curvature_family(n, sig)
        pt = random_point(Rg, free)
        Rm = rend({k_: v.subs(pt) for k_, v in Rg.items()}, n, gi)
        X = [rat() for _ in range(n)]
        Y = [rat() for _ in range(n)]
        h, hp = random_sym_invertible(n), random_sym_invertible(n)
        BX = [sum((X[a] * Rm[(a, i)] for a in range(n)), sp.zeros(n, n)) for i in range(n)]
        BY = [sum((Y[a] * Rm[(a, i)] for a in range(n)), sp.zeros(n, n)) for i in range(n)]
        H, Hp = gi * h, gi * hp
        d = F_value(Rm, n, gi, X, Y, h, hp) - F_value(Rm, n, gi, X, Y, hp, h)
        jac = -sum(gi[i, i] * ((BX[i] * BY[i] - BY[i] * BX[i]) * (Hp * H - H * Hp)).trace() for i in range(n))
        sym = F_value(Rm, n, gi, X, Y, hp, h) - F_value(Rm, n, gi, Y, X, h, hp)
        check(f"{lab} (random rational data): F(X,Y;h,h') - F(X,Y;h',h) = -sum_i eta_i tr([B_i,C_i][H',H]) = {d}"
              f" (nonzero, as expected: the bilinear form is not symmetric in h, h'), and F(X,Y;h',h) = F(Y,X;h,h')",
              sp.expand(d - jac) == 0 and d != 0 and sym == 0)


# ================================================================================================ L2
def L2():
    say("\n=== L2. Corollary L2: Einstein tensors (n = 4) and the contradiction with J5(b) ===")
    n = 4
    for lab, sig in SIG4.items():
        R, free = curvature_family(n, sig, einstein=True)
        g = sp.diag(*sig)
        R1 = const_curv(n, sig)
        W = {k_: sp.expand(v - Lam / 3 * R1[k_]) for k_, v in R.items()}
        Wric = ricci(W, n, g.inv())
        weyl_ok = Wric == sp.zeros(n, n) and all(sp.expand(W[a, b, c, d] + W[a, c, d, b] + W[a, d, b, c]) == 0
                                                 for a, b, c, d in itertools.product(range(n), repeat=4))
        C1, C2s, C2, C2p = contractions(R, n, sig)
        W2 = norm2(W, n, sig)
        ok = (len(free) == 10 and weyl_ok and ricci(R, n, g.inv()) == Lam * g and norm2(R1, n, sig) == 24
              and sp.expand(C1 - 6 * W2 - 16 * Lam ** 2) == 0
              and sp.expand(C2s - sp.Rational(3, 2) * W2 - 12 * Lam ** 2) == 0
              and sp.expand(4 * C2s - C1 - 32 * Lam ** 2) == 0
              and sp.expand(C2 - sp.Rational(3, 4) * W2 - 16 * Lam ** 2) == 0)
        check(f"{lab} Einstein, R = W + (Lambda/3) R_1 with W an arbitrary Weyl tensor (10 parameters), |R_1|^2 = 24:"
              f" C1 = 6|W|^2 + 16 Lambda^2, C2sym = 3/2 |W|^2 + 12 Lambda^2, so 4 C2sym - C1 = 32 Lambda^2"
              f" (side: C2 = 3/4 |W|^2 + 16 Lambda^2)", ok)
        if lab in ('(1,3)', '(2,2)'):
            cc1, cc2s, _, _ = contractions(const_curv(n, sig), n, sig)
            say(f"     constant curvature +1 (Lambda = 3, W = 0): C1 = {cc1}, C2sym = {cc2s}  (nonzero, as in J5(d))")
    say("  So F(X,Y;h,h) = 0 for all X, Y, h forces Lambda = 0 for every Einstein tensor, every signature; the scale")
    say("  of R is irrelevant (F is homogeneous of degree 2 in R), and so is the sign of Lambda.")
    cH, cV = sp.symbols('c_H c_V', real=True, nonzero=True)
    Lg = -cH / cV
    check("J5(b): at x0 = pi(p), Ric(g) = Lambda_g g with Lambda_g = -c_H/c_V, so 4 C2sym - C1 = 8|Ric|^2"
          " = 32 Lambda_g^2 = 32 c_H^2/c_V^2 > 0, whereas the order-eps^2 leaf condition gives C1 = C2sym = 0",
          sp.simplify(32 * Lg ** 2 - 32 * cH ** 2 / cV ** 2) == 0 and (32 * Lg ** 2).is_positive)
    say("  Hence no radial Z2 solution of Upsilon = 0 on a neighbourhood of a point of s(X) in the two settings of J1-J2:")
    say("  (14,0) (definite base) and (7,7) (base (1,3) or (3,1), A4; base (2,2), A5); lambda > 1/4 for (7,7).")
    say("  The leaf argument itself does not depend on the signature; (8,6), (10,4) etc. lack J1-J2, hence J5(b).")


# ================================================================================================ L3
def L3():
    say("\n=== L3. Proposition L3: positivity for a Lorentzian (or definite) base, one value of F ===")
    n = 4
    for lab, (sig, T) in {'(4,0)': (SIG4['(4,0)'], 0), '(1,3)': (SIG4['(1,3)'], 0), '(3,1)': (SIG4['(3,1)'], 3),
                          '(2,2)': (SIG4['(2,2)'], 0)}.items():
        g = sp.diag(*sig)
        gi = g.inv()
        R, free = curvature_family(n, sig)
        Rm = rend(R, n, gi)
        e = [[1 if i == a else 0 for i in range(n)] for a in range(n)]
        hT = sp.Matrix(n, n, lambda i, j: g[i, T] * g[j, T])        # T^b (x) T^b
        FT = F_value(Rm, n, gi, e[T], e[T], hT, hT)
        perp = [i for i in range(n) if i != T]
        # R(T,e_i,T,e_j) := g(R(T,e_i)T, e_j)
        RT = {(i, j): (sp.Matrix([e[j]]) * g * Rm[(T, i)] * sp.Matrix(e[T]))[0, 0] for i in perp for j in perp}
        rhs = 2 * g[T, T] * sum(sig[i] * sig[j] * RT[i, j] ** 2 for i in perp for j in perp)
        definite = len(set(sig[i] for i in perp)) == 1
        check(f"{lab}, T = e_{T} (g(T,T) = {g[T, T]}): F(T,T; T^b (x) T^b) = 2 g(T,T) sum_ij eta_i eta_j R(T,e_i,T,e_j)^2"
              f" over T^perp, generic R (20 components); T^perp definite: {definite}",
              sp.expand(FT - rhs) == 0)
    # R(T, X, T, Y) for all T, X, Y determines R: the linear map R -> coefficients has rank 20
    sig = SIG4['(1,3)']
    g = sp.diag(*sig)
    R, free = curvature_family(n, sig)
    Ts, Xs, Ys = sp.symbols('T0:4'), sp.symbols('X0:4'), sp.symbols('Y0:4')
    poly = sp.expand(sum(R[a, b, c, d] * Ys[a] * Ts[b] * Ts[c] * Xs[d] for a, b, c, d in itertools.product(range(n), repeat=4)))
    cs = sp.Poly(poly, *Ts, *Xs, *Ys).coeffs()
    Mx = sp.Matrix([[sp.diff(cf, s) for s in free] for cf in cs])
    check("R(T,X,T,Y) := g(R(T,X)T, Y) as a polynomial in (T, X, Y) determines R: the linear map from the 20"
          " components has rank 20 (so its vanishing for T in the open timelike cone forces R = 0)", Mx.rank() == 20)
    say("  So for a base of signature (1,3), (3,1) or (4,0): F(X,X;h,h) = 0 for all X, h forces R = 0, for every")
    say("  algebraic curvature tensor (Einstein or not).  For (2,2), T^perp is indefinite and no such value is definite (L4).")


# ================================================================================================ L4
def coeff_matrix(F, outer, inner):
    """F polynomial in outer (X, h) variables; its coefficients are polynomials in inner; matrix of those"""
    cs = sp.Poly(F, *outer).coeffs()
    mons = sorted(set().union(*[set(sp.Poly(c, *inner).monoms()) for c in cs]))
    M = sp.Matrix([[sp.Poly(c, *inner).coeff_monomial(m) for m in mons] for c in cs])
    Mr = M.rref()[0]
    Mr = Mr[[i for i in range(Mr.rows) if any(Mr[i, j] != 0 for j in range(Mr.cols))], :]
    return Mr, mons, len(cs)


def in_span(Mr, mons, inner, p):
    pp = sp.Poly(sp.expand(p), *inner)
    if any(m not in mons for m in pp.monoms()):
        return False
    v = sp.Matrix([[pp.coeff_monomial(m) for m in mons]])
    return Mr.rank() == Mr.col_join(v).rank()


def F_poly(R, n, sig):
    g = sp.diag(*sig)
    gi = g.inv()
    Rm = rend(R, n, gi)
    X = sp.symbols('x0:4')
    h, hs = sym_generic(n, 'h')
    return F_value(Rm, n, gi, X, X, h, h), list(X) + list(hs)


PAIRS = [(a, b) for a in range(4) for b in range(a + 1, 4)]


def curvature_operator(R, sig):
    """6 x 6 matrices on 2-forms (basis e^a ^ e^b, a < b): (R w)_ab = 1/2 R_ab^cd w_cd, (*w)_ab = 1/2 eps_ab^cd w_cd,
    and the inner product <w, w'> = 1/2 w_ab w'^ab"""
    eta = sig
    Rop = sp.Matrix(6, 6, lambda P, Q: R[PAIRS[P][0], PAIRS[P][1], PAIRS[Q][0], PAIRS[Q][1]] * eta[PAIRS[Q][0]] * eta[PAIRS[Q][1]])
    star = sp.Matrix(6, 6, lambda P, Q: sp.LeviCivita(*PAIRS[P], *PAIRS[Q]) * eta[PAIRS[Q][0]] * eta[PAIRS[Q][1]])
    ip = sp.diag(*[eta[a] * eta[b] for a, b in PAIRS])
    return Rop, star, ip


def L4():
    say("\n=== L4. Proposition L4: the coefficients of the quadratic form F(X,X;h,h) (n = 4, Einstein); base (2,2) ===")
    n = 4
    keep = {}
    for lab in ('(1,3)', '(2,2)', '(4,0)'):
        sig = SIG4[lab]
        R, free = curvature_family(n, sig, einstein=True)
        F, outer = F_poly(R, n, sig)
        inner = list(free) + [Lam]
        Mr, mons, ncoef = coeff_matrix(F, outer, inner)
        keep[lab] = (R, free, F, outer, inner, Mr, mons)
        check(f"{lab}: the {ncoef} coefficients of F(X,X;h) (X in R^4, h in Sym^2) span a space of quadratic forms in"
              f" (W, Lambda) of dimension {Mr.rows} (of {len(mons)}), and Lambda^2 lies in it"
              f"  [{time.time() - T0:.0f} s]", Mr.rows == 47 and len(mons) == 66 and in_span(Mr, mons, inner, Lam ** 2))
    say("  (3,1) and (0,4) follow from (1,3) and (4,0) by g -> -g, which leaves the endomorphisms R(X,Y) and G_0")
    say("  unchanged and sends F -> -F (and Lambda -> -Lambda): same span, same rank.")
    # base (2,2): self-dual / anti-self-dual blocks
    lab = '(2,2)'
    sig = SIG4[lab]
    R, free, F, outer, inner, Mr, mons = keep[lab]
    Rop, star, ip = curvature_operator(R, sig)
    Vp = sp.Matrix.hstack(*(star - sp.eye(6)).nullspace())
    Vm = sp.Matrix.hstack(*(star + sp.eye(6)).nullspace())
    Rp = (Vp.T * Vp).inv() * Vp.T * Rop * Vp
    Rmn = (Vm.T * Vm).inv() * Vm.T * Rop * Vm
    blocks_ok = (sp.expand(Rop * Vp - Vp * Rp) == sp.zeros(6, 3) and sp.expand(Rop * Vm - Vm * Rmn) == sp.zeros(6, 3)
                 and sp.expand(Rp.trace() - Lam) == 0 and sp.expand(Rmn.trace() - Lam) == 0)
    check("(2,2): *^2 = 1, [R, *] = 0 for the Einstein curvature operator, R = R+ (+) R- on Lambda+ (+) Lambda-, induced"
          f" forms {(Vp.T * ip * Vp).tolist()} and {(Vm.T * ip * Vm).tolist()} (signature (1,2) each), tr R+ = tr R- ="
          " Lambda (so R+- = W+- + Lambda/3)",
          star ** 2 == sp.eye(6) and sp.expand(Rop * star - star * Rop) == sp.zeros(6, 6) and blocks_ok)
    prods = [Rp[i, j] * Rmn[a, b] for i in range(3) for j in range(3) for a in range(3) for b in range(3)]
    Mp, _, _ = coeff_matrix(sum(p * sp.Symbol(f"pp{i}") for i, p in enumerate(prods)), sp.symbols('pp0:81'), inner)
    check(f"(2,2): all 81 products (R+)_ij (R-)_ab lie in the span of the coefficients of F (they span {Mp.rows} of its"
          f" {Mr.rows} dimensions); hence F = 0 forces R+ = 0 or R- = 0",
          all(in_span(Mr, mons, inner, p) for p in prods) and Mp.rows == 36)
    for nm, Rz, Ro in (('R+ = 0', Rp, Rmn), ('R- = 0', Rmn, Rp)):
        sol = sp.solve(list(Rz), inner, dict=True)[0]
        rest = [v for v in inner if v not in sol]
        Fr = sp.expand(F.subs(sol))
        Mz, monz, _ = coeff_matrix(Fr, outer, rest)
        Ro2 = sp.expand(Ro.subs(sol) ** 2)
        ent = [Ro2[i, j] for i in range(3) for j in range(3)]
        Me, mone, _ = coeff_matrix(sum(e * sp.Symbol(f"zz{i}") for i, e in enumerate(ent)), sp.symbols('zz0:9'), rest)
        back = all(in_span(Me, mone, rest, sum(Mz[i, j] * sp.prod([v ** e_ for v, e_ in zip(rest, monz[j])])
                                                 for j in range(Mz.cols))) for i in range(Mz.rows))
        check(f"(2,2), on {nm}: Lambda = 0 there, {len(rest)} parameters remain (the other half W); the span of the"
              f" coefficients of F equals the span of the entries of (the other block)^2 (both of dimension {Mz.rows})",
              sol.get(Lam) == 0 and Mz.rows == Me.rows and all(in_span(Mz, monz, rest, e) for e in ent) and back)
    say("  So in base (2,2) the Einstein solutions of F = 0 are exactly: Lambda = 0, one of W+, W- zero and the other")
    say("  self-adjoint with square 0, i.e. R = c n (x) n with n a null element of Lambda+ or Lambda- (c real).")
    # the null example and the full condition
    g = sp.diag(*sig)
    gi = g.inv()
    nv = Vm[:, 0] + Vm[:, 2]
    nf = {}
    for P, (a, b) in enumerate(PAIRS):
        nf[(a, b)], nf[(b, a)] = nv[P], -nv[P]
    for a in range(n):
        nf[(a, a)] = 0
    Rn = {(a, b, c, d): nf[(a, b)] * nf[(c, d)] for a, b, c, d in itertools.product(range(n), repeat=4)}
    bian = all(Rn[a, b, c, d] + Rn[a, c, d, b] + Rn[a, d, b, c] == 0 for a, b, c, d in itertools.product(range(n), repeat=4))
    Fn, _ = F_poly(Rn, n, sig)
    ys, _ = sym_generic(n, 'y')
    Xs, Ys = sp.symbols('X0:4'), sp.symbols('Y0:4')
    Rm = rend(Rn, n, gi)
    BX = [sum((Xs[a] * Rm[(a, i)] for a in range(n)), sp.zeros(n, n)) for i in range(n)]
    BY = [sum((Ys[a] * Rm[(a, i)] for a in range(n)), sp.zeros(n, n)) for i in range(n)]
    adj = ys.adjugate()
    full = sp.expand(sum(gi[i, i] * (adj * act(BX[i], ys) * adj * act(BY[i], ys)).trace() for i in range(n)))
    nn = sp.Matrix(n, n, lambda b, d: sum(gi[a, c] * nf[(a, b)] * nf[(c, d)] for a in range(n) for c in range(n)))
    check(f"(2,2) example: n = {list(nv)} (components on e^a^e^b, a<b) is anti-self-dual and null; R = n (x) n satisfies"
          " Bianchi, Ric = 0, n o n = 0, F = 0, and the full leaf condition det(y)^2 Phi_y(X,Y) = 0 identically in"
          " the symmetric y and in X, Y",
          star * nv == -nv and (nv.T * ip * nv)[0, 0] == 0 and bian and ricci(Rn, n, gi) == sp.zeros(n, n)
          and nn == sp.zeros(n, n) and Fn == 0 and full == 0)
    say("  (By hand: R(X,e_i) = n(X,e_i) N with N the endomorphism of n, so every term of Phi_y carries the factor")
    say("   sum_ij g^ij n(X,e_i) n(Y,e_j) = (n o n)(X,Y) = 0.)  These Ricci-flat tensors pass the HH leaf condition to")
    say("   all orders; they are excluded in J5 only because there Lambda_g = -c_H/c_V != 0.")


def main():
    L0()
    L1()
    L2()
    L3()
    L4()
    nf = [lab for lab, ok in RES if not ok]
    say(f"\nSUMMARY: {len(RES)} checks, {len(RES) - len(nf)} True, {len(nf)} False")
    say(f"\nruntime {time.time() - T0:.0f} s")


if __name__ == '__main__':
    main()
