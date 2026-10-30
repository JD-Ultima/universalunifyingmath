"""GU Phase 2(a): signature of natural fiber metrics on the metric bundle Y = Met(X^n)
and of the chimeric metric on C = V (+) H*.

Claims tested: GU-012, GU-013, GU-014, GU-015 (claims/geometric-unity.md).
Write-up: formalization/geometric-unity/A-chimeric-signature.md

Everything is exact (sympy Rationals). Inertia of a real symmetric matrix is computed
by Descartes' rule of signs on the characteristic polynomial, which is exact because a
real symmetric matrix has only real eigenvalues (imported lemma L-A3), and independently
by an exact LDL^T (Sylvester) elimination; the two must agree.
"""
import itertools
import random
import sympy as sp

lam = sp.Symbol('lambda', real=True)


# ---------------------------------------------------------------- linear algebra helpers
def inertia_descartes(M):
    """(n_plus, n_minus, n_zero) of a real symmetric matrix via Descartes on charpoly."""
    x = sp.Symbol('x')
    p = sp.Poly(M.charpoly(x).as_expr(), x)
    coeffs = p.all_coeffs()            # highest degree first
    # strip factors of x (zero eigenvalues)
    n0 = 0
    while coeffs and coeffs[-1] == 0:
        coeffs.pop()
        n0 += 1

    def sign_changes(cs):
        cs = [c for c in cs if c != 0]
        return sum(1 for a, b in zip(cs, cs[1:]) if (a > 0) != (b > 0))
    npos = sign_changes(coeffs)
    deg = len(coeffs) - 1
    neg_coeffs = [c * (-1) ** (deg - i) for i, c in enumerate(coeffs)]
    nneg = sign_changes(neg_coeffs)
    assert npos + nneg + n0 == M.shape[0]
    return npos, nneg, n0


def inertia_ldl(M):
    """Exact Sylvester inertia by symmetric elimination with 2x2 handling of zero pivots."""
    A = sp.Matrix(M)
    n = A.shape[0]
    pos = neg = zero = 0
    idx = list(range(n))
    while idx:
        # choose a nonzero diagonal pivot
        piv = next((i for i in idx if A[i, i] != 0), None)
        if piv is not None:
            d = A[piv, piv]
            pos += int(bool(d > 0))
            neg += int(bool(d < 0))
            rest = [j for j in idx if j != piv]
            for i in rest:
                for j in rest:
                    A[i, j] = A[i, j] - A[i, piv] * A[piv, j] / d
            idx = rest
            continue
        # all diagonal entries zero: find an off-diagonal nonzero; 2x2 block [[0,b],[b,0]] -> (+1,-1)
        pair = next(((i, j) for i in idx for j in idx if i < j and A[i, j] != 0), None)
        if pair is None:
            zero += len(idx)
            break
        i, j = pair
        # change basis e_i -> e_i + e_j makes diagonal entry 2*A[i,j] != 0
        for k in range(n):
            A[i, k] = A[i, k] + A[j, k]
        for k in range(n):
            A[k, i] = A[k, i] + A[k, j]
    return pos, neg, zero


# ---------------------------------------------------------------- the DeWitt family
def sym_basis(n):
    """Basis of Sym^2(R^n): E_ii and E_ij + E_ji (i<j)."""
    B = []
    for i in range(n):
        for j in range(i, n):
            E = sp.zeros(n, n)
            E[i, j] = 1
            E[j, i] = 1
            B.append(E)
    return B


def dewitt_gram(g, lam_val):
    """Gram matrix of G_lam(h,k) = tr(g^-1 h g^-1 k) - lam tr(g^-1 h) tr(g^-1 k) on Sym^2."""
    gi = g.inv()
    B = sym_basis(g.shape[0])
    N = len(B)
    G = sp.zeros(N, N)
    for a in range(N):
        for b in range(a, N):
            v = (gi * B[a] * gi * B[b]).trace() - lam_val * (gi * B[a]).trace() * (gi * B[b]).trace()
            G[a, b] = G[b, a] = sp.nsimplify(v) if not isinstance(lam_val, sp.Symbol) else sp.expand(v)
    return G


def eta(p, q):
    return sp.diag(*([1] * p + [-1] * q))


def closed_form_fiber(p, q, regime):
    """Theorem A1: signature of G_lam for base signature (p,q); regime in {'<', '>'} relative to 1/n."""
    P = p * (p + 1) // 2 + q * (q + 1) // 2
    if regime == '<':
        return (P, p * q)
    return (P - 1, p * q + 1)


def main():
    random.seed(20210401)
    print("=== A1. Fiber signature of G_lambda, symbolic check of the spectrum (diagonal g) ===")
    for (p, q) in [(1, 3), (3, 1), (4, 0), (2, 2)]:
        G = dewitt_gram(eta(p, q), lam)
        # spectrum of the symbolic Gram matrix
        ev = G.eigenvals()
        print(f"base ({p},{q}): eigenvalues of Gram matrix in lambda = "
              + ", ".join(f"{sp.factor(k)} (x{v})" for k, v in ev.items()))

    print("\n=== A2. Exact inertia at sample lambda, diagonal and non-diagonal g, both methods ===")
    lam_samples = [sp.Rational(-3), sp.Rational(0), sp.Rational(1, 8), sp.Rational(1, 4),
                   sp.Rational(3, 10), sp.Rational(1), sp.Rational(7)]
    for (p, q) in [(1, 3), (3, 1), (4, 0)]:
        n = p + q
        for lv in lam_samples:
            res = []
            for trial in range(3):
                if trial == 0:
                    g = eta(p, q)
                else:   # random GL(4,Z) congruence: g = A^T eta A (non-diagonal, same signature)
                    while True:
                        A = sp.Matrix(n, n, lambda i, j: random.randint(-3, 3))
                        if A.det() != 0:
                            break
                    g = A.T * eta(p, q) * A
                G = dewitt_gram(g, lv)
                i1 = inertia_descartes(G)
                i2 = inertia_ldl(G)
                assert i1 == i2, (p, q, lv, i1, i2)
                res.append(i1)
            assert len(set(res)) == 1
            pos, neg, z = res[0]
            if lv == sp.Rational(1, n):
                pred = "degenerate"
                ok = (z == 1)
            else:
                pred = closed_form_fiber(p, q, '<' if lv < sp.Rational(1, n) else '>')
                ok = (pos, neg, z) == (pred[0], pred[1], 0)
            assert ok
            print(f"base ({p},{q}) lambda={str(lv):>5}: inertia (+,-,0) = {res[0]}  closed form {pred}  OK")

    print("\n=== A3. GL(n)-invariance G_lam(A^T h A, A^T k A)|_{A^T g A} = G_lam(h,k)|_g (symbolic) ===")
    n = 4
    hs = sp.Matrix(n, n, lambda i, j: sp.Symbol(f'h{min(i,j)}{max(i,j)}'))
    ks = sp.Matrix(n, n, lambda i, j: sp.Symbol(f'k{min(i,j)}{max(i,j)}'))
    g = eta(1, 3)
    A = sp.Matrix([[1, 2, 0, -1], [0, 1, 3, 0], [2, 0, 1, 1], [0, -1, 0, 2]])
    gA, hA, kA = A.T * g * A, A.T * hs * A, A.T * ks * A

    def Gl(gg, h, k):
        gi = gg.inv()
        return (gi * h * gi * k).trace() - lam * (gi * h).trace() * (gi * k).trace()
    diff = sp.simplify(Gl(gA, hA, kA) - Gl(g, hs, ks))
    print("difference simplifies to:", diff)
    assert diff == 0
    print("invariance under g -> -g:", sp.simplify(Gl(-g, hs, ks) - Gl(g, hs, ks)))

    print("\n=== A4. Chimeric signature on C = V (+) H*: g_C = sV*G_lam (+) sH*g^{-1}, n = 4 ===")
    rows = []
    for (p, q) in [(1, 3), (3, 1), (4, 0)]:
        for regime in ['<', '>']:
            fp, fn = closed_form_fiber(p, q, regime)
            for sV, sH in itertools.product([1, -1], repeat=2):
                Vp, Vn = (fp, fn) if sV == 1 else (fn, fp)
                Hp, Hn = (p, q) if sH == 1 else (q, p)
                tot = (Vp + Hp, Vn + Hn)
                rows.append(((p, q), regime, sV, sH, (Vp, Vn), (Hp, Hn), tot))
                print(f"base ({p},{q}) lambda{regime}1/4 sV={sV:+d} sH={sH:+d}: "
                      f"V {Vp, Vn} + H {Hp, Hn} = {tot}{'   <-- (7,7)' if tot == (7, 7) else ''}")
    # exact cross-check of one (7,7) case by building the 14x14 Gram matrix explicitly
    print("\nexact 14x14 cross-check:")
    for (p, q, lv, sV, sH) in [(1, 3, sp.Rational(1), 1, 1), (3, 1, sp.Rational(1), 1, -1),
                               (1, 3, sp.Rational(0), 1, 1), (4, 0, sp.Rational(1), 1, -1)]:
        g = eta(p, q)
        GC = sp.diag(sV * dewitt_gram(g, lv), sH * g.inv())
        print(f"  base ({p},{q}) lambda={lv} sV={sV:+d} sH={sH:+d}: inertia {inertia_ldl(GC)}")
    sevens = [r for r in rows if r[6] == (7, 7)]
    print("\nall (7,7) cases for n=4 with base in {(1,3),(3,1),(4,0)}:")
    for r in sevens:
        print("  ", r)

    print("\n=== A5. Classification: all base signatures (p,q), n<=40, giving neutral chimeric signature ===")
    neutral = []
    for n in range(1, 41):
        for p in range(n + 1):
            q = n - p
            for regime in ['<', '>']:
                fp, fn = closed_form_fiber(p, q, regime)
                for sV, sH in itertools.product([1, -1], repeat=2):
                    Vp, Vn = (fp, fn) if sV == 1 else (fn, fp)
                    Hp, Hn = (p, q) if sH == 1 else (q, p)
                    if Vp + Hp == Vn + Hn:
                        neutral.append((p, q, regime, sV, sH, Vp + Hp))
    for t in neutral:
        print(f"  (p,q)=({t[0]},{t[1]}) lambda{t[2]}1/n sV={t[3]:+d} sH={t[4]:+d} -> ({t[5]},{t[5]})")
    # exact verification of the closed form for every (p,q) with n<=6 at lambda = 0 and lambda = 1
    print("\nexact check of closed form (Theorem A1) for all (p,q), 2<=n<=6, lambda in {0, 1}:")
    for n in range(2, 7):
        for p in range(n + 1):
            q = n - p
            for lv, regime in [(sp.Rational(0), '<'), (sp.Rational(1), '>')]:
                G = dewitt_gram(eta(p, q), lv)
                assert inertia_ldl(G) == closed_form_fiber(p, q, regime) + (0,), (p, q, lv)
        print(f"  n={n}: all {n+1} signatures OK")

    print("\n=== A6. Source formula (3.9): traceless signature (i n - i^2, (n^2+2i^2-2in+n-2)/2) ===")
    for (i, n) in [(1, 4), (0, 4), (2, 4)]:
        src = (i * n - i * i, (n * n + 2 * i * i - 2 * i * n + n - 2) // 2)
        p, q = i, n - i
        G = dewitt_gram(eta(p, q), sp.Rational(0))
        # traceless subspace: orthogonal complement of g in Frobenius form; Frobenius(g,g)=n>0
        pos, neg, _ = inertia_ldl(G)
        ours_traceless = (pos - 1, neg)     # trace direction is Frobenius-positive
        print(f"  i={i}, n={n}: source (first,second)={src}; ours (+,-) traceless = {ours_traceless};"
              f" source lists (minus, plus) -> match: {src == (ours_traceless[1], ours_traceless[0])}")


if __name__ == '__main__':
    main()
