"""Far-end conjugation and the direct-sum trace (gaps CT-G03, CT-G04).

Companion of formalization/chain-theory/alpha/conjugation.md.  Every check is exact (sympy over the Gaussian
rationals) and asserted; nothing is random.  "g" is u(H1) + u(H2), pairs (X1, X2) of anti-Hermitian operators on
the two ends of a link, acting on the link L = H1 (x) H2 by the tensor-product derivation X1 (x) 1 + 1 (x) X2
(the canonical action).  A reading is a Lie subalgebra h of g (or, equivalently, a relabelling of g together with
the Letter's trace condition); see the .md, Sec. 2.

Sections
  0  the quoted source passages sit at the cited lines (only with --upstream; the files are read as text)
  1  convention-free form: off-diagonal blocks of Lambda^2 W, Sym^2 W, End W for W = H1 + H2 (Prop. 1, 2)
  2  invariant linear trace conditions from natural spaces F: the coefficients (alpha_F, beta_F) (Lemma 7, Thm 8)
  3  relabellings eps in {+1,-1}^2 (conjugation at one or both ends): what they fix and what they move (Thm 3)
  4  the matching algebra: stabilizer, matching closure, generation (Thm 4, Thm 5, Prop. 6)
  5  readouts: acting types, Theorem 1's selection (config_enumeration.py machinery), hypercharge (Thm 8, Cor. 9)

Run from this directory:
  python3 -I conjugation.py --upstream /path/to/chain-theory > conjugation.out
(--upstream is a clone of the author's repository at a28db89; it is used only by Sec. 0 and only read as text.
Without it Sec. 0 is skipped and the rest is unchanged.)  Imports config_enumeration.py, which is not modified;
its own printout is captured and compared with config_enumeration.out.
"""
import argparse
import contextlib
import hashlib
import io
import itertools
import os
import sys
import time

T0 = time.perf_counter()
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)          # python3 -I does not put the script directory on sys.path

import sympy as sp
from sympy import I, Rational, eye, zeros

ap = argparse.ArgumentParser()
ap.add_argument("--upstream", help="clone of the author's chain-theory repository at a28db89 (read as text)")
ARGS = ap.parse_args()

NCHECK = 0


def check(cond, msg):
    global NCHECK
    assert cond, msg
    NCHECK += 1


def head(s):
    print()
    print(s)
    print("-" * len(s))


# ----------------------------------------------------------------------------------------------------------
# Sec. 0  source anchors
# ----------------------------------------------------------------------------------------------------------
head("0. Source anchors (alpha.tex, supplement.tex, ChainAlpha.lean at a28db89)")
ANCHORS = [
    ("alpha.tex", 43, r"a relative phase that neither end owns"),
    ("alpha.tex", 60, r"with the complex-conjugate basis at the far end"),
    ("alpha.tex", 60, r"acts on the pair as $U\otimes\bar{U}$ and leaves it unchanged"),
    ("alpha.tex", 60, r"nothing lines up the bases of a shared qudit's copies"),
    ("alpha.tex", 81, r"the block-diagonal traceless anti-Hermitian operators on $\mathbb{C}^{d_{1}}\oplus\mathbb{C}^{d_{2}}$, the infinitesimal basis changes at its two ends"),
    ("alpha.tex", 81, r"acts on the link as $\xi_{1}\otimes\openone+\openone\otimes\xi_{2}$"),
    ("alpha.tex", 86, r"as the scalar $i(d_{2}-d_{1})$. At $(2,2)$ it vanishes"),
    ("alpha.tex", 88, r"that every matching preserves and that sees no global phase is the count of $\mathcal{Q}$"),
    ("alpha.tex", 88, r"The count sums to zero over all directions of the two ends, which is the tracelessness of $\mathcal{Q}$"),
    ("supplement.tex", 97, r"By the transpose identity $U$ at the near end acts as $\bar{U}$ at the far end"),
    ("supplement.tex", 170, r"under which the far end transforms by $\bar{U}$ where its partner transforms by $U$"),
    ("supplement.tex", 175, r"An element $\xi=\xi_{1}\oplus\xi_{2}$ acts on the pair over the matched directions as $\xi_{1}\otimes\openone+\openone\otimes\xi_{2}$"),
    ("supplement.tex", 182, r"tracelessness gives $d_{1}x+d_{2}y=0$"),
    ("supplement.tex", 187, r"the elements $\zeta\oplus\bar{\zeta}$ with $\zeta$ on the qubit and $\bar{\zeta}$ on the two-dimensional common support of the qutrit"),
    ("supplement.tex", 187, r"at the unbroken generator $\mathcal{U}=i\,\mathrm{diag}(1,1;\,-1,-1,0)$, with the third qutrit level fixed at zero"),
    ("supplement.tex", 200, r"sees no global phase of a pair is, at $(2,3)$, a multiple of the count of $\mathcal{Q}$; at $(2,2)$ it is zero"),
    ("supplement.tex", 203, r"A read sees no global phase of a pair, which gives $2x+3y=0$ at $(2,3)$"),
    ("supplement.tex", 203, r"On the full space of the two ends the count sums over all $d_{1}+d_{2}$ directions and is zero"),
    ("supplement.tex", 233, r"the fourth, $2$R $2$R $2$U, satisfies (a) and (b) but carries no acting charge"),
    ("supplement.tex", 281, r"with $\xi$ acting on the link as $\xi_{1}\otimes\openone+\openone\otimes\xi_{2}$ and not on the direct sum $\mathbb{C}^{d_{1}}\oplus\mathbb{C}^{d_{2}}$"),
    ("supplement.tex", 284, r"that is the space the identification algebra acts on, not the space the readout lives on"),
    ("lean/ChainAlpha.lean", 489, r"the far end transforms by the complex conjugate of a matching"),
    ("lean/ChainAlpha.lean", 515, r"def actingCharge (d1 d2 : Nat) : Prop := d1 ≠ d2"),
    ("lean/ChainAlpha.lean", 541, r"as X1 (x) 1 + 1 (x) X2"),
    ("lean/ChainAlpha.lean", 551, r"def Q : Mat GI := tabulate 5 5"),
    ("lean/ChainAlpha.lean", 561, r"def onLink (x : Mat GI) : Mat GI := addG (kronG (block2 x) (idG 3)) (kronG (idG 2) (block3 x))"),
    ("lean/ChainAlpha.lean", 594, r"def U : Mat GI := tabulate 5 5"),
    ("lean/ChainAlpha.lean", 619, r"theorem acting_charge (d1 d2 : Nat) : ((d2 : Int) - d1 ≠ 0) ↔ d1 ≠ d2"),
    ("lean/ChainAlpha.lean", 897, r"(hq : ∀ a ∈ c.arms, actingCharge c.center a.far)"),
]
if ARGS.upstream:
    base = os.path.join(ARGS.upstream, "papers", "alpha")
    lean_path = os.path.join(base, "lean", "ChainAlpha.lean")
    with open(lean_path, "rb") as fh:
        sha = hashlib.sha256(fh.read()).hexdigest()
    check(sha == "522844e911cf5f6091b603db0d1d042f50544d2ec8d10ed9336b722945c331da", "ChainAlpha.lean pin")
    print("ChainAlpha.lean sha256 matches the pin 522844e9...c331da")
    cache = {}
    for fn, ln, snip in ANCHORS:
        if fn not in cache:
            with open(os.path.join(base, fn), encoding="utf-8") as fh:
                cache[fn] = fh.read().split("\n")
        line = cache[fn][ln - 1]
        check(snip in line, f"anchor {fn}:{ln}")
        short = snip if len(snip) <= 70 else snip[:67] + "..."
        print(f"  {fn}:{ln}: found  \"{short}\"")
    # Table S-inputs (supplement.tex:23-41) does not list the identification algebra, a trace or a conjugation
    tab = "\n".join(cache["supplement.tex"][22:41])
    for w in ["identification algebra", "direct sum", "oplus", "conjugat", "traceless", "acting"]:
        check(w not in tab, f"Table S-inputs mentions {w}")
    print("  supplement.tex:23-41 (Table S-inputs): none of 'identification algebra', 'direct sum', 'oplus', 'conjugat',")
    print("  'traceless', 'acting' occurs (its only 'identification' is the Photon Conjecture's, with alpha(0))")
else:
    print("skipped (no --upstream given)")

# ----------------------------------------------------------------------------------------------------------
# helpers
# ----------------------------------------------------------------------------------------------------------


def Emat(n, a, b):
    m = zeros(n, n)
    m[a, b] = 1
    return m


def u_basis(d):
    B = []
    for j in range(d):
        for k in range(j + 1, d):
            B.append(Emat(d, j, k) - Emat(d, k, j))
            B.append(I * (Emat(d, j, k) + Emat(d, k, j)))
    for j in range(d):
        B.append(I * Emat(d, j, j))
    return B


def su_basis(d):
    B = [b for b in u_basis(d) if not (b.is_diagonal())]
    for j in range(d - 1):
        B.append(I * (Emat(d, j, j) - Emat(d, j + 1, j + 1)))
    return B


def g_basis(d1, d2):
    """basis of g = u(d1) + u(d2), as pairs."""
    return [(b, zeros(d2, d2)) for b in u_basis(d1)] + [(zeros(d1, d1), b) for b in u_basis(d2)]


def kron(A, B):
    return sp.kronecker_product(A, B)


def conjm(A):
    return A.applyfunc(sp.conjugate)


def bar(A, e):
    """A^(e): A for e = +1, its complex conjugate for e = -1."""
    return A if e == 1 else conjm(A)


def on_link(X, eps=(1, 1)):
    X1, X2 = X
    return kron(bar(X1, eps[0]), eye(X2.rows)) + kron(eye(X1.rows), bar(X2, eps[1]))


def tr(A):
    return sp.nsimplify(sp.expand(A.trace()))


def realvec(mats):
    v = []
    for M in mats:
        for z in M:
            z = sp.expand(z)
            v += [sp.re(z), sp.im(z)]
    return v


def rank_of(vectors):
    if not vectors:
        return 0
    return sp.Matrix(vectors).rank()


def perm_sign(seq):
    s = 1
    seq = list(seq)
    for i in range(len(seq)):
        for j in range(i + 1, len(seq)):
            if seq[i] > seq[j]:
                s = -s
    return s


def ext_power(A, k):
    n = A.rows
    basis = list(itertools.combinations(range(n), k))
    idx = {S: i for i, S in enumerate(basis)}
    M = zeros(len(basis), len(basis))
    for j, S in enumerate(basis):
        for p, s in enumerate(S):
            for c in range(n):
                a = A[c, s]
                if a == 0:
                    continue
                T = list(S)
                T[p] = c
                if len(set(T)) < k:
                    continue
                M[idx[tuple(sorted(T))], j] += perm_sign(T) * a
    return M, basis


def sym2(A):
    n = A.rows
    basis = list(itertools.combinations_with_replacement(range(n), 2))
    idx = {S: i for i, S in enumerate(basis)}
    M = zeros(len(basis), len(basis))
    for j, (a, b) in enumerate(basis):
        for c in range(n):
            if A[c, a] != 0:
                M[idx[tuple(sorted((c, b)))], j] += A[c, a]
            if A[c, b] != 0:
                M[idx[tuple(sorted((a, c)))], j] += A[c, b]
    return M, basis


def ad_on_End(Xw):
    """adjoint action Y -> [Xw, Y] on gl(n), basis E_ab (row-major)."""
    n = Xw.rows
    basis = [(a, b) for a in range(n) for b in range(n)]
    idx = {S: i for i, S in enumerate(basis)}
    M = zeros(n * n, n * n)
    for j, (a, b) in enumerate(basis):
        Y = Emat(n, a, b)
        Z = Xw * Y - Y * Xw
        for p in range(n):
            for q in range(n):
                if Z[p, q] != 0:
                    M[idx[(p, q)], j] += Z[p, q]
    return M, basis


def omega_vec(d1, d2):
    m = min(d1, d2)
    w = zeros(d1 * d2, 1)
    for i in range(m):
        w[i * d2 + i] = 1
    return w


def Qgen(d1, d2):
    return (I * d2 * eye(d1), -I * d1 * eye(d2))


TYPES = [(2, 2), (2, 3), (3, 3)]
EXTRA = [(2, 4), (3, 4)]

# ----------------------------------------------------------------------------------------------------------
# Sec. 1  convention-free form
# ----------------------------------------------------------------------------------------------------------
head("1. Off-diagonal blocks of Lambda^2 W, Sym^2 W and End W, W = H1 + H2 (Prop. 1, Prop. 2)")
print("For every basis element X = (X1, X2) of g (checked on a basis; all maps are linear):")
for (d1, d2) in TYPES + EXTRA:
    n = d1 + d2
    nb = 0
    for X in g_basis(d1, d2):
        X1, X2 = X
        Xw = sp.diag(X1, X2)
        A_link = kron(X1, eye(d2)) + kron(eye(d1), X2)            # reading A link action
        As_link = kron(X1, eye(d2)) + kron(eye(d1), conjm(X2))    # reading A* link action
        As2_link = kron(X2, eye(d1)) + kron(eye(d2), conjm(X1))   # the other off-diagonal block of End W
        # Lambda^2 W: block span{e_a ^ f_b} ~ e_a (x) f_b
        M, basis = ext_power(Xw, 2)
        rows = [basis.index((a, d1 + b)) for a in range(d1) for b in range(d2)]
        others = [i for i in range(len(basis)) if i not in rows]
        blk = M.extract(rows, rows)
        check(blk == A_link, "Lambda2 block")
        check(M.extract(others, rows).is_zero_matrix, "Lambda2 block invariant")
        # Sym^2 W
        M, basis = sym2(Xw)
        rows = [basis.index((a, d1 + b)) for a in range(d1) for b in range(d2)]
        others = [i for i in range(len(basis)) if i not in rows]
        check(M.extract(rows, rows) == A_link, "Sym2 block")
        check(M.extract(others, rows).is_zero_matrix, "Sym2 block invariant")
        # End W, adjoint action: Hom(H2,H1) = span{E_{a, d1+b}}, Hom(H1,H2) = span{E_{d1+b, a}}
        M, basis = ad_on_End(Xw)
        rows = [basis.index((a, d1 + b)) for a in range(d1) for b in range(d2)]
        others = [i for i in range(len(basis)) if i not in rows]
        check(M.extract(rows, rows) == As_link, "End block Hom(H2,H1)")
        check(M.extract(others, rows).is_zero_matrix, "End block invariant")
        rows2 = [basis.index((d1 + b, a)) for b in range(d2) for a in range(d1)]
        check(M.extract(rows2, rows2) == As2_link, "End block Hom(H1,H2)")
        nb += 1
    print(f"  ({d1},{d2}): {nb} basis elements: Lambda^2 W and Sym^2 W restrict on H1(x)H2 to X1(x)1+1(x)X2 (reading A);"
          f" ad on Hom(H2,H1) is X1(x)1+1(x)conj(X2) (reading A*), on Hom(H1,H2) X2(x)1+1(x)conj(X1); blocks invariant")

print()
print("Central charge of Q = i diag(d2 1_d1; -d1 1_d2) on each summand (value c means Q acts as i c):")
print("  type   | L2(H1) | H1(x)H2 | L2(H2) | S2(H1) | S2(H2) | Hom(H2,H1) | Hom(H1,H2) | End H1, End H2")
for (d1, d2) in TYPES + EXTRA:
    Q1, Q2 = Qgen(d1, d2)
    Qw = sp.diag(Q1, Q2)
    vals = {}
    M, basis = ext_power(Qw, 2)
    blocks = {"L2(H1)": [S for S in basis if S[1] < d1], "H1(x)H2 [L2]": [S for S in basis if S[0] < d1 <= S[1]],
              "L2(H2)": [S for S in basis if S[0] >= d1]}
    for name, Ss in blocks.items():
        ev = {sp.simplify(M[basis.index(S), basis.index(S)] / I) for S in Ss}
        check(len(ev) <= 1, "scalar on block")
        vals[name] = ev.pop() if ev else None
    M, basis = sym2(Qw)
    blocks = {"S2(H1)": [S for S in basis if S[1] < d1], "H1(x)H2 [S2]": [S for S in basis if S[0] < d1 <= S[1]],
              "S2(H2)": [S for S in basis if S[0] >= d1]}
    for name, Ss in blocks.items():
        ev = {sp.simplify(M[basis.index(S), basis.index(S)] / I) for S in Ss}
        check(len(ev) == 1, "scalar on block")
        vals[name] = ev.pop()
    M, basis = ad_on_End(Qw)
    blocks = {"Hom(H2,H1)": [(a, b) for (a, b) in basis if a < d1 <= b],
              "Hom(H1,H2)": [(a, b) for (a, b) in basis if b < d1 <= a],
              "End": [(a, b) for (a, b) in basis if (a < d1) == (b < d1)]}
    for name, Ss in blocks.items():
        ev = {sp.simplify(M[basis.index(S), basis.index(S)] / I) for S in Ss}
        check(len(ev) == 1, "scalar on block")
        vals[name] = ev.pop()
    # closed forms
    check(vals["H1(x)H2 [L2]"] == d2 - d1 and vals["H1(x)H2 [S2]"] == d2 - d1, "link charge A")
    check(vals["Hom(H2,H1)"] == d1 + d2 and vals["Hom(H1,H2)"] == -(d1 + d2) and vals["End"] == 0, "End charges")
    check(vals["S2(H1)"] == 2 * d2 and vals["S2(H2)"] == -2 * d1, "Sym2 diagonal charges")
    check(vals["L2(H2)"] == -2 * d1 and (d1 == 1 or vals["L2(H1)"] == 2 * d2), "Lambda2 diagonal charges")
    print(f"  ({d1},{d2})  | {str(vals['L2(H1)']):>6} | {str(vals['H1(x)H2 [L2]']):>7} | {str(vals['L2(H2)']):>6} |"
          f" {str(vals['S2(H1)']):>6} | {str(vals['S2(H2)']):>6} | {str(vals['Hom(H2,H1)']):>10} |"
          f" {str(vals['Hom(H1,H2)']):>10} | {vals['End']}")
print("  closed forms (asserted): H1(x)H2: d2-d1; Hom(H2,H1): d1+d2; Hom(H1,H2): -(d1+d2); Lambda^2/Sym^2 of H1: 2 d2,"
      " of H2: -2 d1; End H1, End H2: 0")
print("  acting on the link: H1(x)H2 (reading A) iff d1 != d2; Hom(H2,H1) (reading A*) for every type")

print()
print("Hypercharge Y = Q/(6i) at (d1,d2) = (2,3), the Letter's Y = diag(1/2,1/2;-1/3,-1/3,-1/3):")
d1, d2 = 2, 3
Yv = {"L2(H1)": Rational(2 * d2, 6), "H1(x)H2": Rational(d2 - d1, 6), "L2(H2)": Rational(-2 * d1, 6),
      "S2(H1)": Rational(2 * d2, 6), "S2(H2)": Rational(-2 * d1, 6),
      "Hom(H2,H1)": Rational(d1 + d2, 6), "Hom(H1,H2)": Rational(-(d1 + d2), 6)}
check(Yv["H1(x)H2"] == Rational(1, 6) and Yv["Hom(H2,H1)"] == Rational(5, 6), "Y values")
print(f"  Lambda^2 W = 10 of SU(5): (1,1)_{Yv['L2(H1)']} + (2,3)_{Yv['H1(x)H2']} + (1,3bar)_{Yv['L2(H2)']}"
      "   [SU(2) x SU(3) labels; Lambda^2 C^3 = 3bar]")
print(f"  Sym^2 W   = 15 of SU(5): (3,1)_{Yv['S2(H1)']} + (2,3)_{Yv['H1(x)H2']} + (1,6)_{Yv['S2(H2)']}")
print(f"  End W = 24 + 1: off-diagonal (2,3bar)_{Yv['Hom(H2,H1)']} + (2,3)_{Yv['Hom(H1,H2)']}   [Hom(H2,H1) = H1 (x) H2*]")
check(sp.binomial(5, 2) == 10 and 15 == 5 * 6 // 2, "dims")

print()
print("Reading B (trace on L itself): h_B = {X : tr_L(X1(x)1+1(x)X2) = 0} = (d rho)^{-1}(su(L)).")
for (d1, d2) in TYPES + EXTRA:
    for eps in [(1, 1), (1, -1), (-1, 1), (-1, -1)]:
        # centre (i x 1, i y 1) acts on L as i(e1 x + e2 y); tr_L-condition is d1 d2 (e1 x + e2 y) = 0
        x, y = sp.symbols("x y", real=True)
        Xc = (I * x * eye(d1), I * y * eye(d2))
        act = sp.expand(on_link(Xc, eps)[0, 0] / I)
        trL = sp.expand(on_link(Xc, eps).trace() / I)
        check(sp.simplify(trL - d1 * d2 * act) == 0, "B trace = d1 d2 times the link scalar")
print("  for every type and every eps in {+1,-1}^2: tr_L of the centre (ix,iy) is d1 d2 times its link scalar,")
print("  so the centre of h_B acts on L as 0: no acting charge for any type under any conjugation (asserted).")

# ----------------------------------------------------------------------------------------------------------
# Sec. 2  invariant trace conditions from natural spaces
# ----------------------------------------------------------------------------------------------------------
head("2. Invariant linear trace conditions tr_F(d sigma_F(X)) = alpha_F tr X1 + beta_F tr X2 (Lemma 7, Thm 8)")
for (d1, d2) in TYPES + EXTRA:
    B = g_basis(d1, d2)
    br = []
    for (a, b) in itertools.combinations(range(len(B)), 2):
        X, Y = B[a], B[b]
        br.append(realvec([X[0] * Y[0] - Y[0] * X[0], X[1] * Y[1] - Y[1] * X[1]]))
    r = rank_of(br)
    check(r == d1 * d1 + d2 * d2 - 2, "[g,g] = su + su")
print("  dim_R [g,g] = d1^2 + d2^2 - 2 at", TYPES + EXTRA, "(so [g,g] = su(d1)+su(d2); an Ad-invariant")
print("  linear functional vanishes on it and is a real combination of the two traces)")

H1, H2 = ("H1",), ("H2",)


def Bar(F):
    return ("bar", F)


def Sum(F, G):
    return ("sum", F, G)


def Ten(F, G):
    return ("tens", F, G)


def Ext(k, F):
    return ("ext", k, F)


def Sym2(F):
    return ("sym2", F)


def End(F):
    return Ten(F, Bar(F))


def rep(F, X):
    op = F[0]
    if op == "H1":
        return X[0]
    if op == "H2":
        return X[1]
    if op == "bar":
        return conjm(rep(F[1], X))
    if op == "sum":
        return sp.diag(rep(F[1], X), rep(F[2], X))
    if op == "tens":
        A, B = rep(F[1], X), rep(F[2], X)
        return kron(A, eye(B.rows)) + kron(eye(A.rows), B)
    if op == "ext":
        return ext_power(rep(F[2], X), F[1])[0]
    if op == "sym2":
        return sym2(rep(F[1], X))[0]
    raise ValueError(op)


D1, D2 = sp.symbols("d1 d2", positive=True, integer=True)


def coeffs(F):
    """(dim, alpha, beta) as polynomials in d1, d2, by the trace calculus of the constructions."""
    op = F[0]
    if op == "H1":
        return (D1, sp.Integer(1), sp.Integer(0))
    if op == "H2":
        return (D2, sp.Integer(0), sp.Integer(1))
    if op == "bar":
        n, a, b = coeffs(F[1])
        return (n, -a, -b)
    if op == "sum":
        n, a, b = coeffs(F[1])
        m, c, e = coeffs(F[2])
        return (n + m, a + c, b + e)
    if op == "tens":
        n, a, b = coeffs(F[1])
        m, c, e = coeffs(F[2])
        return (n * m, sp.expand(m * a + n * c), sp.expand(m * b + n * e))
    if op == "ext":
        k = F[1]
        n, a, b = coeffs(F[2])
        return (sp.expand_func(sp.binomial(n, k)), sp.expand(sp.expand_func(sp.binomial(n - 1, k - 1)) * a),
                sp.expand(sp.expand_func(sp.binomial(n - 1, k - 1)) * b))
    if op == "sym2":
        n, a, b = coeffs(F[1])
        return (sp.expand(n * (n + 1) / 2), sp.expand((n + 1) * a), sp.expand((n + 1) * b))
    raise ValueError(op)


W = Sum(H1, H2)
Wp = Sum(H1, Bar(H2))
Lk = Ten(H1, H2)
Lp = Ten(H1, Bar(H2))
SPACES = [
    ("W = H1+H2", W),
    ("Lambda^2 W", Ext(2, W)),
    ("Sym^2 W", Sym2(W)),
    ("W (x) W", Ten(W, W)),
    ("Lambda^3 W", Ext(3, W)),
    ("W' = H1+conj(H2)", Wp),
    ("Lambda^2 W'", Ext(2, Wp)),
    ("End W = W (x) conj W", End(W)),
    ("End H1 + End H2", Sum(End(H1), End(H2))),
    ("L = H1 (x) H2", Lk),
    ("Lambda^2 L", Ext(2, Lk)),
    ("Sym^2 L", Sym2(Lk)),
    ("End L", End(Lk)),
    ("L' = H1 (x) conj(H2)", Lp),
    ("H1 (centre end)", H1),
    ("H2 (far end)", H2),
    ("H1 + H1 + H2", Sum(Sum(H1, H1), H2)),
]
print("  F                      | dim F                | alpha_F (coefficient of tr X1) | beta_F (coefficient of tr X2)")
COEF = {}
for name, F in SPACES:
    n, a, b = coeffs(F)
    COEF[name] = (n, a, b)
    print(f"  {name:22s} | {str(sp.factor(n)):20s} | {str(sp.factor(a)):30s} | {sp.factor(b)}")
# explicit matrices: the trace calculus agrees with the representation matrices on every basis element
print("  explicit check (representation matrices on every basis element of g; skipped where dim F > 40):")
for (d1, d2) in [(2, 2), (2, 3), (3, 3)]:
    done = []
    for name, F in SPACES:
        n, a, b = COEF[name]
        nn = int(n.subs({D1: d1, D2: d2}))
        if nn > 40:
            continue
        av, bv = a.subs({D1: d1, D2: d2}), b.subs({D1: d1, D2: d2})
        for X in g_basis(d1, d2):
            M = rep(F, X)
            check(M.rows == nn, "dimension")
            check(sp.simplify(M.trace() - (av * X[0].trace() + bv * X[1].trace())) == 0, f"trace form {name}")
        done.append(name)
    print(f"    ({d1},{d2}): {len(done)} spaces agree: {', '.join(done)}")

# ----------------------------------------------------------------------------------------------------------
# Sec. 3  relabellings
# ----------------------------------------------------------------------------------------------------------
head("3. Relabellings theta_eps = (conj^e1, conj^e2), eps in {+1,-1}^2 (Thm 3)")
EPS = [(1, 1), (-1, -1), (1, -1), (-1, 1)]
print("(a) the image of g, and of su+su, on L does not depend on eps:")
for (d1, d2) in TYPES:
    imgs = {}
    for eps in EPS:
        imgs[eps] = [realvec([on_link(X, eps)]) for X in g_basis(d1, d2)]
    r0 = rank_of(imgs[(1, 1)])
    for eps in EPS:
        check(rank_of(imgs[eps]) == r0, "image dim")
        check(rank_of(imgs[(1, 1)] + imgs[eps]) == r0, "same image")
    sub = {eps: [realvec([on_link((a, zeros(d2, d2)), eps)]) for a in su_basis(d1)] +
           [realvec([on_link((zeros(d1, d1), b), eps)]) for b in su_basis(d2)] for eps in EPS}
    rs = rank_of(sub[(1, 1)])
    for eps in EPS:
        check(rank_of(sub[(1, 1)] + sub[eps]) == rs, "same su image")
    check(r0 == d1 * d1 + d2 * d2 - 1 and rs == d1 * d1 + d2 * d2 - 2, "image dims")
    print(f"  ({d1},{d2}): d rho_eps(g) = l (dim {r0}, contains i 1_L) and d rho_eps(su+su) = l_0 (dim {rs}) for all four eps")

print("(b) theta_eps maps the trace functional tr1 + tr2 to e1 tr1 + e2 tr2 (tr conj(X) = -tr X on u(d)):")
for d in [2, 3]:
    for b in u_basis(d):
        check(sp.simplify(conjm(b).trace() + b.trace()) == 0, "tr conj = -tr")
print("  checked on bases of u(2), u(3); so theta_eps^{-1}(h_W) = h_W iff e1 = e2")

print("(c) the Letter's h_W = s(u(d1)+u(d2)) under the link action rho_eps; its centre Q = i(d2 1, -d1 1):")
print("   eps     | type  | Q on L (/i) | i 1_L in d rho_eps(h_W)? | acting")
ACT_EPS = {}
for eps in EPS:
    for (d1, d2) in TYPES:
        Q = Qgen(d1, d2)
        q = sp.simplify(on_link(Q, eps)[0, 0] / I)
        check(on_link(Q, eps) == I * q * eye(d1 * d2), "Q scalar on L")
        check(q == eps[0] * d2 - eps[1] * d1, "formula e1 d2 - e2 d1")
        hb = [(a, zeros(d2, d2)) for a in su_basis(d1)] + [(zeros(d1, d1), b) for b in su_basis(d2)] + [Q]
        img = [realvec([on_link(X, eps)]) for X in hb]
        one = realvec([I * eye(d1 * d2)])
        contains = rank_of(img + [one]) == rank_of(img)
        check(contains == (q != 0), "phase in image iff acting")
        ACT_EPS[(eps, (d1, d2))] = q != 0
        print(f"   {str(eps):8s}| ({d1},{d2}) | {str(q):>11} | {str(contains):>24} | {q != 0}")
for (d1, d2) in TYPES:
    check(ACT_EPS[((1, 1), (d1, d2))] == ACT_EPS[((-1, -1), (d1, d2))], "diagonal relabelling invariant")
    check(ACT_EPS[((1, -1), (d1, d2))] == ACT_EPS[((-1, 1), (d1, d2))] == True, "relative sign -1: always")
print("  the acting predicate depends on eps only through e1 e2: e1 e2 = +1 gives reading A (acting iff d1 != d2),")
print("  e1 e2 = -1 gives reading A* (acting at every type).  (asserted)")

print("(d) the automorphism U -> conj(U):")
a, b = sp.symbols("a b")
U = sp.Matrix([[a, -sp.conjugate(b)], [b, sp.conjugate(a)]])         # SU(2) form (det = |a|^2+|b|^2)
sy = sp.Matrix([[0, -I], [I, 0]])
check(sp.simplify(sy * conjm(U) * sy - U) == zeros(2, 2), "sigma_y conj(U) sigma_y = U on SU(2)")
ph = sp.symbols("phi", real=True)
Uph = sp.exp(I * ph) * U
check(sp.simplify(sy * conjm(Uph) * sy - sp.exp(-2 * I * ph) * Uph) == zeros(2, 2), "U(2): central phase flipped")
print("  SU(2): sigma_y conj(U) sigma_y = U for U = [[a,-conj b],[b,conj a]] (inner);")
print("  U(2):  sigma_y conj(e^{i phi} U) sigma_y = e^{-2 i phi} (e^{i phi} U): the central phase is inverted (not inner)")
U3 = sp.diag(I, I, -1)
check(U3.det() == 1 and sp.expand(U3.trace()) == -1 + 2 * I and sp.expand(conjm(U3).trace()) == -1 - 2 * I, "SU(3) witness")
print("  SU(3): U = diag(i, i, -1) has det 1, tr U = -1+2i, tr conj(U) = -1-2i; inner automorphisms keep the trace,")
print("         so conj is outer on SU(3).  By (a) it does not change the image of su+su on L; it relabels 3 <-> 3bar.")

print("(e) transposition at both ends (the covariant exchange E' = E o T of CT-G11 transposes every slot):")
for (d1, d2) in TYPES:
    for X in g_basis(d1, d2):
        check(on_link(X).T == -on_link(X, (-1, -1)), "transpose = -(eps=(-1,-1))")
print("  (X1(x)1+1(x)X2)^T = -(conj X1 (x) 1 + 1 (x) conj X2) on every basis element at (2,2),(2,3),(3,3):")
print("  it is the diagonal relabelling eps = (-1,-1), e1 e2 = +1, which keeps the selection and sends Y to -Y.")

# ----------------------------------------------------------------------------------------------------------
# Sec. 4  matching algebra, matching closure, generation
# ----------------------------------------------------------------------------------------------------------
head("4. The matching algebra: stabilizer, matching closure, generation (Thm 4, Thm 5, Prop. 6)")


def stabilizer(d1, d2, eps=(1, 1)):
    """real basis of {X in g : rho_eps(X) omega = 0}."""
    B = g_basis(d1, d2)
    w = omega_vec(d1, d2)
    cols = [realvec([on_link(X, eps) * w]) for X in B]
    A = sp.Matrix(cols).T
    ns = A.nullspace()
    out = []
    for v in ns:
        X1 = sum((v[k] * B[k][0] for k in range(len(B))), zeros(d1, d1))
        X2 = sum((v[k] * B[k][1] for k in range(len(B))), zeros(d2, d2))
        out.append((sp.simplify(X1), sp.simplify(X2)))
    return out


def m0_basis(d1, d2):
    """matchings with identity off the support, as operator pairs: (zeta, conj(zeta) + 0), d1 <= d2."""
    assert d1 <= d2
    out = []
    for z in u_basis(d1):
        X2 = zeros(d2, d2)
        X2[:d1, :d1] = conjm(z)
        out.append((z, X2))
    return out


def lie_closure(gens, d1, d2):
    basis, vecs = [], []
    queue = list(gens)
    while queue:
        X = queue.pop(0)
        v = realvec([X[0], X[1]])
        if rank_of(vecs + [v]) > len(vecs):
            vecs.append(v)
            basis.append(X)
            for Y in list(basis):
                queue.append((X[0] * Y[0] - Y[0] * X[0], X[1] * Y[1] - Y[1] * X[1]))
    return basis


print("(a) stabilizer of the pair ket in g (canonical action), and the matching algebra m_0:")
for (d1, d2) in TYPES + EXTRA:
    st = stabilizer(d1, d2)
    check(len(st) == d1 * d1 + (d2 - d1) ** 2, "stab dim")
    m0 = m0_basis(d1, d2)
    w = omega_vec(d1, d2)
    for X in m0:
        check((on_link(X) * w).is_zero_matrix, "matchings fix omega")
    check(rank_of([realvec([X[0], X[1]]) for X in st] + [realvec([X[0], X[1]]) for X in m0]) == len(st), "m0 in stab")
    print(f"  ({d1},{d2}): dim stab = {len(st)} = d1^2 + (d2-d1)^2  (pairs (zeta, conj(zeta) + c), c in u(d2-d1));"
          f" m_0 = {{(zeta, conj(zeta) + 0)}} has dim {len(m0)} and fixes omega")

print("(b) generation: Lie algebra generated by su(H1) + su(H2) and m_0, and by su + su and the full stabilizer:")
for (d1, d2) in TYPES + EXTRA:
    su_pairs = [(a, zeros(d2, d2)) for a in su_basis(d1)] + [(zeros(d1, d1), b) for b in su_basis(d2)]
    gen = lie_closure(su_pairs + m0_basis(d1, d2), d1, d2)
    check(len(gen) == d1 * d1 + d2 * d2 - 1, "generated dim")
    for X in gen:
        check(sp.simplify(X[0].trace() + X[1].trace()) == 0, "generated in h_W")
    genfull = lie_closure(su_pairs + stabilizer(d1, d2), d1, d2)
    expect = d1 * d1 + d2 * d2 - (1 if d1 == d2 else 0)
    check(len(genfull) == expect, "generated with full stabilizer")
    print(f"  ({d1},{d2}): <su+su, m_0> has dim {len(gen)} = d1^2+d2^2-1 and lies in ker(tr1+tr2) = s(u(d1)+u(d2));"
          f" <su+su, stab> has dim {len(genfull)}" + (" (= g)" if d1 != d2 else " (= s(u+u))"))

print("(c) matching closure over the invariant hyperplanes h_[alpha:beta] = {alpha tr X1 + beta tr X2 = 0}:")
al, be = sp.symbols("alpha beta", real=True)
for (d1, d2) in TYPES + EXTRA:
    vals = [sp.expand(al * X[0].trace() + be * X[1].trace()) for X in m0_basis(d1, d2)]
    check(all(sp.expand(v.subs(be, al)) == 0 for v in vals), "alpha = beta contains m_0")
    check(any(sp.expand(v.subs({al: 1, be: 0})) != 0 for v in vals), "alpha != beta misses m_0")
    print(f"  ({d1},{d2}): m_0 in h_[alpha:beta]  <=>  alpha = beta   (on the phase generator: {sp.factor(sp.expand(al*I*d1 - be*I*d1))} = 0)")
print("  in eps-labels the matchings are theta_eps^{-1}(m_0) = {(zeta^(e1), (conj zeta + 0)^(e2))}:")
for (d1, d2) in TYPES:
    for eps in EPS:
        labs = [(bar(X[0], eps[0]), bar(X[1], eps[1])) for X in m0_basis(d1, d2)]
        for X, Xl in zip(m0_basis(d1, d2), labs):
            check(on_link(Xl, eps) == on_link(X), "same operators")
        inW = all(sp.simplify(Xl[0].trace() + Xl[1].trace()) == 0 for Xl in labs)
        check(inW == (eps[0] == eps[1]), "MC(eps) iff e1 = e2")
print("  at (2,2),(2,3),(3,3) for all four eps: the labels act on L as the same operators, and lie in the Letter's")
print("  h_W iff e1 = e2, i.e. exactly when the selection is that of reading A  (asserted)")

print("(d) h cap stab at (2,3), as in Thm S-stabilizer; operator pairs, so reading A is h_+ = {tr1+tr2=0}, A* is h_- = {tr1-tr2=0}:")
d1, d2 = 2, 3
st = stabilizer(d1, d2)
stv = [realvec([X[0], X[1]]) for X in st]
m0v = [realvec([X[0], X[1]]) for X in m0_basis(d1, d2)]
x_ = sp.symbols(f"s0:{len(st)}", real=True)
gen = (sum((x_[k] * st[k][0] for k in range(len(st))), zeros(d1, d1)),
       sum((x_[k] * st[k][1] for k in range(len(st))), zeros(d2, d2)))
for sign, name in [(1, "A  (h_+)"), (-1, "A* (h_-)")]:
    cond = sp.expand((gen[0].trace() + sign * gen[1].trace()) / I)
    # basis of the intersection
    lin = sp.Matrix([[sp.diff(cond, s) for s in x_]])
    ns = lin.nullspace()
    inter = [(sum((v[k] * st[k][0] for k in range(len(st))), zeros(d1, d1)),
              sum((v[k] * st[k][1] for k in range(len(st))), zeros(d2, d2))) for v in ns]
    iv = [realvec([X[0], X[1]]) for X in inter]
    equal_m0 = rank_of(iv + m0v) == rank_of(iv) == rank_of(m0v)
    # centre of the intersection: elements commuting with all of it
    cs = sp.symbols(f"c0:{len(inter)}", real=True)
    Z = (sum((cs[k] * inter[k][0] for k in range(len(inter))), zeros(d1, d1)),
         sum((cs[k] * inter[k][1] for k in range(len(inter))), zeros(d2, d2)))
    eqs = []
    for X in inter:
        for M in [Z[0] * X[0] - X[0] * Z[0], Z[1] * X[1] - X[1] * Z[1]]:
            for z in M:
                z = sp.expand(z)
                eqs += [sp.re(z), sp.im(z)]
    csol = sp.Matrix([[sp.diff(e, c) for c in cs] for e in eqs]).nullspace()
    check(len(csol) == 1, "one-dimensional centre")
    Zc = (sum((csol[0][k] * inter[k][0] for k in range(len(inter))), zeros(d1, d1)),
          sum((csol[0][k] * inter[k][1] for k in range(len(inter))), zeros(d2, d2)))
    diag_pair = [sp.simplify(Zc[0][i, i] / I) for i in range(d1)] + [sp.simplify(Zc[1][i, i] / I) for i in range(d2)]
    check(Zc[0].is_diagonal() and Zc[1].is_diagonal(), "diagonal centre")
    nz = [v for v in diag_pair if v != 0][0]
    diag_pair = [sp.nsimplify(v / nz) for v in diag_pair]
    # the Letter's form uses labels: far label = operator (A) or conj(operator) (A*)
    lab = diag_pair[:d1] + ([v for v in diag_pair[d1:]] if sign == 1 else [-v for v in diag_pair[d1:]])
    opL = on_link(Zc)
    check(opL.is_diagonal(), "diag on L")
    Ldiag = [sp.simplify(opL[i, i] / I) for i in range(d1 * d2)]
    nzL = [v for v in Ldiag if v != 0][0]
    Ldiag = [sp.nsimplify(v / nzL) for v in Ldiag]
    check(Ldiag == [0, 0, 1, 0, 0, 1], "centre acts on L as 1 (x) diag(0,0,1)")
    check(len(inter) == 4, "dim 4")
    if sign == 1:
        check(equal_m0 and lab == [1, 1, -1, -1, 0], "A: intersection = m_0, centre = U")
    else:
        check((not equal_m0) and lab == [1, 1, 1, 1, -4], "A*: intersection != m_0, centre i diag(1,1;1,1,-4)")
    print(f"  {name}: dim {len(inter)}; equal to m_0: {equal_m0}; centre, as operator pair i diag{tuple(diag_pair)},"
          f" in the Letter's labels i diag{tuple(lab)}; on L it is proportional to 1 (x) diag(0,0,1)")

print("(e) operator level: is the global phase i 1_L in d rho(su+su) + d rho(m_0)?")
for (d1, d2) in TYPES + EXTRA:
    su_pairs = [(a, zeros(d2, d2)) for a in su_basis(d1)] + [(zeros(d1, d1), b) for b in su_basis(d2)]
    img = [realvec([on_link(X)]) for X in su_pairs + m0_basis(d1, d2)]
    one = realvec([I * eye(d1 * d2)])
    inside = rank_of(img + [one]) == rank_of(img)
    check(inside == (d1 != d2), "phase iff d1 != d2")
    print(f"  ({d1},{d2}): {inside}")
print("  (iff d1 != d2: the subgroup of U(H1)(x)U(H2) generated by SU(H1)(x)SU(H2) and the matchings is the image of")
print("  S(U(d1) x U(d2)) under reading A; it is a statement on L alone, with no labels)")

# ----------------------------------------------------------------------------------------------------------
# Sec. 5  readouts per reading
# ----------------------------------------------------------------------------------------------------------
head("5. Readouts: acting types, Theorem 1's selection, hypercharge (Thm 8, Cor. 9)")
buf = io.StringIO()
with contextlib.redirect_stdout(buf):
    import config_enumeration as CE
with open(os.path.join(HERE, "config_enumeration.out"), encoding="utf-8") as fh:
    committed = fh.read()
check(buf.getvalue() == committed, "config_enumeration.py reproduces config_enumeration.out")
print("config_enumeration.py imported unmodified; its printout equals config_enumeration.out")
check((len(CE.configs), len(CE.admitted), len(CE.readable)) == (76, 42, 23), "76/42/23")


def fmt(k):
    return f"centre {k[0]}: " + " ".join(f"{e}{t}" for e, t in k[1])


def select(acting):
    return [k for k in CE.readable if CE.end_rule(k) and CE.two_point(k) and all(acting(k[0], e) for e, _ in k[1])]


OBJ = (3, ((2, "R"), (2, "R"), (2, "U")))
QUB = (2, ((2, "R"), (2, "R"), (2, "U")))
pre = [k for k in CE.readable if CE.end_rule(k) and CE.two_point(k)]
check(pre == [QUB, OBJ], "end rule + two-point leave the qubit and qutrit configurations")
print("end rule + two-point, before the acting test:", [fmt(k) for k in pre])
print("so Theorem 1's selection is {object} iff the acting test fails at (2,2) and holds at (3,2) (asserted)")


def acting_from(alpha, beta):
    """canonical link action, trace condition alpha tr X_centre + beta tr X_far = 0 (H1 = centre end)."""
    def f(c, e):
        a = alpha.subs({D1: c, D2: e}) if hasattr(alpha, "subs") else alpha
        b = beta.subs({D1: c, D2: e}) if hasattr(beta, "subs") else beta
        if a == 0 and b == 0:
            return True                      # no condition: the whole centre u(1)+u(1) acts
        return b * e - a * c != 0            # generator (b e, -a c); acts on L as x + y
    return f


def mc_from(alpha, beta, c, e):
    a = alpha.subs({D1: c, D2: e}) if hasattr(alpha, "subs") else alpha
    b = beta.subs({D1: c, D2: e}) if hasattr(beta, "subs") else beta
    return a == b


def Y_link(alpha, beta):
    """at the object (centre qutrit H1, far qubit H2): primitive generator, qubit coefficient > 0; Y = link charge / 6
    when |alpha| = |beta| (a direct sum of the two defining blocks, SU(5) normalisation), else None."""
    a = alpha.subs({D1: 3, D2: 2}) if hasattr(alpha, "subs") else alpha
    b = beta.subs({D1: 3, D2: 2}) if hasattr(beta, "subs") else beta
    if a == 0 and b == 0:
        return "no u(1) singled out"
    x, y = 2 * b, -3 * a                 # x on the qutrit (centre), y on the qubit (far)
    g = sp.gcd(x, y)
    x, y = x / g, y / g
    if y < 0 or (y == 0 and x < 0):
        x, y = -x, -y
    if abs(a) != abs(b):
        return f"generator (qutrit {x}, qubit {y}); link charge {x + y}"
    return f"|Y| = {abs(Rational(x + y, 6))}"


print()
print("Readings (rows A, A*, B, C of F-CT-3 and every space F of Sec. 2):")
print("  reading / trace space          | acts (2,2) | acts (3,2) | MC (2,2),(3,2) | selected by Theorem 1             | (2,3) link")
ROWS = [("A   (h_+, trace on H1+H2)", sp.Integer(1), sp.Integer(1)),
        ("A*  (h_-, trace on H1+conj H2)", sp.Integer(1), sp.Integer(-1)),
        ("B   (trace on L)", D2, D1),
        ("C   (no condition)", sp.Integer(0), sp.Integer(0))]
ROWS += [("F = " + name, COEF[name][1], COEF[name][2]) for name, _ in SPACES]
SEL = {}
for name, al_, be_ in ROWS:
    f = acting_from(al_, be_)
    sel = select(f)
    SEL[name] = sel
    mc = (mc_from(al_, be_, 2, 2), mc_from(al_, be_, 3, 2))
    s = [fmt(k).replace("centre ", "c") for k in sel]
    print(f"  {name:32s}| {str(f(2, 2)):10s} | {str(f(3, 2)):10s} | {str(mc):14s} | {', '.join(s) if s else 'none':33s} | {Y_link(al_, be_)}")
check(SEL[ROWS[0][0]] == [OBJ] and SEL[ROWS[1][0]] == [QUB, OBJ] and SEL[ROWS[2][0]] == [] and SEL[ROWS[3][0]] == [QUB, OBJ],
      "A, A*, B, C reproduce config_enumeration.out")
# every MC-compatible nonvacuous row selects the object only; MC-compatible vacuous rows select two
for name, al_, be_ in ROWS:
    mc_all = mc_from(al_, be_, 2, 2) and mc_from(al_, be_, 3, 2)
    vac22 = (al_.subs({D1: 2, D2: 2}) if hasattr(al_, "subs") else al_) == 0 and \
            (be_.subs({D1: 2, D2: 2}) if hasattr(be_, "subs") else be_) == 0
    if mc_all and not vac22:
        check(SEL[name] == [OBJ], f"MC nonvacuous -> object only ({name})")
    if mc_all and vac22:
        check(SEL[name] == [QUB, OBJ], f"MC vacuous -> two ({name})")
    if not mc_from(al_, be_, 2, 2):
        check(SEL[name] != [OBJ], f"MC fails at (2,2) -> not object only ({name})")
print("  asserted: every row that satisfies matching closure at (2,2) and (3,2) with a nonzero condition selects the")
print("  object alone; every row with no condition at (2,2) selects two configurations; no row that fails matching")
print("  closure at (2,2) selects the object alone.")

print()
print("General [alpha:beta] fixed across types (alpha on the centre end, beta on the far end), canonical link:")
print("acting at (2,2) iff alpha != beta; at (3,2) iff 2 beta != 3 alpha")
grid = {}
for a_ in range(-4, 5):
    for b_ in range(-4, 5):
        if a_ == 0 and b_ == 0:
            continue
        sel = tuple(select(acting_from(sp.Integer(a_), sp.Integer(b_))))
        grid.setdefault(sel, []).append((a_, b_))
for sel, pts in grid.items():
    pr = sorted({(a_ // sp.gcd(a_, b_), b_ // sp.gcd(a_, b_)) if (a_, b_) != (0, 0) else (0, 0) for a_, b_ in pts})
    lab = ", ".join(fmt(k) for k in sel) if sel else "none"
    special = [p for p in pr if p in [(1, 1), (-1, -1), (2, 3), (-2, -3)]]
    print(f"  {lab:45s}: {len(pts):2d} of the 80 nonzero (alpha,beta) in [-4,4]^2" +
          (f"; projective classes include {special}" if special else ""))
check(set(grid[(OBJ,)]) == {(k, k) for k in range(-4, 5) if k != 0}, "object only iff alpha = beta")
check(set(grid[(QUB,)]) == {(2, 3), (-2, -3)}, "qubit only iff [2:3]")
print("  asserted: the object alone iff alpha = beta; the qubit configuration alone iff [alpha:beta] = [2:3]")
print("  (the ratio of reading B at (3,2)); otherwise both.  'none' does not occur for a fixed nonzero [alpha:beta].")

print()
print(f"{NCHECK} checks passed.")
print(f"runtime: {time.perf_counter() - T0:.1f} s")
