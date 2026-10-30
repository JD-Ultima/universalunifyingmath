"""Gap CT-G09: is there a uniform definition of the 172 records, with no exception clause?

Companion of formalization/chain-theory/alpha/records-uniform.md.  Everything is exact (integers and rationals;
sympy DomainMatrix over QQ for ranks and kernels), deterministic and asserted.  The only floating-point steps are
the cross-checks against records_independent.py and symmetry.py (imported, not modified) in Sec. 1, which are
rank tests with tolerance 1e-8, and the printing of 1/value (mpmath, 40 digits).

Setting.  H = C^3_{Q1} (x) C^2_{q1} (x) C^3_{Q2} (x) C^2_{q2}.  At each slot the operator basis is the identity
and the real traceless letters of ChainAlpha.lean's tracelessBasis (symmetric, antisymmetric, diag(1,..,1,-k,0..)),
which span the same complex spaces as the Gell-Mann and Pauli matrices.  Every candidate record space is a
subspace of P<=2, the span of the basis operators with at most two letters (dim 192); coordinates are taken in
that basis and the trace pairing is its diagonal Gram matrix.  Exchanges: E (plain slot swap, the Letter's and
the Lean file's), E' = E o T (F-CT-10), and the G_Omega-covariant (eps_Q)^{N_Q} (eps_q)^{N_q} E' of Sec. 5.

Sections
  0  source anchors (only with --upstream; the files are read as text)
  1  exact model; cells; cross-check with records_independent.py and symmetry.py
  2  the record spaces R_f = P2 + (f-odd single-qubit operators); the commutator spans C_L, C, C+; matchings
  3  order-definite composites (Thm 2); anticommuting letters at the qutrit (Prop. 3)
  4  uniform sign rules (Thm 4); isotypic sums (Prop. 5)
  5  identification rules (Thm 6); the covariant exchanges (Thm 7)
  6  collapse X -> X|Omega> under the transpose identity (Prop. 8)
  7  the declared family F (Prop. 9)
  8  single-letter content of Phi(P2) for natural superoperators Phi (Prop. 10)
  9  record-space variants and their values

Run (from the repository root or from this directory):
  python3 -X utf8 computation/chain-theory/alpha/records_uniform.py
  python3 -I -X utf8 records_uniform.py --upstream /path/to/chain-theory > records_uniform.out
(--upstream: the author's repository at a28db89, read as text by Sec. 0 only; without it Sec. 0 is skipped and the
check count is lower by the number of anchors.)
"""
import argparse
import contextlib
import hashlib
import importlib.util
import io
import itertools
import os
import time
from collections import Counter
from fractions import Fraction as Fr

T0 = time.perf_counter()
HERE = os.path.dirname(os.path.abspath(__file__))

import numpy as np
import mpmath as mp
import sympy as sp
from sympy import QQ
from sympy.polys.matrices import DomainMatrix

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


# ==========================================================================================================
# 0. source anchors
# ==========================================================================================================
head("0. Source anchors (alpha.tex, supplement.tex, ChainAlpha.lean at a28db89)")
ANCHORS = [
    ("alpha.tex", 60, r"Each pair has its own matching, so nothing lines up the bases of a shared qudit's copies"),
    ("alpha.tex", 60, r"three at a qubit (the Pauli matrices) and eight at a qutrit (the Gell-Mann matrices)"),
    ("alpha.tex", 62, r"Two reads fix their product jointly, and exchanging them returns the product times the phase by which they fail to commute"),
    ("alpha.tex", 62, r"two anticommuting reads thus fix one thing more than either alone, their product, with their order as its sign"),
    ("alpha.tex", 79, r"the far end can take that basis as its own only if a matching undoes the conjugation, as $\sigma_{y}$ does at a qubit and no unitary does at a qutrit"),
    ("alpha.tex", 109, r"A qutrit end has $3^{2}-1=8$ letters, the Gell-Mann matrices"),
    ("alpha.tex", 109, r"A \emph{measurement record} carries the two insertions as two letters"),
    ("alpha.tex", 109, r"Two letters at the same slot act on the same end, where their product is a multiple of the identity plus one traceless observable, a single insertion"),
    ("alpha.tex", 109, r"So the two letters of a record sit at different slots [Fig.~\ref{fig:link}(b)], except in the commutator records below"),
    ("alpha.tex", 109, r"At a qubit two reads in different bases anticommute (the Weyl phase there is $-1$)"),
    ("alpha.tex", 109, r"At the qutrit no two reads anticommute, since $AB=-BA$ would give $\det(AB)=-\det(BA)$ in odd dimension, and no such record arises"),
    ("alpha.tex", 109, r"of the two combinations only their difference reads the order sign $-1$ under it, so there are three \emph{commutator records}"),
    ("alpha.tex", 109, r"so the exchange sorts that space into its two eigenspaces, keeping both, and does not quotient it"),
    ("alpha.tex", 111, r"These are the $73$ products with both slots at one end and the three commutator records, whose two qubit letters it reorders"),
    ("alpha.tex", 113, r"which is $+1$ for a product, whose two letters sit at distinct slots and commute, and the order sign $-1$ for a commutator record"),
    ("supplement.tex", 35, r"a record: two letters, the two insertions"),
    ("supplement.tex", 70, r"the mutually unbiased eigenbases of the clock and shift operators"),
    ("supplement.tex", 78, r"Two reads that anticommute therefore fix a further record, their product, with their order as its sign."),
    ("supplement.tex", 81, r"the first part fixed whatever their order and the second by it"),
    ("supplement.tex", 81, r"For two invertible reads with $AB=\lambda BA$ at an end of dimension $d$, determinants give $\lambda^{d}=1$"),
    ("supplement.tex", 81, r"a sign is possible only at even $d$, which the qubit's reads realize"),
    ("supplement.tex", 170, r"At a qubit $\sigma_{y}$ does, $\sigma_{y}\bar{U}\sigma_{y}=(\det\bar{U})\,U$, and on traceless letters $A^{T}=-\sigma_{y}A\sigma_{y}$"),
    ("supplement.tex", 170, r"so this route is a second one and not the one the Letter takes"),
    ("supplement.tex", 300, r"A product record is a product of one letter at each of two of the four slots."),
    ("supplement.tex", 303, r"A single letter is one insertion, and three or four letters are more than two."),
    ("supplement.tex", 303, r"Two letters at one slot multiply to a multiple of the identity plus one traceless observable, a single insertion, which is not a record of two insertions, except where anticommuting letters make their commutator one."),
    ("supplement.tex", 306, r"it does not quotient it, and both eigenvalues are kept, the eigenvalue labelling the sector"),
    ("supplement.tex", 310, r"The commutator records are the three differences"),
    ("supplement.tex", 313, r"on their span it reads $+1$ at the sum and $-1$ at the difference. The definition of the known makes the order sign the sign of the record, so the record is the difference, and the sum, which reads $+1$, is not a commutator record."),
    ("supplement.tex", 313, r"Each difference is a single qubit element at each qubit slot, so it is orthogonal to the $169$ product records"),
    ("supplement.tex", 364, r"the three commutator records span a copy of $\Lambda^{2}V_{3}$"),
    ("supplement.tex", 461, r"sums and differences of the commutators both kept"),
    ("supplement.tex", 573, r"product of letters at two slots, or a commutator record"),
    ("supplement.tex", 575, r"$+1$ for a product, the order sign for a commutator record"),
    ("lean/ChainAlpha.lean", 216, r"A read can place a letter, an element of the real traceless basis of an end (the paper's Gell-Mann and"),
    ("lean/ChainAlpha.lean", 217, r"Pauli matrices, up to factors of i)"),
    ("lean/ChainAlpha.lean", 219, r"carries letters at two of the four slots and the identity at the others; a commutator record"),
    ("lean/ChainAlpha.lean", 220, r"is a single letter at the far end of one read minus its image at the far end of the other."),
    ("lean/ChainAlpha.lean", 280, r"def weylIsSign (d : Nat) : Bool := (List.range d).all (fun j => weylPhaseExp d j == d)"),
    ("lean/ChainAlpha.lean", 295, r"def orderSign (d : Nat) : Int := if weylIsSign d then -1 else 1"),
    ("lean/ChainAlpha.lean", 300, r"def commutatorCandidates (dA dB i j : Nat) : List Formal :="),
    ("lean/ChainAlpha.lean", 301, r"let c1 := commAt dA dB i j 1"),
    ("lean/ChainAlpha.lean", 302, r"let c3 := commAt dA dB i j 3"),
    ("lean/ChainAlpha.lean", 324, r"def commutatorRecs (dA dB : Nat) : List Rec :="),
    ("lean/ChainAlpha.lean", 325, r"if weylIsSign dB then"),
    ("lean/ChainAlpha.lean", 329, r"if eqZ (parity fg.1) (orderSign dB) then some ⟨fg.1, orderSign dB, some fg.2⟩ else none)))"),
    ("lean/ChainAlpha.lean", 350, r"odd ++ even ++ fixedP ++ commutatorRecs dA dB"),
    ("lean/ChainAlpha.lean", 638, r"theorem weyl_sign : weylIsSign 2 = true ∧ weylIsSign 3 = false := by decide +kernel"),
    ("lean/ChainAlpha.lean", 665, r"theorem record_count : (records 3 2).length = 172 := by decide +kernel"),
]
if ARGS.upstream:
    base = os.path.join(ARGS.upstream, "papers", "alpha")
    with open(os.path.join(base, "lean", "ChainAlpha.lean"), "rb") as fh:
        sha = hashlib.sha256(fh.read()).hexdigest()
    check(sha == "522844e911cf5f6091b603db0d1d042f50544d2ec8d10ed9336b722945c331da", "ChainAlpha.lean pin")
    print("ChainAlpha.lean sha256 matches the pin 522844e9...c331da")
    cache_src = {}
    for fn, ln, snip in ANCHORS:
        if fn not in cache_src:
            with open(os.path.join(base, fn), encoding="utf-8") as fh:
                cache_src[fn] = fh.read().split("\n")
        check(snip in cache_src[fn][ln - 1], f"anchor {fn}:{ln}")
        short = snip if len(snip) <= 74 else snip[:71] + "..."
        print(f"  {fn}:{ln}: found  \"{short}\"")
else:
    print("skipped (no --upstream given)")

# ==========================================================================================================
# 1. exact model
# ==========================================================================================================
head("1. Exact model: letters, P<=2, cells; cross-check with records_independent.py and symmetry.py")
DIMS = (3, 2, 3, 2)                      # Q1, q1, Q2, q2
QS, qS = (0, 2), (1, 3)


def Emat(n, a, b):
    m = np.zeros((n, n), dtype=np.int64)
    m[a, b] = 1
    return m


def letters(d):
    """identity, then ChainAlpha.lean's tracelessBasis: sym (j<k), anti (j<k), diag(1,..,1,-k,0,..)"""
    B, ts = [np.eye(d, dtype=np.int64)], [1]
    for j in range(d):
        for k in range(j + 1, d):
            B.append(Emat(d, j, k) + Emat(d, k, j)); ts.append(1)
    for j in range(d):
        for k in range(j + 1, d):
            B.append(Emat(d, j, k) - Emat(d, k, j)); ts.append(-1)
    for k in range(1, d):
        m = np.zeros((d, d), dtype=np.int64)
        for a in range(k):
            m[a, a] = 1
        m[k, k] = -k
        B.append(m); ts.append(1)
    return B, ts


LB = {d: letters(d) for d in (2, 3, 4)}
NORM = {d: [int(np.trace(b.T @ b)) for b in LB[d][0]] for d in LB}
TSG = {d: LB[d][1] for d in LB}
for d in LB:
    B = LB[d][0]
    check(len(B) == d * d, f"basis size d={d}")
    for i, a in enumerate(B):
        check((a.T == TSG[d][i] * a).all(), "transpose sign")
        for j, b in enumerate(B):
            if i != j:
                check(int(np.trace(a.T @ b)) == 0, "orthogonal basis")


def struct_consts(d):
    B = LB[d][0]
    n = len(B)
    C = [[[Fr(int(np.trace(B[c].T @ (B[a] @ B[b]))), NORM[d][c]) for c in range(n)] for b in range(n)] for a in range(n)]
    for a in range(n):
        for b in range(n):
            acc = np.zeros((d, d), dtype=object)
            for c in range(n):
                acc = acc + C[a][b][c] * B[c]
            check((acc == (B[a] @ B[b])).all(), "structure constants reconstruct products")
    return C


SC = {d: struct_consts(d) for d in (2, 3, 4)}
print("letters: 8 at a qutrit, 3 at a qubit (15 at d = 4, for Thm 2); structure constants exact and checked")

P2T = [t for t in itertools.product(*[range(d * d) for d in DIMS]) if sum(1 for x in t if x) <= 2]
IDX = {t: i for i, t in enumerate(P2T)}
NP = len(P2T)
GR = [int(np.prod([NORM[DIMS[s]][t[s]] for s in range(4)])) for t in P2T]
check(NP == 192, "dim P<=2 = 192")
print(f"P<=2 (at most two letters): dim {NP} = 1 + 22 + 169")


def occ(t):
    return tuple(s for s in range(4) if t[s])


def cls(t):
    o = occ(t)
    if not o:
        return "id"
    if len(o) == 1:
        return "Q" if o[0] in QS else "q"
    a, b = o
    if a in QS and b in QS:
        return "QQ"
    if a in qS and b in qS:
        return "qq"
    return "mS" if (a, b) in ((0, 1), (2, 3)) else "mC"


CLASSES = ["id", "Q", "q", "QQ", "qq", "mS", "mC"]
PRODS = ["QQ", "qq", "mS", "mC"]


def fsign(kind, t):
    """sign of the exchange on a basis tuple: kind 'E' (plain) or (eps_Q, eps_q) for (eps_Q)^NQ (eps_q)^Nq E'"""
    if kind == "E":
        return 1
    eQ, eq = kind
    sg = 1
    for s in range(4):
        sg *= TSG[DIMS[s]][t[s]]
    return sg * eQ ** sum(1 for s in QS if t[s]) * eq ** sum(1 for s in qS if t[s])


def fimage(t):
    return (t[2], t[3], t[0], t[1])


FR = {"E": "E", "E'": (1, 1), "Fm": (1, -1), "FQ": (-1, 1), "sE'": (-1, -1)}
FRAME_LABEL = {"E": "E (plain; Letter, Lean)", "E'": "E' = E o T (F-CT-10)", "Fm": "Fm = (-1)^Nq E'",
               "FQ": "(-1)^NQ E'", "sE'": "sigma E' = (-1)^N E'"}


def unit(t, c=1):
    v = [0] * NP
    v[IDX[t]] = c
    return v


def add(*vs):
    return [sum(x) for x in zip(*vs)]


def scal(c, v):
    return [c * x for x in v]


def fapply(kind, v):
    out = [0] * NP
    for i, x in enumerate(v):
        if x:
            t = P2T[i]
            out[IDX[fimage(t)]] += fsign(kind, t) * x
    return out


def make_cells(kind):
    res = {}
    seen = set()
    for t in P2T:
        if t in seen:
            continue
        u = fimage(t)
        seen |= {t, u}
        sg = fsign(kind, t)
        c = cls(t)
        if u == t:
            res.setdefault((c, sg), []).append(unit(t))
        else:
            for p in (1, -1):
                v = unit(t)
                v[IDX[u]] = p * sg
                res.setdefault((c, p), []).append(v)
    for (c, p), vs in res.items():
        for v in vs:
            check(fapply(kind, v) == scal(p, v), "cell vector is an eigenvector")
    return res


CELLS = {k: make_cells(v) for k, v in FR.items()}
CELLKEYS = [(c, p) for c in CLASSES for p in (1, -1) if not (c == "id" and p == -1)]
for k in FR:
    check(set(CELLS[k]) == set(CELLKEYS), "13 cells")
CDIM = {ck: len(CELLS["E'"][ck]) for ck in CELLKEYS}
for k in FR:
    check({ck: len(CELLS[k][ck]) for ck in CELLKEYS} == CDIM, "cell dimensions frame-independent")
check(sum(CDIM.values()) == 192, "cells add up")


def ckname(ck):
    return f"{ck[0]}{'+' if ck[1] > 0 else '-'}"


print("cells (class, exchange parity) and their dimensions, the same in every frame:")
print("  " + ", ".join(f"{ckname(ck)}:{CDIM[ck]}" for ck in CELLKEYS))
print("  (id identity; Q, q single letter at a qutrit copy / read qubit; QQ, qq, mS same-link mixed, mC cross mixed)")


# ---------------------------------------------------------------- exact linear algebra (sparse rows)
def _q(x):
    return QQ(x.numerator, x.denominator) if isinstance(x, Fr) else QQ(int(x))


def dmat(rows, ncols):
    d = {}
    for i, r in enumerate(rows):
        e = {j: _q(x) for j, x in enumerate(r) if x}
        if e:
            d[i] = e
    return DomainMatrix(d, (max(len(rows), 1), ncols), QQ)


def rank(vecs):
    vecs = [v for v in vecs if any(v)]
    if not vecs:
        return 0
    return dmat(vecs, len(vecs[0])).rank()


def basis_of(vecs):
    vecs = [v for v in vecs if any(v)]
    if not vecs:
        return []
    _, piv = dmat(vecs, len(vecs[0])).transpose().rref()
    return [vecs[j] for j in piv]


def same_space(U, V):
    r = rank(U + V)
    return rank(U) == r == rank(V)


def contains(U, V):
    return rank(U + V) == rank(U)


def inter_dim(U, V):
    return rank(U) + rank(V) - rank(U + V)


def gram_row(y):
    return [y[i] * GR[i] for i in range(NP)]


def sparse(v):
    return [(i, x) for i, x in enumerate(v) if x]


def apply_rows(rows, Ub):
    """the matrix (rows . u) for u in Ub (rows: functionals on Q^192)"""
    sU = [sparse(u) for u in Ub]
    return [[sum(r[i] * x for i, x in su) for su in sU] for r in rows]


def to_fr(c):
    n, d = sp.fraction(c)
    return Fr(int(n), int(d))


def kernel_in(Ub, rows):
    """basis of {x in span Ub : rows . x = 0} (exact)"""
    Ub = basis_of(Ub)
    if not rows or not Ub:
        return list(Ub)
    A = [r for r in apply_rows(rows, Ub) if any(r)]
    if not A:
        return list(Ub)
    ns = dmat(A, len(Ub)).nullspace().to_Matrix()
    out = []
    sU = [sparse(u) for u in Ub]
    for k in range(ns.rows):
        v = [0] * NP
        for j in range(len(Ub)):
            c = ns[k, j]
            if c != 0:
                cf = to_fr(c)
                for i, x in sU[j]:
                    v[i] += cf * x
        out.append(v)
    return out


def rank_rows_on(rows, Ub):
    if not rows or not Ub:
        return 0
    return rank(apply_rows(rows, Ub))


def cellspace(kind, keys):
    out = []
    for ck in keys:
        out += CELLS[kind][ck]
    return out


def content(kind, W):
    """dimension of W in each cell of the frame, and whether W is the sum of these intersections"""
    cont = {ck: (inter_dim(W, CELLS[kind][ck]) if W else 0) for ck in CELLKEYS}
    return cont, sum(cont.values()) == rank(W)


# ---------------------------------------------------------------- weights
SQ = {1: (1, 0, 0, 0), 2: (0, 1, 0, 0), 3: (0, 0, 1, 0), 6: (0, 0, 0, 1)}
NQc = {"id": 0, "Q": 1, "q": 0, "QQ": 2, "qq": 0, "mS": 1, "mC": 1}
Nqc = {"id": 0, "Q": 0, "q": 1, "QQ": 0, "qq": 2, "mS": 1, "mC": 1}
SIG = {c: (-1) ** (NQc[c] + Nqc[c]) for c in CLASSES}


def w1(ck):
    """operator form M = (1+(r3-1) S Pi_X)(1+(r2-1) S Pi_Y), S = (1 + sigma f)/2, Pi_X = [NQ Nq = 0], Pi_Y = [Nq = 1]
    (operator-origin.md Prop. 2, with the frame's exchange f in place of E')"""
    c, p = ck
    S = 1 if SIG[c] * p == 1 else 0
    return (3 if S and NQc[c] * Nqc[c] == 0 else 1) * (2 if S and Nqc[c] == 1 else 1)


def w2(ck):
    """the Letter's word rule, extended to single-letter composites: X acts on two letters of one type (same end;
    at one slot by reordering, so a single-slot composite with order sign -1 survives it); Y acts on what sits at
    the far ends (it moves a single qubit letter q1 <-> q2; it does not act on the shared qutrit's copies);
    survival = returned times own sign.  None for the identity (no letters)."""
    c, p = ck
    if c == "id":
        return None
    if c in ("QQ", "qq"):
        return 3 if p == 1 else 1
    if c in ("mS", "mC"):
        return 2 if p == 1 else 1
    if c == "q":
        return 6 if p == -1 else 3
    return 3


for ck in CELLKEYS:
    if ck[0] in PRODS:
        check(w1(ck) == w2(ck), "the two weight laws agree on products")
check(w1(("q", -1)) == w2(("q", -1)) == 6, "and on the commutator records")


def numerator(cont, wf):
    num = [0, 0, 0, 0]
    for ck, n in cont.items():
        if n:
            k = wf(ck)
            if k is None:
                return None
            for i in range(4):
                num[i] += n * SQ[k][i]
    return tuple(num)


mp.mp.dps = 40


def inv_value(num):
    a, b, c, d = num
    return mp.mpf(31104) / (a + b * mp.sqrt(2) + c * mp.sqrt(3) + d * mp.sqrt(6))


def fmt_num(num):
    if num is None:
        return "n/a"
    parts = [str(num[0])] if num[0] else []
    for co, r in ((num[1], "r2"), (num[2], "r3"), (num[3], "r6")):
        if co:
            parts.append(f"{co}{r}")
    return "+".join(parts) if parts else "0"


def fmt_inv(num):
    return "n/a" if num is None else mp.nstr(inv_value(num), 9)


LETTER = (79, 48, 42, 3)


def letter_classes(cont):
    """(odd products, even mixed, even same-end, odd single-qubit) and the dimension outside these classes"""
    odd = sum(cont[(c, -1)] for c in PRODS)
    evm = cont[("mS", 1)] + cont[("mC", 1)]
    evs = cont[("QQ", 1)] + cont[("qq", 1)]
    com = cont[("q", -1)]
    rest = cont[("id", 1)] + cont[("Q", 1)] + cont[("Q", -1)] + cont[("q", 1)]
    return (odd, evm, evs, com), rest


# ---------------------------------------------------------------- operators on H (exact integers)
def opmat(t):
    m = np.array([[1]], dtype=np.int64)
    for s in range(4):
        m = np.kron(m, LB[DIMS[s]][0][t[s]])
    return m


OPM = [opmat(t) for t in P2T]
PHI = np.zeros(36, dtype=np.int64)               # |Omega> unnormalized: sum_ij |i j i j>  (order Q1 q1 Q2 q2)
for i in range(3):
    for j in range(2):
        PHI[np.ravel_multi_index((i, j, i, j), DIMS)] = 1
LKET = np.zeros(36, dtype=np.int64)              # |L1 L2> = |w2>_{Q1 q1} |w2>_{Q2 q2}, support levels 0, 1
for i in range(2):
    for j in range(2):
        LKET[np.ravel_multi_index((i, i, j, j), DIMS)] = 1


def vec_to_op(v):
    m = np.zeros((36, 36), dtype=object)
    for i, x in sparse(v):
        m = m + x * OPM[i].astype(object)
    return m


spec = importlib.util.spec_from_file_location("sym", os.path.join(HERE, "symmetry.py"))
sym = importlib.util.module_from_spec(spec)
with contextlib.redirect_stdout(io.StringIO()):
    spec.loader.exec_module(sym)          # rebuilds the Letter's records (records_independent.py) and the covariant ones
rec = sym.rec


def num_rank(mats):
    M = np.array([np.asarray(m, dtype=complex).reshape(-1) for m in mats])
    return np.linalg.matrix_rank(M, tol=1e-8)


RECS = {k: cellspace(k, [(c, p) for c in PRODS for p in (1, -1)] + [("q", -1)]) for k in FR}
mine_E = [vec_to_op(v).astype(float) for v in RECS["E"]]
mine_Ep = [vec_to_op(v).astype(float) for v in RECS["E'"]]
check(num_rank(mine_E) == 172 and num_rank([r[1] for r in rec.records] + mine_E) == 172, "R_E = Letter's span")
check(num_rank(mine_Ep) == 172 and num_rank(list(sym.cov_ops) + mine_Ep) == 172, "R_E' = covariant span")
mine_CL = [vec_to_op(v).astype(float) for v in CELLS["E"][("q", -1)]]
mine_C = [vec_to_op(v).astype(float) for v in CELLS["E'"][("q", -1)]]
check(num_rank(mine_CL + list(rec.comm)) == 3, "C_L = the Letter's commutator records")
check(num_rank(mine_C + list(sym.comm_cov)) == 3, "C = symmetry.py's covariant commutator-type records")
print("cross-check (records_independent.py and symmetry.py, imported unmodified; numerical rank tests, tol 1e-8):")
print("  R_E = span of the Letter's 172 records; R_E' = span of the 172 covariant records of F-CT-10; their")
print("  single-qubit parts are the Letter's three commutator records and symmetry.py's D_x, D_z, S_y.")

# ==========================================================================================================
# 2. record spaces, commutator spans, matchings
# ==========================================================================================================
head("2. Record spaces R_f = P2 + (f-odd single-qubit operators); C_L, C, C+; matchings")
P2 = cellspace("E'", [(c, p) for c in PRODS for p in (1, -1)])
C_L, C_Lp = CELLS["E"][("q", -1)], CELLS["E"][("q", 1)]
C_c, C_p = CELLS["E'"][("q", -1)], CELLS["E'"][("q", 1)]
QSING = cellspace("E'", [("Q", 1), ("Q", -1)])
qSING = cellspace("E'", [("q", 1), ("q", -1)])
IDV = [unit((0, 0, 0, 0))]
check(rank(P2) == 169 and same_space(P2, cellspace("E", [(c, p) for c in PRODS for p in (1, -1)])), "P2")
Dlist = []                     # the Letter's / Lean's commutator records [b_mu, b_nu]_{q1} - [b_mu, b_nu]_{q2}
for mu, nu in itertools.combinations(range(1, 4), 2):
    v = [0] * NP
    for c in range(4):
        co = SC[2][mu][nu][c] - SC[2][nu][mu][c]
        if co:
            v[IDX[(0, c, 0, 0)]] += co
            v[IDX[(0, 0, 0, c)]] -= co
    check(len(sparse(v)) == 2, "one letter at q1 minus the same letter at q2")
    Dlist.append(v)
check(same_space(Dlist, C_L), "commutator differences span C_L = E-odd single-qubit")
check(same_space(CELLS["Fm"][("q", -1)], C_p), "Fm-odd single-qubit = E'-even = C+")
print("C_L = E-odd single-qubit = span of the Letter's [s_mu,s_nu]_q1 - [s_mu,s_nu]_q2;  C = E'-odd;  C+ = E'-even:")
for (n1, X), (n2, Y) in itertools.combinations((("C_L", C_L), ("C", C_c), ("C+", C_p)), 2):
    print(f"  dim({n1} cap {n2}) = {inter_dim(X, Y)}")
check(inter_dim(C_L, C_c) == 2 and inter_dim(C_L, C_p) == 1 and inter_dim(C_c, C_p) == 0, "intersections 2, 1, 0")
check(same_space(C_L + C_Lp, qSING) and same_space(C_c + C_p, qSING), "both are splittings of the 6 single-qubit")

LIE = []          # Lie algebra of G_Omega: K = b^{(slot 1)} - (b^T)^{(slot 2)}, at the join and at the qubit pair
for d, (s1, s2) in ((3, (0, 2)), (2, (1, 3))):
    for a in range(1, d * d):
        v = unit(tuple(a if s == s1 else 0 for s in range(4)))
        v[IDX[tuple(a if s == s2 else 0 for s in range(4))]] = -TSG[d][a]
        LIE.append(v)
LIE3, LIE2 = LIE[:8], LIE[8:]
check(rank(LIE) == 11, "dim Lie(G_Omega) = 11")
check(same_space(LIE, cellspace("E'", [("Q", -1), ("q", -1)])), "Lie(G_Omega) = E'-odd single-letter operators")
check(same_space(LIE2, C_c), "C = matching generators of the qubit pair")
check(inter_dim(C_L, LIE) == 2 and inter_dim(C_Lp, LIE) == 1, "plain frame: two differences and one sum are matchings")
print("Lie(G_Omega) = E'-odd single-letter operators (8 at the join Q1Q2, 3 at the qubit pair q1q2); C is its qubit")
print("  part.  Plain frame: dim(C_L cap Lie) = 2 (the x and z differences), dim(C_L+ cap Lie) = 1 (the y sum).")
for k in FR:
    cont, gr = content(k, RECS[k])
    lc, rest = letter_classes(cont)
    n1, n2 = numerator(cont, w1), numerator(cont, w2)
    check(gr and lc == LETTER and rest == 0 and n1 == n2 == LETTER, f"R_{k}: classes and value")
    print(f"  R_f, f = {FRAME_LABEL[k]:24s}: classes {lc}, numerator {fmt_num(n1)}, 1/value {fmt_inv(n1)}")
check(same_space(RECS["E'"], RECS["FQ"]) and same_space(RECS["Fm"], RECS["sE'"]), "R_f depends on eps_q only")
check(same_space(RECS["Fm"], P2 + C_p), "R_Fm = P2 + C+")
check(not same_space(RECS["E"], RECS["E'"]) and not same_space(RECS["E"], RECS["Fm"])
      and not same_space(RECS["E'"], RECS["Fm"]), "three distinct spaces")
print("  R_E = P2 + C_L,  R_E' = R_(-1)^NQ E' = P2 + C,  R_Fm = R_sigma E' = P2 + C+: three distinct spaces, same value")

# ==========================================================================================================
# 3. order-definite composites; anticommuting letters
# ==========================================================================================================
head("3. Order-definite composites (Thm 2); anticommuting letters at the qutrit (Prop. 3)")
print("at one end of dimension d: dim m(Sym^2 V), dim m(Lambda^2 V), dim of the intersection, dim traceless part of m(Sym^2 V)")
for d in (2, 3, 4):
    n = d * d
    C = SC[d]
    s_ = [[C[a][b][c] + C[b][a][c] for c in range(n)] for a in range(1, n) for b in range(a, n)]
    a_ = [[C[a][b][c] - C[b][a][c] for c in range(n)] for a in range(1, n) for b in range(a + 1, n)]
    rs, ra = rank(s_), rank(a_)
    ri = rs + ra - rank(s_ + a_)
    tl = rank([r[1:] for r in s_])
    print(f"  d = {d}: {rs:2d}, {ra:2d}, {ri:2d}, {tl:2d}")
    check(ra == n - 1, "the commutators span the traceless operators")
    if d == 2:
        check(rs == 1 and ri == 0 and tl == 0, "qubit: the anticommutator of two letters is a scalar")
    else:
        check(rs == n and ri == n - 1 and tl == n - 1, f"d = {d}: every traceless composite is also a symmetric one")

LET = [(s, a) for s in range(4) for a in range(1, DIMS[s] ** 2)]
check(len(LET) == 22, "22 letters")


def letter_vec(x):
    return unit(tuple(x[1] if s == x[0] else 0 for s in range(4)))


def mprod(x, y):
    (s, a), (t, b) = x, y
    v = [0] * NP
    if s != t:
        tt = [0, 0, 0, 0]
        tt[s], tt[t] = a, b
        v[IDX[tuple(tt)]] = 1
    else:
        d = DIMS[s]
        for c in range(d * d):
            co = SC[d][a][b][c]
            if co:
                tt = [0, 0, 0, 0]
                tt[s] = c
                v[IDX[tuple(tt)]] += co
    return v


MALL = [mprod(x, y) for x, y in itertools.product(LET, LET)]
MSYM, MANT = [], []
for i, x in enumerate(LET):
    for y in LET[i:]:
        MSYM.append(add(mprod(x, y), mprod(y, x)))
        if y != x:
            MANT.append(add(mprod(x, y), scal(-1, mprod(y, x))))
check(rank(MALL) == 192, "m(V (x) V) = P<=2")
check(same_space(MSYM, IDV + P2 + QSING), "m(Sym^2 V) = 1 + P2 + V_Q1 + V_Q2")
check(same_space(MANT, QSING + qSING), "m(Lambda^2 V) = the 22 single letters")
check(inter_dim(MSYM, MANT) == 16 and contains(MSYM, QSING) and contains(MANT, QSING), "A = V_Q1 + V_Q2")
D = IDV + P2 + qSING
D0 = P2 + qSING
DP, DM = IDV + P2, qSING
check(contains(MSYM, DP) and contains(MANT, DM) and inter_dim(MSYM, DM) == 0 and inter_dim(MANT, DP) == 0,
      "reversal acts as +1 on D_+ and -1 on D_-")
for ck in CELLKEYS:
    if ck[0] in ("id", "q") or ck[0] in PRODS:
        check(SIG[ck[0]] == (1 if ck[0] != "q" else -1), "on D the own sign (-1)^N is the order sign")
print("m(V(x)V) = P<=2 (192); m(Sym^2 V) = 1 + P2 + V_Q1 + V_Q2 (186); m(Lambda^2 V) = the 22 single letters;")
print("A = m(Sym^2 V) cap m(Lambda^2 V) = V_Q1 + V_Q2 (16).  The reversal of the two insertions descends to")
print("m(V(x)V)/A; with the trace pairing, D = A^perp = 1 + P2 + V_q1 + V_q2 (176), D_+ = 1 + P2, D_- = V_q1 + V_q2,")
print("and on D the order sign equals the own sign (-1)^N.  D_0 = D cap 1^perp = P2 + V_q1 + V_q2 (175).")
ONEPER = [mprod(x, y) for x in LET for y in LET if x[0] in (0, 1) and y[0] in (2, 3)]
ONEPER_CELLS = [(c, p) for c in ("QQ", "qq", "mC") for p in (1, -1)]
check(same_space(ONEPER, cellspace("E'", ONEPER_CELLS)) and rank(ONEPER) == 121, "one letter per read: 121")
SAMEREAD = [mprod(x, y) for x in LET for y in LET if (x[0] in (0, 1)) == (y[0] in (0, 1))]
check(same_space(SAMEREAD, IDV + QSING + qSING + cellspace("E'", [("mS", 1), ("mS", -1)])), "both from one read")
print("one letter per read: QQ + qq + cross-mixed (121); both letters from one read: 1 + singles + same-link mixed (71)")

I_ = sp.I
W3 = sp.Rational(-1, 2) + sp.sqrt(3) * I_ / 2


def gell_mann():
    E_ = lambda a, b: sp.Matrix(3, 3, lambda i, j: 1 if (i, j) == (a, b) else 0)
    return [E_(0, 1) + E_(1, 0), -I_ * E_(0, 1) + I_ * E_(1, 0), E_(0, 0) - E_(1, 1),
            E_(0, 2) + E_(2, 0), -I_ * E_(0, 2) + I_ * E_(2, 0), E_(1, 2) + E_(2, 1),
            -I_ * E_(1, 2) + I_ * E_(2, 1), (E_(0, 0) + E_(1, 1) - 2 * E_(2, 2)) / sp.sqrt(3)]


def pauli():
    return [sp.Matrix([[0, 1], [1, 0]]), sp.Matrix([[0, -I_], [I_, 0]]), sp.Matrix([[1, 0], [0, -1]])]


def weyl(d):
    w = W3 if d == 3 else sp.Integer(-1)
    X = sp.Matrix(d, d, lambda i, j: 1 if i == (j + 1) % d else 0)
    Z = sp.diag(*[w ** j for j in range(d)])
    return [sp.expand(X ** a * Z ** b) for a in range(d) for b in range(d) if (a, b) != (0, 0)]


def cayley_basis(d):
    """a generic orthogonal basis of traceless Hermitian letters: the Hermitian letters S, iA, D mixed by a rational
    Cayley transform O with O^T G O = G (G = their Gram matrix tr(h_a h_b)), so the new letters are Hermitian,
    traceless and orthogonal with the same norms"""
    B = [sp.Matrix(b) * (I_ if TSG[d][a + 1] == -1 else 1) for a, b in enumerate(LB[d][0][1:])]
    n = len(B)
    G = sp.diag(*[(b * b).trace() for b in B])
    A = sp.Matrix(n, n, lambda i, j: 0 if i == j else ((i * 7 + j * 3) % 5 - 2))
    A = A - A.T
    K = G.inv() * A
    O = (sp.eye(n) - K).inv() * (sp.eye(n) + K)
    check((O.T * G * O - G).applyfunc(sp.simplify) == sp.zeros(n, n), "Cayley transform preserves the pairing")
    return [sum((O[a, j] * sp.Matrix(B[a]) for a in range(n)), sp.zeros(d, d)) for j in range(n)]


def anticommuting(L):
    pairs = [(i, j) for i, j in itertools.combinations(range(len(L)), 2)
             if (L[i] * L[j] + L[j] * L[i]).applyfunc(sp.expand) == sp.zeros(*L[i].shape)]
    comms = [(L[i] * L[j] - L[j] * L[i]).applyfunc(sp.expand) for i, j in pairs]
    return pairs, (sp.Matrix([list(c) for c in comms]).rank(simplify=True) if comms else 0)


GM = gell_mann()
for nm, L, expect in (("Gell-Mann (the Letter's letters)", GM, (5, 4)),
                      ("Lean real basis, d = 3", [sp.Matrix(b) for b in LB[3][0][1:]], (5, 4)),
                      ("Weyl X^a Z^b, d = 3", weyl(3), (0, 0)),
                      ("generic rational orthogonal, d = 3", cayley_basis(3), (0, 0)),
                      ("Pauli", pauli(), (3, 3)),
                      ("Lean real basis, d = 2", [sp.Matrix(b) for b in LB[2][0][1:]], (3, 3)),
                      ("Weyl X^a Z^b, d = 2", weyl(2), (3, 3)),
                      ("generic rational orthogonal, d = 2", cayley_basis(2), (3, 3))):
    pairs, ds = anticommuting(L)
    check((len(pairs), ds) == expect, f"anticommuting pairs: {nm}")
    extra = f"; pairs {[(i + 1, j + 1) for i, j in pairs]}" if nm.startswith("Gell") else ""
    print(f"  {nm:36s}: {len(pairs)} anticommuting pairs, their commutators span dim {ds}{extra}")
check([k + 1 for k, m in enumerate(GM) if m.det() == 0] == [1, 2, 3, 4, 5, 6, 7], "lambda_1..lambda_7 singular")
pairs, _ = anticommuting(GM)
cm = sp.Matrix([list((GM[i] * GM[j] - GM[j] * GM[i]).applyfunc(sp.expand)) for i, j in pairs])
span4 = sp.Matrix([list(GM[0]), list(GM[1]), list(GM[2]), list(GM[7])])
check(sp.Matrix.vstack(cm, span4).rank(simplify=True) == 4 == span4.rank(simplify=True), "span lambda_1,2,3,8")
print("  lambda_1..lambda_7 are singular (det 0), lambda_8 is not; the five anticommuting pairs' commutators span")
print("  {lambda_1, lambda_2, lambda_3, lambda_8}.  The det argument of alpha.tex:109 excludes invertible reads only.")

# ==========================================================================================================
# 4. sign rules; isotypic sums
# ==========================================================================================================
head("4. Uniform sign rules (Thm 4); isotypic sums (Prop. 5)")


def coset_sets():
    """subsets of the character group of <sigma, f> = Z2^2 cut out by one sign test (and the whole group)"""
    chars = list(itertools.product((1, -1), repeat=2))
    sets = {frozenset(chars)}
    for g in ((1, 0), (0, 1), (1, 1)):
        for eps in (1, -1):
            sets.add(frozenset(ch for ch in chars if (ch[0] if g[0] else 1) * (ch[1] if g[1] else 1) == eps))
    return sets


COS = coset_sets()
check(len(COS) == 7, "seven uniform sign rules (no test, or one sign test)")
for k in ("E", "E'", "Fm"):
    ucl = {(1, 1): cellspace(k, [(c, 1) for c in PRODS]), (1, -1): cellspace(k, [(c, -1) for c in PRODS]),
           (-1, 1): CELLS[k][("q", 1)], (-1, -1): CELLS[k][("q", -1)]}
    dims = {ch: rank(v) for ch, v in ucl.items()}
    check(dims == {(1, 1): 90, (1, -1): 79, (-1, 1): 3, (-1, -1): 3}, "uniform classes on D_0")
    hits = [frozenset(sub) for r in range(1, 5) for sub in itertools.combinations(sorted(ucl), r)
            if sum(dims[c] for c in sub) == 172]
    rset = frozenset({(1, 1), (1, -1), (-1, -1)})
    check(len(hits) == 2 and rset in hits and all(h not in COS for h in hits), "172 only for two non-coset unions")
    check(same_space(sum((ucl[c] for c in rset), []), RECS[k]), "R_f = {sigma = +1} cup {sigma f = +1}")
    if k == "E":
        print(f"On D_0, classes by (own sign sigma, exchange parity f): (+,+) 90, (+,-) 79, (-,+) 3, (-,-) 3.")
        print(f"  the seven uniform sign rules give dims {sorted(sum(dims[c] for c in s) for s in COS)}")
    print(f"  frame {FRAME_LABEL[k]:24s}: unions of dim 172: {[sorted(h, reverse=True) for h in hits]}; neither is"
          f" one sign test; R_f = (+,+)(+,-)(-,-) = 'sigma = +1 or sigma f = +1'")
for nm, with_id, with_Q in (("P<=2,0", False, True), ("P<=2", True, True), ("D", True, False)):
    for k in ("E", "E'", "Fm"):
        ucl = {(1, 1): cellspace(k, [(c, 1) for c in PRODS]) + (IDV if with_id else []),
               (1, -1): cellspace(k, [(c, -1) for c in PRODS]),
               (-1, 1): CELLS[k][("q", 1)] + (CELLS[k][("Q", 1)] if with_Q else []),
               (-1, -1): CELLS[k][("q", -1)] + (CELLS[k][("Q", -1)] if with_Q else [])}
        dims = {ch: rank(v) for ch, v in ucl.items()}
        alld = {sum(dims[c] for c in sub) for r in range(1, 5) for sub in itertools.combinations(sorted(ucl), r)}
        cosd = {sum(dims[c] for c in s) for s in COS}
        check(172 not in cosd, f"no sign rule of dim 172 on {nm}")
        if nm == "P<=2,0":
            check(172 not in alld, "on P<=2,0 no union of the four classes has dim 172")
        if k == "E'":
            print(f"  on {nm:7s}: classes {[dims[c] for c in sorted(ucl, reverse=True)]}; sign-rule dims {sorted(cosd)};"
                  f" 172 is {'' if 172 in alld else 'not '}a union of classes")


def ad_vec(K, v):
    """[K, X] for K a combination of single-letter tuples and X in P<=2 (coordinates); exact, stays in P<=2"""
    out = [0] * NP
    sK = [(occ(P2T[j])[0], P2T[j], k) for j, k in sparse(K) if occ(P2T[j])]      # the identity commutes
    for i, x in sparse(v):
        t = P2T[i]
        for s, tk, k in sK:
            a = tk[s]
            d = DIMS[s]
            for c in range(d * d):
                co = SC[d][a][t[s]][c] - SC[d][t[s]][a][c]
                if co:
                    tt = list(t)
                    tt[s] = c
                    out[IDX[tuple(tt)]] += k * x * co
    return out


TRSQ = {d: [int(np.trace(b @ b)) for b in LB[d][0]] for d in (2, 3)}


def casimir(gens, d, v):
    out = [0] * NP
    for a, K in enumerate(gens):
        w = ad_vec(K, ad_vec(K, v))
        c = Fr(1, TRSQ[d][a + 1])
        out = [o + c * x for o, x in zip(out, w)]
    return out


EV3 = {Fr(0): "1", Fr(6): "8", Fr(12): "10+10b", Fr(16): "27"}     # trace-form normalization tr(T_a T_b) = delta
EV2 = {Fr(0): "1", Fr(4): "3", Fr(12): "5"}
iso = Counter()
iso_cells = {}
for ck in CELLKEYS:
    V = CELLS["E'"][ck]
    c3 = [casimir(LIE3, 3, v) for v in V]
    c2 = [casimir(LIE2, 2, v) for v in V]
    for l3, n3 in EV3.items():
        for l2, n2 in EV2.items():
            cols = [add(a, scal(-l3, v)) + add(b, scal(-l2, v)) for a, b, v in zip(c3, c2, V)]
            dimk = len(V) - rank(cols)
            if dimk:
                iso[(n3, n2, ck[1])] += dimk
                iso_cells[(ck, n3, n2)] = dimk
check(sum(iso.values()) == 192, "the Casimir eigenspaces exhaust P<=2")
print("isotypic decomposition of P<=2 under SU(3) x SU(2) x <E'> (quadratic Casimirs, exact; 10 and 10b share one):")
comps = []
for (n3, n2, p), dm_ in sorted(iso.items(), key=lambda z: (-z[0][2], z[0][0], z[0][1])):
    where = ", ".join(f"{ckname(ck)} {n}" for (ck, a3, a2), n in iso_cells.items() if (a3, a2) == (n3, n2) and ck[1] == p)
    print(f"  ({n3},{n2})_{'+' if p > 0 else '-'}: dim {dm_:2d}  [{where}]")
    if n3 == "10+10b":
        comps += [("10", n2, p, dm_ // 2), ("10b", n2, p, dm_ // 2)]
    else:
        comps.append((n3, n2, p, dm_))
check(len(comps) == 11, "11 isotypic components")
check(iso_cells[(("q", -1), "1", "3")] == 3 and iso_cells[(("qq", -1), "1", "3")] == 3, "(1,3)_- = qq- + C")
check(iso_cells[(("Q", -1), "8", "1")] == 8 and iso_cells[(("QQ", -1), "8", "1")] == 8, "(8,1)_- = QQ- 8 + Q-")
sums172 = [sub for r in range(1, 12) for sub in itertools.combinations(range(11), r)
           if sum(comps[i][3] for i in sub) == 172]
check(len(sums172) == 1, "exactly one sum of isotypic components has dim 172")
miss = [comps[i][:3] for i in range(11) if i not in sums172[0]]
check(sorted(m[0] for m in miss) == ["10", "10b"], "it omits (10,1) and (10b,1)")
print("  of the 2^11 sums of isotypic components, exactly one has dim 172: all but (10,1)_- and (10b,1)_- (it keeps")
print("  the identity and all 22 single letters and drops 20 odd QQ records).  R_E' meets (8,1)_- (dim 16) in dim 8,")
print("  (8,1)_+ in 8 and (1,1)_+ in 2 of 3: it is not a sum of isotypic components.")
check(inter_dim(RECS["E'"], CELLS["E'"][("Q", -1)]) == 0, "R_E' contains no single qutrit letter")


def iso_of_record_space(kind):
    out = Counter()
    keys = [(c, p) for c in PRODS for p in (1, -1)] + [("q", -1)]
    for ck in keys:
        V = CELLS[kind][ck]
        c3 = [casimir(LIE3, 3, v) for v in V]
        c2 = [casimir(LIE2, 2, v) for v in V]
        for l3, n3 in EV3.items():
            for l2, n2 in EV2.items():
                dimk = len(V) - rank([add(a, scal(-l3, v)) + add(b, scal(-l2, v)) for a, b, v in zip(c3, c2, V)])
                if dimk:
                    out[(n3, n2, ck[1])] += dimk
    return out


ISO_Ep, ISO_Fm = iso_of_record_space("E'"), iso_of_record_space("Fm")
check(ISO_Ep == ISO_Fm and sum(ISO_Fm.values()) == 172, "R_Fm under <Fm> has the isotypic table of R_E' under <E'>")
check(ISO_Ep[("1", "3", -1)] == 6 and ISO_Ep[("8", "3", 1)] == 48 and ISO_Ep[("10+10b", "1", -1)] == 20, "F-CT-10 table")
print("  R_Fm under SU(3) x SU(2) x <Fm> has the same isotypic table as R_E' under <E'> (F-CT-10 item 4): "
      + ", ".join(f"({a},{b})_{'+' if p > 0 else '-'} {n}" for (a, b, p), n in sorted(ISO_Fm.items(), key=lambda z: (-z[0][2], z[0][0], z[0][1]))))

# ==========================================================================================================
# 5. identification rules; covariant exchanges
# ==========================================================================================================
head("5. Identification rules (Thm 6); the covariant exchanges (Thm 7)")
IDEN = IDV + QSING + qSING                         # image of the read links' identification algebras = P<=1
check(rank(IDEN) == 23, "I = P<=1 (dim 23)")
AD = {}


def ad_cols(Xv):
    key = tuple(Xv)
    if key not in AD:
        AD[key] = {i: sparse(ad_vec(Xv, unit(P2T[i]))) for i in range(NP)}
    return AD[key]


def comp_space(kind):
    """Comp_f = {X in I : f o ad_X o f = ad_X on P<=2}; f is an involutive signed permutation"""
    Fs = [(IDX[fimage(P2T[i])], fsign(kind, P2T[i])) for i in range(NP)]
    feats = []
    for Xv in IDEN:
        cols = ad_cols(Xv)
        diff = Counter()
        for i, nz in cols.items():
            for j, c in nz:
                diff[(i, j)] -= c
            fi, si = Fs[i]
            for j, c in cols[fi]:
                fj, sj = Fs[j]
                diff[(i, fj)] += si * sj * c
        feats.append({kk: vv for kk, vv in diff.items() if vv})
    keys = sorted(set().union(*[set(f) for f in feats]))
    if not keys:
        return list(IDEN)
    rows = [[f.get(kk, 0) for f in feats] for kk in keys]
    ns = dmat(rows, len(IDEN)).nullspace().to_Matrix()
    out = []
    for r in range(ns.rows):
        v = [0] * NP
        for j in range(len(IDEN)):
            if ns[r, j] != 0:
                v = add(v, scal(to_fr(ns[r, j]), IDEN[j]))
        out.append(v)
    return out


KAP = [[int(x) for x in (OPM[i] @ PHI)] for i in range(NP)]          # kappa(e_i) = e_i |Omega>
KAPROWS = [[KAP[i][r] for i in range(NP)] for r in range(36)]
NPHI = int(PHI @ PHI)
EXPV = [int(PHI @ OPM[i] @ PHI) for i in range(NP)]
EIGROWS = [[KAP[i][r] * NPHI - EXPV[i] * int(PHI[r]) for i in range(NP)] for r in range(36)]
STAB = kernel_in(IDEN, EIGROWS)
check(same_space(STAB, IDV + LIE), "stab_Omega = {X in P<=1: X|Omega> in C|Omega>} = 1 + Lie(G_Omega)")
COMP = {k: comp_space(FR[k]) for k in FR}
INV = {k: cellspace(k, [("id", 1), ("Q", 1), ("q", 1)]) for k in FR}
check(same_space(COMP["E"], INV["E"]) and not same_space(COMP["E"], STAB), "Comp_E = I^E != stab_Omega")
for k in ("E'", "Fm", "FQ", "sE'"):
    check(same_space(COMP[k], STAB), f"Comp_{k} = stab_Omega")
check(not same_space(INV["E'"], COMP["E'"]), "I^E' != Comp_E'")
check(same_space(INV["Fm"], cellspace("E'", [("id", 1), ("Q", 1), ("q", -1)])), "I^Fm")
print("I = image of the identification algebras of the two read links = P<=1 (identity + 22 single letters, dim 23).")
print("stab_Omega = {X in I : X|Omega> in C|Omega>} = 1 + Lie(G_Omega) (dim 12).")
print("Comp_f = {X in I : f ad_X f = ad_X}:  Comp_E = I^E (the E-even part, 12);  Comp_f = stab_Omega for each of the")
print("  four covariant exchanges;  I^f = f-even part of I, and I^E' != Comp_E'.")
br = ad_vec(C_p[0], C_p[1])
check(any(br) and not contains(INV["E'"], [br]) and contains(LIE, [br]), "I^E' is not closed under brackets")
check(all(contains(IDV + LIE, [ad_vec(x, y)]) for x in LIE for y in LIE), "1 + Lie(G_Omega) is a Lie algebra")
print("  I^E' is not a Lie algebra: the bracket of two E'-even qubit identifications is E'-odd (a matching).")
print("  the part of I^f that D sees, I^f cap D = 1 + (f-even single-qubit), closed under brackets? equal to Comp_f cap D?")
for k in FR:
    part = cellspace(k, [("id", 1), ("q", 1)])
    closed = all(contains(part, [ad_vec(x, y)]) for x in part for y in part)
    eqc = same_space(part, kernel_in(COMP[k], [gram_row(u) for u in QSING]))
    check(closed == eqc == (k in ("E", "Fm", "sE'")), f"I^f cap D is the Lie algebra Comp_f cap D iff eps_q = -1 or E ({k})")
    print(f"    {FRAME_LABEL[k]:24s}: closed {closed}; = Comp_f cap D {eqc}")

RULES = {}
for k in FR:
    RULES[(k, "inv")] = kernel_in(D, [gram_row(y) for y in INV[k]])
    RULES[(k, "comp")] = kernel_in(D, [gram_row(y) for y in COMP[k]])
    check(same_space(RULES[(k, "inv")], RECS[k]), f"D cap (I^f)^perp = R_f ({k})")
RSTAB = kernel_in(D, [gram_row(y) for y in STAB])
check(same_space(RSTAB, P2 + C_p), "D cap stab_Omega^perp = P2 + C+")
check(same_space(RULES[("E", "comp")], RECS["E"]), "D cap Comp_E^perp = R_E (the Letter's 172)")
for k in ("E'", "FQ"):
    check(same_space(RULES[(k, "comp")], RSTAB) and not same_space(RSTAB, RECS[k]), f"{k}: rule gives C+, not C")
for k in ("Fm", "sE'"):
    check(same_space(RULES[(k, "comp")], RECS[k]), f"{k}: rule gives R_f")
check(rank(kernel_in(IDV + P2 + QSING + qSING, [gram_row(y) for y in STAB])) == 180, "without D: dim 180")
print("Rules on D:")
for k in FR:
    for rule in ("inv", "comp"):
        W = RULES[(k, rule)]
        cont, gr = content(k, W)
        lc, rest = letter_classes(cont)
        check(gr, "graded")
        print(f"  f = {k:3s} D cap ({'I^f' if rule == 'inv' else 'Comp_f'})^perp: dim {rank(W)}, classes {lc}"
              f"{' + ' + str(rest) + ' outside' if rest else ''}, = R_f: {same_space(W, RECS[k])};"
              f" 1/value {fmt_inv(numerator(cont, w1))} (M_L), {fmt_inv(numerator(cont, w2))} (word rule)")
print("  D cap stab_Omega^perp = P2 + C+ (every frame); P<=2 cap stab_Omega^perp has dim 180 (D is needed).")
two1 = P2 + kernel_in(STAB, [gram_row(u) for u in IDV + QSING])           # P2 + (D_0 cap stab_Omega)
check(same_space(kernel_in(STAB, [gram_row(u) for u in IDV + QSING]), C_c), "D_0 cap stab_Omega = C")
two2 = kernel_in(D0, [gram_row(u) for u in QSING + C_p])                  # D_0 cap (P2 + Lie(G_Omega))
check(same_space(two1, RECS["E'"]) and same_space(two2, RECS["E'"]) and same_space(
    kernel_in(D0 + LIE, [gram_row(u) for u in QSING + C_p]), RECS["E'"]), "R_E' two-clause forms")
print("  two-clause covariant form for E': R_E' = P2 + (D_0 cap stab_Omega) = D_0 cap (P2 + Lie(G_Omega)).")

print("Exchanges on P<=2: G_Omega-covariance (f ad_K = ad_K f, 11 generators), involution, algebra type on two letters:")
for k in FR:
    cov = all(fapply(FR[k], ad_vec(K, unit(t))) == ad_vec(K, fapply(FR[k], unit(t))) for K in LIE for t in P2T)
    check(cov == (k != "E"), f"covariance of {k}")
    invo = all(fapply(FR[k], fapply(FR[k], unit(t))) == unit(t) for t in P2T)
    check(invo, f"{k} is an involution")
    auto = anti = True
    for x in LET:
        fx = fapply(FR[k], letter_vec(x))
        (ix, sx), = sparse(fx)
        lx = (occ(P2T[ix])[0], P2T[ix][occ(P2T[ix])[0]])
        for y in LET:
            fy = fapply(FR[k], letter_vec(y))
            (iy, sy), = sparse(fy)
            ly = (occ(P2T[iy])[0], P2T[iy][occ(P2T[iy])[0]])
            lhs = fapply(FR[k], mprod(x, y))
            auto = auto and lhs == scal(sx * sy, mprod(lx, ly))
            anti = anti and lhs == scal(sx * sy, mprod(ly, lx))
    kind_ = "automorphism" if auto else ("anti-automorphism" if anti else "neither")
    check(kind_ == {"E": "automorphism", "E'": "anti-automorphism"}.get(k, "neither"), f"algebra type of {k}")
    print(f"  {FRAME_LABEL[k]:24s}: covariant {cov}; involution {invo}; {kind_}")
for d, s1, s2, gens in ((2, 1, 3, LIE2), (3, 0, 2, LIE3)):
    n = d * d - 1
    src = [tuple(a if s == s1 else 0 for s in range(4)) for a in range(1, d * d)]
    dst = [tuple(a if s == s2 else 0 for s in range(4)) for a in range(1, d * d)]
    eqs = []
    for K in gens:
        As = [[ad_vec(K, unit(src[b]))[IDX[src[a]]] for b in range(n)] for a in range(n)]
        Ad = [[ad_vec(K, unit(dst[b]))[IDX[dst[a]]] for b in range(n)] for a in range(n)]
        for a in range(n):
            for b in range(n):
                row = [0] * (n * n)
                for c in range(n):
                    row[a * n + c] += As[c][b]
                    row[c * n + b] -= Ad[a][c]
                eqs.append(row)
    check(n * n - rank(eqs) == 1, f"Schur: intertwiners V_slot1 -> V_slot2 at d = {d} form a line")
print("  Schur: the G_Omega-intertwiners V_q1 -> V_q2 and V_Q1 -> V_Q2 are lines (b -> b^T), so the covariant")
print("  letter-wise exchanges are (eps_Q)^NQ (eps_q)^Nq E'; only E' is an (anti-)automorphism of B(H).")
# the readers remark (alpha.tex:79, supplement.tex:170): at a qubit A^T = -sigma_y A sigma_y, so Fm acts on the read
# qubits as conjugation by the unitary W = (sigma_y (x) sigma_y) SWAP_q; sigma_y = -i J with J = A01 (real)
J2 = LB[2][0][2]
for b in LB[2][0][1:]:
    check((b.T == J2 @ b @ J2).all(), "b^T = J b J, i.e. -b^T = sigma_y b sigma_y")
SWq = np.zeros((36, 36), dtype=np.int64)
for idx in itertools.product(range(3), range(2), range(3), range(2)):
    SWq[np.ravel_multi_index((idx[0], idx[3], idx[2], idx[1]), DIMS), np.ravel_multi_index(idx, DIMS)] = 1
Wq = np.kron(np.kron(np.eye(3, dtype=np.int64), J2), np.kron(np.eye(3, dtype=np.int64), J2)) @ SWq
nqt = 0
for i, t in enumerate(P2T):
    if all(t[s] == 0 for s in QS):
        img = fapply(FR["Fm"], unit(t))
        check((vec_to_op(img).astype(np.int64) == Wq @ OPM[i] @ Wq.T).all(), "Fm = Ad_W on the read qubits")
        nqt += 1
check((Wq @ Wq.T == np.eye(36, dtype=np.int64)).all(), "W unitary")
print(f"  readers remark: -b^T = sigma_y b sigma_y on qubit letters, and Fm = Ad_W, W = (sigma_y (x) sigma_y) SWAP_q, on all")
print(f"  {nqt} basis operators of P<=2 without qutrit letters: Fm identifies the read qubits by the matching that undoes")
print("  the conjugation; at the qutrit copies the only intertwiner is the transpose (no unitary does it).")
# operator form of each frame's record space with the occupation projector [N = 1]
for k in FR:
    rows = []
    for j in range(NP):
        r = [0] * NP
        if len(occ(P2T[j])) == 1:
            r[j] += 1
        fj = IDX[fimage(P2T[j])]
        if len(occ(P2T[fj])) == 1:
            r[fj] += fsign(FR[k], P2T[fj])
        if any(r):
            rows.append(r)
    check(same_space(kernel_in(D0, rows), RECS[k]), f"R_f = D_0 cap ker((1+f)[N=1]) ({k})")
print("  operator form, every frame: R_f = D_0 cap ker((1 + f)[N = 1]); the occupation projector [N = 1] is the clause.")

# ==========================================================================================================
# 6. collapse
# ==========================================================================================================
head("6. Collapse X -> X|Omega> under the transpose identity (Prop. 8)")


def kap_of(vs):
    return [[sum(x * KAP[i][r] for i, x in sparse(v)) for r in range(36)] for v in vs]


check(same_space(kernel_in(IDEN, KAPROWS), LIE), "single-letter operators annihilating |Omega> = Lie(G_Omega)")
check(all(not any(r) for r in kap_of(C_c)), "kappa(C) = 0")
kP2, kCp = kap_of(P2), kap_of(C_p)
check(rank(kP2 + kCp) == rank(kP2), "kappa(C+) inside kappa(P2)")
check(same_space(kap_of(CELLS["E'"][("qq", -1)]), kCp) and rank(kCp) == 3, "kappa(C+) = kappa(E'-odd qq products)")
check(rank(kap_of(C_L)) == 1 and rank(kap_of(C_Lp)) == 2, "plain frame: dims 1 and 2")
check(same_space(P2 + kernel_in(D0, KAPROWS), D0), "P2 + (D_0 cap ker kappa) = D_0")
print(f"rank kappa on P<=2: {rank(kap_of(IDV + P2 + QSING + qSING))} (of 36);  ker kappa cap P<=1 = Lie(G_Omega) (11);  kappa(C) = 0;")
print("kappa(C+) = kappa(E'-odd qq products), dim 3, inside kappa(P2);  plain frame: dim kappa(C_L) = 1, dim kappa(C_L+) = 2;")
print("P2 + (D_0 cap ker kappa) = D_0.")

# ==========================================================================================================
# 7. the declared family
# ==========================================================================================================
head("7. The declared family F (Prop. 9)")
BASES = {
    "P2": PRODS,
    "P<=2": CLASSES,
    "m(Sym2)": ["id", "Q"] + PRODS,
    "m(Lam2)": ["Q", "q"],
    "D": ["id", "q"] + PRODS,
    "D+": ["id"] + PRODS,
    "D-": ["q"],
    "1/read": ["QQ", "qq", "mC"],
}
check(same_space(cellspace("E'", [ck for ck in CELLKEYS if ck[0] in BASES["m(Sym2)"]]), MSYM), "base m(Sym2)")
check(same_space(cellspace("E'", [ck for ck in CELLKEYS if ck[0] in BASES["m(Lam2)"]]), MANT), "base m(Lam2)")


def gram_of_comm(Pm):
    """rows of C^T C for C: X -> [X, P] (exact integers); its kernel on P<=2 is the commutant of P"""
    Cm = np.array([(OPM[i] @ Pm - Pm @ OPM[i]).reshape(-1) for i in range(NP)], dtype=np.int64)
    G_ = Cm @ Cm.T
    check(np.abs(Cm).max() < 2 ** 20 and np.abs(G_).max() < 2 ** 62, "no int64 overflow")
    return [[int(x) for x in row] for row in G_]


P3 = np.zeros((36, 36), dtype=np.int64)
P2q = np.zeros((36, 36), dtype=np.int64)
for i1, j1, i2, j2 in itertools.product(range(3), range(2), range(3), range(2)):
    if i1 == i2:
        for k in range(3):
            P3[np.ravel_multi_index((i1, j1, i2, j2), DIMS), np.ravel_multi_index((k, j1, k, j2), DIMS)] = 1
    if j1 == j2:
        for k in range(2):
            P2q[np.ravel_multi_index((i1, j1, i2, j2), DIMS), np.ravel_multi_index((i1, k, i2, k), DIMS)] = 1
PL = []
for link in (0, 1):
    M = np.zeros((36, 36), dtype=np.int64)
    for idx in itertools.product(range(3), range(2), range(3), range(2)):
        a, b = (idx[0], idx[1]) if link == 0 else (idx[2], idx[3])
        if a == b and a < 2:
            for k in range(2):
                jdx = list(idx)
                jdx[2 * link], jdx[2 * link + 1] = k, k
                M[np.ravel_multi_index(idx, DIMS), np.ravel_multi_index(tuple(jdx), DIMS)] = 1
    PL.append(M)
POM = np.outer(PHI, PHI)
check(((P3 @ P3) == 3 * P3).all() and ((P2q @ P2q) == 2 * P2q).all() and all(((m @ m) == 2 * m).all() for m in PL),
      "scaled projectors (3 P_3, 2 P_2, 2 P_Li)")
check(((P3 @ P2q) == POM).all(), "6 P_Omega = (3 P_3)(2 P_2)")
check((PL[0] @ LKET == 2 * LKET).all() and (PL[1] @ LKET == 2 * LKET).all(), "read-link kets")
TMAP = {1: (1, 1), 2: (4, -1), 3: (7, 1)}     # transport t: S01 -> S01, A01 -> -A01, D1 -> D1 (d = 3 indices)
RLM = []
for sQ, sq in ((0, 1), (2, 3)):
    for b in (1, 2, 3):
        a, sg = TMAP[b]
        v = unit(tuple(a if s == sQ else 0 for s in range(4)), sg)
        v[IDX[tuple(b if s == sq else 0 for s in range(4))]] -= 1
        check(not (vec_to_op(v).astype(np.int64) @ LKET).any(), "read-link matching annihilates |L1 L2>")
        RLM.append(v)
ADROWS = []
for K in LIE:
    cols = ad_cols(K)
    rows = {}
    for i, nz in cols.items():
        for j, c in nz:
            rows.setdefault(j, [0] * NP)[i] += c
    ADROWS += list(rows.values())
INVG = kernel_in(IDV + P2 + QSING + qSING, ADROWS)
check(rank(INVG) == 3, "three G_Omega-invariants in P<=2")
KAPL = [[int(x) for x in (OPM[i] @ LKET)] for i in range(NP)]
EIGT = [[int((OPM[i].T @ PHI)[r]) * NPHI - EXPV[i] * int(PHI[r]) for i in range(NP)] for r in range(36)]
COND = {
    "none": ("cell", lambda k: set(CELLKEYS)),
    "traceless": ("cell", lambda k: set(CELLKEYS) - {("id", 1)}),
    "f-even": ("cell", lambda k: {ck for ck in CELLKEYS if ck[1] == 1}),
    "f-odd": ("cell", lambda k: {ck for ck in CELLKEYS if ck[1] == -1}),
    "own sign +": ("cell", lambda k: {ck for ck in CELLKEYS if SIG[ck[0]] == 1}),
    "own sign -": ("cell", lambda k: {ck for ck in CELLKEYS if SIG[ck[0]] == -1}),
    "survivor": ("cell", lambda k: {ck for ck in CELLKEYS if SIG[ck[0]] * ck[1] == 1}),
    "non-survivor": ("cell", lambda k: {ck for ck in CELLKEYS if SIG[ck[0]] * ck[1] == -1}),
    "perp I": ("cell", lambda k: {ck for ck in CELLKEYS if ck[0] in PRODS}),
    "perp I^f": ("cell", lambda k: {ck for ck in CELLKEYS if not (ck[0] in ("id", "Q", "q") and ck[1] == 1)}),
    "perp Comp_f": ("rows", lambda k: [gram_row(y) for y in COMP[k]]),
    "perp stab": ("rows", lambda k: [gram_row(y) for y in STAB]),
    "perp Lie": ("rows", lambda k: [gram_row(y) for y in LIE]),
    "perp G-inv": ("rows", lambda k: [gram_row(y) for y in INVG]),
    "G-invariant": ("rows", lambda k: ADROWS),
    "perp RL-match": ("rows", lambda k: [gram_row(y) for y in RLM]),
    "perp all-match": ("rows", lambda k: [gram_row(y) for y in LIE + RLM]),
    "ker kappa": ("rows", lambda k: KAPROWS),
    "<Om|X|Om>=0": ("rows", lambda k: [EXPV]),
    "X Om in C Om": ("rows", lambda k: EIGROWS),
    "[X,P_Om]=0": ("rows", lambda k: EIGROWS + EIGT),
    "ker kappa_L": ("rows", lambda k: [[KAPL[i][r] for i in range(NP)] for r in range(36)]),
    "[X,P_Li]=0": ("rows", lambda k: gram_of_comm(PL[0]) + gram_of_comm(PL[1])),
    "[X,P_3]=[X,P_2]=0": ("rows", lambda k: gram_of_comm(P3) + gram_of_comm(P2q)),
}
FRAMEDEP = {"f-even", "f-odd", "survivor", "non-survivor", "perp I^f", "perp Comp_f"}
FAM = ("E", "E'", "Fm")
ROWS = {}
for cname, (ctype, fn) in COND.items():
    if ctype == "rows":
        for k in (FAM if cname in FRAMEDEP else ("E'",)):
            ROWS[(cname, k)] = basis_of(fn(k)) if len(fn(k)) > 60 else fn(k)
CELLRK = {}
results = []
for bname, bcls in BASES.items():
    for tl in (False, True):
        if tl and "id" not in bcls:
            continue
        for cname, (ctype, fn) in COND.items():
            for k in FAM:
                bcells = [ck for ck in CELLKEYS if ck[0] in bcls and not (tl and ck[0] == "id")]
                if ctype == "cell":
                    keep = [ck for ck in bcells if ck in fn(k)]
                    cont = {ck: (CDIM[ck] if ck in keep else 0) for ck in CELLKEYS}
                    dW, gr = sum(cont.values()), True
                else:
                    rk = (cname, k if cname in FRAMEDEP else "E'")
                    R = ROWS[rk]
                    cont = {}
                    for ck in CELLKEYS:
                        if ck in bcells:
                            key = (rk, k, ck)
                            if key not in CELLRK:
                                CELLRK[key] = CDIM[ck] - rank_rows_on(R, CELLS[k][ck])
                            cont[ck] = CELLRK[key]
                        else:
                            cont[ck] = 0
                    Bv = cellspace(k, bcells)
                    dW = sum(CDIM[ck] for ck in bcells) - rank_rows_on(R, Bv)
                    gr = sum(cont.values()) == dW
                lc, rest = letter_classes(cont)
                isR = {}
                if dW == 172:
                    W = cellspace(k, keep) if ctype == "cell" else kernel_in(cellspace(k, bcells), ROWS[rk])
                    isR = {kk: same_space(W, RECS[kk]) for kk in FAM}
                    if not gr:
                        cont_, gr_ = content(k, W)
                        check(not gr_, "gradedness consistent")
                results.append(((bname, tl, cname, k), dW, gr, lc, rest,
                                numerator(cont, w1) if gr else None, numerator(cont, w2) if gr else None, isR))
nm_ = len(results)
print(f"F: {len(BASES)} composite spaces (with and without the identity where it occurs) x {len(COND)} conditions x"
      f" 3 frames (E, E', Fm) = {nm_} members")
print("  composite spaces: P2 (two letters at distinct slots), P<=2 = m(V(x)V), m(Sym2 V), m(Lambda2 V), D, D_+, D_-,")
print("  one letter per read.  Conditions: " + "; ".join(COND))
dimset = Counter(r[1] for r in results)
print(f"  dimensions attained: {sorted(dimset)}")
h172 = [r for r in results if r[1] == 172]
print(f"  members of dimension 172: {len(h172)}")
for (bname, tl, cname, k), dW, gr, lc, rest, n1, n2, isR in h172:
    eq = [kk for kk, v in isR.items() if v]
    print(f"    [{k:3s}] {bname}{'_0' if tl else ''} cap {cname:12s}: graded {str(gr):5s}; classes {lc} + {rest} outside;"
          f" equals R_f for f in {eq}; 1/value {fmt_inv(n1)} (M_L), {fmt_inv(n2)} (word)")
hitsR = {k: sorted({(r[0][0] + ("_0" if r[0][1] else ""), r[0][2]) for r in h172 if r[0][3] == k and r[7].get(k)})
         for k in FAM}
for k in FAM:
    print(f"  members equal to the frame's own record space R_f, f = {k}: {hitsR[k]}")
check({c for _, c in hitsR["E"]} == {"perp I^f", "perp Comp_f"}, "E: the two identification rules")
check({c for _, c in hitsR["E'"]} == {"perp I^f"}, "E': only perp I^E' (not a Lie algebra)")
check({c for _, c in hitsR["Fm"]} == {"perp I^f", "perp Comp_f", "perp stab", "perp Lie"}, "Fm: four rules")
check(all(b.startswith("D") for k in FAM for b, _ in hitsR[k]), "all on D")
okclasses = [r for r in results if r[3] == LETTER and r[4] == 0 and r[2]]
check(all(r[7].get(r[0][3]) for r in okclasses), "every member with the Letter's classes is the frame's R_f")
check(not any(r[1] == 172 and r[0][0] not in ("D",) for r in results if r[7] and any(r[7].values())), "only D")

# ==========================================================================================================
# 8. one-step images
# ==========================================================================================================
head("8. Single-letter content of Phi(P2) for natural superoperators Phi (Prop. 10)")
P1T = [t for t in P2T if len(occ(t)) <= 1]
P1N = [GR[IDX[t]] for t in P1T]
P1M = [OPM[IDX[t]] for t in P1T]
PRIMS = {"join P3": P3, "qubit pair P2": P2q, "joint state P_Om": POM, "read link P_L1": PL[0],
         "read link P_L2": PL[1], "read links P_L1 P_L2": PL[0] @ PL[1]}
p2idx = [i for i in range(NP) if cls(P2T[i]) in PRODS]
anyhit = False
print("Phi(X) for the 169 basis products X, projected onto P<=1 (identity and single letters); dims per cell (E'):")
for pn, Pm in PRIMS.items():
    for side in ("L", "R", "Ad"):
        vs = []
        for i in p2idx:
            X = OPM[i]
            Y = Pm @ X if side == "L" else (X @ Pm if side == "R" else Pm @ X @ Pm)
            v = [0] * NP
            for t, Bm, nb in zip(P1T, P1M, P1N):
                v[IDX[t]] = Fr(int(np.sum(Bm * Y)), nb)
            vs.append(v)
        W = basis_of(vs)
        cont, gr = content("E'", W)
        hit = bool(W) and (same_space(W, C_c) or same_space(W, C_L) or same_space(W, C_p))
        anyhit = anyhit or hit
        cs = ", ".join(f"{ckname(ck)}:{n}" for ck, n in cont.items() if n)
        print(f"  {side + '_' + pn:24s}: dim {rank(W):2d}; graded {str(gr):5s}; in cells: {cs if cs else '-'}")
check(not anyhit, "no Phi(P2) has single-letter part C, C_L or C+")
print("  none of these 18 has single-letter part equal to C, C_L or C+.")

# ==========================================================================================================
# 9. variants
# ==========================================================================================================
head("9. Record-space variants and their values (weights: M_L of the stated frame, and the word rule)")
GML = []                       # lambda_1, lambda_2, lambda_3, lambda_8 differences: S01, A01, D1, D2 at Q1 minus Q2
for a in (1, 4, 7, 8):
    v = unit((a, 0, 0, 0))
    v[IDX[(0, 0, a, 0)]] = -1
    GML.append(v)
check(contains(CELLS["E"][("Q", -1)], GML) and rank(GML) == 4, "Gell-Mann literal records: 4 E-odd qutrit")
VARS = [
    ("the Letter: P2 + C_L, E", "E", RECS["E"]),
    ("F-CT-10: P2 + C, E'", "E'", RECS["E'"]),
    ("D cap stab^perp = P2 + C+, Fm", "Fm", RECS["Fm"]),
    ("D cap stab^perp = P2 + C+, E'", "E'", P2 + C_p),
    ("P2 only (no commutator records)", "E'", P2),
    ("D_0 = P2 + all six single-qubit (both kept)", "E'", D0),
    ("P2 + Lie(G_Omega) (all matchings)", "E'", P2 + LIE),
    ("P<=2,0 (all traceless two-letter composites)", "E'", P2 + QSING + qSING),
    ("survivors only: even products + C", "E'", cellspace("E'", [(c, 1) for c in PRODS] + [("q", -1)])),
    ("Gell-Mann literal: R_E + 4 qutrit differences", "E", RECS["E"] + GML),
]
for nm, k, W in VARS:
    cont, gr = content(k, W)
    check(gr, f"variant graded: {nm}")
    n1, n2 = numerator(cont, w1), numerator(cont, w2)
    print(f"  {nm:46s} dim {rank(W):3d}  M_L {fmt_num(n1):20s} {fmt_inv(n1):>11s}   word {fmt_num(n2):20s} {fmt_inv(n2):>11s}")
check(numerator(content("E", RECS["E"] + GML)[0], w2) == (79, 48, 46, 3), "Gell-Mann literal")
check(numerator(content("E'", P2 + C_p)[0], w1) == (82, 48, 42, 0), "C+ under E', M_L: Table S5 'no own sign'")
check(numerator(content("E'", P2 + C_p)[0], w2) == (79, 48, 45, 0), "C+ under E', word: 'survive X only'")
check(numerator(content("E'", P2)[0], w1) == (79, 48, 42, 0), "no commutator records")
S5_both, own_plus = (79, 48, 42, 6), (79, 51, 42, 3)
check(abs(inv_value(S5_both) - mp.mpf("132.738534414")) < mp.mpf("1e-8"), "Table S5 'sums and differences kept'")
print(f"  Table S5 'sums and differences kept' gives the sums r6: {fmt_num(S5_both)}, 1/value {fmt_inv(S5_both)}; the Letter's")
print(f"  survival rule gives them r3 at own sign -1 (the word rule above, {fmt_inv((79, 48, 45, 3))}) or r2 at own sign +1")
print(f"  ({fmt_num(own_plus)}, 1/value {fmt_inv(own_plus)}).")
check(abs(inv_value(LETTER) - mp.mpf("137.0359991895661550291837440208667995664")) < mp.mpf(10) ** -30, "value")

print()
print(f"{NCHECK} checks passed.")
print(f"runtime: {time.perf_counter() - T0:.1f} s")
