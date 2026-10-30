"""Independent operator-level reconstruction of the 172 records of the
Chain Theory alpha Letter (Wissner-Gross 2026, v1.0.0, Sec. S7).

Written from the Letter's prose only; does not import or copy the author's
enumerate_records.py or ChainAlpha.lean.

Hilbert space of the four slots, in the order Q1, q1, Q2, q2:
    H = C^3 (x) C^2 (x) C^3 (x) C^2,  dim 36.
Records are 36x36 operators. The exchange of the two reads (= exchange of the
two read links) is conjugation by the unitary S that swaps the factor pairs
(Q1,q1) <-> (Q2,q2).

Checks performed (all must pass, else AssertionError):
  1. 169 two-slot product operators; slot-pair sizes 24,24,64,24,24,9.
  2. S permutes the products; it fixes 11 and pairs 158 -> 79 odd + 90 even.
  3. three commutator records [s_mu,s_nu]_{q1} - [s_mu,s_nu]_{q2}, odd under S.
  4. the 172 records are pairwise Hilbert-Schmidt orthogonal (rank 172).
  5. sector sizes 79/48/42/3 under the Letter's weight rule.
  6. closed form (79 + 48 r2 + 42 r3 + 3 r6)/31104 and 1/alpha to 60 digits.

Run:  python3 records_independent.py
"""
import itertools
import numpy as np
import sympy as sp
import mpmath as mp

# ---------- letters ----------
def gell_mann():
    L = []
    for j in range(3):
        for k in range(j + 1, 3):
            m = np.zeros((3, 3), complex); m[j, k] = m[k, j] = 1; L.append(m)       # symmetric
            m = np.zeros((3, 3), complex); m[j, k] = -1j; m[k, j] = 1j; L.append(m)  # antisymmetric
    L.append(np.diag([1, -1, 0]).astype(complex))
    L.append(np.diag([1, 1, -2]).astype(complex) / np.sqrt(3))
    return L

def pauli():
    return [np.array([[0, 1], [1, 0]], complex),
            np.array([[0, -1j], [1j, 0]], complex),
            np.array([[1, 0], [0, -1]], complex)]

GM, PA = gell_mann(), pauli()
assert len(GM) == 8 and len(PA) == 3
for B in (GM, PA):
    for A in B:
        assert abs(np.trace(A)) < 1e-12 and np.allclose(A, A.conj().T)
    G = np.array([[np.trace(A @ C).real for C in B] for A in B])
    assert np.allclose(G, 2 * np.eye(len(B)))          # orthogonal, tr(AB) = 2 delta

DIMS = (3, 2, 3, 2)                 # Q1, q1, Q2, q2
LETTERS = (GM, PA, GM, PA)
END = ('Q', 'q', 'Q', 'q')          # which end (qutrit / qubit) each slot is
READ = (1, 1, 2, 2)                 # which read link each slot belongs to

def kron_all(ms):
    out = np.array([[1.0 + 0j]])
    for m in ms:
        out = np.kron(out, m)
    return out

def slot_op(assign):
    """assign: dict slot -> 2x2/3x3 matrix; identity elsewhere."""
    return kron_all([assign.get(s, np.eye(DIMS[s])) for s in range(4)])

# ---------- the exchange S: (Q1,q1,Q2,q2) -> (Q2,q2,Q1,q1) ----------
def swap_unitary():
    n = 36
    S = np.zeros((n, n))
    for idx in itertools.product(*[range(d) for d in DIMS]):
        a, b, c, d = idx
        src = np.ravel_multi_index(idx, DIMS)
        dst = np.ravel_multi_index((c, d, a, b), DIMS)
        S[dst, src] = 1
    return S

S = swap_unitary()
assert np.allclose(S @ S.T, np.eye(36))
ex = lambda A: S @ A @ S.T

# ---------- 1. products ----------
products = []   # (label, operator, slots)
pair_sizes = {}
for s, t in itertools.combinations(range(4), 2):
    cnt = 0
    for i, A in enumerate(LETTERS[s]):
        for j, B in enumerate(LETTERS[t]):
            products.append(((s, i, t, j), slot_op({s: A, t: B}), (s, t)))
            cnt += 1
    pair_sizes[(s, t)] = cnt
assert len(products) == 169
assert sorted(pair_sizes.values()) == sorted([24, 24, 64, 24, 24, 9])

# ---------- 2. S permutes products ----------
key = {lab: k for k, (lab, _, _) in enumerate(products)}
perm = {}
for k, (lab, A, _) in enumerate(products):
    B = ex(A)
    hits = [m for m, (_, C, _) in enumerate(products) if np.allclose(B, C)]
    assert len(hits) == 1
    perm[k] = hits[0]
fixed = [k for k in perm if perm[k] == k]
assert len(fixed) == 11
orbits = sorted({tuple(sorted((k, perm[k]))) for k in perm})
two_orbits = [o for o in orbits if o[0] != o[1]]
assert len(two_orbits) == 79

def same_end(slots):
    return END[slots[0]] == END[slots[1]]

# eigen-combinations of S
records = []   # (kind, operator, E-eigenvalue, category)
for o in orbits:
    A = products[o[0]][1]
    slots = products[o[0]][2]
    cat = 'same' if same_end(slots) else 'mixed'
    if o[0] == o[1]:
        records.append(('fixed', A, +1, cat))
    else:
        B = products[o[1]][1]
        records.append(('even', A + B, +1, cat))
        records.append(('odd', A - B, -1, cat))
for kind, A, e, _ in records:
    assert np.allclose(ex(A), e * A)

# ---------- 3. commutator records ----------
comm = []
for mu, nu in itertools.combinations(range(3), 2):
    C = PA[mu] @ PA[nu] - PA[nu] @ PA[mu]
    D = slot_op({1: C}) - slot_op({3: C})
    assert np.allclose(ex(D), -D)
    comm.append(D)
    records.append(('comm', D, -1, 'comm'))
assert len(records) == 172

# ---------- 4. orthogonality ----------
M = np.array([r[1].reshape(-1) for r in records])
G = M.conj() @ M.T
off = G - np.diag(np.diag(G))
assert np.max(np.abs(off)) < 1e-9
assert np.linalg.matrix_rank(M) == 172

# ---------- 5. sectors under the Letter's rule ----------
# Letter, Eq. (S-exchanges): X acts on same-end products and commutators,
# Y on mixed products and commutators; both act by the SAME involution S on
# products. A product "survives" (reads -1) when S returns it times +1; a
# commutator record survives when S returns it times its own sign -1.
def XY(rec):
    kind, _, e, cat = rec
    if cat == 'comm':
        return (-1, -1)
    survive = (e == +1)
    if cat == 'same':
        return (-1 if survive else +1, +1)
    return (+1, -1 if survive else +1)

from collections import Counter
sectors = Counter(XY(r) for r in records)
assert sectors == Counter({(1, 1): 79, (1, -1): 48, (-1, 1): 42, (-1, -1): 3}), sectors

# ---------- 6. value ----------
r2, r3, r6 = sp.sqrt(2), sp.sqrt(3), sp.sqrt(6)
def k_of(xy):
    X, Y = xy
    return sp.Integer(3) ** int(X == -1) * sp.Integer(2) ** int(Y == -1)
num = sum(sp.sqrt(k_of(XY(r))) for r in records)
N = 24 * 36 * 36
alpha = sp.nsimplify(num) / N
target = (79 + 48 * r2 + 42 * r3 + 3 * r6) / 31104
assert sp.simplify(alpha - target) == 0
mp.mp.dps = 80
inv = mp.mpf(31104) / (79 + 48 * mp.sqrt(2) + 42 * mp.sqrt(3) + 3 * mp.sqrt(6))

print("products:", len(products), "slot-pair sizes:", [pair_sizes[p] for p in sorted(pair_sizes)])
print("S-fixed products:", len(fixed), " exchange pairs:", len(two_orbits))
print("records:", len(records), " Gram matrix diagonal, rank", np.linalg.matrix_rank(M))
print("sectors (X,Y):count", dict(sectors))
print("numerator:", sp.nsimplify(num), " N =", N)
print("alpha_cand =", mp.nstr(1 / inv, 40))
print("1/alpha_cand =", mp.nstr(inv, 60))
print("ALL CHECKS PASSED")
