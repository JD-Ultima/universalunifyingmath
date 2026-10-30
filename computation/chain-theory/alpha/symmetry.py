"""Symmetry analysis of the record space (audit F-CT-10).

G_Omega: the symmetry group of the joint read state
    |Omega> = |omega_3>_{Q1Q2} (x) |omega_2>_{q1q2}
(verified in weight_functionals.py), i.e. U (x) Ubar on (Q1,Q2) and V (x) Vbar on (q1,q2),
U in SU(3), V in SU(2), acting on operators by conjugation.

Questions:
 1. Does the Letter's slot exchange E (plain swap Q1<->Q2, q1<->q2) commute with G_Omega?
    Does E' = E composed with the transpose at every slot?
 2. Is the span of the Letter's 172 records G_Omega-invariant?  Is the span of the
    'covariant' records (built with E' instead of E) invariant?
 3. Sector counts and value with E' (do they change?).
 4. Isotypic decomposition of the covariant record space under SU(3) x SU(2) x <E'>,
    dimension of the commutant, and the number of free parameters in tr(M) for an
    invariant dressing M.
Run: python3 symmetry.py
"""
import importlib.util, io, contextlib, itertools
import numpy as np
from collections import Counter
spec = importlib.util.spec_from_file_location("rec", __file__.replace("symmetry.py", "records_independent.py"))
rec = importlib.util.module_from_spec(spec)
with contextlib.redirect_stdout(io.StringIO()):
    spec.loader.exec_module(rec)
rng = np.random.default_rng(7)
DIMS, GM, PA, slot_op, kron_all = rec.DIMS, rec.GM, rec.PA, rec.slot_op, rec.kron_all

def haar_su(n):
    Z = (rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))) / np.sqrt(2)
    Q, R = np.linalg.qr(Z); Q = Q @ np.diag(np.diag(R) / abs(np.diag(R)))
    return Q / np.linalg.det(Q) ** (1 / n)

def g_unitary(U, V):
    return kron_all([U, V, U.conj(), V.conj()])

# --- transpose at every slot (partial transposes on all four factors = full transpose) ---
def T_all(A):
    return A.T
S = rec.S
E = lambda A: S @ A @ S.T
Ep = lambda A: T_all(E(A))

# 1. commutation with G_Omega
U, V = haar_su(3), haar_su(2)
W = g_unitary(U, V)
g = lambda A: W @ A @ W.conj().T
A = rng.normal(size=(36, 36)) + 1j * rng.normal(size=(36, 36))
print("1. ||E g - g E|| / ||A|| =", f"{np.linalg.norm(E(g(A)) - g(E(A))) / np.linalg.norm(A):.3e}",
      "   ||E' g - g E'|| / ||A|| =", f"{np.linalg.norm(Ep(g(A)) - g(Ep(A))) / np.linalg.norm(A):.3e}")

def span_basis(ops):
    M = np.array([o.reshape(-1) for o in ops]).T
    Q, R = np.linalg.qr(M)
    r = np.linalg.matrix_rank(M, tol=1e-9)
    return Q[:, :r]

def invariance_defect(Bq, gfun):
    defect = 0.0
    for k in range(Bq.shape[1]):
        o = Bq[:, k].reshape(36, 36)
        v = gfun(o).reshape(-1)
        defect = max(defect, np.linalg.norm(v - Bq @ (Bq.conj().T @ v)) / np.linalg.norm(v))
    return defect

# 2a. Letter's records
Bq_letter = span_basis([r[1] for r in rec.records])
print("2a. Letter's 172-record span: dim", Bq_letter.shape[1],
      " max relative leakage under a random g in G_Omega:", f"{invariance_defect(Bq_letter, g):.3e}")
# which part leaks: the commutator records
Bq_comm = span_basis(rec.comm)
Bq_prod = span_basis([p[1] for p in rec.products])
print("    product span (169): leakage", f"{invariance_defect(Bq_prod, g):.3e}",
      "; E-even product span leakage", f"{invariance_defect(span_basis([r[1] for r in rec.records if r[0] in ('even','fixed')]), g):.3e}",
      "; commutator span (3): leakage", f"{invariance_defect(Bq_comm, g):.3e}")

# 3. covariant records with E'
prods = [p[1] for p in rec.products]
key = lambda A: tuple(np.round(A.reshape(-1), 9))
covariant = []
used = set()
for i, A in enumerate(prods):
    if i in used: continue
    B = Ep(A)
    j = next(k for k, C in enumerate(prods) if np.allclose(C, B) or np.allclose(C, -B))
    sign = 1 if np.allclose(prods[j], B) else -1
    used |= {i, j}
    slots = rec.products[i][2]
    cat = 'same' if rec.END[slots[0]] == rec.END[slots[1]] else 'mixed'
    if j == i:
        covariant.append(('fixed' if sign == 1 else 'antifixed', A, sign, cat))
    else:
        covariant.append(('even', A + sign * prods[j], +1, cat))
        covariant.append(('odd', A - sign * prods[j], -1, cat))
for kind, A_, e, cat in covariant:
    assert np.allclose(Ep(A_), e * A_)
# commutator-type records: single-letter commutators at q1 and q2, E'-odd combinations
singles = []
for mu, nu in itertools.combinations(range(3), 2):
    Cm = PA[mu] @ PA[nu] - PA[nu] @ PA[mu]
    singles.append((slot_op({1: Cm}), slot_op({3: Cm})))
comm_cov = []
for c1, c3 in singles:
    for s in (+1, -1):
        D = c1 + s * c3
        if np.allclose(Ep(D), -D):
            comm_cov.append(D)
assert len(comm_cov) == 3
cov_ops = [c[1] for c in covariant] + comm_cov
Bq_cov = span_basis(cov_ops)
print("2b. covariant record span: dim", Bq_cov.shape[1], " leakage under g:", f"{invariance_defect(Bq_cov, g):.3e}",
      " under E':", f"{invariance_defect(Bq_cov, Ep):.3e}")

def sector(kind, e, cat):
    if cat == 'comm': return (-1, -1)
    surv = (e == +1)
    return ((-1 if surv else 1), 1) if cat == 'same' else (1, (-1 if surv else 1))
sec = Counter([sector(k, e, c) for k, _, e, c in covariant] + [(-1, -1)] * 3)
print("3. sectors with E':", dict(sec), "-> same as Letter:",
      sec == Counter({(1, 1): 79, (1, -1): 48, (-1, 1): 42, (-1, -1): 3}))
print("   E'-antifixed products (E' p = -p):", sum(1 for c in covariant if c[0] == 'antifixed'))

# 4. isotypic decomposition on the covariant span
def ad_matrix(K, Bq):
    cols = []
    for k in range(Bq.shape[1]):
        o = Bq[:, k].reshape(36, 36)
        cols.append(Bq.conj().T @ (K @ o - o @ K).reshape(-1))
    return np.array(cols).T
T3 = [l / 2 for l in GM]; T2 = [p / 2 for p in PA]
K3 = [slot_op({0: t}) - slot_op({2: t.T}) for t in T3]   # generator of U (x) Ubar on (Q1,Q2)
K2 = [slot_op({1: t}) - slot_op({3: t.T}) for t in T2]   # generator of V (x) Vbar on (q1,q2)
C3 = sum(ad_matrix(K, Bq_cov) @ ad_matrix(K, Bq_cov) for K in K3)
C2 = sum(ad_matrix(K, Bq_cov) @ ad_matrix(K, Bq_cov) for K in K2)
EP = np.array([Bq_cov.conj().T @ Ep(Bq_cov[:, k].reshape(36, 36)).reshape(-1) for k in range(Bq_cov.shape[1])]).T
print("4. [C3,C2], [C3,E'], [C2,E'] norms:",
      [f"{np.linalg.norm(X @ Y - Y @ X):.1e}" for X, Y in [(C3, C2), (C3, EP), (C2, EP)]])
Hmix = C3 + np.pi * C2 + np.e * 100 * EP          # generic combination for joint diagonalization
w, Vv = np.linalg.eigh((Hmix + Hmix.conj().T) / 2)
labels = Counter()
for k in range(len(w)):
    v = Vv[:, k]
    c3 = (v.conj() @ C3 @ v / (v.conj() @ v)).real
    c2 = (v.conj() @ C2 @ v / (v.conj() @ v)).real
    ep = (v.conj() @ EP @ v / (v.conj() @ v)).real
    labels[(round(c3, 4), round(c2, 4), int(round(ep)))] += 1
irrep3 = {0.0: '1', 3.0: '8', 6.0: '10/10bar', 8.0: '27'}
irrep2 = {0.0: '1', 2.0: '3', 6.0: '5'}
dims3 = {'1': 1, '8': 8, '10/10bar': 10, '27': 27}; dims2 = {'1': 1, '3': 3, '5': 5}
rows = []
for (c3, c2, ep), n in sorted(labels.items()):
    r3, r2 = irrep3[c3], irrep2[c2]
    mult = n / (dims3[r3] * dims2[r2])
    rows.append((r3, r2, ep, n, mult))
    print(f"   SU(3) {r3:9s} SU(2) {r2}  E'={ep:+d}: dim {n:3d}, multiplicity {mult:g}")
comm_dim = sum(m * m for *_, m in rows)
print("   isotypic components:", len(rows), "(10/10bar counted together; separate them -> +1)",
      "  commutant dimension (sum of multiplicity^2, 10 and 10bar separate):",
      comm_dim if all(r[0] != '10/10bar' for r in rows) else comm_dim - sum(m*m for r0,_,_,_,m in rows if r0=='10/10bar') + 2*sum((m/2)**2 for r0,_,_,_,m in rows if r0=='10/10bar'))
