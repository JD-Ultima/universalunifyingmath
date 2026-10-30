"""Gap CT-G01 / CT-G12: is there a natural operator origin of the weights sqrt(k_r)?

Bounded search; write-up in formalization/chain-theory/alpha/operator-origin.md.

The Letter's dressed insertion M is DEFINED to carry sqrt(k_r), k in {1,2,3,6}, on
covariant record r (E'-odd products 1; E'-even mixed 2; E'-even same-end 3;
commutator-type 6).  We ask whether an operator built from the model's primitives
alone, without the record list, reproduces these numbers.

Parts (every printed number is computed here):
  0. covariant record basis (172), classes, Letter's M_L = diag sqrt(k), K_L = diag k
  1. primitives: G_Omega-covariance, Q|q factorization, rationality (Galois lemma)
  2. isotypic decomposition, commutant (22), per-isotypic traces of M_L
  3. exact identities: M_L as a polynomial in natural projectors; Schmidt rank of
     the mixed block (factorization obstruction)
  4. family F1: all words of length <= 3 in a 38-letter alphabet of primitives,
     deduplicated up to scale; tests (a) invariance, (b) class-scalar / diagonal
     pattern vs sqrt(k) and vs k (also after twirl, and for |C| and sqrt(C)),
     (c) per-isotypic traces vs the Letter's point
  5. family F2: operator-level "survival x attribution" dressings
     M = prod_x (1 + (sqrt(d_Q)-1) S_x A_x^Q + (sqrt(d_q)-1) S_x A_x^q)
     (7 survival operators x 81 attributions per exchange, two exchanges)

Reuses records_independent.py and symmetry.py (imported, not modified).
Run: python3 operator_origin.py > operator_origin.out     (~16 min, 4 cores)
"""
import importlib.util, io, contextlib, itertools, os, sys, time, hashlib
from fractions import Fraction
from collections import Counter, defaultdict
import numpy as np
import sympy as sp

T_START = time.time()
HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("sym", os.path.join(HERE, "symmetry.py"))
sym = importlib.util.module_from_spec(spec)
with contextlib.redirect_stdout(io.StringIO()):
    spec.loader.exec_module(sym)          # builds the covariant records (and loads records_independent)
rec = sym.rec
DIMS = (3, 2, 3, 2)                       # Q1, q1, Q2, q2
SH = DIMS + DIMS
QS, qS = (0, 2), (1, 3)
NAME = ('Q1', 'q1', 'Q2', 'q2')
TOL = 1e-9
rng = np.random.default_rng(20261030)
nrm = np.linalg.norm


def hdr(s):
    print("\n" + "=" * 78 + "\n" + s + "\n" + "=" * 78)


# ---------------------------------------------------------------- batch helpers
def as8(Y):  return Y.reshape(SH + (Y.shape[-1],))
def as2(T):  return T.reshape(36, 36, T.shape[-1])


def tau_b(s):
    """scalar projection at slot s:  A -> (1_s/d_s) (x) Tr_s A   (partial trace, reinsert)"""
    d = DIMS[s]
    def f(Y):
        t = as8(Y)
        tr = np.trace(t, axis1=s, axis2=4 + s)
        out = np.multiply.outer(tr, np.eye(d) / d)
        return as2(np.ascontiguousarray(np.moveaxis(out, [7, 8], [s, 4 + s])))
    return f


def perm_b(p):
    def f(Y):
        return as2(np.ascontiguousarray(as8(Y).transpose(tuple(p) + (8,))))
    return f


def compose(*fs):
    def f(Y):
        for g in reversed(fs):
            Y = g(Y)
        return Y
    return f


def lin(*terms):
    """linear combination sum c_i f_i"""
    def f(Y):
        return sum(c * g(Y) for c, g in terms)
    return f


ID = lambda Y: Y
def L_b(X):  return lambda Y: np.tensordot(X, Y, axes=(1, 0))
def R_b(X):  return lambda Y: np.ascontiguousarray(np.tensordot(Y, X, axes=(1, 0)).transpose(0, 2, 1))
def Ad_b(X):
    l, r = L_b(X), R_b(X.conj().T)
    return lambda Y: r(l(Y))


def T_slot(s):                            # partial transpose at slot s
    p = list(range(8)); p[s], p[4 + s] = p[4 + s], p[s]; return perm_b(p)


T_full = perm_b((4, 5, 6, 7, 0, 1, 2, 3))
F_Q = perm_b((2, 1, 0, 3, 6, 5, 4, 7))    # swap Q1<->Q2 (both indices)
F_q = perm_b((0, 3, 2, 1, 4, 7, 6, 5))    # swap q1<->q2
E_plain = perm_b((2, 3, 0, 1, 6, 7, 4, 5))  # the Letter's exchange E = Ad_S
Ep = compose(T_full, E_plain)             # E' = T o E   (covariant exchange, F-CT-10)
EpQ = compose(T_slot(0), T_slot(2), F_Q)  # E'_Q: exchange + transpose at the qutrit slots only
Epq = compose(T_slot(1), T_slot(3), F_q)  # E'_q: same at the qubit slots only
TAU = [tau_b(s) for s in range(4)]
COTAU = [lin((1, ID), (-1, TAU[s])) for s in range(4)]
SIGMA = compose(*[lin((2, TAU[s]), (-1, ID)) for s in range(4)])   # (-1)^{N}, N = number of occupied slots


def occ_proj(occ):
    """projector onto operators that are traceless exactly at the slots in occ, scalar elsewhere"""
    return compose(*[(COTAU[s] if s in occ else TAU[s]) for s in range(4)])


def type_proj(nQ, nq):
    occs = [o for r in range(5) for o in itertools.combinations(range(4), r)
            if sum(s in QS for s in o) == nQ and sum(s in qS for s in o) == nq]
    fs = [occ_proj(o) for o in occs]
    return lambda Y: sum(f(Y) for f in fs)


# ---------------------------------------------------------------- pair operators
def ket(idx):
    v = np.zeros(36); v[np.ravel_multi_index(idx, DIMS)] = 1; return v


def pair_projector(s, t, m):
    """|omega_m><omega_m| on slots (s,t), common support = first m levels, identity on the other slots"""
    P = np.zeros((36, 36))
    others = [u for u in range(4) if u not in (s, t)]
    for o in itertools.product(*[range(DIMS[u]) for u in others]):
        for i in range(m):
            for j in range(m):
                a = [0] * 4; b = [0] * 4
                for u, val in zip(others, o): a[u] = b[u] = val
                a[s] = a[t] = i; b[s] = b[t] = j
                P[np.ravel_multi_index(a, DIMS), np.ravel_multi_index(b, DIMS)] += 1.0 / m
    return P


P3 = pair_projector(0, 2, 3)              # the join |omega_3><omega_3| on (Q1,Q2)
P2 = pair_projector(1, 3, 2)              # the qubit pair of the joint state on (q1,q2)
PL1 = pair_projector(0, 1, 2)             # read link 1: |omega_2> on (Q1,q1), support levels 0,1 of Q1
PL2 = pair_projector(2, 3, 2)             # read link 2 on (Q2,q2)
POm = P3 @ P2                             # |Omega><Omega|
for P in (P3, P2, PL1, PL2, POm):
    assert np.allclose(P @ P, P) and np.allclose(P, P.T)
Om = np.zeros(36)
for i in range(3):
    for a in range(2):
        Om[np.ravel_multi_index((i, a, i, a), DIMS)] = 1 / np.sqrt(6)
assert np.allclose(POm, np.outer(Om, Om))


def twirl_pair(slots_pair, d):
    """HS projection onto span{1, P_d} on the pair (= twirl over U (x) Ubar there), identity on the rest"""
    s, t = slots_pair
    u1 = np.eye(d * d).reshape(d, d, d, d) / d                        # (r_s, r_t, c_s, c_t)
    Pd = np.zeros((d, d, d, d))
    for i in range(d):
        for j in range(d):
            Pd[i, i, j, j] = 1 / d
    w = Pd - np.einsum('abcd,abcd', u1, Pd) * u1                      # Gram-Schmidt against u1
    u2 = w / np.sqrt(np.einsum('abcd,abcd', w, w))
    U = np.array([u1, u2])
    other = [x for x in range(4) if x not in (s, t)]
    letters = 'abcdefgh'
    # contraction string: row axes 0..3, col axes 4..7
    def f(Y):
        T = as8(Y)
        rs, rt, cs, ct = letters[s], letters[t], letters[4 + s], letters[4 + t]
        full = letters + 'm'
        keep = ''.join(letters[x] for x in other) + ''.join(letters[4 + x] for x in other) + 'm'
        coef = np.einsum('k' + rs + rt + cs + ct + ',' + full + '->k' + keep, U.conj(), T)
        out = np.einsum('k' + rs + rt + cs + ct + ',k' + keep + '->' + full, U, coef)
        return as2(out)
    return f


TW_Q = twirl_pair((0, 2), 3)
TW_q = twirl_pair((1, 3), 2)
TW = compose(TW_Q, TW_q)                  # G_Omega twirl: projection onto span{1,P3} (x) span{1,P2}

# ---------------------------------------------------------------- 0. records
hdr("0. Covariant record basis (from symmetry.py) and the Letter's weights")
ops, kinds = [], []
for kind, A, e, cat in sym.covariant:
    ops.append(A); kinds.append((kind, e, cat))
for D in sym.comm_cov:
    ops.append(D); kinds.append(('comm', -1, 'comm'))
NR = len(ops)
assert NR == 172
B3 = np.stack([A / nrm(A) for A in ops], axis=-1)        # (36,36,172), orthonormal columns
B = B3.reshape(1296, NR)
Bh = B.conj().T
assert np.allclose(Bh @ B, np.eye(NR), atol=1e-12)

# occupation type of each record from the slot projectors (no record list used for the operators)
TYPE_OPS = {'QQ': type_proj(2, 0), 'qq': type_proj(0, 2), 'mix': type_proj(1, 1), 'comm': type_proj(0, 1)}
rtype = []
for r in range(NR):
    hits = [t for t, f in TYPE_OPS.items() if nrm(f(B3[:, :, r:r + 1])[:, :, 0] - B3[:, :, r]) < 1e-12]
    assert len(hits) == 1
    rtype.append(hits[0])
rtype = np.array(rtype)
epar = np.array([e for _, e, _ in kinds])
cls = np.array(['comm' if t == 'comm' else ('odd' if e == -1 else ('evsame' if t in ('QQ', 'qq') else 'evmix'))
                for t, e in zip(rtype, epar)])
KOF = {'odd': 1, 'evmix': 2, 'evsame': 3, 'comm': 6}
CLASSES = ['odd', 'evmix', 'evsame', 'comm']
k = np.array([KOF[c] for c in cls], float)
print("records:", NR, " class sizes:", {c: int((cls == c).sum()) for c in CLASSES},
      " occupation types:", dict(Counter(rtype)))
assert {c: int((cls == c).sum()) for c in CLASSES} == {'odd': 79, 'evmix': 48, 'evsame': 42, 'comm': 3}
ML = np.diag(np.sqrt(k)).astype(complex)
KL = np.diag(k).astype(complex)
CIDX = {c: np.flatnonzero(cls == c) for c in CLASSES}
PCL = {c: np.diag((cls == c).astype(float)) for c in CLASSES}


def compress(f):
    Y = f(B3)
    Yf = Y.reshape(1296, NR)
    C = Bh @ Yf
    return C, nrm(Yf - B @ C) / max(nrm(Yf), 1e-300)


# ---------------------------------------------------------------- 1. primitives
hdr("1. Primitive set: covariance, Q|q factorization, rationality")
LETTERS = [
    # name, function, group
    ("E'", Ep, 'cov'), ("E'_Q", EpQ, 'cov'), ("E'_q", Epq, 'cov'),
    *[("tau_" + NAME[s], TAU[s], 'cov') for s in range(4)],
    *[("1-tau_" + NAME[s], COTAU[s], 'cov') for s in range(4)],
    ("sigma", SIGMA, 'cov'),
    ("Pi+", lin((0.5, ID), (0.5, Ep)), 'cov'), ("Pi-", lin((0.5, ID), (-0.5, Ep)), 'cov'),
    ("L_P3", L_b(P3), 'cov'), ("R_P3", R_b(P3), 'cov'), ("Ad_P3", Ad_b(P3), 'cov'),
    ("L_P2", L_b(P2), 'cov'), ("R_P2", R_b(P2), 'cov'), ("Ad_P2", Ad_b(P2), 'cov'),
    ("Ad_POm", Ad_b(POm), 'cov'),
    ("Tw", TW, 'cov'), ("Tw_Q", TW_Q, 'cov'), ("Tw_q", TW_q, 'cov'),
    ("E", E_plain, 'non'), ("F_Q", F_Q, 'non'), ("F_q", F_q, 'non'),
    ("T", T_full, 'non'), *[("T_" + NAME[s], T_slot(s), 'non') for s in range(4)],
    ("L_PL1", L_b(PL1), 'non'), ("R_PL1", R_b(PL1), 'non'), ("Ad_PL1", Ad_b(PL1), 'non'),
    ("L_PL2", L_b(PL2), 'non'), ("R_PL2", R_b(PL2), 'non'), ("Ad_PL2", Ad_b(PL2), 'non'),
]
NL = len(LETTERS)


def haar_su(n):
    Z = (rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))) / np.sqrt(2)
    Q, R = np.linalg.qr(Z); Q = Q @ np.diag(np.diag(R) / abs(np.diag(R)))
    return Q / np.linalg.det(Q) ** (1 / n)


U, V = haar_su(3), haar_su(2)
W = rec.kron_all([U, V, U.conj(), V.conj()])
gW = Ad_b(W)
Atest = rng.normal(size=(36, 36, 3)) + 1j * rng.normal(size=(36, 36, 3))
ID1296 = np.eye(1296).reshape(36, 36, 1296)


def full_matrix(f):
    return f(ID1296.astype(complex)).reshape(1296, 1296)


def schmidt_rank_Qq(Mx):
    """operator-Schmidt rank of a superoperator across B(Q1 Q2) | B(q1 q2)"""
    T = Mx.reshape(SH + SH)                                   # out(8 axes), in(8 axes)
    Qax = [0, 2, 4, 6, 8, 10, 12, 14]; qax = [1, 3, 5, 7, 9, 11, 13, 15]
    Rm = T.transpose(Qax + qax).reshape(3 ** 8, 2 ** 8)
    s = np.linalg.svd(Rm, compute_uv=False)
    return int((s > 1e-9 * s[0]).sum())


def min_denominator(Mx):
    """least common denominator of the entries (real and imaginary parts), if all are rationals 2^a 3^b-denominated"""
    vals = np.unique(np.round(np.concatenate([Mx.real.ravel(), Mx.imag.ravel()]), 12))
    D = 1
    for x in vals:
        fr = Fraction(float(x)).limit_denominator(10 ** 6)
        if abs(float(fr) - x) > 1e-10:
            return None
        D = D * fr.denominator // np.gcd(D, fr.denominator)
    m = D
    for p in (2, 3):
        while m % p == 0:
            m //= p
    return D if m == 1 else None


FACT, DEN = {}, {}
print(f"{'letter':9s} {'G_Om-cov defect':>15s} {'E\'-cov':>8s} {'Schmidt(Q|q)':>12s} {'denominator':>11s}  group")
for name, f, grp in LETTERS:
    cov = nrm(f(gW(Atest)) - gW(f(Atest))) / nrm(f(Atest)) if nrm(f(Atest)) > 0 else 0.0
    epc = nrm(f(Ep(Atest)) - Ep(f(Atest))) / max(nrm(f(Atest)), 1e-300)
    Mx = full_matrix(f)
    sr = schmidt_rank_Qq(Mx)
    den = min_denominator(Mx)
    FACT[name] = (sr == 1); DEN[name] = den
    assert (cov < 1e-10) == (grp == 'cov'), (name, cov)
    print(f"{name:9s} {cov:15.1e} {epc:8.1e} {sr:12d} {str(den):>11s}  {grp}")
print("all letters have entries in (1/(2^a 3^b)) Z[i]  (Galois lemma premise):",
      all(v is not None for v in DEN.values()))
print("Q|q-factorizing letters:", sum(FACT.values()), "of", NL,
      "; non-factorizing:", [n for n, v in FACT.items() if not v])
# the amplitude 'summed over matched directions'
for d in (2, 3):
    Phi = np.eye(d).reshape(-1); om = Phi / np.sqrt(d)
    rho = np.eye(d) / d
    print(f"d={d}: <Phi_d|omega_d> = {Phi @ om:.12f} = sqrt(d); Tr sqrt(rho_d) = {np.trace(np.sqrt(rho)):.12f};"
          f" |Phi_d><omega_d| = sqrt(d) P_d: {np.allclose(np.outer(Phi, om), np.sqrt(d) * np.outer(om, om))}")

# ---------------------------------------------------------------- 2. symmetry structure on R
hdr("2. Isotypic decomposition of R (record basis), commutant, Letter's point")
GEN = []
for K in sym.K3 + sym.K2:
    GEN.append(compress(lambda Y, K=K: np.einsum('ij,jkr->ikr', K, Y) - np.einsum('ijr,jk->ikr', Y, K))[0])
EPR, lk = compress(Ep)
assert lk < 1e-12
C3 = sum(G @ G for G in GEN[:8]); C2 = sum(G @ G for G in GEN[8:])
Hmix = C3 + np.pi * C2 + 100 * np.e * EPR
w, Vv = np.linalg.eigh((Hmix + Hmix.conj().T) / 2)
groups = defaultdict(list)
for j in range(NR):
    v = Vv[:, j]
    lab = (round((v.conj() @ C3 @ v).real, 4), round((v.conj() @ C2 @ v).real, 4), int(round((v.conj() @ EPR @ v).real)))
    groups[lab].append(j)
IR3 = {0.0: '1', 3.0: '8', 6.0: '10', 8.0: '27'}; IR2 = {0.0: '1', 2.0: '3', 6.0: '5'}
D3 = {'1': 1, '8': 8, '10': 10, '27': 27}; D2 = {'1': 1, '3': 3, '5': 5}
gens3 = [sum(rng.normal() * G for G in GEN) for _ in range(3)]


def commutant_basis(Vb):
    n = Vb.shape[1]
    rows = []
    for g in gens3:
        h = Vb.conj().T @ g @ Vb
        rows.append(np.kron(h, np.eye(n)) - np.kron(np.eye(n), h.T))
    Mt = np.vstack(rows)
    _, s, vh = np.linalg.svd(Mt, full_matrices=False)
    null = vh[(s > 1e-8).sum():].conj()
    return [x.reshape(n, n) for x in null]


ISO = []                                  # (label, Vb, multiplicity)
for lab in sorted(groups):
    Vb = Vv[:, groups[lab]]
    r3, r2 = IR3[lab[0]], IR2[lab[1]]
    if r3 == '10':                        # split 10 and 10bar with a generic element of the commutant
        Zs = commutant_basis(Vb)
        assert len(Zs) == 2
        Z = sum(rng.normal() * z for z in Zs); Z = Z + Z.conj().T
        wz, vz = np.linalg.eigh(Z)
        cut = np.argmax(np.diff(wz)) + 1
        for nm, sel in (('10', slice(0, cut)), ('10b', slice(cut, None))):
            Vs = Vb @ vz[:, sel]
            ISO.append(((nm, r2, lab[2]), Vs, Vs.shape[1] // (10 * D2[r2])))
    else:
        ISO.append(((r3, r2, lab[2]), Vb, Vb.shape[1] // (D3[r3] * D2[r2])))
COMM = []
for lab, Vb, m in ISO:
    Zs = commutant_basis(Vb)
    assert len(Zs) == m * m, (lab, len(Zs), m)
    COMM += [Vb @ z @ Vb.conj().T for z in Zs]
Zmat = np.array([z.reshape(-1) for z in COMM])
q, _ = np.linalg.qr(Zmat.T)
Zmat = q.T                                 # orthonormal basis of the commutant (rows)
for z in Zmat:
    zz = z.reshape(NR, NR)
    assert max(nrm(zz @ G - G @ zz) for G in GEN + [EPR]) < 1e-8
print("isotypic components (SU3, SU2, E'):", len(ISO), " commutant dim:", Zmat.shape[0])


def twirl(C):
    return (Zmat.T @ (Zmat.conj() @ C.reshape(-1))).reshape(NR, NR)


ISO_LABELS = [f"({a},{b}){'+' if e > 0 else '-'}" for (a, b, e), _, _ in ISO]
ISO_P = [Vb @ Vb.conj().T for _, Vb, _ in ISO]


ISO_PT = np.array([P.T for P in ISO_P])


def iso_traces(C):
    return np.einsum('kij,ij->k', ISO_PT, C)


ncl = {c: [int(round(np.trace(P @ PCL[c]).real)) for P in ISO_P] for c in CLASSES}
r2_, r3_, r6_ = sp.sqrt(2), sp.sqrt(3), sp.sqrt(6)
SQK = {'odd': 1, 'evmix': r2_, 'evsame': r3_, 'comm': r6_}
print(f"{'component':12s} {'mult':>4s} {'dim':>4s}  class content (odd,evmix,evsame,comm)   tr(M_L) exact      tr(K_L)")
tM_exact = []
for i, (lab, Vb, m) in enumerate(ISO):
    content = tuple(ncl[c][i] for c in CLASSES)
    tex = sum(content[j] * SQK[c] for j, c in enumerate(CLASSES))
    tM_exact.append(tex)
    print(f"{ISO_LABELS[i]:12s} {m:4d} {Vb.shape[1]:4d}  {str(content):38s} {str(sp.nsimplify(tex)):18s} "
          f"{sum(content[j] * KOF[c] for j, c in enumerate(CLASSES))}")
tM = iso_traces(ML); tK = iso_traces(KL)
assert np.allclose(tM, [float(x) for x in tM_exact])
print("M_L in the commutant of G_Omega x <E'>: defect", f"{nrm(ML - twirl(ML)) / nrm(ML):.1e}")
print("sum of per-isotypic traces of M_L:", sp.nsimplify(sum(tM_exact)))
print("(1,3)- component: Letter weight 1 on the q1q2-odd copy, sqrt6 on the commutator copy; the two copies are\n"
      "  separated by the slot-occupation number N (2 vs 1), which is built from partial traces.")

# ---------------------------------------------------------------- 3. exact identities
hdr("3. Exact identities: M_L from natural projectors; factorization obstruction")
S_glob = lin((0.5, ID), (0.5, compose(SIGMA, Ep)))       # survival: E' returns A times its own sign (-1)^N
PI_X = lambda Y: TYPE_OPS['QQ'](Y) + TYPE_OPS['qq'](Y) + TYPE_OPS['comm'](Y)   # letters at one end type only
PI_Y = lambda Y: TYPE_OPS['mix'](Y) + TYPE_OPS['comm'](Y)                       # exactly one occupied qubit slot
SX = compose(S_glob, PI_X); SY = compose(S_glob, PI_Y)
M_nat = compose(lin((1, ID), (np.sqrt(3) - 1, SX)), lin((1, ID), (np.sqrt(2) - 1, SY)))
for nm, f in (('S', S_glob), ('Pi_X', PI_X), ('Pi_Y', PI_Y), ('M_nat', M_nat)):
    C, lk = compress(f)
    if nm == 'M_nat':
        print(f"M_nat = (1+(sqrt3-1) S Pi_X)(1+(sqrt2-1) S Pi_Y) on R: leakage {lk:.1e}, "
              f"||M_nat - M_L|| = {nrm(C - ML):.1e}")
    else:
        print(f"{nm:5s}: maps R into R (leakage {lk:.1e}), diagonal on records: {nrm(C - np.diag(np.diag(C))) < 1e-12}")
print("  S = (1 + (-1)^N E')/2 with N = sum_s (1 - tau_s);  Pi_X = [N_Q N_q = 0];  Pi_Y = [N_q = 1].")
print("  On R the projector Pi_X contains the 9 q1q2 records and the 3 commutator records (no qutrit letter).")
# Schmidt rank of the mixed block across W_Q | W_q
mixidx = np.flatnonzero(rtype == 'mix')
Bmix = B[:, mixidx]


def mixed_block_schmidt(C172):
    """Schmidt rank across (single-Q-letter ops) x (single-q-letter ops) of the mixed block, via the full operator"""
    X = Bmix @ C172[np.ix_(mixidx, mixidx)] @ Bmix.conj().T     # superoperator on B(H) supported on the mixed space
    return schmidt_rank_Qq(X)


for nm, C in (('M_L', ML), ('K_L', KL), ("E'", EPR), ('identity', np.eye(NR))):
    print(f"Schmidt rank (Q|q) of the mixed block of {nm:8s}: {mixed_block_schmidt(C)}")
print("  A Q|q-product superoperator has a product mixed block (rank 1); class-scalar products therefore have\n"
      "  |mu(E'-even mixed)| = |mu(E'-odd mixed)|, but M_L needs sqrt2 : 1 and K_L needs 2 : 1.")

# ---------------------------------------------------------------- 4. family F1: words
hdr("4. Family F1: words of length <= 3 in the 38 letters, acting on R, compressed to R")


def canon_key(Y):
    flat = Y.reshape(-1)
    a = np.abs(flat); m = a.max()
    if m < 1e-10:
        return None
    ref = flat[np.flatnonzero(a > m * (1 - 1e-6))[0]]
    z = np.round(flat / ref, 7)
    return (np.round(z.real, 7) + 0.0).tobytes() + (np.round(z.imag, 7) + 0.0).tobytes()


def proj_resid(x, y):
    """min_c ||x - c y|| / ||x||"""
    nx = nrm(x)
    if nx < 1e-14:
        return 1.0
    c = (y.conj() @ x) / (y.conj() @ y)
    return nrm(x - c * y) / nx


CLS_I = np.array([CLASSES.index(c) for c in cls])


def class_scalar(C):
    """returns (residual, values) for the best class-scalar operator sum_c v_c Pi_c"""
    d = np.diag(C)
    v = np.array([d[CIDX[c]].mean() for c in CLASSES])
    return nrm(C - np.diag(v[CLS_I])) / max(nrm(C), 1e-300), v


sqk = np.sqrt(np.array([1, 2, 3, 6.0])); kk = np.array([1, 2, 3, 6.0])
MIXE = np.flatnonzero((rtype == 'mix') & (epar == 1)); MIXO = np.flatnonzero((rtype == 'mix') & (epar == -1))
GEN_TEST = gens3 + [EPR]


def analyse(C):
    nC = nrm(C)
    out = {}
    out['inv'] = max(nrm(C @ g - g @ C) / (nC * nrm(g)) for g in gens3)
    out['epinv'] = nrm(C @ EPR - EPR @ C) / (nC * nrm(EPR))
    d = np.diag(C)
    out['offdiag'] = nrm(C - np.diag(d)) / nC
    dm = np.array([d[CIDX[c]].mean() for c in CLASSES])
    out['dspread'] = max(np.abs(d[CIDX[c]] - dm[i]).max() for i, c in enumerate(CLASSES)) / max(np.abs(d).max(), 1e-300)
    out['dmeans'] = dm
    dvec_M = np.sqrt(k); dvec_K = k
    out['rM'] = proj_resid(d, dvec_M.astype(complex)); out['rK'] = proj_resid(d, dvec_K.astype(complex))
    t = iso_traces(C)
    out['rtM'] = proj_resid(t, tM); out['rtK'] = proj_resid(t, tK)
    out['cs'], out['csv'] = class_scalar(C)
    Ct = twirl(C)
    if nrm(Ct) > 1e-10 * nC:
        out['cs_tw'], out['csv_tw'] = class_scalar(Ct)
    else:
        out['cs_tw'], out['csv_tw'] = 1.0, np.zeros(4)
    # functional calculus: |C| = (C^dag C)^(1/2) for every word; sqrt(C) when C is Hermitian PSD
    wv, vv = np.linalg.eigh(C.conj().T @ C)
    Cabs = (vv * np.sqrt(np.clip(wv, 0, None))) @ vv.conj().T
    out['abs_cs'], out['abs_csv'] = class_scalar(Cabs)
    out['abs_rM'] = proj_resid(np.diag(Cabs), np.sqrt(k).astype(complex))
    out['abs_rtM'] = proj_resid(iso_traces(Cabs), tM)
    out['abs_dmeans'] = np.array([np.diag(Cabs)[CIDX[c]].mean() for c in CLASSES])
    herm = nrm(C - C.conj().T) / nC < 1e-9
    out['psd'] = False
    if herm:
        wh, vh = np.linalg.eigh((C + C.conj().T) / 2)
        if wh.min() > -1e-9 * np.abs(wh).max():
            out['psd'] = True
            Cs = (vh * np.sqrt(np.clip(wh, 0, None))) @ vh.conj().T
            out['sq_cs'], out['sq_csv'] = class_scalar(Cs)
            out['sq_rM'] = proj_resid(np.diag(Cs), np.sqrt(k).astype(complex))
            out['sq_rtM'] = proj_resid(iso_traces(Cs), tM)
            out['sq_dmeans'] = np.array([np.diag(Cs)[CIDX[c]].mean() for c in CLASSES])
    # mixed block test (factorization theorem): E'-scalar on the mixed records after twirl?
    Cm = Ct[np.ix_(mixidx, mixidx)]
    dCt = np.diag(Ct)
    mu_e = dCt[MIXE].mean(); mu_o = dCt[MIXO].mean()
    Sm = np.diag(np.where(epar[mixidx] == 1, mu_e, mu_o))
    out['mix_scalar'] = nrm(Cm - Sm) / max(nrm(Cm), 1e-300) if nrm(Cm) > 1e-10 * nC else 1.0
    out['mu'] = (mu_e, mu_o)
    return out


t_words = time.time()
alpha_names = [n for n, _, _ in LETTERS]
seen = set()
nwords = Counter(); nnew = Counter()
results = []                                 # (word tuple, analysis or None)


def visit(w, Z):
    """record a new image Z of R under word w (letters applied right to left); return True if new"""
    key = canon_key(Z)
    if key is None:
        return False
    h = hashlib.sha1(key).digest()
    if h in seen:
        return False
    seen.add(h)
    nnew[len(w)] += 1
    Yf = Z.reshape(1296, NR)
    C = Bh @ Yf
    if nrm(C) < 1e-10 * nrm(Yf):
        results.append((w, None)); return True
    a = analyse(C)
    a['leak'] = nrm(Yf - B @ C) / nrm(Yf)
    a['fact'] = all(FACT[alpha_names[i]] for i in w)
    a['cov'] = all(LETTERS[i][2] == 'cov' for i in w)
    results.append((w, a))
    return True


visit((), B3)                                # the empty word: identity of R (the bare insertion)
nwords[0] = 1
lev1, lev2 = [], []
for li, (_, f, _) in enumerate(LETTERS):
    nwords[1] += 1
    Z = f(B3)
    if visit((li,), Z):
        lev1.append(((li,), Z))
for w1, Y1 in lev1:                          # words of length 1 that were not new are scalar multiples of earlier ones
    for li, (_, f, _) in enumerate(LETTERS):
        nwords[2] += 1
        if visit((li,) + w1, f(Y1)):
            lev2.append((li,) + w1)
for w2 in lev2:                              # recompute each new length-2 image, then extend it on the left
    Y2 = LETTERS[w2[0]][1](LETTERS[w2[1]][1](B3))
    for li, (_, f, _) in enumerate(LETTERS):
        nwords[3] += 1
        visit((li,) + w2, f(Y2))
# words whose length-(L-1) right factor duplicates an earlier image are skipped: their images coincide (up to
# scale) with left-extensions of that earlier image, which are enumerated
for L in (0, 1, 2, 3):
    print(f"length {L}: {nwords[L]:6d} words formed, {nnew[L]:6d} new distinct nonzero images of R (up to scale)")
print("distinct images in F1:", len(results), "  (words of length <= 3 over 38 letters, with the empty word:",
      1 + 38 + 38 ** 2 + 38 ** 3, ")")
print(f"word analysis time {time.time() - t_words:.0f} s", file=sys.stderr)
wname = lambda w: ' . '.join(alpha_names[i] for i in w) if w else '1 (identity)'
good = [(w, a) for w, a in results if a is not None]
print("images with zero compression onto R (word leaves R entirely):", sum(a is None for _, a in results))
print("analysed words (distinct, nonzero on R):", len(good))
print("  all letters covariant:", sum(a['cov'] for _, a in good),
      "| compression G_Omega-invariant (defect<1e-9):", sum(a['inv'] < TOL for _, a in good),
      "| also E'-invariant:", sum((a['inv'] < TOL) and (a['epinv'] < TOL) for _, a in good))
print("  maps R into R (leakage<1e-9):", sum(a['leak'] < TOL for _, a in good),
      "| diagonal in the record basis:", sum(a['offdiag'] < TOL for _, a in good))


def rat_pattern(v):
    """normalize a class-value vector and identify it as rationals"""
    j = int(np.argmax(np.abs(v) > 1e-9 * np.abs(v).max()))
    u = v / v[j]
    if np.abs(u.imag).max() > 1e-9:
        return None
    fr = tuple(Fraction(float(x)).limit_denominator(10000) for x in u.real)
    if max(abs(float(f) - x) for f, x in zip(fr, u.real)) > 1e-9:
        return ('irrational',) + tuple(np.round(u.real, 6))
    return fr


def fmt_pat(p):
    return '(' + ', '.join(str(x) for x in p) + ')' if p is not None else 'None'


for tag, cskey, vkey in (('raw C', 'cs', 'csv'), ("twirled over G x <E'>", 'cs_tw', 'csv_tw'),
                         ('|C| = (C^dag C)^1/2', 'abs_cs', 'abs_csv')):
    pats = Counter()
    for w, a in good:
        if a[cskey] < TOL and np.abs(a[vkey]).max() > 1e-9:
            pats[rat_pattern(a[vkey])] += 1
    print(f"\n[{tag}] class-scalar operators (C = sum_c v_c Pi_c on R): {sum(pats.values())} words, "
          f"{len(pats)} distinct patterns (v_odd, v_evmix, v_evsame, v_comm) up to scale:")
    for p, n in sorted(pats.items(), key=lambda x: -x[1]):
        print("   ", fmt_pat(p), "x", n)
    print("   pattern prop. to (1, sqrt2, sqrt3, sqrt6):",
          sum(n for p, n in pats.items() if p is not None and p[0] == 'irrational'
              and np.allclose(np.array(p[1:], float), sqk / sqk[0], atol=1e-6)),
          "  prop. to (1, 2, 3, 6):", sum(n for p, n in pats.items()
                                          if p is not None and p[0] != 'irrational' and tuple(p) == (1, 2, 3, 6)))
sqs = [(w, a) for w, a in good if a['psd'] and a['sq_cs'] < TOL and np.abs(a['sq_csv']).max() > 1e-9]
print(f"\n[sqrt(C), C Hermitian PSD] {sum(a['psd'] for _, a in good)} PSD words; class-scalar sqrt: {len(sqs)};"
      f" patterns: {Counter(fmt_pat(rat_pattern(a['sq_csv'])) for _, a in sqs).most_common()}")

# diagonal-only criterion (b'): diagonal constant on the four classes
dconst = [(w, a) for w, a in good if a['dspread'] < TOL and np.abs(a['dmeans']).max() > 1e-9]
print(f"\n(b') diagonal constant on each class (operator need not be diagonal): {len(dconst)} words; patterns:")
for p, n in Counter(fmt_pat(rat_pattern(a['dmeans'])) for _, a in dconst).most_common():
    print("   ", p, "x", n)

# factorization theorem check on the mixed block
fact_mix = [(w, a) for w, a in good if a['fact'] and a['mix_scalar'] < TOL and max(abs(a['mu'][0]), abs(a['mu'][1])) > 1e-9]
viol = [(w, a) for w, a in fact_mix if abs(abs(a['mu'][0]) - abs(a['mu'][1])) > 1e-9 * max(abs(a['mu'][0]), abs(a['mu'][1]))]
nonfact_mix = [(w, a) for w, a in good if (not a['fact']) and a['mix_scalar'] < TOL and max(abs(a['mu'][0]), abs(a['mu'][1])) > 1e-9]
print(f"\nfactorization check: Q|q-factorizing words whose twirled mixed block is E'-scalar: {len(fact_mix)};"
      f" with |mu_even| != |mu_odd|: {len(viol)}")
rat = Counter()
for w, a in nonfact_mix:
    me, mo = a['mu']
    rat['inf' if abs(mo) < 1e-12 else str(Fraction(float((me / mo).real)).limit_denominator(1000))] += 1
print(f"non-factorizing words with E'-scalar twirled mixed block: {len(nonfact_mix)}; ratios mu_even/mu_odd: {dict(rat)}")
exotic = [(w, a) for w, a in nonfact_mix if abs(a['mu'][1]) > 1e-12 and
          min(abs(a['mu'][0] / a['mu'][1] - x) for x in (0, 1, -1)) > 1e-9]
print(f"  the {len(exotic)} words with mu_even/mu_odd not in {{0, 1, -1, inf}}: every one uses a read link:",
      all(any(alpha_names[i].endswith(('PL1', 'PL2')) for i in w) for w, _ in exotic),
      "; class-scalar after twirl:", sum(a['cs_tw'] < TOL for _, a in exotic))
for w, a in sorted(exotic, key=lambda x: (abs((x[1]['mu'][0] / x[1]['mu'][1]).real), len(x[0])))[::9]:
    v = a['csv_tw'] / a['csv_tw'][0] if abs(a['csv_tw'][0]) > 1e-12 else a['csv_tw']
    print(f"   e.g. [{' . '.join(alpha_names[i] for i in w)}]  mu_e/mu_o = {(a['mu'][0] / a['mu'][1]).real:+.4f};"
          f" twirled class means (odd, evmix, evsame, comm) rel. to odd {np.round(v.real, 4)};"
          f" class-scalar residual after twirl {a['cs_tw']:.3f}")

# closest misses
A_ID = [a for w, a in good if w == ()][0]


def show(title, key, n=6, filt=lambda a: True, dkey='dmeans'):
    print(f"\n{title}   [identity: {A_ID[key]:.4f}]")
    lst = sorted([(a[key], w, a) for w, a in good if key in a and filt(a)], key=lambda x: x[0])[:n]
    for val, w, a in lst:
        dm = a[dkey] / a[dkey][0] if abs(a[dkey][0]) > 1e-12 else a[dkey]
        print(f"   resid {val:.4f}  word [{wname(w)}]  G-inv defect {a['inv']:.0e}  class means of diagonal"
              f" (odd, evmix, evsame, comm), rel. to odd: {np.round(dm.real, 4)}")
    return lst


show("closest misses, (b) diagonal vs c*sqrt(k_r)  [Letter: (1, 1.4142, 1.7321, 2.4495)]", 'rM')
show("closest misses, (b) diagonal vs c*k_r  [K-route: (1, 2, 3, 6)]", 'rK')
show("closest misses, (c) per-isotypic traces vs the Letter's point", 'rtM')
show("closest misses, (c) via |C|: per-isotypic traces of |C| vs the Letter's point", 'abs_rtM', dkey='abs_dmeans')
show("closest misses, (b) diagonal of |C| vs c*sqrt(k_r)", 'abs_rM', dkey='abs_dmeans')
show("closest misses, (b) diagonal of sqrt(C) (C Hermitian PSD) vs c*sqrt(k_r)", 'sq_rM', dkey='sq_dmeans')
show("closest misses, (c) per-isotypic traces of sqrt(C) (C Hermitian PSD) vs the Letter's point", 'sq_rtM', dkey='sq_dmeans')
show("closest misses, (b)+(c) among G_Omega-invariant words: per-isotypic traces", 'rtM', filt=lambda a: a['inv'] < TOL)
best_t = show("best (c) miss, detail:", 'rtM', n=1)[0]
Cb = Bh @ (compose(*[LETTERS[i][1] for i in best_t[1]])(B3)).reshape(1296, NR)
tb = iso_traces(Cb); cfit = (tM.conj() @ tb) / (tM.conj() @ tM)
print("   per-isotypic traces, Letter vs best miss (rescaled):")
for lab, x, y in zip(ISO_LABELS, tM, tb / cfit):
    print(f"     {lab:9s} {x.real:9.4f} {y.real:9.4f}")
KEYS = ('rM', 'rK', 'rtM', 'abs_rM', 'abs_rtM', 'sq_rM', 'sq_rtM')
print("\nminimum residual over F1 / residual of the identity (0 = exact reproduction):")
for key in KEYS:
    print(f"   {key:8s} {min(a[key] for _, a in good if key in a):.4f} / {A_ID[key]:.4f}")
print("   (rM, rK: diagonal vs sqrt(k), k; rtM: per-isotypic traces vs the Letter's point; abs_/sq_: same for |C|, sqrt(C))")

# ---------------------------------------------------------------- 5. family F2: survival x attribution
hdr("5. Family F2: survival x attribution dressings (operator form of the Letter's mechanism)")
EQR, l1 = compress(EpQ); EqR, l2 = compress(Epq); SGR, l3 = compress(SIGMA)
TPR = {t: compress(f) for t, f in TYPE_OPS.items()}
assert max(l1, l2, l3, *[v[1] for v in TPR.values()]) < 1e-12
assert nrm(EQR @ EqR - EPR) < 1e-12 and nrm(EQR @ EqR - EqR @ EQR) < 1e-12
print("E'_Q, E'_q, sigma and the four occupation projectors map R into R; E' = E'_Q E'_q on R:",
      nrm(EQR @ EqR - EPR) < 1e-12)
Hj = EPR + np.pi * EQR + np.e * EqR + 7.1 * SGR + sum((i + 1) * 11.3 * TPR[t][0] for i, t in enumerate(TYPE_OPS))
wj, Vj = np.linalg.eigh((Hj + Hj.conj().T) / 2)
labs = []
for j in range(NR):
    v = Vj[:, j]
    rq = lambda X: int(round((v.conj() @ X @ v).real))
    t = [tt for tt in TYPE_OPS if rq(TPR[tt][0]) == 1][0]
    labs.append((t, rq(EPR), rq(EQR), rq(EqR), rq(SGR)))
LABC = Counter(labs)
print("joint eigen-labels (type, E', E'_Q, E'_q, sigma): dim")
for l, n in sorted(LABC.items()):
    print("   ", l, n)
ULAB = sorted(LABC)
OVL = np.abs(Vj) ** 2                      # |<record r | joint eigvec j>|^2
SURV = {
    '1': lambda l: 1,
    "Pi+(E')": lambda l: int(l[1] == 1), "Pi+(sigma E')": lambda l: int(l[1] * l[4] == 1),
    "Pi+(E'_Q)": lambda l: int(l[2] == 1), "Pi+(sigma E'_Q)": lambda l: int(l[2] * l[4] == 1),
    "Pi+(E'_q)": lambda l: int(l[3] == 1), "Pi+(sigma E'_q)": lambda l: int(l[3] * l[4] == 1),
}
TYPES = ['QQ', 'qq', 'mix', 'comm']
HAS = {'Q': {'QQ', 'mix'}, 'q': {'qq', 'mix', 'comm'}}
DIM = {'Q': 3, 'q': 2}
ATTR = list(itertools.product([None, 'Q', 'q'], repeat=4))


def factor_vec(sname, attr):
    """K-factor (d or 1) on each label class for one exchange"""
    s = SURV[sname]
    return tuple((DIM[attr[TYPES.index(l[0])]] if (attr[TYPES.index(l[0])] is not None and s(l)) else 1) for l in ULAB)


def is_local(attr):
    return all(a is None or TYPES[i] in HAS[a] for i, a in enumerate(attr))


def letter_class(l):
    if l[0] == 'comm':
        return 'comm'
    return 'odd' if l[1] == -1 else ('evsame' if l[0] in ('QQ', 'qq') else 'evmix')


ML_lab = tuple(KOF[letter_class(l)] for l in ULAB)
singles = {}
for sname in SURV:
    for attr in ATTR:
        singles.setdefault(factor_vec(sname, attr), []).append((sname, attr))
print("distinct single-exchange factors:", len(singles), " (7 survival operators x 81 attributions =", 7 * 81, ")")
dress = {}
for f1, c1 in singles.items():
    for f2, c2 in singles.items():
        kv = tuple(a * b for a, b in zip(f1, f2))
        if kv not in dress:
            dress[kv] = (c1[0], c2[0], any(is_local(a) for _, a in c1) and any(is_local(a) for _, a in c2))
        elif not dress[kv][2] and any(is_local(a) for _, a in c1) and any(is_local(a) for _, a in c2):
            dress[kv] = (c1[0], c2[0], True)
print("ordered pairs of exchanges:", (7 * 81) ** 2, " -> distinct dressings K on R:", len(dress))
dimv = np.array([LABC[l] for l in ULAB])


SQRT_COORD = {1: (1, 0, 0, 0), 2: (0, 1, 0, 0), 3: (0, 0, 1, 0), 4: (2, 0, 0, 0), 6: (0, 0, 0, 1), 9: (3, 0, 0, 0)}


def numerator(kv):
    """exact sum of sqrt(K) over R, as integer coordinates on (1, sqrt2, sqrt3, sqrt6)"""
    c = np.zeros(4, int)
    for kval, n in zip(kv, dimv):
        c += n * np.array(SQRT_COORD[kval])
    return tuple(int(x) for x in c)


def num_str(c):
    return f"{c[0]}+{c[1]}r2+{c[2]}r3+{c[3]}r6"


num_float = lambda c: c[0] + c[1] * 2 ** 0.5 + c[2] * 3 ** 0.5 + c[3] * 6 ** 0.5


AINV = 137.035999177
N = 31104
rows = []
for kv, (c1, c2, loc) in dress.items():
    num = numerator(kv)
    inv = N / num_float(num)
    diag_on_records = all(len({kv[i] for i, l in enumerate(ULAB) if (l[0], l[1]) == key}) == 1
                          for key in {(l[0], l[1]) for l in ULAB})
    rows.append((kv, num, inv, loc, diag_on_records, c1, c2))
Mw_of = lambda kv: (Vj * np.sqrt(np.array([kv[ULAB.index(l)] for l in labs], float))) @ Vj.conj().T
exact = [r for r in rows if r[0] == ML_lab]
print("dressings equal to M_L on R:", len(exact))
nloc = sum(r[3] for r in rows)
print("dressings with a local attribution (d_e credited only to records with a letter at an e-slot):", nloc,
      "; of these equal to M_L:", sum(r[3] for r in exact))
print("number of (survival, attribution) pairs giving M_L, and two of them:")
cnt = 0; ex_list = []
for f1, c1 in singles.items():
    for f2, c2 in singles.items():
        if tuple(a * b for a, b in zip(f1, f2)) == ML_lab:
            cnt += len(c1) * len(c2); ex_list.append((c1[0], c2[0]))
print("   ", cnt, "ordered choices;", ex_list[:2])
print("   every one credits d_Q = 3 to the E'-even q1q2 records and to the commutator records,"
      " which carry no qutrit letter:",
      all(not is_local(a1) or not is_local(a2)
          for f1, c1 in singles.items() for f2, c2 in singles.items()
          if tuple(a * b for a, b in zip(f1, f2)) == ML_lab for (_, a1) in c1 for (_, a2) in c2))
# Letter's Table S5 rows inside F2
S5 = {                                       # Supplement Table S5, weight-law rows (numerators over 31104)
    'computed (Letter)': (79, 48, 42, 3),
    'two ends swapped': (79, 42, 48, 3),
    'both at qubit end': (85, 90, 0, 0),
    'both at qutrit end': (88, 0, 90, 0),
    'X where letters sit': (85, 54, 36, 0),
    'qq moved by Y': (79, 48, 36, 9),
    'copies separate': (79, 0, 42, 51),
    'no own sign': (82, 48, 42, 0),
    'comm survive X only': (79, 48, 45, 0),
    'no weights': (172, 0, 0, 0),
}
print("Letter's Table S5 weight-law rows realized as F2 dressings (numerator found; local?):")
for nm, val in S5.items():
    m = [r for r in rows if r[1] == val]
    print(f"    {nm:22s} {num_str(val):22s} 1/value {N / num_float(val):9.3f}  in F2: {len(m) > 0}  local: {any(r[3] for r in m)}")
loc_rows = [r for r in rows if r[3]]
print("local dressings: distinct values", len({r[1] for r in loc_rows}),
      "; 1/value range [%.3f, %.3f]" % (min(r[2] for r in loc_rows), max(r[2] for r in loc_rows)))
print("local dressings closest to M_L as operators (||M - M_L||_HS / ||M_L||_HS):")
MLn = nrm(ML)
dl = sorted([(nrm(Mw_of(r[0]) - ML) / MLn, r) for r in loc_rows], key=lambda x: x[0])[:5]
for dist, r in dl:
    print(f"    dist {dist:.4f}  numerator {num_str(r[1]):22s} 1/value {r[2]:9.3f}  diagonal on records: {r[4]}  "
          f"X: {r[5]}  Y: {r[6]}")
print("local dressings closest in value to 1/alpha(CODATA 2022):")
for r in sorted(loc_rows, key=lambda r: abs(r[2] - AINV))[:5]:
    print(f"    1/value {r[2]:11.5f}  numerator {num_str(r[1]):22s} diagonal on records: {r[4]}  X: {r[5]}  Y: {r[6]}")
print("all F2 dressings: distinct values", len({r[1] for r in rows}),
      "; within 1e-3 relative of 1/alpha:", sum(abs(r[2] / AINV - 1) < 1e-3 for r in rows),
      "; of these local:", sum(abs(r[2] / AINV - 1) < 1e-3 and r[3] for r in rows))

lv = sorted({r[1] for r in loc_rows})
dens = sum(120 < N / num_float(c) < 155 for c in lv) / 35.0
print(f"local distinct values with 1/value in (120,155): {sum(120 < N / num_float(c) < 155 for c in lv)}"
      f" -> density {dens:.2f} per unit; density-window reference count in the +-1e-3 relative window: {dens * 2e-3 * AINV:.2f}")
print("\nSUMMARY")
print("  F1 (", 1 + 38 + 38 ** 2 + 38 ** 3, "words;", sum(nwords.values()), "formed explicitly, the rest are left extensions of",
      "duplicate images;", len(good), "distinct nonzero compressions): no word, twirled word, |word| or sqrt(word)",
      "reproduces sqrt(k_r) or k_r as an operator, its diagonal, or its per-isotypic traces.")
allpats = [rat_pattern(a[v]) for _, a in good for c, v in (('cs', 'csv'), ('cs_tw', 'csv_tw'), ('abs_cs', 'abs_csv'))
           if a[c] < TOL and np.abs(a[v]).max() > 1e-9]
allpats += [rat_pattern(a['sq_csv']) for _, a in good if a['psd'] and a['sq_cs'] < TOL]
print("  every class-scalar pattern found in F1 (raw, twirled, |C|, sqrt C;", len(allpats), "cases) has entries in {0, +1, -1}:",
      all(p is not None and p[0] != 'irrational' and all(x in (0, 1, -1) for x in p) for p in allpats))
print("  F2 (", (7 * 81) ** 2, "dressings,", len(dress), "distinct): M_L is reproduced exactly, but only with a non-local",
      "attribution (qutrit factor on records without qutrit letters).")
print(f"runtime {time.time() - T_START:.0f} s", file=sys.stderr)
