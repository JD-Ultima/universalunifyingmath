"""Independent coordinator cross-check of GU claim B (spinor_branching.py):
"at the SM level the 144 contains a fourth generation-equivalent".

Method (independent of spinor_branching.py, which uses weight multisets):
explicit Cl(10) on C^32, Weyl spinors S+ (16) and S- (16), vector V = C^10.
144 := kernel of Clifford multiplication  V (x) S+ -> S-  (dim 160 - 16).
SM embedding: so(10) > so(6)+so(4); su(3) = commutant in so(6) of the complex
structure J = e12+e34+e56 (traceless part); B-L = -(2/3) * (spinor generator of J)/i
up to sign fixed by requiring 16+ to carry one left-handed generation;
su(2)_L + su(2)_R = so(4); Y = (B-L)/2 + T3R.
We then decompose 16+ and the 144 into (SU(3) Casimir, SU(2)_L Casimir, Y) blocks.
Run: python3 xcheck_144.py
"""
import numpy as np
from collections import Counter
from itertools import combinations
sx = np.array([[0, 1], [1, 0]], complex); sy = np.array([[0, -1j], [1j, 0]]); sz = np.diag([1., -1]).astype(complex); I2 = np.eye(2)
def kron(*ms):
    out = np.array([[1.0 + 0j]])
    for m in ms: out = np.kron(out, m)
    return out
# Cl(10) generators, Hermitian, {g_i,g_j} = 2 delta
G = []
for k in range(5):
    pre = [sz] * k; post = [I2] * (4 - k)
    G.append(kron(*pre, sx, *post)); G.append(kron(*pre, sy, *post))
for i in range(10):
    for j in range(10):
        assert np.allclose(G[i] @ G[j] + G[j] @ G[i], 2 * (i == j) * np.eye(32))
chi = np.eye(32, dtype=complex)
for g in G: chi = chi @ g
chi = 1j * chi                      # (g0...g9)^2 = (-1)^45 = -1
assert np.allclose(chi @ chi, np.eye(32))
w, Vc = np.linalg.eigh((chi + chi.conj().T) / 2)
Pp = Vc[:, w > 0]; Pm = Vc[:, w < 0]           # bases of S+ and S-
assert Pp.shape[1] == 16 and Pm.shape[1] == 16
# so(10) generators: spinor M_ij = (i/4)[g_i,g_j]  (Hermitian), vector (M_ij)_{ab} = -i(d_ia d_jb - d_ib d_ja)
def Msp(i, j): return 0.25j * (G[i] @ G[j] - G[j] @ G[i])
def Mv(i, j):
    m = np.zeros((10, 10), complex); m[i, j] = -1j; m[j, i] = 1j; return m
# so(6) on indices 0..5, so(4) on 6..9.  J = M01 + M23 + M45.
pairs6 = [(0, 1), (2, 3), (4, 5)]
def gen(i, j): return Msp(i, j), Mv(i, j)
# su(3) inside so(6): commutant of J; build so(6) basis, project
so6 = [(i, j) for i, j in combinations(range(6), 2)]
Jv = sum(Mv(i, j) for i, j in pairs6)
basis6 = [Mv(i, j) for i, j in so6]
# find combinations commuting with Jv
A = np.array([(Jv @ b - b @ Jv).reshape(-1) for b in basis6]).T
_, s, vh = np.linalg.svd(A)
null = vh[np.sum(s > 1e-9):].conj()
u3 = [sum(c * gen(*p)[0] for c, p in zip(v, so6)) for v in null]  # spinor reps of u(3)
u3v = [sum(c * gen(*p)[1] for c, p in zip(v, so6)) for v in null]
assert len(u3) == 9
# remove the J direction (u(1)) to get su(3): orthogonalize w.r.t. Jv under trace form on vector rep
Jsp = sum(Msp(i, j) for i, j in pairs6)
su3 = []
for X, Xv in zip(u3, u3v):
    c = np.trace(Xv @ Jv) / np.trace(Jv @ Jv)
    su3.append(X - c * Jsp)
# su(2)_L, su(2)_R from so(4) on 6..9: T_L = (1/2)(M67+M89, M68-M79? ...) use standard self-dual combos
a, b, c_, d = 6, 7, 8, 9
L = [0.5 * (Msp(a, b) + Msp(c_, d)), 0.5 * (Msp(a, c_) - Msp(b, d)), 0.5 * (Msp(a, d) + Msp(b, c_))]
R = [0.5 * (Msp(a, b) - Msp(c_, d)), 0.5 * (Msp(a, c_) + Msp(b, d)), 0.5 * (Msp(a, d) - Msp(b, c_))]
for X in L:
    for Y in R: assert np.allclose(X @ Y, Y @ X)
def cas(gens, P):
    Ms = [P.conj().T @ X @ P for X in gens]
    return sum(M @ M for M in Ms)
# normalise su(3) so that on the vector-rep triplet the Casimir is 4/3 (fundamental)
def norm_su3(gens):
    # vector rep restricted: 6 = 3 + 3bar; Casimir sum X^2 on C^10
    Cv = sum(Xv @ Xv for Xv in [sum(cc * gen(*p)[1] for cc, p in zip(v, so6)) - np.trace(sum(cc * gen(*p)[1] for cc, p in zip(v, so6)) @ Jv) / np.trace(Jv @ Jv) * Jv for v in null])
    ev = np.linalg.eigvalsh((Cv + Cv.conj().T) / 2)
    top = max(ev)
    return [X * np.sqrt((4 / 3) / top) for X in gens]
su3 = norm_su3(su3)
def joint(C3, C2, Y, BL):
    H = 1.0 * C3 + 0.731 * C2 + 0.377 * Y + 0.0913 * BL
    H = (H + H.conj().T) / 2
    w, V = np.linalg.eigh(H)
    cnt = Counter()
    for k in range(len(w)):
        v = V[:, k]
        q = lambda M: round(float((v.conj() @ M @ v).real), 4)
        cnt[(q(C3), q(C2), q(Y), q(BL))] += 1
    return cnt
def blocks(P, BL_sign):
    C3 = cas(su3, P); C2 = cas(L, P)
    Jp = P.conj().T @ Jsp @ P
    T3R = P.conj().T @ R[0] @ P
    BL = BL_sign * (-(2 / 3) * Jp)
    Y = BL / 2 + T3R
    return joint(C3, C2, Y, BL)
names3 = {0.0: '1', 1.3333: '3', 3.0: '8', 3.3333: '6', 6.0: '10', 4.0: '?'}
names2 = {0.0: '1', 0.75: '2', 2.0: '3', 3.75: '4'}
def irreps(cnt):
    # SU(3) label from Casimir; 3 vs 3bar (6 vs 6bar, 10 vs 10bar) from triality = 3(B-L) mod 3
    out = Counter()
    for (c3, c2, y, bl), n in cnt.items():
        base = names3[c3]
        t = int(round(3 * bl)) % 3
        lab = base if base in ('1', '8') else (base if t == 1 else base + 'bar') if base in ('3', '6') else (base if t == 0 else base + '?')
        if base == '3' and t == 0 or base == '6' and t == 0: lab = base + '?'
        d3 = {'1': 1, '3': 3, '8': 8, '6': 6, '10': 10}[base]; d2 = {'1': 1, '2': 2, '3': 3, '4': 4}[names2[c2]]
        assert n % (d3 * d2) == 0, (c3, c2, y, bl, n)
        out[(lab, names2[c2], y)] += n // (d3 * d2)
    return out
# fix the sign of B-L: 16+ must contain (1,1)_{Y=1} (e^c) and (1,2)_{-1/2} (L), the left-handed generation
for BL_sign in (+1, -1):
    gen16 = irreps(blocks(Pp, BL_sign))
    if gen16.get(('1', '1', 1.0), 0) == 1 and gen16.get(('1', '2', -0.5), 0) == 1:
        break
print("B-L sign chosen:", BL_sign)
print("16+ content (SU3 Casimir label, SU2 label, Y):", dict(gen16))
# 144_eps = kernel of Clifford multiplication V (x) S_eps -> S_-eps, for eps = +, -
dim3 = lambda l: {'1': 1, '3': 3, '8': 8, '6': 6, '10': 10}[l.replace('bar', '').replace('?', '')]
conj_lab = lambda l: l if l in ('1', '8') else (l[:-3] if l.endswith('bar') else l + 'bar')
gen16_conj = Counter({(conj_lab(k[0]), k[1], round(-k[2], 4) + 0.0): n for k, n in gen16.items()})
def content_144(Pin, Pout):
    Cm = np.zeros((16, 160), complex)
    for a_ in range(10):
        Cm[:, a_ * 16:(a_ + 1) * 16] = Pout.conj().T @ G[a_] @ Pin
    _, sv, vh = np.linalg.svd(Cm)
    K = vh[np.sum(sv > 1e-9):].conj().T
    assert K.shape[1] == 144
    def tensor_gen(Xsp, Xv):
        # vector index transforms by -Mv so that Clifford multiplication is equivariant ([M,G_a] = sum_b G_b (-Mv)_{ba})
        return np.kron(-Xv, np.eye(16)) + np.kron(np.eye(10), Pin.conj().T @ Xsp @ Pin)
    su3_vs = []
    for v in null:
        Xsp = sum(cc * gen(*p)[0] for cc, p in zip(v, so6)); Xv = sum(cc * gen(*p)[1] for cc, p in zip(v, so6))
        c = np.trace(Xv @ Jv) / np.trace(Jv @ Jv)
        su3_vs.append(tensor_gen(Xsp - c * Jsp, Xv - c * Jv))
    v0 = null[0]
    raw0 = sum(cc * gen(*p)[0] for cc, p in zip(v0, so6)) - np.trace(sum(cc * gen(*p)[1] for cc, p in zip(v0, so6)) @ Jv) / np.trace(Jv @ Jv) * Jsp
    scale = np.linalg.norm(su3[0]) / np.linalg.norm(raw0)
    su3_vs = [X * scale for X in su3_vs]
    vgen = lambda i, j: tensor_gen(Msp(i, j), Mv(i, j))
    L_vs = [0.5 * (vgen(a, b) + vgen(c_, d)), 0.5 * (vgen(a, c_) - vgen(b, d)), 0.5 * (vgen(a, d) + vgen(b, c_))]
    R0_vs = 0.5 * (vgen(a, b) - vgen(c_, d))
    J_vs = sum(vgen(i, j) for i, j in pairs6)
    for X in su3_vs + L_vs + [R0_vs, J_vs]:
        assert np.linalg.norm(X @ K - K @ (K.conj().T @ X @ K)) < 1e-9     # the 144 is invariant
    proj = lambda M: K.conj().T @ M @ K
    C3 = sum(proj(X) @ proj(X) for X in su3_vs); C2 = sum(proj(X) @ proj(X) for X in L_vs)
    BLp = BL_sign * (-(2 / 3) * proj(J_vs))
    Y = BLp / 2 + proj(R0_vs)
    cnt = joint(C3, C2, Y, BLp)
    return irreps(cnt), cnt
for name, Pin, Pout in [("144_+ = V(x)S+ - S-", Pp, Pm), ("144_- = V(x)S- - S+", Pm, Pp)]:
    c144, cnt = content_144(Pin, Pout)
    assert sum(n * dim3(k[0]) * {'1': 1, '2': 2, '3': 3, '4': 4}[k[1]] for k, n in c144.items()) == 144
    has_gen = all(c144.get(k, 0) >= n for k, n in gen16.items())
    has_conj = all(c144.get(k, 0) >= n for k, n in gen16_conj.items())
    print(f"{name}: {len(c144)} SM irreps; contains SM content of one 16+ generation: {has_gen}; of its conjugate: {has_conj}")
    print("   ", dict(sorted(c144.items())))
    charges = sorted({round(y + t, 4) for (c3, c2, y, bl) in cnt for t in ([-0.5, 0.5] if c2 == 0.75 else [-1, 0, 1] if c2 == 2.0 else [-1.5, -0.5, 0.5, 1.5] if c2 == 3.75 else [0])})
    print("    electric charges Q = T3L + Y:", charges)
