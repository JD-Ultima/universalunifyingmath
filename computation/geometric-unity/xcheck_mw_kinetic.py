"""Independent coordinator cross-check of round-3 claim G (fermions): in signature (7,7) a
single Majorana-Weyl spinor admits no (Dirac-type) kinetic term.

Own construction (not clifford77.py): real Cl(1,1) on R^2 with e+ = sx (e+^2 = +1),
e- = [[0,1],[-1,0]] (e-^2 = -1), G = e+ e-; Cl(7,7) as the graded tensor power on R^128.
Spin-invariant bilinears: C+ = g1...g7 (positive gammas), C- = g8...g14 (negative gammas).
A kinetic term for a Grassmann spinor psi in S+ is int psi^T M d psi with M = C gamma^mu; it is
a total derivative unless the S+ x S+ block of M has a nonzero symmetric part.
Run: python3 xcheck_mw_kinetic.py
"""
import numpy as np
from functools import reduce
ep = np.array([[0, 1], [1, 0]], float); em = np.array([[0, 1], [-1, 0]], float); G = ep @ em; I = np.eye(2)
kr = lambda ms: reduce(np.kron, ms)
pos, neg = [], []
for k in range(7):
    pos.append(kr([G] * k + [ep] + [I] * (6 - k)))
    neg.append(kr([G] * k + [em] + [I] * (6 - k)))
gam = pos + neg
eta = [1] * 7 + [-1] * 7
for a in range(14):
    for b in range(14):
        assert np.allclose(gam[a] @ gam[b] + gam[b] @ gam[a], 2 * eta[a] * (a == b) * np.eye(128))
chi = reduce(np.matmul, gam)
assert np.allclose(chi @ chi, np.eye(128))            # real chirality: Majorana-Weyl spinors exist
w, V = np.linalg.eigh((chi + chi.T) / 2)
Sp, Sm = V[:, w > 0], V[:, w < 0]
assert Sp.shape[1] == 64
Cp, Cm = reduce(np.matmul, pos), reduce(np.matmul, neg)
sig = [(a, b) for a in range(14) for b in range(a + 1, 14)]
for name, C in [("C+", Cp), ("C-", Cm)]:
    inv = max(np.abs(C @ gam[a] @ gam[b] + (gam[a] @ gam[b]).T @ C).max() for a, b in sig)
    sym = "symmetric" if np.allclose(C, C.T) else "antisymmetric" if np.allclose(C, -C.T) else "neither"
    pp = np.abs(Sp.T @ C @ Sp).max()                   # pairing S+ with S+
    pm = np.abs(Sp.T @ C @ Sm).max()
    kin_sym = max(np.abs((Sp.T @ C @ g @ Sp) + (Sp.T @ C @ g @ Sp).T).max() for g in gam)
    kin_any = max(np.abs(Sp.T @ C @ g @ Sp).max() for g in gam)
    print(f"{name}: Spin-invariant (max defect {inv:.1e}); {sym}; |S+ x S+| = {pp:.1e}, |S+ x S-| = {pm:.1e}; "
          f"kinetic block on S+: max |M| = {kin_any:.1e}, max |M + M^T| = {kin_sym:.1e}")
print("=> for a single MW spinor the kinetic matrices C gamma^mu are antisymmetric on S+:"
      " psi^T M d psi is a total derivative for Grassmann psi (no Majorana-Weyl kinetic term in (7,7)).")
# two MW copies (same chirality): psi1^T M d psi2 - psi2^T M d psi1 is allowed (complex Weyl)

# Rarita-Schwinger kinetic term for zeta_mu in S- (vector index mu):  zeta_mu^T C gamma^{mu nu rho} d_nu zeta_rho.
# Total matrix M_{(mu,i),(rho,j)} = (C gamma^{mu nu rho})_{ij}; it survives for Grassmann zeta iff it has a nonzero
# symmetric part under the combined exchange (mu,i) <-> (rho,j).  gamma^{mu nu rho} is antisymmetric in mu<->rho.
import itertools
worst_sym, worst_any = 0.0, 0.0
for C in (Cp, Cm):
    for m, n, r in itertools.permutations(range(14), 3):
        g3 = gam[m] @ gam[n] @ gam[r]           # distinct indices: totally antisymmetric product
        A = Sm.T @ C @ g3 @ Sm                  # spinor block on S- x S-
        g3r = gam[r] @ gam[n] @ gam[m]
        B = Sm.T @ C @ g3r @ Sm                 # the (rho, mu) block
        worst_sym = max(worst_sym, np.abs(A + B.T).max())   # symmetric part of the full matrix
        worst_any = max(worst_any, np.abs(A).max())
print(f"RS kinetic block (zeta in S-): max |block| = {worst_any:.1e}, max |symmetric part of full matrix| = {worst_sym:.1e}")
print("=> the Rarita-Schwinger kinetic term of a single Majorana-Weyl zeta is also a total derivative in (7,7)"
      if worst_sym < 1e-12 else "=> RS kinetic term survives")
