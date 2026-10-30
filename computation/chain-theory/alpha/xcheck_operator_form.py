"""Coordinator cross-check of operator-origin.md, Prop. 2 (CT-G01).

Written separately from operator_origin.py. It applies the superoperators
S = (1 + (-1)^N E')/2, Pi_X = [N_Q N_q = 0], Pi_Y = [N_q = 1] directly to each of the
172 covariant records (built by symmetry.py), with N_s = 1 - tau_s and
tau_s(A) = 1_s/d_s (x) Tr_s A computed by tensor contraction. It checks that every
record is a joint eigenvector, that M = (1+(sqrt3-1) S Pi_X)(1+(sqrt2-1) S Pi_Y)
has eigenvalue sqrt(k_r) with the Letter's k_r, and that the trace reproduces
(79 + 48 sqrt2 + 42 sqrt3 + 3 sqrt6)/31104. It also counts, per class, the records
that carry a letter at a qutrit slot (the locality question of Prop. 5).
"""
import contextlib, io, importlib.util
from collections import Counter
import numpy as np
import sympy as sp

spec = importlib.util.spec_from_file_location("sym", __file__.replace("xcheck_operator_form.py", "symmetry.py"))
sym = importlib.util.module_from_spec(spec)
with contextlib.redirect_stdout(io.StringIO()):
    spec.loader.exec_module(sym)

DIMS = (3, 2, 3, 2)          # Q1, q1, Q2, q2
END = ('Q', 'q', 'Q', 'q')

def tau(A, s):
    """1_s/d_s (x) Tr_s A on B(C^3 (x) C^2 (x) C^3 (x) C^2)."""
    T = A.reshape(DIMS + DIMS)
    tr = np.trace(T, axis1=s, axis2=s + 4)                     # remove slot s
    out = np.expand_dims(np.expand_dims(tr, s), s + 4)        # re-insert identity at slot s
    out = out * np.eye(DIMS[s]).reshape([DIMS[s] if k in (s, s + 4) else 1 for k in range(8)]) / DIMS[s]
    return out.reshape(36, 36)

def occupation(A):
    """Return the eigenvalues (N_Q, N_q) of A, asserting A is a joint eigenvector.
    (Single-slot occupations are not diagonal on E'-combinations; the sums over an end type are.)"""
    out = []
    for end in ('Q', 'q'):
        occ = sum(A - tau(A, s) for s in range(4) if END[s] == end)
        c = np.vdot(A, occ) / np.vdot(A, A)
        assert np.allclose(occ, c * A, atol=1e-12)
        assert abs(c - round(c.real)) < 1e-12
        out.append(round(c.real))
    return tuple(out)

def letter_at_Q(A):
    """Does A have a nonzero component with a traceless factor at some qutrit slot?"""
    return any(np.linalg.norm(A - tau(A, s)) > 1e-12 for s in (0, 2))

rows = []
cov = [(k, A, e, c) for k, A, e, c in sym.covariant] + [('comm', D, -1, 'comm') for D in sym.comm_cov]
for kind, A, e, cat in cov:
    assert np.allclose(sym.Ep(A), e * A)
    nQ, nq = occupation(A)
    N = nQ + nq
    sigma = (-1) ** N
    S = (1 + sigma * e) // 2                      # eigenvalue of S on A
    PX = int(nQ * nq == 0)
    PY = int(nq == 1)
    m = (1 + (sp.sqrt(3) - 1) * S * PX) * (1 + (sp.sqrt(2) - 1) * S * PY)
    klass = 'comm' if cat == 'comm' else ('odd' if e == -1 else ('evmix' if cat == 'mixed' else 'evsame'))
    rows.append((klass, nQ, nq, sp.nsimplify(sp.expand(m)), letter_at_Q(A)))

print("records:", len(rows))
print("class, (N_Q, N_q), eigenvalue of M_nat, letter at a qutrit slot? -> count")
for key, n in sorted(Counter(rows).items(), key=lambda t: str(t[0])):
    print("  ", key, "->", n)

k_letter = {'odd': 1, 'evmix': sp.sqrt(2), 'evsame': sp.sqrt(3), 'comm': sp.sqrt(6)}
assert all(sp.simplify(m - k_letter[c]) == 0 for c, _, _, m, _ in rows)
total = sp.nsimplify(sum(m for *_, m, _ in rows))
print("sum of eigenvalues:", total)
assert sp.simplify(total - (79 + 48 * sp.sqrt(2) + 42 * sp.sqrt(3) + 3 * sp.sqrt(6))) == 0
print("1/value:", sp.N(31104 / total, 30))
noQ = Counter(c for c, _, _, _, q in rows if not q)
print("records with no qutrit letter, by class:", dict(noQ), " total:", sum(noQ.values()))
print("M_nat equals the Letter's sqrt(k_r) on every record: True")
