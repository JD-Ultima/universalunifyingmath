"""Numerical cross-check of the imported lemmas behind the read cap (CT-005).

For a linear map Phi: M_m -> M_n with Choi matrix C = sum_ij E_ij (x) Phi(E_ij):
  positive      <=> <x (x) y| C |x (x) y> >= 0 for all x, y   (block-positivity)
  decomposable  <=> C = P + Q^{Gamma} with P, Q >= 0          (Gamma: partial transpose)
We compute lambda* = max { lambda : C - lambda*1 = P + Q^Gamma, P, Q >= 0 } by SDP.
lambda* >= 0 (within solver tolerance) iff C is decomposable.

Checks:
  1. Choi's map on M_3 (Choi 1975), as the generalized Choi map Phi[2,0,1]:
     positive (sampled block-positivity) and NOT decomposable (lambda* < 0).
  2. Transpose map on M_3: positive, decomposable (lambda* = 0) -- solver sanity.
  3. Choi's map compressed/embedded into M_3 -> M_4 and M_4 -> M_4: still
     non-decomposable (the embedding step of the read-cap argument).
Run: python3 readcap_lemmas.py      (cvxpy 1.6.0, Clarabel/SCS solvers)
"""
import numpy as np, cvxpy as cp
rng = np.random.default_rng(20261030)

def E(n, i, j):
    m = np.zeros((n, n), complex); m[i, j] = 1; return m

def choi_matrix(phi, m):
    return sum(np.kron(E(m, i, j), phi(E(m, i, j))) for i in range(m) for j in range(m))

def choi_map(X, a=2, b=0, c=1):
    d = np.diag([a*X[0,0]+b*X[1,1]+c*X[2,2], c*X[0,0]+a*X[1,1]+b*X[2,2], b*X[0,0]+c*X[1,1]+a*X[2,2]])
    return d - X

def ptranspose(M, m, n):
    T = M.reshape(m, n, m, n)
    return T.transpose(0, 3, 2, 1).reshape(m*n, m*n)

def block_positivity(C, m, n, samples=200000):
    worst = np.inf
    X = rng.normal(size=(samples, m)) + 1j*rng.normal(size=(samples, m))
    Y = rng.normal(size=(samples, n)) + 1j*rng.normal(size=(samples, n))
    X /= np.linalg.norm(X, axis=1, keepdims=True); Y /= np.linalg.norm(Y, axis=1, keepdims=True)
    V = np.einsum('si,sj->sij', X, Y).reshape(samples, m*n)
    vals = np.real(np.einsum('si,ij,sj->s', V.conj(), C, V))
    return vals.min()

def lambda_star(C, m, n):
    N = m*n
    P = cp.Variable((N, N), hermitian=True); Q = cp.Variable((N, N), hermitian=True)
    lam = cp.Variable()
    # partial transpose of Q as a linear map: build via index permutation
    idx = np.arange(N*N).reshape(m, n, m, n).transpose(0, 3, 2, 1).reshape(N*N)
    Qv = cp.reshape(Q, (N*N, 1), order='C')
    QG = cp.reshape(Qv[idx], (N, N), order='C')
    cons = [P >> 0, Q >> 0, C - lam*np.eye(N) == P + QG]
    prob = cp.Problem(cp.Maximize(lam), cons)
    prob.solve(solver=cp.CLARABEL)
    return lam.value, prob.status

def embed(phi, m_in, m, n_out, n):
    """Phi on M_m_in -> M_n_out, extended to M_m -> M_n by compressing the input and zero-padding output."""
    def big(X):
        out = np.zeros((n, n), complex)
        out[:n_out, :n_out] = phi(X[:m_in, :m_in])
        return out
    return big

cases = [
    ("Choi map M3->M3", choi_map, 3, 3),
    ("transpose M3->M3", lambda X: X.T, 3, 3),
    ("Choi map embedded M3->M4", embed(choi_map, 3, 3, 3, 4), 3, 4),
    ("Choi map embedded M4->M4", embed(choi_map, 3, 4, 3, 4), 4, 4),
]
for name, phi, m, n in cases:
    C = choi_matrix(phi, m)
    bp = block_positivity(C, m, n)
    ls, st = lambda_star(C, m, n)
    print(f"{name:28s} min sampled <xy|C|xy> = {bp:+.6f}   lambda* = {ls:+.6f} ({st})  ->",
          "positive" if bp > -1e-9 else "NOT positive", "/",
          "decomposable" if ls > -1e-6 else "NOT decomposable")
