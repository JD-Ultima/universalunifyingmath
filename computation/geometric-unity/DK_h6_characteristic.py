"""Exact real characteristic and transverse rank for the three DK gradients."""

from itertools import combinations
from math import lcm
from time import perf_counter

import numpy as np
from scipy.sparse import bmat, csr_matrix, eye

P = 65537
IMAGINARY = 256
D = 286
SHIFTS = (20, 4, -6, -10, -20, -28)


def spin_generators():
    I = np.eye(2, dtype=np.complex128)
    x = np.array([[0, 1], [1, 0]])
    y = np.array([[0, -1j], [1j, 0]])
    z = np.diag([1, -1])
    gamma = [np.kron(a, b) for a, b in [(x, I), (y, I), (z, x), (z, y), (z, z)]]
    for i, j in combinations(range(5), 2):
        assert not np.any(gamma[i] @ gamma[j] + gamma[j] @ gamma[i])
    assert all(np.array_equal(g @ g, np.eye(4)) for g in gamma)
    monomials = [(a, b, c, 10-a-b-c) for a in range(11)
                 for b in range(11-a) for c in range(11-a-b)]
    assert len(monomials) == D
    where = {a: i for i, a in enumerate(monomials)}
    gen = []
    for i, j in combinations(range(5), 2):
        q = gamma[i] @ gamma[j]
        A = np.zeros((D, D), dtype=np.complex128)
        for a, b in zip(*np.nonzero(q)):
            for n, alpha in enumerate(monomials):
                if alpha[b]:
                    v = list(alpha)
                    v[b] -= 1
                    v[a] += 1
                    A[where[tuple(v)], n] += alpha[b]*q[a, b]
        gen.append(csr_matrix(A))
    return gen, monomials


def tangent_generators():
    vectors = []
    for i, j in combinations(range(5), 2):
        q = np.zeros((5, 5), dtype=np.int64)
        q[i, j] = q[j, i] = 1
        vectors.append(q)
    vectors.extend(np.diag([1]*k + [-k] + [0]*(4-k)) for k in range(1, 5))
    T = np.array(vectors)
    inner = T.reshape(14, -1)
    g = np.sum(inner*inner, axis=1)
    assert np.array_equal(inner @ inner.T, np.diag(g))
    scale = lcm(*map(int, g))
    gen = []
    for i, j in combinations(range(5), 2):
        Q = np.zeros((5, 5), dtype=np.int64)
        Q[i, j], Q[j, i] = 1, -1
        raw = scale * (inner @ (Q @ T - T @ Q).reshape(14, -1).T)
        assert not np.any(raw % g[:, None])
        gen.append(raw // g[:, None])
    return scale, gen


def polynomial_symbol():
    """Return exact Gaussian integer matrices for a common scaled cubic."""
    ga, monomials = spin_generators()
    s, va = tangent_generators()
    blocks = {}
    for V, J in zip(va, ga):
        for a, b in zip(*np.nonzero(V)):
            blocks[a,b] = blocks.get((a,b), 0) - V[a,b]*J
    M = bmat([[blocks.get((a,b)) for b in range(14)] for a in range(14)], format="csr")
    assert M.shape == (14*D, 14*D)
    assert np.all(M.data == np.rint(M.data.real) + 1j*np.rint(M.data.imag))
    # Check a sufficient exact-integer bound at each of the three sparse
    # products. Use the complex integer coefficient l1 row sum to bound every
    # component product and every partial sum under ordinary binary64.
    ab = abs(M.real) + abs(M.imag)
    m = int(max(ab.sum(axis=1).A.ravel()))
    need = eye(14*D, dtype=np.complex128, format="csr")[:, 10*D:]
    for u in SHIFTS[:3]:
        f = int(np.max(abs(need.data.real) + abs(need.data.imag), initial=0))
        assert (m + abs(u*s))*f < 1 << 53
        need = M @ need - u*s*need
        assert np.all(need.data == np.rint(need.data.real) + 1j*np.rint(need.data.imag))
    return need.toarray().reshape(14*D, 4, D), monomials


def orthogonal_projection_columns():
    """Exact Gaussian-integer numerators of the Cartan projector columns."""
    from fractions import Fraction as F
    from math import prod

    gs, monos = spin_generators()
    s, vg = tangent_generators()
    blocks = {}
    for V,J in zip(vg,gs):
        for a,b in zip(*np.nonzero(V)):
            blocks[a,b]=blocks.get((a,b),0)-V[a,b]*J
    M=bmat([[blocks.get((a,b)) for b in range(14)] for a in range(14)],format="csr")
    column=eye(14*D,dtype=np.complex128,format="csr")[:,10*D:]
    # Interpolate the spectral projector on all three forbidden eigenvalues.
    coefficients=[F(0)]*6
    for r in SHIFTS[3:]:
        others=[u for u in SHIFTS if u!=r]
        q=[1]
        for u in others:
            q=[-u*q[0], *(q[j-1]-u*q[j] for j in range(1,len(q))), q[-1]]
        coefficients=[a+F(b,prod(r-u for u in others)) for a,b in zip(coefficients,q)]
    assert all(sum(q*u**k for k,q in enumerate(coefficients))==int(j>=3)
               for j,u in enumerate(SHIFTS))
    # Six Newton coefficients and roots give the same projection with
    # smaller intermediate products than dense powers in the monomial basis.
    nodes = SHIFTS
    newton = [F(int(j>=3)) for j in range(6)]
    for j in range(1,6):
        for k in range(5,j-1,-1):
            newton[k]=(newton[k]-newton[k-1])/(nodes[k]-nodes[k-j])
    denom=lcm(*(q.denominator*s**k for k,q in enumerate(newton)))
    # Store each power as an exact Gaussian integer by an a priori row bound.
    m=int(max((abs(M.real)+abs(M.imag)).sum(axis=1).A.ravel()))
    num=np.zeros((4*D,4*D),dtype=object)
    imag=num.copy()
    for k,v in enumerate(newton):
        c=F(denom,s**k)*v
        assert c.denominator==1
        sub=column[10*D:].tocoo()
        for i,j,a in zip(sub.row,sub.col,sub.data):
            num[i,j]+=int(a.real)*c.numerator
            imag[i,j]+=int(a.imag)*c.numerator
        if k<5:
            root=nodes[k]*s
            # Split each Gaussian coefficient into low and high 24-bit
            # integer pieces. Each multiplication and their signed sum are
            # still bounded exactly; the checked output itself stays <2^53.
            bound=int(np.max(abs(column.data.real)+abs(column.data.imag),initial=0))
            assert bound<1<<53
            if (m+abs(root))*bound >= 1<<53:
                split=1<<24
                high=column.copy()
                high.data=np.floor(column.data.real/split)+1j*np.floor(column.data.imag/split)
                lo=column-split*high
                assert (m+abs(root))*max(split*2,bound//split+4)<1<<53
                low_image=M@lo-root*lo
                high_image=M@high-root*high
                high_image,low_image=high_image.tocoo(),low_image.tocoo()
                result={(int(i),int(j)):complex(a) for i,j,a in zip(low_image.row,low_image.col,low_image.data)}
                for i,j,a in zip(high_image.row,high_image.col,high_image.data):
                    old=result.get((int(i),int(j)),0)
                    re=int(old.real)+split*int(a.real)
                    im=int(old.imag)+split*int(a.imag)
                    assert max(abs(re),abs(im))<1<<53
                    result[int(i),int(j)]=complex(re,im)
                i,j=zip(*result)
                column=csr_matrix((list(result.values()),(i,j)),shape=column.shape)
            else:
                column=M@column-root*column
            assert np.all(column.data==np.rint(column.data.real)+1j*np.rint(column.data.imag))
    return num,imag,denom,monos


def rank_mod(A):
    A = A.copy() % P
    m,n = A.shape
    r = 0
    for j in range(n):
        found = np.flatnonzero(A[r:,j])
        if len(found):
            k = r + int(found[0])
            A[[r,k]] = A[[k,r]]
            choices = np.flatnonzero(A[r+1:,j]) + r + 1
            e = (A[r,j:]*pow(int(A[r,j]),-1,P)) % P
            A[choices,j:] = (A[choices,j:] - A[choices,j,None]*e) % P
            r += 1
        if r == m:
            break
    return r


def verify():
    t = perf_counter()
    E, monomials = polynomial_symbol()
    # Cartan basis rows are diag(1,-1,0,0,0), diag(1,1,-2,0,0),
    # diag(1,1,1,-3,0), and diag(1,1,1,1,-4).
    # 4 times the second + 2 times the third gives 6*(1,1,-1,-1,0).
    A = 4*E[:,1,:] + 2*E[:,2,:]
    pure = [j for j, p in enumerate(monomials) if max(p) == 10]
    assert len(pure) == 4 and not np.any(A[:,pure])
    remaining = [j for j in range(D) if j not in pure]
    # e1,e2,e4 are linearly independent modulo the nonzero axial direction.
    # Their components on all four axial kernel directions are tested at once.
    joined = np.column_stack((A[:,remaining], *(E[:,a,pure] for a in (0,1,3))))
    assert IMAGINARY**2 % P == P-1 and joined.shape == (14*D,294)
    assert max(np.max(abs(joined.real)), np.max(abs(joined.imag))) * IMAGINARY < 1 << 62
    assert rank_mod((joined.real.astype(np.int64) + IMAGINARY*joined.imag.astype(np.int64)) % P) == 294
    print("Exact joint characteristic at diagonal (1,1,-1,-1,0), dimension 4, transverse rank 294.")
    print(f"Runtime {perf_counter()-t:.3f} seconds.")


if __name__ == "__main__":
    verify()
