"""Exact classification of real characteristic rays of the three DK gradients."""

from collections import deque
from fractions import Fraction as F
import json
from math import factorial, isqrt, lcm
from pathlib import Path
from time import perf_counter

import numpy as np
from scipy.linalg import cholesky, solve_triangular
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components

import DK_h6_characteristic as DK

HERE = Path(__file__).resolve().parent


def gram_parts(r,g):
    h=(2,6,12,20)
    return [[h[a]*g[:,None]*r[a*286:(a+1)*286,b*286:(b+1)*286] for b in range(4)] for a in range(4)]


def combine(P,x,y):
    return sum((x[i]*y[j]*P[i][j] for i in range(4) for j in range(4) if x[i] and y[j]),
               np.zeros(P[0][0].shape,dtype=object))


def exact_pos_scaled(C,D,roundbits=36):
    # Multiply by diagonal inverse square-root bounds to remove the
    # factorial weights before exact lower bound verification.
    n=len(D)
    u=1<<roundbits
    scale=np.asarray([isqrt(int(u*u*F(1,int(z)))) for z in D],dtype=object)
    assert np.all(scale>0)
    M=np.asarray(C,dtype=object)
    assert np.array_equal(M,M.T)
    B=scale[:,None]*M*scale[None,:]
    diag=max(max(abs(int(v)) for v in s) for s in B)
    Bf=np.asarray(B,dtype=float)/diag
    Lf=cholesky(Bf,lower=True)
    def integers(array,bits):
        src=np.ldexp(array,bits)
        assert np.max(abs(src))<1<<62
        return np.asarray(np.rint(src),dtype=np.int64).astype(object)
    L=integers(Lf,34)
    Q=L@L.T
    res=B*(1<<68)-diag*Q
    err=max(sum(abs(v) for v in row) for row in res)
    # Gershgorin on B is insufficient. To show approximate factor has safe
    # smallest singular use its approximate inverse and an exact norm deficit.
    inv=solve_triangular(Lf,np.eye(n),lower=True,check_finite=False)
    I=integers(inv,27)
    W=(1<<61)*np.eye(n,dtype=object)-L@I
    def l1(W):
        return max(sum(abs(int(v)) for v in row) for row in W)
    wi,w1=l1(W),l1(W.T)
    assert 4*wi*w1 < (1<<122)
    ii,ij=l1(I),l1(I.T)
    assert 4*int(err)*ii*ij < int(diag)*(1<<(68+54))
    return F(int(err),diag*(1<<68))


def exact_local_and_lipschitz(raw):
    r,q,monos,g=raw
    C=gram_parts(r,g)
    G=np.diag(g)
    h=[2,6,12,20]
    block=np.asarray(r,dtype=object)*(np.repeat(h,286)*np.tile(g,4))[:,None]
    assert np.array_equal(block,block.T)
    major=np.diag(np.repeat(h,286)*np.tile(g,4)).astype(object)
    # Prove the norm bound ||Pbad(x⊗b)||²≤(5/8)||x⊗b||²
    # by certifying all connected Cartan-fiber projector blocks.
    F5=5*q*major-8*block
    n,L=connected_components(csr_matrix(F5!=0),directed=False)
    proof=[]
    for i in range(n):
        use=np.flatnonzero(L==i)
        V=exact_pos_scaled(F5[np.ix_(use,use)],major.diagonal()[use])
        proof.append((len(use),float(V)))
    assert sum(k for k,v in proof)==1144
    print("Exact Cartan projection norm bound squared < 5/8.",flush=True)

    # Unsorted projective center ξ0=diag(1,1,-1,-1,0)/2 has rational
    # basis coordinates (0,1/3,1/6,0).  Use primitive orthogonal axes
    # e1, (0,1,-1,0), e4, with respective squared norms 2,18,20.
    x=[0,2,1,0]
    T=[[1,0,0,0],[0,1,-1,0],[0,0,0,1]]
    tmetric=[2,18,20]
    K=[j for j,p in enumerate(monos) if max(p)==10]
    E=[j for j in range(286) if j not in K]
    A0=combine(C,x,x)
    assert not np.any(A0[:,K])
    A=A0[np.ix_(E,E)]
    assert len(K)==4
    # Base gap σ²>1/16 at unit ξ0. x is six times that covector.
    e0=exact_pos_scaled(16*A-36*q*G[np.ix_(E,E)],g[E])
    H=np.vstack([np.column_stack([combine(C,u,v)[np.ix_(K,K)] for v in T]) for u in T])
    B=np.column_stack([combine(C,x,u)[np.ix_(E,K)] for u in T])
    # Schur lower > 9/25 for the full tensor product of the three
    # orthogonal transverse directions with the characteristic kernel.
    Gnew=(np.asarray(tmetric,dtype=object)[:,None]*g[K][None,:]).ravel()
    J=np.block([[25*A,25*B],[25*B.T,25*H-9*q*np.diag(Gnew)]])
    exact_pos_scaled(J,np.r_[g[E],Gnew])
    # sqrt(5/8) < 4/5 and the near-ball Euclidean radius is 1/25.
    # The image equation gives
    # ||e|| <= c*r / ((1-r)/4 - c*r) ||k||.
    # The quotient equation is incompatible when c times this ratio <3/5.
    c,radius=F(4,5),F(1,25)
    d=(1-radius)/4-c*radius
    assert d>0 and c*c*radius < F(3,5)*d
    print("Exact local rank gap >1/4 and tensor Schur lower >3/5; base defect",float(e0),flush=True)


def verify_cover(raw,proof):
    r,q,monos,g=raw
    C=gram_parts(r,g)
    T=np.array([[1]*k+[-k]+[0]*(4-k) for k in range(1,5)],dtype=np.int64)
    m=np.array([2,6,12,20],dtype=np.int64)
    h=np.array([[(5-k) if j<k else -k for j in range(5)] for k in (1,2,3,4)],dtype=np.int64)
    H=h@h.T
    cl=[[(F(int((h@T.T)[i,j]),5*int(m[j]))) for j in range(4)] for i in range(4)]
    p=np.array([0,F(1,2),F(1,2),0],dtype=object)
    v=np.eye(4,dtype=object)
    todo=deque([np.vstack((p,v[0],v[i],v[3])) for i in (1,2)])
    n=0
    near=0
    exterior=0
    t=perf_counter()

    def sqnorm(z):
        return F(z@H@z,25)

    def transformed(z):
        a=z@np.array(cl,dtype=object)
        d=lcm(*(F(c).denominator for c in a))
        b=[int(c*d) for c in a]
        return d,b

    while todo:
        W=todo.popleft()
        decision=proof[n]
        n+=1
        if isinstance(decision,list):
            a,b=decision
            assert 0<=a<b<4
            X=W.copy();Y=W.copy()
            mid=np.array([F(z,2) for z in W[a]+W[b]],dtype=object)
            X[b]=Y[a]=mid
            todo.extend((X,Y))
        elif decision==1:
            assert max(sqnorm(z-p) for z in W)<F(1,625)
            near+=1
        else:
            assert decision==2
            exterior+=1
            z=np.array([F(c,4) for c in sum(W)],dtype=object)
            radius=max(sqnorm(v-z) for v in W)
            d,a=transformed(z)
            u=radius*F(5,8)*d*d*q
            M=u.denominator*combine(C,a,a)-u.numerator*np.diag(g.astype(object))
            nums,labels=connected_components(csr_matrix(M!=0),directed=False)
            for i in range(nums):
                ix=np.flatnonzero(labels==i)
                exact_pos_scaled(M[np.ix_(ix,ix)],g[ix])

    assert n==len(proof)==700 and (near,exterior)==(27,324)
    print("Exact Weyl cover:",n,"nodes,",near,"near the characteristic ray,",exterior,
          "positive elsewhere.",round(perf_counter()-t,3),"seconds.",flush=True)


def main():
    t=perf_counter()
    R,I,q,monomials=DK.orthogonal_projection_columns()
    assert not np.any(I)
    g=np.array([np.prod([factorial(k) for k in a]) for a in monomials],dtype=np.int64)
    raw=(R,q,monomials,g)
    exact_local_and_lipschitz(raw)
    decisions=json.loads((HERE/"DK_h6_characteristic_cover.json").read_text())["decisions"]
    verify_cover(raw,decisions)
    print(f"Exact real projective classification completed in {perf_counter()-t:.3f}s.")


if __name__=="__main__":
    main()
