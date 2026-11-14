"""CV: exact all-frequency H6 obstruction on the trivial axial sign character."""

import argparse
from fractions import Fraction as F
from functools import reduce
from itertools import combinations,product
import gzip
import hashlib
import heapq
import json
from math import comb,gcd,isqrt,prod
from pathlib import Path
import sys
from time import perf_counter

import numpy as np
import sympy as s

HERE=Path(__file__).resolve().parent
import BC_affine_hooks as BC
import BG_h6_curvature as BG
import BO_h6_full_kernel as BO
import BR_h6_an_closed as BR
import BV_h6_an_positive as BV

V=range(14)
G=tuple(BC.METRIC)
P=65537
D=43
MODULI=(65521,65537,65557,65581,65609,65617,65629,65633,65657,65677,65701)
DEGREES=[0]*2+[1]*23+[2]*13+[3]*5


def archive(name):
    raw=gzip.decompress((HERE/name).read_bytes())
    return json.loads(raw),hashlib.sha256(raw).hexdigest()


def physical(a,m):
    return G[a]*prod(G[j] for j in BC.bit_indices(m))


def build_pencil():
    coord,_=BR.stack(0)
    row=dict()
    B=dict()
    for j,(a,m) in enumerate(coord):
        w=BR.wedge(a,m)
        if w:
            row.setdefault((0,0,w[0]),{})[j]=BR.den*w[1]
        else:
            m0=m^(1<<a)
            h=F(BR.den,G[a])
            assert h.denominator==1
            row.setdefault((1,0,m0),{})[j]=h.numerator*BC.sign(1<<a,m0)
        for (b,k),c in BR.closure_column(a,m).items():
            row.setdefault((2,b,k),{})[j]=c
        w=BR.wedge(13,m)
        if w:
            B.setdefault((2,a,w[0]),{})[j]=BR.den*w[1]
    labels=sorted(row.keys()|B.keys())
    return coord,[row.get(k,{}) for k in labels],[B.get(k,{}) for k in labels]


def read_chains(data,coords):
    n=len(coords)
    loc={key:i for i,key in enumerate(coords)}
    assert data["character"]==0 and len(data["chains"])==D
    cols=[]
    for item,d in zip(data["chains"],DEGREES):
        assert item["degree"]==d
        z=[{} for _ in range(4)]
        for k,a,m,w in item["entries"]:
            assert 0<=k<=d and (a,m) in loc and loc[a,m] not in z[k]
            z[k][loc[a,m]]=F(w)
        cols.append(z)
    return cols


def closure(B,A,E):
    def row_dot(r,X):
        return sum((v*X.get(k,0) for k,v in r.items()),F(0))
    for z,d in zip(B,DEGREES):
        for a,b in zip(A,E):
            assert not row_dot(a,z[0])
            for j in range(1,d+1):
                assert not (row_dot(a,z[j])+row_dot(b,z[j-1]))
            assert not row_dot(b,z[d])
    print("All 43 rational polynomial columns satisfy the full hook/closure pencil and every endpoint.",flush=True)


def finite(c):
    x=F(c)
    return x.numerator*pow(x.denominator,-1,P)%P


def modular_sparse_kernel(rows,n):
    A=[dict(row) for row in rows]
    m=len(A)
    active=[True]*m
    heap=[(len(a),i) for i,a in enumerate(A)]
    heapq.heapify(heap)
    index=[set() for _ in range(n)]
    for i,a in enumerate(A):
        for j in a:
            index[j].add(i)
    pivots=[]
    while heap:
        d,i=heapq.heappop(heap)
        a=A[i]
        if not active[i] or len(a)!=d:
            continue
        active[i]=False
        for j in a:
            index[j].remove(i)
        if not a:
            continue
        j=min(a,key=lambda k:(len(index[k]),k))
        den=pow(a[j],-1,P)
        nrm={k:c*den%P for k,c in a.items()}
        pivots.append((j,nrm))
        for k in list(index[j]):
            b=A[k];c=b[j]
            for u,v in nrm.items():
                w=(b.get(u,0)-c*v)%P
                if w:
                    if u not in b:
                        index[u].add(k)
                    b[u]=w
                elif u in b:
                    del b[u]
                    index[u].remove(k)
            heapq.heappush(heap,(len(b),k))
    free=sorted(set(range(n))-{j for j,r in pivots})
    ans=np.zeros((n,len(free)),dtype=np.int64)
    ans[free]=np.eye(len(free),dtype=np.int64)
    for j,a in reversed(pivots):
        t=[(k,c) for k,c in a.items() if k!=j]
        if t:
            q,v=zip(*t)
            ans[j]=-np.asarray(v,dtype=np.int64)@ans[list(q)]%P
    assert all(not np.any(np.asarray(list(a.values()),dtype=np.int64)@ans[list(a)]%P) for a in rows)
    return ans


def fixed_dimension():
    # Kostant multiplicities in B2. Double all weights to remove half-integral rho.
    source=Path(BO.__file__).with_name("BO_h6_full_kernel_witnesses.json")
    data=json.loads(source.read_text(encoding="utf-8"))["extremal_vectors"]

    def part(a,b):
        return sum(b+x-y>=0 for x in range(a+1) for y in range(a-x+1)) if a>=0 else 0

    def weight(l,m,r,t):
        v=0
        for exchange in (False,True):
            C=(2*m+1,2*l+3) if exchange else (2*l+3,2*m+1)
            for a,b in product((-1,1),repeat=2):
                v+=(-1 if exchange else 1)*a*b*part((a*C[0]-(2*r+3))//2,(b*C[1]-(2*t+1))//2)
        return v

    dimension=0
    fixed=0
    for x in data:
        l,m=x["weight"]; copies=len(x["vectors"])
        tab={(a,b):weight(l,m,a,b) for a in range(-l,l+1) for b in range(-l,l+1)}
        assert min(tab.values())>=0
        z=sum(tab.values())
        assert z==(l-m+1)*(l+m+2)*(2*l+3)*(2*m+1)//6
        n2=sum((1-2*(a%2))*v for (a,b),v in tab.items())
        n4=sum((1-2*((a+b)%2))*v for (a,b),v in tab.items())
        assert (z+10*n2+5*n4)%16==0
        fixed+=copies*((z+10*n2+5*n4)//16)
        dimension+=copies*z
    assert (dimension,fixed)==(6834,389)
    return fixed


def rank_proof(coords,A,E,B):
    R=BG.curvature_model()
    cov=[[(i,j,b,finite(-v)) for (i,j),matrix in R.items() for b,v in enumerate(matrix[a]) if v] for a in V]
    rows={}
    for i,(a,m) in enumerate(coords):
        if (w:=BR.wedge(a,m)):
            rows.setdefault((0,0,w[0]),{})[i]=w[1]%P
        else:
            w=m^(1<<a)
            rows.setdefault((1,0,w),{})[i]=finite(F(BC.sign(1<<a,w),G[a]))
        for j,k,b,v in cov[a]:
            t=(1<<j)|(1<<k)
            if not t&m:
                coefficient=v*BC.sign(t,m)%P
                if coefficient:
                    rows.setdefault((2,b,t|m),{})[i]=coefficient
    K=modular_sparse_kernel(list(rows.values()),len(coords))
    assert K.shape==(2672,fixed_dimension())

    def dense(source,n):
        d=np.zeros((len(source),n),dtype=np.int64)
        for r,row in enumerate(source):
            for i,x in row.items():
                d[r,i]=x%P
        return d

    A,E=dense(A,len(coords)),dense(E,len(coords))
    X=np.zeros((len(coords),107),dtype=np.int64)
    I=[]
    j=0
    for z,d in zip(B,DEGREES):
        for k in range(d+1):
            for i,c in z[k].items():
                X[i,j+k]=finite(c)
        I.extend(range(j,j+d));j+=d+1
    assert j==107 and len(I)==64
    T=E@X[:,I]%P
    direct=BO.modular_rank(X.copy()),BO.modular_rank(T.copy()),BO.modular_rank(np.column_stack((T,A@K%P,E@K%P)))
    assert direct==(107,64,628), direct
    print("Modular coefficient/image/joint ranks 107,64,628; exact sign-fixed K dimension 389.",flush=True)


def dyadic_proxies(B,coords,precision):
    unit=1<<precision; grid=1<<90
    out=np.zeros((len(coords),D,4),dtype=np.int64)
    roots=[F(isqrt(q.denominator*grid*grid//q.numerator),grid) for a,m in coords if (q:=physical(a,m))]
    dist=[]
    for i,z in enumerate(B):
        loss=F(0)
        for p,x in enumerate(z):
            for k,v in x.items():
                r=roots[k]
                integer=round(unit*v*r)
                assert abs(integer)<1<<60
                out[k,i,p]=integer
                loss+=(abs(F(integer,unit)-v*r)+abs(v)/grid)**2
        dist.append(loss)
    for k in range(4):
        total=sum(int(z)**2 for z in out[:,:,k].ravel())
        assert total<49*unit*unit
    return out,dist


def gram_CRT(orth,receiver):
    cols,Y,bits=receiver
    n0=orth.shape[1]
    W=orth.reshape(len(cols),-1)
    h=max(sum(abs(int(v)) for v in W[:,j]) for j in range(W.shape[1]))
    y=sum(abs(v) for a,m,v in Y)
    bound=120*h*h*y
    assert prod(MODULI)>2*bound>0
    H=None
    mult=1
    for p in MODULI:
        assert p%4==1 and s.isprime(p)
        M=BV.Gradient(0,cols,Y,p)
        transformed=np.column_stack([M.core(M.matrix(x)) for x in W.T])
        out=((W%p).T@transformed)%p
        H=out.astype(object) if H is None else H+mult*(((out-H%p)*pow(mult,-1,p))%p)
        mult*=p
        print("Exact proxy Gram modulo",p,flush=True)
    H=np.where(H>mult//2,H-mult,H)
    assert np.array_equal(H,H.T) and max(abs(c) for c in H.ravel())<=bound
    q=sum(v for a,m,v in Y if m==1<<a)
    matrix=(-7*H+26*q*(W.astype(object).T@W.astype(object))).reshape(n0,4,n0,4)
    gram_denom=21*2**(bits*2+20)
    result=[]
    for k in range(7):
        a=np.zeros((D,D),dtype=object);b=a.copy();c=a.copy();e=a.copy()
        for i,j in product(range(D),repeat=2):
            for p in range(4):
                r=k-p
                if not 0<=r<4:
                    continue
                plane=(r-p)%4
                v=int(matrix[i,p,j,r])
                if plane%2:
                    b[i,j]+=(1 if plane==1 else -1)*v
                else:
                    a[i,j]+=(1 if plane==0 else -1)*v
                a0,a1=DEGREES[i]-p,DEGREES[j]-r
                if a0<0 or a1<0:
                    continue
                v=int(matrix[i,a0,j,a1])
                if plane%2:
                    e[i,j]+=(1 if plane==1 else -1)*v
                else:
                    c[i,j]+=(1 if plane==0 else -1)*v
        result.append(dict(power=k,real=a.tolist(),imaginary=b.tolist(),revreal=c.tolist(),revimaginary=e.tolist()))
    return result,gram_denom,bound


def Bernstein(M,lo,hi):
    n=len(M)-1
    coeff=[]
    for i in range(n+1):
        coeff.append(sum((M[k]*F(comb(k,i))*lo**(k-i)*(hi-lo)**i for k in range(i,n+1)),np.zeros_like(M[0])))
    return [sum((coeff[j]*F(comb(i,j),comb(n,j)) for j in range(i+1)),np.zeros_like(M[0])) for i in range(n+1)]


def hermitian_plot(polynomial,tag,denom):
    out=[]
    for p in polynomial:
        R=np.asarray(p[f"{tag}real"],dtype=object)
        I=np.asarray(p[f"{tag}imaginary"],dtype=object)
        assert np.array_equal(R,R.T) and np.array_equal(I,-I.T)
        out.append(np.block([[R,-I],[I,R]])/F(denom))
    return np.array(out)


def exact_positive(M,shift,margin):
    # The floating Cholesky only proposes a dyadic factor. The exact
    # rational max-row defect certifies the bound independently.
    q=np.linalg.cholesky(np.asarray(M,dtype=float)-float(shift)*np.eye(M.shape[0]))
    d=1<<32
    proposed=np.rint(d*q).astype(np.int64).astype(object)
    certified=proposed@proposed.T
    loss=max(sum(abs(v-F(c,d*d)-(shift if i==j else 0)) for j,(v,c) in enumerate(zip(row,b)))
             for i,(row,b) in enumerate(zip(M,certified)))
    assert loss<shift-margin
    return shift-loss-margin


def certify(gram,dist,Y,record):
    U=F(16)*sum(abs(F(v,2**20)) for a,m,v in Y)
    assert U < 60000
    M=hermitian_plot(gram["polynomial_matrices"],"",gram["gram_denominator"])
    N=hermitian_plot(gram["polynomial_matrices"],"rev",gram["gram_denominator"])
    d=sum(dist,F(0));k=10**16
    e=F(isqrt(d.numerator*k*k//d.denominator)+1,k)
    assert e*e>d and e<F(1,10**10)
    minimum=F(1)
    for i,raw_shift in enumerate(record["shifts"]):
        if i<42:
            lo,hi=F(i,2),F(i+1,2)
            c=[1/(1+hi)**d for d in DEGREES]*2
            G0=M*np.outer(c,c)[None,...]
            K=F(4)
        else:
            assert i==42
            lo,hi=F(0),F(1,21)
            G0=N
            K=sum(hi**j for j in range(4))
        uncertainty=F(60000)*K*K*(14*e+e*e)
        for B in Bernstein(G0,lo,hi):
            minimum=min(minimum,exact_positive(B,F(raw_shift),uncertainty))
    assert len(record["shifts"])==43 and minimum>F(1,10**6)
    assert 4*7+4*e<29
    print("Exact rational positivity in all 301 Bernstein matrices, minimum Gram bound",float(minimum),flush=True)
    print("Uniform normalized receiver margin > 1/10^9 for every real axial frequency.",flush=True)


def run(full):
    t=perf_counter()
    record=json.loads((HERE/"CV_h6_axial_certificate.json").read_text(encoding="utf-8"))
    chain,h=archive("CV_h6_axial_chains.json.gz")
    assert h==record["chain_sha256"]
    result,h=archive("CV_h6_axial_proxy_gram.json.gz")
    assert h==record["proxy_gram_sha256"]
    coords,A,E=build_pencil()
    B=read_chains(chain,coords)
    closure(B,A,E)
    rank_proof(coords,A,E,B)
    bits=result["dyadic_inputs_bits"]
    proxy,dist=dyadic_proxies(B,coords,bits)
    assert dist==list(map(F,result["operator_input_errors"]))
    encoded=json.loads((HERE/"CV_h6_axial_receiver.json").read_text(encoding="utf-8"))
    Y=[tuple(x) for x in encoded["rows"]]
    assert encoded["denominator"]==2**20 and len(Y)==len(set((a,m) for a,m,v in Y))
    labels=set(BV.RECEIVERS)
    assert all((a,m) in labels for a,m,v in Y)
    if full:
        polys,g,bound=gram_CRT(proxy,(coords,Y,bits))
        assert polys==result["polynomial_matrices"] and g==result["gram_denominator"]
        assert bound==int(result["integer_coefficient_bound"])
    certify(result,dist,Y,record)
    print(f"Runtime {perf_counter()-t:.3f} seconds.",flush=True)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--full",action="store_true",help="also rebuild the exact Gaussian-integer proxy Gram in all eleven primes")
    run(p.parse_args().full)
