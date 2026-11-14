"""CZ: exact pointwise homogeneous separator for the 107D axial H6 value space."""

import argparse
from fractions import Fraction as F
import gzip
import hashlib
import json
from itertools import combinations
from math import isqrt,prod
from pathlib import Path
from time import perf_counter
import CV_h6_axial as CV
import numpy as np
import sympy as s

HERE=Path(__file__).resolve().parent
BITS=48
DYADIC=2**20
N=107


def raw_data(name):
    raw=gzip.decompress((HERE/name).read_bytes())
    return json.loads(raw),hashlib.sha256(raw).hexdigest()


def arrays():
    coords,A,D=CV.build_pencil()
    raw=gzip.decompress((HERE/"CV_h6_axial_chains.json.gz").read_bytes())
    chain=CV.read_chains(json.loads(raw),coords)
    return coords,A,D,chain


def raw_chains(coords):
    source,_=raw_data("CZ_h6_axial_local_raw_chains.json.gz")
    data=dict(character=source["ch"],chains=[dict(degree=d,entries=r) for d,r in source["chains"]])
    chains=CV.read_chains(data,coords)
    B=np.zeros((len(coords),N),dtype=np.int64)
    out=[]
    for n,c in zip(CV.DEGREES,chains):
        out.extend(c[:n+1])
    for j,c in enumerate(out):
        for k,v in c.items():
            B[k,j]=CV.finite(v)
    assert CV.BO.modular_rank(B)==N
    return chains,out


def reconstruct_basis(rows,pos,coords,expected):
    assert len(rows)==N and all(len(x)==N for x in rows)
    basis=[]
    for j,src in enumerate(rows):
        out={}
        for k,val in enumerate(src):
            c=F(val)
            if c:
                for i,v in pos[k].items():
                    out[i]=out.get(i,F(0))+c*v
        basis.append({i:v for i,v in out.items() if v})
    packed=json.dumps(dict(basis=[[[*coords[i],str(v)] for i,v in sorted(z.items())] for z in basis],character=0),
                      separators=(",",":")).encode("ascii")
    assert hashlib.sha256(packed).hexdigest()==expected
    print("Every rational source vector is reconstructed in the exact 107D polynomial coefficient span.",flush=True)
    return basis


def unit_proxies(coords,rows,proof):
    assert len(rows)==N
    points={i:CV.physical(a,m) for i,(a,m) in enumerate(coords)}
    unit=1<<BITS;grid=1<<95
    root={}
    for i,G in points.items():
        h=F(isqrt(G.denominator*grid*grid//G.numerator),grid)
        assert h*h<=1/G<(h+F(1,grid))**2
        root[i]=h
    W=np.zeros((len(coords),N),dtype=np.int64)
    err=F(0)
    for j,r in enumerate(rows):
        for i,x in r.items():
            f=x*root[i]
            w=round(unit*f)
            assert abs(w)<1<<60
            W[i,j]=w
            err+=(abs(F(w,unit)-f)+abs(x)/grid)**2
    assert err==F(proof["basis_error_squared_exact"])<F(1,10**22)
    assert sum(int(x)**2 for x in W.ravel())<121*unit*unit
    return W


def exact_core(coords,W,Y,p):
    M=CV.BV.Gradient(0,coords,Y,p)
    return (W%p).T@np.column_stack([M.core(M.matrix(x)) for x in W.T])%p


def verify_gram(coords,W,proof,full):
    Y=proof["receiver"]
    assert proof["receiver_denominator"]==DYADIC and len(set((a,m) for a,m,v in Y))==len(Y)
    possible=set(CV.BV.RECEIVERS)
    assert all((a,m) in possible and isinstance(v,int) for a,m,v in Y)
    G=np.array(proof["matrix"],dtype=object)
    assert G.shape==(N,N) and np.array_equal(G,G.T)
    t=sum(v for a,m,v in Y if m==1<<a)
    core=26*t*(W.astype(object).T@W.astype(object))-G
    assert all(int(x)%7==0 for x in core.ravel())
    core//=7
    y=sum(abs(v) for a,m,v in Y)
    h=max(sum(int(abs(v)) for v in W[:,i]) for i in range(N))
    cap=120*y*h*h
    primes=proof["primes"]
    assert all(s.isprime(p) and p%4==1 for p in primes) and len(set(primes))==len(primes)
    assert prod(primes)>2*cap and cap==int(proof["integer_bound"])
    assert max(abs(x) for x in core.ravel()) <=cap
    for p in (primes if full else primes[:1]):
        Z=exact_core(coords,W,Y,p)
        assert np.array_equal(core%p,Z)
        print("Exact 107 by 107 polarization modulo",p,flush=True)
    return G,F(16*y,DYADIC)


def definiteness(G,L,proof):
    unit=proof["matrix_denominator"]
    assert unit==21*(1<<(2*BITS+20))
    H=np.asarray(proof["factor"],dtype=object)
    scale=proof["factor_denominator"]
    assert scale==1<<36
    assert H.shape==(N,N)
    shift=F(proof["shift"])
    q=H@H.T
    numerical=F(1,10**11)
    bound=max(sum(abs(F(int(G[i,j]),unit)-F(int(q[i,j]),scale**2)-(shift if i==j else 0)) for j in range(N))
              for i in range(N))
    assert bound==F(proof["exact_defect"])
    error=L*(22*numerical+numerical*numerical)
    assert error==F(proof["exact_propagation_error"])
    assert L<30000 and shift-bound-error>F(39,100)
    assert (shift-bound-error)/121==F(proof["relative_margin"])>F(3,1000)
    print("Uniform exact receiver margin > 3/1000 on the entire real 107D point-value space.")


def main(full):
    t=perf_counter()
    proof,_=raw_data("CZ_h6_axial_local_proof.json.gz")
    coords,A,D,chain=arrays()
    CV.closure(chain,A,D)
    CV.rank_proof(coords,A,D,chain)
    alternative,Q=raw_chains(coords)
    CV.closure(alternative,A,D)
    # Any polynomial chain with these bounded degrees consists of that
    # same CV coefficient span: every evaluation lies in R, and finitely
    # many distinct evaluations recover every coefficient by Vandermonde.
    change,_=raw_data("CZ_h6_axial_local_coordinates.json.gz")
    assert change["point_basis_sha256"]==proof["basis_sha256"]
    rows=reconstruct_basis(change["rows"],Q,coords,proof["basis_sha256"])
    W=unit_proxies(coords,rows,proof)
    G,L=verify_gram(coords,W,proof,full)
    definiteness(G,L,proof)
    print(f"Runtime {perf_counter()-t:.3f} seconds.")


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--full",action="store_true",help="verify exact proxy Gram in every CRT prime")
    main(p.parse_args().full)
