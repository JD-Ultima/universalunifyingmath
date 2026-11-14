"""Verify the CR compact-orbit first-jet test on the complete BS box."""

import argparse
from fractions import Fraction as F
import gzip
import hashlib
import json
from math import isqrt, lcm, prod
from pathlib import Path
from time import perf_counter

HERE = Path(__file__).resolve().parent

import BC_affine_hooks as BC
import BR_h6_an_closed as BR
import BG_h6_curvature as BG
import BO_h6_full_kernel as BO

G=tuple(BC.METRIC)
V=range(14)
D=140


def bit(m, i):
    return (m>>i)&1


def blade_gram(m):
    return prod(G[i] for i in BC.bit_indices(m))


def orbit_generators(t):
    # Six BC SO4 rotations have unit plane convention; the four
    # rotations in a plane containing its fifth axis need a factor of two.
    return [{key:c*(1 if a<6 else 2) for key,c in BC.action(a,t).items()} for a in range(10)]


def rank_modular(rows, derivatives, p):
    import numpy as np

    length=len(rows)
    assert length==D and BO.PRIME==p
    A=np.zeros((D,D),dtype=np.int64)
    for r,(a,m) in enumerate(rows):
        assert m.bit_count()==7 and a in V
        for i in BC.bit_indices(m):
            src=m^(1<<i)
            s=BC.sign(1<<i,src)
            for j,z in enumerate(derivatives):
                # Branch also verifies covariant slot and exterior sign order.
                u=F(z.get((a,src),0),1)*s
                A[r,i*10+j]=u.numerator*pow(u.denominator,-1,p)%p
    determinant=1
    for j in range(D):
        hits=np.flatnonzero(A[j:,j])
        assert len(hits)
        i=j+int(hits[0])
        if j!=i:
            A[[i,j]]=A[[j,i]]
            determinant=(-determinant)%p
        c=int(A[j,j])
        determinant=determinant*c%p
        A[j,j:]=A[j,j:]*pow(c,-1,p)%p
        A[j+1:,j:]=(A[j+1:,j:]-A[j+1:,j,None]*A[j,None,j:])%p
    return determinant


def exact_gram(t, derivatives):
    # Accumulate integer numerator weights instead of floating Gram entries.
    T={}
    for j,z in enumerate(derivatives):
        for (a,m),c in z.items():
            for i in V:
                if not bit(m,i):
                    target=m|(1<<i)
                    T.setdefault((a,target),{})[i*10+j]=BC.sign(1<<i,m)*c
    weights={key:G[key[0]]*blade_gram(key[1]) for key in T}
    g=lcm(*(c.numerator for c in weights.values()))
    q=lcm(*(c.denominator for z in derivatives for c in z.values()))
    acc=[[0]*D for _ in range(D)]
    for key,row in T.items():
        k=F(g,weights[key])
        assert k.denominator==1
        terms=sorted((i,int(c*q)) for i,c in row.items())
        for i,(j,c) in enumerate(terms):
            u=c*k.numerator
            for b,z in terms[i:]:
                acc[j][b]+=u*z
    return [[F(acc[min(i,j)][max(i,j)],g*q*q) for j in range(D)] for i in range(D)]


def check_positive(matrix, lower):
    U=[x.copy() for x in matrix]
    for i in range(D):
        U[i][i]-=lower/G[i//10]
    if __debug__:
        assert all(U[i][j]==U[j][i] for i in range(D) for j in range(i))
    for i in range(D):
        d=U[i][i]
        assert d>0
        for j in range(i+1,D):
            if U[i][j]:
                v=U[i][j]/d
                for k in range(j,D):
                    U[j][k]-=v*U[i][k]
    print("All 140 pivots of the exact shifted Gram are positive.")


def verify_generators():
    # Compare the entire tangent covariant action for the canonical unit
    # physical plane rotations with BC's documented generator normalization.
    from itertools import combinations

    Q=BG.integer_phi()
    for j,(a,b) in enumerate(combinations(range(5),2)):
        plane=tuple(tuple(int((i,k)==(a,b))-int((i,k)==(b,a)) for k in range(5)) for i in range(5))
        action=[]
        for ext in V:
            x=tuple(-F(BG.trace_product(BG.bracket(plane,q),Q[ext]),25*G[ext]) for q in Q)
            action.append(x)
        candidates=[]
        for k,rs in enumerate(BC.COV):
            multiplier=2 if k>=6 else 1
            source=[tuple(multiplier*dict(r).get(i,0) for i in V) for r in rs]
            if source==action or source==[tuple(-v for v in r) for r in action]:
                candidates.append(k)
        assert len(candidates)==1
        yield candidates[0]


def counterexample(R):
    data=json.loads((HERE/"CL_h6_carrier_witness.json").read_text(encoding="utf-8"))
    source=data["fields"][0]
    assert source["character"]==24
    b={(a,m):F(c) for a,m,c in source["entries"]}
    assert BC.hook_traces(b)==({}, {}) and not BG.full_sparse_curvature(R,b) and not BR.check_sparse(b)
    J=BC.action(0,b)
    assert J[(4,8433)]==-48 and J[(7,10642)]==1368
    assert all(bit(m,13) for a,m in J)
    assert BC.hook_traces(J)==({}, {}) and not BG.full_sparse_curvature(R,J)
    print("CL value violates universal orbit-jet injectivity: (J_01 b)_(4,8433)=-48, e13 wedge J_01 b=0.")


def main(full):
    t=perf_counter()
    record=json.loads((HERE/"CR_h6_orbit_certificate.json").read_text(encoding="utf-8"))
    raw=gzip.decompress((HERE/"partial_fields.json.gz").read_bytes())
    assert hashlib.sha256(raw).hexdigest()==record["anchor_hash"]
    data=json.loads(raw)
    assert tuple(map(F,data["metric_gram"]))==G
    fields=[{(a,m):F(x,record["denominator"]) for a,m,x in record["entries"]} for record in data["fields"]]
    assert len(fields)==9
    assert sorted(verify_generators())==list(range(10))
    b=fields[0]
    R=BG.curvature_model()
    assert BC.hook_traces(b)==({}, {}) and not BG.full_sparse_curvature(R,b)
    derivatives=orbit_generators(b)
    # Selected minor was obtained in the original BC generator convention,
    # whereas shifted positivity uses the 4 doubled unit generators.
    scaled=[{key:F(c,1 if j<6 else 2) for key,c in z.items()} for j,z in enumerate(derivatives)]
    determinant=rank_modular(record["selected_rows"],scaled,record["prime"])
    assert determinant==record["determinant"]!=0
    print("Nonzero exact 140x140 rational minor:",determinant,"mod",record["prime"])
    counterexample(R)
    radius=F(record["box_radius"])
    assert radius==F(1,10**8)
    ns=[sum((c*c/(G[a]*blade_gram(m)) for (a,m),c in z.items()),F(0)) for z in fields[1:]]
    assert ns==list(map(F,record["axis_norms_squared"]))
    integers=10**6
    roots=[F(isqrt(z.numerator*integers**2//z.denominator)+1,integers) for z in ns]
    assert all(r*r>z for r,z in zip(roots,ns))
    perturb=120*radius*sum(roots,F(0))
    assert perturb==F(record["perturbation_bound"])
    assert F(1,2)>(F(3,5)+perturb)**2
    if full:
        check_positive(exact_gram(b,derivatives),F(1,2))
        print("Uniform lower bound sigma_min > 3/5 throughout the complete BS box.")
    else:
        print("Use --full to independently check the shifted Gram for the entire box.")
    print(f"Runtime {perf_counter()-t:.3f} seconds.")


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--full",action="store_true",help="also reconstruct the complete rational shifted Gram")
    main(p.parse_args().full)
