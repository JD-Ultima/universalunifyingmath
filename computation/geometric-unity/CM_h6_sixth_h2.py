"""Exact H2 projection guard, including verification of the 14-vector module."""
from functools import lru_cache
from fractions import Fraction as Q
import json
import CF_h6_fourth_jet as v
BO=v.BO;HERE=v.HERE

@lru_cache(None)
def sharp_form(mask):
    out={0:Q(1)}
    for i in BO.axes(mask):
        nxt={}
        for m,c in out.items():
            for j in range(14):
                q=Q(v.inverse[i][j])
                if q and not m>>j&1:BO.add(nxt,m|(1<<j),c*q*(-1 if (m>>(j+1)).bit_count()%2 else 1))
        out=nxt
    return out

def sharp(tensor):
    out={}
    for (a,m),c in tensor.items():
        for b in range(14):
            q=Q(v.inverse[a][b])
            if q:
                for n,z in sharp_form(m).items():BO.add(out,(b,n),c*q*z)
    return out

def build():
    import sympy as s
    import gzip
    data=json.loads(gzip.decompress((HERE/'CM_h6_sixth_aux.json.gz').read_bytes()))['h2_module'];blocks={};basis=[]
    for block in data['blocks']:
        wt=tuple(block['weight']);vectors=[]
        for entries in block['vectors']:
            t={(a,m):Q(c) for a,m,c in entries};p=min(t);t={k:c/t[p] for k,c in t.items()}
            assert all(tuple(v.weights[a][z]+sum(v.weights[i][z] for i in BO.axes(m)) for z in range(2))==wt for a,m in t)
            vectors.append(t);basis.append(t)
        blocks[wt]=vectors
    assert len(basis)==14
    highest=blocks[2,0][0];assert all(not BO.root_action(highest,a) for a in v.raising)
    for tensor in basis:
        residual={}
        for (a,m),c in tensor.items():
            for key,z in BO.constraint_entries(a,m,v.inverse,v.curvature):BO.add(residual,key,c*z)
        assert not residual,('H2 certificate outside K',next(iter(residual.items())))
    tau=[4,6,5,7,0,2,1,3,10,11,8,9,12,13];actions=list(v.raising)
    for action in v.raising:
        opposite=[[] for i in range(14)]
        for a,row in enumerate(action):opposite[tau[a]]=[(tau[b],c) for b,c in row]
        actions.append(opposite)
    for wt,vectors in blocks.items():
        for t in vectors:
            for action,alpha in zip(actions,[(1,-1),(0,1),(-1,1),(0,-1)]):
                out=BO.root_action(t,action);target=tuple(wt[z]+alpha[z] for z in range(2))
                for prior in blocks.get(target,[]):
                    c=out.get(min(prior),0)
                    for key,z in prior.items():BO.add(out,key,-c*z)
                assert not out,(wt,alpha)
    dual=[sharp(t) for t in basis]
    gram=s.Matrix([[sum(c*t.get(k,0) for k,c in d.items()) for t in basis] for d in dual]);assert gram.det()!=0
    bykey={}
    for i,d in enumerate(dual):
        for k,c in d.items():bykey.setdefault(k,[]).append((i,c))
    return bykey

def verify(poly):
    dual=build();count=0
    for mon,tensor in poly.items():
        pair=[Q(0)]*14
        for key,c in tensor.items():
            for i,z in dual.get(key,[]):pair[i]+=c*z
        assert not any(pair),('H2 component',mon,pair)
        count+=1
    return count
