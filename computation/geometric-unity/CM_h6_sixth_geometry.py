"""Direct characteristic-zero quintic closure right hand side."""
from collections import defaultdict
from itertools import combinations
from functools import lru_cache
from fractions import Fraction as Q
import json
import CF_h6_fourth_jet as v
BO=v.BO;HERE=v.HERE;union=v.union
import gzip
AUX=json.loads(gzip.decompress((HERE/'CM_h6_sixth_aux.json.gz').read_bytes()))

def subtract(a,b):
    z=list(a)
    for i in b:z.remove(i)
    return tuple(z)
@lru_cache(None)
def submons(mon,n):return tuple(sorted(set(combinations(mon,n))))

T=[defaultdict(dict) for i in range(14)];U=[defaultdict(dict) for i in range(14)]
for r in range(14):
    for i in range(14):
        for d,row in enumerate(v.cov(r,i)):
            for s,c in row:BO.add(T[i][union((r,),(s,))],d,-c)
for i in range(14):
    for mon,ds in T[i].items():
        for e,c in ds.items():
            for r in range(14):
                for d,row in enumerate(v.cov(r,e)):
                    for s,z in row:BO.add(U[i][union(mon,(r,s))],d,-c*z)
def curve(poly):
    out=[defaultdict(dict) for i in range(14)]
    for i in range(14):
        for mon,ds in poly[i].items():
            for e,c in ds.items():
                for r in range(14):
                    if r!=e:BO.add(out[i][union(mon,(r,))],(r,e),c)
    return out
J=curve(T);L=curve(U)
@lru_cache(None)
def curvature_preimages(key,r,i):
    transposed=[[] for _ in range(14)]
    for a,row in enumerate(v.cov(r,i)):
        for b,c in row:transposed[b].append((a,c))
    return tuple(BO.root_action({key:Q(1)},transposed).items())

class RHS:
    def __init__(self):
        _,self.b2=v.load(2);_,self.b4=v.load(4)
        source=v.PACKED['source']
        self.b={(a,m):Q(c) for a,m,c in source}
    def row(self,label):
        mon,a,mask=label;E=H=F=JJ=LL=Q(0)
        for i in BO.axes(mask):
            key=(a,mask^(1<<i));sgn=BO.wedge_sign(i,mask)
            for q in submons(mon,2):
                rest=subtract(mon,q)
                for d,c in T[i].get(q,{}).items():
                    source=union(rest,(d,));E+=sgn*c*source.count(d)*self.b4.get(source,{}).get(key,0)
            for r in set(mon):
                source=subtract(mon,(r,));tensor=self.b4.get(source,{})
                for pre,c in curvature_preimages(key,r,i):H+=sgn*c*tensor.get(pre,0)
            for q in submons(mon,4):
                rest=subtract(mon,q)
                for d,c in U[i].get(q,{}).items():
                    source=union(rest,(d,));F+=sgn*c*source.count(d)*self.b2.get(source,{}).get(key,0)
            for q in submons(mon,3):
                source=subtract(mon,q);tensor=self.b2.get(source,{})
                for (r,e),c in J[i].get(q,{}).items():
                    for pre,z in curvature_preimages(key,r,e):JJ+=sgn*c*z*tensor.get(pre,0)
            for (r,e),c in L[i].get(mon,{}).items():
                for pre,z in curvature_preimages(key,r,e):LL+=sgn*c*z*self.b.get(pre,0)
        return E/6-H/2-7*F/360+JJ/24-LL/240

def selected_rhs():
    rhs=RHS()
    return [rhs.row((tuple(mon),a,m)) for mon,a,m in AUX['closure_labels']]
