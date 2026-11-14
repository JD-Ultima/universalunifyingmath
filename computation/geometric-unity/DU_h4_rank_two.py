"""Exact O2 x O12 sign-block spectrum for derivative rank at most two.

Only differentiated hooks and coclosure are imposed; curvature
restrictions can only increase the minimum.
"""
from pathlib import Path
from itertools import product
import sympy as s,json
from math import comb
records=[];t=s.Symbol('t')
def sign(i,m):return (-1)**((m&((1<<i)-1)).bit_count())
def matrix(columns):
    keys=sorted(set().union(*map(set,columns)))
    return s.Matrix([[col.get(key,0) for col in columns] for key in keys]) if keys else s.zeros(0,len(columns))
for inside,total in product(range(3),[2,4,6]):
    outside=total-inside;sigma=sum(1<<i for i in range(inside))+sum(1<<i for i in range(2,2+outside))
    columns=[(i,a,sigma^(1<<i)^(1<<a)) for i in range(2) for a in range(14) if (sigma^(1<<i)^(1<<a)).bit_count()==4]
    if not columns:continue
    loc={key:j for j,key in enumerate(columns)};constraints=[];deltas=[];X=s.zeros(len(columns))
    for j,(i,a,m) in enumerate(columns):
        out={};delta={}
        if not m>>a&1:out[0,i,m|(1<<a)]=sign(a,m)
        if m>>a&1:out[1,i,m^(1<<a)]=sign(a,m)
        if m>>i&1:out[2,a,m^(1<<i)]=sign(i,m)
        if not m>>i&1:delta[a,m|(1<<i)]=sign(i,m)
        if (a,i,m) in loc:X[loc[a,i,m],j]=1
        constraints.append(out);deltas.append(delta)
    M=matrix(constraints);h=M.nullspace()
    if not h:continue
    R=s.Matrix.hstack(*h);G=R.T*R;D=matrix(deltas);Q=R.T*(2*X-s.eye(len(columns))+D.T*D)*R
    polynomial=s.factor((Q-t*G).det()/G.det());roots=s.polys.polytools.all_roots(s.Poly(polynomial,t));minimum=min(roots);maximum=max(roots)
    rec=dict(inside=inside,outside=outside,character_multiplicity=comb(2,inside)*comb(12,outside),raw_columns=len(columns),kernel_dimension=len(h),minimum=str(minimum),maximum=str(maximum),normalized_characteristic=str(polynomial));records.append(rec);print(rec,flush=True)
result=dict(records=records,sharp_ambient_minimum=str(min(s.sympify(r['minimum']) for r in records)),sharp_ambient_maximum=str(max(s.sympify(r['maximum']) for r in records)),full_dimension=sum(r['character_multiplicity']*r['kernel_dimension'] for r in records),scope='All real rank-at-most-two derivative tensors satisfying differentiated H4 hooks and closure; no curvature or nonlinear tangency imposed.')
assert result['sharp_ambient_minimum']=='-7/5' and result['sharp_ambient_maximum']=='7/5';print('Exact rank-two endpoints',result['sharp_ambient_minimum'],result['sharp_ambient_maximum'])


"""Raw physical exact H4 rank-two witnesses for both stress endpoints."""
from pathlib import Path
from fractions import Fraction as F
from math import prod
import sys,json,time
HERE=Path(__file__).resolve().parent
import argparse
parser=argparse.ArgumentParser();parser.add_argument('--source-dir',type=Path);args=parser.parse_args()
SOURCE=args.source_dir or (HERE if (HERE/'BC_affine_hooks.py').exists() else HERE.parents[2]/'universalunifyingmath/computation/geometric-unity')
sys.path.insert(0,str(SOURCE))
import BC_affine_hooks as BC
import BG_h6_curvature as BG
tick=time.time();g=BC.METRIC;R=BG.curvature_model();full=(1<<14)-1
def wedge(i,m):
    if m>>i&1:return None
    return m|(1<<i),(-1)**((m&((1<<i)-1)).bit_count())
def star(m):
    n=full^m;sign=(-1)**sum((n&((1<<i)-1)).bit_count() for i in BC.bit_indices(m))
    return n,F(sign,prod(g[i] for i in BC.bit_indices(m)))
omega=1871;assert omega.bit_count()==8 and omega>>12==0
B={}
for i in [12,13]:
    for a in range(14):
        if a==i:continue
        first=wedge(a,omega)
        if first is None:continue
        m,sg=wedge(i,first[0])
        B[i,a,m]=g[i]*g[a]*(F(1) if a in [12,13] else F(-1,4))*first[1]*sg
J={}
for (i,a,m),v in B.items():
    n,c=star(m);J[i,a,n]=v*c
def norms(T):
    norm=F(0);exchange=F(0);D={};delta={}
    for (i,a,m),v in T.items():
        den=g[i]*g[a]*prod(g[j] for j in BC.bit_indices(m))
        norm+=v*v/den;exchange+=v*T.get((a,i,m),0)/den
        hit=wedge(i,m)
        if hit:BC.plus(D,(a,hit[0]),v*hit[1])
        if m>>i&1:BC.plus(delta,(a,m^(1<<i)),v*(-1)**((m&((1<<i)-1)).bit_count())/g[i])
    assert not delta
    dn=sum((v*v/(g[a]*prod(g[j] for j in BC.bit_indices(m))) for (a,m),v in D.items()),F(0))
    return norm,exchange,dn,2*exchange-norm+dn
records=[]
for row_sign in [1,-1]:
    T={key:v*(row_sign if key[0]==13 else 1) for key,v in J.items()}
    for i in [12,13]:
        C={(a,m):v for (j,a,m),v in T.items() if j==i}
        assert BC.hook_traces(C)==({}, {})
        dual={}
        for (a,m),v in C.items():
            n,c=star(m);BC.plus(dual,(a,n),v*c)
        assert not BG.full_sparse_curvature(R,dual)
    n,x,d,q=norms(T);assert q/n==F(-7*row_sign,5)
    assert T and {i for i,a,m in T}=={12,13}
    # The two derivative slices have disjoint coefficient/form support.
    assert not ({(a,m) for i,a,m in T if i==12}&{(a,m) for i,a,m in T if i==13})
    records.append(dict(row_sign=row_sign,norm=str(n),exchange=str(x),wedge_norm=str(d),stress=str(q),ratio=str(q/n),entries=[[i,a,m,str(v)] for (i,a,m),v in sorted(T.items())]))
result=dict(records=records,rank=2,omega_dual_mask=omega,seconds=time.time()-tick,scope='Exact first-symbol membership and stress equality; differentiated nonlinear tangency and higher jet realization are not asserted.')
print({k:v for k,v in result.items() if k!='records'});print([(r['ratio'],len(r['entries'])) for r in records])
