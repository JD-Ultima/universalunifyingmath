"""Exact ordinary S5 matrix-conjugation action and character counts."""
from fractions import Fraction as F
from functools import lru_cache
from itertools import permutations,product,combinations_with_replacement
from collections import Counter
from math import lcm
import numpy as np
import sympy as s
import BC_affine_hooks as BC
import BG_h6_curvature as BG
def permutation_actions():
 phi=[s.Matrix(x)/5 for x in BG.integer_phi()];maps=[]
 for i in range(4):
  P=s.eye(5);P.row_swap(i,i+1)
  U=s.Matrix(14,14,lambda a,b:s.trace(phi[a]*P*phi[b]*P.T)/s.Rational(BC.METRIC[a]));assert U*U==s.eye(14)
  metric=s.diag(*map(s.Rational,BC.METRIC));assert U.T*metric*U==metric
  rows=[[(j,F(v)) for j,v in enumerate(row) if v] for row in U.tolist()]
  @lru_cache(None)
  def wedge(mask,rows=rows):
   out={0:F(1)}
   for a in BC.bit_indices(mask):
    new={}
    for m,c in out.items():
     for b,v in rows[a]:
      if not m>>b&1:BC.plus(new,m|(1<<b),c*v*(-1)**((m>>(b+1)).bit_count()))
    out=new
   return tuple(out.items())
  maps.append((rows,wedge))
 return maps
def apply(t,operation):
 rows,wedge=operation;out={}
 for (a,m),c in t.items():
  for b,v in rows[a]:
   for n,d in wedge(m):BC.plus(out,(b,n),c*v*d)
 return out
@lru_cache(None)
def partition(a,b):
 if a<0:return 0
 return sum(b+r1-r4>=0 for r1 in range(a+1) for r4 in range(a-r1+1))
def character(l,m):
 out={};highest=(2*l+3,2*m+1)
 for a,b in product(range(-l,l+1),repeat=2):
  v=0
  for perm in permutations((0,1)):
   for signs in product((-1,1),repeat=2):
    x=tuple(signs[i]*highest[perm[i]] for i in range(2));x=((x[0]-(2*a+3))//2,(x[1]-(2*b+1))//2)
    v+=(1 if perm==(0,1) else -1)*signs[0]*signs[1]*partition(*x)
  assert v>=0
  if v:out[a,b]=v
 assert sum(out.values())==(2*l+3)*(2*m+1)*(l-m+1)*(l+m+2)//6
 return out
def restriction(l,m):
 weights=character(l,m);x=s.Symbol('x');classes=[(1,0,0,1,1),(2,1,1,10,-1),(2,1,0,15,1),(3,1,0,20,1),(6,3,1,20,-1),(4,2,1,30,-1),(5,1,2,24,1)];values=[]
 for n,a,b,size,sign in classes:
  poly=sum(c*x**((a*i+b*j)%n) for (i,j),c in weights.items());v=s.rem(poly,s.cyclotomic_poly(n,x),x);assert v.is_Integer;values.append(int(v))
 totals=[sum(v*cl[3]*(cl[4] if sign else 1) for v,cl in zip(values,classes)) for sign in (False,True)];assert all(t%120==0 for t in totals)
 return values,[t//120 for t in totals]

def rank(A,p=65521):
 A=np.array(A,dtype=np.int64).copy()%p;r=0
 for j in range(A.shape[1]):
  hits=np.flatnonzero(A[r:,j])
  if not len(hits):continue
  k=r+int(hits[0]);A[[r,k]]=A[[k,r]];A[r,j:]=A[r,j:]*pow(int(A[r,j]),-1,p)%p
  A[r+1:,j:]=(A[r+1:,j:]-A[r+1:,j,None]*A[r,None,j:])%p;r+=1
  if r==A.shape[0]:break
 return r
def check_basis(w):
 tensors=[{(a,m):F(v) for a,m,v in t} for t in w['tensors']];R=BG.curvature_model();ops=permutation_actions()
 for t in tensors:
  assert all(m.bit_count()==4 for a,m in t)
  assert BC.hook_traces(t)==({},{});assert not BG.full_sparse_curvature(R,BC.hodge_scaled(t))
  assert all(apply(t,op)==t for op in ops)
 slots=sorted(set().union(*(t.keys() for t in tensors)));p=65521
 assert rank([[int(t.get(slot,0))%p for slot in slots] for t in tensors])==68
 counts={tuple(weight):copies*restriction(*weight)[1][0] for weight,copies in w['inventory']}
 assert sum(counts.values())==68
 assert all(b['dimension']==counts[tuple(b['weight'])] for b in w['blocks'])
 assert rank([[int(v)%p for v in row] for row in w['preparation']['basis_coefficients']])==68
 return {'hooks_curvature_four_S5_generators':True,'source_basis_rank':68,'chosen_basis_change_rank':68,'full_fixed_dimension_from_CA_inventory':68}

def output_count():
 classes=[([1,1,1,1,1],1),([2,1,1,1],10),([2,2,1],15),([3,1,1],20),([3,2],20),([4,1],30),([5],24)]
 counts=[];symmetric=[];x=s.Symbol('x')
 for cycles,size in classes:
  n=lcm(*cycles);v=[j*n//k for k in cycles for j in range(k)];w=[(a+b)%n for a,b in combinations_with_replacement(v,2)];w.remove(0);char=[]
  E=[Counter() for _ in range(15)];E[0][0]=1
  for value in w:
   for k in range(14,0,-1):
    for e,c in list(E[k-1].items()):E[k][(e+value)%n]+=c
  for k in range(15):
   z=s.rem(sum(c*x**e for e,c in E[k].items()),s.cyclotomic_poly(n,x),x);assert z.is_Integer;char.append(int(z))
  sym=Counter((a+b)%n for a,b in combinations_with_replacement(w,2));z=s.rem(sum(c*x**e for e,c in sym.items()),s.cyclotomic_poly(n,x),x);assert z.is_Integer;symmetric.append((size,int(z)))
  counts.append((size,char));print(cycles,char[1],char[5],char[6],char[8])
 print('VtensorLambda5 invariant',sum(n*c[1]*c[5] for n,c in counts)//120)
 print('Lambda6 invariant',sum(n*c[6] for n,c in counts)//120)
 print('Lambda8 invariant',sum(n*c[8] for n,c in counts)//120)
 result={'symmetric_square_invariants':sum(n*c for n,c in symmetric)//120,'V_tensor_Lambda5_invariants':sum(n*c[1]*c[5] for n,c in counts)//120,'Lambda6_invariants':sum(n*c[6] for n,c in counts)//120,'Lambda8_invariants':sum(n*c[8] for n,c in counts)//120}
 assert result=={'symmetric_square_invariants':5,'V_tensor_Lambda5_invariants':223,'Lambda6_invariants':28,'Lambda8_invariants':28}
 result['total_homogeneous_output_bound']=5-1+223-28+28
 return result
