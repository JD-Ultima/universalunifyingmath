"""Explicit rank6 real H6 hook with full scalarity and zero top grade."""
from pathlib import Path
from fractions import Fraction as F
from itertools import combinations
import sys,json,time
HERE=Path(__file__).resolve().parent;import CD_clifford_hook_identities as C
start=time.perf_counter();R=list(range(6));T=list(range(6,14));fullT=sum(1<<i for i in T)
# Standard SU(3) Cartan3-form on8 Euclidean axes. Its norm squared is4
# and every contraction has norm squared3/2. Split rational/sqrt3 parts.
f0={(1,2,3):F(1),(1,4,7):F(1,2),(2,4,6):F(1,2),(2,5,7):F(1,2),(3,4,5):F(1,2),(1,5,6):-F(1,2),(3,6,7):-F(1,2)};f3={(4,5,8):F(1,2),(6,7,8):F(1,2)}
def hodge_t(form):
 out={}
 for inds,v in form.items():
  mask=sum(1<<T[i-1] for i in inds);complement=fullT^mask;sign=(-1)**sum((complement&((1<<i)-1)).bit_count() for i in C.BC.bit_indices(mask));out[complement]=sign*v
 return out
eta0,eta3=hodge_t(f0),hodge_t(f3);d=[1,1,1,-1,-1,-1];U=[C.scale(C.wedge(a,eta0),d[a]) if a<6 else {} for a in range(14)];V=[C.scale(C.wedge(a,eta3),d[a]) if a<6 else {} for a in range(14)];W=[{} for _ in range(14)]
# Six distinct4-blades, as three complementary pairs, make an isotropic
# aggregate4-form density. A pair of3-cycles gives balanced hook columns.
halves=[(0,1,2,3),(0,1,4,5),(0,2,4,6)];omega=[]
for half in halves:
 mask=sum(1<<T[i] for i in half);omega.extend([mask,fullT^mask])
matrix=[[0]*6 for _ in range(6)]
for group in ((0,1,2),(3,4,5)):
 for i,a in enumerate(group):matrix[a][group[(i+1)%3]]=1;matrix[a][group[(i-1)%3]]=-1
for a in R:
 for b in R:
  if matrix[a][b]:W[a]=C.add(W[a],C.scale(C.wedge(a,C.wedge(b,{omega[b]:1})),matrix[a][b]))
parts=[U,V,W]
for B in parts:
 assert not C.add(*(C.wedge(a,b) for a,b in enumerate(B))) and not C.add(*(C.interior(a,b) for a,b in enumerate(B)));assert all(all(m>>a&1 for m in b) for a,b in enumerate(B))
def quadratic(B):return C.BC.Clifford([1]*14).quadratic([C.mul(b,{C.FULL:1}) for b in B])
Qs=[quadratic(B) for B in parts];cross={}
for i,j in combinations(range(3),2):
 Q=quadratic([C.add(a,b) for a,b in zip(parts[i],parts[j])]);cross[i,j]=[C.add(q,C.scale(x,-1),C.scale(y,-1)) for q,x,y in zip(Q,Qs[i],Qs[j])]
# B=sqrt2 U+sqrt6 V+sqrt3 W. The four independent radical channels are
# 1, sqrt2, sqrt3, sqrt6; their quadratic coefficients follow below.
channels={'rational':[C.add(C.scale(a,2),C.scale(b,6),C.scale(c,3)) for a,b,c in zip(*Qs)],'sqrt2':[C.scale(q,3) for q in cross[1,2]],'sqrt3':[C.scale(q,2) for q in cross[0,1]],'sqrt6':cross[0,2]}
inner=lambda x,y:sum(v*y.get(m,0) for m,v in x.items());N=sum(scale*sum(inner(b,b) for b in B) for scale,B in zip((2,6,3),parts));assert N==84
E=[[sum(scale*inner(B[a],B[b]) for scale,B in zip((2,6,3),parts)) for b in range(14)] for a in range(14)];assert E==[[14 if a==b and a<6 else 0 for b in range(14)] for a in range(14)]
for label,Q in channels.items():
 for a in range(14):
  assert C.grade(Q[a],1)==({1<<a:-104} if label=='rational' else {});assert not C.grade(Q[a],13)
result={'exact':True,'field':'sqrt(2) U+sqrt(6) V+sqrt(3) W','norm_squared':84,'external_Gram':'14 I6 direct_sum 0_8','external_rank':6,'full_hooks_scalarity_and_grade13':True,'curvature_not_imposed':True,'R_degree1_energy_fraction':'4/7','R_degree2_energy_fraction':'3/7','higher_nonzero_coefficients':{label:{str(g):sum(len(C.grade(q,g)) for q in Q) for g in (5,9)} for label,Q in channels.items()},'parts':[[[a,m,str(v)] for a,b in enumerate(B) for m,v in b.items()] for B in parts],'seconds':time.perf_counter()-start};import argparse
parser=argparse.ArgumentParser();parser.add_argument('--receipt',type=Path);args=parser.parse_args()
if args.receipt:args.receipt.write_text(json.dumps(result,indent=2))
print({k:v for k,v in result.items() if k!='parts'})
