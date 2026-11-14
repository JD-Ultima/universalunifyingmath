"""Exact exterior, contraction, and operator-square identities for H4/H6.

The reduced Clifford calculation is checked against BC's independent nested
commutator formula on complete sparse rational test tensors in both grades.
No numerical rank or positivity test is used. Run with --receipt FILE to save
the exact validation record.
"""
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations
from pathlib import Path
import argparse,json,random,time
import BC_affine_hooks as BC
N=14;FULL=(1<<N)-1
def add(*terms):
 out={}
 for term in terms:
  for m,v in term.items():out[m]=out.get(m,0)+v
 return {m:v for m,v in out.items() if v}
def scale(A,c):return {m:v*c for m,v in A.items() if v*c}
@lru_cache(None)
def factor(a,b):
 sign=1
 # Append the right word in increasing generator order.
 for j in range(N):
  if b>>j&1:
   if (a>>(j+1)).bit_count()%2:sign=-sign
   a^=1<<j
 return sign
def mul(A,B):
 out={}
 for a,x in A.items():
  for b,y in B.items():
   m=a^b;out[m]=out.get(m,0)+x*y*factor(a,b)
 return {m:v for m,v in out.items() if v}
def dagger(A):return {m:v*(-1)**(m.bit_count()*(m.bit_count()-1)//2) for m,v in A.items()}
def grade(A,r):return {m:v for m,v in A.items() if m.bit_count()==r}
def interior(a,A):return {m^(1<<a):v*(-1)**((m&((1<<a)-1)).bit_count()) for m,v in A.items() if m>>a&1}
def wedge(a,A):return {m|(1<<a):v*(-1)**((m&((1<<a)-1)).bit_count()) for m,v in A.items() if not m>>a&1}
def exterior_product(A,B):return {m:v for m,v in mul(A,B).items() if m.bit_count()==next(iter(A)).bit_count()+next(iter(B)).bit_count()} if A and B else {}
GAMMA=[{1<<a:1} for a in range(N)]
def hook_sample(degree,seed,terms=12):
 rng=random.Random(seed);B=[{} for _ in range(N)]
 for j in range(terms):
  indices=rng.sample(range(N),degree+1);a,b=indices[:2];tail=sum(1<<i for i in indices[2:]);v=F(rng.choice([-3,-2,-1,1,2,3]),rng.choice([1,2,3]))
  x,y,s=(a,b,-1) if j%2==0 else (b,a,1)
  B[a]=add(B[a],scale(wedge(x,{tail:1}),v));B[b]=add(B[b],scale(wedge(y,{tail:1}),s*v))
 if degree==6:
  # A deliberately nonzero top exterior square, independent of random luck.
  for b,tail in [(0,sum(1<<i for i in range(1,6))),(6,sum(1<<i for i in range(7,12)))]:
   B[12]=add(B[12],wedge(b,{tail:1}));B[b]=add(B[b],wedge(12,{tail:1}))
 assert not add(*(wedge(a,B[a]) for a in range(N)))
 assert not add(*(interior(a,B[a]) for a in range(N)))
 return B
def reduced(C):
 products=[[mul(a,b) for b in C] for a in C];S=add(*(products[a][a] for a in range(N)));Q=[]
 for a in range(N):
  U=add(*(mul(GAMMA[b],products[a][b]) for b in range(N)))
  V=add(*(mul(products[b][a],GAMMA[b]) for b in range(N)))
  W=add(*(mul(mul(C[b],GAMMA[a]),C[b]) for b in range(N)))
  Q.append(scale(add(mul(GAMMA[a],S),mul(S,GAMMA[a]),scale(U,-1),scale(V,-1),scale(W,-1)),-F(2,3)))
 V=add(*(mul(mul(GAMMA[a],products[b][a]),GAMMA[b]) for a in range(N) for b in range(N)))
 return S,V,Q
def verify(degree,seed):
 B=hook_sample(degree,seed);C=B if degree==4 else [mul(b,{FULL:1}) for b in B];k=4 if degree==4 else 8
 assert all(dagger(c)==c for c in C)
 assert not add(*(mul(GAMMA[a],C[a]) for a in range(N)))
 assert not add(*(mul(C[a],GAMMA[a]) for a in range(N)))
 S,V,Q=reduced(C);reference=BC.Clifford([1]*N).quadratic(C);assert Q==reference
 assert all(dagger(q)==q and all(m.bit_count()%4==1 for m in q) for q in Q)
 FS=add(*(mul(mul(g,S),g) for g in GAMMA));left=add(*(mul(GAMMA[a],Q[a]) for a in range(N)));right=add(*(mul(Q[a],GAMMA[a]) for a in range(N)))
 assert left==right==dagger(left)
 assert left==scale(add(scale(S,2*k-2),FS,scale(V,-1)),-F(2,3))
 A=[[add(interior(a,C[b]),scale(interior(b,C[a]),-1)) for b in range(N)] for a in range(N)]
 assert all(dagger(A[a][b])==scale(A[a][b],-1) for a in range(N) for b in range(N))
 P=scale(add(*(mul(A[a][b],A[a][b]) for a in range(N) for b in range(N))),-1)
 assert left==scale(add(scale(S,N-2*k),scale(P,2)),-F(2,3))
 assert P.get(0,0)==2*(k-1)*S[0]
 for a in range(N):assert add(*(mul(GAMMA[b],A[a][b]) for b in range(N)))==scale(C[a],-(k-1))
 square={}
 for a in range(N):
  for b in range(N):
   remainder=add(A[a][b],scale(add(mul(GAMMA[b],C[a]),scale(mul(GAMMA[a],C[b]),-1)),F(k-1,N-2)))
   square=add(square,mul(dagger(remainder),remainder))
 assert square==add(P,scale(S,-F(2*(k-1)**2,N-2)))
 components={};contractions={}
 for r in (1,5,9,13):
  qr=[grade(q,r) for q in Q];alt=add(*(wedge(a,qr[a]) for a in range(N)));ctr=add(*(interior(a,qr[a]) for a in range(N)));assert not alt
  expected=scale(add(scale(grade(S,r-1),N+2*k-2-2*(r-1)),scale(grade(V,r-1),-1)),-F(2,3));assert ctr==expected
  if r==1:assert ctr=={0:-F(2,3)*(N+2*k-4)*S[0]}
  trace=[scale(wedge(a,ctr),F(1,N-r+1)) for a in range(N)];hook=[add(q,scale(t,-1)) for q,t in zip(qr,trace)]
  assert not add(*(interior(a,hook[a]) for a in range(N)))
  components[r]={'trace_nonzero':bool(ctr),'hook_nonzero':any(hook),'coefficients':sum(map(len,qr))};contractions[r]=len(ctr)
 top=2*degree+1;omega=add(*(exterior_product(b,b) for b in B));coefficient=-F(2,3) if degree==4 else F(2)
 assert all(grade(Q[a],top)==scale(wedge(a,omega),coefficient) for a in range(N))
 assert grade(V,top-1)==scale(grade(S,top-1),-2)
 assert add(*(interior(a,grade(Q[a],top)) for a in range(N)))==scale(grade(S,top-1),-4)
 independence={}
 for r in ([4] if degree==4 else [4,8]):
  pairs=[(m,S.get(m,0),V.get(m,0)) for m in sorted(S.keys()|V.keys()) if m.bit_count()==r]
  witness=next(((a,b,x*w-y*z) for (a,x,y),(b,z,w) in combinations(pairs,2) if x*w-y*z),None)
  assert witness;independence[r]=[witness[0],witness[1],str(witness[2])]
 return {'degree':degree,'clifford_degree':k,'seed':seed,'input_coefficients':sum(map(len,C)),'complete_independent_reference_coefficients':sum(map(len,Q)),'all_exterior_identities':True,'left_and_right_contraction_identity':True,'positive_square_identity':True,'top_trace_identity':True,'output_components':components,'contraction_coefficients':contractions,'middle_S_V_independence_minors':independence}
def main():
 parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--receipt');args=parser.parse_args();start=time.perf_counter();records=[verify(d,s) for d,s in [(4,1404),(4,2404),(6,1406),(6,2406)]]
 for d in (4,6):
  part=[r for r in records if r['degree']==d]
  assert any(r['output_components'][1]['hook_nonzero'] for r in part)
  assert any(r['output_components'][5]['trace_nonzero'] and r['output_components'][5]['hook_nonzero'] for r in part)
  assert any(r['output_components'][2*d+1]['trace_nonzero'] for r in part)
  if d==6:assert any(r['output_components'][9]['trace_nonzero'] and r['output_components'][9]['hook_nonzero'] for r in part)
 result={'exact':True,'records':records,'h4_operator_cap':'2 N','h6_operator_cap':'78 N /43','h4_rank_lower_bound':64,'h6_rank_lower_bound':72,'seconds':time.perf_counter()-start}
 if args.receipt:Path(args.receipt).write_text(json.dumps(result,indent=2))
 print(json.dumps(result,indent=2))
if __name__=='__main__':main()
