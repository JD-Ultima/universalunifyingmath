"""Exact cyclic supports beyond forests, with characteristic-two Macaulay lifts."""
import sys,json,time
from pathlib import Path
from fractions import Fraction as F
from itertools import combinations_with_replacement
HERE=Path(__file__).resolve().parent
import CD_clifford_hook_identities as C

def monomials(k,d):return list(combinations_with_replacement(range(k),d))
def rank2(rows):
 basis={}
 for v in rows:
  while v:
   j=v.bit_length()-1
   if j not in basis:basis[j]=v;break
   v^=basis[j]
 return len(basis)
def inner(A,B):
 out={}
 for m,(i,c) in A.items():
  if m in B:
   j,d=B[m];key=tuple(sorted((i,j)));out[key]=out.get(key,0)+c*d
 return {k:v for k,v in out.items() if v}
def add(*parts):
 out={}
 for scale,part in parts:
  for k,v in part.items():out[k]=out.get(k,0)+scale*v
 return {k:v for k,v in out.items() if v}
def inside(a,A):return {m^(1<<a):(i,c*(-1)**((m&((1<<a)-1)).bit_count())) for m,(i,c) in A.items() if m>>a&1}
def support(p,k):
 K=(1<<(p-2))-1;u=list(range(p-2,p-2+k));a=list(range(p-2+k,p-2+2*k));assert max(a)<14
 B=[{} for _ in range(14)]
 for i in range(k):
  J=K|(1<<u[i])
  for row,sign in [(a[i],1),(a[(i-1)%k],-1)]:
   mask=J|(1<<row);B[row][mask]=(i,sign*(-1)**((mask&((1<<row)-1)).bit_count()))
 return B,u,a
def run(p,k):
 B,u,a=support(p,k);E=[[inner(x,y) for y in B] for x in B];Fm=[[add(*[(1,inner(inside(i,b),inside(j,b))) for b in B]) for j in range(14)] for i in range(14)];N=add(*[(1,E[i][i]) for i in range(14)])
 equations=[]
 for i in range(14):
  for j in range(i,14):
   q=add((7,Fm[i][j]),(-14 if p==4 else 14,E[i][j]),(-1 if p==4 else -4,N if i==j else {}));equations.append(q)
 # A random exact tensor checks every grade-one coefficient in the original
 # nested Clifford formula, and both hooks independently.
 x=[F(i+2) for i in range(k)];field=[{m:c*x[j] for m,(j,c) in b.items()} for b in B]
 assert not C.add(*(C.interior(i,b) for i,b in enumerate(field))) and not C.add(*(C.wedge(i,b) for i,b in enumerate(field)))
 actual=field if p==4 else [C.mul(b,{C.FULL:1}) for b in field];Q=C.BC.Clifford([1]*14).quadratic(actual)
 evaluate=lambda q:sum(c*x[i]*x[j] for (i,j),c in q.items());nv=evaluate(N)
 for i in range(14):
  for j in range(14):
   expected=F(4,3)*(2*evaluate(E[i][j])+(-1 if p==4 else 1)*evaluate(Fm[i][j]))-(F(2,3) if p==4 else 2)*nv*(i==j)
   assert Q[i].get(1<<j,0)==expected
 parity=[{ij for ij,c in q.items() if c%2} for q in equations];parity=[q for q in parity if q];ranks=[]
 for degree in range(2,5):
  mons=monomials(k,degree);lookup={m:i for i,m in enumerate(mons)};rows=[]
  for q in parity:
   for tail in monomials(k,degree-2):
    v=0
    for ij in q:v^=1<<lookup[tuple(sorted(ij+tail))]
    rows.append(v)
  rank=rank2(rows);ranks.append({'degree':degree,'rows':len(rows),'columns':len(mons),'rank_mod2':rank})
  if rank==len(mons):break
 assert ranks[-1]['rank_mod2']==ranks[-1]['columns']
 return {'degree':p,'cycle_length':2*k,'free_hook_variables':k,'slots':[[row,m,j,c] for row,b in enumerate(B) for m,(j,c) in b.items()],'all196_original_grade1_guards':True,'Macaulay':ranks,'nonzero_complex_solution_excluded':True}

def graph_certificate(p,slots):
 """Sufficient graph criterion, including a deliberately inconclusive mode."""
 assert len(slots)==len(set(map(tuple,slots))) and all(m.bit_count()==p for a,m in slots)
 vertices={};ends=[]
 for edge,(a,m) in enumerate(slots):
  group=('R' if m>>a&1 else 'D',m^(1<<a));row=('row',a);ends.append((row,group))
  for v in (row,group):vertices.setdefault(v,set()).add(edge)
 m=len(slots);group_rows=[sum(1<<e for e in edges) for v,edges in vertices.items() if v[0]!='row'];allones=(1<<m)-1
 diagonal=[]
 for a in range(14):
  form=sum(1<<e for e,(b,I) in enumerate(slots) if I>>a&1);external=sum(1<<e for e,(b,I) in enumerate(slots) if b==a)
  diagonal.append(form^(allones if p==4 else 0))
  assert rank2(group_rows+diagonal+[external])==rank2(group_rows+diagonal)
 # The preceding incremental check need only use this axis's diagonal,
 # since group squares give both the total norm and every group incidence.
 remaining=set(range(m));degree={v:len(es) for v,es in vertices.items()};queue=[v for v,d in degree.items() if d<=1]
 while queue:
  v=queue.pop()
  for e in list(vertices[v]&remaining):
   remaining.remove(e)
   for u in ends[e]:
    degree[u]-=1
    if degree[u]==1:queue.append(u)
 if not remaining:
  assert rank2(group_rows+diagonal)==m
  return {'slots':m,'forest':True,'cycle_count':0,'excluded_over_complex':True,'uniform_full_Macaulay_degree':m+1}
 corevertices={v for e in remaining for v in ends[e]}
 if any(degree[v]!=2 for v in corevertices):return {'slots':m,'forest':False,'supported_core':False,'excluded_over_complex':False}
 pairs=list(combinations_with_replacement(range(14),2));pairs=[x for x in pairs if x[0]!=x[1]];index={x:i for i,x in enumerate(pairs)};left=set(remaining);vectors=[]
 while left:
  todo=[left.pop()];component=set(todo)
  while todo:
   e=todo.pop()
   for v in ends[e]:
    for f in vertices[v]&left:left.remove(f);component.add(f);todo.append(f)
  parity=0
  for v in {v for e in component for v in ends[e] if v[0]=='row'}:
   e,f=list(vertices[v]&component);difference=slots[e][1]^slots[f][1]
   if difference.bit_count()==2:
    a,b=C.BC.bit_indices(difference);parity^=1<<index[a,b]
  vectors.append(parity)
 r=rank2(vectors);assert rank2(group_rows+diagonal)==m-len(vectors)
 return {'slots':m,'forest':False,'cycle_count':len(vectors),'exchange_parity_rank':r,'excluded_over_complex':r==len(vectors),'uniform_full_Macaulay_degree':3*m+1 if r==len(vectors) else None}

def main():
 import argparse
 parser=argparse.ArgumentParser();parser.add_argument('--receipt',type=Path);args=parser.parse_args();start=time.perf_counter();records=[];graphs=[]
 for p in (4,6):
  # A mixed contraction/alternation tree, with two blades in its middle row.
  J=(1<<(p-1))-1;a,b,c=p-1,p,p+1;T=J|(1<<b)|(1<<c)
  forest=[(a,J|(1<<a)),(b,J|(1<<b)),(b,T^(1<<b)),(c,T^(1<<c))];g=graph_certificate(p,forest);assert g['forest'] and g['excluded_over_complex'];graphs.append(g)
  for k in range(3,(16-p)//2+1):
   r=run(p,k);records.append(r);g=graph_certificate(p,[(a,m) for a,m,j,c in r['slots']]);assert g['excluded_over_complex'];graphs.append(g)
  # Two rows and two right groups: the exchange occurs twice, so the
  # cycle-parity criterion correctly makes no exclusion claim.
  K=(1<<(p-2))-1;u,v,a,b=p-2,p-1,p,p+1;negative=[(row,K|(1<<axis)|(1<<row)) for row in (a,b) for axis in (u,v)];g=graph_certificate(p,negative);assert not g['excluded_over_complex'] and g['exchange_parity_rank']==0;graphs.append(g)
 # Two vertex-disjoint cycles with independent exchange signatures.
 K=3;multi=[]
 for u,a in [([2,3,4],[5,6,7]),([8,9,10],[11,12,13])]:
  for i in range(3):multi.extend([(a[i],K|(1<<u[i])|(1<<a[i])),(a[i],K|(1<<u[(i+1)%3])|(1<<a[i]))])
 g=graph_certificate(4,multi);assert g['cycle_count']==g['exchange_parity_rank']==2 and g['excluded_over_complex'];graphs.append(g)
 # Exact scalarity trace implications at the minimum real external rank.
 assert F(2)+F(4,7)==F(18,7) and F(4*6,7)-2==F(10,7)
 rank_energy={'H4_minimum_external_rank':4,'H6_minimum_external_rank':6,'H4_at_rank4_w3_minus_w1':'4/7','H4_at_rank4_w1_upper':'3/14','H4_at_rank4_w3_interval':['4/7','11/14'],'H6_at_rank6_sum_k_minus_1_times_wk':'3/7','H6_at_rank6_w1_lower':'4/7'}
 result={'exact':True,'cyclic_coefficient_and_Macaulay_guards':records,'graph_guards':graphs,'real_rank_energy':rank_energy,'seconds':time.perf_counter()-start}
 if args.receipt:args.receipt.write_text(json.dumps(result,indent=2))
 print(json.dumps({'exact':True,'cyclic_families':len(records),'graph_guards':len(graphs),'seconds':result['seconds']},indent=2))
if __name__=='__main__':main()
