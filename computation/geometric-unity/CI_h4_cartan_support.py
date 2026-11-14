"""Exact grade-five Cartan obstruction for 488 H4 dimensions, then BY for512.

Only BC/BG/X geometry libraries and SymPy are used. No numerical rank or
precomputed quadratic coefficient is trusted.
"""
from collections import defaultdict
from fractions import Fraction as F
from itertools import combinations,combinations_with_replacement
from math import prod
from pathlib import Path
import argparse,gzip,json,time
import sympy as s
import BC_affine_hooks as BC
import BG_h6_curvature as BG
import BC_quadratic_xcheck as X
HERE=Path(__file__).resolve().parent

def combine(terms):
 out=defaultdict(F)
 for c,t in terms:
  for k,v in t.items():out[k]+=c*v
 return {k:v for k,v in out.items() if v}

def add(t,u):return tuple(combine(((1,t[i]),(1,u[i]))) for i in (0,1))

def root_actions():
 C=s.Matrix([[1,1,0,0,0],[s.I,-s.I,0,0,0],[0,0,1,1,0],[0,0,s.I,-s.I,0],[0,0,0,0,1]])
 roots=[]
 for entries in (((2,0,1),(1,3,-1)),((3,4,1),(4,2,-2))):
  r=s.zeros(5)
  for i,j,c in entries:r[i,j]=c
  r=C*r*C.inv();assert r.T==-r
  terms=[]
  for k,(i,j) in enumerate(combinations(range(4),2)):
   if r[i,j]:terms.append((k,F(s.re(r[i,j])),F(s.im(r[i,j]))))
  for i in range(4):
   if r[i,4]:terms.append((6+i,F(2*s.re(r[i,4])),F(2*s.im(r[i,4]))))
  roots.append(terms)
 return roots

def highest(t,w,roots):
 re,im=t
 for k,l in ((0,w[0]),(5,w[1])):
  assert BC.action(k,re)==combine(((l,im),))
  assert BC.action(k,im)==combine(((-l,re),))
 for root in roots:
  real=[];imag=[]
  for k,a,b in root:
   ar,ai=BC.action(k,re),BC.action(k,im)
   real.extend(((a,ar),(-b,ai)));imag.extend(((b,ar),(a,ai)))
  assert not combine(real) and not combine(imag)

def hermitian_inner(t,u):
 def inner(a,b):return sum(v*b.get((i,m),0)/(BC.METRIC[i]*prod(BC.METRIC[j] for j in range(14) if m>>j&1)) for (i,m),v in a.items())
 return s.Rational(inner(t[0],u[0])+inner(t[1],u[1]))+s.I*s.Rational(inner(t[0],u[1])-inner(t[1],u[0]))

def field(t):
 out=[{} for _ in range(14)]
 for (a,m),v in t.items():out[a][m]=v
 return out

def coefficient(t,a,m,engine):
 re,im=t;both=combine(((1,re),(1,im)))
 values=[engine(field(v),a,m) for v in (re,im,both)]
 return s.Rational(values[0]-values[1])+s.I*s.Rational(values[2]-values[0]-values[1])

def exponents(n,d):
 out=[]
 for indices in combinations_with_replacement(range(n),d):
  powers=[0]*n
  for i in indices:powers[i]+=1
  out.append(tuple(powers))
 return out

def modular_rank(matrix,p):
 pivots={};selected=[]
 for rowindex,row in enumerate(matrix):
  row=[int(x)%p for x in row]
  for j,base in sorted(pivots.items()):
   if row[j]:
    c=row[j];row=[(v-c*w)%p for v,w in zip(row,base)]
  if not any(row):continue
  j=next(j for j,v in enumerate(row) if v);c=pow(row[j],-1,p)
  pivots[j]=[(c*v)%p for v in row];selected.append(rowindex)
 return len(pivots),selected

def macaulay(rows,n,d):
 qmon=exponents(n,2);mon=exponents(n,d);pos={e:i for i,e in enumerate(mon)};matrix=[]
 for row in rows:
  for multiplier in exponents(n,d-2):
   out=[0]*len(mon)
   for c,e in zip(row,qmon):out[pos[tuple(x+y for x,y in zip(e,multiplier))]]+=c
   matrix.append(out)
 p=65521;rank,selected=modular_rank(matrix,p);assert rank==len(mon)
 minor=s.Matrix([matrix[i] for i in selected]);det=minor.det(method='domain-ge');assert det and int(det)%p
 return {'degree':d,'rows':len(matrix),'columns':len(mon),'rank_mod_65521':rank,'selected_rows':selected,'exact_minor_determinant':str(det)}

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--full-independent',action='store_true');ap.add_argument('--receipt');args=ap.parse_args()
 start=time.perf_counter();data=json.loads(gzip.decompress((HERE/'CI_h4_cartan_witness.json.gz').read_bytes()))
 roots=root_actions();curvature=BG.curvature_model();independent=X.BladeAlgebra(BC.METRIC)
 character=json.loads((HERE/'BX_data/h4_character.json').read_text(encoding='utf-8'))
 known={tuple(row['weight']):(row['copies'],row['dimension']) for row in character['components']}
 def primary(f,a,m):return BC.quadratic_coefficient(f,BC.METRIC,a,m)
 records=[]
 for group in data['groups']:
  w=group['type'];tensors=[]
  for entries in group['tensors']:
   t=tuple({(a,m):F(v[i]) for a,m,*v in entries if F(v[i])} for i in (0,1))
   assert all(m.bit_count()==4 for a,m,*v in entries)
   highest(t,w,roots)
   for part in t:
    assert BC.hook_traces(part)==({}, {})
    assert not BG.full_sparse_curvature(curvature,BC.hodge_scaled(part))
   tensors.append(t)
  n=len(tensors);assert n==group['multiplicity']
  assert (n,group['per_copy'])==known[tuple(w)]
  gram=s.Matrix(n,n,lambda i,j:hermitian_inner(tensors[i],tensors[j]));assert gram.det()!=0
  pairs=list(combinations_with_replacement(range(n),2));rows=[];guards=0
  for r,((a,m),stored) in enumerate(zip(group['receivers'],group['rows'])):
   assert m.bit_count()==5
   diag=[coefficient(t,a,m,primary) for t in tensors];row=[]
   for index,(i,j) in enumerate(pairs):
    c=diag[i] if i==j else coefficient(add(tensors[i],tensors[j]),a,m,primary)-diag[i]-diag[j]
    want=s.sympify(stored['phase'])*s.Rational(stored['scale'])*stored['integers'][index];assert c==want
    if args.full_independent or (r==0 and index in (0,1,len(pairs)-1)):
     e=independent.quadratic_coefficient
     q=coefficient(tensors[i],a,m,e) if i==j else coefficient(add(tensors[i],tensors[j]),a,m,e)-coefficient(tensors[i],a,m,e)-coefficient(tensors[j],a,m,e)
     assert q==c;guards+=1
    row.append(stored['integers'][index])
   rows.append(row)
  rank=macaulay(rows,n,group['macaulay_degree'])
  rec={'type':w,'multiplicity':n,'per_copy':group['per_copy'],'highest_weight_hooks_and_full_dual_curvature':True,'highest_gram_determinant':str(gram.det()),'all_selected_coefficients_reconstructed':len(rows)*len(pairs),'independent_guards':guards,'macaulay':rank}
  print(json.dumps(rec),flush=True);records.append(rec)
 assert sum(g['multiplicity']*g['per_copy'] for g in data['groups'])==488
 chain=[tuple(g['type']) for g in data['groups']]+[(2,0),(1,0)]
 for a,b in zip(chain,chain[1:]):assert a[0]>=b[0] and sum(a)>=sum(b) and a!=b
 result={'exact':True,'records':records,'Cartan_support_dimension':488,'with_BY_low_support_dimension':512,'complex_exclusion':True,'scope':'only the indicated complete isotypic sum; not the full H4 kernel','seconds':time.perf_counter()-start}
 if args.receipt:Path(args.receipt).write_text(json.dumps(result,indent=2))
 print(json.dumps(result,indent=2))
if __name__=='__main__':main()
