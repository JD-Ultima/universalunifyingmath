"""Exact L5 multiplicity forms on the entire10388D H4 curvature kernel."""
from collections import defaultdict
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations,product
from math import prod,lcm
from pathlib import Path
import argparse,gzip,json,time
import sympy as s
import BC_affine_hooks as BC
import BG_h6_curvature as BG
import BF_h6_casimir as BF
import BQ_h6_invariant_forms as BQ
import BC_quadratic_xcheck as X
HERE=Path(__file__).resolve().parent;G=tuple(BC.METRIC);FULL=(1<<14)-1;DET=prod(G)

def combine(terms):
 out=defaultdict(F)
 for c,t in terms:
  for key,v in t.items():out[key]+=c*v
 return {key:v for key,v in out.items() if v}

def roots():
 C=s.Matrix([[1,1,0,0,0],[s.I,-s.I,0,0,0],[0,0,1,1,0],[0,0,s.I,-s.I,0],[0,0,0,0,1]])
 answer=[]
 for entries in (((2,0,1),(1,3,-1)),((3,4,1),(4,2,-2))):
  r=s.zeros(5)
  for i,j,c in entries:r[i,j]=c
  r=C*r*C.inv();assert r.T==-r
  terms=[(k,F(s.re(r[i,j])),F(s.im(r[i,j]))) for k,(i,j) in enumerate(combinations(range(4),2)) if r[i,j]]
  terms.extend((6+i,F(2*s.re(r[i,4])),F(2*s.im(r[i,4]))) for i in range(4) if r[i,4]);answer.append(terms)
 return answer

def highest(re,w,raising):
 im=combine(((F(1,w[0]),BC.action(0,re)),))
 for k,l in ((0,w[0]),(5,w[1])):
  assert BC.action(k,re)==combine(((l,im),))
  assert BC.action(k,im)==combine(((-l,re),))
 for root in raising:
  real=[];imag=[]
  for k,a,b in root:
   ar,ai=BC.action(k,re),BC.action(k,im)
   real.extend(((a,ar),(-b,ai)));imag.extend(((b,ar),(a,ai)))
  assert not combine(real) and not combine(imag)
 # Physical conjugation is the same SO5 Weyl element on every copy.
 for tensor,sign in ((re,1),(im,-1)):
  assert all((-1)**((BC.char_of(a,m)&10).bit_count())==sign for a,m in tensor)

PREF=[]
for mask in range(1<<14):
 z=mask
 for k in (1,2,4,8):z^=z<<k
 PREF.append((z<<1)&FULL)
PAR=[1-(-1)**xy+(-1)**mx+(-1)**(mx+my)-(-1)**(xy+my)-(-1)**(mx+my+xy) for mx,my,xy in product(range(2),repeat=3)]
WEIGHT=[int(5*prod(G[j] for j in range(14) if not mask>>j&1)) for mask in range(1<<14)]
TRIPLES={a:[(b,c,(1<<a)|(1<<b)|(1<<c),BC.sign(1<<a,1<<b)*BC.sign((1<<a)|(1<<b),1<<c)) for b,c in combinations(range(14),2) if a not in (b,c)] for a in range(14)}

def value(t,receiver,independent,full,matrix_checker=None):
 den=lcm(*(v.denominator for v in t.values()));fields=[{} for _ in range(14)]
 for (a,m),v in t.items():fields[a][m]=int(v*den)
 coefficients=[];answer=F(0)
 for a,target,weight in receiver:
  total=0
  for b,c,m,ts in TRIPLES[a]:
   for x,vx in fields[b].items():
    y=m^x^target;vy=fields[c].get(y,0)
    if not vy:continue
    parity=PAR[4*((m&x).bit_count()&1)+2*((m&y).bit_count()&1)+((x&y).bit_count()&1)]
    if not parity:continue
    sign=1-2*((((m&PREF[x])^((m^x)&PREF[y])).bit_count())&1)
    total+=ts*parity*sign*WEIGHT[(m&x)|(m&y)|(x&y)]*vx*vy
  q=-F(2,15)*G[a]*F(total,den**2)/DET;answer+=weight*q;coefficients.append((a,target,q))
 active=[r for r in coefficients if r[2]];guards=coefficients if full else active[:2]+active[-2:]
 rational=[{m:F(v,den) for m,v in f.items()} for f in fields]
 for a,m,q in guards:assert independent.quadratic_coefficient(rational,a,m)==q
 if matrix_checker:
  check=matrix_checker.coefficients(t,receiver)
  for a,m,q in coefficients:assert check[a,m]==matrix_checker.scalar(q)
 return answer,len(guards)

@lru_cache(None)
def gram(mask):return prod(G[i] for i in BC.bit_indices(mask))
def interior(mask,i):
 if mask>>i&1:
  n=mask^(1<<i);return n,BC.sign(1<<i,n)
def act(t,maps):
 out={}
 for (a,m),c in t.items():
  for b,v in maps[a]:BC.plus(out,(b,m),c*v)
  for i in BC.bit_indices(m):
   n,e=interior(m,i)
   for b,v in maps[i]:
    if not n>>b&1:BC.plus(out,(a,n|(1<<b)),c*v*e*BC.sign(1<<b,n))
 return out
def curvature_source(t,R):
 linear=defaultdict(F)
 for (i,j),A in R.items():
  maps=tuple(tuple((b,-c) for b,c in enumerate(row) if c) for row in A)
  for (a,m),coefficient in act(t,maps).items():
   for p,q,ordering in ((i,j,-1),(j,i,1)):
    c=coefficient*ordering/(G[i]*G[j])
    if a==q:linear[p,m]+=2*c/gram(m)
    hit=interior(m,q)
    if hit is None:continue
    n,e=hit
    if not n>>p&1:linear[a,n|(1<<p)]-=c*e*BC.sign(1<<p,n)/(G[a]*gram(n))
 return sum(c*linear.get(key,0) for key,c in t.items())

def main():
 parser=argparse.ArgumentParser();parser.add_argument('--all-independent',action='store_true',help='slow: check every output coefficient again with nested rational Clifford products');parser.add_argument('--matrix-independent',action='store_true',help='check every output coefficient by independent128D matrices modulo65521');parser.add_argument('--receipt');args=parser.parse_args();start=time.perf_counter()
 matrix_checker=None
 if args.matrix_independent:import CK_modular_clifford as matrix_checker
 data=json.loads(gzip.decompress((HERE/'CK_h4_invariant_witness.json.gz').read_bytes()));R=BG.curvature_model();raising=roots();independent=X.BladeAlgebra(G)
 character=json.loads((HERE/'BX_data/h4_character.json').read_text(encoding='utf-8'))
 known={tuple(r['weight']):(r['copies'],r['dimension']) for r in character['components']}
 raw=BQ.receiver(5);assert len(raw)==1143 and all(not BC.action(i,raw) for i in range(10))
 receiver=[(a,m,v/(G[a]*gram(m))) for (a,m),v in raw.items()]
 records=[];lookup={};signature=[0,0];whole=[0,0];copies=0;coefficients=0;guards=0
 for group in data['groups']:
  w=group['weight'];n=group['multiplicity'];tensors=[]
  assert (n,group['per_copy'])==known.pop(tuple(w))
  for v in group['vectors']:
   t={(a,m):F(c,v['denominator']) for a,m,c in v['entries']};assert len(t)==len(v['entries']) and all(m.bit_count()==4 for a,m in t)
   highest(t,w,raising);assert BC.hook_traces(t)==({},{});assert not BG.full_sparse_curvature(R,BC.hodge_scaled(t));tensors.append(t)
  assert len(tensors)==n;copies+=n
  g=s.Matrix(n,n,lambda i,j:s.Rational(BF.inner(tensors[i],tensors[j])));assert g==s.Matrix(group['gram']) and all(g[:i,:i].det()>0 for i in range(1,n+1))
  a=s.zeros(n);polars=0
  for i,t in enumerate(tensors):
   a[i,i],count=value(t,receiver,independent,args.all_independent,matrix_checker);guards+=count;polars+=1
  for i,j in combinations(range(n),2):
   z,count=value(combine(((1,tensors[i]),(1,tensors[j]))),receiver,independent,args.all_independent,matrix_checker);guards+=count;polars+=1;a[i,j]=a[j,i]=(s.Rational(z)-a[i,i]-a[j,j])/2
  coefficients+=1143*polars;assert a==s.Matrix(group['L5'])
  lower,diagonal=a.LDLdecomposition(hermitian=False);assert lower*diagonal*lower.T==a
  pivots=list(diagonal.diagonal());assert all(pivots);inertia=[sum(bool(v>0) for v in pivots),sum(bool(v<0) for v in pivots)]
  # Exact generalized-eigenvalue bounds and a gap on both sides of zero.
  for shifted,expected_positive in ((a+s.Rational(17,3)*g,n),(6*g-a,n),(a-g/100,inertia[0]),(a+g/15,inertia[0])):
   ll,dd=shifted.LDLdecomposition(hermitian=False);assert ll*dd*ll.T==shifted and all(dd.diagonal())
   assert sum(bool(v>0) for v in dd.diagonal())==expected_positive
  transform=lower.T.inv();assert transform.T*a*transform==diagonal
  signature=[signature[i]+group['per_copy']*inertia[i] for i in (0,1)]
  for i in (0,1):
   if inertia[i]==n:whole[i]+=group['per_copy']*n
  assert group['per_copy']==(2*w[0]+3)*(2*w[1]+1)*(w[0]-w[1]+1)*(w[0]+w[1]+2)//6
  record={'weight':w,'multiplicity':n,'per_copy':group['per_copy'],'inertia':inertia,'exact_LDL_pivots':list(map(str,pivots)),'rational_congruence':[[str(v) for v in row] for row in transform.tolist()],'all_receiver_coefficients':1143*polars,'four_exact_spectral_shift_inertias':True}
  records.append(record);lookup[tuple(w)]=(group,tensors,g,a);print(w,inertia,'seconds',round(time.perf_counter()-start,2),flush=True)
 assert not known and len(records)==22 and copies==57 and signature==[5893,4495] and whole==[2097,853]
 frontier=[]
 for item in data['H_L5_frontier']:
  weights=list(map(F,item['weights']));assert all(v>0 for v in weights) and sum(weights)==1
  L5=H=F(0)
  for weight,ray in zip(weights,item['rays']):
   group,tensors,g,a=lookup[tuple(ray['weight'])];i=ray['copy'];norm=F(g[i,i]);source=curvature_source(tensors[i],R)
   assert source==F(group['curvature_form'][i][i])
   assert source/norm==F(ray['H']) and F(a[i,i])/norm==F(ray['L5'])
   H+=weight*source/norm;L5+=weight*F(a[i,i])/norm
  assert L5==0 and H==F(item['H']);assert item['rays'][0]['weight']!=item['rays'][1]['weight'];frontier.append({'norm_squared':'1','L5':'0','H':str(H)})
 assert F(frontier[0]['H'])<0<F(frontier[1]['H'])
 assert F(1,100)/(F(1,100)+F(17,3))==F(3,1703)
 assert F(1,15)/(F(1,15)+6)==F(1,91)
 result={'exact':True,'types':22,'copies':57,'dimension':10388,'L5_signature_positive_negative_zero':signature+[0],'maximal_whole_isotypic_definite_dimensions':whole,'positive_spectrum_open_interval':['1/100','6'],'negative_spectrum_open_interval':['-17/3','-1/15'],'strict_energy_lower_bounds':{'negative_spectral_subspace':'3/1703','positive_spectral_subspace':'1/91'},'all_actual_receiver_coefficients':coefficients,'independent_nested_coefficient_checks':guards,'independent_modular_matrix_coefficients':coefficients if matrix_checker else 0,'blocks':records,'actual_real_zero_L5_H_frontier':frontier,'seconds':time.perf_counter()-start}
 if args.receipt:Path(args.receipt).write_text(json.dumps(result,indent=2))
 print(json.dumps({k:v for k,v in result.items() if k!='blocks'},indent=2))
if __name__=='__main__':main()
