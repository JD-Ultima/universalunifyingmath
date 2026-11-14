"""Exact invariant-output counts and positive Haar moments in H4 and H6.

Verifies two small curvature-compatible H4 tensors with opposite invariant
grade-five values. Their positively weighted Haar covariances annihilate every
homogeneous necessary quadratic output. Hence no linear combination of those
quadratics can be positive definite on the full H4 curvature kernel.
The exact coefficients from BQ give the analogous H6 certificate
on the submodule already admitting closed fields.
"""
from collections import Counter,defaultdict
from fractions import Fraction as F
from itertools import combinations,combinations_with_replacement,permutations,product
from math import gcd,lcm,prod
from pathlib import Path
import argparse,json,time
import BC_affine_hooks as BC
import BG_h6_curvature as BG
import BO_h6_full_kernel as BO
import BQ_h6_invariant_forms as BQ
import BC_quadratic_xcheck as X
HERE=Path(__file__).resolve().parent
WEIGHTS=[(2,0),(-2,0),(0,2),(0,-2),(1,1),(1,-1),(-1,1),(-1,-1),(1,0),(-1,0),(0,1),(0,-1),(0,0),(0,0)]
def invariant_counts():
 exterior=[Counter() for _ in range(15)];exterior[0][0,0]=1
 for w in WEIGHTS:
  for k in range(14,0,-1):
   for v,n in list(exterior[k-1].items()):exterior[k][v[0]+w[0],v[1]+w[1]]+=n
 terms=[]
 for p in permutations((0,1)):
  for signs in product((-1,1),repeat=2):
   image=[signs[i]*(3,1)[p[i]] for i in range(2)]
   terms.append((((3-image[0])//2,(1-image[1])//2),(1 if p==(0,1) else -1)*signs[0]*signs[1]))
 def multiplicity(character):return sum(sign*character[w] for w,sign in terms)
 def tensor_character(degree):
  tensor=Counter()
  for v,a in Counter(WEIGHTS).items():
   for w,b in exterior[degree].items():tensor[v[0]+w[0],v[1]+w[1]]+=a*b
  return tensor
 symmetric=Counter((WEIGHTS[i][0]+WEIGHTS[j][0],WEIGHTS[i][1]+WEIGHTS[j][1]) for i,j in combinations_with_replacement(range(14),2))
 ext=[multiplicity(c) for c in exterior];assert ext==[1,0,0,0,0,1,0,0,0,1,0,0,0,0,1]
 assert multiplicity(symmetric)==1
 assert multiplicity(tensor_character(5))==multiplicity(tensor_character(9))==1
 return {'method':'Exact B2 Weyl denominator: sum det(w) m(rho-w rho), 2rho=(3,1)',
         'exterior_invariants':ext,'symmetric_tracefree_invariants':0,
         'V_tensor_Lambda5_invariants':1,'V_tensor_Lambda9_invariants':1,
         'complete_effective_H4_output_invariants':1,'complete_effective_H6_output_invariants':2}
def mask_character(mask):
 value=0
 for i in range(14):
  if mask>>i&1:value^=BC.CHARS[i]
 return value
def sparse_rank(rows,columns,target,p=65521):
 pivots={}
 for original in sorted(rows,key=len):
  row={i:int(v)%p for i,v in original.items() if int(v)%p}
  while row:
   j=min(row)
   if j not in pivots:
    inv=pow(row[j],-1,p);pivots[j]={i:v*inv%p for i,v in row.items()};break
   factor=row[j]
   for i,v in pivots[j].items():
    value=(row.get(i,0)-factor*v)%p
    if value:row[i]=value
    elif i in row:del row[i]
  if len(pivots)==target:return target
 return len(pivots)
def lie_count(degree,tensor=False):
 masks=[sum(1<<i for i in inds) for inds in combinations(range(14),degree)]
 columns=[(a,m) for m in masks for a in range(14) if BC.char_of(a,m)==0] if tensor else [m for m in masks if mask_character(m)==0]
 denominator=lcm(*(v.denominator for moves in BC.COV for row in moves for b,v in row));equations=defaultdict(dict)
 for j,column in enumerate(columns):
  for which in range(10):
   action=BC.action_basis(which,*column).items() if tensor else BC.form_action(which,column)
   for target,value in action:equations[which,target][j]=int(denominator*value)
 expected=1 if tensor else 0;rank=sparse_rank(list(equations.values()),len(columns),len(columns)-expected);assert rank==len(columns)-expected
 return {'degree':degree,'tensor_external_slot':tensor,'M_fixed_columns':len(columns),'action_rows':len(equations),'prime':65521,'certified_rank':rank,'invariant_dimension':expected}


def h6_moment(full):
 selected=((3,0),(5,0),(6,3))
 source=json.loads(Path(BO.__file__).with_name('BO_h6_full_kernel_witnesses.json').read_text(encoding='utf-8'))
 vectors={tuple(r['weight']):r['vectors'] for r in source['extremal_vectors']}
 data=json.loads((HERE/'BQ_h6_invariant_forms.json').read_text(encoding='utf-8'))
 blocks={tuple(r['weight']):r for r in data['blocks']}
 ratios=[]
 for typ in selected:
  record=blocks[typ]
  assert record['copies']==len(vectors[typ])==1
  g=F(record['gram'][0][0])
  u=F(record['forms']['5'][0][0]);v=F(record['forms']['9'][0][0])
  ratios.append((u/g,384*v/g))
 if full:
  transform=BQ.BP.split_to_physical_covectors()
  receivers=BQ.weighted_receivers()
  engine=X.BladeAlgebra(BC.METRIC)
  for typ in selected:
   g,u,v=BQ.exact_block(vectors[typ],transform,receivers,engine)
   ratio=ratios[selected.index(typ)]
   assert (F(u[0,0]/g[0,0]),384*F(v[0,0]/g[0,0]))==ratio
 assert ratios==[(F(-355,186),F(184,93)),(F(35,8),F(7,3)),(F(-5,6),F(-3))]
 p=(F(4991,14135),F(3292,14135),F(532,1285))
 assert all(v>0 for v in p) and sum(p)==1
 assert all(sum(p[i]*ratios[i][j] for i in range(3))==0 for j in range(2))
 a,b,c=ratios
 determinant=(b[0]-a[0])*(c[1]-a[1])-(b[1]-a[1])*(c[0]-a[0])
 assert determinant==F(-70675,2232)
 return {'weights':selected,'normalized_moments':[list(map(str,r)) for r in ratios],
         'probabilities':list(map(str,p)),'affine_determinant':str(determinant),
         'full_H6_recalibration':full}


def main():
 parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--full-lie',action='store_true',help='independently prove the1717-column tensor-action rank');parser.add_argument('--all-coefficients',action='store_true',help='cross-check all1143 grade-five receiver coordinates by sorted-word nested arithmetic');parser.add_argument('--receipt');args=parser.parse_args();start=time.perf_counter();counts=invariant_counts();receiver=BQ.receiver(5)
 assert receiver and all(not BC.action(j,receiver) for j in range(10));assert BC.hook_traces(receiver)==({},{});lie=[lie_count(8)]
 receiver9=BQ.BN.invariant_receiver()
 assert receiver9 and all(not BC.action(j,receiver9) for j in range(10))
 if args.full_lie:lie.extend([lie_count(6),lie_count(5,True),lie_count(9,True),lie_count(10),lie_count(12)])
 R=BG.curvature_model();data=json.loads((HERE/'CG_h4_invariant_moment.json').read_text());records=[];algorithm=BC.Clifford(BC.METRIC);independent=X.BladeAlgebra(BC.METRIC)
 for item in data['witnesses']:
  tensor={(a,m):F(v) for a,m,v in item['entries']};assert BC.hook_traces(tensor)==({},{});assert not BG.full_sparse_curvature(R,BC.hodge_scaled(tensor));field=[{} for _ in range(14)]
  for (a,m),v in tensor.items():field[a][m]=v
  output=algorithm.quadratic(field);norm=sum(v*v/(BC.METRIC[a]*prod(BC.METRIC[i] for i in range(14) if m>>i&1)) for (a,m),v in tensor.items());value=sum(v*output[a].get(m,0)/(BC.METRIC[a]*prod(BC.METRIC[i] for i in range(14) if m>>i&1)) for (a,m),v in receiver.items());assert norm==F(item['norm_squared']) and value==F(item['L5'])
  entries=list(receiver);guards=entries if args.all_coefficients else [key for key in entries if output[key[0]].get(key[1],0)][:4]+entries[:2]
  for a,m in guards:assert independent.quadratic_coefficient(field,a,m)==output[a].get(m,0)
  records.append({'entries':len(tensor),'norm_squared':str(norm),'L5':str(value),'normalized_L5':str(value/norm),'full_hooks_and_dual_curvature':True,'independent_quadratic_guards':len(guards)})
 weights=list(map(F,data['haar_covariance_weights']));assert all(w>0 for w in weights);assert sum(w*F(r['norm_squared']) for w,r in zip(weights,records))==1;assert sum(w*F(r['L5']) for w,r in zip(weights,records))==0
 result={'exact':True,'invariant_counts':counts,'independent_lie_counts':lie,
         'H4_witnesses':records,'H4_positive_haar_covariance_weights':[str(w) for w in weights],
         'H6_Haar_moment':h6_moment(args.all_coefficients),
         'unit_trace_zero_quadratic_moments':True,'seconds':time.perf_counter()-start}
 if args.receipt:Path(args.receipt).write_text(json.dumps(result,indent=2))
 print(json.dumps(result,indent=2))
if __name__=='__main__':main()
