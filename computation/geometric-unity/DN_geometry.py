"""Independent compact-frame rational geometry audit for the high Cartan sector."""
import os,sys
os.environ['OPENBLAS_NUM_THREADS']='1';sys.dont_write_bytecode=True
from pathlib import Path
from fractions import Fraction as F
from collections import defaultdict,Counter
import gzip,json,time
import sympy as s
HERE=Path(__file__).resolve().parent
import BC_affine_hooks as BC
import BG_h6_curvature as BG
PUBLIC=Path(BG.__file__).resolve().parent
from CN_h4_s5_cubic import weight_multiplicities
G=tuple(map(F,BC.METRIC));phi=[s.Matrix(x)/5 for x in BG.integer_phi()]
def operator(a,b):
 J=s.zeros(5);J[a,b]=1;J[b,a]=-1;A=s.Matrix(14,14,lambda i,j:s.trace(phi[i]*(J*phi[j]-phi[j]*J))/s.Rational(G[i]));D=-A.T
 assert A.T*s.diag(*map(s.Rational,G))+s.diag(*map(s.Rational,G))*A==s.zeros(14)
 return [[(i,F(D[i,j])) for i in range(14) if D[i,j]] for j in range(14)]
def act(v,D):
 out=defaultdict(F)
 for (a,m),c in v.items():
  for b,x in D[a]:out[b,m]+=c*x
  for i in BC.bit_indices(m):
   tail=m^(1<<i);sgn=(-1)**((m&((1<<i)-1)).bit_count())
   for b,x in D[i]:
    if not tail>>b&1:out[a,tail|(1<<b)]+=c*x*sgn*(-1)**((tail&((1<<b)-1)).bit_count())
 return {k:v for k,v in out.items() if v}
def combine(*pairs):
 out=defaultdict(F)
 for c,v in pairs:
  for k,x in v.items():out[k]+=c*x
 return {k:v for k,v in out.items() if v}
def audit(basis,rec,full=True):
 start=time.perf_counter();D01,D23,D24,D34=[operator(*ij) for ij in [(0,1),(2,3),(2,4),(3,4)]]
 receiver={(a,m):F(v) for a,m,v in rec['receiver']}
 for D in (D01,D23,D24,D34):assert not act(receiver,D)
 assert sum(receiver.get((a,1<<a),0)/G[a] for a in range(14))==0
 R=BG.curvature_model();checks=[]
 for block in basis['blocks']:
  w,l=block['w1'],block['SO3_spin']
  for j,v in enumerate(block['highest_vectors']):
   re={(a,m):F(x) for a,m,x in v['real']};im={(a,m):F(x) for a,m,x in v['imag']}
   for tensor in (re,im):
    assert BC.hook_traces(tensor)==({},{})
    if full or j==0:assert not BG.full_sparse_curvature(R,BC.hodge_scaled(tensor))
   assert act(re,D01)==combine((w,im)) and act(im,D01)==combine((-w,re))
   assert act(re,D23)==combine((l,im)) and act(im,D23)==combine((-l,re))
   assert not combine((1,act(re,D24)),(1,act(im,D34)))
   assert not combine((1,act(im,D24)),(-1,act(re,D34)))
  checks.append({'weight':w,'spin':l,'multiplicity':len(block['highest_vectors'])});print('exact geometry',w,l,flush=True)
 inventory=json.loads((PUBLIC/'BX_data/h4_character.json').read_text())['components'];counts=Counter()
 for component in inventory:
  for weight,mult in weight_multiplicities(*component['weight']).items():counts[weight]+=component['copies']*mult
 expected={(w,l):counts[w,l]-counts[w,l+1] for w,l in counts if w>=4 and l>=0 and counts[w,l]>counts[w,l+1]};actual={(b['w1'],b['SO3_spin']):len(b['highest_vectors']) for b in basis['blocks']};assert expected==actual;dim=sum(m*(2*l+1) for (w,l),m in actual.items());assert dim==480
 return {'exact':True,'real_dimension':2*dim,'highest_columns':sum(actual.values()),'all_four_receiver_Lie_actions_zero':True,'receiver_trace_zero':True,'source_hooks_curvature_weights_raising':True,'complete_by_exact_CA_character':True,'blocks':checks,'seconds':time.perf_counter()-start}
