"""DT: exact real exclusion on the coherent H4 Cartan (3,3) slice.

Uses published DN Clifford arithmetic and compact geometry helpers.
Default checks one complete modular image; --full checks all ten.
"""
import os,sys
os.environ['OPENBLAS_NUM_THREADS']='1';sys.dont_write_bytecode=True
from pathlib import Path
from fractions import Fraction as F
from math import lcm,prod,isqrt
import argparse,gzip,json,time
import numpy as np
import DN_clifford as H
import DN_geometry as G
from DT_original_guards import guard
HERE=Path(__file__).resolve().parent
def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--full',action='store_true');ap.add_argument('--receipt');args=ap.parse_args();start=time.perf_counter()
 with gzip.open(HERE/'DT_coherent_witness.json.gz','rt') as f:data=json.load(f)
 certificate=data['certificate'];b=data['basis']['blocks'][0];assert (b['w1'],b['SO3_spin'],b['multiplicity'])==(3,3,28);receiver={(a,m):F(v) for a,m,v in certificate['receiver']};D01,D23,D24,D34=[G.operator(*ij) for ij in [(0,1),(2,3),(2,4),(3,4)]]
 assert all(m.bit_count() in (1,5,9) and H.parity(a,m)==1 for a,m in receiver);assert not G.act(receiver,D01) and not G.act(receiver,D23);assert sum(receiver.get((a,1<<a),0)/H.G[a] for a in range(14))==0
 R=G.BG.curvature_model();vectors=[]
 for v in b['highest_vectors']:
  re={(a,m):F(x) for a,m,x in v['real']};im={(a,m):F(x) for a,m,x in v['imag']};vectors.append(re)
  for tensor in (re,im):assert G.BC.hook_traces(tensor)==({},{}) and not G.BG.full_sparse_curvature(R,G.BC.hodge_scaled(tensor))
  assert all(H.parity(a,m)==1 for a,m in re) and all(H.parity(a,m)==-1 for a,m in im)
  for D in (D01,D23):assert G.act(re,D)==G.combine((3,im)) and G.act(im,D)==G.combine((-3,re))
  assert not G.combine((1,G.act(re,D24)),(1,G.act(im,D34))) and not G.combine((1,G.act(im,D24)),(-1,G.act(re,D34)))
 inventory=json.loads((G.PUBLIC/'BX_data/h4_character.json').read_text())['components'];dimension=0
 for c in inventory:
  weights=G.weight_multiplicities(*c['weight']);dimension+=c['copies']*(weights.get((3,3),0)-weights.get((3,4),0))
 assert dimension==28;b['vectors']=vectors
 dy=lcm(*(v.denominator for v in receiver.values()));db=lcm(*(v.denominator for vec in vectors for v in vec.values()));den=3*dy*db*db*prod(g.numerator**3 for g in H.G);assert den==certificate['denominator'];ynorm=sum(v*v/H.metric(a,m) for (a,m),v in receiver.items());Bnorm=max(sum(v*v/H.metric(a,m) for (a,m),v in vec.items()) for vec in vectors);bound=150*(isqrt(ynorm.numerator//ynorm.denominator)+1)*(Bnorm.numerator//Bnorm.denominator+1)*den;assert bound==certificate['absolute_entry_bound']
 primes=certificate['primes'];modulus=prod(primes);assert modulus==certificate['modulus'] and modulus>2*bound and len(primes)==len(set(primes))
 for p in primes:assert p%120==1 and all(p%d for d in range(2,isqrt(p)+1)) and 2*896*(p-1)**2<2**53
 record=certificate['blocks'][0];assert record['tag']=='3_3';A=np.array(record['numerator'],dtype=object);assert np.array_equal(A,A.T) and max(abs(int(x)) for x in A.flat)<bound;assert H.exact_pd(A)==record['positive']
 Gram=[[sum(x*vectors[j].get(slot,F(0))/H.metric(*slot) for slot,x in vectors[i].items()) for j in range(28)] for i in range(28)];gd=lcm(*(v.denominator for row in Gram for v in row));Gnum=np.array([[int(v*gd) for v in row] for row in Gram],dtype=object);assert H.exact_pd(20*gd*A-19*den*Gnum)==certificate['gram_gap_certificate'];assert ynorm<F(9,2)
 original=guard(receiver,[b],certificate);checked=[]
 for p in (primes if args.full else primes[:1]):
  actual=H.modular(receiver,[b],p)['3_3'];assert np.array_equal((actual.astype(object)*(den%p)%p).astype(np.int64),(A%p).astype(np.int64));checked.append(p);print('verified modular image',len(checked),p,flush=True)
 rec={'exact':True,'full':args.full,'complex_highest_dimension':28,'real_coherent_dimension':56,'source_hooks_curvature_two_weights_and_raising':True,'receiver_two_torus_invariant_tracefree_and_conjugation_even':True,'source_dimension_complete_by_CA_character':True,'Gram_gap':'19/20','receiver_norm_squared':str(ynorm),'full_residual_squared_gap':'1/5','integer_entry_bound_bits':bound.bit_length(),'CRT_modulus_bits':modulus.bit_length(),'verified_primes':checked,'all_CRT_primes':primes,'original128_guards':original,'seconds':time.perf_counter()-start}
 if args.receipt:Path(args.receipt).write_text(json.dumps(rec,indent=2))
 print(json.dumps(rec,indent=2))
if __name__=='__main__':main()
