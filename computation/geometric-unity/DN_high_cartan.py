"""DN: an exact real high-Cartan exclusion in the full H4 curvature kernel.

Default audits all exact geometry/positivity and one full modular image.
--full recomputes every Gram entry modulo all ten reconstruction primes.
"""
import os,sys
os.environ['OPENBLAS_NUM_THREADS']='1';sys.dont_write_bytecode=True
from pathlib import Path
from fractions import Fraction as F
from math import lcm,prod,isqrt
import argparse,gzip,json,time
import numpy as np
import DN_clifford as H
import DN_geometry as geometry
from DN_original_guards import guard
HERE=Path(__file__).resolve().parent
def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--full',action='store_true');ap.add_argument('--receipt');args=ap.parse_args();start=time.perf_counter()
 with gzip.open(HERE/'DN_high_cartan_witness.json.gz','rt') as f:data=json.load(f)
 basis=data['basis'];certificate=data['certificate'];receiver={(a,m):F(v) for a,m,v in certificate['receiver']};blocks=basis['blocks'];records={b['tag']:b for b in certificate['blocks']}
 assert all(H.parity(a,m)==1 for a,m in receiver)
 for b in blocks:
  b['vectors']=[{(a,m):F(v) for a,m,v in r['real']} for r in b['highest_vectors']];assert all(H.parity(a,m)==1 for v in b['vectors'] for a,m in v)
  assert all(H.parity(a,m)==-1 for r in b['highest_vectors'] for a,m,v in r['imag'])
 geom=geometry.audit(basis,certificate,full=True)
 dy=lcm(*(v.denominator for v in receiver.values()));db=lcm(*(v.denominator for b in blocks for v in b['vectors'] for v in v.values()));den=3*dy*db*db*prod(g.numerator**3 for g in H.G);assert den==certificate['denominator']
 ynorm=sum(v*v/H.metric(a,m) for (a,m),v in receiver.items());Bnorm=max(sum(v*v/H.metric(a,m) for (a,m),v in vec.items()) for b in blocks for vec in b['vectors']);bound=150*(isqrt(ynorm.numerator//ynorm.denominator)+1)*(Bnorm.numerator//Bnorm.denominator+1)*den;assert bound==certificate['absolute_entry_bound']
 primes=certificate['primes'];modulus=prod(primes);assert modulus==certificate['modulus'] and modulus>2*bound and len(primes)==len(set(primes))
 for p in primes:assert p%120==1 and all(p%d for d in range(2,isqrt(p)+1)) and 2*896*(p-1)**2<2**53
 positivity=[]
 for b in blocks:
  tag=f"{b['w1']}_{b['SO3_spin']}";A=np.array(records[tag]['numerator'],dtype=object);assert np.array_equal(A,A.T) and max(abs(int(v)) for v in A.flat)<bound;proof=H.exact_pd(A);assert proof==records[tag]['positive']
  vectors=b['vectors'];Gram=[[sum(x*vectors[j].get(slot,F(0))/H.metric(*slot) for slot,x in vectors[i].items()) for j in range(len(vectors))] for i in range(len(vectors))];gd=lcm(*(v.denominator for row in Gram for v in row));Gnum=np.array([[int(v*gd) for v in row] for row in Gram],dtype=object)
  lower=H.exact_pd(20*gd*A-den*Gnum);positivity.append({'tag':tag,'dimension':len(vectors),'strict_positive_receiver':proof,'receiver_minus_norm_over20_positive':lower})
 original=guard(receiver,blocks,certificate);checked=[]
 for p in (primes if args.full else primes[:1]):
  actual=H.modular(receiver,blocks,p)
  for tag,A in actual.items():assert np.array_equal((A.astype(object)*(den%p)%p).astype(np.int64),(np.array(records[tag]['numerator'],dtype=object)%p).astype(np.int64)),(p,tag)
  checked.append(p);print('verified complete modular image',len(checked),p,flush=True)
 assert ynorm<F(3,400)
 rec={'exact':True,'full':args.full,'real_domain_dimension':960,'complete_complex_positive_weight_dimension':480,'multiplicity_blocks':len(blocks),'highest_columns':sum(len(b['vectors']) for b in blocks),'geometry':geom,'positive_blocks':positivity,'exact_receiver_norm_squared':str(ynorm),'homogeneous_receiver_lower_bound':'1/20','complete_residual_squared_lower_bound':str(F(1,400)/ynorm),'clean_complete_residual_squared_lower_bound':'1/3','integer_entry_bound_bits':bound.bit_length(),'CRT_modulus_bits':modulus.bit_length(),'recomputed_modular_images':checked,'all_stored_modular_images':primes,'original128_guards':original,'seconds':time.perf_counter()-start}
 if args.receipt:Path(args.receipt).write_text(json.dumps(rec,indent=2))
 print(json.dumps(rec,indent=2))
if __name__=='__main__':main()
