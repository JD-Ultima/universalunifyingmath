"""DC: a certified real algebraic four-atom all-grade H4 covariance.

Default: all exact bounds/geometry, six sample guards and six quartic pairs.
--full: independently regenerate every cached integer Clifford value.
"""
import os,sys
os.environ['OPENBLAS_NUM_THREADS']='1';sys.dont_write_bytecode=True
from pathlib import Path
from fractions import Fraction as F
from math import prod,isqrt
import argparse,json,gzip,time,hashlib
import numpy as np
import CW_geometry as G
import CW_h4_s5_moment as CW
from DC_integer_clifford import Engine
from DC_rank4_quartic import verify as quartic_verify
HERE=Path(__file__).resolve().parent
def source_hash(w):
 p=w['preparation'];a={'tensors':w['tensors'],'C':p['basis_coefficients'],'scales':p['source_column_powers'],'receivers':w['receivers']}
 return hashlib.sha256(json.dumps(a,separators=(',',':')).encode()).hexdigest()
def prepare(w,p,E):
 old=w['preparation'];S=old['common_source_power'];cp=old['basis_coefficient_power'];xp=p['factor_power'];dp=p['unit_power'];ip=80;n=68
 assert (cp,xp,dp,p['rightinverse_power'])==(48,48,56,36)
 V=np.zeros((14014,n),dtype=object)
 for j,t in enumerate(w['tensors']):
  for a,m,v in t:V[a*1001+E.obj.frame.location[m],j]=int(v)*(1<<(S-old['source_column_powers'][j]))
 B=V@np.array(old['basis_coefficients'],dtype=object);Xi=np.array(p['factor'],dtype=object)
 metrics=[F(G.BC.METRIC[a])*prod(F(G.BC.METRIC[i]) for i in range(14) if m>>i&1) for a in range(14) for m in E.obj.frame.masks];inv=[isqrt((g.denominator<<(2*ip))//g.numerator) for g in metrics]
 def rounded(A,power,label):
  cols=A.shape[1];shift=power+ip-dp;denpow=power+ip;D=np.empty(A.shape,dtype=np.int64);errors=[0]*cols
  for i,row in enumerate(A):
   lo=row*inv[i];hi=row*(inv[i]+1);di=(lo+(1<<(shift-1)))//(1<<shift);D[i]=np.asarray(di,dtype=np.int64);err=np.maximum(abs(lo-di*(1<<shift)),abs(hi-di*(1<<shift)))
   for j,e in enumerate(err):errors[j]+=int(e)**2
  nums=[isqrt(v)+1 for v in errors];norms=[sum(int(v)**2 for v in D[:,j]) for j in range(cols)];r=p[label+'_rounding']
  assert nums==list(map(int,r['numerators'])) and norms==list(map(int,r['norm_squared_integer_sums'])) and denpow==r['denominator_power'];assert max(nums)*10**15<1<<denpow
  return D,norms
 DB,bnorm=rounded(B,S+cp,'basis');DF,fnorm=rounded(B@Xi,S+cp+xp,'anchor')
 assert max(bnorm)*10000<11881*(1<<(2*dp));assert max(fnorm)*100<81*(1<<(2*dp))
 factorproof=CW.exact_pd(5*(Xi.T@Xi)-np.eye(4,dtype=object)*(1<<(2*xp)))
 # Exact14014-term Gram products via14-bit limbs, bounded below2^53.
 assert 14014*((1<<14)-1)**2<2**53
 base=1<<14;Z=DB.copy();limbs=[]
 for k in range(3):lo=Z%base;limbs.append(lo.astype(float));Z=(Z-lo)//base
 assert np.max(abs(Z))<base;limbs.append(Z.astype(float));Gram=np.zeros((68,68),dtype=object)
 for i,L in enumerate(limbs):
  for j,R in enumerate(limbs):
   dot=L.T@R;assert np.all(dot==np.rint(dot));Gram+=dot.astype(np.int64).astype(object)*(1<<(14*(i+j)))
 for i,j in [(0,0),(0,67),(35,35)]:assert sum(int(a)*int(b) for a,b in zip(DB[:,i],DB[:,j]))==Gram[i,j]
 gramproof=CW.exact_pd(3*(1<<(2*dp))*np.eye(68,dtype=object)-2*Gram)
 return DB,DF,{'factor_Gram_lower':'1/5','factor_Gram_proof':factorproof,'basis_dyadic_Gram_upper':'3/2','basis_Gram_proof':gramproof}
def contraction(p,J,v):
 Y=np.array(p['rightinverse'],dtype=object);den=21*(1<<(2*p['unit_power']));ys=1<<p['rightinverse_power'];root=lambda v:isqrt(int(v))+1
 E=den*ys*np.eye(228,dtype=object)-Y@J
 one=max(sum(abs(int(x)) for x in E[:,j]) for j in range(228));inf=max(sum(abs(int(x)) for x in E[i,:]) for i in range(228));defect=F(root(one*inf),den*ys);ynorm=F(root(sum(int(x)**2 for x in Y.flat)),ys);assert ynorm<1800
 delta=F(1,10**15);true_defect=defect+ynorm*16*840*delta;residual_d=F(root(sum(int(x)**2 for x in v)),den);residual=residual_d+1600*delta;beta=ynorm*residual
 rho=F(1,10**8);constant=true_defect+802*ynorm*rho
 assert true_defect<F(1,1000) and beta<rho/2 and constant<F(1,50)
 return {'exact':True,'factor_rank':4,'affine_dimension':228,'domain_dimension':68,'center_displacement_radius':str(rho),'basis_and_anchor_rounding_bound':str(delta),'inverse_Frobenius_bound':str(ynorm),'dyadic_inverse_defect':str(defect),'exact_center_inverse_defect':str(true_defect),'dyadic_center_residual':str(residual_d),'exact_center_residual_bound':str(residual),'Newton_displacement_at_zero_bound':str(beta),'contraction_bound':str(constant),'factor_minimum_singular_value_lower':'2/5','numeric_readout':{'inverse_norm':float(ynorm),'center_residual':float(residual),'displacement_at_zero':float(beta),'inverse_defect':float(true_defect),'contraction':float(constant)}}
def main():
 parser=argparse.ArgumentParser();parser.add_argument('--full',action='store_true');parser.add_argument('--receipt',type=Path);args=parser.parse_args();start=time.perf_counter()
 w=json.loads(gzip.decompress((Path(CW.__file__).resolve().parent/'CW_h4_s5_moment_witness.json.gz').read_bytes()));bundle=json.loads(gzip.decompress((HERE/'DC_h4_rank4_witness.json.gz').read_bytes()));assert source_hash(w)==bundle['CW_source_sha256'];p=bundle['preparation'];cache=bundle['integer_values'];geometry=G.check_basis(w);assert G.output_count()['total_homogeneous_output_bound']==227
 E=Engine(w['receivers']);DB,DF,positive=prepare(w,p,E);guards=CW.original_guard(E);print('exact geometry and preparation passed',round(time.perf_counter()-start,2),flush=True);Xi=np.array(p['factor'],dtype=object);cols=p['derivative_columns'];assert len(cols)==len(set(cols))==228 and all(0<=c<272 for c in cols);den=21*(1<<(2*p['unit_power']))
 keys=['f'+str(i) for i in range(4)]+['b'+str(i) for i in sorted({c//4 for c in cols})]+['p'+str(c) for c in cols];assert set(keys)==set(cache);ps=[k for k in keys if k[0]=='p'];bs=[k for k in keys if k[0]=='b'];chosen=keys if args.full else ['f0','f3',bs[0],bs[-1],ps[0],ps[-1]]
 for count,key in enumerate(chosen):
  j=int(key[1:])
  if key[0]=='f':D=DF[:,j]
  elif key[0]=='b':D=DB[:,j]
  else:i,r=divmod(j,4);D=DF[:,r]+DB[:,i]
  assert E.value(D)==list(map(int,cache[key])),key
  if (count+1)%30==0:print('integer values',count+1,'/',len(chosen),'seconds',round(time.perf_counter()-start,2),flush=True)
 val=lambda key:np.array(list(map(int,cache[key])),dtype=object);J=[]
 for c in cols:
  i,r=divmod(c,4);top=val('p'+str(c))-val('f'+str(r))-val('b'+str(i));J.append(np.r_[top,2*int(Xi[i,r])*(den//(1<<p['factor_power']))])
 J=np.array(J,dtype=object).T;v=np.r_[sum((val('f'+str(i)) for i in range(4)),np.zeros(227,dtype=object)),(sum(int(x)**2 for x in Xi.flat)-(1<<(2*p['factor_power'])))*(den//(1<<(2*p['factor_power'])))];result=contraction(p,J,v);result.update({'full':args.full,'geometry':geometry,'positive_Gram_bounds':positive,'original_Clifford_guards':guards,'recomputed_integer_values':len(chosen),'total_integer_values':len(keys),'seconds':time.perf_counter()-start})
 result['quartic_span']=quartic_verify(p,bundle['quartic'],cache,DF,E)
 added=result['quartic_span']['additional_exact_anchor_pair_values']
 result['recomputed_integer_values']+=added
 result['total_integer_values']+=added
 assert result['total_integer_values']==306
 result['seconds']=time.perf_counter()-start
 if args.receipt:args.receipt.write_text(json.dumps(result,indent=2))
 print(json.dumps(result,indent=2))
if __name__=='__main__':main()
