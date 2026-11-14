"""CW: exact positive all-grade covariance on the full68D S5-fixed H4 kernel.

Default checks geometry, exact rounding, positivity and correction bounds,
and six representative integer Clifford values. --full recomputes all363.
"""
import os,sys
os.environ['OPENBLAS_NUM_THREADS']='1';sys.dont_write_bytecode=True
from pathlib import Path
from fractions import Fraction as F
from itertools import combinations_with_replacement
from math import prod,isqrt,lcm
import argparse,gzip,json,time,hashlib
import numpy as np
import CW_geometry as G
from CW_integer_clifford import Engine
HERE=Path(__file__).resolve().parent

def exact_pd(A):
 A=np.array(A,dtype=object);old=1
 for k in range(len(A)):
  pivot=int(A[k,k]);assert pivot>0,(k,pivot)
  if k+1<len(A):
   T=A[k+1:,k+1:]*pivot-A[k+1:,k,None]*A[k,None,k+1:]
   assert all(int(v)%old==0 for v in T.flat)
   A[k+1:,k+1:]=T//old
  old=pivot
 return {'positive':True,'dimension':len(A),'last_minor_sha256':hashlib.sha256(str(pivot).encode()).hexdigest()}

def prepare(w,E):
 p=w['preparation'];n=68;S=p['common_source_power'];powers=p['source_column_powers'];cp=p['basis_coefficient_power'];lp=p['covariance_factor_power'];dp=p['unit_power'];ip=80
 assert (cp,lp,dp,p['rightinverse_power'])==(48,40,36,36)
 masks=E.obj.frame.masks;location=E.obj.frame.location;V=np.zeros((14014,n),dtype=object)
 for j,t in enumerate(w['tensors']):
  for a,m,v in t:V[a*1001+location[m],j]=int(v)*(1<<(S-powers[j]))
 C=np.array(p['basis_coefficients'],dtype=object);L=np.array(p['covariance_factor'],dtype=object);B=V@C;Fmat=B@L
 metrics=[F(G.BC.METRIC[a])*prod(F(G.BC.METRIC[i]) for i in range(14) if m>>i&1) for a in range(14) for m in masks]
 inv=[isqrt((g.denominator<<(2*ip))//g.numerator) for g in metrics]
 def rounded(A,power,label):
  shift=power+ip-dp;denpow=power+ip;D=np.empty(A.shape,dtype=np.int64);errors=[0]*n
  for i,row in enumerate(A):
   lo=row*inv[i];hi=row*(inv[i]+1);di=(lo+(1<<(shift-1)))//(1<<shift);D[i]=np.asarray(di,dtype=np.int64)
   error=np.maximum(abs(lo-di*(1<<shift)),abs(hi-di*(1<<shift)))
   for j,v in enumerate(error):errors[j]+=int(v)**2
  nums=[isqrt(v)+1 for v in errors];norms=[sum(int(v)**2 for v in D[:,j]) for j in range(n)];rec=p[label+'_rounding']
  assert nums==list(map(int,rec['rounding_error_numerators'])) and norms==list(map(int,rec['norm_squared_integer_sums'])) and denpow==rec['rounding_error_denominator_power']
  assert max(nums)*2*10**9<1<<denpow
  return D,norms
 D,dnorms=rounded(B,S+cp,'basis');DF,fnorms=rounded(Fmat,S+cp+lp,'covariance')
 assert max(dnorms)*10000<11881*(1<<(2*dp));assert sum(fnorms)<80*(1<<(2*dp))
 M=L@L.T;proof=exact_pd(M-np.eye(n,dtype=object)*(1<<(2*lp-2)))
 return D,DF,proof

def correction(w,A,v):
 p=w['preparation'];Y=np.array(p['rightinverse'],dtype=object);den=21*(1<<(2*p['unit_power']));ys=1<<p['rightinverse_power'];n=227;root=lambda x:isqrt(int(x))+1
 error=den*ys*np.eye(n,dtype=object)-A@Y
 one=max(sum(abs(int(x)) for x in error[:,j]) for j in range(n));inf=max(sum(abs(int(x)) for x in error[i,:]) for i in range(n));e=F(root(one*inf),den*ys)
 ynorm=F(root(sum(int(x)**2 for x in Y.flat)),ys);assert ynorm<450
 delta=F(1,2*10**9);neumann=e+n*880*delta*ynorm;assert neumann<F(1,20)
 residual_d=F(root(sum(int(x)**2 for x in v)),den);residual_b=residual_d+200*delta*(148+68*delta)
 correction=2*ynorm*residual_b/(1-neumann);assert correction<F(1,64)
 return {'exact':True,'dimension':68,'output_rank':227,'basis_and_ray_rounding_bound':str(delta),'dyadic_inverse_defect_bound':str(e),'exact_basis_neumann_defect_bound':str(neumann),'dyadic_covariance_residual_bound':str(residual_d),'moment_correction_operator_bound':str(correction),'final_moment_lower_bound':'15/64','numeric_readout':{'neumann':float(neumann),'dyadic_residual':float(residual_d),'correction':float(correction)}}

def original_guard(E):
 # Independent sparse Clifford multiplication from BC, on a generic exact hook.
 import CD_clifford_hook_identities as CD
 B=CD.hook_sample(4,7814,40);den=lcm(*(v.denominator for b in B for v in b.values()));D=np.zeros((14,1001),dtype=np.int64)
 for a,b in enumerate(B):
  for m,v in b.items():D[a,E.obj.frame.location[m]]=int(v*den)
 field=[{m:int(v*den) for m,v in b.items()} for b in B];Q=G.BC.Clifford([1]*14).quadratic(field);norm=sum(int(v)**2 for v in D.ravel());J=E.value(D)
 for (a,m),j in zip(E.receivers,J):assert j==21*Q[a].get(m,0)+(18*norm if m==1<<a else 0)
 return len(J)

def main():
 parser=argparse.ArgumentParser();parser.add_argument('--full',action='store_true');parser.add_argument('--receipt',type=Path);args=parser.parse_args();start=time.perf_counter()
 w=json.loads(gzip.decompress((HERE/'CW_h4_s5_moment_witness.json.gz').read_bytes()));assert w['dimension']==68 and len(w['receivers'])==227 and len(set(map(tuple,w['receivers'])))==227
 geometry=G.check_basis(w);counts=G.output_count();print('exact geometry passed',round(time.perf_counter()-start,2),flush=True)
 E=Engine(w['receivers']);D,DF,positive=prepare(w,E);guard=original_guard(E);cache=w['integer_values'];pairs=list(combinations_with_replacement(range(68),2));cols=w['preparation']['selected_quadratic_columns'];assert len(cols)==227 and len(set(cols))==227
 keys=['b'+str(i) for i in range(68)]+['p'+str(c) for c in cols if pairs[c][0]!=pairs[c][1]]+['f'+str(i) for i in range(68)];assert set(keys)==set(cache)
 pkeys=[k for k in keys if k[0]=='p'];chosen=keys if args.full else ['b0','b67',pkeys[0],pkeys[-1],'f0','f67']
 for i,key in enumerate(chosen):
  j=int(key[1:])
  if key[0]=='b':x=D[:,j]
  elif key[0]=='f':x=DF[:,j]
  else:a,b=pairs[j];x=D[:,a]+D[:,b]
  assert E.value(x)==list(map(int,cache[key])),key
  if (i+1)%40==0:print('exact Clifford values',i+1,'/',len(chosen),'seconds',round(time.perf_counter()-start,2),flush=True)
 val=lambda key:np.array(list(map(int,cache[key])),dtype=object)
 A=[]
 for col in cols:
  i,j=pairs[col];A.append(val('b'+str(i)) if i==j else val('p'+str(col))-val('b'+str(i))-val('b'+str(j)))
 A=np.array(A,dtype=object).T;v=sum((val('f'+str(i)) for i in range(68)),np.zeros(227,dtype=object));result=correction(w,A,v)
 result.update({'geometry':geometry,'output_counts':counts,'initial_covariance_pd':positive,'independent_original_Clifford_guards':guard,'recomputed_integer_values':len(chosen),'total_integer_values':len(keys),'full':args.full,'seconds':time.perf_counter()-start})
 if args.receipt:args.receipt.write_text(json.dumps(result,indent=2))
 print(json.dumps(result,indent=2))
if __name__=='__main__':main()
