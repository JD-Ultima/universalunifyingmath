"""Independent original128-spin nested-product guards for high Cartan Gram data."""
import DN_clifford as H
from fractions import Fraction as F
import numpy as np
import sympy as s
import json,time
def guard(receiver,blocks,certificate):
 start=time.perf_counter();d={b['tag']:np.array(b['numerator'],dtype=object) for b in certificate['blocks']};p=int(certificate['primes'][0]);den=int(certificate['denominator'])%p;roots=[int(s.sqrt_mod(H.residue(1/g,p),p)) for g in H.G];rr=np.arange(128);par=np.array([[1-2*((z&r).bit_count()%2) for r in range(128)] for z in range(128)])
 def mod(x):return np.remainder(x.real,p)+1j*np.remainder(x.imag,p)
 def matrix(v):
  out=np.zeros((14,128,128),complex)
  for (a,m),v in v.items():
   x,z,phase=H.previous.blade(m);coef=H.residue(v,p)*roots[a]*H.prod(roots[i] for i in range(14) if m>>i&1)%p;out[a,rr^x,rr]+=coef*phase*par[z]
  return mod(out)
 gamma=H.previous.GAMMA.copy();Y=matrix(receiver);checks=[]
 for tag,selected in [('4_0',[0]),('4_1',[0,1])]:
  b=next(b for b in blocks if f"{b['w1']}_{b['SO3_spin']}"==tag);v={}
  for j in selected:
   for key,c in b['vectors'][j].items():v[key]=v.get(key,F(0))+c
  C=matrix(v);S=mod(np.sum(mod(C@C),axis=0));Q=[]
  for a in range(14):
   U=mod(sum(mod(mod(gamma[b]@C[a])@C[b]) for b in range(14)));V=mod(sum(mod(mod(C[b]@C[a])@gamma[b]) for b in range(14)));W=mod(sum(mod(mod(C[b]@gamma[a])@C[b]) for b in range(14)));Q.append(mod(gamma[a]@S+S@gamma[a]-U-V-W))
  Q=np.array(Q);pair=int(np.sum(Y.real.astype(np.int64)*Q.real.astype(np.int64)+Y.imag.astype(np.int64)*Q.imag.astype(np.int64)))%p;actual=(-2*pow(3,-1,p)*pow(128,-1,p)*pair)%p;expected=int(sum(int(d[tag][i,j]) for i in selected for j in selected))*pow(den,-1,p)%p;assert actual==expected,(tag,actual,expected);checks.append({'tag':tag,'indices':selected,'residue':actual})
 return {'exact_modular_original128':True,'prime':p,'checks':checks,'seconds':time.perf_counter()-start}
