"""Exact raw-rational high-Cartan receiver matrices by modular Clifford adjoints."""
import os,sys
os.environ['OPENBLAS_NUM_THREADS']='1';sys.dont_write_bytecode=True
from pathlib import Path
from fractions import Fraction as F
from itertools import combinations
from math import prod,lcm,isqrt
import gzip,json,time,argparse
import numpy as np
import sympy as s
HERE=Path(__file__).resolve().parent
import BC_affine_hooks as BC
import BG_h6_curvature as BG
from CW_integer_clifford import Clifford,previous
G=tuple(map(F,BC.METRIC));N=128;M4=[sum(1<<i for i in c) for c in combinations(range(14),4)];loc={m:i for i,m in enumerate(M4)};cc=Clifford(M4,dual=False,target=6/7)
SG=[];sg=np.array([1,-1,1,-1,1]);phi=np.array(BG.integer_phi())
for P in phi:
 plus=bool(np.all(P*sg[:,None]*sg[None,:]==P));minus=bool(np.all(P*sg[:,None]*sg[None,:]==-P));assert plus!=minus;SG.append(1 if plus else -1)
def parity(a,m):return SG[a]*prod(SG[i] for i in range(14) if m>>i&1)
def metric(a,m):return G[a]*prod(G[i] for i in range(14) if m>>i&1)
def field(v,modulus=None,sqrts=None):
 out=np.zeros((14,1001))
 for (a,m),v in v.items():
  if modulus is None:out[a,loc[m]]=float(v)/np.sqrt(float(metric(a,m)))
  else:out[a,loc[m]]=residue(v,modulus)*sqrts[a]*prod(int(sqrts[i]) for i in range(14) if m>>i&1)%modulus
 return out
def residue(v,p):return v.numerator*pow(v.denominator,-1,p)%p
def exact_pd(A):
 import hashlib
 A=np.array(A,dtype=object);old=1
 for k in range(len(A)):
  pivot=int(A[k,k]);assert pivot>0,(k,pivot)
  if k+1<len(A):
   T=A[k+1:,k+1:]*pivot-A[k+1:,k,None]*A[k,None,k+1:];assert all(int(v)%old==0 for v in T.flat);A[k+1:,k+1:]=T//old
  old=pivot
 return {'positive':True,'dimension':len(A),'last_minor_sha256':hashlib.sha256(str(pivot).encode()).hexdigest()}
def receptor(v,modulus=None,sqrts=None):
 grid=np.zeros((14,N,N),complex)
 for (a,m),v in v.items():
  x,z,phase=previous.blade(m)
  coef=float(v)/np.sqrt(float(metric(a,m))) if modulus is None else residue(v,modulus)*sqrts[a]*prod(int(sqrts[i]) for i in range(14) if m>>i&1)%modulus
  grid[a,x,z]=coef*phase
 previous.fwht(grid);full=grid[:,previous.ROWS[:,None]^previous.ROWS[None,:],previous.ROWS[None,:]];return cc.split(full,odd=True)
def core(C,Y,modulus=None):
 def mod(x):
  if modulus is None:return x
  assert np.max(abs(x.real-np.rint(x.real)))==0 and np.max(abs(x.imag-np.rint(x.imag)))==0
  return np.remainder(x.real,modulus)+1j*np.remainder(x.imag,modulus)
 C=mod(C);Y=mod(Y);h=mod(np.sum(cc.left(slice(None),Y)+cc.right(Y,slice(None),odd=True),axis=0));flip=Y[:,::-1];out=[]
 for a in range(14):
  t=mod(cc.sum_product(cc.left(slice(None),Y[a]),C));u=mod(cc.sum_product(cc.right(Y,a,odd=True),C));v=mod(cc.sum_product(cc.left(slice(None),C[a]),flip));out.append(mod(C[a]@h+h@C[a]-cc.anti(t+u+v)))
 return np.array(out)
def coefficients(A,p=None):
 full=cc.merge(A);grid=np.take_along_axis(full,cc.xor[None,...],axis=-2);previous.fwht(grid);out=(grid[:,cc.x,cc.z]*cc.phase.conj()).real
 return out/128 if p is None else (np.remainder(out,p)*pow(128,-1,p))%p
def prime_list():
 n=1_000_081
 while True:
  if s.isprime(n):yield n
  n+=120
def modular(receiver,blocks,p):
 assert p%120==1 and 2*896*(p-1)**2<2**53
 roots=[int(s.sqrt_mod(residue(1/g,p),p)) for g in G];Y=receptor(receiver,p,roots);answers={}
 for b in blocks:
  tag=f"{b['w1']}_{b['SO3_spin']}";D=np.column_stack([field(v,p,roots).ravel() for v in b['vectors']]);cores=np.column_stack([coefficients(core(cc.from_tensor(D[:,j].reshape(14,1001)),Y,p),p).ravel() for j in range(D.shape[1])]);A=np.zeros((D.shape[1],D.shape[1]))
  for lo in range(0,len(D),64):A=np.remainder(A-D[lo:lo+64].T@cores[lo:lo+64],p)
  A=np.remainder(A*pow(3,-1,p),p).astype(np.int64);assert np.array_equal(A,A.T),(p,tag);answers[tag]=A
 return answers
