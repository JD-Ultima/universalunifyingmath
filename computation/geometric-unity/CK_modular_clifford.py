"""Independent128D finite-field Clifford check for CK's rational coefficients."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1'
from fractions import Fraction as F
import numpy as np
import sympy as s
import BC_affine_hooks as BC
P=65521;N=128;ROWS=np.arange(N);INV128=pow(N,-1,P);FACTOR=-2*pow(3,-1,P)%P
assert s.isprime(P) and P%120==1 and 14*N*(P-1)**2<2**53
IM=int(s.sqrt_mod(-1,P))
def scalar(value):
 value=F(value);return value.numerator*pow(value.denominator,-1,P)%P
SQ=[int(s.sqrt_mod(scalar(g),P)) for g in BC.METRIC];ISQ=[pow(v,-1,P) for v in SQ]
GENERATORS=[(1<<(j//2),((1<<(j//2))-1)|((1<<(j//2)) if j%2 else 0),IM if j%2 else 1) for j in range(14)]
def blade(mask):
 x=z=0;phase=1
 for j,(a,b,c) in enumerate(GENERATORS):
  if mask>>j&1:phase=phase*c*(-1 if (z&a).bit_count()%2 else 1)%P;x^=a;z^=b
 return x,z,phase
PARITY=np.array([[1-2*((z&j).bit_count()&1) for j in range(N)] for z in range(N)])
LROWS=np.array([ROWS^x for x,z,c in GENERATORS])
RPH=np.remainder(np.array([c*PARITY[z] for x,z,c in GENERATORS],dtype=float),P)
LPH=np.remainder(np.array([c*PARITY[z,ROWS^x] for x,z,c in GENERATORS],dtype=float),P)
def mod(A):return np.remainder(A,P)
def leftall(A):return mod(A[LROWS,:]*LPH[:,:,None])
def left(A,a):return mod(A[...,LROWS[a],:]*LPH[a,:,None])
def right(A,a):return mod(A[...,LROWS[a]]*RPH[a])
def fwht(A):
 A=A.copy();h=1
 while h<N:
  part=A.reshape(*A.shape[:-1],N//(2*h),2,h);u=part[...,0,:].copy();v=part[...,1,:].copy();part[...,0,:]=u+v;part[...,1,:]=u-v;h*=2
 return mod(A)
def matrix(t):
 grid=np.zeros((14,N,N))
 for (a,m),c in t.items():
  value=scalar(c)*ISQ[a]%P
  for j in range(14):
   if m>>j&1:value=value*ISQ[j]%P
  x,z,phase=blade(m);grid[a,x,z]=(grid[a,x,z]+value*phase)%P
 grid=fwht(grid);A=np.empty_like(grid)
 for x in range(N):A[:,ROWS^x,ROWS]=grid[:,x,:]
 return A
def coefficients(t,receivers):
 A=matrix(t);S=mod(np.sum(A@A,axis=0));Q=np.empty_like(A)
 assert np.array_equal(mod(A[0]@A[0]),np.remainder(A[0].astype(np.int64)@A[0].astype(np.int64),P))
 assert np.max(mod(sum(left(A[a],a) for a in range(14))))==0
 for a in range(14):
  U=mod(np.sum(leftall(A[a])@A,axis=0));temp=mod(A@A[a]);V=mod(sum(right(temp[b],b) for b in range(14)));W=mod(np.sum(right(A,a)@A,axis=0))
  Q[a]=mod(FACTOR*mod(left(S,a)+right(S,a)-U-V-W))
 norm=int(np.trace(S))*INV128%P
 assert sum(int(np.trace(left(Q[a],a))) for a in range(14))*INV128%P==-12*norm%P
 grid=np.empty_like(Q)
 for x in range(N):grid[:,x,:]=Q[:,ROWS^x,ROWS]
 grid=fwht(grid);out={}
 for a,m,weight in receivers:
  x,z,phase=blade(m);value=int(grid[a,x,z])*INV128%P*pow(phase,-1,P)%P*SQ[a]%P
  for j in range(14):
   if m>>j&1:value=value*SQ[j]%P
  out[a,m]=value
 return out
