"""Exact three21-bit limb Clifford values with integer Walsh transform."""
import os,sys
os.environ['OPENBLAS_NUM_THREADS']='1';sys.dont_write_bytecode=True
from pathlib import Path
HERE=Path(__file__).resolve().parent
import numpy as np
from CW_integer_clifford import Engine as Base,fwht,N,HALF
class Engine(Base):
 def value(self,D):
  D=np.asarray(D,dtype=np.int64).reshape(self.obj.frame.shape);assert max(sum(abs(int(v)) for v in row) for row in D)<2**62
  cc=self.cc;parts=[]
  for phase in (cc.phase.real,cc.phase.imag):
   grid=np.zeros((14,N,N),dtype=np.int64);grid[:,cc.x,cc.z]=D*phase.astype(np.int64);fwht(grid)
   parts.append(np.stack([np.take_along_axis(grid[...,cc.order[k:k+HALF]],cc.even_xor[:,parity],axis=-2) for parity,k in enumerate((0,HALF))],axis=1))
  assert max(int(np.max(abs(v))) for v in parts)<2**62
  limbs=[]
  for i in range(3):
   if i<2:
    lo=[v%self.base for v in parts];parts=[(v-w)//self.base for v,w in zip(parts,lo)]
   else:lo=parts
   assert all(np.max(abs(v))<self.base for v in lo);limbs.append(lo[0].astype(float)+1j*lo[1].astype(float))
  accumulated=[[0,0] for _ in self.receivers]
  for i,L in enumerate(limbs):
   for j,R in enumerate(limbs):
    value=self.core(L,R);scale=1<<(self.bits*(i+j))
    for out,(re,im) in zip(accumulated,value):out[0]+=re*scale;out[1]+=im*scale
  norm=sum(int(v)**2 for v in D.ravel());result=[]
  for (a,m),(re,im) in zip(self.receivers,accumulated):
   assert im==0 and re%128==0;result.append(-14*(re//128)+(18*norm if m==1<<a else 0))
  return result
