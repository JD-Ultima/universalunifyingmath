"""CW: exact integer-limb H4 Clifford contractions (NumPy only)."""
import os,sys
os.environ['OPENBLAS_NUM_THREADS']='1';sys.dont_write_bytecode=True
from itertools import combinations
from types import SimpleNamespace
import numpy as np
N=128;HALF=64;ROWS=np.arange(N)
PAULI=[]
for j in range(14):
 k=j//2;x=1<<k;z=(1<<k)-1
 if j%2:z|=x;phase=1j
 else:phase=1.+0j
 PAULI.append((x,z,phase))
def blade(mask):
 x=z=0;phase=1.+0j
 for j,(a,b,c) in enumerate(PAULI):
  if mask>>j&1:phase*=c*(-1 if (z&a).bit_count()%2 else 1);x^=a;z^=b
 return x,z,phase
PARITY=np.array([[1-2*((z&c).bit_count()&1) for c in range(N)] for z in range(N)],dtype=np.int8)
GAMMA=np.zeros((14,N,N),dtype=np.complex128)
for a,(x,z,p) in enumerate(PAULI):GAMMA[a,ROWS^x,ROWS]=p*PARITY[z]
def fwht(array):
 h=1
 while h<N:
  view=array.reshape(*array.shape[:-1],N//(2*h),2,h);u=view[...,0,:].copy();v=view[...,1,:].copy();view[...,0,:]=u+v;view[...,1,:]=u-v;h*=2
 return array
previous=SimpleNamespace(blade=blade,ROWS=ROWS,PARITY=PARITY,GAMMA=GAMMA,fwht=fwht)
class Clifford:
    """Exploit the parity block decomposition of the same Jordan–Wigner matrices."""

    def __init__(self, masks, *, dual=False, target=6/7):
        self.order = np.array(sorted(range(N), key=lambda i: i.bit_count() % 2))
        self.inverse = np.argsort(self.order)
        self.rows = np.arange(HALF)
        self.gamma = self.split(previous.GAMMA, odd=True)
        row, col = np.argmax(abs(self.gamma), axis=-1), np.argmax(abs(self.gamma), axis=-2)
        phase = np.take_along_axis(self.gamma, row[..., None], axis=-1)[..., 0]
        right_phase = np.take_along_axis(self.gamma, col[..., None, :], axis=-2)[..., 0, :]
        self.left_rows = row
        self.left_phase = phase[..., None]
        self.right_columns = col
        self.right_phase = right_phase[..., None, :]
        self.target = target
        pauli_blade = previous.blade_volume if dual else previous.blade
        self.x, self.z, self.phase = map(np.asarray, zip(*map(pauli_blade, map(int, masks))))
        self.xor = previous.ROWS[:, None] ^ previous.ROWS[None, :]
        self.perm_xor = self.xor[np.ix_(self.order, self.order)]
        self.even_xor = np.stack([self.perm_xor[:HALF, :HALF], self.perm_xor[HALF:, HALF:]])
        self.even_xor = self.even_xor[None, ...]

    def split(self, full, odd=False):
        matrix = full[..., self.order[:, None], self.order[None, :]]
        a = HALF if odd else 0
        return np.stack((matrix[..., :HALF, a:a + HALF], matrix[..., HALF:, HALF - a:N - a]), axis=-3)

    def merge(self, piece, odd=False):
        matrix = np.zeros((*piece.shape[:-3], N, N), dtype=piece.dtype)
        a = HALF if odd else 0
        matrix[..., :HALF, a:a + HALF] = piece[..., 0, :, :]
        matrix[..., HALF:, HALF - a:N - a] = piece[..., 1, :, :]
        return matrix[..., self.inverse[:, None], self.inverse[None, :]]

    def from_tensor(self, b):
        grid = np.zeros((14, N, N), dtype=complex)
        grid[:, self.x, self.z] = b * self.phase
        previous.fwht(grid)
        return np.stack([np.take_along_axis(grid[..., self.order[k:k + HALF]], self.even_xor[:, parity], axis=-2)
                         for parity, k in enumerate((0, HALF))], axis=1)

    def left(self, a, m):
        inds = self.left_rows[a, ..., None]
        mat = m[..., ::-1, :, :]
        ndim = max(inds.ndim, mat.ndim)
        return np.take_along_axis(mat.reshape((1,) * (ndim - mat.ndim) + mat.shape),
                                  inds.reshape((1,) * (ndim - inds.ndim) + inds.shape), axis=-2) * self.left_phase[a]

    def right(self, m, a, *, odd=False):
        k = slice(None, None, -1) if odd else slice(None)
        inds = self.right_columns[a, k, None, :]
        ndim = max(inds.ndim, m.ndim)
        return np.take_along_axis(m.reshape((1,) * (ndim - m.ndim) + m.shape),
                                  inds.reshape((1,) * (ndim - inds.ndim) + inds.shape), axis=-1) * self.right_phase[a, k]

    @staticmethod
    def anti(m, *, odd=False):
        k = slice(None, None, -1) if odd else slice(None)
        return m + m[k, :, :].swapaxes(-1, -2).conj()

    @staticmethod
    def sum_product(l, r):
        # Combine all 14 contracted terms into two larger matrix multiplies.
        return l.transpose(1, 2, 0, 3).reshape(2, HALF, -1) @ r.transpose(1, 0, 2, 3).reshape(2, -1, HALF)


class Engine:
 def __init__(self,receivers):
  masks=[sum(1<<i for i in p) for p in combinations(range(14),4)];self.obj=SimpleNamespace(frame=SimpleNamespace(masks=masks,shape=(14,1001),location={m:i for i,m in enumerate(masks)}));self.cc=Clifford(masks);self.receivers=[tuple(x) for x in receivers];self.bits=21;self.base=1<<21;self.groups=[[] for _ in range(14)]
  assert 2*896*(self.base-1)**2<2**53
  self.lre=np.rint(self.cc.left_phase.real).astype(np.int64);self.lim=np.rint(self.cc.left_phase.imag).astype(np.int64);self.rre=np.rint(self.cc.right_phase.real).astype(np.int64);self.rim=np.rint(self.cc.right_phase.imag).astype(np.int64)
  for k,(a,m) in enumerate(self.receivers):
   x,z,p=previous.blade(m);row=self.cc.inverse[previous.ROWS^x];col=self.cc.inverse[previous.ROWS];self.groups[a].append((k,row//64,row%64,col%64,previous.PARITY[z].astype(np.int64),int(p.real),int(p.imag)))
 def integer(self,A):
  assert np.all(A.real==np.rint(A.real)) and np.all(A.imag==np.rint(A.imag))
  return A.real.astype(np.int64),A.imag.astype(np.int64)
 def sumprod(self,A,B):
  q=self.cc.sum_product(A,B)
  return self.integer(q)
 def left(self,Z,a):
  inds=self.cc.left_rows[a,...,None];out=[]
  for X in Z:
   mat=X[...,::-1,:,:];ndim=max(inds.ndim,mat.ndim);out.append(np.take_along_axis(mat.reshape((1,)*(ndim-mat.ndim)+mat.shape),inds.reshape((1,)*(ndim-inds.ndim)+inds.shape),axis=-2))
  re,im=out;return re*self.lre[a]-im*self.lim[a],im*self.lre[a]+re*self.lim[a]
 def right(self,Z,a):
  inds=self.cc.right_columns[a][...,None,:];out=[]
  for X in Z:
   ndim=max(inds.ndim,X.ndim);out.append(np.take_along_axis(X.reshape((1,)*(ndim-X.ndim)+X.shape),inds.reshape((1,)*(ndim-inds.ndim)+inds.shape),axis=-1))
  re,im=out;return re*self.rre[a]-im*self.rim[a],im*self.rre[a]+re*self.rim[a]
 def core(self,X,Y):
  S=self.sumprod(X,Y);out=[None]*len(self.receivers);flip=Y[:,::-1,:,:]
  # An explicit int64 dot guards the actual BLAS integer product.
  left=X[:,0,0,:].ravel();right=Y[:,0,:,0].ravel()
  er=int(left.real.astype(np.int64)@right.real.astype(np.int64)-left.imag.astype(np.int64)@right.imag.astype(np.int64));ei=int(left.real.astype(np.int64)@right.imag.astype(np.int64)+left.imag.astype(np.int64)@right.real.astype(np.int64));assert (int(S[0][0,0,0]),int(S[1][0,0,0]))==(er,ei)
  for a in range(14):
   if not self.groups[a]:continue
   U=self.sumprod(self.cc.left(slice(None),X[a]),flip)
   temp=self.integer(X@Y[a]);right=self.right(temp,slice(None));V=tuple(np.sum(v,axis=0,dtype=np.int64) for v in right)
   W=self.sumprod(self.cc.right(X,a),flip);ls=self.left(S,a);rs=self.right(S,a)
   Q=tuple(ls[i]+rs[i]-U[i]-V[i]-W[i] for i in (0,1))
   assert max(int(np.max(abs(v))) for v in Q)<2**61
   for k,parity,row,col,sign,pr,pi in self.groups[a]:
    re=sum(int(v)*int(e) for v,e in zip(Q[0][parity,row,col],sign));im=sum(int(v)*int(e) for v,e in zip(Q[1][parity,row,col],sign));out[k]=(re*pr+im*pi,im*pr-re*pi)
  return out
 def value(self,D):
  D=np.asarray(D,dtype=np.int64).reshape(self.obj.frame.shape);assert max(sum(abs(int(v)) for v in row) for row in D)<2**52
  C=self.cc.from_tensor(D.astype(float));re,im=self.integer(C);assert max(int(np.max(abs(re))),int(np.max(abs(im))))<self.base**2
  lo=[v%self.base for v in (re,im)];hi=[(v-w)//self.base for v,w in zip((re,im),lo)]
  assert all(np.max(abs(v))<self.base for v in lo+hi)
  L=lo[0].astype(float)+1j*lo[1].astype(float);H=hi[0].astype(float)+1j*hi[1].astype(float)
  q00=self.core(L,L);q01=self.core(L,H);q10=self.core(H,L);q11=self.core(H,H)
  norm=sum(int(v)*int(v) for v in D.ravel());result=[]
  for (a,m),v,w,x,y in zip(self.receivers,q00,q01,q10,q11):
   real=v[0]+self.base*(w[0]+x[0])+self.base**2*y[0];imag=v[1]+self.base*(w[1]+x[1])+self.base**2*y[1]
   assert imag==0 and real%128==0
   q3=-2*(real//128);result.append(7*q3+(18*norm if m==1<<a else 0))
  return result
