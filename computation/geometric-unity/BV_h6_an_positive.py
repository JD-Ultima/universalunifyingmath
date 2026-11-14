"""Exact positive receiver on the entire 952D real AN-closed H6 family.

Default reconstructs and fully verifies character0 (43 dimensions), while
checking the stored rational positivity matrices for all16 characters.
--full reconstructs and exactly evaluates every character: about15 minutes.
Only public BR/BC dependencies, NumPy, and the bundled rational data are used.
No floating rank, eigenvalue, or optimizer result enters the proof.
"""
import os,sys
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1'
sys.dont_write_bytecode=True
from pathlib import Path
from itertools import combinations
import numpy as np
from types import SimpleNamespace
def original_blade(mask):
 x=z=0;phase=1+0j
 for j in range(14):
  if mask>>j&1:
   k=j//2;a=1<<k;b=(1<<k)-1
   if j%2:b|=a;c=1j
   else:c=1
   phase*=c*(-1 if (z&a).bit_count()%2 else 1);x^=a;z^=b
 return x,z,phase
def original_volume(mask):
 x,z,p=original_blade(mask);a,b,c=original_blade((1<<14)-1)
 return x^a,z^b,p*c*(-1 if (z&a).bit_count()%2 else 1)
T=SimpleNamespace(blade=original_blade,blade_volume=original_volume)
FULL=(1<<14)-1;N=32;ROWS=np.arange(N)
CHARS=[(1<<i)|(1<<4) for i in range(4)]+[(1<<i)|(1<<j) for i,j in combinations(range(4),2)]+[0]*4
def symp(u,v):return (((u&127)&(v>>7)).bit_count()+((u>>7)&(v&127)).bit_count())&1
def hermitian(v):
 x,z=v&127,v>>7;return x,z,(1j)**((x&z).bit_count())
def multiply(a,b):
 x,z,p=a;y,w,q=b;return x^y,z^w,p*q*(-1 if (z&y).bit_count()%2 else 1)
def project(v,pairs):
 for e,f in pairs:
  if symp(v,f):v^=e
  if symp(v,e):v^=f
 return v
def complete(candidates,pairs):
 while any(candidates):
  e=next(v for v in candidates if v);f=next(v for v in candidates if symp(e,v))
  pairs.append((e,f));candidates=[project(v,[(e,f)]) for v in candidates]
 return pairs
generators=[]
for i in range(4):
 flip=(1<<i)|(1<<4);mask=sum(1<<a for a,ch in enumerate(CHARS) if (ch&flip).bit_count()%2);assert mask.bit_count()==6
 x,z,p=T.blade(mask);generators.append(x|(z<<7))
PAIRS=complete(generators,[]);assert len(PAIRS)==2
PAIRS=complete([project(1<<j,PAIRS) for j in range(14)],PAIRS);assert len(PAIRS)==7
assert all(symp(e,f)==1 for e,f in PAIRS)
assert all(not symp(v,w) for i,pair in enumerate(PAIRS) for v in pair for other in PAIRS[i+1:] for w in other)
def remap(old):
 x,z,phase=old;v=x|(z<<7);nx=nz=0;product=(0,0,1.+0j)
 for j,(e,f) in enumerate(PAIRS):
  if symp(v,f):nx|=1<<j;product=multiply(product,hermitian(e))
  if symp(e,v):nz|=1<<j;product=multiply(product,hermitian(f))
 assert product[:2]==(x,z)
 newphase=phase/product[2];low=(nx&3)|((nz&3)<<2);lowphase=(1j)**(((nx&3)&(nz&3)).bit_count())
 return low,nx>>2,nz>>2,newphase/lowphase
def low_mu(a,b):
 x,z=a&3,a>>2;y,w=b&3,b>>2;u,v=x^y,z^w
 return (1j)**((x&z).bit_count()+(y&w).bit_count()-(u&v).bit_count())*(-1 if (z&y).bit_count()%2 else 1)
def triple(a,b,c):return low_mu(a,b)*low_mu(a^b,c)
GDATA=[remap(T.blade(1<<a)) for a in range(14)]
LABEL={0:0}
for a,ch in enumerate(CHARS):
 for old,label in list(LABEL.items()):
  target=old^ch;value=label^GDATA[a][0]
  if target in LABEL:assert LABEL[target]==value
  LABEL[target]=value
assert len(LABEL)==16
PARITY=np.array([[1-2*((z&c).bit_count()&1) for c in range(N)] for z in range(N)])
def fwht(array):
 h=1
 while h<N:
  view=array.reshape(*array.shape[:-1],N//(2*h),2,h);u=view[...,0,:].copy();v=view[...,1,:].copy();view[...,0,:]=u+v;view[...,1,:]=u-v;h*=2
 return array
GAMMA=np.zeros((14,N,N),dtype=complex)
for a,(label,x,z,p) in enumerate(GDATA):GAMMA[a,ROWS^x,ROWS]=p*PARITY[z]
LROWS=np.array([ROWS^x for label,x,z,p in GDATA]);RCOLS=LROWS.copy()
RPHASE=np.array([p*PARITY[z] for label,x,z,p in GDATA]);LPHASE=np.array([p*PARITY[z,ROWS^x] for label,x,z,p in GDATA])
def leftall(M):return M[LROWS,:]*LPHASE[:,:,None]
def left(M,a):return M[...,LROWS[a],:]*LPHASE[a,:,None]
def right(M,a):return M[...,RCOLS[a]]*RPHASE[a]
def sum_product(A,B):return A.transpose(1,0,2).reshape(N,-1)@B.reshape(-1,N)
def char(mask):
 result=0
 for a in range(14):
  if mask>>a&1:result^=CHARS[a]
 return result
RECEIVERS=[(a,m) for m in range(1<<14) if m.bit_count() in (1,5,9,13) for a in range(14) if CHARS[a]==char(m)]
EXTRACT=[remap(T.blade(m)) for a,m in RECEIVERS]
assert all(label==GDATA[a][0] for (a,m),(label,x,z,p) in zip(RECEIVERS,EXTRACT))
class Character:
 def __init__(self,ch,coords):
  self.ch=ch;self.coords=coords;self.labels=[LABEL[ch]^gd[0] for gd in GDATA]
  self.data=[remap(T.blade_volume(m)) for a,m in coords]
  assert all(label==self.labels[a] for (a,m),(label,x,z,p) in zip(coords,self.data))
  self.phase=np.array([[triple(GDATA[b][0],self.labels[a],self.labels[b]) for b in range(14)] for a in range(14)])
  self.sandwich=np.array([[triple(self.labels[b],GDATA[a][0],self.labels[b]) for b in range(14)] for a in range(14)])
  assert np.max(abs(self.sandwich.imag))==0
 def matrix(self,coeff):
  grid=np.zeros((14,N,N),dtype=complex)
  for value,(a,m),(label,x,z,p) in zip(coeff,self.coords,self.data):grid[a,x,z]+=value*p
  fwht(grid);C=np.empty_like(grid)
  for x in range(N):C[:,ROWS^x,ROWS]=grid[:,x,:]
  return C
 def residual(self,C):
  S=sum_product(C,C);norm=float(np.trace(S).real/N);q=np.empty_like(C)
  for a in range(14):
   U=sum_product(self.phase[a,:,None,None]*leftall(C[a]),C)
   W=sum_product(self.sandwich[a,:,None,None]*right(C,a),C)
   q[a]=(-2/3)*(left(S,a)+right(S,a)-U-U.conj().T-W)+(26/21)*norm*GAMMA[a]
  return q
 def coefficients(self,q):
  grid=np.empty_like(q)
  for x in range(N):grid[:,x,:]=q[:,ROWS^x,ROWS]
  fwht(grid)
  values=np.array([grid[a,x,z]/(N*p) for (a,m),(label,x,z,p) in zip(RECEIVERS,EXTRACT)])
  assert np.max(abs(values.imag))<1e-9
  return values.real
 def value(self,coeff):return self.coefficients(self.residual(self.matrix(coeff)))

C=sys.modules[__name__]
class Gradient:
 def __init__(self,ch,coords,receiver,p=None):
  self.p=p;self.coords=coords;self.ch=ch;self.labels=[C.LABEL[ch]^d[0] for d in C.GDATA]
  self.data=[C.remap(C.T.blade_volume(m)) for a,m in coords]
  if p:
   self.im=next(pow(a,(p-1)//4,p) for a in range(2,p) if pow(a,(p-1)//2,p)==p-1)
   self.inv32=pow(32,-1,p)
  def phase(z):return (int(z.real)+self.im*int(z.imag))%p if p else z
  self.phase=phase;self.dtype=float if p else complex
  self.rph=np.array([[phase(z)*v for v in C.PARITY[b]] for label,x,b,z in C.GDATA],dtype=self.dtype)
  self.lph=np.array([[phase(z)*v for v in C.PARITY[b,C.ROWS^x]] for label,x,b,z in C.GDATA],dtype=self.dtype)
  self.rph=self.mod(self.rph);self.lph=self.mod(self.lph)
  self.forward=np.array([[phase(C.triple(C.GDATA[b][0],self.labels[a],self.labels[b])) for b in range(14)] for a in range(14)],dtype=self.dtype)
  self.reverse=np.array([[phase(C.triple(self.labels[b],self.labels[a],C.GDATA[b][0])) for b in range(14)] for a in range(14)],dtype=self.dtype)
  self.sandwich=np.array([[phase(C.triple(self.labels[b],C.GDATA[a][0],self.labels[b])) for b in range(14)] for a in range(14)],dtype=self.dtype)
  self.R=self.matrix_entries(receiver,False)
  self.H=self.mod(sum(self.right(self.R[a],a)+self.left(self.R[a],a) for a in range(14)))
  self.Rright=[self.right(self.R,c) for c in range(14)]
  self.Rleft=[self.leftall(self.R[c]) for c in range(14)]
 def mod(self,a):return np.remainder(a,self.p) if self.p else a
 def mul(self,a,b):return self.mod(a@b)
 def leftall(self,A):return self.mod(A[C.LROWS,:]*self.lph[:,:,None])
 def left(self,A,a):return self.mod(A[...,C.LROWS[a],:]*self.lph[a,:,None])
 def right(self,A,a):return self.mod(A[...,C.RCOLS[a]]*self.rph[a])
 def matrix_entries(self,entries,volume):
  grid=np.zeros((14,32,32),dtype=self.dtype)
  for a,m,v in entries:
   label,x,z,phase=C.remap(C.T.blade_volume(m) if volume else C.T.blade(m))
   grid[a,x,z]+=self.phase(phase)*(int(v)%self.p if self.p else v)
  grid=self.mod(grid);C.fwht(grid);grid=self.mod(grid);A=np.empty_like(grid)
  for x in range(32):A[:,C.ROWS^x,C.ROWS]=grid[:,x,:]
  return A
 def matrix(self,v):return self.matrix_entries([(a,m,x) for (a,m),x in zip(self.coords,v)],True)
 def weighted(self,weights,products):return self.mod(np.sum(self.mod(weights[:,None,None]*products),axis=0))
 def core(self,A):
  result=np.empty_like(A)
  for c in range(14):
   value=self.mod(self.mul(A[c],self.H)+self.mul(self.H,A[c]))
   temp=self.mul(A,self.R[c]);term=self.mod(sum(self.mod(self.forward[c,b]*self.right(temp[b],b)) for b in range(14)))
   value=self.mod(value-term-self.weighted(self.forward[:,c],self.mul(self.Rright[c],A)))
   value=self.mod(value-self.weighted(self.reverse[c],self.mul(self.Rleft[c],A)))
   value=self.mod(value-self.weighted(self.reverse[:,c],self.mul(self.right(A,c),self.R)))
   value=self.mod(value-self.weighted(self.sandwich[:,c],self.mul(self.leftall(A[c]),self.R)))
   temp=self.mul(self.R,A[c]);term=self.mod(sum(self.mod(self.sandwich[a,c]*self.right(temp[a],a)) for a in range(14)))
   result[c]=self.mod(value-term)
  grid=np.empty_like(result)
  for x in range(32):grid[:,x,:]=result[:,C.ROWS^x,C.ROWS]
  C.fwht(grid);grid=self.mod(grid)
  if self.p:return np.array([int(grid[a,x,z])*self.inv32%self.p*pow(self.phase(phase),-1,self.p)%self.p for (a,m),(label,x,z,phase) in zip(self.coords,self.data)],dtype=np.int64)
  out=np.array([grid[a,x,z]/(32*phase) for (a,m),(label,x,z,phase) in zip(self.coords,self.data)])
  assert np.max(abs(out.imag))<1e-7
  return out.real

def exact_pd(A):
 """Sylvester criterion via exact fraction-free Schur complements."""
 A=np.array(A,dtype=object);previous=1;pivots=[]
 for k in range(len(A)):
  pivot=int(A[k,k]);assert pivot>0,(k,pivot);pivots.append(pivot)
  if k+1<len(A):
   numer=A[k+1:,k+1:]*pivot-A[k+1:,k,None]*A[k,None,k+1:]
   assert not any(int(v)%previous for v in numer.flat)
   A[k+1:,k+1:]=numer//previous
  previous=pivot
 return {'all_leading_minors_positive':True,'dimension':len(A),'leading_minor_bit_lengths':[v.bit_length() for v in pivots],'last_minor_sha256':hashlib.sha256(str(pivots[-1]).encode()).hexdigest()}

# Appended to the public verifier by build_public_an952.py.
import argparse,gzip,hashlib,io,json,math,time
from fractions import Fraction as F
import BR_h6_an_closed as BR
HERE=Path(__file__).resolve().parent
ALLCHARS=[0,3,5,6,9,10,12,15,17,18,20,23,24,27,29,30]
def independent128(receiver):
 """Original128-dimensional ordered products guard the compressed gradient."""
 p=65521;im=pow(17,(p-1)//4,p)
 if im*im%p!=p-1:im=next(pow(a,(p-1)//4,p) for a in range(2,p) if pow(a,(p-1)//2,p)==p-1)
 phase=lambda z:(int(z.real)+im*int(z.imag))%p
 rr=np.arange(128);par=np.array([[1-2*((z&r).bit_count()%2) for r in range(128)] for z in range(128)],dtype=float)
 gd=[original_blade(1<<a) for a in range(14)];lr=np.array([rr^x for x,z,v in gd]);rp=np.array([phase(v)*par[z]%p for x,z,v in gd]);lp=np.array([phase(v)*par[z,rr^x]%p for x,z,v in gd])
 mod=lambda x:np.remainder(x,p)
 def left(A,a):return mod(A[...,lr[a],:]*lp[a,:,None])
 def leftall(A):return mod(A[lr,:]*lp[:,:,None])
 def right(A,a):return mod(A[...,lr[a]]*rp[a])
 def matrices(entries,volume):
  A=np.zeros((14,128,128))
  for a,m,v in entries:
   x,z,ph=original_volume(m) if volume else original_blade(m)
   A[a,rr^x,rr]=mod(A[a,rr^x,rr]+(int(v)%p)*phase(ph)*par[z])
  return A
 R=matrices(receiver['rows'],False);checks=[]
 for ch in (0,15,27):
  coords=[(a,sum(1<<i for i in inds)) for inds in combinations(range(14),6) for a in range(14) if CHARS[a]^char(sum(1<<i for i in inds))==ch]
  v=np.zeros(len(coords),dtype=np.int64)
  for j in range(7):v[(j*397)%len(v)]=j-3
  A=matrices([(a,m,int(c)) for (a,m),c in zip(coords,v) if c],True);S=mod(np.sum(A@A,axis=0));q=np.empty_like(A)
  for a in range(14):
   U=mod(np.sum(leftall(A[a])@A,axis=0));temp=mod(A@A[a]);V=mod(sum(right(temp[b],b) for b in range(14)));W=mod(np.sum(right(A,a)@A,axis=0))
   q[a]=mod((-2*pow(3,-1,p)%p)*mod(left(S,a)+right(S,a)-U-V-W))
  original=int(np.sum(R*q.transpose(0,2,1)))%p*pow(128,-1,p)%p
  engine=Gradient(ch,coords,receiver['rows'],p);core=engine.core(engine.matrix(v));compressed=(-sum(int(x)*int(y) for x,y in zip(v,core))*pow(3,-1,p))%p
  assert original==compressed,(ch,original,compressed);checks.append(ch)
 return checks
def reconstruct(record):
 ch=record['character'];coords,rows=BR.stack(ch);rank,basis=BR.exact_kernel(rows,len(coords));m=len(basis)
 assert rank==record['rank'] and m==record['dimension']
 assert hashlib.sha256(json.dumps(basis,separators=(',',':')).encode()).hexdigest()==record['basis_sha256']
 coeff=record['coefficient_rows'];den=math.lcm(*(r['denominator'] for r in coeff));Cnum=np.array([[v*(den//r['denominator']) for v in r['numerators']] for r in coeff],dtype=object);Bnum=np.array(basis,dtype=object).T@Cnum
 sp=record['sqrt_interval_power'];dp=record['unit_tensor_power'];metric=list(map(F,BR.g));unitden=den<<sp
 D=np.empty(Bnum.shape,dtype=np.int64);maximum=F(0)
 for i,(a,mask) in enumerate(coords):
  w=1/metric[a]
  for j in range(14):
   if mask>>j&1:w/=metric[j]
  lo=math.isqrt((w.numerator<<(2*sp))//w.denominator)
  assert F(lo*lo,1<<(2*sp))<=w<F((lo+1)**2,1<<(2*sp))
  for j,v in enumerate(Bnum[i]):
   v=int(v);low=v*lo;high=v*(lo+1)
   if low>high:low,high=high,low
   q=((low+high)*(1<<dp)+unitden)//(2*unitden);D[i,j]=q
   maximum=max(maximum,abs(F(low,unitden)-F(q,1<<dp)),abs(F(high,unitden)-F(q,1<<dp)))
 stream=io.BytesIO();np.save(stream,D,allow_pickle=False)
 assert hashlib.sha256(stream.getvalue()).hexdigest()==record['unit_tensor_sha256']
 error2=maximum*maximum*D.size;assert error2==F(*record['frobenius_error_squared_bound']) and error2<F(1,10**16)
 gram=D.astype(object).T@D.astype(object);scale=1<<(2*dp)
 lower=min(int(gram[i,i])-sum(abs(int(gram[i,j])) for j in range(m) if i!=j) for i in range(m));upper=max(sum(abs(int(v)) for v in row) for row in gram)
 assert 2*lower>scale and 10000*upper<10001*scale
 return coords,D,gram
def check_contraction(record,cert,receiver):
 coords,D,gram=reconstruct(record);m=D.shape[1];primes=cert['crt_primes'];modulus=math.prod(primes);J=np.array(cert['matrix_integer'],dtype=object)
 assert cert['matrix_denominator']==3*receiver['denominator']*(1<<(2*record['unit_tensor_power']))
 norm2=sum(v*v for a,mask,v in receiver['rows']);bound=450*(math.isqrt(norm2)+1)*max(int(gram[i,i]) for i in range(m))
 assert len(primes)==len(set(primes))==6
 assert J.shape==(m,m)
 assert bound==cert['absolute_integer_entry_bound'] and modulus>2*bound and max(abs(int(v)) for v in J.flat)<=bound
 for p in primes:
  # Trial division is sufficient for the fixed23-bit primes.
  assert p%4==1 and all(p%d for d in range(2,math.isqrt(p)+1))
  assert 64*(p-1)**2<2**53
  engine=Gradient(record['character'],coords,receiver['rows'],p);Dmod=np.remainder(D,p).astype(float)
  cores=np.column_stack([engine.core(engine.matrix(D[:,j])) for j in range(m)]).astype(float);actual=np.zeros((m,m))
  for start in range(0,len(D),64):actual=np.remainder(actual-Dmod[start:start+64].T@cores[start:start+64],p)
  assert np.array_equal(actual.astype(np.int64),np.array(J%p,dtype=np.int64))
  for i,j in [(0,0),(0,m-1),(m-1,m-1)]:assert int(actual[i,j])==(-sum(int(D[r,i])%p*int(cores[r,j]) for r in range(len(D))))%p
 return {'character':record['character'],'dimension':m,'exact_kernel':True,'dyadic_sqrt_enclosure':True,'six_prime_full_contraction':True,'closed_basis_lower_bound':[399,10000]}
def main():
 parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--full',action='store_true');parser.add_argument('--character',type=int,default=0,choices=ALLCHARS);parser.add_argument('--receipt');args=parser.parse_args();start=time.perf_counter()
 data=json.loads(gzip.decompress((HERE/'BV_basis_coefficients.json.gz').read_bytes()));receiver=json.loads((HERE/'BV_receiver.json').read_text());basis={r['character']:r for r in data['basis']};certs={r['character']:r for r in data['certificates']}
 assert sorted(basis)==ALLCHARS and sorted(certs)==ALLCHARS and sum(r['dimension'] for r in basis.values())==952
 ys=receiver['denominator'];assert ys==1<<36
 assert len({(a,m) for a,m,v in receiver['rows']})==len(receiver['rows'])
 assert all(0<=a<14 and 0<m<=FULL and m.bit_count() in (1,5,9,13)
            and CHARS[a]==char(m) for a,m,v in receiver['rows'])
 assert sum(v for a,m,v in receiver['rows'] if m==1<<a)==0
 assert sum(v*v for a,m,v in receiver['rows'])<121*ys*ys
 assert F(150*11,10**8)*(4+F(1,10**8))<F(1,10**4)
 guards=independent128(receiver)
 for ch,cert in certs.items():
  J=np.array(cert['matrix_integer'],dtype=object);assert np.array_equal(J,J.T)
  exact_pd(25*J-cert['matrix_denominator']*np.eye(len(J),dtype=object))
 print('All16 rational Gram matrices exceed1/25 I; exact receiver trace and sign checks passed.',flush=True)
 selected=ALLCHARS if args.full else [args.character];records=[]
 for ch in selected:
  record=check_contraction(basis[ch],certs[ch],receiver);records.append(record);print('Verified',ch,'dimension',record['dimension'],'seconds',round(time.perf_counter()-start,3),flush=True)
 result={'full_theorem_verified':args.full,'verified_dimension':sum(r['dimension'] for r in records),'all_stored_matrices_exactly_positive':True,'independent_original128_guards':guards,'blocks':records,'seconds':time.perf_counter()-start}
 if args.receipt:Path(args.receipt).write_text(json.dumps(result,indent=2))
 print(json.dumps(result,indent=2))
if __name__=='__main__':main()
