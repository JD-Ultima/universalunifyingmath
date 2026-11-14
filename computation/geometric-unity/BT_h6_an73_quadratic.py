"""Exact pointwise quadratic exclusion of the complete 73D AN character15.

Default reconstructs BR's exact closed kernel and checks representative
quadratic rows with independent rational nested Clifford arithmetic.
--full generates every homogeneous necessary quadratic row and proves rank
2701 modulo65521 (roughly15 minutes on one CPU). No private data is used.
All BLAS products are exact integer sums:128*(p-1)^2<2^53. Each operand is
reduced first; an independent int64 guard checks the runtime arithmetic.
"""
import os,sys,json,time,argparse,hashlib
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1'
sys.dont_write_bytecode=True
from pathlib import Path
from fractions import Fraction as F
from itertools import combinations
import numpy as np
import sympy as s
import BC_affine_hooks as BC
import BC_quadratic_xcheck as X
import BR_h6_an_closed as BR
P=65521;N=128;FULL=(1<<14)-1;START=time.perf_counter()
assert s.isprime(P) and P%120==1 and N*(P-1)**2<2**53
IM=int(s.sqrt_mod(-1,P));INV128=pow(128,-1,P);FACTOR=(-2*pow(3,-1,P))%P
ROWS=np.arange(N);PAULI=[]
for j in range(14):
 k=j//2;x=1<<k;z=(1<<k)-1
 if j%2:z|=x;phase=IM
 else:phase=1
 PAULI.append((x,z,phase))
def blade(mask):
 x=z=0;phase=1
 for j,(a,b,c) in enumerate(PAULI):
  if mask>>j&1:phase=phase*c*(-1 if (z&a).bit_count()%2 else 1)%P;x^=a;z^=b
 return x,z,phase
VOLUME=blade(FULL)
def blade_volume(mask):
 x,z,p=blade(mask);a,b,c=VOLUME
 return x^a,z^b,p*c*(-1 if (z&a).bit_count()%2 else 1)%P
PARITY=np.array([[1-2*((z&c).bit_count()&1) for c in range(N)] for z in range(N)])
LROWS=np.array([ROWS^x for x,z,p in PAULI]);RCOLS=LROWS.copy()
def scalar(x):
 x=F(x);return x.numerator*pow(x.denominator,-1,P)%P
G=list(map(F,BC.METRIC));SQ=[int(s.sqrt_mod(scalar(g),P)) for g in G];ISQ=[pow(x,-1,P) for x in SQ]
RPH=np.array([(p*PARITY[z])%P for x,z,p in PAULI],dtype=float)
LPH=np.array([(p*PARITY[z,ROWS^x])%P for x,z,p in PAULI],dtype=float)
def mod(x):return np.remainder(x,P)
def left_all(M):return mod(M[LROWS,:]*LPH[:,:,None])
def left(M,a):return mod(M[...,LROWS[a],:]*LPH[a,:,None])
def right(M,a):return mod(M[...,RCOLS[a]]*RPH[a])
def transform(array):
 array=array.copy();h=1
 while h<N:
  shape=array.reshape(*array.shape[:-1],N//(2*h),2,h);u=shape[...,0,:].copy();v=shape[...,1,:].copy();shape[...,0,:]=u+v;shape[...,1,:]=u-v;h*=2
 return mod(array)
def tensor_matrix(entries):
 grid=np.zeros((14,N,N))
 for a,m,c in entries:
  v=scalar(c)*ISQ[a]%P
  for j in range(14):
   if m>>j&1:v=v*ISQ[j]%P
  x,z,phase=blade_volume(m)
  grid[a,x,z]=(grid[a,x,z]+v*phase)%P
 grid=transform(grid);C=np.empty_like(grid)
 for x in range(N):C[:,ROWS^x,ROWS]=grid[:,x,:]
 return C
def qvalue(C,guards=False):
 S=mod(np.sum(C@C,axis=0));q=np.empty_like(C)
 for a in range(14):
  U=mod(np.sum(left_all(C[a])@C,axis=0))
  temp=mod(C@C[a]);V=mod(sum(right(temp[b],b) for b in range(14)))
  W=mod(np.sum(right(C,a)@C,axis=0))
  q[a]=mod(FACTOR*mod(left(S,a)+right(S,a)-U-V-W))
 norm=int(np.trace(S))%P*INV128%P
 if guards:
  assert np.max(mod(sum(left(C[a],a) for a in range(14))))==0
  assert np.max(mod(sum(right(C[a],a) for a in range(14))))==0
  trace=sum(int(np.trace(left(q[a],a))) for a in range(14))*INV128%P
  assert trace==(-52*pow(3,-1,P)*norm)%P,(trace,norm)
 return q,norm
CHARS=[(1<<i)|(1<<4) for i in range(4)]+[(1<<i)|(1<<j) for i,j in combinations(range(4),2)]+[0]*4
def char(m):
 x=0
 for i in range(14):
  if m>>i&1:x^=CHARS[i]
 return x
RECEIVERS=[]
for m in range(1<<14):
 if m.bit_count() not in (1,5,9,13):continue
 for a in range(14):
  if char(m)==CHARS[a]:RECEIVERS.append((a,m))
RECEIVER_PAULI=[blade(m) for a,m in RECEIVERS]
def coefficients(q,norm):
 grid=np.empty_like(q)
 for x in range(N):grid[:,x,:]=q[:,ROWS^x,ROWS]
 grid=transform(grid);out=[]
 for (a,m),(x,z,ph) in zip(RECEIVERS,RECEIVER_PAULI):
  v=int(grid[a,x,z])*INV128%P*pow(ph,-1,P)%P
  if m==(1<<a):v=(v+26*pow(21,-1,P)*norm)%P
  out.append(v)
 return np.array(out,dtype=np.int64)
def validate(entries):
 rng=np.random.default_rng(1);A=rng.integers(0,P,size=(N,N),dtype=np.int64);B=rng.integers(0,P,size=(N,N),dtype=np.int64)
 assert np.array_equal(A@B,(A.astype(float)@B.astype(float)).astype(np.int64))
 q,n=qvalue(tensor_matrix(entries),True);out=coefficients(q,n)
 # Every computed entry is an integer (including all intermediate reduced matrices).
 assert np.array_equal(q,np.round(q))
 alg=X.BladeAlgebra(G);fields=[{} for _ in range(14)];det=np.prod(np.array(G,dtype=object));checked=0
 for a,m,v in entries:
  target,c=alg.multiply(m,FULL);fields[a][target]=F(v)*c
 for degree in (1,5,9,13):
  indices=[i for i,(a,m) in enumerate(RECEIVERS) if m.bit_count()==degree]
  active=[i for i in indices if out[i]];zeros=[i for i in indices if not out[i]]
  for index in active[:2]+zeros[:1]:
   a,m=RECEIVERS[index];exact=det*alg.quadratic_coefficient(fields,a,m);expected=scalar(exact)*ISQ[a]%P
   for j in range(14):
    if m>>j&1:expected=expected*ISQ[j]%P
   if m==(1<<a):expected=(expected+26*pow(21,-1,P)*n)%P
   assert int(out[index])==expected;checked+=1
 return checked

def panel_rank(A):
 A=A.astype(np.int64).copy();n,m=A.shape;r=0;rows=np.arange(n);cols=[]
 for j in range(m):
  hit=np.flatnonzero(A[r:,j])
  if not len(hit):continue
  k=r+int(hit[0]);A[[r,k]]=A[[k,r]];rows[[r,k]]=rows[[k,r]]
  A[r,j:]=A[r,j:]*pow(int(A[r,j]),-1,P)%P
  if r+1<n:A[r+1:,j:]=(A[r+1:,j:]-A[r+1:,j,None]*A[r,None,j:])%P
  cols.append(j);r+=1
  if r==n:break
 return rows[:r],cols
def inverse(A):
 n=len(A);R=np.c_[A.astype(np.int64),np.eye(n,dtype=np.int64)];det=1
 for j in range(n):
  k=j+int(np.flatnonzero(R[j:,j])[0])
  if k!=j:R[[j,k]]=R[[k,j]];det=-det
  v=int(R[j,j]);det=det*v%P;R[j]=R[j]*pow(v,-1,P)%P
  inds=np.r_[np.arange(j),np.arange(j+1,n)];R[inds]=(R[inds]-R[inds,j,None]*R[j,None,:])%P
 assert np.array_equal((A.astype(np.int64)@R[:,n:])%P,np.eye(n,dtype=np.int64))
 return R[:,n:].astype(float),det%P
def rank(matrix,blocksize=64,reverse=False):
 A=np.asarray(matrix,dtype=float).copy();rowids=np.arange(A.shape[0]);offset=0;total=0;panels=[];det=1;start=time.perf_counter()
 if reverse:A=A[::-1].copy();rowids=rowids[::-1]
 while A.shape[0] and A.shape[1]:
  width=min(blocksize,A.shape[1]);sel,cols=panel_rank(A[:,:width]);r=len(sel)
  keep=np.ones(A.shape[0],dtype=bool);keep[sel]=False
  if r:
   inv,pdet=inverse(A[np.ix_(sel,cols)]);coeff=np.remainder(A[keep][:,cols]@inv,P)
   A=np.remainder(A[keep,width:]-coeff@A[sel,width:],P)
   panels.append({'rows':rowids[sel].tolist(),'columns':[offset+c for c in cols],'determinant':pdet});det=det*pdet%P
  else:A=A[:,width:]
  rowids=rowids[keep];offset+=width;total+=r
  if len(panels)%5==0:print('rank',total,'consumed',offset,'shape',A.shape,'seconds',round(time.perf_counter()-start,2),flush=True)
 return {'rank':total,'panels':panels,'nonzero_schur_determinant_product':det,'seconds':time.perf_counter()-start}

def main():
 parser=argparse.ArgumentParser(description=__doc__)
 parser.add_argument('--full',action='store_true',help='Generate all2701 columns and prove complete modular rank.')
 parser.add_argument('--receipt',type=Path,help='Optionally write full minor metadata to this JSON file.')
 args=parser.parse_args()
 # Longest Weyl reversal preserves the kernel, transferring15 to30.
 phi=[s.Matrix(q)/5 for q in BR.Q];metric=s.diag(*map(s.Rational,BR.g));P5=s.zeros(5)
 for i in range(5):P5[i,4-i]=1
 U=s.Matrix(14,14,lambda a,b:s.trace(phi[a]*P5*phi[b]*P5.T)/metric[a,a]);connection=list(map(s.Matrix,BR.L))
 assert P5.det()==U.det()==1 and U*U==s.eye(14) and U.T*metric*U==metric
 for i in range(14):
  assert sum((U[j,i]*connection[j] for j in range(14)),s.zeros(14))+U*connection[i]*U==s.zeros(14)
  reverse=sum(((BR.chars[i]>>j)&1)<<(4-j) for j in range(5))
  assert all(not U[j,i] or BR.chars[j]==reverse for j in range(14))
 columns,rows=BR.stack(15);linear_rank,basis=BR.exact_kernel(rows,len(columns))
 assert len(basis)==73 and linear_rank==len(columns)-73
 entries=[[[a,m,v] for (a,m),v in zip(columns,vec) if v] for vec in basis]
 checks=sum(validate(entries[i]) for i in (0,12,40,72))
 matrices=np.array([tensor_matrix(v) for v in entries]);k=len(matrices)
 pairs=[(i,i) for i in range(k)]+list(combinations(range(k),2)) if args.full else [(i,i) for i in range(8)]+[(0,j) for j in range(1,9)]
 array=np.zeros((len(RECEIVERS),len(pairs)),dtype=np.int64);diagonal={}
 print('AN73 exact basis reconstructed; independent all-grade guards',checks,'columns to test',len(pairs),flush=True)
 for col,(i,j) in enumerate(pairs):
  for z in (i,j):
   if z not in diagonal:
    q,n=qvalue(matrices[z]);diagonal[z]=coefficients(q,n)
  if i==j:value=diagonal[i]
  else:
   q,n=qvalue(mod(matrices[i]+matrices[j]));value=(coefficients(q,n)-diagonal[i]-diagonal[j])%P
  array[:,col]=value
  if (col+1)%100==0:print('Quadratic columns',col+1,'/',len(pairs),'seconds',round(time.perf_counter()-START,2),flush=True)
 traces=[i for i,(a,m) in enumerate(RECEIVERS) if m==(1<<a)]
 zero13=[i for i,(a,m) in enumerate(RECEIVERS) if m==(FULL^(1<<a))]
 assert len(traces)==len(zero13)==14
 assert not np.any(np.sum(array[traces],axis=0)%P)
 assert not np.any(array[zero13])
 proof=rank(array);assert proof['rank']==len(pairs)
 if args.full:assert proof['rank']==2701
 result={'full_rank_proof':args.full,'character':15,'input_dimension':73,'linear_closed_kernel_rank':linear_rank,
  'receiver_count':len(RECEIVERS),'columns_checked':len(pairs),'prime':P,'rank':proof['rank'],
  'corrected_grade_one_trace_exactly_zero':True,'fourteen_structural_grade13_rows_exactly_zero':True,
  'independent_exact_rational_guards':checks,'exact_integer_BLAS_runtime_guard':True,
  'exact_connection_reversal_transfers_exclusion_to_character30':True,
  'basis_sha256':hashlib.sha256(json.dumps(basis,separators=(',',':')).encode()).hexdigest(),
  'matrix_sha256':hashlib.sha256(array.tobytes()).hexdigest(),'schur_minor_product':proof['nonzero_schur_determinant_product'],
  'seconds':time.perf_counter()-START}
 if args.receipt:
  detail=dict(result,panels=proof['panels'],receivers=RECEIVERS,monomials=pairs)
  args.receipt.write_text(json.dumps(detail,indent=2))
 print(json.dumps(result,indent=2))
if __name__=='__main__':main()
