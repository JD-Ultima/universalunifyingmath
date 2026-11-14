"""Independent exact decomposition and norm eigenvalues of the pure55 symbol."""
from pathlib import Path
from fractions import Fraction as F
from itertools import combinations
from functools import lru_cache
import gzip,json,sys,time
HERE=Path(__file__).resolve().parent
import argparse
parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--source-dir',type=Path);parser.add_argument('--write-receipt',action='store_true');args=parser.parse_args()
SOURCE=args.source_dir or (HERE if (HERE/'BO_h6_full_kernel.py').exists() else HERE.parents[2]/'universalunifyingmath/computation/geometric-unity')
sys.path.insert(0,str(SOURCE))
import BO_h6_full_kernel as BO
s=BO.s;weights,inv,curv,roots=BO.build_model();tick=time.time();tau=[4,6,5,7,0,2,1,3,10,11,8,9,12,13]
metric=s.zeros(5);metric[0,1]=metric[1,0]=metric[2,3]=metric[3,2]=2;metric[4,4]=1
p=[]
for i in range(5):
    for j in range(i,5):
        if (i,j) in [(0,1),(2,3),(4,4)]:continue
        e=BO.matrix_units(i,j);p.append((e if i==j else e+e.T)*metric)
for i,j in [(0,1),(2,3)]:p.append((BO.matrix_units(i,j)+BO.matrix_units(j,i)-4*BO.matrix_units(4,4))*metric)
dual=[sum((s.Rational(inv[a][j])*p[j] for j in range(14)),s.zeros(5)) for a in range(14)]
k=[(BO.matrix_units(i,j)-BO.matrix_units(j,i))*metric for i,j in combinations(range(5),2)]
g=s.Matrix([[-s.trace(a*b)/2 for b in k] for a in k]).inv()
actions=[[[ (b,F(s.trace(p[b]*(K*Q-Q*K)))) for b in range(14) if s.trace(p[b]*(K*Q-Q*K))] for Q in dual] for K in k]
def plus(out,key,v):BO.add(out,key,v)
def casimir(vec):
    first=[BO.root_action(vec,a) for a in actions];out={}
    for i in range(10):
        for j in range(10):
            if g[i,j]:
                for key,v in BO.root_action(first[j],actions[i]).items():plus(out,key,-F(g[i,j])*v)
    return out
def raising(T,action):
    out={}
    for (i,a,m),v in T.items():
        for j,c in action[i]:plus(out,(j,a,m),c*v)
        for (b,n),c in BO.root_action({(a,m):v},action).items():plus(out,(i,b,n),c)
    return out
H=[[inv[i][tau[j]] for j in range(14)] for i in range(14)];hnz=[[(j,v) for j,v in enumerate(row) if v] for row in H]
@lru_cache(None)
def form_rows(m):
    candidates=[m]
    if bool(m>>12&1)!=bool(m>>13&1):candidates.append(m^(1<<12)^(1<<13))
    rows=[]
    for n in candidates:
        value=F(s.Matrix([[s.Rational(H[i][j]) for j in BO.axes(n)] for i in BO.axes(m)]).det())
        if value:rows.append((n,value))
    return rows
def norm(T):
    total=F(0)
    for key,v in T.items():
        if len(key)==3:
            i,a,m=key
            for j,c in hnz[i]:
                for b,d in hnz[a]:
                    for n,e in form_rows(m):total+=v*c*d*e*T.get((j,b,n),0)
        else:
            a,m=key
            for b,d in hnz[a]:
                for n,e in form_rows(m):total+=v*d*e*T.get((b,n),0)
    return total
pack=json.loads(gzip.decompress((HERE/'DE_h6_symbol_witness.json.gz').read_bytes()));records=[]
for rec in pack['records']:
    T={(i,a,m):F(v) for i,a,m,v in rec['entries']};lam=tuple(rec['weight']);assert T
    assert all(tuple(weights[i][j]+weights[a][j]+sum(weights[z][j] for z in BO.axes(m)) for j in range(2))==lam for i,a,m in T)
    assert all(not raising(T,r) for r in roots)
    for i in range(14):
        vec={(a,m):v for (j,a,m),v in T.items() if j==i}
        if not vec:continue
        con={}
        for (a,m),v in vec.items():
            for key,c in BO.constraint_entries(a,m,inv,curv):plus(con,key,c*v)
        assert not con
        assert casimir(vec)=={key:70*v for key,v in vec.items()}
    out={}
    contraction={}
    for (i,a,m),v in T.items():
        if not m>>i&1:plus(out,(a,m|(1<<i)),BO.wedge_sign(i,m)*v)
        for j in BO.axes(m):
            if inv[i][j]:plus(contraction,(a,m^(1<<j)),inv[i][j]*BO.wedge_sign(j,m)*v)
    source_norm=norm(T);target_norm=norm(out);assert source_norm>0 and target_norm>=0
    l,m=lam;dim=(2*l+3)*(2*m+1)*(l-m+1)*(l+m+2)//6;diff=l*(l+3)+m*(m+1)-80
    if not out:assert not contraction
    row=dict(weight=list(lam),dimension=dim,cross_casimir=diff,symbol_squared_eigenvalue=str(target_norm/source_norm),source_entries=len(T),image_entries=len(out),contraction_entries=len(contraction));records.append(row);print(row,flush=True)
assert sum(x['dimension'] for x in records)==14*286
assert all((F(x['symbol_squared_eigenvalue'])==0)==(x['cross_casimir']>0) for x in records)
minimum=min(F(x['symbol_squared_eigenvalue']) for x in records if F(x['symbol_squared_eigenvalue'])>0)
result=dict(records=records,domain_dimension=4004,kernel_dimension=2024,complement_dimension=1980,real_ellipticity_constant=str(minimum/8),seconds=time.time()-tick,scope='Exact real ellipticity of pure55 exterior closure; decomposition completeness uses BO unique Casimir70 type and the six verified highest vectors with dimensions summing4004.')
b0={(0,811):F(1),(1,4875):F(-1),(1,8971):F(1),(8,783):F(1)};qR={}
for i,j in combinations(range(14),2):
    rb=BO.root_action(b0,[curv[a,i,j] for a in range(14)])
    for (a,m),v in rb.items():
        for k0 in BO.axes(m):
            reduced=m^(1<<k0);sgn=BO.wedge_sign(k0,m)
            for ii,jj,sign in [(i,j,-1),(j,i,1)]:
                if inv[jj][k0] and not reduced>>ii&1:plus(qR,(a,reduced|(1<<ii)),sign*inv[jj][k0]*sgn*BO.wedge_sign(ii,reduced)*v)
assert qR=={key:-7*v for key,v in b0.items()}
result['automatic_coclosure']=True
result['weitzenbock_curvature_scalar']='-7'
result['rough_laplacian_eigenvalue']='7'
C=286*(1+s.sqrt(s.Rational(3,2)));a0=1/(8*C+1)
assert s.simplify(2*C*a0/(1-a0)-s.Rational(1,4))==0
result['periodic_evaluation_error_bound']='1/4'
receipt=HERE/'DE_h6_elliptic_periodic.json'
if args.write_receipt:receipt.write_text(json.dumps(result,indent=2))
else:
    known=json.loads(receipt.read_text());assert {k:v for k,v in result.items() if k!='seconds'}=={k:v for k,v in known.items() if k!='seconds'}
print(json.dumps(result,indent=2))
