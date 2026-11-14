"""Exact first-order equations forced on the pure55 part of mixed fields.

DE supplies independently checked geometry and Hermitian norm routines.
No private factor caches are needed: every highest tensor is included.
"""
from pathlib import Path
from fractions import Fraction as F
from collections import Counter
import argparse,contextlib,io,runpy,sys,json,gzip,time
HERE=Path(__file__).resolve().parent
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--source-dir',type=Path)
parser.add_argument('--write-receipt',action='store_true')
args=parser.parse_args();tick=time.time()
old=sys.argv;sys.argv=[str(HERE/'DE_h6_elliptic_periodic.py')]
if args.source_dir:sys.argv+=['--source-dir',str(args.source_dir)]
with contextlib.redirect_stdout(io.StringIO()):
    geom=runpy.run_path(str(HERE/'DE_h6_elliptic_periodic.py'))
sys.argv=old
BO=geom['BO'];s=BO.s;weights=geom['weights'];inv=geom['inv'];curv=geom['curv']
roots=geom['roots'];norm=geom['norm'];casimir=geom['casimir'];raising=geom['raising']
def add(out,key,v):BO.add(out,key,v)
def dim(w):
    l,m=w;return (2*l+3)*(2*m+1)*(l-m+1)*(l+m+2)//6
def C(w):return w[0]*(w[0]+3)+w[1]*(w[1]+1)
V=Counter({(2,0):1,(-2,0):1,(0,2):1,(0,-2):1,(1,1):1,(1,-1):1,(-1,1):1,(-1,-1):1,(1,0):1,(-1,0):1,(0,1):1,(0,-1):1,(0,0):2})
K={(7,3):1,(7,2):1,(7,1):1,(6,3):1,(6,2):1,(6,1):2,(5,5):1,(5,2):2,(5,1):1,(5,0):1,(4,1):2,(3,2):1,(3,0):1,(2,1):1}
def tensor(w):
    out=Counter()
    for mu,n in V.items():
        a,b=2*(w[0]+mu[0])+3,2*(w[1]+mu[1])+1
        if a==0 or b==0 or abs(a)==abs(b):continue
        sign=1 if a*b>0 else -1;a,b=abs(a),abs(b)
        if a<b:a,b=b,a;sign=-sign
        out[(a-3)//2,(b-1)//2]+=sign*n
    out=Counter({w:n for w,n in out.items() if n})
    assert all(n>0 for n in out.values())
    assert sum(dim(t)*n for t,n in out.items())==14*dim(w)
    return out
decompositions={w:tensor(w) for w in K}
assert sum(dim(w)*n for w,n in K.items())==6820
pack=json.loads(gzip.decompress((HERE/'DK_h6_mixed_gradient_witness.json.gz').read_bytes()))
grouped={};ratios=[]
for rec in pack['records']:
    lam=tuple(rec['output']);fiber=tuple(rec['fiber']);copy=rec['copy']
    T={(i,a,m):F(v) for i,a,m,v in rec['entries']}
    assert T and all(not raising(T,r) for r in roots)
    assert all(tuple(weights[i][j]+weights[a][j]+sum(weights[z][j] for z in BO.axes(m)) for j in range(2))==lam for i,a,m in T)
    for i in range(14):
        vec={(a,m):v for (j,a,m),v in T.items() if i==j}
        if not vec:continue
        constraints={}
        for (a,m),v in vec.items():
            for key,c in BO.constraint_entries(a,m,inv,curv):add(constraints,key,c*v)
        assert not constraints
        assert casimir(vec)=={key:C(fiber)*v for key,v in vec.items()}
    image={};contraction={}
    for (i,a,m),v in T.items():
        if not m>>i&1:add(image,(a,m|(1<<i)),BO.wedge_sign(i,m)*v)
        for j in BO.axes(m):
            if inv[i][j]:add(contraction,(a,m^(1<<j)),inv[i][j]*BO.wedge_sign(j,m)*v)
    grouped.setdefault(lam,[]).append((fiber,copy,T,image))
    if lam==(6,4):ratios.append(dict(fiber=list(fiber),wedge_squared=str(norm(image)/norm(T)),contraction_squared=str(norm(contraction)/norm(T))))
records=[]
for lam,columns in grouped.items():
    expected={(w,c) for w,n in K.items() for c in range(n) if decompositions[w][lam]}
    assert all(decompositions[w][lam]==1 for w,c in expected)
    assert {(w,c) for w,c,T,im in columns}==expected
    labels=[[ *w,c] for w,c,T,im in columns]
    keys=sorted(set().union(*(set(T) for w,c,T,im in columns)))
    source=s.Matrix([[T.get(k,0) for w,c,T,im in columns] for k in keys])
    assert source.rank()==len(columns)
    keys=sorted(set().union(*(set(im) for w,c,T,im in columns)))
    matrix=s.Matrix([[im.get(k,0) for w,c,T,im in columns] for k in keys])
    rank=matrix.rank();lower=matrix[:,1:].rank()
    assert columns[0][0]==(5,5)
    records.append(dict(output=list(lam),labels=labels,rank=rank,lower_only_rank=lower,can_cancel_pure55=(rank==lower)))
assert [r['can_cancel_pure55'] for r in records]==[True,False,False,False]
assert ratios[0]['wedge_squared']=='26/27' and ratios[0]['contraction_squared']=='13/27'
spectrum=[dict(p=p,dimension=(p+1)*(11-p),Q_eigenvalue=35-5*(p-5)**2) for p in range(11)]
assert sum(r['dimension'] for r in spectrum)==286
assert sum(r['dimension']*r['Q_eigenvalue'] for r in spectrum)==0
assert F(14*sum(r['dimension']*r['Q_eigenvalue']**2 for r in spectrum),286*20)==945
result=dict(records=records,allowed_55_derivative_types=[[7,5],[6,5],[6,4]],channel64_norm_ratios=ratios,axial_Q_spectrum=spectrum,TTstar=945,delta_to_d_squared_ratio='1/2',scope='Exact mixed first-order closure projection; no nonlinear existence conclusion.',seconds=time.time()-tick)
receipt=HERE/'DK_h6_mixed_gradient.json'
if args.write_receipt:receipt.write_text(json.dumps(result,indent=2))
else:
    known=json.loads(receipt.read_text());assert {k:v for k,v in known.items() if k!='seconds'}=={k:v for k,v in result.items() if k!='seconds'}
print(json.dumps(result,indent=2))
