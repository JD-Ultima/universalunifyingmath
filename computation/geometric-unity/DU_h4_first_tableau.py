"""Independent exact raw H4 tableau verification, without module caches."""
from pathlib import Path
from fractions import Fraction as F
from functools import lru_cache
from itertools import product
import sys,json,gzip,time,argparse
HERE=Path(__file__).resolve().parent
parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--block');parser.add_argument('--full',action='store_true');parser.add_argument('--source-dir',type=Path);args=parser.parse_args()
SOURCE=args.source_dir or (HERE if (HERE/'BO_h6_full_kernel.py').exists() else HERE.parents[2]/'universalunifyingmath/computation/geometric-unity')
sys.path.insert(0,str(SOURCE))
import BO_h6_full_kernel as BO
from sympy.polys.matrices import DomainMatrix
s=BO.s;weights,inv,curv,roots=BO.build_model();tick=time.time()
receipt=json.loads((HERE/'DU_h4_first_tableau.json').read_text());graph=receipt['character_graph'];expected={tuple(r['output']):r['source_copies'] for r in graph['records']}
# Independently reconstruct tensor multiplicities from BX's complete value character.
from collections import Counter
V=Counter(weights);assert sum(V.values())==14
calc=Counter()
for row in json.loads((SOURCE/'BX_data/h4_character.json').read_text())['components']:
    l,m=row['weight']
    for (u,v),n in V.items():
        a,b=2*(l+u)+3,2*(m+v)+1
        if a==0 or b==0 or abs(a)==abs(b):continue
        sign=1 if a*b>0 else -1;a,b=abs(a),abs(b)
        if a<b:a,b=b,a;sign=-sign
        calc[(a-3)//2,(b-1)//2]+=sign*n*row['copies']
assert {w:n for w,n in calc.items() if n}==expected
pack=json.loads(gzip.decompress((HERE/'DU_h4_first_tableau_witness.json.gz').read_bytes()))
selected={tuple(map(int,args.block.split(',')))} if args.block else {(0,0),(2,1),(4,1),(6,4),(9,0)}
expected_receipts={tuple(r['output']):r for r in receipt['blocks']}
def raising(T,action):
    out={}
    for (i,a,m),v in T.items():
        for j,c in action[i]:BO.add(out,(j,a,m),c*v)
        for (b,n),c in BO.root_action({(a,m):v},action).items():BO.add(out,(i,b,n),c)
    return out
@lru_cache(None)
def star(m):
    inds=list(BO.axes(m));assert len(inds)==4
    rows={tuple(sorted(choice)) for choice in product(*[[j for j in range(14) if inv[i][j]] for i in inds]) if len(set(choice))==4}
    out=[]
    for row in rows:
        d=F(s.Matrix([[inv[i][j] for j in inds] for i in row]).det())
        n=((1<<14)-1)^sum(1<<i for i in row)
        sign=(-1)**sum((n&((1<<i)-1)).bit_count() for i in row)
        if d:out.append((n,d*sign))
    return out
@lru_cache(None)
def constraints(a,m):
    out={}
    for key,v in BO.constraint_entries(a,m,inv,curv):
        if key[0]<2:BO.add(out,key,v)
    for n,d in star(m):
        for key,v in BO.constraint_entries(a,n,inv,curv):
            if key[0]==2:BO.add(out,key,d*v)
    return out
def matrix(columns):
    keys=sorted(set().union(*map(set,columns)))
    return s.Matrix([[col.get(k,0) for col in columns] for k in keys]) if keys else s.zeros(0,len(columns))
def independent(columns):
    # Full column rank at one good prime is a characteristic-zero lower bound.
    p=65537;rows={};pivots={}
    for j,col in enumerate(columns):
        for key,v in col.items():
            v=F(v);z=v.numerator*pow(v.denominator,-1,p)%p
            if z:rows.setdefault(key,{})[j]=z
    for row in rows.values():
        row=dict(row)
        while row:
            j=min(row);v=row[j]
            if j not in pivots:
                vi=pow(v,-1,p);pivots[j]={k:z*vi%p for k,z in row.items()};break
            for k,z in pivots[j].items():
                row[k]=(row.get(k,0)-v*z)%p
                if not row[k]:del row[k]
        if len(pivots)==len(columns):return True
    return False
results=[]
for block in pack['blocks']:
    known=block['result'];lam=tuple(known['output'])
    if not args.full and lam not in selected:continue
    sources=[];images=[];wedges=[];bt=time.time()
    for rec in block['witnesses']:
        T={(i,a,m):F(z) for i,a,m,z in rec['entries']};assert tuple(rec['output'])==lam
        assert T and all(not raising(T,r) for r in roots)
        assert all(tuple(weights[i][j]+weights[a][j]+sum(weights[z][j] for z in BO.axes(m)) for j in range(2))==lam for i,a,m in T)
        cons={}
        for (i,a,m),z in T.items():
            for key,c in constraints(a,m).items():BO.add(cons,(i,*key),z*c)
        assert not cons,(lam,rec['fiber'],'raw hook/dual-curvature')
        contraction={};wedge={}
        for (i,a,m),z in T.items():
            if not m>>i&1:BO.add(wedge,(a,m|(1<<i)),BO.wedge_sign(i,m)*z)
            for j in BO.axes(m):
                if inv[i][j]:BO.add(contraction,(a,m^(1<<j)),inv[i][j]*BO.wedge_sign(j,m)*z)
        sources.append(T);images.append(contraction);wedges.append(wedge)
    assert len(sources)==expected[lam] and independent(sources)
    delta=matrix(images);d=matrix(wedges);rank=DomainMatrix.from_Matrix(delta).rank();assert rank==known['coclosed_rank']
    relations=s.Matrix.hstack(*[s.Matrix([s.Rational(x) for x in r]) for r in known['relations']]) if known['relations'] else s.zeros(len(sources),0)
    assert relations.cols==len(sources)-rank and DomainMatrix.from_Matrix(relations).rank()==relations.cols
    assert delta*relations==s.zeros(delta.rows,relations.cols)
    drank=DomainMatrix.from_Matrix(d*relations).rank();assert drank==known['wedge_rank_on_coclosed_kernel']
    rec=dict(output=list(lam),source_copies=len(sources),coclosed_rank=rank,coclosed_kernel_copies=relations.cols,wedge_image_copies=drank,raw_source_entries=sum(map(len,sources)),seconds=time.time()-bt)
    results.append(rec);print(rec,flush=True)
    assert {k:v for k,v in rec.items() if k!='seconds'}==expected_receipts[lam]
result=dict(blocks=results,seconds=time.time()-tick)
if args.full:
    assert len(results)==36
    dim=lambda w:(2*w[0]+3)*(2*w[1]+1)*(w[0]-w[1]+1)*(w[0]+w[1]+2)//6
    result.update(source_dimension=sum(dim(r['output'])*r['source_copies'] for r in results),coclosed_rank=sum(dim(r['output'])*r['coclosed_rank'] for r in results),coclosed_kernel_dimension=sum(dim(r['output'])*r['coclosed_kernel_copies'] for r in results),wedge_image_dimension_on_coclosed_kernel=sum(dim(r['output'])*r['wedge_image_copies'] for r in results))
    assert result['source_dimension']==145432
    assert all(result[k]==receipt[k] for k in result if k not in ['blocks','seconds'])
else:
    assert len(results)==len(selected)
    result['scope']='Selected blocks only; use --full for all36.'
print({k:v for k,v in result.items() if k!='blocks'},flush=True)
