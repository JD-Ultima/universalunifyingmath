"""Independent exact H6 first tableau: selected default or --full verification."""
from pathlib import Path
from fractions import Fraction as F
import contextlib,io,runpy,sys,json,gzip,time,argparse
HERE=Path(__file__).resolve().parent
parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--block',help='One output weight, e.g. 4,1');parser.add_argument('--full',action='store_true');parser.add_argument('--source-dir',type=Path);args=parser.parse_args()
old=sys.argv;sys.argv=[str(HERE/'DK_h6_mixed_gradient.py')]
if args.source_dir:sys.argv+=['--source-dir',str(args.source_dir)]
with contextlib.redirect_stdout(io.StringIO()):g=runpy.run_path(sys.argv[0])
sys.argv=old;BO=g['BO'];s=g['s'];weights=g['weights'];inv=g['inv'];curv=g['curv'];roots=g['roots'];C=g['C'];dim=g['dim'];casimir=g['casimir'];raising=g['raising'];K=g['K'];tensor=g['tensor'];decompositions=g['decompositions'];tick=time.time()
results=[];receipt=json.loads((HERE/'DM_h6_first_tableau.json').read_text());expected_receipts={tuple(r['output']):r for r in receipt['blocks']}
pack=json.loads(gzip.decompress((HERE/'DM_h6_first_tableau_witness.json.gz').read_bytes()))
selected={tuple(map(int,args.block.split(',')))} if args.block else {(1,0),(4,1),(5,3),(6,4),(9,3)}
assert len(pack['blocks'])==38 and len(expected_receipts)==38
for block in pack['blocks']:
    known=block['result'];lam=tuple(known['output'])
    if not args.full and lam not in selected:continue
    sources=[];images=[];contractions=[];labels=[];bt=time.time()
    for rec in block['witnesses']:
        fiber=tuple(rec['fiber']);T={(i,a,m):F(z) for i,a,m,z in rec['entries']};assert tuple(rec['output'])==lam
        assert T and all(not raising(T,r) for r in roots)
        assert all(tuple(weights[i][j]+weights[a][j]+sum(weights[z][j] for z in BO.axes(m)) for j in range(2))==lam for i,a,m in T)
        for i in range(14):
            vec={(a,m):z for (j,a,m),z in T.items() if i==j}
            if not vec:continue
            constraints={}
            for (a,m),z in vec.items():
                for key,c in BO.constraint_entries(a,m,inv,curv):BO.add(constraints,key,z*c)
            assert not constraints,(lam,fiber,i,'constraints')
            assert casimir(vec)=={key:C(fiber)*z for key,z in vec.items()},(lam,fiber,i,'Casimir')
        image={};contraction={}
        for (i,a,m),z in T.items():
            if not m>>i&1:BO.add(image,(a,m|(1<<i)),BO.wedge_sign(i,m)*z)
            for j in BO.axes(m):
                if inv[i][j]:BO.add(contraction,(a,m^(1<<j)),inv[i][j]*BO.wedge_sign(j,m)*z)
        labels.append([*fiber,rec['copy'],rec['highest']]);sources.append(T);images.append(image);contractions.append(contraction)
    expected={(w,c,h) for w,n in K.items() for c in range(n) for h in range(decompositions[w][lam])}
    assert {(tuple(x[:2]),x[2],x[3]) for x in labels}==expected
    assert labels==known['labels']
    def matrix(columns):
        keys=sorted(set().union(*map(set,columns)))
        return s.Matrix([[col.get(k,0) for col in columns] for k in keys]) if keys else s.zeros(0,len(columns))
    assert matrix(sources).rank()==len(sources)
    dm=matrix(images);delta=matrix(contractions);rank=dm.rank()
    assert rank==known['closure_rank']
    relations=s.Matrix.hstack(*[s.Matrix([s.Rational(x) for x in r]) for r in known['relations']]) if known['relations'] else s.zeros(len(sources),0)
    assert relations.cols==len(sources)-rank and relations.rank()==relations.cols
    assert dm*relations==s.zeros(dm.rows,relations.cols)
    delta_rank=(delta*relations).rank();assert delta_rank==known['delta_rank_on_closure_kernel']
    rec=dict(output=list(lam),source_copies=len(sources),closure_rank=rank,closure_kernel_copies=relations.cols,delta_image_copies=delta_rank,raw_source_entries=sum(map(len,sources)),seconds=time.time()-bt)
    assert {k:v for k,v in rec.items() if k!='seconds'}==expected_receipts[lam]
    results.append(rec);print(rec,flush=True)
result=dict(blocks=results,seconds=time.time()-tick)
if args.full:
    assert len(results)==38
    result.update(source_dimension=sum(dim(r['output'])*r['source_copies'] for r in results),closure_rank=sum(dim(r['output'])*r['closure_rank'] for r in results),closure_kernel_dimension=sum(dim(r['output'])*r['closure_kernel_copies'] for r in results),delta_image_dimension_on_closure_kernel=sum(dim(r['output'])*r['delta_image_copies'] for r in results))
    assert result['source_dimension']==95480
    assert all(result[k]==receipt[k] for k in ['source_dimension','closure_rank','closure_kernel_dimension','delta_image_dimension_on_closure_kernel'])
else:
    assert len(results)==len(selected)
    result['scope']='Selected exact blocks passed; use --full to verify all38 blocks and total ranks.'
print({k:v for k,v in result.items() if k!='blocks'},flush=True)
