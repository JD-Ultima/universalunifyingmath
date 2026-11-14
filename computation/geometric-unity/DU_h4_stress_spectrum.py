"""Exact multiplicity Gram matrices for the CJ derivative quadratic on H4 coclosed g1."""
from pathlib import Path
from fractions import Fraction as F
from functools import lru_cache
import contextlib,io,runpy,sys,json,gzip,time
HERE=Path(__file__).resolve().parent
import argparse
parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--source-dir',type=Path);args=parser.parse_args()
old=sys.argv;sys.argv=[str(HERE/'DE_h6_elliptic_periodic.py')]
if args.source_dir:sys.argv+=['--source-dir',str(args.source_dir)]
with contextlib.redirect_stdout(io.StringIO()):g=runpy.run_path(sys.argv[0])
sys.argv=old;BO=g['BO'];s=g['s'];H=g['H'];hnz=g['hnz'];inv=g['inv'];tick=time.time()
form_rows=g['form_rows']

def lower(T):
    out={}
    for key,v in T.items():
        if len(key)==3:
            i,a,m=key
            for j,c in hnz[i]:
                for b,d in hnz[a]:
                    for n,e in form_rows(m):BO.add(out,(j,b,n),v*c*d*e)
        else:
            a,m=key
            for b,d in hnz[a]:
                for n,e in form_rows(m):BO.add(out,(b,n),v*d*e)
    return out
def dot(T,U):return sum((v*U.get(k,0) for k,v in T.items()),F(0))
def psd(A):
    A=A.copy();piv=[]
    for i in range(A.rows):
        d=A[i,i];assert d>=0;piv.append(str(d))
        if d==0:
            assert all(A[i,j]==0 for j in range(i,A.rows));continue
        for j in range(i+1,A.rows):
            for k in range(j,A.rows):A[j,k]=A[k,j]=A[j,k]-A[j,i]*A[i,k]/d
    return piv
pack=json.loads(gzip.decompress((HERE/'DU_h4_first_tableau_witness.json.gz').read_bytes()))
receipt=json.loads((HERE/'DU_h4_stress_spectrum.json').read_text());known_receipts={tuple(r['output']):r for r in receipt['records']};records=[]
for block in pack['blocks']:
    known=block['result'];lam=known['output'];h=known['coclosed_kernel_copies']
    Ts=[];Ds=[];bt=time.time()
    for rec in block['witnesses']:
        T={(i,a,m):F(z) for i,a,m,z in rec['entries']};Ts.append(T);delta={}
        for (i,a,m),v in T.items():
            if not m>>i&1:BO.add(delta,(a,m|(1<<i)),BO.wedge_sign(i,m)*v)
        Ds.append(delta)
    lows=list(map(lower,Ts));dlows=list(map(lower,Ds));swaps=[{(a,i,m):v for (i,a,m),v in T.items()} for T in Ts]
    M=s.Matrix([[dot(T,U) for U in Ts] for T in lows]);X=s.Matrix([[dot(T,U) for U in swaps] for T in lows]);D=s.Matrix([[dot(T,U) for U in Ds] for T in dlows])
    assert M==M.T and X==X.T and D==D.T
    R=s.Matrix.hstack(*[s.Matrix([s.Rational(x) for x in col]) for col in known['relations']])
    gram=R.T*M*R;quadratic=R.T*(2*X-M+D)*R
    gp=psd(gram);assert all(s.Rational(x)>0 for x in gp)
    lp=psd(quadratic+3*gram);up=psd(s.Rational(17,5)*gram-quadratic)
    rec=dict(output=lam,multiplicity=h,positive_Gram_pivots=gp,lower_PSD_pivots=lp,upper_PSD_pivots=up,lower_equality_copies=lp.count('0'),upper_equality_copies=up.count('0'))
    assert rec==known_receipts[tuple(lam)]
    records.append(rec);print(lam,h,rec['lower_equality_copies'],rec['upper_equality_copies'],flush=True)
dim=lambda w:(2*w[0]+3)*(2*w[1]+1)*(w[0]-w[1]+1)*(w[0]+w[1]+2)//6
assert len(records)==36
assert sum(dim(r['output'])*r['multiplicity'] for r in records)==140622
assert sum(dim(r['output'])*r['lower_equality_copies'] for r in records)==42537
assert sum(dim(r['output'])*r['upper_equality_copies'] for r in records)==19019
print('Exact endpoints -3 and17/5; equality dimensions42537 and19019; seconds',time.time()-tick)
