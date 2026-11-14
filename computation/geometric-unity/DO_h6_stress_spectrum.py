"""Exact sharp CJ derivative bound on the complete DM first tableau."""
from pathlib import Path
from fractions import Fraction as F
import argparse,contextlib,io,runpy,sys,json,gzip,time
HERE=Path(__file__).resolve().parent
parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--source-dir',type=Path);parser.add_argument('--verify-tableau',action='store_true');parser.add_argument('--write-receipt',action='store_true');args=parser.parse_args()
old=sys.argv
if args.verify_tableau:
    sys.argv=[str(HERE/'DM_h6_first_tableau.py'),'--full']
    if args.source_dir:sys.argv+=['--source-dir',str(args.source_dir)]
    runpy.run_path(sys.argv[0])
sys.argv=[str(HERE/'DE_h6_elliptic_periodic.py')]
if args.source_dir:sys.argv+=['--source-dir',str(args.source_dir)]
with contextlib.redirect_stdout(io.StringIO()):g=runpy.run_path(sys.argv[0])
sys.argv=old;BO=g['BO'];s=g['s'];hnz=g['hnz'];inv=g['inv'];form_rows=g['form_rows'];tick=time.time()
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
def pivots(A):
    A=A.copy();out=[]
    for i in range(A.rows):
        d=A[i,i];assert d>0;out.append(str(d))
        for j in range(i+1,A.rows):
            for k in range(j,A.rows):A[j,k]=A[k,j]=A[j,k]-A[j,i]*A[i,k]/d
    return out
def semidefinite_pivots(A):
    A=A.copy();out=[]
    for i in range(A.rows):
        d=A[i,i];assert d>=0;out.append(str(d))
        if d==0:
            assert all(A[i,j]==0 for j in range(i,A.rows))
            continue
        for j in range(i+1,A.rows):
            for k in range(j,A.rows):A[j,k]=A[k,j]=A[j,k]-A[j,i]*A[i,k]/d
    return out
pack=json.loads(gzip.decompress((HERE/'DM_h6_first_tableau_witness.json.gz').read_bytes()));records=[];t=s.Symbol('t');dimensions=0
for block in pack['blocks']:
    known=block['result'];lam=known['output'];h=known['closure_kernel_copies']
    if not h:continue
    Ts=[];Ds=[]
    for rec in block['witnesses']:
        T={(i,a,m):F(z) for i,a,m,z in rec['entries']};Ts.append(T);delta={}
        for (i,a,m),v in T.items():
            for j in BO.axes(m):
                if inv[i][j]:BO.add(delta,(a,m^(1<<j)),inv[i][j]*BO.wedge_sign(j,m)*v)
        Ds.append(delta)
    lows=list(map(lower,Ts));dlows=list(map(lower,Ds));swaps=[{(a,i,m):v for (i,a,m),v in T.items()} for T in Ts]
    M=s.Matrix([[dot(T,U) for U in Ts] for T in lows]);X=s.Matrix([[dot(T,U) for U in swaps] for T in lows]);D=s.Matrix([[dot(T,U) for U in Ds] for T in dlows])
    assert M==M.T and X==X.T and D==D.T
    R=s.Matrix.hstack(*[s.Matrix([s.Rational(x) for x in col]) for col in known['relations']])
    G=R.T*M*R;Q=R.T*(2*X-M+D)*R
    rec=dict(output=lam,multiplicity=h,gram=[[str(x) for x in row] for row in G.tolist()],quadratic=[[str(x) for x in row] for row in Q.tolist()],positive_Gram_pivots=pivots(G))
    if lam in [[8,1],[8,3]]:
        polynomial=s.factor((Q-t*G).det()/G.det());expected=(t-1)**(2 if lam==[8,1] else 1)*(5*t*t+20*t+19)/5
        if lam==[8,3]:expected=-expected
        assert s.expand(polynomial-expected)==0
        rec['normalized_characteristic']=str(polynomial)
    else:rec['positive_Q_plus_12_over_5_Gram_pivots']=pivots(Q+s.Rational(12,5)*G)
    if lam in [[7,5],[6,5],[6,4]]:rec['nonnegative_Q_plus_2_Gram_pivots']=semidefinite_pivots(Q+2*G)
    records.append(rec)
    l,m=lam;dimensions+=h*(2*l+3)*(2*m+1)*(l-m+1)*(l+m+2)//6
assert len(records)==34 and dimensions==78229
result=dict(records=records,tableau_dimension=dimensions,sharp_lower_bound='-2-sqrt(5)/5',saturation_types=[[8,1],[8,3]],saturation_dimension=2565,full_solution_gradient_bound='5193/(25000*(2+sqrt(5)/5))',pure55_gradient_penalty='sqrt(5)/5',rank_one_quadratic=0,scope='Exact derivative quadratic on DM g1; the second-order consequence uses CJ stress divergence and its curvature-source bound.',seconds=time.time()-tick)
receipt=HERE/'DO_h6_stress_spectrum.json'
if args.write_receipt:receipt.write_text(json.dumps(result,indent=2))
else:
    known=json.loads(receipt.read_text());assert {k:v for k,v in known.items() if k!='seconds'}=={k:v for k,v in result.items() if k!='seconds'}
print(json.dumps({k:v for k,v in result.items() if k!='records'},indent=2))
