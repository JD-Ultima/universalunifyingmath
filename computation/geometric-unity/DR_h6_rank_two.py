"""Exact rank-two stress bound and universal curvature nonattainment."""
from pathlib import Path
from itertools import combinations,product
from fractions import Fraction as F
from math import comb
import argparse,json,sys,time
import sympy as s
HERE=Path(__file__).resolve().parent
parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--source-dir',type=Path);parser.add_argument('--write-receipt',action='store_true');args=parser.parse_args()
SOURCE=args.source_dir or (HERE if (HERE/'BG_h6_curvature.py').exists() else HERE.parents[2]/'universalunifyingmath/computation/geometric-unity')
sys.path.insert(0,str(SOURCE));import BG_h6_curvature as BG;import BC_affine_hooks as BC
tick=time.time();t=s.Symbol('t');records=[]
def sign(i,m):return (-1)**((m&((1<<i)-1)).bit_count())
def matrix(columns):
    keys=sorted(set().union(*map(set,columns)))
    return s.Matrix([[col.get(key,0) for col in columns] for key in keys]) if keys else s.zeros(0,len(columns))
for inside,total in product(range(3),[4,6,8]):
    outside=total-inside;sigma=sum(1<<i for i in range(inside))+sum(1<<i for i in range(2,2+outside))
    columns=[(i,a,sigma^(1<<i)^(1<<a)) for i in range(2) for a in range(14) if (sigma^(1<<i)^(1<<a)).bit_count()==6]
    if not columns:continue
    loc={key:j for j,key in enumerate(columns)};constraints=[];deltas=[];X=s.zeros(len(columns))
    for j,(i,a,m) in enumerate(columns):
        out={};delta={}
        if not m>>a&1:out[0,i,m|(1<<a)]=sign(a,m)
        if m>>a&1:out[1,i,m^(1<<a)]=sign(a,m)
        if not m>>i&1:out[2,a,m|(1<<i)]=sign(i,m)
        if m>>i&1:delta[a,m^(1<<i)]=sign(i,m)
        if (a,i,m) in loc:X[loc[a,i,m],j]=1
        constraints.append(out);deltas.append(delta)
    h=matrix(constraints).nullspace()
    if not h:continue
    R=s.Matrix.hstack(*h);G=R.T*R;D=matrix(deltas);Q=R.T*(2*X-s.eye(len(columns))+D.T*D)*R
    polynomial=s.factor((Q-t*G).det()/G.det());roots=s.polys.polytools.all_roots(s.Poly(polynomial,t));minimum=min(roots)
    assert minimum>=-s.Rational(5,3)
    records.append(dict(inside=inside,outside=outside,character_multiplicity=comb(2,inside)*comb(12,outside),kernel_dimension=len(h),minimum=str(minimum),normalized_characteristic=str(polynomial)))
assert len(records)==6 and sum(r['character_multiplicity']*r['kernel_dimension'] for r in records)==38764
assert [r['character_multiplicity'] for r in records if r['minimum']=='-5/3']==[495]
# Explicit equality formula in an orthonormal frame, omega=e2 e3 e4 e5.
omega=sum(1<<i for i in [2,3,4,5]);J={}
for i in [0,1]:
    for a in range(14):
        if a==i or omega>>a&1:continue
        m=omega|(1<<i)|(1<<a);v=(F(1) if a in [0,1] else F(-1,8))*sign(a,omega)*sign(i,omega|(1<<a));J[i,a,m]=v
norm=sum(z*z for z in J.values());exchange=sum(z*J.get((a,i,m),0) for (i,a,m),z in J.items());delta={};constraints={}
for (i,a,m),z in J.items():
    if not m>>a&1:BC.plus(constraints,(0,i,m|(1<<a)),sign(a,m)*z)
    if m>>a&1:BC.plus(constraints,(1,i,m^(1<<a)),sign(a,m)*z)
    if not m>>i&1:BC.plus(constraints,(2,a,m|(1<<i)),sign(i,m)*z)
    if m>>i&1:BC.plus(delta,(a,m^(1<<i)),sign(i,m)*z)
assert not constraints and norm==F(9,4) and exchange==-2
assert sum(z*z for z in delta.values())==F(5,2)
assert 2*exchange-norm+sum(z*z for z in delta.values())==F(-15,4)
# The compact moment forms are assembled directly from the exact5 by5 model.
prime=65537;P=BG.integer_phi();chars=[(1<<i)|(1<<4) for i in range(4)]+[(1<<i)|(1<<j) for i,j in combinations(range(4),2)]+[0]*4
def char(m):
    c=0
    for i in BC.bit_indices(m):c^=chars[i]
    return c
def mod(v):return v.numerator*pow(v.denominator,-1,prime)%prime
all_forms={}
for a,b in combinations(range(5),2):
    K=[[0]*5 for _ in range(5)];K[a][b]=1;K[b][a]=-1
    all_forms[a,b]=[(i,j,F(BG.trace_product(BG.bracket(K,P[i]),P[j]),25)) for i,j in combinations(range(14),2)]
    all_forms[a,b]=[(i,j,z) for i,j,z in all_forms[a,b] if z]
groups={}
for inds in combinations(range(14),6):
    m=sum(1<<i for i in inds);groups.setdefault(char(m),[]).append(m)
annihilators=[]
for kind in ['so3','torus']:
    forms=[v for (a,b),v in all_forms.items() if not (a<3 and b<3 if kind=='so3' else (a,b) in [(0,1),(2,3)])]
    blocks=[]
    for ch,ms in sorted(groups.items()):
        rows={}
        for c,m in enumerate(ms):
            for k,form in enumerate(forms):
                for i,j,z in form:
                    if m>>i&1 or m>>j&1:continue
                    rows.setdefault((k,m|(1<<i)|(1<<j)),{})[c]=mod(z*sign(i,m)*sign(j,m))
        pivots={}
        for row in rows.values():
            row=dict(row)
            while row:
                j=min(row);v=row[j]
                if j not in pivots:
                    vi=pow(v,-1,prime);pivots[j]={k:z*vi%prime for k,z in row.items()};break
                for k,z in pivots[j].items():
                    row[k]=(row.get(k,0)-v*z)%prime
                    if not row[k]:del row[k]
            if len(pivots)==len(ms):break
        assert len(pivots)==len(ms)
        blocks.append(dict(character=ch,columns=len(ms),rank_mod_prime=len(pivots)))
    assert sum(b['columns'] for b in blocks)==3003
    annihilators.append(dict(removed_subalgebra=kind,number_of_twoforms=len(forms),blocks=blocks,exact_annihilator_dimension=0))
result=dict(ambient_blocks=records,ambient_rank_two_tableau_dimension=38764,ambient_sharp_lower_bound='-5/3',ambient_equality_dimension=495,curvature_nonattainment='all real derivative planes',moment_form_prime=prime,moment_form_annihilators=annihilators,rank_two_full_solution_gradient_bound='15579/125000',explicit_uniform_gap_below_5_over_3=None,scope='The universal centralizer classification and curvature reduction are proved in DR_h6_rank_two_proof.txt. Modular full column rank proves the two exact injectivity statements.',seconds=time.time()-tick)
receipt=HERE/'DR_h6_rank_two.json'
if args.write_receipt:receipt.write_text(json.dumps(result,indent=2))
else:
    known=json.loads(receipt.read_text());assert {k:v for k,v in result.items() if k!='seconds'}=={k:v for k,v in known.items() if k!='seconds'}
print(json.dumps({k:v for k,v in result.items() if k not in ['ambient_blocks','moment_form_annihilators']},indent=2))
