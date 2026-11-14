"""Exact second-jet obstruction to every H2 component of a closed H6 field.

Reuses the public BO full-kernel witnesses and BR exact sparse elimination.
It reconstructs the complete seven-dimensional equivariant symmetric-jet
space, independently evaluates eight closure/curvature rows, and checks a
nonzero cokernel functional. No numerical rank inference is used.
"""
from pathlib import Path
from fractions import Fraction as F
from functools import lru_cache,reduce
from collections import Counter,deque
from itertools import product
from math import gcd,lcm
from time import perf_counter
import json
import sympy as s

HERE=Path(__file__).resolve().parent
import BO_h6_full_kernel as BO
import BR_h6_an_closed as BR

TYPES=((5,1),(3,2),(2,1),(2,0))
ROWS=((0,0,254),(0,0,1468),(0,0,2001),(0,0,2020),(0,0,3892),(0,0,4572),(0,0,5976),(0,1,253))
EXPECTED=[
 [0,60,-90,-22890,60,942,0],
 [0,20,-23,-3033,20,471,0],
 [-2,-16,-12,-4796,-16,628,-9],
 [0,-18,-182,-1800,-18,360,0],
 [0,-1,7,-782,-1,-314,0],
 [0,-20,30,-11520,-16,360,0],
 [0,0,-14,-7346,5,1438,0],
 [-5,-80,-105,-13700,-80,364,-20],
]
EXPECTED_RHS=[F(-30),F(-10),F(8),F(9),F(1,2),F(8),-F(3,2),F(40)]
DUAL=[915357947,-4246099842,0,-228528579,-717216798,-1258489440,1006791552,0]

def dimension(w):
    l,m=w;return (l-m+1)*(l+m+2)*(2*l+3)*(2*m+1)//6

@lru_cache(None)
def partition(a,b):
    if a<0:return 0
    return sum(b+n1-n4>=0 for n1 in range(a+1) for n4 in range(a-n1+1))

def multiplicity(lam,mu):
    x=(2*lam[0]+3,2*lam[1]+1);y=(2*mu[0]+3,2*mu[1]+1);value=0
    for perm in ((0,1),(1,0)):
        for signs in product((-1,1),repeat=2):
            a=signs[0]*x[perm[0]]-y[0];b=signs[1]*x[perm[1]]-y[1]
            if a%2 or b%2:continue
            value+=(1 if perm==(0,1) else -1)*signs[0]*signs[1]*partition(a//2,b//2)
    return value

def tensor_product_decomposition(weights,kernel):
    character=Counter((weights[i][0]+weights[j][0]+weights[k][0],weights[i][1]+weights[j][1]+weights[k][1]) for i in range(14) for j in range(i,14) for k in range(14))
    assert sum(character.values())==1470
    result=[]
    while any(character.values()):
        lam=max(w for w,n in character.items() if n);assert lam[0]>=lam[1]>=0
        copies=character[lam];total=0
        for w in character:
            mult=multiplicity(lam,w);assert mult>=0
            character[w]-=copies*mult;assert character[w]>=0;total+=mult
        assert total==dimension(lam);result.append((lam,copies))
    assert sum(dimension(w)*n for w,n in result)==1470
    present={tuple(item['weight']):len(item['vectors']) for item in kernel}
    intersection={w:n for w,n in result if w in present}
    assert intersection=={(5,1):1,(3,2):1,(2,1):1,(2,0):4}
    assert all(present[w]==1 for w in intersection)
    print('Complete S²V tensor V character:',result,flush=True)
    print('Exactly7 intertwiners from K types',intersection,flush=True)
    return result

def primitive(tensor):
    denominator=lcm(*(v.denominator for v in tensor.values()));values={key:int(v*denominator) for key,v in tensor.items()};divisor=reduce(gcd,values.values())
    return {key:F(v//divisor) for key,v in values.items()}

def generate_modules(weights,raising,kernel):
    tau=[4,6,5,7,0,2,1,3,10,11,8,9,12,13];lowering=[]
    for action in raising:
        opposite=[[] for _ in range(14)]
        for a,row in enumerate(action):opposite[tau[a]]=[(tau[b],c) for b,c in row]
        lowering.append(opposite)
    for action,alpha in zip(lowering,[(-1,1),(0,-1)]):
        for a,row in enumerate(action):
            for b,c in row:assert tuple(weights[b][k]-weights[a][k] for k in range(2))==alpha
    source={tuple(d['weight']):d['vectors'] for d in kernel};modules={}
    for typ in TYPES:
        assert len(source[typ])==1
        seed={(a,m):F(v) for a,m,v in source[typ][0]};basis={};queue=deque()
        def insert(wt,tensor):
            tensor=dict(tensor)
            for key,previous in basis.get(wt,[]):
                factor=tensor.get(key,0)
                if factor:
                    for label,value in previous.items():BO.add(tensor,label,-factor*value)
            if not tensor:return
            key=min(tensor);value=tensor[key];tensor={label:v/value for label,v in tensor.items()}
            basis.setdefault(wt,[]).append((key,tensor));queue.append((wt,tensor))
        insert(typ,seed)
        while queue:
            wt,tensor=queue.popleft()
            for action,alpha in zip(lowering,[(-1,1),(0,-1)]):
                new=BO.root_action(tensor,action)
                if new:insert(tuple(wt[k]+alpha[k] for k in range(2)),new)
        assert sum(map(len,basis.values()))==dimension(typ)
        modules[typ]={wt:[primitive(t) for key,t in vectors] for wt,vectors in basis.items()}
        print('Generated exact K module',typ,'dimension',dimension(typ),flush=True)
    assert sum(dimension(w) for w in modules)==414
    return modules

def symmetric_highest_space(weights,raising,modules):
    alphas=((1,-1),(0,1));lam=(2,0)
    def inventory(wt):
        columns=[]
        for i in range(14):
            for j in range(i,14):
                mw=tuple(wt[k]-weights[i][k]-weights[j][k] for k in range(2))
                for typ,blocks in modules.items():
                    for z in range(len(blocks.get(mw,[]))):columns.append((i,j,typ,mw,z))
        return columns
    @lru_cache(None)
    def root_coordinates(typ,wt,z,which):
        target=BO.root_action(modules[typ][wt][z],raising[which]);out=[];tw=tuple(wt[k]+alphas[which][k] for k in range(2))
        for n,vec in enumerate(modules[typ].get(tw,[])):
            key=min(vec);factor=target.get(key,F(0))/vec[key]
            if factor:
                out.append((n,factor))
                for key,value in vec.items():BO.add(target,key,-factor*value)
        assert not target
        return out
    columns=inventory(lam);locations=[];offset=0
    for alpha in alphas:
        target=inventory(tuple(lam[k]+alpha[k] for k in range(2)));locations.append({key:offset+i for i,key in enumerate(target)});offset+=len(target)
    rows=[{} for _ in range(offset)]
    for col,(i,j,typ,wt,z) in enumerate(columns):
        for which,action in enumerate(raising):
            tw=tuple(wt[k]+alphas[which][k] for k in range(2))
            for n,factor in root_coordinates(typ,wt,z,which):BO.add(rows[locations[which][i,j,typ,tw,n]],col,factor)
            pairs={}
            for a,b in ([(i,j),(j,i)] if i!=j else [(i,i)]):
                for d,value in action[a]:BO.add(pairs,(d,b),value)
                for d,value in action[b]:BO.add(pairs,(a,d),value)
            for (a,b),value in pairs.items():
                if a<=b:BO.add(rows[locations[which][a,b,typ,wt,z]],col,value)
    integer=[]
    for row in rows:
        if not row:continue
        denominator=lcm(*(v.denominator for v in row.values()));d={j:int(v*denominator) for j,v in row.items()};divisor=reduce(gcd,d.values());integer.append({j:v//divisor for j,v in d.items()})
    rank,basis=BR.exact_kernel(integer,len(columns))
    assert len(columns)==862 and rank==855 and len(basis)==7
    print('Complete equivariant symmetric-jet space:862 columns, raising rank855, dimension7.',flush=True)
    return columns,basis

def closure_minor(columns,intertwiners,modules):
    output=[]
    for outer,a,mask in ROWS:
        row=[]
        for i,j,typ,wt,z in columns:
            value=F(0);tensor=modules[typ][wt][z]
            if outer==i and mask>>j&1:value+=BO.wedge_sign(j,mask)*tensor.get((a,mask^(1<<j)),0)
            if i!=j and outer==j and mask>>i&1:value+=BO.wedge_sign(i,mask)*tensor.get((a,mask^(1<<i)),0)
            row.append(value)
        output.append([sum(c*v for c,v in zip(row,vector)) for vector in intertwiners])
    matrix=s.Matrix(output)
    for col in range(7):
        denominator=lcm(*(F(v).denominator for v in matrix[:,col]));divisor=reduce(gcd,[int(v*denominator) for v in matrix[:,col]])
        matrix[:,col]/=s.Rational(divisor,denominator)
    assert matrix.tolist()==EXPECTED
    return matrix

def curvature_rhs(kernel,curvature):
    entries=next(item for item in kernel if item['weight']==[2,0])['vectors'][0];tensor={(a,m):F(v) for a,m,v in entries}
    @lru_cache(None)
    def action(i,j):
        if i==j:return {}
        sign=1 if i<j else -1;i,j=sorted((i,j))
        cov=tuple(tuple((b,sign*c) for b,c in curvature[a,i,j]) for a in range(14))
        return BO.root_action(tensor,cov)
    rhs=[]
    for i,a,mask in ROWS:
        value=-sum((F(1,2)*BO.wedge_sign(j,mask)*action(i,j).get((a,mask^(1<<j)),0) for j in BO.axes(mask)),F(0))
        rhs.append(value)
    assert rhs==EXPECTED_RHS
    return s.Matrix(rhs)

def main():
    tick=perf_counter();kernel=json.loads((HERE/'BO_h6_full_kernel_witnesses.json').read_text())['extremal_vectors']
    weights,inverse,curvature,raising=BO.build_model();columns,rows,slot=BO.inventory(weights)
    assert BO.kernel_lower_bound(kernel,slot,inverse,curvature,raising)==6834
    tensor_product_decomposition(weights,kernel)
    modules=generate_modules(weights,raising,kernel)
    columns,intertwiners=symmetric_highest_space(weights,raising,modules)
    matrix=closure_minor(columns,intertwiners,modules);rhs=curvature_rhs(kernel,curvature)
    dual=s.Matrix([DUAL]);assert dual*matrix==s.zeros(1,7)
    defect=(dual*rhs)[0];assert defect==1006791552
    determinant=matrix.row_join(rhs).det();assert determinant==241629972480
    print('Exact dual obstruction:',DUAL,flush=True)
    print('Dual curvature defect',defect,'augmented determinant',determinant,flush=True)
    print('No H2 value can occur in any smooth covariantly closed H6 field.',flush=True)
    print('Runtime',round(perf_counter()-tick,3),'seconds.',flush=True)

if __name__=='__main__':main()
