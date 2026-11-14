"""Exact AN closed H6 values and a 6534-dimensional closed-field lift.

Default: verify selected global-field/orbit witnesses and reconstruct one exact
character block. --full-exact reconstructs all952 basis directions. No private
research files or floating nullspaces are used. Reuses public BO/BP conventions.
"""
from pathlib import Path
from fractions import Fraction as F
from itertools import combinations
from functools import reduce
from math import gcd,lcm,prod
from time import perf_counter
import argparse
import hashlib
import heapq
import json
import numpy as np
HERE=Path(__file__).resolve().parent
import BC_affine_hooks as BC
import BG_h6_curvature as BG
import BO_h6_full_kernel as BO
import BP_h6_cartan_grade9 as BP
Q=BG.integer_phi();g=list(BC.METRIC)
def matadd(a,b,scale=1):return tuple(tuple(x+scale*y for x,y in zip(ar,br)) for ar,br in zip(a,b))
def matmul(a,b):return tuple(tuple(sum(a[i][k]*b[k][j] for k in range(len(b))) for j in range(len(b[0]))) for i in range(len(a)))
L=[]
for i in range(14):
    K=tuple(tuple(Q[i][r][s]*(1 if r<s else -1 if r>s else 0) for s in range(5)) for r in range(5))
    columns=[BG.bracket(K,Q[a]) for a in range(14)]
    L.append(tuple(tuple(F(BG.trace_product(columns[a],Q[d]),250*g[d]) for a in range(14)) for d in range(14)))
assert all(g[d]*L[i][d][a]+g[a]*L[i][a][d]==0 for i in range(14) for d in range(14) for a in range(14))
brackets={(i,j):tuple(L[i][d][j]-L[j][d][i] for d in range(14)) for i,j in combinations(range(14),2)}
R={}
for i,j in combinations(range(14),2):
    Rij=matadd(matmul(L[i],L[j]),matmul(L[j],L[i]),-1)
    for d,c in enumerate(brackets[i,j]):
        if c:Rij=matadd(Rij,L[d],-c)
    comm=BG.bracket(Q[i],Q[j])
    expected=tuple(tuple(-F(BG.trace_product(BG.bracket(comm,Q[a]),Q[d]),2500*g[d]) for a in range(14)) for d in range(14))
    assert Rij==expected
    R[i,j]=Rij
cov=[{a:[(b,-Li[a][b]) for b in range(14) if Li[a][b]] for a in range(14)} for Li in L]
den=lcm(*(v.denominator for Li in cov for cs in Li.values() for b,v in cs),*(F(1,x).denominator for x in g))
chars=[(1<<i)|(1<<4) for i in range(4)]+[(1<<i)|(1<<j) for i,j in combinations(range(4),2)]+[0]*4
def char(mask):
    c=0
    for i in range(14):
        if mask>>i&1:c^=chars[i]
    return c
def wedge(i,mask):
    if mask>>i&1:return None
    return mask|(1<<i),1-2*((mask&((1<<i)-1)).bit_count()%2)
def add(out,key,c):
    out[key]=out.get(key,0)+c
    if not out[key]:del out[key]
def closure_column(a,mask):
    result={}
    for i in range(14):
        if not any(cov[i].values()):continue
        wed=wedge(i,mask)
        if wed:
            for b,c in cov[i][a]:add(result,(b,wed[0]),wed[1]*int(den*c))
        for r in range(14):
            if not mask>>r&1:continue
            reduced=mask^(1<<r);sgn=1-2*((mask&((1<<r)-1)).bit_count()%2)
            for t,c in cov[i][r]:
                first=wedge(t,reduced)
                if not first:continue
                second=wedge(i,first[0])
                if second:add(result,(a,second[0]),sgn*first[1]*second[1]*int(den*c))
    return result
def check_sparse(tensor):
    out={}
    for (a,m),v in tensor.items():
        for key,c in closure_column(a,m).items():add(out,key,c*v)
    return out

def stack(character):
    columns=[]
    for inds in combinations(range(14),6):
        mask=sum(1<<i for i in inds)
        for a in range(14):
            if chars[a]^char(mask)==character:columns.append((a,mask))
    labels={};rows=[]
    def put(label,column,value):
        if label not in labels:labels[label]=len(rows);rows.append({})
        add(rows[labels[label]],column,int(value))
    for column,(a,mask) in enumerate(columns):
        wed=wedge(a,mask)
        if wed:put((0,0,wed[0]),column,den*wed[1])
        if mask>>a&1:put((1,0,mask^(1<<a)),column,(-1)**(mask&((1<<a)-1)).bit_count()*int(den/g[a]))
        for (b,target),value in closure_column(a,mask).items():
            assert chars[b]^char(target)==character
            put((2,b,target),column,value)
    return columns,rows

def exact_kernel(input_rows,ncols):
    rows=[];colrows=[set() for _ in range(ncols)];heap=[]
    for i,source in enumerate(input_rows):
        divisor=reduce(gcd,source.values()) if source else 1
        row={j:F(v//divisor) for j,v in source.items()};rows.append(row)
        for j in row:colrows[j].add(i)
        heapq.heappush(heap,(len(row),i))
    active=set(range(len(rows)));pivots=[]
    while heap:
        size,i=heapq.heappop(heap)
        if i not in active or len(rows[i])!=size:continue
        active.remove(i);row=rows[i]
        for j in row:colrows[j].remove(i)
        if not row:continue
        j=min(row,key=lambda j:(len(colrows[j]),row[j].denominator,abs(row[j].numerator)))
        value=row[j];row={k:c/value for k,c in row.items()};rows[i]=row
        affected=list(colrows[j]);pivots.append((j,row))
        for r in affected:
            target=rows[r];factor=target[j]
            for k,c in row.items():
                old=target.get(k,F(0));new=old-factor*c
                if new:
                    target[k]=new
                    if not old:colrows[k].add(r)
                elif old:del target[k];colrows[k].remove(r)
            heapq.heappush(heap,(len(target),r))
    pivotcols={j for j,row in pivots};free=[j for j in range(ncols) if j not in pivotcols]
    values=[{} for _ in range(ncols)]
    for z,j in enumerate(free):values[j][z]=F(1)
    for j,row in reversed(pivots):
        out={}
        for k,c in row.items():
            if k!=j:
                for z,v in values[k].items():out[z]=out.get(z,F(0))-c*v
        values[j]={z:v for z,v in out.items() if v}
    basis=[]
    for z in range(len(free)):
        rational=[d.get(z,F(0)) for d in values]
        denominator=lcm(*(v.denominator for v in rational));vec=[int(v*denominator) for v in rational]
        divisor=reduce(gcd,vec);vec=[v//divisor for v in vec]
        assert all(sum(c*vec[j] for j,c in row.items())==0 for row in input_rows)
        basis.append(vec)
    return len(pivots),basis

def verify_block(record):
    tick=perf_counter();character=record['character'];columns,rows=stack(character)
    assert len(columns)==record['columns']
    matrix=np.zeros((len(rows),len(columns)),dtype=np.int64)
    for i,row in enumerate(rows):
        for j,value in row.items():matrix[i,j]=value%BO.PRIME
    modular=BO.modular_rank(matrix)
    exact,basis=exact_kernel(rows,len(columns))
    assert modular==exact==record['rank'] and len(basis)==record['nullity']
    digest=hashlib.sha256(json.dumps(basis,separators=(',',':')).encode()).hexdigest()
    print('Exact AN character',character,'columns',len(columns),'rank',exact,'nullity',len(basis),'seconds',round(perf_counter()-tick,3),'basis SHA256',digest,flush=True)
    return dict(character=character,rank=exact,nullity=len(basis),basis_sha256=digest)

def fp(value):
    value=F(value);assert value.denominator%BO.PRIME
    return value.numerator*pow(value.denominator,-1,BO.PRIME)%BO.PRIME

def verify_witnesses(data):
    selected={tuple(v['label']):{(a,m):F(c) for a,m,c in v['entries']} for v in data['selected_vectors']}
    for label,tensor in selected.items():
        assert BC.hook_traces(tensor)==({}, {}) and not check_sparse(tensor)
        assert all(chars[a]^char(m)==label[0] for a,m in tensor)
    print('All',len(selected),'selected real global AN field values pass both hooks and full closure.',flush=True)
    source=json.loads((HERE/'BO_h6_full_kernel_witnesses.json').read_text())['extremal_vectors']
    weights,inverse,curvature,roots=BO.build_model();columns,rows,slot=BO.inventory(weights)
    assert BO.kernel_lower_bound(source,slot,inverse,curvature,roots)==6834
    byweight={tuple(d['weight']):d['vectors'] for d in source}
    transform=BP.split_to_physical_covectors();total=0;imaginary=256
    assert imaginary*imaginary%BO.PRIME==BO.PRIME-1
    for record in data['orbit_minors']:
        w=tuple(record['weight']);vectors=byweight[w];labels=list(map(tuple,record['pivot_AN_values']));matrix=[]
        assert len(vectors)==len(labels)==record['certified_multiplicity_lower_bound']
        for vector in vectors:
            real,imag=BP.transport(vector,transform);keys=set(real)|set(imag);row=[]
            for label in labels:
                tensor=selected[label];value=0
                for a,m in keys & tensor.keys():
                    norm=g[a]*prod(g[i] for i in range(14) if m>>i&1)
                    value+=(fp(real.get((a,m),0))-imaginary*fp(imag.get((a,m),0)))*fp(tensor[a,m]/norm)
                row.append(value%BO.PRIME)
            matrix.append(row)
        determinant=matrix[0][0] if len(matrix)==1 else (matrix[0][0]*matrix[1][1]-matrix[0][1]*matrix[1][0])%BO.PRIME
        assert matrix==record['minor_mod_p'] and determinant==record['determinant_mod_p'] and determinant
        l,m=w;dimension=(l-m+1)*(l+m+2)*(2*l+3)*(2*m+1)//6
        assert dimension==record['dimension_per_copy'];total+=len(vectors)*dimension
        print('Complex highest-pairing minor',w,'copies',len(vectors),'det mod65537',determinant,flush=True)
    assert total==data['orbit_dimension_lower_bound']==6534
    return dict(selected_global_closed_fields=len(selected),orbit_dimension_lower_bound=total)

def verify_parabolic_characters():
    data=json.loads((HERE/'BG_h6_curvature_kernel.json').read_text())
    basis=[{(a,m):F(v) for a,m,v in item['entries']} for item in data['kernel_tensors']]
    derivative=[check_sparse(t) for t in basis];twist=[]
    for tensor in basis:
        out={}
        for (a,m),v in tensor.items():
            wed=wedge(13,m)
            if wed:add(out,(a,wed[0]),den*wed[1]*v)
        twist.append(out)
    assert [[d.get((0,8413),0) for d in derivative],[d.get((0,8413),0) for d in twist]]==[[0,0,-768],[0,0,0]]
    assert [[d.get((0,8615),0) for d in derivative],[d.get((0,8615),0) for d in twist]]==[[11520,1152,0],[1920,-960,0]]
    for terms in [[(-6,derivative[0]),(60,derivative[1])],[(5,derivative[0]),(10,derivative[1]),(-6,twist[0]),(60,twist[1])],[(5,twist[0]),(10,twist[1])]]:
        out={}
        for coefficient,tensor in terms:
            for key,v in tensor.items():add(out,key,coefficient*v)
        assert not out
    # Row1 forces K2 amplitude0. Row2 is nonzero for every complex z. The
    # nonvanishing polynomial vector (5z-6,10z+60,0) is the complete kernel.
    print('All affine-character parameters: exact rank2; kernel (5z-6,10z+60,0); no H2 component.',flush=True)
    return dict(rank_for_every_complex_parameter=2,kernel='span((5z-6)K0+(10z+60)K1)',H2_component=False)

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--full-exact',action='store_true',help='Reconstruct exact bases/ranks of all16 blocks.')
    parser.add_argument('--block',type=int,action='append',help='Exact character block(s); default15.')
    parser.add_argument('--witnesses-only',action='store_true',help='Skip block reconstruction; verify closure and orbit witnesses only.')
    parser.add_argument('--parabolic-characters',action='store_true',help='Also verify the complete affine-parabolic character family on the SO4-fixed curvature kernel.')
    args=parser.parse_args();tick=perf_counter()
    data=json.loads(HERE.joinpath('BR_h6_an_closed_witnesses.json').read_text())
    assert sum(d['columns'] for d in data['rank_table'])==42042
    assert sum(d['rank'] for d in data['rank_table'])==41090
    assert sum(d['nullity'] for d in data['rank_table'])==952
    result=verify_witnesses(data)
    if args.parabolic_characters:result['parabolic_characters']=verify_parabolic_characters()
    requested=[d['character'] for d in data['rank_table']] if args.full_exact else args.block or [15]
    result['verified_blocks']=[] if args.witnesses_only else [verify_block(d) for d in data['rank_table'] if d['character'] in requested]
    assert args.witnesses_only or len(result['verified_blocks'])==len(set(requested))
    if args.full_exact:assert sum(d['nullity'] for d in result['verified_blocks'])==952
    result['seconds']=perf_counter()-tick
    print(json.dumps(result,indent=2))

if __name__=='__main__':main()
