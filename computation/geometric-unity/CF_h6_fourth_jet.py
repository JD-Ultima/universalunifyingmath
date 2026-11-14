"""Exact characteristic-zero four-jet at a nonzero H6 type(5,5) value.

The only dependency is the public BO_h6_full_kernel geometry builder.
All witness coefficients are checked directly over QQ; CRT provenance has
no role in acceptance. The gzip stores common denominators and integers.

Convention: B(X)=b0+b2(X)+b4(X), ordinary monomial coefficients, in the
radial parallel frame s(X)=exp(X/2) on SL5/SO5. The frame and connection are
E=(ad_X/2)/sinh(ad_X/2), omega=-tanh(ad_X/4)/2. Therefore E is even and
omega odd. For this even jet closure has no degree0 or degree2 term;
the two explicitly verified identities exhaust closure through degree3.
No actual closed field or nonlinear solution follows from a finite jet.
See the CF formalization for formulas and the reality argument.
"""
from pathlib import Path
from fractions import Fraction as F
from collections import Counter,defaultdict
from functools import lru_cache
import sys,json,time
sys.dont_write_bytecode=True
import gzip
HERE=Path(__file__).resolve().parent
PACKED=json.loads(gzip.decompress((HERE/'CF_h6_fourth_jet_witness.json.gz').read_bytes()))
import BO_h6_full_kernel as BO
weights,inverse,curvature,raising=BO.build_model();tick=time.time()
def wedge(i,m):
    if m>>i&1:return None
    return m|(1<<i),BO.wedge_sign(i,m)
def union(a,b):return tuple(sorted(a+b))
def add(d,k,v):BO.add(d,k,v)
def load(degree):
    jet=next(d for d in PACKED['jets'] if d['degree']==degree);den=int(jet['denominator']);out={}
    for entry in jet['entries']:
        key=entry[:-1];mon=tuple(key[:degree]);a,m=key[degree:];v=F(int(entry[-1]),den)
        assert len(BO.axes(m))==6 and 0<=a<14 and 0<m<1<<14
        assert len(mon)==degree and tuple(sorted(mon))==mon and all(0<=i<14 for i in mon)
        assert (a,m) not in out.get(mon,{}) and v
        assert tuple(sum(weights[i][k] for i in mon)+weights[a][k]+sum(weights[i][k] for i in BO.axes(m)) for k in range(2))==(5,5)
        if v:out.setdefault(mon,{})[a,m]=v
    return dict(primes=PACKED['primes_used_only_for_proposals']),out
@lru_cache(None)
def cov(r,i):
    if r==i:return tuple(() for _ in range(14))
    lo,hi=sorted((r,i));sign=1 if r<i else -1
    return tuple(tuple((b,sign*c) for b,c in curvature[a,lo,hi]) for a in range(14))
def R(tensor,r,i):return BO.root_action(tensor,cov(r,i))
def constraints(poly,degree):
    columns=0
    for mon,tensor in poly.items():
        out={}
        for (a,m),v in tensor.items():
            for label,c in BO.constraint_entries(a,m,inverse,curvature):add(out,label,c*v)
        assert not out,(degree,mon,next(iter(out.items())))
        columns+=len(tensor)
    print('Exact hook/curvature membership',degree,'coefficients',columns,'seconds',time.time()-tick,flush=True)
def highest(poly,degree):
    for action in raising:
        out={}
        for mon,tensor in poly.items():
            for key,v in BO.root_action(tensor,action).items():add(out,(mon,*key),v)
            for a,count in Counter(mon).items():
                for b,c in action[a]:
                    new=list(mon);new.remove(a);new.append(b);new=tuple(sorted(new))
                    for key,v in tensor.items():add(out,(new,*key),count*c*v)
        assert not out,(degree,next(iter(out.items())))
    print('Exact two raising equations',degree,'seconds',time.time()-tick,flush=True)
def delta(poly):
    out={}
    for mon,tensor in poly.items():
        for i,count in Counter(mon).items():
            new=list(mon);new.remove(i);new=tuple(new)
            for (a,m),v in tensor.items():
                wed=wedge(i,m)
                if wed:add(out,(new,a,wed[0]),count*wed[1]*v)
    return out

def verify_real_part(b):
    # Conjugation in the compact real form exchanges the standard isotropic
    # axes0<->1 and2<->3; its induced coframe permutation is tau below.
    tau=[4,6,5,7,0,2,1,3,10,11,8,9,12,13];conjugate={}
    for (a,m),v in b.items():
        sequence=[tau[i] for i in BO.axes(m)]
        sign=(-1)**sum(sequence[i]>sequence[j] for i in range(len(sequence)) for j in range(i+1,len(sequence)))
        add(conjugate,(tau[a],sum(1<<i for i in sequence)),sign*v)
    real=dict(b)
    for k,v in conjugate.items():add(real,k,v)
    real={k:v/2 for k,v in real.items()};assert real
    norm=F(0)
    for (a,m),v in real.items():
        for (c,n),w in real.items():
            if inverse[a][c]:
                gram=BO.s.Matrix([[BO.s.Rational(inverse[i][j]) for j in BO.axes(n)] for i in BO.axes(m)])
                norm+=v*w*inverse[a][c]*F(gram.det())
    assert norm>0,norm
    return norm

def verify_no_H2():
    # Weyl denominator extraction of [type(5,5):Sym^d(V) tensor V].
    # It is zero for d=0,2,4, so the highest jets have no H2-valued part.
    from itertools import product
    chars=[Counter() for _ in range(5)];chars[0][0,0]=1
    for w in weights:
        for d in range(1,5):
            for mu,n in list(chars[d-1].items()):chars[d][mu[0]+w[0],mu[1]+w[1]]+=n
    roots=[(1,-1),(0,1),(1,0),(1,1)]
    for d in [0,2,4]:
        multiplicity=0
        for bits in product((0,1),repeat=4):
            target=tuple(5+sum(b*r[k] for b,r in zip(bits,roots)) for k in range(2))
            multiplicity+=(-1)**sum(bits)*sum(chars[d].get((target[0]-w[0],target[1]-w[1]),0) for w in weights)
        assert multiplicity==0,(d,multiplicity)

def main():
    d2,b2=load(2);d4,b4=load(4)
    source=PACKED['source']
    b={(a,m):F(v) for a,m,v in source}
    constraints({():b},0);highest({():b},0)
    real_norm=verify_real_part(b)
    assert real_norm==F(3,131072)
    verify_no_H2()
    for degree,poly in [(2,b2),(4,b4)]:constraints(poly,degree);highest(poly,degree)
    rhs1={}
    for r in range(14):
        for i in range(14):
            for (a,m),v in R(b,r,i).items():
                wed=wedge(i,m)
                if wed:add(rhs1,((r,),a,wed[0]),-F(1,2)*wed[1]*v)
    assert delta(b2)==rhs1
    print('Exact degree-one closure passed','seconds',time.time()-tick,flush=True)
    T=[defaultdict(dict) for _ in range(14)]
    for r in range(14):
        for i in range(14):
            for d,row in enumerate(cov(r,i)):
                for s,c in row:add(T[i][union((r,),(s,))],d,-c)
    rhs3={}
    for mon,tensor in b2.items():
        for r in range(14):
            new=union(mon,(r,))
            for i in range(14):
                for (a,m),v in R(tensor,r,i).items():
                    wed=wedge(i,m)
                    if wed:add(rhs3,(new,a,wed[0]),-F(1,2)*wed[1]*v)
        for d,count in Counter(mon).items():
            reduced=list(mon);reduced.remove(d);reduced=tuple(reduced)
            for i in range(14):
                for q,ds in T[i].items():
                    c=ds.get(d,F(0))
                    if not c:continue
                    new=union(reduced,q)
                    for (a,m),v in tensor.items():
                        wed=wedge(i,m)
                        if wed:add(rhs3,(new,a,wed[0]),F(1,6)*count*c*wed[1]*v)
    for i in range(14):
        for q,ds in T[i].items():
            for e,c in ds.items():
                for r in range(14):
                    new=union(q,(r,))
                    for (a,m),v in R(b,r,e).items():
                        wed=wedge(i,m)
                        if wed:add(rhs3,(new,a,wed[0]),F(1,24)*c*wed[1]*v)
    lhs3=delta(b4);assert lhs3==rhs3,('cubic mismatch',len(lhs3),len(rhs3),next((k for k in set(lhs3)|set(rhs3) if lhs3.get(k,0)!=rhs3.get(k,0)),None))
    result=dict(exact_rational_fourth_jet=True,highest_value=[5,5],value_module='H2-free K6820',primes_used_for_proposal=d4['primes'],degree2_entries=sum(map(len,b2.values())),degree4_entries=sum(map(len,b4.values())),nonzero_cubic_closure_coordinates=len(lhs3),exact_checks=['Every polynomial coefficient: Alt,Ctr,coefficient-curvature wedge','Both complete highest raising identities at degrees2 and4','Every degree-one closure coordinate','Every degree-three closure coordinate'],seconds=time.time()-tick,scope='A genuine characteristic-zero even four-jet of the linear closure system at a nonzero (5,5) value. Taking compact real part gives real finite Taylor data. No actual closed field or nonlinear solution is claimed.')
    result['compact_real_value_norm_squared']=str(real_norm)
    assert result['degree2_entries']==PACKED['receipt']['degree2_entries'] and result['degree4_entries']==PACKED['receipt']['degree4_entries']
    print(json.dumps(result,indent=2),flush=True)

if __name__=='__main__':main()
