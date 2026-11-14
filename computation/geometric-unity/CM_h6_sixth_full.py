"""Direct raw-tensor verification of the sixth-order rational proposal."""
from collections import Counter,defaultdict
from fractions import Fraction as Q
from functools import lru_cache
import json,time,sys
import CF_h6_fourth_jet as v
import CM_h6_sixth_geometry as g
import CM_h6_sixth_h2 as h2
BO=v.BO;HERE=v.HERE;tick=time.time()

def configure_arithmetic(use_fraction=False,flint_dir=None):
    """FLINT accelerates exact arithmetic; --stdlib keeps Fraction only."""
    dual=h2.build()
    if use_fraction:return dual,'fractions.Fraction'
    if flint_dir:sys.path.insert(0,str(flint_dir))
    try:
        from flint import fmpq
    except ImportError:
        print('python-flint unavailable; using exact Fraction arithmetic',flush=True)
        return dual,'fractions.Fraction'
    global Q
    oldQ=Q
    def convert(x):return fmpq(x.numerator,x.denominator) if isinstance(x,oldQ) else fmpq(x)
    original_cov=v.cov
    @lru_cache(None)
    def cov(r,i):return tuple(tuple((b,convert(c)) for b,c in row) for row in original_cov(r,i))
    v.cov=cov;v.raising=[[[ (b,convert(c)) for b,c in row] for row in a] for a in v.raising]
    for vectors in [g.T,g.U,g.J,g.L]:
        for terms in vectors:
            for ds in terms.values():
                for key,c in ds.items():ds[key]=convert(c)
    for key,items in dual.items():dual[key]=[(i,convert(c)) for i,c in items]
    original_constraints=BO.constraint_entries
    @lru_cache(None)
    def constraints_column(a,m):return tuple((key,convert(c)) for key,c in original_constraints(a,m,v.inverse,v.curvature))
    def constraints(poly,degree):
        for mon,tensor in poly.items():
            out={}
            for (a,m),c in tensor.items():
                for key,z in constraints_column(a,m):BO.add(out,key,c*z)
            assert not out,(degree,mon,next(iter(out.items())))
        print('Exact hook/curvature membership',degree,'coefficients',sum(map(len,poly.values())),'seconds',time.time()-tick,flush=True)
    v.constraints=constraints;Q=fmpq
    return dual,'FLINT fmpq'

def derivative_term(out,poly,vector,scale):
    byder=defaultdict(list)
    for i,terms in enumerate(vector):
        for mon,ds in terms.items():
            for d,c in ds.items():byder[d].append((i,mon,scale*c))
    for mon,tensor in poly.items():
        for d,n in Counter(mon).items():
            reduced=g.subtract(mon,(d,))
            for i,q,c in byder[d]:
                new=g.union(reduced,q)
                for (a,m),z in tensor.items():
                    if not m>>i&1:BO.add(out,(new,a,m|(1<<i)),n*c*BO.wedge_sign(i,m)*z)

def curvature_term(out,poly,piece,scale):
    bycurve=defaultdict(list)
    for i,terms in enumerate(piece):
        for mon,ds in terms.items():
            for (r,e),c in ds.items():bycurve[r,e].append((i,mon,scale*c))
    for n,(mon,tensor) in enumerate(poly.items()):
        for (r,e),targets in bycurve.items():
            image=v.R(tensor,r,e)
            if not image:continue
            for i,q,c in targets:
                new=g.union(mon,q)
                for (a,m),z in image.items():
                    if not m>>i&1:BO.add(out,(new,a,m|(1<<i)),c*BO.wedge_sign(i,m)*z)

def main(certificate,use_fraction=False,flint_dir=None):
    selected_rhs=g.selected_rhs()
    dual,backend=configure_arithmetic(use_fraction,flint_dir)
    b6=certificate.expand(Q)
    assert len(b6) == 13817 and sum(map(len, b6.values())) == 1936021
    v.constraints(b6,6);v.highest(b6,6)
    for mon,tensor in b6.items():
        pair=[Q(0)]*14
        for key,c in tensor.items():
            for i,z in dual.get(key,[]):pair[i]+=c*z
        assert not any(pair),('H2 component',mon,pair)
    print('Direct H2 projection',len(b6),'monomials; seconds',time.time()-tick,flush=True)
    rhs=g.RHS()
    for poly in [rhs.b2,rhs.b4,{():rhs.b}]:
        for tensor in poly.values():
            for key,c in tensor.items():tensor[key]=Q(c.numerator,c.denominator)
    out={}
    derivative_term(out,rhs.b4,g.T,Q(1,6));print('Quintic T derivative complete',len(out),'seconds',time.time()-tick,flush=True)
    linear=[{(r,):{(r,i):Q(1)} for r in range(14) if r!=i} for i in range(14)]
    curvature_term(out,rhs.b4,linear,Q(-1,2));print('Quintic R b4 complete',len(out),'seconds',time.time()-tick,flush=True)
    derivative_term(out,rhs.b2,g.U,Q(-7,360));print('Quintic U derivative complete',len(out),'seconds',time.time()-tick,flush=True)
    curvature_term(out,rhs.b2,g.J,Q(1,24));print('Quintic R T b2 complete',len(out),'seconds',time.time()-tick,flush=True)
    curvature_term(out,{():rhs.b},g.L,Q(-1,240));print('Quintic RHS complete',len(out),'seconds',time.time()-tick,flush=True)
    # Independent backwards-coordinate RHS at the full fixed pivot set.
    for (mon,a,m),value in zip(g.AUX['closure_labels'],selected_rhs):
        assert out.get((tuple(mon),a,m),0)==Q(value.numerator,value.denominator)
    lhs=v.delta(b6)
    assert len(out) == 604246
    assert lhs==out,('quintic mismatch',len(lhs),len(out),next((k for k in set(lhs)|set(out) if lhs.get(k,0)!=out.get(k,0)),None))
    result=dict(exact_rational_sixth_jet=True,highest_value=[5,5],value_module='H2-free K6820',arithmetic_backend=backend,degree6_entries=sum(map(len,b6.values())),degree6_monomials=len(b6),nonzero_quintic_closure_coordinates=len(out),exact_checks=['Every degree6 coefficient: Alt,Ctr,coefficient-curvature wedge','Both complete highest raising identities at degree6','Direct invariant pairing with verified full14D H2 module at every coefficient','All quintic closure coefficients, forward operator','Independent backwards operator at986 selected coordinates'],seconds=time.time()-tick,scope='Characteristic-zero even six-jet of linear closure at a nonzero (5,5) value, extending the exact CF four-jet. No actual local field or nonlinear solution is claimed.')
    print(json.dumps(result,indent=2),flush=True)
