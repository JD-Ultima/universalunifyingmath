"""Exact grade-nine obstruction for the certified BL mixed point value.

The complete six-dimensional polynomial family is lifted using rational
SO5 plane rotations and signed coordinate-permutation averaging. The
selected higher quadratic receiver is separated from zero on BL's real
certificate box. All tensors are rebuilt in memory; adjacent modules are
libraries and no earlier main is called. Running writes only stdout.
"""
from fractions import Fraction as F
from functools import lru_cache
from itertools import product, combinations_with_replacement, combinations, permutations
from math import prod, factorial
from pathlib import Path
from time import perf_counter
import json
import sys
import sympy as s

HERE=Path(__file__).resolve().parent
PUBLIC=HERE
sys.path.insert(0,str(PUBLIC))
import BC_affine_hooks as BC
import BF_h6_casimir as BF
import BG_h6_curvature as BG
import BC_quadratic_xcheck as X

G=tuple(BC.METRIC)
FULL=(1<<14)-1
DET=prod(G)
PHI=[s.Matrix(m)/5 for m in BG.integer_phi()]
GRADE9_TARGET=(0,1279)
GRADE9_ROW=(0,570343710720,0,0,-51352963200,-3736212480,
            35548651520,49612780800,5095537920,1183248000,0,
            0,0,11396602800,157894110,0,337252410,-8742492,1884938121,0,0)
GRADE9_SCALE=s.Rational(5,1672704)
AMPLITUDE_PAIRS=tuple(combinations_with_replacement(range(6),2))


def rotation_coefficients(degree):
    c,t=s.Rational(4,5),s.Rational(3,5)
    cos=[s.Integer(1)];sin=[s.Integer(0)]
    for _ in range(degree):
        cos.append(s.expand(c*cos[-1]-t*sin[-1]))
        sin.append(s.expand(t*cos[-2]+c*sin[-1]))
    even=s.Matrix([[(-k*k)**m for m in range(degree+1)] for k in range(degree+1)]).inv()*s.Matrix(cos)
    odd=s.Matrix([[(-k*k)**m for m in range(degree)] for k in range(1,degree+1)]).inv()*s.Matrix([sin[k]/k for k in range(1,degree+1)])
    coeff=[F(0)]*(2*degree+1)
    for m,v in enumerate(even):coeff[2*m]=F(v)
    for m,v in enumerate(odd):coeff[2*m+1]=F(v)
    for k in range(-degree,degree+1):
        expected=(c+s.I*t)**k
        assert s.simplify(sum(s.Rational(v)*(s.I*k)**j for j,v in enumerate(coeff))-expected)==0
    return coeff


def rotate(tensor,degree,axis):
    coeff=rotation_coefficients(degree)
    powers=tensor
    result={}
    for order,c in enumerate(coeff):
        for key,v in powers.items():BC.plus(result,key,c*v)
        if order<len(coeff)-1:powers={key:2*v for key,v in BC.action(6+axis,powers).items()}
    return result


@lru_cache(None)
def swap_rows(i,j):
    p=list(range(5));p[i],p[j]=p[j],p[i]
    images=[]
    for a in range(14):
        matrix=s.Matrix(5,5,lambda r,c:PHI[a][p[r],p[c]])
        images.append([s.trace(matrix*PHI[d])/s.Rational(G[d]) for d in range(14)])
    tangent=s.Matrix.hstack(*(s.Matrix(v) for v in images))
    assert tangent*tangent==s.eye(14)
    assert tangent.T*s.diag(*map(s.Rational,G))*tangent==s.diag(*map(s.Rational,G))
    # Source covector a maps to sum_b X^{-1}[a,b] e^b. Here X^{-1}=X.
    return tuple(tuple((b,F(tangent[a,b])) for b in range(14) if tangent[a,b]) for a in range(14))


@lru_cache(None)
def swap_form(i,j,mask):
    rows=swap_rows(i,j)
    out={0:F(1)}
    for source in BC.bit_indices(mask):
        nxt={}
        for existing,value in out.items():
            for target,c in rows[source]:
                if not existing>>target&1:
                    sign=(-1)**((existing>>(target+1)).bit_count())
                    BC.plus(nxt,existing|(1<<target),value*c*sign)
        out=nxt
    return tuple(out.items())


def swap(tensor,i,j):
    rows=swap_rows(i,j)
    out={}
    for (a,mask),v in tensor.items():
        for b,c in rows[a]:
            for target,d in swap_form(i,j,mask):BC.plus(out,(b,target),v*c*d)
    return out


def reynolds(tensor):
    # Signed ambient S5 projection = natural S5 projection under the
    # determinant-twisted odd-harmonic SO5 intertwiner.
    result=tensor
    for n in range(2,6):
        result=BF.combine((F(1,n),result),*((-F(1,n),swap(result,j,n-1)) for j in range(n-1)))
    return result


def full_casimir(tensor):
    out=BF.casimir(tensor)
    for j in range(6):out=BF.combine((1,out),(-1,BC.action(j,BC.action(j,tensor))))
    return out


def polynomial_coordinates():
    """Derive the small Reynolds matrices from the actual rotation directions."""
    import BL_symmetric_mixed as FAMILY
    x=s.symbols('x0:5')
    radius=sum(v*v for v in x)
    parts={1:((1,0,0,0,0),),3:FAMILY.PARTS3,5:FAMILY.PARTS5}
    orbits={d:[FAMILY.orbit(part,x) for part in group] for d,group in parts.items()}
    directions=[]
    for length in range(5):
        n=[s.Integer(0)]*5
        n[4]=s.Integer(1)
        for i in range(length):
            n[i]=s.Rational(3,5)*n[4]
            n[4]*=s.Rational(4,5)
        assert sum(v*v for v in n)==1
        directions.append(tuple(n))
    def average_power(n,degree):
        images=tuple(permutations(n))
        result=0
        for part,polynomial in zip(parts[degree],orbits[degree]):
            c=s.Rational(factorial(degree),prod(factorial(p) for p in part))
            c*=sum(prod(v**p for v,p in zip(image,part)) for image in images)/120
            result+=c*polynomial
        return s.expand(result)
    columns={3:[],5:[]}
    bases={3:s.Matrix(FAMILY.BASIS3).T,5:s.Matrix(FAMILY.BASIS5).T}
    for n in directions:
        p1,p3,p5=(average_power(n,d) for d in (1,3,5))
        harmonics={3:s.expand(p3-s.Rational(3,7)*radius*p1),
                   5:s.expand(p5-s.Rational(10,11)*radius*p3+s.Rational(5,33)*radius**2*p1)}
        for d,value in harmonics.items():
            assert s.expand(sum(s.diff(value,v,2) for v in x))==0
            p=s.Poly(value,*x,domain=s.QQ)
            coordinate=s.Matrix([p.coeff_monomial(part) for part in parts[d]])
            solution,free=bases[d].gauss_jordan_solve(coordinate)
            assert not free and bases[d]*solution==coordinate
            columns[d].append(solution)
    selected={3:(0,2),5:(0,1,2,4)}
    determinants={3:s.Rational(1728,546875),
                  5:s.Rational(71230958862336,7212162017822265625)}
    inverses={}
    print('Rational unit directions:',directions)
    for d in (3,5):
        matrix=s.Matrix.hstack(*columns[d])
        assert tuple(matrix.rref()[1])==selected[d]
        square=matrix[:,selected[d]]
        assert square.det()==determinants[d]
        inverse=square.inv()
        assert square*inverse==s.eye(bases[d].cols)
        inverses[d]=inverse
        print('Degree',d,'Reynolds columns:',matrix.tolist())
        print('Degree',d,'selected directions:',selected[d],'; determinant:',square.det())
    return inverses,selected


def reconstruct():
    receipt=json.loads((PUBLIC/'BC_h6_closed_generator.json').read_text(encoding='utf-8'))
    assert receipt['degree']==6 and tuple(map(F,receipt['metric_gram']))==G
    original={(a,m):F(v) for a,m,v in receipt['entries']}
    assert len(original)==len(receipt['entries'])==1358
    casimir=BF.casimir(original)
    seeds={3:BF.combine((F(40,22),original),(-F(1,22),casimir)),
           5:BF.combine((F(1,22),casimir),(-F(18,22),original))}
    inverses,indices=polynomial_coordinates()
    basis=[]
    for degree in (3,5):
        prefix=seeds[degree]
        averaged={}
        for length in range(max(indices[degree])+1):
            if length:
                prefix=rotate(prefix,degree,length-1)
            if length in indices[degree]:
                averaged[length]=reynolds(prefix)
        vectors=[averaged[j] for j in indices[degree]]
        for column in range(len(vectors)):
            value=BF.combine(*((F(inverses[degree][row,column]),vectors[row]) for row in range(len(vectors))))
            assert BC.hook_traces(value)==({},{})
            assert not BF.combine((1,full_casimir(value)),(-degree*(degree+3),value))
            for j in range(4):assert not BF.combine((1,value),(1,swap(value,j,j+1)))
            basis.append(value)
            print("Basis",len(basis)-1,"hooks, Casimir and signed-S5 action verified; support",len(value))
    # Compare the complete norm Gram, including every cubic/quintic cross.
    expected=s.zeros(6)
    expected[:2,:2]=s.diag(s.Rational(76899375,88),s.Rational(3661875,352))
    expected[2:,2:]=s.Rational(1366875,32768)*s.Matrix([[380,20,0,0],[20,8,0,0],[0,0,33,0],[0,0,0,s.Rational(1,120)]])
    actual=s.Matrix(6,6,lambda i,j:s.Rational(BF.inner(basis[i],basis[j])))
    assert actual==expected
    print('Full 6x6 norm Gram matches the BL/BI/BK normalization exactly.')
    return basis


def target_probe(basis):
    import BI_h3_stress as CUBIC
    import BK_quintic_gram as QUINTIC
    import BL_mixed_stress as MIXED
    import BL_symmetric_mixed as FAMILY
    import BL_real_certificate as CERTIFICATE
    algebra=BC.Clifford(G)
    def convert(tensor):
        return [algebra.mul({m:v for (b,m),v in tensor.items() if b==a},{FULL:F(1)}) for a in range(14)]
    values=list(map(convert,basis))
    assert all(mask.bit_count()==8 for value in values for item in value for mask in item)
    def combine(weights):
        return [BC.add(*(BC.scale(value[a],weight) for weight,value in zip(weights,values))) for a in range(14)]
    def q(field,a,mask):return DET*BC.quadratic_coefficient(field,G,a,mask)
    x=CUBIC.x
    poly3=sum(c*FAMILY.orbit(part,x) for c,part in zip(FAMILY.BASIS3[0],FAMILY.PARTS3))
    poly5=sum(c*FAMILY.orbit(part,x) for c,part in zip(FAMILY.BASIS5[0],FAMILY.PARTS5))
    predictions=(F((s.Matrix(FAMILY.PURE3).T*CUBIC.patterns(CUBIC.tensors(poly3),PHI[0]))[0]),
                 F((s.Matrix(FAMILY.PURE5).T*QUINTIC.patterns(QUINTIC.tensor_blocks(poly5),PHI[0]))[0]))
    observed=(q(values[0],0,1),q(values[2],0,1))
    assert predictions==observed
    crossdata=MIXED.cross_data(MIXED.polynomial_tensor(poly3,3),MIXED.polynomial_tensor(poly5,5))
    prediction=F((s.Matrix(MIXED.CROSS_COEFFICIENTS).T*MIXED.patterns(crossdata,PHI[0],PHI[0]))[0])
    actual=q(combine((1,0,1,0,0,0)),0,1)-sum(observed)
    assert actual==prediction and actual!=0
    print('Actual grade-one pure/cross guards agree with BI/BK/BL:',observed,actual,flush=True)
    selected=GRADE9_TARGET
    a,mask=selected
    assert mask.bit_count()==9
    diagonal=[q(value,a,mask) for value in values]
    pairs=AMPLITUDE_PAIRS
    coeff=[]
    for i,j in pairs:
        if i==j:coeff.append(diagonal[i])
        else:
            weights=[0]*6;weights[i]=weights[j]=1
            coeff.append(q(combine(weights),a,mask)-diagonal[i]-diagonal[j])
    den=s.ilcm(*(v.denominator for v in coeff))
    ints=[int(v*den) for v in coeff]
    common=s.igcd(*ints)
    row=[v//common for v in ints]
    scale=s.Rational(common,den)
    assert tuple(row)==GRADE9_ROW and scale==GRADE9_SCALE
    amplitudes=s.symbols('u0:6')
    polynomial=sum(c*amplitudes[i]*amplitudes[j] for c,(i,j) in zip(row,pairs))
    center=[s.Rational(n,CERTIFICATE.CENTER_DENOMINATOR) for n in CERTIFICATE.CENTER_NUMERATORS]
    restricted=polynomial.subs({amplitudes[4]:-s.Rational(1,100),amplitudes[5]:1})
    atcenter=restricted.subs(dict(zip(amplitudes[:4],center)))
    radius=s.Rational(1,10**10)
    gradient=[s.diff(restricted,v).subs(dict(zip(amplitudes[:4],center))) for v in amplitudes[:4]]
    hessian=s.hessian(restricted,amplitudes[:4])
    bound=radius*sum(abs(v) for v in gradient)+radius**2*sum(abs(v) for v in hessian)/2
    assert not polynomial.atoms(s.Float)
    independent=X.BladeAlgebra(G)
    control=combine((2,-3,1,0,1,-1))
    expected=scale*polynomial.subs(dict(zip(amplitudes,(2,-3,1,0,1,-1))))
    assert DET*independent.quadratic_coefficient(control,a,mask)==expected
    assert abs(atcenter)>bound and atcenter < -10**7 and bound < 2
    print('Selected receiver:',selected,'grade',mask.bit_count(),flush=True)
    print('Amplitude monomials:',pairs)
    print('Primitive 21-coefficient row:',row)
    print('Physical coordinate quadratic = scale * primitive polynomial; scale:',scale)
    print('Primitive polynomial at the BL certificate center:',atcenter)
    print('Exact certificate-box variation bound:',bound)
    print('Nonzero on the entire certified box:',abs(atcenter)>bound)
    print('Exact coarse guards: center value < -10000000; variation bound < 2.')
    print('One independent nested Clifford target agrees exactly.')
    print('The certified scalar-stress point fails this grade-nine equation.')
    print('No exclusion of all mixed scalar zeros or unrestricted fields is inferred.')


def main():
    started=perf_counter()
    basis=reconstruct()
    target_probe(basis)
    print(f'Runtime {perf_counter()-started:.3f} seconds.')


if __name__=='__main__':
    main()
