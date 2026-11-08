"""Exact scalarity reduction on the complete coordinate-S5 mixed family.

BI/BK and BL's calibrated cross formula are inputs. All 196 ordered receiver
entries and all 21 amplitude coefficients are recomputed. The four integer
rows and monomial order are exported for BL_real_certificate.py. No earlier
main or numerical solver is run; this helper writes stdout only.
"""
from pathlib import Path
from itertools import product, permutations, combinations_with_replacement
from math import factorial, prod, gcd, lcm
from functools import reduce
from time import perf_counter
import hashlib
import json
import sys
import sympy as s

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import BG_h6_curvature as BG
import BC_affine_hooks as BC
import BL_mixed_stress as MIXED

PARTS3=((3,0,0,0,0),(2,1,0,0,0),(1,1,1,0,0))
PARTS5=((5,0,0,0,0),(4,1,0,0,0),(3,2,0,0,0),(3,1,1,0,0),(2,2,1,0,0),(2,1,1,1,0),(1,1,1,1,1))
BASIS3=((4,-3,0),(0,0,1))
BASIS5=((4,5,-10,0,0,0,0),(0,1,0,0,-2,0,0),(0,0,0,2,0,-9,0),(0,0,0,0,0,0,1))
PURE3=tuple(map(s.Rational,('-9124875/968','4095875/484','-268625/484','1031625/484')))
PURE5=tuple(map(s.Rational,('-2311875/32768','2244375/32768','1873125/32768','860625/16384')))
CROSS=MIXED.CROSS_COEFFICIENTS
NORM_RATIOS=MIXED.NORM_RATIOS
AMPLITUDE_PAIRS=tuple(combinations_with_replacement(range(6),2))
QUADRATIC_ROWS=(
    (252019998720,0,81589939200,-10878658560,0,0,-2000158720,0,0,1359832320,0,32566762800,3674198880,0,0,267214464,0,0,-1440452970,0,-695871),
    (84608778240,92408709120,60513868800,11714250240,335301120,0,-8700805120,1942617600,-475960320,-2821737600,19831680,-4225211100,-2692138680,1002054240,0,-302615676,-1177060896,0,-1406139669,49044204,0),
    (0,6000476160,0,0,-437088960,48565440,-250019840,0,194261760,437088960,0,0,0,-1377824580,0,0,-175359492,4175226,-231725043,0,0),
    (28082503680,12721192960,5348851200,1134524160,2259290880,0,-893419520,-323769600,-141166080,-521928000,5881920,-1064388600,-258746400,-167009040,0,-20307672,-86797656,0,-78094368,3616569,0),
)
ROW_JSON_SHA256='8d365144715aac0307b26b0051d5c59610bc6ead63a903f66d727d005788d6bb'


def canonical_rows_bytes():
    """Stable integer data only; independent of source formatting or timing."""
    return json.dumps(QUADRATIC_ROWS,separators=(',',':')).encode('ascii')


def quadrics(amplitudes):
    assert len(amplitudes)==6
    return s.Matrix([sum(c*amplitudes[i]*amplitudes[j] for c,(i,j) in zip(row,AMPLITUDE_PAIRS))
                     for row in QUADRATIC_ROWS])


def check_normalizations():
    """Retain the exact pure-module constants and physical norm conventions."""
    assert PURE3==(-s.Rational(9124875,968),s.Rational(4095875,484),
                   -s.Rational(268625,484),s.Rational(1031625,484))
    assert PURE5==(-s.Rational(2311875,32768),s.Rational(2244375,32768),
                   s.Rational(1873125,32768),s.Rational(860625,16384))
    assert CROSS==(-s.Rational(2131875,2816),s.Rational(2874375,11264),s.Rational(2874375,11264))
    assert NORM_RATIOS==(s.Rational(2197125,352),s.Rational(1366875,32768))
    for (alpha,beta,gamma,delta),ratio in zip((PURE3,PURE5),NORM_RATIOS):
        scalar=14*alpha+s.Rational(14,5)*beta+gamma+s.Rational(3,10)*delta
        assert scalar==-s.Rational(52,3)*ratio
    assert hashlib.sha256(canonical_rows_bytes()).hexdigest()==ROW_JSON_SHA256


def orbit(part,x):
    return sum(prod(v**n for v,n in zip(x,p)) for p in sorted(set(permutations(part))))


def components(poly,degree,x):
    p=s.Poly(poly,*x,domain=s.QQ)
    return {indices:p.coeff_monomial(tuple(indices.count(i) for i in range(5)))
            *prod(factorial(indices.count(i)) for i in range(5))/factorial(degree)
            for indices in product(range(5),repeat=degree)}


def primitive(row):
    den=lcm(*(v.q for v in row))
    nums=[int(v*den) for v in row]
    g=reduce(gcd,nums)
    nums=[v//g for v in nums]
    return tuple(-v for v in nums) if next(v for v in nums if v)<0 else tuple(nums)


def main():
    start=perf_counter()
    check_normalizations()
    x=s.symbols('x0:5')
    groups=[]
    for degree,parts,columns in ((3,PARTS3,BASIS3),(5,PARTS5,BASIS5)):
        orbits=[orbit(p,x) for p in parts]
        laps=[s.Poly(sum(s.diff(p,v,2) for v in x),*x) for p in orbits]
        powers=sorted(set().union(*(set(p.monoms()) for p in laps)))
        matrix=s.Matrix([[p.coeff_monomial(m) for p in laps] for m in powers])
        assert matrix.rank()==len(parts)-len(columns)
        assert s.Matrix(columns).rank()==len(columns)
        basis=[s.expand(sum(c*p for c,p in zip(column,orbits))) for column in columns]
        tensors=[components(p,degree,x) for p in basis]
        for p,T in zip(basis,tensors):
            assert sum(s.diff(p,v,2) for v in x).expand()==0
            assert not p.atoms(s.Float)
            assert all(sum(T[(i,i)+tail] for i in range(5))==0
                       for tail in product(range(5),repeat=degree-2))
        groups.append(tensors)
    phi=[s.Matrix(m)/5 for m in BG.integer_phi()]
    sparse=[[(i,j,m[i,j]) for i,j in product(range(5),repeat=2) if m[i,j]] for m in phi]
    G=tuple(map(s.Rational,BC.METRIC))
    assert all(s.trace(phi[a]*phi[b])==(G[a] if a==b else 0)
               for a,b in product(range(14),repeat=2))
    receivers=tuple(product(range(14),repeat=2))
    pairs=AMPLITUDE_PAIRS
    coefficients={receiver:[] for receiver in receivers}
    norm_coefficients=[]
    for i,j in pairs:
        if (i<2)==(j<2):
            group=0 if i<2 else 1
            T,U=groups[group][i if group==0 else i-2],groups[group][j if group==0 else j-2]
            degree=3 if group==0 else 5
            tails=tuple(product(range(5),repeat=degree-2))
            K={(a,b,c,d):sum(T[(a,b)+tail]*U[(c,d)+tail]
                            +(T[(c,d)+tail]*U[(a,b)+tail] if i!=j else 0)
                            for tail in tails) for a,b,c,d in product(range(5),repeat=4)}
            A=s.Matrix(5,5,lambda a,c:sum(K[a,b,c,b] for b in range(5)))
            N=s.trace(A)
            alpha,beta,gamma,delta=(PURE3,PURE5)[group]
            norm_coefficients.append(NORM_RATIOS[group]*N)
            for a,b in receivers:
                first=alpha*N*(G[a] if a==b else 0)+beta*s.trace(A*phi[a]*phi[b])
                rest=sum(v*w*(gamma*K[i,j,k,l]+delta*K[i,l,j,k])
                         for i,j,v in sparse[a] for k,l,w in sparse[b])
                coefficients[a,b].append(s.expand(first+rest))
        else:
            T,U=groups[0][i],groups[1][j-2]
            Z={(a,b,c,d):sum(T[i,j,a]*U[i,j,b,c,d] for i,j in product(range(5),repeat=2))
               for a,b,c,d in product(range(5),repeat=4)}
            C=s.Matrix(5,5,lambda c,d:sum(Z[a,a,c,d] for a in range(5)))
            assert C==C.T and C.trace()==0
            norm_coefficients.append(s.Integer(0))
            for a,b in receivers:
                value=CROSS[0]*s.trace(C*phi[a]*phi[b])+sum(
                    v*w*(CROSS[1]*Z[i,j,k,l]+CROSS[2]*Z[k,l,i,j])
                    for i,j,v in sparse[a] for k,l,w in sparse[b])
                coefficients[a,b].append(s.expand(value))
    assert len(norm_coefficients)==21
    assert all(not value.atoms(s.Float) for values in coefficients.values() for value in values)
    assert all(not value.atoms(s.Float) for value in norm_coefficients)
    assert all(coefficients[a,b]==coefficients[b,a] for a,b in receivers)
    trace=[sum(coefficients[a,a][k]/G[a] for a in range(14)) for k in range(21)]
    assert trace==[-s.Rational(52,3)*n for n in norm_coefficients]
    rows={}
    for a,b in receivers:
        row=[s.expand(v-(trace[k]*G[a]/14 if a==b else 0))
             for k,v in enumerate(coefficients[a,b])]
        if any(row):rows.setdefault(primitive(row),(a,b))
    M=s.Matrix(list(rows))
    _,pivot=M.T.rref()
    selected=tuple(list(rows)[i] for i in pivot)
    assert tuple(rows)==QUADRATIC_ROWS and selected==QUADRATIC_ROWS and M.rank()==4
    assert not M.atoms(s.Float)
    amplitudes=s.symbols('a0:6',real=True)
    equations=quadrics(amplitudes)
    norm_gram=s.zeros(6)
    for value,(i,j) in zip(norm_coefficients,pairs):
        norm_gram[i,j]=norm_gram[j,i]=value/(1 if i==j else 2)
    assert all(norm_gram[:size,:size].det()>0 for size in range(1,7))
    print('Complete coordinate-S5 harmonic spaces: cubic dimension 2, quintic dimension 4.',flush=True)
    print('Cubic columns in m3,m21,m111:',BASIS3)
    print('Quintic columns in degree-five partition order:',BASIS5)
    print('All 196 receiver coefficient vectors computed; full stress symmetric; weighted trace=-52/3*norm.')
    print('Amplitude order: a0,a1 cubic; a2,a3,a4,a5 quintic.')
    print('Amplitude monomials:',pairs)
    print('Intertwiner norm coefficients:',tuple(norm_coefficients))
    print('Distinct primitive nonscalar rows:',len(rows),'; exact rank:',M.rank())
    print('Independent receiver labels:',[rows[tuple(v)] for v in selected])
    print('Exact independent coefficient rows:',selected)
    for index,value in enumerate(equations):print('f'+str(index)+' =',value)
    print('Every one of 196 scalarity equations lies in this exact row span.',flush=True)
    print('Published BI/BK constants and norm ratios agree; the six-amplitude norm Gram is positive definite.')
    print('Canonical integer-row JSON SHA256:',ROW_JSON_SHA256)
    print('The real-zero certificate consumes these exported rows and their exact monomial order.')
    print('This is a pointwise scalarity reduction; higher grades and field reconstruction remain separate.')
    print(f'Runtime {perf_counter()-start:.3f} seconds.')


if __name__=='__main__':main()
