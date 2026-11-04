"""Exact Clifford controls for the three complementary even hook degrees.

The full derivative, its lower-grade cancellation and its contraction are
checked on a spanning SO14 orbit cover. The hand proof and field-sector
limits are in BA-hook-closure.md. Running this helper writes no files.
"""
from itertools import combinations
from math import comb
from time import perf_counter

N = 14
FULL = (1 << N)-1


def add(a, b, factor=1):
    out = dict(a)
    for k, v in b.items():
        out[k] = out.get(k, 0)+factor*v
    return {k:v for k,v in out.items() if v}


def scale(a, factor):
    return {k:v*factor for k,v in a.items() if v*factor}


def mul(a, b):
    out = {}
    for x,u in a.items():
        for y,v in b.items():
            parity = sum((y & ((1 << j)-1)).bit_count()
                         for j in range(N) if x & (1 << j))
            out[x^y] = out.get(x^y,0)+u*v*(-1)**parity
    return {k:v for k,v in out.items() if v}


def wedge(a, b):
    out = {}
    for x,u in a.items():
        for y,v in b.items():
            if not x & y:
                out = add(out,mul({x:u},{y:v}))
    return out


def interior(i, a):
    return {k^(1 << i):v*(-1)**((k & ((1 << i)-1)).bit_count())
            for k,v in a.items() if k & (1 << i)}


def gen(k, kind):
    K = {sum(1 << j for j in range(2,k+1)):1}
    B = [{} for _ in range(N)]
    if kind == 'D':
        B[0],B[1] = wedge({2:1},K),wedge({1:1},K)
    else:
        B[0],B[1] = wedge({1:1},K),scale(wedge({2:1},K),-1)
    alt,ctr = {},{}
    for a in range(N):
        alt = add(alt,wedge({1 << a:1},B[a]))
        ctr = add(ctr,interior(a,B[a]))
    assert not alt and not ctr
    return B


def control(k, kind, xi):
    B = gen(k,kind)
    D = B[xi]
    direct = []
    for a in range(N):
        raw,lower = {},{}
        for b,c in combinations(range(N),2):
            H = add(B[c] if b == xi else {},B[b] if c == xi else {},-1)
            if not H:
                continue
            lower = add(lower,scale(interior(a,interior(b,interior(c,H))),2))
            if a in (b,c):
                continue
            triple = mul(mul({1 << a:1},{1 << b:1}),{1 << c:1})
            raw = add(raw,add(mul(triple,H),mul(H,triple),-1))
        target = scale(add(wedge({1 << xi:1},B[a]),
                           wedge({1 << a:1},D),-1),2)
        assert raw == target
        assert not lower
        assert all(m.bit_count() == k+1 for m in raw)
        actual = scale(mul(raw,{FULL:1}),-1)
        assert all(m.bit_count() == 13-k for m in actual)
        direct.append(raw)
    contraction = {}
    for a in range(N):
        contraction = add(contraction,interior(a,direct[a]))
    assert contraction == scale(D,-2*(13-k))


def main():
    started = perf_counter()
    assert mul({FULL:1},{FULL:1}) == {0:-1}
    count = 0
    for k in (2,6,10):
        dimension = k*comb(N,k+1)+(N-k)*comb(N,k-1)
        assert dimension == N*comb(N,k)-comb(N,k+1)-comb(N,k-1)
        assert dimension == {2:896,6:36608,10:11648}[k]
        for kind in ('D','R'):
            for xi in range(N):
                control(k,kind,xi);count+=1
            print('PASS k',k,'type',kind,': all14 covectors, every output slot, lower-grade cancellation and contraction')
        # The triple-contraction term exists before imposing the hook jet.
        if k >= 6:
            blade={sum(1 << j for j in range(k)):1}
            triple={7:1}
            raw=add(mul(triple,blade),mul(blade,triple),-1)
            assert {m.bit_count() for m in raw} == {k-3}
            assert raw == scale(interior(0,interior(1,interior(2,blade))),2)
        print('Hk dimension',dimension,'; converted grade',k+1,
             '; actual derivative grade',13-k,'; contraction coefficient',-2*(13-k))
    assert count == 84
    print('84 exact orbit controls passed; generator span and SO14 covariance supply universality.')
    print('Lower-grade terms on unconstrained inputs are not silently dropped.')
    print('No curvature rank, quadratic coefficient or full-field existence claim.')
    print(f'Runtime {perf_counter()-started:.3f} seconds.')


if __name__ == '__main__':
    main()
