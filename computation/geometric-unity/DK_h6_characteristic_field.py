"""Exact local mixed-field lift of a balanced real characteristic value."""
import argparse
from fractions import Fraction as F
from itertools import combinations
import json
from pathlib import Path
import time

import sympy as s

import BO_h6_full_kernel as BO
from CU_h6_local_field import tangent_matrices


HERE=Path(__file__).resolve().parent


def geometry():
    _, inverse, curvature, roots = BO.build_model()
    metric = s.zeros(5)
    metric[0, 1] = metric[1, 0] = metric[2, 3] = metric[3, 2] = 2
    metric[4, 4] = 1
    tangent = tangent_matrices(metric)
    compact = [(BO.matrix_units(i, j) - BO.matrix_units(j, i)) * metric
               for i, j in combinations(range(5), 2)]
    conjugation = [4, 6, 5, 7, 0, 2, 1, 3, 10, 11, 8, 9, 12, 13]
    return tangent, compact, inverse, curvature, roots, conjugation, metric


def verify():
    p, k, inv, curv, roots, tau, metric = geometry()
    result, v, leaf, stabilizer = verify_leaf(p, k, inv, curv, roots, tau)
    return result, p, metric, v, leaf, stabilizer


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-receipt", action="store_true")
    args = parser.parse_args()
    start = time.perf_counter()
    result, *_ = verify()
    receipt = HERE / "DK_h6_characteristic_field.json"
    if args.write_receipt:
        receipt.write_text(json.dumps(result, indent=2), encoding="utf-8")
    else:
        saved = json.loads(receipt.read_text(encoding="utf-8"))
        assert result == {key: value for key, value in saved.items() if key != "seconds"}
    print(json.dumps({**result, "seconds": time.perf_counter() - start}, indent=2))


def verify_leaf(p, k, inv, curv, roots, tau):
    v=sum((s.Rational(inv[1][j])*p[j] for j in range(14)),s.zeros(5))
    D=[p[i] for i in [2,4,5,6,7,10,11]]+[p[12]+p[13]]
    T=[p[i] for i in [0,1,3,8,9]]+[(p[12]-p[13])/2]
    def flat(a):return s.Matrix(list(a))
    def span(items):
        A=s.Matrix.hstack(*map(flat,items));return [items[j] for j in A.rref()[1]]
    def contained(a,items):
        A=s.Matrix.hstack(*map(flat,items));return A.row_join(flat(a)).rank()==A.cols
    assert s.Matrix.hstack(*map(flat,D+T)).rank()==14
    comm=span([a*b-b*a for a in D for b in D if a*b!=b*a])
    conditions=s.Matrix.hstack(*[flat(a*v-v*a) for a in k])
    kc=[sum((c[i]*k[i] for i in range(10)),s.zeros(5)) for c in conditions.nullspace()]
    assert len(comm)==len(kc)==4
    assert all(a*v==v*a for a in comm)
    assert all(contained(a*b-b*a,D) for a in kc for b in D)
    assert all(contained(a*b-b*a,kc) for a in kc for b in kc)
    beta={4875:F(1),8971:F(-1)};B={(1,m):z for m,z in beta.items()}
    constraints={}
    for (a,m),z in B.items():
        for key,c in BO.constraint_entries(a,m,inv,curv):BO.add(constraints,key,z*c)
    assert not constraints
    def wedge(vec):
        out={}
        for (a,m),z in vec.items():
            for i,c in [(12,1),(13,-1)]:
                if not m>>i&1:BO.add(out,(a,m|(1<<i)),c*BO.wedge_sign(i,m)*z)
        return out
    assert not wedge(B)
    source=HERE;data=json.loads((source/'BO_h6_full_kernel_witnesses.json').read_text())['extremal_vectors']
    def seed(w):return {(a,m):F(z) for rec in data if tuple(rec['weight'])==w for a,m,z in rec['vectors'][0]}
    lower=[[] for _ in range(14)]
    for i,row in enumerate(roots[0]):lower[tau[i]]=[(tau[j],z) for j,z in row]
    b55=seed((5,5));b73=BO.root_action(BO.root_action(seed((7,3)),lower),lower)
    combo=dict(b73)
    for key,z in b55.items():BO.add(combo,key,-2*z)
    assert combo=={key:6*z for key,z in B.items()}
    # Compact conjugation sends the leaf direction p4 to p0, outside D.
    assert not contained(p[0],D)
    real={key:z/2 for key,z in B.items()}
    for (a,m),z in B.items():
        seq=[tau[i] for i in BO.axes(m)];sign=(-1)**sum(seq[i]>seq[j] for i in range(len(seq)) for j in range(i+1,len(seq)))
        BO.add(real,(tau[a],sum(1<<i for i in seq)),sign*z/2)
    N=F(0)
    for (a,m),z in real.items():
        for (b,n),w in real.items():
            if inv[a][b]:N+=z*w*inv[a][b]*F(s.Matrix([[s.Rational(inv[i][j]) for j in BO.axes(n)] for i in BO.axes(m)]).det())
    assert N==F(1,131072)
    result=dict(leaf_dimension=8,compact_coefficient_centralizer_dimension=4,complex_h_dimension=12,transverse_dimension=6,real_form_preserves_h=False,origin_entries=[[a,m,str(z)] for (a,m),z in sorted(B.items())],compact_real_origin_norm_squared=str(N),origin_decomposition='(F_alpha1^2 b73 - 2 b55)/6',scope='Local holomorphic field restricted to the real symmetric space; taking its real part preserves all linear equations. No global extension is asserted.')
    return result, v, D, kc


if __name__ == "__main__":
    main()
