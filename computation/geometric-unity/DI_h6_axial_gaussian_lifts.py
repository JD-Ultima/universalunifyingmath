"""DI: exact compact orbit witnesses for Gaussian axial H6 profiles."""

from fractions import Fraction as F
from itertools import combinations
import json
from math import prod
from pathlib import Path
from time import perf_counter

import numpy as np

from DS_h6_axial_multisign import CV, archive, build_pencil, check_chain

HERE = Path(__file__).resolve().parent
P = 65537
IMAGINARY = 256


def compact_rotation(moves):
    Z = np.eye(5, dtype=np.int64)
    for i, j, t in moves:
        Q = np.eye(5, dtype=np.int64)
        c, s = F(1-t*t, 1+t*t), F(2*t, 1+t*t)
        assert c*c+s*s == 1 and 0 <= i < j < 5
        Q[i,i] = Q[j,j] = CV.finite(c)
        Q[i,j], Q[j,i] = CV.finite(s), CV.finite(-s)
        Z = Q @ Z % P
    assert np.array_equal(Z.T @ Z % P, np.eye(5, dtype=np.int64))
    return Z


def wedge_minor_matrix(T, degree):
    previous = {0: np.ones(1, dtype=np.int64)}
    for d in range(1, degree+1):
        masks = [sum(1 << i for i in b) for b in combinations(range(14), d)]
        loc = {m:i for i,m in enumerate(previous)}
        rows = np.array([[loc[m ^ (1 << i)] for i in range(14) if m >> i & 1] for m in masks])
        axes = np.array([[i for i in range(14) if m >> i & 1] for m in masks])
        earlier = np.column_stack(tuple(previous.values()))
        signs = (-1)**np.arange(d-1, -1, -1)
        previous = {}
        for m, row in zip(masks, axes):
            a = int(row[-1])
            sub = earlier[:, loc[m ^ (1 << a)]]
            previous[m] = np.einsum("ij,ij->i", sub[rows], T[axes, a]*signs) % P
    return masks, np.column_stack(tuple(previous.values()))


def compact_pair_rows(refs, rotation):
    # Hermitian rows already include the inverse Gram of all seven slots.
    # Pull a linear functional back to input tensors; hence rotate the
    # vector coefficients by the inverse, represented by Z transpose.
    Q = np.asarray(CV.BG.integer_phi(), dtype=np.int64)
    X = ((rotation.T @ Q % P) @ rotation) % P
    gram_inverse = np.array([CV.finite(1/(25*g)) for g in CV.G], dtype=np.int64)
    U = (Q.reshape(14, -1) @ X.reshape(14, -1).T % P * gram_inverse[:,None]) % P
    assert np.array_equal(U.T @ np.diag([CV.finite(g) for g in CV.G]) @ U % P,
                          np.diag([CV.finite(g) for g in CV.G]))
    masks, W = wedge_minor_matrix(U, 6)
    loc = {m:i for i,m in enumerate(masks)}
    # Audit independent literal determinants in one compact wedge block.
    axes = [(0,44),(401,45),(907,1000),(2098,2842),(3002,1378)]
    from itertools import permutations
    perms = [(w,(-1)**sum(w[i]>w[j] for i,j in combinations(range(6),2))) for w in permutations(range(6))]
    for i,j in axes:
        r = list(CV.BC.bit_indices(masks[i]))
        c = list(CV.BC.bit_indices(masks[j]))
        expected = sum(sign*prod(int(U[r[k], c[perm[k]]]) for k in range(6)) for perm,sign in perms)
        assert expected % P == W[i,j]
    out = []
    # Gaussian coefficients reduced mod P are small integers. The longest
    # matrix inner dimension is 3003 and 3003(P-1)^2 < 2^53. Binary64 BLAS
    # therefore performs these integer products and sums exactly.
    assert len(masks)*(P-1)**2 < 1<<53
    W = W.astype(float)
    Uf = U.astype(float)
    for sparse in refs:
        coeff = np.zeros((14,len(masks)))
        for (a,m),v in sparse.items():
            coeff[a,loc[m]] = v
        out.append(((Uf @ (coeff @ W.T % P)) % P).astype(np.int64))
    return out, loc


def chosen_columns(data):
    need = {tuple(label) for r in data["weights"] for label in r["labels"]}
    columns = {}
    for ch in sorted({label[1] for label in need}):
        filename = "CZ_h6_axial_local_raw_chains.json.gz" if ch == 0 else f"DS_chains_{ch:02d}.json.gz"
        bundle = archive(filename)
        coords,A,D = build_pencil(ch)
        loc = {x:i for i,x in enumerate(coords)}
        requested = sorted({chain for rot, c, chain, grade in need if c == ch})
        for chain in requested:
            d, z = bundle["chains"][chain]
            B=[{} for _ in range(d+1)]
            for k,a,m,x in z:
                assert 0 <= k <= d and loc[a,m] not in B[k]
                B[k][loc[a,m]]=F(x)
            check_chain(B,A,D)
            for rot, c, j, k in need:
                if c == ch and j == chain:
                    columns[rot,c,j,k] = {coords[i]:CV.finite(v) for i,v in B[k].items()}
    assert set(columns)==need
    return columns


def main():
    t = perf_counter()
    data = json.loads((HERE / "DI_axial_orbit_witnesses.json").read_text(encoding="utf-8"))
    assert P == CV.P and IMAGINARY*IMAGINARY%P == P-1
    highest = json.loads(Path(CV.BO.__file__).with_name("BO_h6_full_kernel_witnesses.json")
                         .read_text(encoding="utf-8"))["extremal_vectors"]
    assert sum(len(r["vectors"])*(l-m+1)*(l+m+2)*(2*l+3)*(2*m+1)//6
               for r in highest for l,m in [r["weight"]]) == 6834
    selected = {(tuple(r["weight"])):r for r in data["weights"]}
    assert set(selected) == {tuple(r["weight"]) for r in highest} - {(5,5),(2,0)}
    columns = chosen_columns(data)
    transform = CV.BR.BP.split_to_physical_covectors()
    refs = []
    for item in highest:
        if tuple(item["weight"]) in selected:
            for h in item["vectors"]:
                R,I = CV.BR.BP.transport(h, transform)
                refs.append({k:(CV.finite(R.get(k,0)/CV.physical(*k))
                                - IMAGINARY*CV.finite(I.get(k,0)/CV.physical(*k))) % P
                             for k in R.keys() | I.keys()})
    rotated, loc = compact_pair_rows(refs, compact_rotation(data["rotation"]))
    j = 0
    total = 0
    for h in highest:
        if tuple(h["weight"]) not in selected:
            continue
        labels = [tuple(x) for x in selected[tuple(h["weight"])]["labels"]]
        assert len(labels) == len(h["vectors"]) in (1,2)
        mat = []
        for i in range(len(labels)):
            row = []
            for label in labels:
                b = columns[label]
                if label[0] == 0:
                    row.append(sum(refs[j+i].get(k,0)*v for k,v in b.items()) % P)
                else:
                    assert label[0] == 1
                    row.append(sum(int(rotated[j+i][a,loc[m]])*v for (a,m),v in b.items()) % P)
            mat.append(row)
        rec = selected[tuple(h["weight"])]
        assert mat == rec["matrix"]
        det = mat[0][0] if len(mat)==1 else (mat[0][0]*mat[1][1]-mat[0][1]*mat[1][0])%P
        assert det == rec["determinant"] and det
        j += len(labels)
        l,m = h["weight"]
        total += len(labels)*(l-m+1)*(l+m+2)*(2*l+3)*(2*m+1)//6
        print(f"Weight {(l,m)}: {len(labels)} copies, exact Gaussian minor {det} mod {P}.", flush=True)
    assert j==len(refs) and total==6534
    print(f"Gaussian axial coefficient orbit contains the 6534-dimensional BR complement. {perf_counter()-t:.3f}s.")


if __name__ == "__main__":
    main()
