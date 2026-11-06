"""Exact complete SO4-fixed H6 coefficient-curvature kernel.

BG-h6-fixed-values.md supplies the completeness and nonlinear argument.
All tensor slots are covariant. BC and BF are reused as libraries without
running their mains. The adjacent certificate is reconstructed and checked;
running this helper writes stdout only.
"""
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path
from time import perf_counter
import json
import sys
import sympy as s

PUBLIC = Path(__file__).resolve().parent
sys.path.insert(0, str(PUBLIC))
import BC_affine_hooks as BC
import BF_h6_casimir as BF

N = 14
PAIRS = tuple(combinations(range(N), 2))


def mm(a, b):
    return tuple(tuple(sum(a[i][k]*b[k][j] for k in range(5))
                       for j in range(5)) for i in range(5))


def bracket(a, b):
    ab, ba = mm(a, b), mm(b, a)
    return tuple(tuple(ab[i][j]-ba[i][j] for j in range(5)) for i in range(5))


def trace_product(a, b):
    return sum(a[i][j]*b[j][i] for i in range(5) for j in range(5))


def integer_phi():
    """Return Q_a=5 Phi_a in the BC Gram-orthogonal tangent basis."""
    sym = []
    for i, j in combinations(range(4), 2):
        mat = [[0]*4 for _ in range(4)]
        mat[i][j] = mat[j][i] = 1
        sym.append(mat)
    for k in range(1, 4):
        diagonal = [1]*k+[-k]+[0]*(3-k)
        sym.append([[diagonal[i] if i == j else 0 for j in range(4)] for i in range(4)])
    sym.append([[int(i == j) for j in range(4)] for i in range(4)])
    phi = []
    for a in range(N):
        mat = [[0]*5 for _ in range(5)]
        if a < 4:
            mat[a][4] = mat[4][a] = -5
        else:
            S = sym[a-4]
            tr = sum(S[i][i] for i in range(4))
            for i in range(4):
                for j in range(4):
                    mat[i][j] = -10*S[i][j]+(2*tr if i == j else 0)
            mat[4][4] = 2*tr
        phi.append(tuple(map(tuple, mat)))
    assert all(sum(mat[i][i] for i in range(5)) == 0 and
               mat == tuple(zip(*mat)) for mat in phi)
    assert all(trace_product(phi[a], phi[b]) ==
               (25*BC.METRIC[a] if a == b else 0)
               for a in range(N) for b in range(N))
    # Exact intertwining with all six BC stabilizer generators.
    for which, (i, j) in enumerate(combinations(range(4), 2)):
        J = [[0]*5 for _ in range(5)]
        J[i][j], J[j][i] = 1, -1
        for a in range(N):
            rhs = bracket(J, phi[a])
            for row in range(5):
                for col in range(5):
                    lhs = sum(-dict(BC.COV[which][d]).get(a, 0)*phi[d][row][col]
                              for d in range(N))
                    assert lhs == rhs[row][col]
    return phi


def curvature_model():
    """R_ij[d,a] is the vector curvature coefficient, full ambient model.

    R(Phi_i,Phi_j)Phi_a=-1/4 [[Phi_i,Phi_j],Phi_a].
    With Q=5Phi the coefficient is -tr([[Qi,Qj],Qa]Qd)/(2500*g_d).
    """
    phi = integer_phi()
    R = {}
    for i, j in PAIRS:
        K = bracket(phi[i], phi[j])
        columns = [bracket(K, phi[a]) for a in range(N)]
        R[i, j] = tuple(tuple(-F(trace_product(columns[a], phi[d]),
                                2500*BC.METRIC[d])
                            for a in range(N)) for d in range(N))

    def value(i, j, d, a):
        if i == j:
            return F(0)
        return R[i, j][d][a] if i < j else -R[j, i][d][a]

    assert all(BC.METRIC[d]*R[i, j][d][a] + BC.METRIC[a]*R[i, j][a][d] == 0
               for i, j in PAIRS for d in range(N) for a in range(N))
    assert all(value(i, j, d, k)+value(j, k, d, i)+value(k, i, d, j) == 0
               for i, j, k in combinations(range(N), 3) for d in range(N))
    assert all(sum(value(i, j, i, k) for i in range(N)) ==
               (-F(5, 4)*BC.METRIC[j] if j == k else 0)
               for j in range(N) for k in range(N))
    print('Full ambient Phi curvature: metric skewness, Bianchi and Ric=-5g/4 verified.')
    print('All six SO4 tangent actions intertwine Phi; all HH, HV and VV channels retained.')
    return R


def output_coefficient(R, tensor, output):
    """One actual coefficient of F^(T*) wedge B; coefficient covaction -R^T."""
    a, mask = output
    indices = BC.bit_indices(mask)
    assert len(indices) == 8
    result = F(0)
    for u, v in combinations(range(8), 2):
        i, j = indices[u], indices[v]
        reduced = mask ^ (1 << i) ^ (1 << j)
        sign = (-1)**(u+v-1)
        result -= sign*sum(R[i, j][d][a]*tensor.get((d, reduced), 0)
                           for d in range(N))
    return result


def full_sparse_curvature(R, tensor):
    """All slots, used only to validate the newly found three nullvectors."""
    out = {}
    for (d, mask), value in tensor.items():
        for i, j in PAIRS:
            pairmask = (1 << i) | (1 << j)
            if mask & pairmask:
                continue
            sign = (-1)**((mask & ((1 << i)-1)).bit_count()
                          + (mask & ((1 << j)-1)).bit_count())
            for a, curvature in enumerate(R[i, j][d]):
                if curvature:
                    BC.plus(out, (a, mask | pairmask), -sign*curvature*value)
    return out


def matrix_at(basis, rows):
    return s.Matrix([[s.Rational(tensor.get(key, 0)) for tensor in basis] for key in rows])


def dense_coefficients(vector, dimension):
    return s.Matrix([s.Rational(vector.get((j, 0), 0)) for j in range(dimension)])


def main():
    start = perf_counter()
    print('BG: flat positive Z1, lambda=1/5, cH=2, cV=1; SO4-fixed H6 point values.')
    print('Gram:', tuple(map(str, BC.METRIC)))
    invariant = BC.invariant_basis(6)
    hooks = BC.hook_basis(invariant)
    assert (len(invariant), len(hooks)) == (30, 23)
    target = [BC.hodge_scaled(tensor) for tensor in invariant]
    assert all(not BC.action(j, tensor) for j in range(6) for tensor in target)
    trank, trows, tcols, tdet, tkernel = BC.rank_witness(target)
    assert trank == 30 and len(trows) == 30 and not tkernel
    T = matrix_at(target, trows)
    assert T.det() == s.Rational(BC.fraction(tdet)) != 0
    print('Complete target: Hodge(raw degree-six invariants), dimension', trank)
    print('Target coordinate-isomorphism rows (external, eight-form mask):', trows)
    print('Target coordinate-isomorphism determinant:', tdet)

    R = curvature_model()
    images = [{key: output_coefficient(R, tensor, key) for key in trows}
              for tensor in hooks]
    images = [{key: v for key, v in image.items() if v} for image in images]
    rank, rows, cols, determinant, kernel = BC.rank_witness(images)
    M = matrix_at(images, trows)
    assert M.shape == (30, 23) and not M.atoms(s.Float)
    assert M.rank() == rank == 20 and len(kernel) == 23-rank == 3
    assert all(M*dense_coefficients(vector, 23) == s.zeros(30, 1) for vector in kernel)
    print('Full invariant coefficient-curvature rank:', rank)
    print('Full invariant coefficient-curvature nullity:', len(kernel))
    print('Nonzero curvature minor rows:', rows)
    print('Nonzero curvature minor columns:', cols)
    print('Nonzero curvature minor determinant:', determinant)
    print('Curvature kernel amplitude vectors in reconstructed hook basis:',
          [tuple((j, str(c)) for (j, _), c in sorted(vector.items())) for vector in kernel])
    kernel_tensors = [BC.combine(vector, hooks) for vector in kernel]
    for tensor in kernel_tensors:
        assert BC.hook_traces(tensor) == ({}, {})
        assert all(not BC.action(j, tensor) for j in range(6))
        assert not full_sparse_curvature(R, tensor)
    print('Every slot of the full sparse curvature output vanishes on all', len(kernel), 'nullvectors.')

    # Basis coordinate reconstruction is separately certified, not inferred
    # from selected coefficients without checking the full sparse tensor.
    srank, srows, scols, sdet, skernel = BC.rank_witness(hooks)
    assert srank == 23 and not skernel
    S = matrix_at(hooks, srows)
    assert S.det() == s.Rational(BC.fraction(sdet)) != 0
    Sinv = S.inv()
    data = json.loads((PUBLIC/'BC_h6_closed_generator.json').read_text(encoding='utf-8'))
    assert data['degree'] == 6 and tuple(map(F, data['metric_gram'])) == tuple(BC.METRIC)
    B = {(a, mask): F(value) for a, mask, value in data['entries']}
    CB = BF.casimir(B)
    known = []
    for name, tensor in (('Bstar', B), ('Casimir Bstar', CB)):
        assert BC.hook_traces(tensor) == ({}, {})
        assert all(not BC.action(j, tensor) for j in range(6))
        coordinates = Sinv*s.Matrix([s.Rational(tensor.get(key, 0)) for key in srows])
        vector = {(j, 0): BC.fraction(value) for j, value in enumerate(coordinates) if value}
        assert BC.combine(vector, hooks) == tensor
        assert M*coordinates == s.zeros(30, 1)
        known.append(coordinates)
        print(name+' hook amplitudes:', tuple(map(str, coordinates)))
    K = s.Matrix.hstack(*known)
    assert K.rank() == 2
    KN = s.Matrix.hstack(*(dense_coefficients(vector, 23) for vector in kernel))
    assert s.Matrix.hstack(K, KN).rank() == len(kernel)
    known_change, known_parameters = KN.gauss_jordan_solve(K)
    assert known_parameters.rows == 0 and KN*known_change == K
    print('Bstar and Casimir Bstar in the complete curvature-kernel basis:', known_change.tolist())
    print('Known span dimension:', K.rank())
    assert len(kernel)-K.rank() == 1
    print('Additional invariant curvature directions beyond the known two:', len(kernel)-K.rank())

    # The full geometric curvature map is SO5-equivariant.  Its SO4-fixed
    # kernel is therefore stable under the positive SO5 Casimir, even though
    # this fixed subspace is not itself an SO5 representation.
    casimir_columns = []
    for tensor in kernel_tensors:
        image = BF.casimir(tensor)
        assert BC.hook_traces(image) == ({}, {})
        assert all(not BC.action(j, image) for j in range(6))
        coordinates = Sinv*s.Matrix([s.Rational(image.get(key, 0)) for key in srows])
        vector = {(j, 0): BC.fraction(value) for j, value in enumerate(coordinates) if value}
        assert BC.combine(vector, hooks) == image
        assert M*coordinates == s.zeros(30, 1)
        column, parameters = KN.gauss_jordan_solve(coordinates)
        assert parameters.rows == 0 and KN*column == coordinates
        casimir_columns.append(column)
    casimir = s.Matrix.hstack(*casimir_columns)
    z = s.Symbol('z')
    polynomial = s.factor(casimir.charpoly(z).as_expr())
    print('Positive SO5 Casimir matrix in the complete invariant curvature kernel:', casimir.tolist())
    assert s.expand(polynomial-(z-10)*(z-18)*(z-40)) == 0
    print('Casimir characteristic polynomial:', polynomial)
    spectra = []
    for eigenvalue, multiplicity, eigenspace in casimir.eigenvects():
        eigenvectors = []
        for column in eigenspace:
            vector = {(j, 0): BC.fraction(value) for j, value in enumerate(column) if value}
            tensor = BC.combine(vector, kernel_tensors)
            assert BF.casimir(tensor) == {key: eigenvalue*value for key, value in tensor.items()}
            assert not full_sparse_curvature(R, tensor)
            eigenvectors.append(tuple(map(str, column)))
        print('Casimir eigenvalue', eigenvalue, 'multiplicity', multiplicity,
              'kernel-amplitude eigenspace:', eigenvectors)
        spectra.append({'eigenvalue': str(eigenvalue), 'multiplicity': multiplicity,
                        'eigenvectors': eigenvectors})
    assert sum(multiplicity for _, multiplicity, _ in casimir.eigenvects()) == len(kernel)
    receipt = {
        'scope': 'Full SO4-fixed H6 coefficient-curvature kernel at flat positive Z1 lambda=1/5,cH=2,cV=1; point values only.',
        'degree': 6,
        'basis_order': data['basis_order'],
        'metric_gram': list(map(str, BC.METRIC)),
        'rank': rank,
        'kernel_dimension': len(kernel),
        'curvature_convention': 'R_ij Z=-1/4 [[Phi_i,Phi_j],Z]; coefficient curvature acts by -R^T.',
        'kernel_tensors': [
            {'name': 'K'+str(j),
             'hook_amplitudes': [(i, str(v)) for (i, _), v in sorted(vector.items())],
             'entries': [(a, mask, str(v)) for (a, mask), v in sorted(tensor.items())]}
            for j, (vector, tensor) in enumerate(zip(kernel, kernel_tensors))],
        'known_Bstar_CasimirBstar_coordinates': [[str(v) for v in row] for row in known_change.tolist()],
        'casimir_matrix': [[str(v) for v in row] for row in casimir.tolist()],
        'casimir_characteristic_polynomial': str(polynomial),
        'casimir_spectra': spectra,
    }
    recorded = json.loads((PUBLIC/'BG_h6_curvature_kernel.json').read_text(encoding='utf-8'))
    assert json.loads(json.dumps(receipt)) == recorded
    print('Complete three-tensor and Casimir certificate reconstructed and matched.')
    print('Receipt tensor term counts:', tuple(map(len, kernel_tensors)))
    print('Source coordinate-isomorphism rows:', srows)
    print('Source coordinate-isomorphism determinant:', sdet)
    print('Boundary: invariant point VALUES only; no closure restriction, no full SO5 kernel or GU solution claim.')
    print(f'Runtime {perf_counter()-start:.3f} seconds.')


if __name__ == '__main__':
    main()
