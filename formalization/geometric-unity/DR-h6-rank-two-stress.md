# DR. Stronger H6 stress obstruction at derivative rank two

Let \(J\) be a real H6 first derivative of rank at most two, where
rank is the rank of the map from tangent vectors to field derivatives.
If \(J\) satisfies the differentiated hooks and exterior closure,
the [DO](DO-h6-sharp-tableau-stress.md) derivative quadratic obeys
\[
 F(J)\ge-\tfrac53\|J\|^2 .
\]
This is sharp before imposing curvature. If every derivative
slice also satisfies coefficient-curvature, equality forces
\(J=0\) on every real derivative plane. Consequently there is
a uniform coefficient \(\beta_2<5/3\) valid for all rank-at-most-two
jets in the [DM](DM-h6-first-tableau.md) tableau.

At a rank-two point of a full H6 field, the [CJ](CJ-h6-universal-stress-gradient.md)
stress identity therefore gives
\[
 \|\nabla B\|^2>\tfrac35H(*B)>\tfrac{15579}{125000}\|B\|^2 .
\]
The last coefficient is \(0.124632\).
On the BS unit partial box the squared gradient must exceed
\(1.002583584\) at such a point if norm and scalarity persist
through second order. The uniform strict gap below \(5/3\) is
proved by compactness.

## DR1. The Euclidean two-plane spectrum

Fix an orthonormal derivative plane \(L=\operatorname{span}(e_0,e_1)\).
The differentiated hooks, exterior symbol and \(F\) are orthogonal
invariants before imposing curvature. Set \(J_i=0\) for derivative
indices outside \(L\). A coordinate \(J_{i,a;M}\) transforms under
diagonal signs by the mask \(\{i\}\mathbin\triangle\{a\}\mathbin\triangle M\).
Its size is four, six or eight. Permutation symmetry within the
derivative plane and its orthogonal complement reduces the
calculation to the number of mask bits inside and outside \(L\).

The verifier computes the exact kernel of the two hooks and closure
in each block and restricts \(F=2X-\|J\|^2+\|\delta J\|^2\) to it.
The six nonzero blocks have the following normalized determinants
\(\det(Q-tG)/\det G\).

| Inside bits | Outside bits | Normalized determinant |
| ---: | ---: | --- |
| 0 | 4 | \((t-1)^7(t+1)^7(3t-5)(3t+5)/9\) |
| 0 | 6 | \((t-1)^5(t+1)^5\) |
| 1 | 3 | \(t^8\) |
| 1 | 5 | \(t^4(t-1)(19t+5)/19\) |
| 2 | 4 | \((t+1)^7(3t-5)/3\) |
| 2 | 6 | \(-(t+1)^5\) |

Their copies sum to dimension 38,764. The minimum \(-5/3\)
is simple in each of the \(\binom{12}{4}=495\) masks of the first row.
Write \(\phi=e^0\wedge e^1\wedge\omega\) for
\(\omega\in\Lambda^4(L^\perp)^*\). The equality space consists of
\[
 J_{0,1}=\phi,\qquad J_{1,0}=-\phi,\qquad
 J_{i,a}=-\tfrac18 e^i\wedge e^a\wedge\omega
 \quad(i\in\{0,1\},\ a\perp L),
\]
with other components zero. The contraction hook follows from
\(\sum_{a\perp L}\iota_a(e^a\wedge\omega)=8\omega\).
The alternating hook and exterior symbol vanish directly. These
derivatives have squared norm \(9\|\omega\|^2/4\),
\(X=-2\|\omega\|^2\), and squared contraction norm
\(5\|\omega\|^2/2\). Their image has dimension 495, exhausting
the lowest eigenspace.

## DR2. Curvature removes all equality tensors

Identify coefficient covectors with vectors and use the solder
form \(P=\sum_a e_a\otimes e^a\). The first Bianchi identity gives
\(R\wedge P=0\). The slice with derivative index zero equals
\[
 J_0=e^0\wedge
 \left(\tfrac98 e_1\otimes(e^1\wedge\omega)
       -\tfrac18 P\wedge\omega\right).
\]
The other slice exchanges zero and one. Consequently
\(R\wedge J_0=(9/8)R(e_1)\wedge\phi\) and
\(R\wedge J_1=-(9/8)R(e_0)\wedge\phi\).
For \(K\in\mathfrak{so}(5)\), let
\(\omega_K(X,Y)=\operatorname{tr}(K[X,Y])\).
The five-by-five curvature normalization of BO gives
\(\langle R(X,Y)z,a\rangle=-\omega_{[z,a]}(X,Y)/4\).
Thus the two curvature equations require
\(\omega_K\wedge\phi=0\) for every
\(K\in[L,\mathfrak p]\), the trace-orthogonal complement of the
common compact centralizer of \(L\).

The common centralizer of two independent real trace-free symmetric
five-by-five matrices is contained, up to compact conjugation, in
a standard \(\mathfrak{so}(3)\) on three axes or a maximal rotation
torus \(\mathfrak{so}(2)\oplus\mathfrak{so}(2)\). To check this,
diagonalize one matrix. Its repeated eigenspaces give a product
of rotation algebras, already contained in a stated subalgebra except
for the partitions \(4+1\) and \(3+2\).

For \(4+1\), a nonzero cross block of the second matrix forces the
common algebra to fix a vector in the four-space, placing it in
\(\mathfrak{so}(3)\). With zero cross block, the symmetric four-block
is nonscalar by independence. Its repeated eigenvalues then give
the required containment. For \(3+2\) the common centralizer is
a proper subalgebra of \(\mathfrak{so}(3)\oplus\mathfrak{so}(2)\).
The full product would force the second matrix to be scalar in
each block and proportional to the first. The projection to
\(\mathfrak{so}(3)\) has dimension zero, one or three. The first
two cases lie in an abelian torus. In the last case the kernel
is zero by properness, and the graph homomorphism from the
perfect algebra \(\mathfrak{so}(3)\) to \(\mathfrak{so}(2)\)
vanishes.

For either stated standard subalgebra \(\mathfrak h_0\), the map
\[
 \Lambda^6\mathfrak p^*\longrightarrow
 \bigoplus_{K\ \text{basis of}\ \mathfrak h_0^\perp}\Lambda^8\mathfrak p^*,
 \qquad \phi\longmapsto(\omega_K\wedge\phi)_K
\]
has full column rank 3,003. The verifier constructs the seven
or eight moment forms from the exact ambient matrices. In each
of 16 diagonal sign blocks all columns are independent modulo
65,537, whose denominators are units. This proves the corresponding
characteristic-zero rank. The inclusion of centralizers puts the
conjugated \(\mathfrak h_0^\perp\) inside \([L,\mathfrak p]\).
Hence the curvature equations force \(\phi=0\) and \(J=0\).

## DR3. Uniformity and verification

The unit sphere of the DM tableau intersected with derivative
rank at most two is compact. The minimum of \(F\) is attained
and is greater than \(-5/3\). The real rank-two CU examples of
DO give \(\beta_2\ge1\). The additional tangent restrictions of the
nonlinear equations could strengthen this interval.

    python -X utf8 computation/geometric-unity/DR_h6_rank_two.py

The exact verifier checks all six ambient spectra, the equality
formula on a representative and both moment-form rank certificates.
