# CX. External rank and sparse support in H4 and H6

For a nonzero real hook in a positive 14-dimensional orthonormal
frame, scalar grade one requires external rank at least four in
H4 and six in H6. An independent sparse-support criterion excludes
every forest support and certain cyclic supports even for complex
hooks. Both results are pointwise and hold before imposing
coefficient curvature.

## CX1. Minimum external rank

Use \(p=4\) for the H4 actual input, or \(p=6\) for the H6
precursor to multiplication by \(\Gamma\). Write
\(E_{ab}=\langle B_a,B_b\rangle\),
\(I_{ab}=\sum_c\langle\iota_aB_c,\iota_bB_c\rangle\) and
\(N=\operatorname{tr}E\). The full scalarity equations are
$$
 \mathrm{H4}:\ I=2E+Ng/7,\qquad
 \mathrm{H6}:\ I+2E=4Ng/7 .
$$
These are the identities of [CB](CB-h4-certified-partial-value.md)
and [BS](BS-h6-certified-partial-value.md).

Let \(R\) be the external range of \(E\), \(r=\dim R\), and
decompose \(B\) by the number \(k\) of its form slots in \(R\).
All external components in \(R^\perp\) vanish. The alternating
hook map is injective when \(k=0\), and the contraction hook map
on \(R^*\otimes\Lambda^r R^*\) is an isomorphism when \(k=r\).
Each map preserves the \(R^\perp\) form degree, so those two
components vanish. Let \(N_k\) be the squared norm in the
remaining components. Then
$$
 N\le\operatorname{tr}_R I=\sum_{k=1}^{r-1}kN_k\le(r-1)N .
$$
For H4 the middle trace is \((2+r/7)N\), giving \(r\ge4\).
For H6 it is \((4r/7-2)N\), giving \(r\ge6\).

At the H4 minimum, \(w_k=N_k/N\) satisfies
\(w_3-w_1=4/7\), with \(0\le w_1\le3/14\)
and \(4/7\le w_3\le11/14\). At the H6 minimum,
\(\sum_{k=2}^5(k-1)w_k=3/7\), so \(w_1\ge4/7\)
and \(\sum_{k\ge j}w_k\le3/[7(j-1)]\) for \(2\le j\le5\).

The H4 rank bound is attained with scalarity and the complete
grade-nine equation before curvature is imposed. The exact
example \(B=B_1+\sqrt{22}B_3\) uses one and three form slots
in its four-dimensional external space. It has
\(E=84I_4\oplus0_{10}\) and \(N=336\).
Every displayed form in row \(a\) contains \(e^a\), so its
exterior square vanishes. The standalone construction also
checks the full grade-nine equation using the nested Clifford
formula. Its grade-five output is nonzero.

## CX2. Forest and cycle support criteria

Fix a support of \(m\) distinct covariant slots \((a,A)\),
where \(A\) is a \(p\)-blade. Form a bipartite graph with
external indices as left vertices. A slot has right vertex
of contraction kind with label \(A\setminus\{a\}\) if
\(a\in A\). Otherwise it has right vertex of alternating kind
with label \(A\cup\{a\}\). Each slot is one edge with variable
\(x_e\). A right-vertex equation is a signed hook sum.

Take the homogeneous integer ideal comprising these hook sums
and the cleared scalarity equations. Reduce it modulo two.
Squaring every hook sum puts the sum of edge squares at its
right vertex in the ideal and in particular \(N\) in the ideal.
For a right vertex labeled by \(S\), its incident form blade
is \(S\mathbin{\triangle}\{a\}\). Thus modulo the right-vertex
square sums \(I_{ii}=E_{ii}\) for each \(i\). The diagonal
scalarity equations now put the sum of squares at every left
vertex in the ideal as well.

Peel all tree edges by successively using a degree-one vertex.
The square of each peeled edge lies in the reduced ideal. For
a forest this gives all squares. Every degree \(m+1\) monomial
contains a square. The degree \(m+1\) Macaulay matrix therefore
has a nonzero maximal minor modulo two, and the same minor
is nonzero over the rationals. Full rank at that degree
forces the common complex zero set to be the origin.

Suppose the remaining core consists of vertex-disjoint simple
cycles. Define a row for every unordered pair of distinct
physical axes \(i,j\), and a column for every cycle.
Its entry is the parity of the external rows in that cycle
whose two cycle blades differ by exchanging \(i,j\).
Vertex equations identify all edge squares in each cycle with
one element \(t_C\). The off-diagonal scalarity equations
give \(I_{ij}=0\) modulo two, and their squares give linear
equations on the \(t_C^2\) with precisely this parity matrix.
If the matrix has full column rank, every edge fourth
power is in the reduced ideal. The same degree argument at
\(3m+1\) excludes every nonzero complex value. For a single
cycle one odd exchange parity suffices. The proof lifts the
full-degree rank, which is the quantity preserved by the
nonzero minor.

## CX3. Exact examples

Choose disjoint axes \(u_1,\ldots,u_k\), external axes
\(a_1,\ldots,a_k\), and a common \((p-2)\)-blade \(K\)
disjoint from both. In row \(a_i\) use the blades
\(K\cup\{u_i,a_i\}\) and \(K\cup\{u_{i+1},a_i\}\)
with cyclic indexing. The graph has one cycle of length
\(2k\), with odd parity for each consecutive pair of
\(u\)-axes when \(k\ge3\). This fits \(k=3,\ldots,6\)
in H4 and \(k=3,\ldots,5\) in H6.

    python -X utf8 computation/geometric-unity/CX_sparse_hook_support.py
    python -X utf8 computation/geometric-unity/CX_h4_minimum_external_rank.py

The verifier checks the actual Clifford grade-one coefficients
and exact Macaulay ranks on all seven families. It also checks
12 graph controls, including forest and mixed-cycle cases.
These criteria do not decide all cyclic supports. The
positive baseline excludes every full field value on the
stated supports.

[Review108](../../reports/reviews/REVIEW108-GU-CW-CX-2026-11-14.md)
records the checks. External peer review remains outstanding.
