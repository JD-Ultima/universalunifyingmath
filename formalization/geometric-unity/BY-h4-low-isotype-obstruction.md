# BY. Grade-five exclusion on a 24-dimensional H4 subspace

On [BX](BX-h4-curvature-and-closed-fields.md)'s pure H4 branch,
consider the two copies of type $(1,0)$ and one copy of type $(2,0)$
in the curvature kernel. Their complete real isotypic sum has
dimension 24. Its grade-five quadratic equation has only the zero
value, even after complexification. A full real field cannot take a
value in this subspace at any point because the zero value fails the
positive grade-one baseline.

## BY1. The complete sum and its receivers

Use the physical frame of [BC](BC-affine-hooks.md) at
$(c_H,c_V)=(2,1)$ and the actual degree-four Clifford input $C$.
The [certificate](../../computation/geometric-unity/BY_data/h4_low_witness.json.gz)
contains 24 integer covariant tensors $C_0,\dots,C_{23}$.
The [verifier](../../computation/geometric-unity/BY_data/verify_h4_low.py)
checks both hook traces and the full coefficient-curvature wedge on
the dual ten-form. A modular minor proves that the tensors are
independent.

The Casimir is the negative sum of squared actions of the ten
unit coordinate plane rotations on all five slots. It equals
$a(a+3)+b(b+1)$ on weight type $(a,b)$ in BX's convention.
Ten of the certificate tensors have eigenvalue 4 and fourteen
have eigenvalue 10, verified exactly in the Gram-orthogonal frame.
The BX character table shows that no other type in the kernel has
either eigenvalue. The dimension and independence give the entire sum.

Write $C=\sum_i x_i C_i$ and $z_{ij}=x_ix_j$ for $i\le j$.
The certificate specifies 297 coordinates of the true grade-five
quadratic in the external and exterior slots. Applying the complete
BC differentiated cubic gives their coefficients as a $297\times300$
matrix $A$ on the monomials $z_{ij}$. The ordered products use
Clifford generator squares $g_i^{-1}$, external factor $g_a$, and
the mixed coefficients for both input orders. The verifier compares
selected diagonal and polarized values to native rational products.

## BY2. Rank-one identities

A pivot minor of $A$ is nonzero modulo 65,537. Its size is 297,
so the rational nullspace has dimension three. Solve in the three
free coordinates using that pivot minor. Its determinant is a unit
in the localization at 65,537. Thus the reduced modular null basis
with identity on the free coordinates is the reduction of the
corresponding rational basis.

Extend each $z_{ij}$ to a symmetric $24\times24$ matrix. The
selected six principal minors of order two are homogeneous
quadratics in the three free coordinates on the nullspace.
The verifier computes their $6\times6$ coefficient matrix using
that modular null basis. The matrix has rank six modulo 65,537.
The defined rational coefficient matrix therefore has rank six.
An independent repetition over 65,521 gives the same two ranks.

An actual tensor has $z=xx^T$, so all these minors vanish.
Since they span the six-dimensional space of quadratics in the
free coordinates, the square of each free coordinate vanishes.
This forces $z=0$ and then $x=0$ over both $\mathbb R$ and
$\mathbb C$. This argument uses the rank-one identities in
addition to the linear receiver rank.

On the declared branch the derivative term has grade three,
so it does not enter this calculation. At $C=0$, [BD](BD-affine-h4.md)
gives baseline $15/(2c_V)+7\kappa^2/12>0$. Constant positive
scales transfer by the same dilation and homothety.

## Reproduction and scope

    python -X utf8 computation/geometric-unity/BY_data/verify_h4_low.py --exhaustive

The [internal review](../../reports/reviews/REVIEW89-GU-BY-2026-11-14.md)
records an independent reconstruction of the selected coefficients.
The restriction to these isotypic components is essential to
the stated proof. Full pointwise feasibility outside this sum
and general field solutions remain open. External peer review
remains outstanding.
