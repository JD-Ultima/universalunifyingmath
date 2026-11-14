# CA. Analytic closed lifts of all H4 curvature values

On the flat positive observer at $\lambda=1/5$, every real value in
the 10,388-dimensional [BX](BX-h4-curvature-and-closed-fields.md)
H4 dual-curvature kernel has a global real analytic closed extension.
The extension can be taken as a finite real sum of compact translates
of 28 explicit closed fields invariant under the positive upper
triangular subgroup. The claim concerns the linear derivative and
hook equations.

## CA1. Compact translates and evaluation

Use $Y=G/K$ with $G=SL(5,\mathbb R)$ and $K=SO(5)$. Write $S=AN$
for the positive upper triangular subgroup. As in BX, a closed
H4 value at $o=K$ for the left-invariant $S$ frame gives a global
analytic invariant field by left transport. In particular, the
28 rational certificate values satisfy both hook traces and
$\delta_\nabla C=0$, equivalently $d_\nabla(*C)=0$.

Let $W_0$ be the span of those values and $M$ the real linear
span of $K\cdot W_0$. The compact group fixes $o$, so the
evaluation of finite real sums of the corresponding translated
fields gives all of $M$. The hook and derivative equations are
equivariant under the isometries. Every such value obeys
dual curvature, hence $M$ is a submodule of the BX kernel.

## CA2. Highest-weight pairings

Complexify and use the invariant Hermitian product induced
by the positive physical metric. If $M$ were proper, its
orthogonal complement in the BX kernel would contain an
irreducible highest vector. The
[verifier](../../computation/geometric-unity/BX_data/verify_h4_lift.py)
checks the complete set of 57 rational extremals across the
22 highest-weight types in BX. In the split coordinates
of BO and BX it applies both hook traces, the full curvature
wedge of the dual form, and both simple root actions to
each witness. The additive covariant slot labels use the
reversed weight convention of BO.

The verifier transports each split extremal to the physical
coframe by the inverse tangent basis change on every covariant
slot. It uses the four-fold exterior minors on the form slots.
Within each weight type, pair the extremals with the selected
real $S$-invariant witnesses using the Hermitian product.
The certificate and verifier give a nonzero determinant for
all 22 resulting square matrices over $\mathbb F_{65537}$,
using $i\mapsto256$. The entries arise from Gaussian rational
numbers with denominators invertible at that prime. The
Gaussian rational determinants are therefore nonzero.
The second prime 65,521 independently reproduces nonzero
determinants in all types.

This verifies independence within each highest-weight space.
The exact multiplicities of BX show that the supplied sets
exhaust those spaces. Each such space intersects the
annihilator of $W_0$ trivially, by its invertible pairings.
It follows that the orthogonal complement of $M$ is zero.
Therefore $M$ is the full 10,388-dimensional kernel.

## CA3. Closed lifts and scope

A finite set of compact translates spans at $o$ because the
evaluation space is finite-dimensional. Choose a basis from
their values and linearly extend the corresponding fields.
This gives a real linear right inverse of evaluation into
global analytic closed H4 fields. Transitivity gives the same
result at any point. Positive scales transfer by the dilation
and homothety of BC, which preserve the linear equations.
The grade-one, grade-five and grade-nine quadratic equations
throughout an open set remain a separate condition.

Run the full verification with

    python -X utf8 computation/geometric-unity/BX_data/verify_h4_lift.py --exhaustive

[Review91](../../reports/reviews/REVIEW91-GU-CA-2026-11-14.md)
records separate exact constructions of the closure and pairing
checks. External peer review remains outstanding.
