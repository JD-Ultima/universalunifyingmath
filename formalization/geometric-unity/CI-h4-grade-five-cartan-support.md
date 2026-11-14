# CI. Grade five excludes a 512-dimensional H4 support

On the pure real H4 branch of [BX](BX-h4-curvature-and-closed-fields.md),
the full grade-five quadratic has only the zero value on the
complete isotypic sum
$$
 3V_{(3,1)}\oplus3V_{(2,2)}\oplus4V_{(2,1)}
 \oplus V_{(2,0)}\oplus2V_{(1,0)} .
$$
The sum has dimension 512. The obstruction holds over the complex
numbers, and the positive grade-one baseline excludes the forced
zero on the real branch.

## CI1. Highest-weight quadratic certificates

The [verifier](../../computation/geometric-unity/CI_h4_cartan_support.py)
uses ten complex highest vectors stored as exact real and imaginary
components in the physical coframe of BC. They have multiplicities
three, three and four as in BX. Each vector passes both hook traces,
the full curvature map on the Hodge dual, the two Cartan frequency
equations and the two raising equations. Each block has invertible
Hermitian Gram matrix. Thus the vectors span the corresponding
highest spaces.

In a block of multiplicity $m$, let $h(x)$ be its highest vector with
coefficients $x_1,\ldots,x_m$. Evaluate actual external and five-form
coordinates of $q_5(h(x))$ from the six ordered Clifford terms.
The exact coefficients include the mixed terms and all physical
Gram factors. The independent complete nested engine agrees on all
136 selected coefficients.

For $(2,1)$ the ten selected quadratics span every quadratic in the
four variables. Their coefficient matrix has an exact nonzero
determinant. For each of $(3,1)$ and $(2,2)$, take the three selected
quadratics and multiply each by all six quadratic monomials.
The resulting 18 quartics span the 15 quartic monomials by a
nonzero exact integer minor found modulo 65,521. In either case
the ideal contains a power of each generator $x_i$. Thus the
necessary grade-five equations force $x=0$ for those highest values.

## CI2. Descending weight elimination

With simple roots $(1,-1),(0,1)$, the highest weights are ordered by
the first coordinate and the sum of the two coordinates. The chain is
$$
 (3,1)>(2,2)>(2,1)>(2,0)>(1,0).
$$
In particular, every weight in a lower module lies strictly below
the highest weight $\lambda$ in the maximal remaining module.
The output weight $2\lambda$ of the quadratic can only come from
the product of the weight-$\lambda$ input components. Consequently
its vanishing forces the highest projection in that isotypic space
to vanish by CI1.

Full grade-five vanishing is equivariant under the complexified
group. Apply the same argument to every group translate. The
orbit of a highest coordinate functional spans the dual of the
irreducible representation, so the maximal isotypic component
vanishes in every multiplicity copy. Remove the first three
types in succession. The remaining $(2,0)$ and $(1,0)$ sum is
exactly the 24-dimensional space excluded in
[BY](BY-h4-low-isotype-obstruction.md). The multiplicity argument
does not require any choice of one distinguished copy.

## Reproduction and scope

    python -X utf8 computation/geometric-unity/CI_h4_cartan_support.py --full-independent

The verifier checks the three exact minors and only uses the
selected necessary receiver coefficients for CI1.
[Review98](../../reports/reviews/REVIEW98-GU-CI-2026-11-14.md) records
the verification. The statement is confined to the displayed
isotypic sum. External peer review remains outstanding.
