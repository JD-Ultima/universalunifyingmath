# BT. Full quadratic exclusion on a 73-dimensional closed H6 space

The upper triangular closure space in
[BR](BR-h6-closed-field-lifts.md) splits into sixteen sign-character
blocks. The two 73-dimensional blocks with character masks 15 and 30
are each excluded from the full pointwise equation on the stated branch.
Each complexified block has only the zero solution of the homogeneous
scalarity and higher-grade equations together.

The [verifier](../../computation/geometric-unity/BT_h6_an73_quadratic.py)
reconstructs the exact closed basis and checks the modular coefficient
map. Use `--full` to recompute the complete rank certificate.

## BT1. The homogeneous map

Let $b_1,\ldots,b_{73}$ be a rational basis of the exact closed
character block. Each field from this space extends by the simply
transitive triangular action and has constant norm. The source for
the following calculation is a point value $B=\sum x_i b_i$.
The grade-one equation requires the scalarity residual
$$
q_1(B)+(26/21)\|B\|^2g
$$
to vanish, as in [BS](BS-h6-certified-partial-value.md). Use
this homogeneous tracefree expression together with all the
grade-five, grade-nine, and grade-thirteen rows. The unadjusted
grade-one trace and the unit-norm equation do not enter the map.

The sign group acts on this input space by one character, whose
square is one. Thus its quadratic image is invariant, including
on the complexification. In the physical coframe, the only
possibly nonzero entries have equal external and form-blade
characters. There are 3,488 such entries in grades 1, 5, 9, 13.
The rank calculation retains structurally zero entries as controls;
the sum of the corrected grade-one diagonal rows is identically
zero. There are 2,701 distinct quadratic monomials in 73 variables.

## BT2. Exact coefficient map and modular rank

Use the 128-dimensional Clifford matrices and the complete
BS1 quadratic. The verifier reduces the Gram square roots and
$i$ in the field of 65,521 elements, after checking that they
exist. For each basis vector it evaluates the quadratic, and for
each distinct pair it subtracts the two pure squares from the
sum's value. This yields the matrix on all 2,701 monomials.
The field reduction respects the coordinate expression because
the required Gram factors are units. The square roots implement
an invertible change to an orthonormal coframe. All original
hook and closure data are rational.

Each dense product used for field arithmetic has exact integer
operands in the available binary floating representation. The
bound on the absolute sum of 128 products is below $2^{53}$,
so every such sum is an exactly representable integer. Operands
are reduced between products. A separate integer multiplication
checks the implementation. The verifier also checks exact rational
nested-Clifford receivers in all four grades. Modular elimination
uses exact field inverses and Schur complements on successively
selected nonzero minors.

The complete coefficient matrix has column rank 2,701 modulo
65,521, certifying full column rank in characteristic zero.
Consequently every monomial $x_i x_j$ must vanish at a common
zero, including each $x_i^2$. Over the complex numbers this
forces all coordinates to vanish. The positive degree-one
baseline excludes the zero real field value on the declared
branch. Thus no point in this block can satisfy the full equation.

## BT3. The reversed character

Let $P$ reverse the five defining axes and let $U$ be its induced
tangent map. Conjugation by $P$ interchanges strict upper and lower
parts. BR1 therefore gives
$\Lambda_{UX}=-U\Lambda_X U^{-1}$. The helper checks the complete
connection identity, positive metric and orientation. On all seven
slots, this intertwines the closure map with its negative. The two
hook traces are natural and the homogeneous Clifford map is
equivariant. Five-bit reversal sends mask 15 to mask 30, so it
identifies the two closed blocks and transfers the quadratic exclusion.

Mixtures of the two blocks or of other triangular characters
require the mixed equations. External peer review remains outstanding.
