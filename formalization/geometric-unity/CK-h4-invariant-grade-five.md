# CK. The full invariant grade-five form on H4

On the 10,388-dimensional [BX](BX-h4-curvature-and-closed-fields.md)
curvature kernel, the invariant grade-five form $L_5$ has exact
signature $(5893,4495,0)$. It gives explicit invariant subspaces of
dimensions 5,893 and 4,495 excluded from the full equation by its
positive and negative signs, respectively. A unit value satisfying
$L_5=0$ must carry positive energy in both orthogonal sign spaces.

## CK1. Exact multiplicity forms

Use the primitive grade-five matrix trace receiver of
[BQ](BQ-h6-invariant-moments.md), now on the actual real degree-four
Clifford input of [CB](CB-h4-certified-partial-value.md).
The [verifier](../../computation/geometric-unity/CK_h4_invariant_forms.py)
checks all 57 supplied compact highest vectors across the 22 BX
types. It verifies both hook traces, the full coefficient curvature
on the output of the form-slot Hodge star, both Cartan frequencies
and both raising equations. Positive exact Gram determinants
certify independence and the BX multiplicities give completeness.

The same compact element
$\operatorname{diag}(1,-1,1,-1,1)$ takes every highest vector to
its compact conjugate. Taking real parts thus uses one common
nonzero representation direction on each multiplicity space.
Schur orthogonality for the compact representation identifies the
computed norm and invariant quadratic matrices with that
multiplicity space's forms up to the same positive factor.

The helper rebuilds all 1,143 receiver coefficients for each of
129 diagonal and polarized evaluations. The rational six-word
Clifford arithmetic includes all external Gram and inverse
contraction factors. There are 147,447 exact calculations and
516 independent nested coefficient comparisons by default.
The full independent 128-dimensional matrix check reproduces all
147,447 coefficients modulo 65,521, using exact integer modular
arithmetic with explicit binary64 dot-product bounds.

## CK2. Exclusions and energy fractions

Exact $LDL^T$ congruence gives no zero diagonal in any multiplicity
matrix. Repeating its sign counts by the irreducible dimension gives
5,893 positive and 4,495 negative directions. The verifier gives
the rational congruence in each type. Its positive and negative
columns extend by compact transport to real invariant definite
subspaces of those dimensions. Neither contains a nonzero $L_5$
zero. The positive baseline prevents the zero value of the full
equation. These are maximal possible dimensions for a definite
subspace of this one form.

For sums of entire isotypic components, the largest positive
sum has dimension 2,097 and consists of
$$
 (6,3),\ (6,2),\ 2(5,1),\ 2(3,3),\ (2,0).
$$
The largest negative sum has dimension 853 and consists of
$(7,0),(5,4),(4,4),(4,0)$. All remaining isotypic components
are indefinite for $L_5$.

Let $A$ and $G$ be the multiplier and Gram matrices of a
given type. Exact inertia of $A+(17/3)G$, $6G-A$,
$A-G/100$ and $A+G/15$ puts the positive spectrum in
$(1/100,6)$ and the negative spectrum in $(-17/3,-1/15)$.
Write $E_+$ and $E_-$ for the squared norms in the orthogonal
positive and negative spectral subspaces of $L_5$, and
$N=E_++E_-$. At a nonzero $L_5$ zero the common magnitude
$q$ of its two sign contributions obeys
$$
 \frac1{100}E_+ <q<\frac{17}{3}E_-,
 \qquad
 \frac1{15}E_-<q<6E_+ .
$$
Thus $E_->(3/1703)N$ and $E_+>N/91$. These are estimates
in the orthogonal spectral splitting of the invariant form.

## CK3. The H4 second-order source on the invariant zero set

The curvature source $H$ of
[CH](CH-h4-second-order-stress.md) has both signs even on the
real unit set $L_5=0$. Exact witnesses use two distinct types
in each case, so invariant cross terms vanish. For normalized
vectors in type $(4,1)$ and $(3,0)$, squared amplitudes
$1963168/4792085$ and $2828917/4792085$ give
$H=-1847830827/1878497320$. For normalized vectors in types
$(5,4)$ and $(6,3)$, squared amplitudes $2/3$ and $1/3$
give $H=7/4$. Each pair has $L_5=0$ exactly. These are real
algebraic unit values obtained by taking positive square
roots of the amplitudes.

    python -X utf8 computation/geometric-unity/CK_h4_invariant_forms.py
    python -X utf8 computation/geometric-unity/CK_h4_invariant_forms.py --matrix-independent

[Review100](../../reports/reviews/REVIEW100-GU-CK-2026-11-14.md)
records the checks. The full grade-five and general field
equations remain open. External peer review remains outstanding.
