# BX. H4 curvature and triangular closure dimensions

On the flat positive Z1 observer at $\lambda=1/5$, the complete
real $H_4$ hook space has dimension 11,648. Its values satisfying
the necessary coefficient-curvature condition have dimension 10,388
and decompose into 57 irreducible copies across 22 $SO(5)$ types.
The positive upper triangular subgroup gives a 6,643-dimensional
space of invariant fields that solve the linear derivative equation
everywhere. General closed fields and the nonlinear field equation
require further analysis.

Use [BD](BD-affine-h4.md)'s compact real branch with $C$ the actual
degree-four hook and $T_o=i\kappa\Phi^1/24$. By
[BA](BA-hook-closure.md), the derivative equation is
$\delta_\nabla C=0$, equivalent to $d_\nabla(*C)=0$ with Hodge
star on the form slots. A necessary condition on the value is
$R^{T^*Y}\wedge(*C)=0$. All three linear equations below use the
trace-metric model $V=\operatorname{Sym}^2_0(\mathbb R^5)$.

## BX1. The pointwise curvature kernel

On $V^*\otimes\Lambda^4 V^*$, of dimension
$14\binom{14}{4}=14{,}014$, stack the alternating trace,
the contraction trace, and $R^{T^*Y}\wedge(*C)$. The first two
targets are $\Lambda^5 V^*$ and $\Lambda^3 V^*$. They give
independent summands by the orthogonal hook decomposition in BA.
Their ranks sum to 2,366. The curvature output is in
$V^*\otimes\Lambda^{12}V^*$. It uses the coefficient-covector
action of $R_{ij}P=-[[P_i,P_j],P]/4$.

Use [BO](BO-h6-full-curvature-kernel.md)'s rational split Cartan
basis, in which the five-dimensional symmetric Gram form has
entries $K_{01}=K_{10}=K_{23}=K_{32}=2$ and $K_{44}=1$.
The tangent trace metric has its full inverse in this basis.
The Hodge map uses the ordered complementary blade with each
corresponding inverse-metric minor. The common volume factor
is omitted because it is nonzero. Opposite Cartan weights have
equal ranks by the Weyl group. Complexification preserves the
dimension of a real kernel.

The [verifier](../../computation/geometric-unity/BX_data/verify_h4_dual_curvature.py)
builds the zero weight and 74 representatives of opposite pairs.
It clears denominators in each rational matrix and computes rank
over each of $\mathbb F_{65537}$ and $\mathbb F_{65521}$.
A nonzero modular minor determines pivot rows and columns.
Writing the full matrix as $A$ and its pivot minor as $D$, it
solves $DX=A_{\text{pivot rows},:}$ exactly and checks the
entire identity $A_{:,\text{pivot columns}}X=A$ over $\mathbb Q$.
This supplies the matching characteristic-zero upper bound.
The resulting full stacked rank is 3,626, so the restriction of
the curvature operator to the hook space has rank 1,260.

The [character check](../../computation/geometric-unity/BX_data/verify_h4_character.py)
takes each weight multiplicity from these exact nullities and
subtracts irreducible characters in descending highest-weight order.
It uses Kostant's multiplicity formula for the positive $B_2$
roots $e_1-e_2,e_2,e_1,e_1+e_2$ and checks every remaining
multiplicity after each subtraction. The resulting
[type table](../../computation/geometric-unity/BX_data/h4_character.json)
contains 57 copies of 22 types and sums to 10,388.

## BX2. The triangular closed fields

Let $S=AN$ be the positive upper triangular subgroup of $SL(5,\mathbb R)$.
It acts simply transitively on $Y=SL(5,\mathbb R)/SO(5)$ as in
[BR](BR-h6-closed-field-lifts.md). In the left-invariant coframe,
an invariant field is determined by its value at the identity.
Use the physical Gram-orthogonal matrices $\Phi_i$ of
[BC](BC-affine-hooks.md). Set $Q_i=5\Phi_i$ and let $K_i$ be the
strict upper triangular part of $Q_i$ minus its transpose.
The Nomizu connection matrix on vectors is
$$
 (L_i)_{da}=\frac{\operatorname{tr}([K_i,Q_a]Q_d)}{250\,g_d},
 \qquad g_d=\operatorname{tr}(\Phi_d^2).
$$
Act by the negative transpose on each covariant slot. The
codifferential is $-\sum_i g_i^{-1}\iota_{\Phi_i}\nabla_i C$,
with interior product on the four form slots. This differentiates
all five slots. The subgroup's invariance carries closure at
the identity to every point of $Y$.

Stack this codifferential with the two hook traces and decompose
under the sixteen determinant-one diagonal sign characters. An
off-diagonal coframe label is the pair of its matrix indices.
An input coefficient has the exclusive-or of the labels of its
five covectors. The diagonal coframe labels are trivial.

| Subset size | Blocks | Columns per block | Rank per block | Nullity per block |
| ---: | ---: | ---: | ---: | ---: |
| 0 | 1 | 814 | 436 | 378 |
| 2 | 10 | 859 | 456 | 403 |
| 4 | 5 | 922 | 475 | 447 |

The [verifier](../../computation/geometric-unity/BX_data/verify_h4_an.py)
rebuilds all sixteen integer matrices and checks the two modular
ranks. All 6,643 archived integer null vectors pass every row.
Within each block each vector has its own nonzero free coordinate
and is zero in the other designated free coordinates. This proves
independence and the matching exact rank.

## Reproduction and scope

Run the complete rank and character checks with

    python -X utf8 computation/geometric-unity/BX_h4_kernels.py

The independently assembled Hodge-closure and double-contraction
checks are recorded in [Review88](../../reports/reviews/REVIEW88-GU-BX-2026-11-14.md).
Constant positive observer scales transfer by BC's dilation and
homothety because these are homogeneous linear conditions.
The 10,388-dimensional kernel bounds all smooth closed point values.
The 6,643-dimensional space is realized by the specified invariant
fields. Nonlinear solutions in that space and arbitrary smooth
closed values remain open. External peer review remains outstanding.
