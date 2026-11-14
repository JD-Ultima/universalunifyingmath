# BR. Analytic closed H6 fields from the upper triangular group

On the flat positive observer of [BO](BO-h6-full-curvature-kernel.md),
the $AN$-invariant space of closed real $H_6$ fields has dimension 952.
Each such field has constant norm. Finite sums of rotated fields give
global analytic closed extensions of every point value in a
6,534-dimensional submodule of the curvature kernel. These are closure
results. The nonlinear field equations remain to be tested.

The [verifier](../../computation/geometric-unity/BR_h6_an_closed.py)
reconstructs the connection and sign blocks. The adjacent
[rank and pairing certificate](../../computation/geometric-unity/BR_h6_an_closed_witnesses.json)
contains only the seven real witnesses used for the orbit conclusion,
along with the rank table and modular pairing minors.
[Review82](../../reports/reviews/REVIEW82-GU-BR-2026-11-14.md)
records the scope and independent checks.

## BR1. A simply transitive family

On $SL(5,\mathbb R)/SO(5)$ use the subgroup $AN$ of upper triangular
matrices with positive diagonal and determinant one. It acts simply
transitively by the triangular decomposition. Let $P\in V
=\operatorname{Sym}_0^2(\mathbb R^5)$ and write
$J(P)=P_{\mathrm{strict\ upper}}-P_{\mathrm{strict\ lower}}$.
With BG's tangent normalization, the triangular Lie algebra lift
of $P$ is $(P+J(P))/2$. The Levi-Civita connection in the transported
tangent frame is
$$
\Lambda_P Q=\tfrac12[J(P),Q].
$$
It is metric skew. The difference $\Lambda_PQ-\Lambda_QP$ is the
transported triangular Lie bracket, giving zero torsion. The
verifier checks that its curvature agrees in every coefficient
with $R_{ij}P=-[[P_i,P_j],P]/4$.

For any value $b$ at the identity, transport by the simply
transitive action gives one analytic, constant-norm invariant tensor
field. In its moving frame the components are constant, so its
covariant derivative acts by $-\Lambda$ in each of the seven
covariant slots. Hence the field is a closed hook exactly when
$$
\operatorname{Alt}b=\operatorname{Ctr}b=0,\qquad
(D_{AN}b)_a=\sum_i e^i\wedge
  \big(\rho_{\mathrm{cov},\mathrm{total}}(\Lambda_i)b\big)_a=0.
$$
Here $\rho_{\mathrm{cov},\mathrm{total}}$ is the induced covariant
action applied to all slots. The old BC generator passes the new
closure map exactly.

## BR2. Exact dimension of the invariant fields

The even coordinate-sign subgroup of $SO(5)$ normalizes $AN$. Its
characters split the raw tensor space into sixteen blocks indexed
by the masks of even size on five axes. The exact nullities are

| Characters | Count | Columns per block | Rank per block | Nullity per block |
|:---|---:|---:|---:|---:|
| Empty | 1 | 2672 | 2629 | 43 |
| All two-axis masks | 10 | 2639 | 2584 | 55 |
| Four-axis masks 15 and 30 | 2 | 2596 | 2523 | 73 |
| Four-axis masks 23, 27 and 29 | 3 | 2596 | 2525 | 71 |

The verifier rebuilds the full stack in every selected block.
Modular elimination gives each lower rank bound. Sparse rational
elimination computes a basis of the matching nullity and checks
every vector against all original integer rows. Together the ranks
sum to 41,090 in 42,042 columns, giving dimension 952.
Use `--full-exact` to verify all sixteen blocks.

## BR3. Closed extensions by rotation

Let $W$ denote the 952-dimensional value space. For each $b\in W$,
translate its closed field by $k\in SO(5)$, which fixes the identity
and acts there on the value. Thus evaluation of finite real sums of
these fields contains $M=\operatorname{span}_{\mathbb R}(SO(5)W)$.
Curvature compatibility puts $M$ inside the full kernel of BO.

For each tested type, pair the closed-field values with BO's complex
highest covectors using the Hermitian form induced by the positive
physical metric. The compact complexification has the same
multiplicities. If any copy of a tested type were missing in $M$,
some nonzero combination of its highest vectors would lie in
$M^\perp$ and annihilate the listed evaluations. Full-row-rank
pairing minors exclude this possibility. The minors reduce Gaussian
rationals modulo 65537 using $i=256$ and have nonzero determinant.
Seven exact real closed-field values suffice to verify the
following entire types:
$$
(7,3),(7,2),(7,1),(6,3),(6,2),2(6,1),2(5,2),(5,1),(5,0),
2(4,1),(3,2),(3,0),(2,1).
$$
Their total dimension is 6,534. Selecting a finite spanning subset
of rotated evaluations gives a linear right inverse to evaluation
on this direct sum. Each resulting field is a finite sum of global
analytic closed fields. The proof gives a lower bound on the
rotated span. [BU](BU-h6-second-jet-obstruction.md) excludes the
type $(2,0)$ from all smooth closed field values. The type $(5,5)$
remains unresolved.

## BR4. The larger affine symmetry and its character modes

For the full affine parabolic used in [BC](BC-affine-hooks.md),
allow multiplication by a scalar character whose logarithmic
differential is $z e^{13}$ in the given coframe. Its origin
value is $SO(4)$-fixed, hence is in BG's three-dimensional fixed
kernel with ordered basis $K_0,K_1,K_2$. Here $K_2$ is the type
$(2,0)$ vector. The twisted closure equation is
$(D_{AN}+z e^{13}\wedge)b=0$. In the integer closure normalization,
two rows of the resulting matrix are
$$
-768,b_2,\qquad (1920 z+11520)b_0+(1152-960z)b_1.
$$
They are independent for every complex $z$. Direct substitution in
all rows verifies the kernel vector
$$
(5z-6)K_0+(10z+60)K_1,
$$
which never vanishes and always has zero $K_2$ component. Thus this
character extension of the full affine ansatz cannot obtain that
missing type. The independent exact check is included with
`--parabolic-characters`.

## Field scope

An $AN$-invariant value solving all remaining algebraic branch
equations at the identity would transport them globally, because
its closure is exact and the background is invariant. The fields
of BR3 generally lose the constant-norm property on taking sums.
They solve the linear closure equation with all coefficient and
form slots included. [BV](BV-h6-triangular-positive-obstruction.md)
subsequently excludes the full 952-dimensional invariant family.
General nonlinear feasibility for rotated sums and arbitrary fields
remains open here.
Positive observer scales transfer through BC's dilation and
homothety. External peer review remains outstanding.
