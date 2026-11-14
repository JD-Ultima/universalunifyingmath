# CG. Limits of quadratic separation on the full hook kernels

Every linear combination of the homogeneous necessary quadratic
equations that is nonnegative on the entire H4 curvature kernel
vanishes identically there. The same holds on the 6,534-dimensional
H6 submodule of linearly liftable values and on its invariant
superspaces inside the H6 curvature kernel. In each case the
quadratic equations have a strictly positive definite
unit-trace moment matrix. This conclusion concerns quadratic
separation and the linear moment relaxation.

## CG1. Effective invariant receivers

Use the compact group $K=SO(5)$ on
$V=\operatorname{Sym}_0^2(\mathbb R^5)$. For both branches remove
the fixed scalar trace from degree one using [CD](CD-hook-clifford-identities.md).
The degree-one output is then in $\operatorname{Sym}_0^2 V$.
For H4 the top degree-nine output lies in the exterior-product
copy of $\Lambda^8 V$. For H6 the top degree-thirteen output
lies in the exterior-product copy of $\Lambda^{12}V$.

The invariant multiplicities in $\operatorname{Sym}_0^2V$,
$\Lambda^6 V$, $\Lambda^8 V$, $\Lambda^{10}V$ and $\Lambda^{12}V$
are zero. Those in $V\otimes\Lambda^5V$ and
$V\otimes\Lambda^9V$ are each one. In particular there is one
invariant receiver for H4 and two for H6. They are the primitive
matrix trace receivers $L_5,L_9$ of
[BQ](BQ-h6-invariant-moments.md), evaluated with the actual
Clifford input on each branch.

The [verifier](../../computation/geometric-unity/CG_hook_no_separator.py)
computes the multiplicities with integer Weyl character
extraction. Its independent Lie-action count uses the
determinant-one diagonal sign-fixed slots and checks all ten
infinitesimal compact generators. The modular lower ranks
exclude any additional invariant vectors. The nonzero trace
receivers provide the matching lower dimension.

## CG2. Exact witness moments

For H4 the two rational tensors in the certificate satisfy
both hook traces and the full coefficient curvature on the
Hodge dual. Their squared norms and invariant grade-five
moments are

| Witness | Squared norm | $L_5$ | $L_5$ per unit squared norm |
| --- | ---: | ---: | ---: |
| First | $105/256$ | $-25/512$ | $-5/42$ |
| Second | $4985/1024$ | $10475/6144$ | $2095/5982$ |

The complete nested Clifford evaluation and the separate
sorted-word calculation agree on all 2,286 selected receiver
coefficients. Multiplying the two unnormalized Haar averaged
rank-one covariances by $53632/29475$ and $512/9825$ makes
the trace one and $L_5$ zero. Equivalently the probabilities
on the two unit witness orbits are $2933/3930$ and $997/3930$.

For H6 the unit moment pairs $(L_5,384L_9)$ for the three
single-copy types in [BR](BR-h6-closed-field-lifts.md) are

| Type | $L_5$ | $384L_9$ | Haar orbit probability |
| --- | ---: | ---: | ---: |
| $(3,0)$ | $-355/186$ | $184/93$ | $4991/14135$ |
| $(5,0)$ | $35/8$ | $7/3$ | $3292/14135$ |
| $(6,3)$ | $-5/6$ | $-3$ | $532/1285$ |

The probabilities sum to one and annihilate both moments.
The exact affine determinant of the three points is
$-70675/2232$, so they span the invariant output plane
affinely. The selected verifier recalibrates all three types
from their exact highest vectors. Each type is included in
BR's 6,534-dimensional liftable submodule.

## CG3. Consequences for separation

Let $W$ be one of the declared invariant domains. For a real
linear combination $F$ of all homogeneous necessary quadratic
outputs, assume $F\ge0$ on $W$. Its normalized Haar average is
again nonnegative. By CG1 it is one of the indicated invariant
linear combinations. Its values at the unit witnesses are
nonnegative and have weighted sum zero with the strictly
positive weights in CG2. All these values are therefore zero.
The witnesses span the invariant output affinely, so the
averaged form vanishes everywhere on $W$. On every orbit,
the continuous nonnegative function $k\mapsto F(kb)$ then has
zero integral. Hence it is zero at the identity for every
$b\in W$, proving the claimed vanishing.

There is also a strictly feasible covariance for the quadratic
moment relaxation. Start with the exact unit-trace Haar moments
above and assign a small positive weight $\epsilon$ to the
unit-trace identity covariance on $W$. Its output is invariant.
Adjust the witness probabilities by solving the same invertible
affine system, including the equation that the total trace is
one. At $\epsilon=0$ the probabilities are strictly positive.
They remain positive for all sufficiently small $\epsilon>0$.
The identity term makes the covariance positive definite, and
all its required quadratic outputs vanish. The covariance
averages coefficient products. The nonlinear equations constrain
the products of a single field value.

The covariance uses compact orbits. These are contained in the
stated invariant domains, while the fixed triangular spaces
need not contain them. Thus the result is consistent with the
positive receiver on the fixed H6 triangular space in BV.

    python -X utf8 computation/geometric-unity/CG_hook_no_separator.py --full-lie --all-coefficients

[Review96](../../reports/reviews/REVIEW96-GU-CG-2026-11-14.md)
records the independent exact checks. External peer review remains
outstanding.
