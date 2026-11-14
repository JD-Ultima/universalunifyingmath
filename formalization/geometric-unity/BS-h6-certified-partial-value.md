# BS. A certified real simultaneous partial zero in the H6 kernel

There is a real value of norm one in the full curvature kernel of
[BO](BO-h6-full-curvature-kernel.md) satisfying complete grade-one
scalarity, complete grade-thirteen vanishing, and the two invariant
projections of grades five and nine from
[BQ](BQ-h6-invariant-moments.md). An exact nonzero grade-nine
coefficient excludes this value from the full equation.

The [verifier](../../computation/geometric-unity/BS_h6_partial_value.py)
checks the rational affine generators, their quadratic polynomials,
the contraction certificate, and the separating grade-nine coordinate.

## BS1. Low and high grade identities for every real hook

Use a positive orthonormal frame, Hermitian Clifford generators
$\gamma_a$, and the unit volume $\Gamma^2=-1$. Write $C_a=B_a\Gamma$,
with $B$ a real six-form hook. Its two traces give
$\sum_b\gamma_bC_b=\sum_b C_b\gamma_b=0$. With $S=\sum_b C_b^2$,
the full quadratic of [BC](BC-affine-hooks.md) becomes
$$
q_a=-\frac23\left(\{\gamma_a,S\}
 -\sum_b(\gamma_bC_aC_b+C_b\gamma_a C_b+C_b C_a\gamma_b)\right).
$$
To see this, expand the nested products and use antisymmetry to
sum over all ordered $b,c$ with factor $-2/3$. This gives the
terms $A_{abc}C_bC_c+C_bC_c A_{abc}+C_b A_{abc}C_c$.
Insert $A_{abc}=\gamma_a\gamma_b\gamma_c-\delta_{ab}\gamma_c
+\delta_{ac}\gamma_b-\delta_{bc}\gamma_a$, and use the two hook
sums. In the first term the leading triple reduces to
$2\gamma_a S$ by commuting $\gamma_b,\gamma_c$. In the second
and third terms order the triple as $-\gamma_c\gamma_b\gamma_a$
and $-\gamma_b\gamma_a\gamma_c$, respectively. The leading
triples there vanish by the hook sums. The listed corrections
give the displayed formula.

Set
$$
N=\sum_c\|B_c\|^2,\qquad
E_{ab}=\langle B_a,B_b\rangle,\qquad
I_{ab}=\sum_c\langle\iota_a B_c,\iota_bB_c\rangle.
$$
The degree-one expression is
$$
q_1=-2Ng+\tfrac83E+\tfrac43 I.
$$
Here is a contraction proof of completeness. The map is
$O(14)$-equivariant, since the two volume signs cancel in the
quadratic. The orthogonal invariant pairing theorem spans the
scalar polarization by pairings of the sixteen slots, fourteen
input and two receiver. See [Lehrer and Zhang](https://annals.math.princeton.edu/2012/176-3/p12)
for the pairing formulation. Internal pairings of a single $B$
vanish by alternation or the contraction trace. The two receiver
slots either pair together or attach to different inputs.
A receiver pair leaves the squared norm or the exchanged
internal-external contraction. Attachment to different inputs
gives the external pairing $E$, the mixed contraction
$K_{ab}=\sum_c\langle\iota_cB_a,\iota_bB_c\rangle$, the internal
pairing $I$, or
$J_{ab}=\sum_{c,d}\langle\iota_a\iota_d B_c,\iota_c\iota_bB_d\rangle$.
The alternating trace gives $B_a=\sum_c e^c\wedge\iota_aB_c$,
whence $K=E$ and $I=E-J$. Thus the span reduces to $Ng,E,I$.
The simple two-component hook
$B_0=e^1\wedge K_0$, $B_1=e^0\wedge K_0$ with
$K_0=e^2\wedge\cdots\wedge e^6$ gives the three ordered
diagonal values $q_{77}=-4,q_{22}=-4/3,q_{00}=0$.
The corresponding contraction rows $(2,0,0),(2,0,2),(2,1,1)$
are independent, yielding the formula. Since $\operatorname{tr}E=N$
and $\operatorname{tr}I=6N$, complete scalarity is equivalent to
$2E+I=(4/7)Ng$, with $q_1=-(26/21)Ng$.

For degree thirteen, substitute $C_b=B_b\Gamma$ in the full
formula. One has $C_aC_b=-B_aB_b$ and
$C_b\gamma_aC_b=B_b\gamma_a B_b$. The highest exterior terms
of the two mixed external contractions vanish by the alternating
hook identity. The anticommutator and middle term contribute
$2e^a\wedge S_{12}$ and $e^a\wedge S_{12}$ inside the common
factor $2/3$, where $S_{12}=\sum_b B_b\wedge B_b$. Thus
$q_{13,a}=2e^a\wedge S_{12}$. This vanishes for all $a$ exactly
when $S_{12}=0$.

## BS2. A finite rational slice and complete detectors

The data give nine rational tensors $b_0,\ldots,b_8$ in the
physical BC frame. Every one is a hook and passes the full
coefficient-curvature map. Set $B(x)=b_0+\sum_{j=1}^8 x_jb_j$.
Each generator transforms by the sign of the coordinate permutation
group $S_5$ on the five-dimensional defining space. In particular
the quadratic outputs are invariant.

The character of $V$ at a permutation $p$ is
$(f(p)^2+f(p^2))/2-1$, where $f$ counts fixed axes. Character
averaging gives five symmetric invariant receiver dimensions,
four after removing the metric, and one alternating dimension.
For the four symmetric tracefree basis matrices in the certificate,
the coordinate readings $(0,0),(0,1),(0,7),(0,10)$ have determinant
one. Hence these readings detect the entire tracefree stress.
All adjacent permutations preserve the tangent orientation. The
Hodge dual of the alternating invariant line gives the invariant
twelve-form line, on which mask 15358 is nonzero. Reading that
coefficient of $S_{12}$ detects the entire grade-thirteen equation.

The other three equations are $L_5=L_9=0$ and $N=1$. The two
moment polynomials follow from BQ's calibrated norm ratios on
single-copy Casimir types. The certificate supplies the exact
Casimir components of the affine anchor and the type of each
direction. Within each type, Gram polarization computes the
norm contributions. This reconstructs all eight symmetric
rational polynomial matrices without an approximate optimizer.

## BS3. Contraction and the full-equation obstruction

Let $F(x)$ be the eight polynomials above and $M$ the invertible
rational preconditioner in the data. Set $G(x)=x-MF(x)$. On the
closed sup-norm box of radius $10^{-8}$, the exact derivative bound
is less than $7\cdot10^{-7}$ and the image norm bound is less
than $1.05\cdot10^{-9}$. Thus $G$ is a contraction into the box.
Its fixed point has $F(x)=0$.

Grade nine at external index 3 and increasing blade mask 13135
is a rational quadratic on the same affine slice. Its value at
the center exceeds 100 in the BC coframe. Exact coefficient
bounds give total variation less than $5.84\cdot10^{-6}$ on the
box. The helper reconstructs all 45 polynomial coefficients by
polarization, and the separate nested Clifford engine checks the
center value. The coefficient is positive at the exact partial root.
No approximation to the root is needed for this exclusion.

At norm one the grade-one value is $-(26/21)g$, so a positive
amplitude normalization gives the full grade-one equation on
the branch. The nonzero grade-nine component survives that
normalization, and BA places the derivative in grade seven.
The construction is a pointwise partial solution and supplies
no full field. External peer review remains outstanding.
