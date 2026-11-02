# BB. One declared nonzero-Weyl direction

The obstruction covers every amplitude in one specified Weyl direction, with arbitrary cubic, quartic and smooth higher metric jets. Two fibre witnesses remove the roots of a single compatibility certificate. A cubic-invariant subsystem and its four-dimensional kernel preserve coefficient rank, so the result excludes fibre-open solution patches rather than only individual points.

Code: `computation/geometric-unity/BB_weyl_direction.py` and `BB_cubic_extension.py`, with their `.out` files; record: `metrics/geometric-unity/BB-weyl-direction.yaml`. Both use U's accepted relaxed-system assembly and the earlier referee curvature engine. [Independent internal Review21](../../reports/reviews/REVIEW21-GU-BB-2026-11-02.md) accepts the quadratic certificates BB1–BB2. [Independent internal Review22](../../reports/reviews/REVIEW22-GU-BB-JETS-2026-11-02.md) separately accepts BB3–BB6 and reproduces all 41 extension checks in 33.2 s. These are internal reviews, not external peer review.

## The declared curvature and metric

Let $R_0=\tfrac12(e_0\otimes e_0-\delta/6)\owedge\delta$, as in U. Define the algebraic Weyl tensor $\mathcal W$ by its nonzero sectional components
$$\mathcal W_{0101}=\mathcal W_{2323}=1,\qquad
\mathcal W_{0202}=\mathcal W_{1313}=-1,$$
with $\mathcal W_{0303}=\mathcal W_{1212}=0$, all cross-pair components zero, and the usual curvature symmetries. The first Bianchi identity holds, and each diagonal Ricci component is the sum $1-1+0=0$; off-diagonal Ricci components vanish. Thus $\mathcal W$ is a genuine zero-Ricci Weyl tensor.

For $\rho\ne0$, $\omega\in\mathbb R$, consider
$$g_{ij}(x)=\delta_{ij}-\tfrac13(\rho R_0+\omega\mathcal W)_{iajb}x^ax^b
+Q_{3;ij}(x)+Q_{4;ij}(x)+R_{5;ij}(x),\tag{BB1}$$
where $Q_3,Q_4$ are arbitrary symmetric homogeneous cubic and quartic matrices, and the smooth symmetric remainder $R_5$ has all derivatives through order four zero at $0$. There are 200 free cubic and 350 free quartic coefficients. The metric is positive definite on some neighbourhood of $0$, with no uniform radius over the coefficients. At the origin $\mathrm{Ric}=\rho e_0\otimes e_0$ and the Weyl tensor is $\omega\mathcal W$.

Use Z1 at $\lambda=1/5$, $c_H,c_V>0$. U's naturality and homothety maps normalize
$$\rho=\sigma\in\{1,-1\},\qquad c_H=2,\quad c_V=3,\qquad
w=\omega/|\rho|\in\mathbb R.\tag{BB2}$$
The normalized cubic and quartic tensors are $Q_3/|\rho|^{3/2}$ and $Q_4/\rho^2$; both remain arbitrary. The analysis covers every real normalized $w$ in this declared curvature direction and its rotated orthonormal-frame copies. The full Weyl tensor is not arbitrary.

## A degree bound before interpolation

First set $Q_3=Q_4=R_5=0$. This auxiliary complete quadratic metric supplies the exact polynomial certificates; higher-jet invariance is proved separately below.

**Lemma BB1.** For this complete quadratic metric, at a fixed fibre matrix $y$, the full relaxed coefficient matrix $A(w)$ has polynomial degree at most two in $w$, and its RHS $b(w)$ has degree at most three. Denominators are independent of $w$.

*Proof.* At $x=0$, the base metric is $I$, its first derivative vanishes, and its quadratic coefficient is affine in $w$. For the complete quadratic base metric, the connection's linear term is affine in $w$, its quadratic term is zero, and its cubic term has degree at most two in $w$. The observer metric at the point, and its inverse, are independent of $w$. Its first derivatives have degree at most one, and its second and third derivatives degree at most two: use $h_{HV}=-c_VG_y(\Psi,\cdot)$ and $h_{HH}=c_Hy+c_VG_y(\Psi,\Psi)$ with $\Psi(0)=0$.

Riemann and Ricci are linear in the second observer derivatives plus quadratic terms in the first derivatives. Their degree is at most two. The covariant derivative of Ricci contains third derivatives, products of first and second derivatives, and cubic products of first derivatives; its degree is at most three. Raising indices uses the fixed observer inverse metric. These are exactly the tensors entering the relaxed matrix and RHS. ∎

Accordingly, three exact nodes $w=-1,0,1$ determine every coefficient polynomial, and four nodes $w=-1,0,1,2$ determine every RHS polynomial. These bounds are established before interpolation. The fourth node also checks the entire reconstructed coefficient matrix, entry by entry. Coefficient determinants are then calculated over $\mathbb Q[w]$ from the reconstructed entries; their degree is not guessed from determinant samples.

## The two three-row witnesses

Use fibre matrices
$$y_1=\operatorname{diag}(2,3,4,5),\qquad y_2=\operatorname{diag}(2,3,5,7).$$
The three rows are again $I=(168,337,492)$, or $(a,b,d)=(0,13,0),(1,13,1),(2,13,2)$. At either witness they have only the common scalar $L_y(n)=h(n,\partial_{y_{33}})$ and $\xi_{13}$ as unknown combinations:
$$L_y(n)+a_i(w)\xi_{13}=b_i(w),\qquad i=0,1,2.\tag{BB3}$$
The unwanted coefficient polynomials have degree at most two and vanish at three distinct nodes, so their absence is certified for all $w$. The $n$ coefficients are independent of curvature directly from the definition of the relaxed system.

At $y_1$ the exact row polynomials are
$$\begin{aligned}
a_0&=-\frac{2925w^2-750\sigma w+18649-10800\sigma}{43200},&
b_0&=\frac1{500},\\
a_1&=-\frac{1968w^2+1056\sigma w+24163}{57600},&
b_1&=\frac{(6w+\sigma)^2}{4500},\\
a_2&=-\frac{4212w^2-2754\sigma w+36467}{86400},&
b_2&=\frac{(6w-\sigma)^2}{12000}.
\end{aligned}$$
At $y_2$ they are
$$\begin{aligned}
a_0&=-\frac{46746w^2-16464\sigma w+231319-132300\sigma}{529200},&
b_0&=\frac5{4116},\\
a_1&=-\frac{8475w^2+4050\sigma w+74107}{176400},&
b_1&=\frac{(6w+\sigma)^2}{6174},\\
a_2&=-\frac{37017w^2-24246\sigma w+224701}{529200},&
b_2&=\frac{(6w-\sigma)^2}{20580}.
\end{aligned}$$
The reconstructed cubic coefficients of these RHS polynomials are zero. This is an exact consequence of the four-node calculation within BB1's bound, rather than an additional assumed degree reduction.

Let $P_{\sigma,j}$ denote the numerator of $\det[(a_i,1,b_i)]_{i=0}^2$ at $y_j$, with denominators $414720000$ and $1555848000$ respectively. The numerators are

| sign | witness | $P_{\sigma,j}(w)$ |
|---|---|---|
| $+$ | $y_1$ | $21168w^4+127440w^3-487728w^2-403540w-11799$ |
| $+$ | $y_2$ | $57132w^4+372708w^3-1482021w^2-1051750w-33734$ |
| $-$ | $y_1$ | $21168w^4-127440w^3+549072w^2-356780w+17001$ |
| $-$ | $y_2$ | $57132w^4-372708w^3+1693179w^2-913850w+54466$ |

For each sign, exact polynomial gcd gives
$$\gcd(P_{\sigma,1},P_{\sigma,2})=1.\tag{BB4}$$
Thus at least one three-row witness is inconsistent at every real $w$. The roots of a single compatibility quartic are exceptions of that certificate, not solutions of the field equation.

## The full-rank step and its finite exceptions

Pointwise inconsistency alone would not establish fibre-open exclusion if the coefficient matrix lost rank at the witness. The coefficient audit is therefore separate.

Use the 38-row list recorded in Z's output. Append either
$$J_5=(61,229,383,523),\qquad J_6=(76,244,398,538),$$
the four rows $(i,5,5)$ or $(i,6,6)$, $i=0,1,2,3$. At each $y_j$, let $D_{\sigma,j,5}(w)$ and $D_{\sigma,j,6}(w)$ be the resulting coefficient determinants. The code computes these exactly over $\mathbb Q[w]$ and prints their complete rational factorizations. The rank-38 kernel of Z is not assumed to remain valid throughout this Weyl family.

If any product $D_{\sigma,j,k}P_{\sigma,j}$ is nonzero, the same witness has coefficient rank 42 and augmented rank 43. The square-free gcd of the four products, up to a nonzero constant, is
$$\begin{aligned}
G_+&=(3w-1)(3w+1)(6w-1)(6w+1)(21w+13)(39w+7),\\
G_-&=(3w-1)(3w+1)(6w-1)(6w+1)(21w-1)(39w-19).
\end{aligned}\tag{BB5}$$
All remaining possible failures of this combined minor argument are therefore rational and finite:
$$\begin{aligned}
E_+&=\{-13/21,-1/3,-7/39,-1/6,1/6,1/3\},\\
E_-&=\{-1/3,-1/6,1/21,1/6,1/3,19/39\}.
\end{aligned}$$
For every member of $E_\sigma$, exact evaluation of the full reconstructed coefficient matrix at $y_1$ gives rank 42, while $P_{\sigma,1}\ne0$. The output records all twelve checks. These are evaluations of coefficient polynomials certified by BB1, rather than extrapolations from nearby amplitudes. No exceptional real amplitude remains.

## The auxiliary quadratic result

**Theorem BB2.** For (BB1) with $Q_3=Q_4=R_5=0$, $\rho\ne0$, $\omega\in\mathbb R$, and positive $c_H,c_V$ at $\lambda=1/5$, the relaxed system is consistent on no nonempty fibre-open patch over $0$. In the slaved vector sector of S.0, the twisted cone has no $\Upsilon=0$ solution on an open set meeting that fibre.

*Proof.* Outside $E_\sigma$, (BB5) supplies a witness with a nonzero compatibility determinant and a nonzero coefficient 42-minor. At every exceptional value the exact rank checks give the same two ranks. Hence some augmented 43-minor is nonzero for every real $w$. S4 excludes fibre-open consistency by the rational-minor argument. The parameter maps in (BB2) cover every nonzero $\rho$, every $\omega$ and positive observer scales.

For the field equation, fix the actual parameters and one nonzero augmented minor at its positive coordinate fibre matrix. That fixed minor persists over a base neighbourhood. The base curvature is nonzero at the origin because its Ricci tensor is nonzero, so shrink the neighbourhood to preserve this as well. U5's argument then applies: S5 excludes Einstein pieces, while the Q5 reduction supplies an open variable patch with $H\ne0$ and therefore relaxed consistency, a contradiction. No continuation of a three-row weight relation is assumed. ∎

The recorded quadratic run takes 52.3 s. It comprises four signed witness analyses, polynomial gcd calculations and twelve remaining coefficient-rank evaluations. Review21 independently reproduced these calculations in 54.1 s. The following argument supplies the additional proof needed for arbitrary higher jets.

## Higher-jet invariance

Write $H=\{0,1,2,3\}$, $V=\{4,\ldots,13\}$ and $D=\{4,8,11,13\}$ for horizontal, vertical and diagonal-vertical coordinates. The row convention remains $(a,b,d)$.

**Lemma BB3.** At a diagonal positive fibre matrix, the three coefficient rows and RHS rows in (BB3) are independent of every $Q_3,Q_4,R_5$. All coefficient rows in the following 722-row subsystem are also independent of them:

| row class | count |
|---|---:|
| $V,V,V$ | 450 |
| $H,D,H$ | 64 |
| $H,V,D$ | 160 |
| $H,H,D$ | 24 |
| $D,D,H$ | 24 |

*Proof.* The jet-order, support and parity arguments of Z1–Z4 and W1–W2 do not require the fixed quadratic curvature to have zero Weyl tensor. Here is their dependence on that hypothesis explicitly. For any fixed quadratic metric jet, the connection has fixed linear part, degree-two part linear in $Q_3$, and degree-three part equal to a fixed curvature contribution plus a term linear in $Q_4$. Thus the observer first jet is fixed, its second jet is affine in $Q_3$, and its third jet is affine jointly in $Q_3,Q_4$. The coefficient matrix is affine in $Q_3$ and independent of $Q_4$; the RHS is affine jointly in $Q_3,Q_4$. The zero-four-jet remainder affects neither.

Only mixed $HV$ second metric derivatives vary with $Q_3$. For a diagonal vertical slot $A=y_{aa}$, the degree-two mixed metric variation is a base gradient,
$$\delta h_{jA}=\partial_jF_a,\qquad
F_a=-\frac{c_V}{d_a}(Q_{3;aa}-\lambda\operatorname{tr}Q_3).$$
Commuting third derivatives eliminates the Riemann variation with three horizontal slots and this diagonal vertical slot. Hence the only cubic curvature variation has three horizontal slots and one off-diagonal vertical slot; the Ricci variation has one horizontal and one off-diagonal vertical slot. Raising indices preserves this support at diagonal $y$. Substitution in the relaxed system fixes all 722 listed coefficient rows and the three selected rows.

Base reflection $(x,y)\mapsto(-x,y)$ reverses $Q_3$ and fixes the quadratic and quartic jets. Each selected RHS row has two horizontal indices and therefore even parity. Affine cubic dependence then makes its cubic variation zero. For the quartic variation, replacing $Q_3$ by $Q_4$ gives the leading degree-three mixed metric variation. With lower observer jets fixed, its only contribution is the principal third-derivative Ricci term
$$\delta(\nabla_i\mathrm{Ric})_{iA}
=\tfrac12\sum_{c=0}^3h^{cc}
\left(\partial_c\partial_i^2\partial_cF_a-
\partial_c^2\partial_i\partial_iF_a\right)=0.$$
The companion term $\delta(\nabla_A\mathrm{Ric})_{ii}$ vanishes by base order. This argument is independent of the fixed curvature, including its Weyl component. There are no mixed cubic-quartic terms. ∎

## The invariant kernel and rank completion

Let $M_{\sigma,j}(w)$ be the 722-row coefficient subsystem at $y_j$ for the auxiliary quadratic metric. For $y=\operatorname{diag}(d_0,d_1,d_2,d_3)$ define four independent sparse columns
$$K_i=e_{28+i}+\frac{24d_i}{5}e_{14+i},\qquad i=0,1,2,3.\tag{BB6}$$
The columns are in the ordering $(v^0,\ldots,v^{13},\xi_0,\ldots,\xi_{13},n^0,\ldots,n^{13})$. Thus they describe $v=0$, $\xi_V=n_V=0$, and $\xi_i=c h_{ii}n^i$ with $c=12/5$ and $h_{ii}=2d_i$.

**Proposition BB4.** At the witness supplied by the quadratic certificates for each real $w$, the full coefficient matrix has rank 42 for every $Q_3,Q_4,R_5$.

*Polynomial kernel certificate.* At both $y_j$ and both signs, exact evaluation gives $M_{\sigma,j}(w)K=0$ at $w=-1,0,1$. Every entry has degree at most two by BB1, so
$$M_{\sigma,j}(w)K=0\quad\text{for all }w.\tag{BB7}$$
This is a polynomial identity from a proved degree bound, rather than a sampled rank assertion. Since the four columns are independent, $\operatorname{rank}M\le38$.

The exact restrictions of the two four-row completions are
$$J_5K=-(3w+\sigma)^2\Delta_{j,5},\qquad
J_6K=-(3w-\sigma)^2\Delta_{j,6},\tag{BB8}$$
where all diagonal entries below are positive:

| witness | $\Delta_{j,5}$ | $\Delta_{j,6}$ |
|---|---|---|
| $y_1$ | $\operatorname{diag}(1/45,1/30,2/45,1/18)$ | $\operatorname{diag}(1/20,3/40,1/10,1/8)$ |
| $y_2$ | $\operatorname{diag}(1/45,1/30,1/18,7/90)$ | $\operatorname{diag}(9/125,27/250,9/50,63/250)$ |

Again three exact nodes determine these identities within BB1's degree bound. Since $\sigma=\pm1$, their two scalar factors never vanish simultaneously.

*Why cubics preserve a completed minor.* In a row $(i,O,O)$ with off-diagonal vertical coordinate $O=5$ or $6$, cubic curvature support excludes Riemann variation. The only variable Ricci coefficient is that of $\xi_O$. Every column in $K$ has $\xi_O=0$, so the row variation annihilates $K$. Quartics change no coefficients. Thus (BB8) is unchanged by all higher jets.

For $w\notin E_\sigma$, the quadratic product-gcd argument supplies a nonzero $D_{\sigma,j,k}P_{\sigma,j}$. Its first 38 invariant rows are independent. By (BB7), their kernel is exactly $\operatorname{span}K$. A cubic change of an appended row annihilates $K$ and therefore belongs to the span of those 38 rows. Row operations leave that same 42-minor unchanged. No new polynomial determinant interpolation is needed.

At every one of the twelve values $w\in E_\sigma$, exact rational evaluation at $y_1$ gives $\operatorname{rank}M_{\sigma,1}=38$; all twelve ranks are recorded by the extension script. Thus its kernel is exactly $\operatorname{span}K$, even if the originally recorded 38-row list becomes dependent. One of (BB8)'s two completions is invertible on this kernel, and remains so for every higher jet. Appending it to any 38 independent invariant rows gives rank 42. The original quadratic computation already checks $P_{\sigma,1}\ne0$ at these values. ∎

## Smooth and intrinsic consequences

**Theorem BB5.** For every smooth metric (BB1), arbitrary $Q_3,Q_4,R_5$, $\rho\ne0$, $\omega\in\mathbb R$, and positive $c_H,c_V$ at $\lambda=1/5$, the relaxed system is consistent on no nonempty fibre-open patch over $0$. Under the slaved-vector hypotheses of S.0, the twisted cone has no $\Upsilon=0$ solution on an open set meeting that fibre.

*Proof.* Outside $E_\sigma$, use the witness with $D_{\sigma,j,k}P_{\sigma,j}\ne0$. At the twelve remaining amplitudes use $y_1$. Lemma BB3 preserves the nonzero compatibility determinant, while BB4 gives coefficient rank 42 for every higher jet. Consequently the augmented rank is 43. S4's rational fibre-minor argument excludes any fibre-open consistency patch. Naturality and homothety transfer the ranks to all original parameters.

For a fixed actual metric, select one nonzero augmented 43-minor at its positive coordinate fibre witness. Its nonvanishing persists over a base neighbourhood. The nonzero base Ricci tensor at $0$ allows that neighbourhood to be chosen with nonzero curvature. S5 then excludes Einstein pieces; the Q5 reduction on its dense variable patch with $H\ne0$ supplies relaxed consistency, contradicting the fibrewise obstruction as in U5. The continued object is a fixed matrix minor, rather than the original three-row weights. The neighbourhood may depend on the entire metric jet; no uniform radius is asserted. ∎

**Corollary BB6 (intrinsic scope).** Let a smooth Riemannian four-metric have an orthonormal frame at $x_0$ in which
$$\mathrm{Ric}_{x_0}=\rho e_0\otimes e_0,\qquad
W_{x_0}=\omega\mathcal W,\qquad \rho\ne0.$$
Then BB5 holds over $x_0$. Rotated copies rotate the distinguished Ricci axis and the specified Weyl tensor together.

*Proof.* The curvature decomposition gives $R_{x_0}=\rho R_0+\omega\mathcal W$. Riemann normal coordinates have precisely the quadratic term in (BB1). Their actual cubic and quartic terms belong to the arbitrary class just proved, and the smooth Taylor remainder has zero four-jet. Observer naturality completes the transfer. ∎

The extension script checks the exact polynomial kernel and completion identities, all twelve exceptional invariant ranks, and mixed nonzero-Weyl cubic/quartic controls. The controls change other coefficients and RHS rows while preserving the proved subsystem, selected rows and completed minor; reflection and affine dependence are checked on the full control system. Universality follows from BB3–BB4, not these controls. No repeated polynomial 42-determinant computation is required by the extension script.

The remaining boundary is arbitrary Weyl curvature outside this declared simultaneous orthonormal-frame orbit, and the Ricci-flat nonflat case $\rho=0$. The result places no additional higher-jet restriction inside its stated curvature orbit.
