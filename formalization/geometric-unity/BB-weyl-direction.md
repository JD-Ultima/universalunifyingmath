# BB. One declared nonzero-Weyl direction

The complete normal quadratic metrics considered here are obstructed for every amplitude in the specified Weyl direction. Two fibre witnesses remove the exceptional roots of a single three-row certificate. Coefficient-rank exceptions are checked separately, so the conclusion excludes fibre-open patches rather than only individual points.

Code: `computation/geometric-unity/BB_weyl_direction.py`; exact output: `BB_weyl_direction.out`; record: `metrics/geometric-unity/BB-weyl-direction.yaml`. The computation uses U's accepted relaxed-system assembly and the earlier referee curvature engine. [Independent internal review](../../reports/reviews/REVIEW21-GU-BB-2026-11-02.md) accepts BB1–BB2 and reproduces the exact calculations. Arbitrary cubic and quartic jets are not included in the present theorem.

## The declared curvature and metric

Let $R_0=\tfrac12(e_0\otimes e_0-\delta/6)\owedge\delta$, as in U. Define the algebraic Weyl tensor $\mathcal W$ by its nonzero sectional components
$$\mathcal W_{0101}=\mathcal W_{2323}=1,\qquad
\mathcal W_{0202}=\mathcal W_{1313}=-1,$$
with $\mathcal W_{0303}=\mathcal W_{1212}=0$, all cross-pair components zero, and the usual curvature symmetries. The first Bianchi identity holds, and each diagonal Ricci component is the sum $1-1+0=0$; off-diagonal Ricci components vanish. Thus $\mathcal W$ is a genuine zero-Ricci Weyl tensor.

For $\rho\ne0$, $\omega\in\mathbb R$, take the complete local metric
$$g_{ij}(x)=\delta_{ij}-\tfrac13(\rho R_0+\omega\mathcal W)_{iajb}x^ax^b.\tag{BB1}$$
It is positive definite near $0$. At the origin $\mathrm{Ric}=\rho e_0\otimes e_0$ and the Weyl tensor is $\omega\mathcal W$. Cubic and quartic metric coefficients are fixed to zero; the connection's higher derivatives are those of this complete polynomial metric.

Use Z1 at $\lambda=1/5$, $c_H,c_V>0$. U's naturality and homothety maps normalize
$$\rho=\sigma\in\{1,-1\},\qquad c_H=2,\quad c_V=3,\qquad
w=\omega/|\rho|\in\mathbb R.\tag{BB2}$$
The analysis below covers every real normalized $w$. It concerns this declared direction and its rotated orthonormal-frame copies, rather than arbitrary Weyl tensors.

## A degree bound before interpolation

**Lemma BB1.** At a fixed fibre matrix $y$, the full relaxed coefficient matrix $A(w)$ has polynomial degree at most two in $w$, and its RHS $b(w)$ has degree at most three. Denominators are independent of $w$.

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

## Result and scope

**Theorem BB2.** For every complete local metric (BB1), with $\rho\ne0$, $\omega\in\mathbb R$, and positive $c_H,c_V$ at $\lambda=1/5$, the relaxed system is consistent on no nonempty fibre-open patch over $0$. In the slaved vector sector of S.0, the twisted cone has no $\Upsilon=0$ solution on an open set meeting that fibre.

*Proof.* Outside $E_\sigma$, (BB5) supplies a witness with a nonzero compatibility determinant and a nonzero coefficient 42-minor. At every exceptional value the exact rank checks give the same two ranks. Hence some augmented 43-minor is nonzero for every real $w$. S4 excludes fibre-open consistency by the rational-minor argument. The parameter maps in (BB2) cover every nonzero $\rho$, every $\omega$ and positive observer scales.

For the field equation, fix the actual parameters and one nonzero augmented minor at its positive coordinate fibre matrix. That fixed minor persists over a base neighbourhood. The base curvature is nonzero at the origin because its Ricci tensor is nonzero, so shrink the neighbourhood to preserve this as well. U5's argument then applies: S5 excludes Einstein pieces, while the Q5 reduction supplies an open variable patch with $H\ne0$ and therefore relaxed consistency, a contradiction. No continuation of a three-row weight relation is assumed. ∎

The recorded full run takes 52.3 s and passes all exact driver assertions. The work consists of four signed witness analyses, polynomial gcd calculations, and twelve remaining coefficient-rank evaluations.

The result is bounded by the declared Weyl direction and complete quadratic metric. It does not yet allow arbitrary cubic or quartic metric coefficients in the nonzero-Weyl cases, and it does not cover arbitrary Weyl tensors or $\rho=0$. The zero-Weyl case $w=0$ has the separate, stronger higher-jet result Z; that does not automatically extend the present theorem to $w\ne0$.

A useful next question is whether Z's cubic-invariant subsystem and its kernel-completion argument can be audited uniformly for this Weyl pencil, including the finite coefficient-minor exceptions. That additional proof is required before passing from the complete quadratic family to arbitrary higher jets.
