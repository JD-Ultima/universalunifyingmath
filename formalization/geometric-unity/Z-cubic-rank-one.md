# Z. The cubic extension and an invariant curvature condition

**Broader curvature class.** [CC](CC-diagonal-weyl.md) treats the complete diagonal Weyl plane, equivalently vanishing mixed Weyl components relative to the nonzero Ricci line. This note retains the simpler certificate for its subcase.

The obstruction survives arbitrary cubic metric coefficients as well as arbitrary quartic coefficients. The prescribed quadratic jet still has nonzero rank-one Ricci tensor and zero Weyl tensor. A fixed coefficient minor proves full rank for every cubic tensor, so the result excludes fibre-open solution patches rather than only the evaluated point.

Code: `computation/geometric-unity/Z_cubic_rank_one.py`; output: `Z_cubic_rank_one.out`; record: `metrics/geometric-unity/Z-cubic-rank-one.yaml`. The code uses the accepted U assembly and the earlier referee engine. The rank certificates are exact rational calculations; the support, parity and normalization arguments establish the parameter ranges. It was independently accepted in `reports/reviews/REVIEW19-GU-Z-2026-11-02.md`, including the coordinate-invariant corollary and a 35-assertion rerun. This is internal review, not external peer review.

## Setting

Near $0\in\mathbb R^4$, let
$$g_{ij}(x)=\delta_{ij}-\frac{\rho}{3}(R_p)_{iajb}x^ax^b
+Q_{3;ij}(x)+Q_{4;ij}(x)+R_{5;ij}(x),\qquad
R_p=\tfrac12(p\otimes p-\delta/6)\owedge\delta,\tag{Z1}$$
where $|p|=1$, $\rho\ne0$, $Q_3,Q_4$ are arbitrary symmetric homogeneous cubic and quartic matrices, and all derivatives of the smooth symmetric remainder $R_5$ through order four vanish at $0$. There are 200 free cubic and 350 free quartic coefficients. The metric is positive definite on a sufficiently small neighbourhood, with no uniform radius over those coefficients.

Use Z1 at $\lambda=1/5$ and $c_H,c_V>0$, with the relaxed system (U2). Rotations, the lifted dilation and the observer homothety of U reduce the problem to
$$\rho=\sigma\in\{1,-1\},\qquad p=e_0,\qquad c_H=2,\quad c_V=3.$$
The normalized cubic and quartic tensors are $Q_3/|\rho|^{3/2}$ and $Q_4/\rho^2$. Both remain arbitrary. Higher remainders do not affect the observer 3-jet at the fibre over $0$.

All coefficient certificates below are at $y_* =\operatorname{diag}(2,3,4,5)$. Write $H=\{0,1,2,3\}$ for horizontal coordinate indices, $V=\{4,\ldots,13\}$ for vertical indices, and $D=\{4,8,11,13\}$ for the diagonal vertical coordinates $y_{00},y_{11},y_{22},y_{33}$.

## Which cubic variations can occur

**Lemma Z1 (jet dependence).** The observer metric and its first jet at $(0,y)$ are fixed. Its second jet is affine in $Q_3$ and independent of $Q_4$. Its third jet is affine jointly in $Q_3,Q_4$. Consequently the relaxed coefficient matrix is affine in $Q_3$, while its RHS is affine jointly in $Q_3,Q_4$.

*Proof.* The base connection has fixed linear part, degree-two part linear in $Q_3$, and degree-three part equal to a fixed quadratic-curvature contribution plus a term linear in $Q_4$. A cubic metric term cannot enter the degree-three connection coefficient: the inverse-metric corrections begin at degree two and the first base derivative vanishes at the origin. In the observer metric, $h_{HV}=-c_VG_y(\Psi,\cdot)$ and $h_{HH}=c_Hy+c_VG_y(\Psi,\Psi)$ with $\Psi(0)=0$. Cubics first change $h_{HV}$ at order two and $h_{HH}$ at order three. Quartics first change them at orders three and four. No product of two variable cubic contributions occurs below order four. The standard curvature formulas then give the stated affine dependence. ∎

**Lemma Z2 (cubic support).** At a diagonal fibre point, the cubic variation of the covariant Riemann tensor has only components with three horizontal indices and one off-diagonal vertical index. The cubic Ricci variation has only horizontal/off-diagonal-vertical components. These statements also hold for the corresponding raised components.

*Proof.* Only the mixed $HV$ second metric derivatives vary, and both derivative indices are horizontal. Thus the principal Riemann variation has exactly three horizontal and one vertical index. For a diagonal vertical coordinate $A=y_{aa}$, its mixed metric variation is a base gradient:
$$\mathcal H_{jA}=\partial_jF_a,\qquad
F_a=-\frac{c_V}{d_a}\big(Q_{3;aa}-\lambda\operatorname{tr}Q_3\big).\tag{Z2}$$
This is the cubic version of W's diagonal-vertical calculation. The principal components with three horizontal indices and this vertical slot are
$$\delta R_{ijkA}=\tfrac12\big(\partial_j\partial_k\mathcal H_{iA}
-\partial_i\partial_k\mathcal H_{jA}\big)=0,$$
because third derivatives of $F_a$ commute. Riemann symmetries put the diagonal vertical slot anywhere. Contraction with the block-diagonal inverse metric leaves only $H/V$ Ricci components, and a diagonal vertical slot again gives zero. The inverse fibre metric at diagonal $y$ preserves its diagonal and off-diagonal subspaces, so raising indices preserves this support. ∎

**Lemma Z3 (the certificate RHS).** The RHS rows $(a,b,d)=(i,13,i)$, $i=0,1,2$, are independent of all cubic and quartic coefficients.

*Proof.* Base reflection $(x,y)\mapsto(-x,y)$ sends $Q_3$ to $-Q_3$ while fixing the quadratic and quartic jets. The tensor $(\nabla_a\widehat{\mathrm{Ric}})^d{}_b-(\nabla_b\widehat{\mathrm{Ric}})^d{}_a$ has even reflection parity in each selected row, since it has two horizontal indices. Its dependence on $Q_3$ is affine by Z1, so its linear cubic variation vanishes. W2 proves that its quartic variation in the same rows vanishes. There are no mixed cubic-quartic terms by Z1. ∎

The coefficients in these three rows are also fixed. Their Riemann terms have either two horizontal/two vertical slots, or a diagonal vertical slot; Z2 excludes both variations. Their Ricci terms are $HH$ or $H/D$, also fixed, and the metric coefficients of $n$ are unchanged. Thus W's three-row contradiction remains valid for every $Q_3,Q_4$.

## Full rank requires a separate certificate

The preceding contradiction by itself proves inconsistency at the witness. Since the full coefficient matrix changes with $Q_3$, full rank needs its own argument before applying the fibre-minor continuation.

**Proposition Z4 (a cubic-invariant subsystem).** Every coefficient row in the following disjoint classes is fixed by all cubic tensors:

| row type $(a,b,d)$ | count |
|---|---:|
| $V,V,V$ | 450 |
| $H,D,H$ | 64 |
| $H,V,D$ | 160 |
| $H,H,D$ | 24 |
| $D,D,H$ | 24 |
| total | 722 |

*Proof.* Apply Z2 to the Riemann and Ricci terms of (U2). A Riemann variation needs exactly three horizontal slots and an off-diagonal vertical slot. Each listed class prevents that combination or puts its sole vertical slot in $D$. A Ricci variation requires an $H$ and an off-diagonal vertical slot; the listed output/input pairs prevent it. The metric coefficients of $n$ are fixed. Quartics do not change any coefficient row by W1. ∎

**Proposition Z5 (exact rank and kernel).** For each canonical sign, the fixed 722-row subsystem has rank 38. Its kernel is exactly
$$v=0,\qquad\xi_V=n_V=0,\qquad
\xi_i=c\,h_{ii}n^i\quad(i\in H),\qquad c=12/5.\tag{Z3}$$

*Certificate.* The output records 38 independent row indices from the subsystem. Rational elimination gives rank 38. The four sparse vectors specified by $n^i=1$ separately in (Z3), with all other $n^j=0$, annihilate all 722 rows and are independent. Their nonzero $(\xi_i,n^i)$ entries are
$$ (48/5,1),\quad(72/5,1),\quad(96/5,1),\quad(24,1).$$
Rank-nullity proves that these four vectors span the entire kernel. The coefficient rows are fixed by Z4, so this kernel is the same for every cubic tensor. ∎

**Proposition Z6 (four extra rows give rank 42).** Append rows
$$J=(61,229,383,523),\qquad(a,b,d)=(i,5,5),\ i=0,1,2,3.$$
Their restriction to the kernel basis in Z5 is, for either sign,
$$\operatorname{diag}(-1/45,-1/30,-2/45,-1/18).\tag{Z4}$$
Consequently the full coefficient matrix has rank 42 for every cubic and quartic tensor.

*Proof and certificate.* In each extra row, the Riemann coefficient variation is zero by its horizontal/vertical count. The only variable Ricci coefficient is that of $\xi_5$; all other coefficients are fixed. Since $\xi_V=0$ in (Z3), this variation vanishes on the entire kernel. At zero cubic tensor, exact restriction gives (Z4). Equivalently, the fixed vertical Einstein defect is
$$\left(\widehat{\mathrm{Ric}}+I/c\right)^5{}_5=1/432,$$
and the restriction is $-c h_{ii}/432=-d_i/90$. Thus the four extra rows annihilate no nonzero vector of the kernel. Together with the rank-38 subsystem they have rank 42 for all coefficients.

The recorded 38 independent rows followed by $J$ define a fixed square coefficient minor. Every cubic change in an appended row annihilates the subsystem kernel, so it lies in the span of the first 38 rows. Row operations therefore leave this determinant unchanged. The output records its exact nonzero rational value for each sign, independently of cubic and quartic coefficients. ∎

## Consequences

**Theorem Z7 (fibre-open exclusion).** For every metric jet (Z1), the relaxed system is consistent on no nonempty open subset of the fibre over $0$. Under the slaved-vector hypotheses of S.0, the twisted cone has no $\Upsilon=0$ solution on an open set meeting that fibre.

*Proof.* Z3 and the unchanged coefficient rows retain W's exact three-row contradictions, with RHS defects $3933/800$ and $5667/800$ for the two canonical signs. Z6 gives coefficient rank 42, so the augmented rank is 43 at the witness. U's coordinate and homothety maps transfer both ranks to every $\rho\ne0$, direction $p$ and positive $c_H,c_V$. S4's rational-minor argument then excludes fibre-open consistency.

For the field equation, fix the actual metric jet and one nonzero augmented 43-minor at the witness. Its nonvanishing persists over a base neighbourhood at the same positive coordinate fibre matrix. Shrink the neighbourhood so that the nonzero base curvature also persists. U5's argument applies there: S5 excludes Einstein pieces, and Q5's dense variable patch supplies an open line-equation domain with $\mu\ne0$, hence relaxed consistency. This contradicts the fibre-open obstruction. It is the fixed minor, rather than the three-row weights, that is continued to the neighbourhood. ∎

**Corollary Z8 (coordinate-invariant scope).** Let a smooth Riemannian four-metric have, at a point $x_0$, a nonzero rank-one Ricci tensor and zero Weyl tensor. For Z1 at $\lambda=1/5$, $c_H,c_V>0$, the conclusion of Z7 holds over $x_0$.

*Proof.* Choose an orthonormal frame in which $\mathrm{Ric}_{x_0}=\rho p\otimes p$, $|p|=1$, $\rho\ne0$. The four-dimensional curvature decomposition with zero Weyl tensor gives $R_{x_0}=\rho R_p$. Riemann normal coordinates have the quadratic expansion in (Z1); their actual cubic and quartic terms belong to the arbitrary class covered by Z7. The remaining smooth Taylor remainder has zero four-jet. Naturality of the observer metric and relaxed system completes the transfer. ∎

This removes the need to assume the special complete polynomial metric of U or the zero cubic jet of W. It still does not cover nonzero Weyl curvature or the Ricci-flat nonflat case.

## Verification and remaining boundary

The recorded run passes 35 assertions in 13.8 s. The symbolic calculation includes all 200 cubic coefficients and verifies every diagonal-vertical gradient identity in (Z2), including the vanishing three-horizontal Riemann variation. Exact controls with mixed cubic tensors change many other coefficients and RHS rows while preserving the 722-row subsystem, the three certificate rows and the fixed coefficient determinant. Reflection and affine-dependence identities are checked on the entire mixed control system. A cubic-plus-quartic control checks the combined scope.

The unresolved rank-$\mathrm{Ric}\le1$ frontier can still include points with nonzero Weyl tensor and Ricci-flat nonflat points. Cubic and quartic freedom at a nonzero rank-one-Ricci, zero-Weyl point is no longer an obstruction to this particular proof. A further extension would need a declared nonzero Weyl family, or a separate argument for the Ricci-flat case; neither follows from these certificates.
