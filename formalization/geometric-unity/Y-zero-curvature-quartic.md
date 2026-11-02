# Y. A fibre obstruction at a zero-curvature base point

A nonzero quartic conformal term can obstruct the relaxed system even when the base curvature vanishes at the origin. This note proves an obstruction throughout the positive metric fibre over that point. It does not infer a local full-field $\Upsilon$ theorem from the rank-deficient witness.

Code: `computation/geometric-unity/Y_zero_curvature_quartic.py`; exact output: `Y_zero_curvature_quartic.out`; record: `metrics/geometric-unity/Y-zero-curvature-quartic.yaml`. The calculation uses W's point helper, U's relaxed-system assembly and the accepted referee curvature engine. Independently reviewed in `reports/reviews/REVIEW17-GU-Y-2026-11-02.md`; all 12 assertions were reproduced. This is internal review, not external peer review.

## The declared four-jet

For $|p|=1$ and $\beta\ne0$, take a smooth local base metric

$$g_{ij}(x)=\big(1+\beta(p\cdot x)^4\big)\delta_{ij}+R_{5;ij}(x),\tag{Y1}$$

where the symmetric remainder has zero derivatives through order four at $0$. Use Z1 with $\lambda=1/5$ and $c_H,c_V>0$. The metric is positive definite on some neighbourhood of $0$. Its Riemann tensor vanishes at $0$; this is not an assertion that the metric is flat or Ricci-flat on a neighbourhood.

A rotation reduces $p$ to $e_0$. At every $(0,y)$, the observer 2-jet is exactly that of the flat-base Z1 metric: the quartic term first changes the base connection at cubic order. Thus the relaxed coefficient matrix satisfies $A(\beta,y)=A(0,y)$, and the RHS is affine in $\beta$. The zero-four-jet remainder does not contribute. These are W1's jet-order arguments, which do not require nonzero quadratic curvature.

The flat control is Einstein at $\lambda=1/5$, so $b(0,y)=0$. Therefore

$$b(\beta,y)=\beta b(1,y).\tag{Y2}$$

W's principal-jet calculation of the six nonzero slopes applies without a quadratic curvature term. At the canonical parameters $c_H=2,c_V=3$ and $y_*=\operatorname{diag}(2,3,4,5)$, those slopes are $3/2,3,9/2,1/4,9/32,27/100$ at rows $5,20,35,57,72,87$ respectively. This dependence is proved from the metric jets, rather than fitted to samples.

## An exact three-row contradiction

Use U's 42 unknowns $(v^i,\xi_i,n^i)$ and its zero-based row ordering $(a,b,d)$. Rows $5,20,191$ are $(0,1,5),(0,2,6),(1,2,9)$. Their equations are

$$\frac14v^4-\frac16v^8=\frac32\beta,\qquad
\frac13v^4-\frac16v^{11}=3\beta,\qquad
\frac13v^8-\frac14v^{11}=0.\tag{Y3}$$

All other coefficients in these rows vanish. Four times the first equation, minus three times the second, plus twice the third gives $0=-3\beta$. Hence every nonzero real $\beta$ is inconsistent at this witness, for either sign.

The exact coefficient and augmented ranks are $(28,29)$. The output records 29 rows and columns of an augmented minor whose determinant at $\beta=1$ is

$$D=\frac{30517578125}{13249474533898474022240256}\ne0.\tag{Y4}$$

Its determinant is $\beta D$ for arbitrary $\beta$, since only its final column changes. This supplies a full recorded minor as well as the smaller row contradiction.

## Why no fibre-open consistent patch exists

**Proposition Y1.** The relaxed coefficient rank is 28 throughout the positive fibre over $0$.

*Proof.* Its coefficient matrix is that of flat-base Z1 at every $y$. For the flat connection, the maps

$$\Phi_B(x,y)=(Bx,B^{-\top}yB^{-1}),\qquad B\in GL(4),$$

are isometries of Z1. The horizontal term $c_Hy(dx,dx)$ is invariant, and the vertical DeWitt pairing is invariant under congruence. These maps fix the base origin and act transitively on positive-definite fibre matrices. Tensor covariance therefore preserves the relaxed coefficient rank throughout that fibre. Its exact rank at $y_*$ is 28. ∎

**Proposition Y2.** For every metric (Y1), the relaxed system is consistent on no nonempty open subset of the fibre over $0$.

*Proof.* At the canonical scales, Proposition Y1 and the nonzero augmented 29-minor give the result by S4's rational-minor argument. Consistency on a fibre-open patch would make this minor vanish there, because the coefficient rank is 28 everywhere. Its entries are rational in the fibre coordinates; clearing their nonzero denominators would make the numerator polynomial identically zero, contrary to (Y4).

Rotations cover every unit $p$. U's observer homothety and positive fibre scaling cover every $c_H,c_V>0$: with $\alpha=3/c_V$ and $q=2c_V/(3c_H)$, the map $(x,y)\mapsto(x,qy)$ pulls $\alpha h(c_H,c_V)$ back to $h(2,3)$. These maps preserve relaxed consistency and rank. ∎

In particular, a line equation with its coefficient $\mu$ nonzero everywhere cannot hold on an open set meeting this fibre, since it would supply a fibre-open consistent patch.

## Limit of the conclusion

This is not W's full-rank inconsistent witness. The coefficient rank may rise above 28 away from the origin. A nonzero augmented 29-minor therefore does not preserve inconsistency under arbitrary base-point or metric-jet perturbations. No base neighbourhood of inconsistent fibres, or local full-field $\Upsilon=0$ exclusion, is claimed here. A reduction whose line coefficient vanishes on the origin fibre is not excluded by the preceding line-equation consequence.

The result concerns only the conformal quartic direction (Y1), with its fixed zero quadratic and cubic jets. It does not classify arbitrary zero-curvature four-jets, arbitrary quartic tensors at $\rho=0$, or Ricci-flat nonflat metrics.
