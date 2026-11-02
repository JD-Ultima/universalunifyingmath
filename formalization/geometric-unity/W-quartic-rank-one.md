# W. An obstruction with arbitrary quartic base terms

The three-row obstruction of U survives every symmetric quartic term when the quadratic jet is fixed and the cubic jet is zero. Arbitrary higher terms are allowed. The conclusion concerns a declared class of metric four-jets; it does not classify all metrics whose Ricci tensor has rank one.

Code: `computation/geometric-unity/W_quartic_rank_one.py`; output: `W_quartic_rank_one.out`; record: `metrics/geometric-unity/W-quartic-rank-one.yaml`. The code reuses U's independently reviewed relaxed-system assembly and the earlier referee curvature engine `xcheck_round10.py`. Its symbolic check includes all 350 quartic coefficients. Exact engine controls check the conventions; the universality follows from the jet and gradient arguments below. It was independently reviewed in `reports/reviews/REVIEW15-GU-W-2026-11-02.md`; all 33 assertions were reproduced. This is internal review, not external peer review.

## The declared four-jets

Use U's notation. For a unit vector $p$, $\rho\ne0$, and any symmetric homogeneous quartic matrix $Q_4(x)$, let
$$g_{ij}(x)=\delta_{ij}-\frac{\rho}{3}(R_p)_{iajb}x^ax^b+Q_{4;ij}(x)+R_{5;ij}(x),\qquad
R_p=\tfrac12(p\otimes p-\delta/6)\owedge\delta.\tag{W1}$$
Here $R_5$ is any smooth symmetric remainder whose derivatives through order four vanish at $0$. The cubic coefficients are zero. The metric is positive definite on some neighbourhood of $0$, without a uniform neighbourhood over all coefficients. At $0$, $\mathrm{Ric}=\rho p\otimes p$ and the Weyl tensor vanishes. A symmetric quartic matrix in four variables has $10\binom74=350$ independent real coefficients.

Use Z1 at $\lambda=1/5$ with $c_H,c_V>0$. Rotations take $p$ to $e_0$ and preserve the class of arbitrary $Q_4$. Under $F_s(x)=sx$, the parameters transform as
$$(\rho,Q_4)\longmapsto(\rho s^2,s^4Q_4).$$
Thus U's lifted dilation and observer-metric homothety reduce the problem to
$$\rho=\sigma\in\{1,-1\},\qquad c_H=2,\quad c_V=3,\quad Q_4\text{ arbitrary}.\tag{W2}$$
The normalized quartic tensor is $Q_4/\rho^2$. U1–U2 prove consistency and rank invariance for the relaxed system, including its free covector and vector unknowns.

## Jet order and a gradient cancellation

**Lemma W1 (jet order).** At $(0,y)$, for every positive-definite $y$, the relaxed coefficient matrix $A$ is independent of $Q_4$ and $R_5$. Its right-hand side $b$ is affine in the quartic coefficients and independent of $R_5$.

*Proof.* At a Euclidean normal base point, the leading connection variation is
$$\delta\Gamma_\mu{}^a{}_i=\tfrac12\big(\partial_\mu Q_{4;ai}+\partial_iQ_{4;a\mu}-\partial_aQ_{4;\mu i}\big)+O(|x|^4).\tag{W3}$$
Its degree-three part is exactly linear in $Q_4$. The remainder $R_5$ first affects degree four of the connection. For Z1,
$$h_{HH}=c_Hy+c_VG_y(\Psi,\Psi),\qquad h_{HV}=-c_VG_y(\Psi,\cdot),\qquad
h_{VV}=c_VG_y,\qquad \Psi_\mu=\Gamma_\mu^\top y+y\Gamma_\mu.$$
Since $\Psi(0)=0$, the quartic variation of $h_{HV}$ begins at base order three, that of $h_{HH}$ at order four, and $h_{VV}$ is fixed. Fibre derivatives cannot lower this base order. Consequently the observer 2-jet is fixed and its 3-jet is affine in $Q_4$; $R_5$ affects neither. Riemann and Ricci depend only on the 2-jet, while $\nabla\mathrm{Ric}$ is linear in the third metric derivatives when the lower jets are fixed. This proves the asserted dependence. ∎

**Lemma W2 (the three diagonal-vertical rows).** At a diagonal fibre point $y=\operatorname{diag}(d_0,d_1,d_2,d_3)$, the RHS rows $(a,b,d)=(i,A,i)$ are unchanged by every $Q_4$, where $i=0,1,2$ and $A$ is the vertical coordinate $y_{33}$.

*Proof.* The fibre coordinate basis vector for $A$ is $E_{33}$. Let $\mathcal H_j$ denote the cubic variation of $h_{jA}$. Directly from the DeWitt pairing,
$$\mathcal H_j=-\frac{2c_V}{d_3}\left(\delta\Gamma_j{}^3{}_3-\lambda\sum_a\delta\Gamma_j{}^a{}_a\right).
$$
Equation (W3) gives $\delta\Gamma_j{}^a{}_a=\tfrac12\partial_jQ_{4;aa}$, separately for each $a$. Hence
$$\mathcal H_j=\partial_jF,\qquad
F=-\frac{c_V}{d_3}\big(Q_{4;33}-\lambda\operatorname{tr}Q_4\big).\tag{W4}$$
Thus this diagonal-vertical metric variation is a gradient for every quartic tensor, including off-diagonal coefficients.

Because the lower observer jets are fixed, only the principal third-derivative part of the covariant derivative of Ricci varies:
$$\delta(\nabla_e\mathrm{Ric})_{ab}=\tfrac12h^{cd}\delta
\big(h_{bd,cae}+h_{ad,cbe}-h_{ab,cde}-h_{cd,abe}\big).$$
For $e=a=i$, $b=A$, terms containing vertical derivatives vanish at $0$ by the base orders in W1. The remaining expression is
$$\delta(\nabla_i\mathrm{Ric})_{iA}
=\tfrac12\sum_{c=0}^3h^{cc}\big(\partial_c\partial_i^2\mathcal H_c-
\partial_c^2\partial_i\mathcal H_i\big)=0,\tag{W5}$$
where (W4) and commutation of fourth derivatives give the last equality. Also $\delta(\nabla_A\mathrm{Ric})_{ii}=0$: the Ricci tensor and the observer connection on the entire fibre are already independent of $Q_4$. Raising the fixed Ricci output index therefore gives $\delta b_{(i,A),i}=0$. ∎

The symbolic calculation constructs all 350 quartic coefficients, verifies (W4), and checks the vanishing principal expressions in (W5). Separate mixed-monomial engine controls change other RHS rows but preserve the coefficient matrix and the three selected rows. No degree or support bound is inferred from samples.

## The exact three-row certificate

Use the canonical witness $y_*=\operatorname{diag}(2,3,4,5)$. Rows are ordered as in U: $(a,b,d)$ lexicographically, $0\le a<b<14$, $0\le d<14$. Vertical coordinate 13 is $y_{33}$, so
$$I=(168,337,492)$$
corresponds to $(0,13,0)$, $(1,13,1)$ and $(2,13,2)$. W1–W2 make both their coefficient rows and RHS independent of the arbitrary quartic tensor.

**Proposition W3 (exact certificates).** For the two canonical signs, the following combinations annihilate all 42 coefficient columns but not the RHS:

| $\rho$ | weights on rows $(168,337,492)$ | coefficient combination | RHS combination |
|---|---|---:|---:|
| $+1$ | $(-445,41538,-41093)$ | $0$ | $3933/800$ |
| $-1$ | $(445,44862,-45307)$ | $0$ | $5667/800$ |

*Certificate.* Evaluate the three rows at $Q_4=0$ over $\mathbb Q$ using the earlier independent referee curvature engine and U's assembly. Every coefficient of the stated combination is zero and the RHS is the displayed nonzero rational. W1–W2 preserve this contradiction for every $Q_4$. The coefficient rank is 42 by U and remains fixed; the augmented rank is consequently 43. ∎

The scalar equations are small. In the coordinate order $(x^0,\ldots,x^3,y_{00},y_{01},\ldots,y_{33})$, put
$$L(n)=-\frac3{50}n^4-\frac1{25}n^8-\frac3{100}n^{11}+\frac{12}{125}n^{13}.$$
All omitted coefficients vanish. The equations are $L(n)+a\xi_{13}=r$ with

| row | $a$ for $\rho=+1$ | $a$ for $\rho=-1$ | $r$ |
|---|---:|---:|---:|
| 168 | $-7849/43200$ | $-29449/43200$ | $1/500$ |
| 337 | $-24163/57600$ | $-24163/57600$ | $1/4500$ |
| 492 | $-36467/86400$ | $-36467/86400$ | $1/12000$ |

For either sign the integer weights sum to zero and annihilate the $a$ column, while their $r$ combination is nonzero. Thus the obstruction can be inspected as three equations in two unknown scalar combinations.

**Theorem W4 (the declared four-jet family).** For every metric (W1), with $\rho\ne0$, arbitrary $Q_4,R_5$ and $c_H,c_V>0$, the relaxed system is consistent on no nonempty open subset of the fibre over $0$. In the slaved vector sector of S.0, the twisted cone admits no $\Upsilon=0$ solution on an open set meeting that fibre.

*Proof.* W3 excludes both signed canonical representatives for every quartic tensor. The normalization (W2) transfers the obstruction to all nonzero amplitudes, directions and positive observer-metric scales. Augmented rank 43 at the witness gives a nonzero augmented minor; S4 then excludes a fibre-open consistent patch.

For the full-field conclusion, fix the metric, including its quartic coefficients and remainder. Coefficient rank 42 and the three-row contradiction give augmented rank 43 at the origin. Choose an actual nonzero augmented 43-minor there and its fixed positive coordinate fibre matrix. Continuity of that fixed minor gives a base neighbourhood on which the same fibre-open obstruction holds. Shrink this neighbourhood so that the base curvature stays nonzero. The integer row weights certify the origin; it is the fixed minor that supplies the neighbourhood.

U5's argument now applies. Einstein pieces are impossible by S5, and the dense variable patch in the Q5 reduction supplies an open set with $H(t)\ne0$, where the line equation implies relaxed consistency. This is a contradiction. All sector and slaving hypotheses of S.0 remain in force. No uniform base neighbourhood over quartic coefficients or remainders is asserted. ∎

## A one-parameter control

The earlier direction $Q_4=\beta(x^0)^4I$ provides a useful check that the engine retains fourth base derivatives. At a diagonal fibre point its only nonzero variations of $\nabla\mathrm{Ric}$ are, with Ricci-slot symmetry,
$$\partial_\beta(\nabla_0\mathrm{Ric})_{i,A_{0i}}=
\frac{12c_V(d_i-d_0)}{c_Hd_0^2d_i},\qquad i=1,2,3.$$
At the canonical witness the six nonzero relaxed RHS slopes are

| row | 5 | 20 | 35 | 57 | 72 | 87 |
|---|---:|---:|---:|---:|---:|---:|
| $\partial_\beta b$ | $3/2$ | $3$ | $9/2$ | $1/4$ | $9/32$ | $27/100$ |

The computed full slope matches this formula for both signs, and $\beta=-2$ is a separate affine control. Nonunit controls at $(\rho,\beta,c_H,c_V)=(4,48/5,5,7)$ and $(-9,-162/7,3/2,5/2)$ verify the lifted tensor identity with normalized coefficients $\beta/\rho^2=3/5$ and $-2/7$. These controls supplement the arbitrary-quartic proof.

## Remaining boundary

The quadratic curvature tensor is fixed, so the Weyl tensor at the origin is zero; the cubic metric coefficients are also fixed to zero. W does not classify arbitrary cubic jets, arbitrary Weyl components, or the full rank-$\mathrm{Ric}\le1$ frontier.

The case $\rho=0$ is outside this result. A nonzero quartic tensor can produce a curved metric near the origin even though its curvature vanishes at the origin; this case must not be identified with U's flat exception. A next bounded question is whether the certificate survives a declared cubic metric direction, which can change the relaxed coefficient matrix itself.
