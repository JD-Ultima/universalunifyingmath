# AO. Global constant invariant fields

[Q15](Q-z1-remaining-members.md) already gives the constant invariant solutions over a flat positive base. This note realizes the whole family on [AM](AM-global-flat-models.md)'s complete observer and torus quotients, and replaces Q15's computed vector-singularity argument with an explicit hand factorization. Neither the metric nor the local solutions are claimed as new.

Use AM's oriented positive Z1 over standard affine $\mathbb R^4$ or a translation torus, with $\lambda=1/5$ and constant $c_H,c_V>0$. Keep the true reconstructed bosonic equation, compact $\mathfrak u(128)$ and real constant $\kappa\ne0$. For [J.0](J-slaved-sector.md)'s four-parameter self-adjoint family,
$$
\circledast=xW_A+yW'_A+zW_C\operatorname{Ad}_\Gamma+wW'_C\operatorname{Ad}_\Gamma,
$$
write its compact-real coefficients as
$$
(x,y,z,w)=-i(X,Y,Z,W),\qquad X,Y,Z,W\in\mathbb R,\qquad n=X^2+Y^2>0.
$$
Consider precisely the constant invariant ansatz
$$
T_o=i(S\Phi^1+R\Phi^1\Gamma),\qquad T_e=0,\qquad S,R\in\mathbb R.
$$
The condition $\det M\ne0$ is imposed only when specifying the regular slaving sector below. Extra even fields are not classified.

## AO1. The full constant-ansatz solutions

Put
$$
C=\frac{c_V\kappa^2}{n}.
$$
Real solutions exist exactly when $C\ge9360$. Their complete list is
$$
\zeta_\pm=\frac{-1\pm\sqrt{1-9360/C}}{624},\qquad
(S_\pm,R_\pm)=\frac{\kappa\zeta_\pm}{n}(X,Y).
$$
There is one double root at equality and two distinct coefficient pairs above it. Both signs of $\kappa$ and purely axial members $X=0$ are included. These counts do not assert gauge inequivalence.

*Proof.* Let $\sigma=XS+YR$ and $\delta=YS-XR$. Taking the difference and sum of [I6](I-elliptic-reduction.md)'s two real invariant channels with coefficients $(Y,-X)$ and $(X,Y)$ gives equivalently
$$
\delta(\kappa+208\sigma)=0,\qquad
312\sigma^2+\kappa\sigma+104\delta^2+\frac{15n}{2c_V}=0.
$$
The change of channels is invertible because $n>0$. The alternative $\sigma=-\kappa/208$ makes the second expression
$$
\frac{\kappa^2}{416}+104\delta^2+\frac{15n}{2c_V}>0.
$$
Thus $\delta=0$, so $(S,R)=\sigma(X,Y)/n$. The remaining quadratic gives the displayed roots. Twisted coefficients $Z,W$ do not enter these channels: their curvature contraction and their quadratic contribution on odd fields vanish, while the constant field has no derivative term. The two channels therefore exhaust the full equation in this ansatz. This uses I6 directly, without assuming that M1 is a symmetry of arbitrary full fields. ∎

## AO2. Exact vector singularities and regular counts

Define
$$
a=\frac{XW-YZ}{n},\qquad A=a^2,\qquad
L_\pm=384\pm32\sqrt{25+144A}.
$$
On the solution line $\delta=0$, J1 gives the full four-by-four determinant
$$
\det M=\kappa^2\left[(\kappa+224\sigma)(\kappa+544\sigma)
-147456A\sigma^2\right]
=\kappa^4(1+L_+\zeta)(1+L_-\zeta),\qquad \zeta=\sigma/\kappa.
$$
Exactly one scale can make a solution root singular:
$$
C_*(A)=\frac{15L_+^2}{2(L_+-312)}
=\frac{960(12+\sqrt{25+144A})^2}{9+4\sqrt{25+144A}}.
$$
It satisfies $C_*\ge9360$, with equality exactly at $A=125/576$.

*Proof.* Every allowed root has $-1/312<\zeta<0$. Since $L_-\le224$, its factor is strictly positive throughout this interval. Only $\zeta_*=-1/L_+$ can be singular. Set $\eta=9360/C$; the solution quadratic gives $\eta=-1248\zeta(1+312\zeta)$. Consequently
$$
\eta_*=\frac{1248(L_+-312)}{L_+^2}>0,\qquad
1-\eta_*=\frac{(L_+-624)^2}{L_+^2}.
$$
This proves the formula, uniqueness and equality case. At $A=0$, it recovers Q15's $C_*=277440/29$. For $A<125/576$, the normalized root $\zeta_-$ is singular at $C_*$; for $A>125/576$, it is $\zeta_+$. At equality the roots coincide. ∎

The regular-sector count is therefore:

| Parameters | Coefficient pairs with $\det M\ne0$ |
|---|---:|
| $C<9360$ | 0 |
| $C=9360$, $A\ne125/576$ | 1 |
| $C=9360$, $A=125/576$ | 0 |
| $C>9360$, $C=C_*$ | 1 |
| $C>9360$, $C\ne C_*$ | 2 |

At the threshold,
$$
\det M=\kappa^4\frac{125-576A}{1521}.
$$
Thus regular existence requires strict inequality only for $A=125/576$. A singular root still gives a smooth solution of the underlying constant-ansatz equation: $dt=du=T_e=0$ uses no inverse of $M$. It is excluded from the regular sector, not from the equation or the positive metric. This says nothing about additional even fields at that boundary.

## AO3. Global scalar and axial descent

Every AO1 field and its fixed operator are globally defined on AM's complete cover and on every translation-torus quotient, for all sixteen spin choices. The regularity condition is constant everywhere; descent adds no parameter or spin restriction.

*Proof.* AM2 supplies a translation-invariant oriented frame. Both $\Phi^1$ and the oriented Clifford volume $\Gamma$ are global and parallel. Central spin signs act trivially by conjugation, so $\Phi^1\Gamma$ descends just as $\Phi^1$ does. With positive Hermitian generators,
$$
\Gamma^\dagger=-\Gamma,\qquad
\Gamma\gamma_a=-\gamma_a\Gamma,\qquad
(\gamma_a\Gamma)^\dagger=\gamma_a\Gamma.
$$
Thus $T_o$ is anti-Hermitian. The Clifford, volume, Hodge and contraction ingredients of the fixed J.0 operator also descend. AO1's full residual vanishes in every chart and hence globally. Adding $T_o$ to the Levi-Civita spin connection gives a unitary connection; the identity gauge transformation is global. ∎

These fields act on the full 128-dimensional Dirac bundle: the odd lifts exchange its chiral halves. No parallel spinor or endomorphism of a single Weyl half is asserted. Orientation-reversing quotients, indefinite forms, $\kappa=0$, arbitrary fields, fermions, joint base variation, finite action and physical interpretation are outside this result. The $n=0$ case is the existing Q3 zero-field sector.

The proof uses I6, J1, Q15 and AM, with no new executable certificate or rerun. [Independent internal Review45](../../reports/reviews/REVIEW45-GU-AO-2026-11-03.md) accepts the argument and its stated scope; external peer review remains outstanding.
