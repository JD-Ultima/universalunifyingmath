# AN. No opposite Lambda5 constant-endpoint domains

The two endpoint cases left open by [AL](AL-opposite-lambda5-system.md) admit no smooth solution on a nonempty open positive-Z1 observer domain. Contracted Bianchi forces the five-form norm to be constant. The resulting scalar metric identity contradicts the radial geometry of a positive fibre. Together with [AH9](AH-lambda5-gram.md#ah9-the-first-pair-constant-endpoint-domains-are-empty), this excludes all four declared constant-odd Lambda5 endpoints.

This is a hand proof using the previously established AL identities and S4/F4 geometry. No computation or existing check was rerun for this result. [Internal Review44](../../reports/reviews/REVIEW44-GU-AN-2026-11-03.md) records the independent derivation and scope review.

## Setting

Retain AL's hypotheses: a smooth positive four-base; reconstructed Z1 with constant $c_H,c_V>0$ and $\lambda<1/4$; positive observer dimension fourteen; compact $\mathfrak u(128)$; the true dualized bosonic equation; fixed parallel real $W_A$; real fields; constant $\kappa\ne0$; and $\det M\ne0$. On an entire connected open observer domain, fix $u=0$ and either

| $a=t/\kappa$ | even field $k$ | $(q,\ell)$ | $b/\kappa^2$ |
|---|---|---|---|
| $-3/704$ | $\mathsf w(\omega)$ | $(112/3,320/3)$ | $87/61952$ |
| $-1/160$ | $\mathsf c(\omega)\Gamma$ | $(40,-32)$ | $19/3200$ |

Only the declared single non-vector kernel is present, away from coincident singular factors. Absorb every chart weight into the smooth real five-form $\omega$. The slaved vector field and its derivative vanish identically. Write $N=|\omega|^2$, $S_{ab}=\langle\iota_a\omega,\iota_b\omega\rangle$ and $J_a=\nabla_a\omega$, where $\nabla$ is the observer Levi-Civita connection. Clifford generators have positive square, and $\Gamma^2=-1$.

## AN1. Bianchi forces constant norm

Every full solution in this setting has $dN=0$.

*Proof.* For an arbitrary smooth five-form, metric compatibility and the exterior/contraction adjunction give
$$
(\operatorname{div}S)_b
=\langle\operatorname{Ctr}J,\iota_b\omega\rangle
+\tfrac12(dN)_b-\langle\iota_b\operatorname{Alt}J,\omega\rangle.\tag{AN1a}
$$
Indeed $\iota_b\operatorname{Alt}J=J_b-\sum_a e^a\wedge\iota_bJ_a$; pairing with $\omega$ identifies the second differentiated factor of $S$.

Use AL's notation $F_j=P_j(\Omega^2)$ and $Y_{d,a}=P_d(\Omega\gamma_a\Omega)$, with $\Omega$ the Clifford image of $\omega$. Reversion gives only grades $0,4,8$ in $\Omega^2$ and $1,5,9$ in the sandwich. The identity $\sum_a\gamma_a\Omega\gamma_a=-4\Omega$ therefore gives
$$
\operatorname{Alt}Y_5=0,\quad\operatorname{Ctr}Y_5=-4F_4,
\qquad \operatorname{Alt}(Y_9\Gamma)=-4F_8\Gamma,
\quad\operatorname{Ctr}(Y_9\Gamma)=0.
$$
AL3's determined jets, $\operatorname{Ctr}(\mathsf wF_4)=10F_4$ and $\operatorname{Alt}(\mathsf c(F_8\Gamma))=6F_8\Gamma$ now give the following table. Put $\xi=[\Omega^3]_1$ and $\eta=[\Omega^3\Gamma]_1$.

| copy | $\operatorname{Alt}J$ | $\operatorname{Ctr}J$ | $dN$ | $\operatorname{div}S$ |
|---|---|---|---|---|
| $\mathsf w\omega$ | $0$ | $-17F_4/9$ | $157\xi/18$ | $89\xi/36$ |
| $\mathsf c\omega\Gamma$ | $-22F_8\Gamma/15$ | $0$ | $8\eta/5$ | $-2\eta/3$ |

For the last column, the same scalar-trace adjunctions used in AL3 give $\langle F_4,\iota_b\omega\rangle=\xi_b$ and $\langle\iota_b(F_8\Gamma),\omega\rangle=-\eta_b$. Thus the coefficients are $157/36-17/9=89/36$ and $4/5-22/15=-2/3$, respectively. The axial sign is retained.

AL2a requires the actual raised Einstein tensor to satisfy
$$G=-(b+qN)I-\ell S.$$
Since $b,q,\ell$ are constant, contracted Bianchi imposes $q\,dN+\ell\operatorname{div}S=0$. In the two rows this becomes
$$\tfrac{1768}{3}\xi=0,\qquad \tfrac{256}{3}\eta=0,$$
respectively. Both coefficients are nonzero, so $dN=0$. ∎

## AN2. Constant norm contradicts the fibre trace

Neither endpoint has a full solution on a nonempty open observer domain.

*Proof.* AH6 gives $N>0$. By AN1 it is a constant $N_0$ on the connected domain. Fix a base point $x$ whose fibre meets the domain, and a nonempty relative-open fibre patch $\mathcal O$ in that intersection. Since $\operatorname{tr}S=5N$,
$$\tau(y):=\operatorname{tr}(G(y)+bI)=-(14q+5\ell)N_0\quad(y\in\mathcal O).$$
The coefficient $D=14q+5\ell$ is respectively $1056$ and $400$, so this constant is negative.

By [S4(a)](S-z1-off-section.md), $\tau(y)$ is rational in the fibre matrix $y$, with denominators powers of $\det y$. Clearing denominators in $\tau(y)+DN_0=0$ gives a polynomial vanishing on the nonempty real open set $\mathcal O$. It vanishes identically. Consequently the **metric equality** $\tau(y)=-DN_0$ holds on the whole positive fibre. No field, norm function or inequality is continued outside its original domain.

Put $c=4c_V(1-4\lambda)>0$. Along $y=\rho g(x)$, the exact F4/K trace formula recorded in AH7 is
$$\tau(\rho)=14b+\frac{30}{c}+\frac{135}{2c_V}
-\frac{6R_g(x)}{\rho c_H}.$$
Its limit as $\rho\to\infty$ is strictly positive, contradicting the continued negative constant. ∎

## AN3. The four constant endpoints

Combining AN2 with AH9, none of $a=-1/96,-3/832,-3/704,-1/160$ admits a full solution on an entire constant-odd, single-kernel domain under the standing positive-Z1 hypotheses. No curvature-trace sign, real analyticity of the base or fixed normalized Gram spectrum is assumed.

AL's nonzero-current first jets remain valid controls of its declared algebraic and derivative equations. They violate the additional Bianchi requirement for an actual metric solution. AL's reduction and conditional fixed-background uniqueness are also unchanged; they supply no endpoint existence theorem.

This does not exclude isolated zero-$u$ points of moving odd fields, nonzero-$u$ Lambda5 conics, coincident kernels, $\det M=0$, other operator families or indefinite real forms. It is a bounded exclusion within the reconstruction, with no physical or literature-first claim.
