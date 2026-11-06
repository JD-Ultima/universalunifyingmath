# BE. The affine-invariant actual H2 branch is empty

There is no affine-invariant full bosonic field on the pure actual $H_2$ branch $t=3\kappa/56,u=0$ over a flat positive four-base with reconstructed positive Z1, for any $\lambda<1/4$ and positive observer scales. The complete invariant ansatz has two real amplitudes. A grade-five coefficient eliminates the first; the remaining grade-one equations contradict positivity for the second.

The actual even field is a real bivector-valued $B$, with no factor $i$ or $\Gamma$. Evidence: [helper](../../computation/geometric-unity/BE_h2_affine.py), [output](../../computation/geometric-unity/BE_h2_affine.out) and [record](../../metrics/geometric-unity/BE-affine-h2.yaml). Internal [Review66](../../reports/reviews/REVIEW66-GU-BE-2026-11-05.md) and [Review67](../../reports/reviews/REVIEW67-GU-BE-2026-11-05.md) accept BE1–BE3 and disclose prior contributions; external peer review remains outstanding.

## Setting and positive-scale normalization

Retain the true reconstructed dualized bosonic equation, fixed parallel real $W_A$, positive observer signature $(14,0)$, compact gauge real form $\mathfrak u(128)$, real smooth fields and constant real $\kappa\ne0$. The four-base is flat and positive; $c_H,c_V>0$ and $\lambda<1/4$ are constant. The even tensor is the restriction of a tensor invariant under
$$G=GL^+(4,\mathbb R)\ltimes\mathbb R^4,\qquad Y=G/SO(4),$$
including its spin/Clifford identification. This invariant ansatz is an essential extra assumption.

Fix $T_o=i(3\kappa/56)\Phi^1$. R2/V leave exactly $T_e=B\in H_2$, with a regular vector block and no other even kernel. Real Clifford bivectors are anti-Hermitian in this compact convention.

Use [AZ1](AZ-first-hook-classification.md)'s affine model and write
$$\mu=\frac1{1-4\lambda}>0.$$
This covers the entire positive range, including negative $\lambda$. Only $\mu=5$ is the $SL(5)/SO(5)$ metric.

The dilation $F_s(x,y)=(sx,y)$, $s=\sqrt{c_H/(2c_V)}$, gives
$$g(c_H,c_V)=c_VF_s^*g(2,1)$$
and conjugates the affine action, preserving the ansatz. For a fixed three-covariant $B_2$ under $g\mapsto hg$, the metric-relative derivative scales by $h^{-2}$ and the metric-relative grade-one quadratic by $h^{-3}$. These powers follow from five, respectively seven, Clifford covectors, their inverse-metric contractions and the external Gram factor. Grade five of the quadratic scales by $h^0$.

After this dilation, set
$$\widehat B=B/c_V,\qquad \rho=c_V\kappa^2>0.$$
Multiplying the grade-one equation by $c_V$ gives exactly the normalized metric's equation for $\widehat B$, with odd-field shift $93\rho/98$. Its grade-five zero condition is unchanged by this nonzero rescaling. Thus it suffices to work at $c_H=2,c_V=1$ with arbitrary $\mu,\rho>0$. Below write $\widehat B=xP_0+yP_1$.

## BE1. Complete invariant space and retained derivative

Public BC functions reconstruct the metric-independent $SO(4)$ tensor invariants. Exact finite-subgroup reduction and six infinitesimal equations leave three vector-valued two-forms. Their coefficients have no external index $13$, the scalar vertical direction. The helper checks this for **all three raw invariants**.

In BC's basis the only changing Gram entry is $g_{13}=16/\mu$; all others are fixed. The contraction weight $g_{13}^{-1}$ therefore never enters on this entire invariant space. Alternation is metric independent, so both hook-constraint maps are unchanged for every $\mu>0$. Their rank-one reduction leaves the same complete two-dimensional hook space, including at special parameter values. Persistence of two example tensors alone would not establish this completeness. Both symbolic hook conditions are checked again on the resulting basis, whose coefficients are printed.

The covariant derivative is minus the Nomizu action on the coefficient covector and both form slots. In the general AZ1 model, only the scalar coordinate $(\Lambda_{H_r})_{13,r}=\mu/8$ changes from the $\mu=5$ matrices. The helper constructs its own parameter-dependent maps, checks metric compatibility and does not mutate BC globals or cached actions.

With $H_{bc}=\nabla_bB_c-\nabla_cB_b$ and the alternating Clifford triple $A_{abc}$, the orthonormal formula is
$$\mathcal D(B)_a/i=\sum_{b<c}\{A_{abc},H_{bc}\}.$$
In the non-unit covector basis this has the external factor $g_a$, with no additional inverse $b,c$ weights outside the Clifford product. Both derivative grades are evaluated first. Its apparent grade five vanishes by the differentiated alternating-hook condition:
$$
2e^a\wedge\sum_{b,c}e^b\wedge e^c\wedge\nabla_bB_c
=2e^a\wedge\sum_b e^b\wedge\nabla_b(\operatorname{Alt}B)=0.
$$
This normal-frame identity is covariant. The helper verifies it on both complete invariant basis directions. It does not set the remaining derivative to zero or impose covariant closure.

The real-field quadratic phase is also retained: BC evaluates $Q(iC)/i$, so $Q(B)/i$ is its negative. Complete outputs on $P_0,P_1,P_0+P_1$ determine the entire quadratic by polarization. Only grades one and five occur.

## BE2. Exact coefficient equations

The three irreducible tangent summands are horizontal $H$ of dimension four, tracefree symmetric $S$ of dimension nine and scalar vertical $I$ of dimension one. The normalized metric-relative derivative entries are
$$
(d_H,d_S,d_I)=
\left(\frac{-9x+27\mu y}{16},\ \frac x4-\mu y,\ \frac{9\mu y}{4}\right),
$$
and the quadratic entries are
$$
(q_H,q_S,q_I)=
\left(-\frac{3x^2+51\mu y^2}{16},\
-\frac{11x^2+109\mu y^2}{24},\
-\frac{3x^2+117\mu y^2}{8}\right).
\tag{BE2a}
$$
Every grade-one diagonal and off-diagonal coefficient is checked symbolically. Direct inverse-Gram norm and trace checks give
$$
N=\|\widehat B\|^2=\frac9{16}(x^2+13\mu y^2),\qquad
\operatorname{tr}d=0,\qquad\operatorname{tr}q=-\frac{28}{3}N.
$$

From AZ1 the normalized Ricci eigenvalues on $H,S,I$ are $-\mu/4,-5/4,-\mu/4$. Therefore the Einstein tensor eigenvalues are
$$
(G_H,G_S,G_I)=
\left(\frac{3(\mu+15)}8,\ \frac{5(\mu+7)}8,\ \frac{3(\mu+15)}8\right).
$$
The generic metric is not assumed Einstein. The complete grade-one equations are
$$G_j+\frac{93}{98}\rho+d_j+q_j=0,\qquad j=H,S,I,$$
since $\kappa t+312t^2=93\kappa^2/98$.

Neither curvature nor the constant invariant-odd baseline contributes grade five. Its derivative contribution has already vanished. A genuine grade-five coefficient at external $0$ and form indices $(0,1,2,4,5)$ is
$$-\frac43x^2.$$
The full equation thus forces $x=0$. One necessary output suffices; no classification of the whole grade-five zero locus is assumed.

## BE3. Positive residual excludes both remaining cases

Since $G_H=G_I$, their difference after $x=0$ is
$$\frac{3\mu}{16}y(61y-3)=0.$$
Positive $\mu$ gives exactly $y=0$ or $y=3/61$. The horizontal residual is
$$
\mathcal R_H=\frac{3(\mu+15)}8+\frac{93}{98}\rho
+\frac{\mu}{16}(27y-51y^2).
$$
At $y=0$ its first two terms are strictly positive. At $y=3/61$ the last term is
$$\frac{2241}{29768}\mu>0.$$
Thus neither case solves the equation. This proves nonexistence for every $\mu,\rho>0$, hence for every positive Z1 parameter and positive scale pair, with either sign of nonzero real $\kappa$. No root sampling, general solver or inequality continuation is used.

## Reproduction and boundary

Run

    python -X utf8 -I computation/geometric-unity/BE_h2_affine.py

The author run takes 0.464 seconds with pinned Python 3.12/SymPy 1.14.0. It reconstructs the complete invariant space, proves the parameter-independent hook constraints, computes every derivative slot and complete quadratic output, and checks the exact symbolic tables, norm, traces and both positive residuals. Coefficients are asserted free of floating-point numbers; $\mu=5$ recovers the earlier normalized tables.

Three new complete-output comparisons with the separately authored T/xcheck implementation verify two non-unit derivative jets and the real-field quadratic phase, including nonzero grade-five controls. These are independent convention checks; the complete parameter-dependent homogeneous residual uses BC's tensor/Clifford library and is not recomputed by a second engine. No old helper main or private input is used.

The result is restricted to affine-invariant actual $H_2$ fields over a flat positive base. Arbitrary smooth $H_2$ fields, curved bases, other hook copies and coincident kernels retain their separate scope. No field construction, literature-priority claim or physical conclusion is asserted.

Prior component work: [K2375](https://github.com/disruptionjoe/gu-formalization/blob/20fc3665c6690b68544b0b8a26fec29295416384/explorations/conditional-build/k2375-axial-hook-stationarity-2026-10-31.md#L126) contains related H2 contraction identities in a different real form and selected-action convention. The [comparison](../../reports/RELATED-WORK.md) records the overlap; its statement does not supply the all-positive-parameter affine exclusion proved here.
