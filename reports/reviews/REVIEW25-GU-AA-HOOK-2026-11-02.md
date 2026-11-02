# Internal review of the Lambda-three hook obstruction

November 2, 2026. Review of AA6 in `formalization/geometric-unity/AA-conic-traces.md`, `AA_lambda3_hook.py` and its output, together with the separately constructed scratch cross-check. AA6 is valid as a necessary pointwise condition on the stated single-kernel branch. Review20 and Review23 retain their original scopes. This is internal review, not external peer review.

## Universal quadratic patterns

The field is $k=r\mathsf c(\omega)+s\mathsf w(\omega)\Gamma$, for a real three-form and real kernel-chart weights. Both input grades, two and ten, are anti-Hermitian without an additional phase. The real $W_A$ cubic normalization and normalized trace are those of T.

I checked the metric-pairing count for the full grade-five output and the converted grade-nine output. Each quadratic map has two alternating three-tensor inputs and output slots consisting of one vector index plus five alternating indices: twelve slots in total. Clifford normal ordering supplies metric contractions after paired volume elements are simplified. Any remaining orientation-dependent term would require a fourteen-index epsilon tensor; contracting its unused indices with the symmetric metric gives zero. It cannot introduce another rank-twelve pattern.

Pairing two indices within an input three-form or within the five alternating output slots vanishes. With no cross-input contraction, the free output vector index must pair with an input index, leaving exactly

$$Z_a=(\iota_a\omega)\wedge\omega.$$

With one cross-input contraction, the free vector index must pair with an alternating output index, giving $\mathsf w(F_4)$, where $F_4=\sum_a(\iota_a\omega)\wedge(\iota_a\omega)$. Two cross-input contractions force an alternating output-output trace and vanish. Thus the two displayed patterns exhaust the actual Clifford contractions.

The nondecomposable calibrations distinguish these patterns. The disjoint-triple form has $F_4=0$ and nonzero $Z$; a shared-index calibration has nonzero $F_4$ and nonzero hook part. Exact full-output calibration with symbolic $r,s$ therefore fixes the coefficient polynomials universally. Decomposable controls alone would not suffice, since both patterns vanish there.

## Signs and projection normalization

The checked full identities are

$$Q_5=\frac{8i}{3}[(43s^2-3r^2)Z+(4r^2-49s^2)\mathsf w(F_4)],$$

$$\widetilde Q_9=Q_9\Gamma
=\frac{32i}{3}rs[5Z-7\mathsf w(F_4)].$$

The conversion is right multiplication by $\Gamma$. Since $\Gamma^2=-1$, recovery is $Q_9=-\widetilde Q_9\Gamma$; replacing the defined conversion by its inverse would reverse the second formula's sign. Both computations use the compact real phases and the true differentiated cubic, including its factor of one third.

On one-form-valued five-forms, $\operatorname{Alt}\mathsf c=6$ and $\operatorname{Ctr}\mathsf w=10$. The cross compositions vanish, so

$$\mathcal H_5=I-\mathsf c\operatorname{Alt}/6-\mathsf w\operatorname{Ctr}/10$$

has the stated normalization. Here $\operatorname{Alt}Z=3\omega\wedge\omega=0$ and $\operatorname{Ctr}Z=F_4$. Hence

$$\mathcal H_5Z=Z-\mathsf w(F_4)/10,$$

and the two projected coefficients are precisely $8i(43s^2-3r^2)/3$ and $160irs/3$.

## Why the full equation supplies this condition

The slaved baseline need not solve its own odd equation. Its invariant terms have grades one and thirteen, and its curvature term has grade one. The vector derivative and vector quadratic terms have types contained in $V\otimes V$, which contains no $\mathcal H_5$. This remains true for multiple vector copies.

The vector/kernel mixed term has input type $V\otimes\Lambda^3=\Lambda^4\oplus\Lambda^2\oplus\mathcal H_3$, so it also has no $\mathcal H_5$ component in the grade-five target or its converted grade-nine partner. AA4 places those components of the full covariant derivative in exterior lifts of four-forms. Their hook projection vanishes, including all derivatives of the weights and of $\omega$.

The full odd equation therefore requires

$$ (43s^2-3r^2)\mathcal H_5Z=0,\qquad rs\,\mathcal H_5Z=0.$$

If the hook tensor were nonzero, $rs=0$ would force one weight to vanish, and the other equation would force the second to vanish too. A valid kernel chart excludes $(r,s)=(0,0)$. Thus $\mathcal H_5Z=0$ is necessary at every point.

## Charts, example and limits

The result uses the single $\Lambda^3:\mathsf c,\mathsf w\Gamma$ kernel, away from every other non-vector singular factor, with $\kappa\ne0$ and $\det M\ne0$. Either smooth nonzero chart from AA4 can be used at a zero-$u$ endpoint. No division by $u$ or representation of the weights as smooth functions of $a$ is required. On chart overlaps the three-form rescales by a nonzero scalar; the hook condition is homogeneous and is unchanged. The background remains compact $(14,0)$ with fixed parallel real $W_A$.

For $\omega=e^{012}+e^{034}$, direct exterior calculation gives $F_4=2e^{1234}$ and $Z_5=0$, while $\mathsf w(F_4)_5=2e^{12345}$. Therefore

$$ (\mathcal H_5Z)_5=-e^{12345}/5\ne0.$$

Every nonzero chart-weight realization of this form is excluded on this branch, regardless of the derivative jet. No curvature-trace sign assumption is used.

The zero locus has not been classified. The condition is not asserted to characterize decomposability or to suffice for a full solution. Other conics, sums of coincident kernels and the remaining odd tensor/PDE equations are not settled. In particular the derivative term is retained; only the specified hook projections vanish universally.

## Reproduction

`python -X utf8 -I computation/geometric-unity/AA_lambda3_hook.py` passed all **26 checks** in **6.3 seconds** with Python 3.12 and SymPy 1.14.0. Deterministic output matched the recorded file. A separate workspace scratch calculation, `Lambda3_independent_hook_scratch.py`, also passed in 8.1 seconds, reconstructing wedge products, contractions and the projector locally. It imports neither the new helper nor its candidate implementation. The repository's reproduction does not depend on that scratch file. Both checks use the accepted T Clifford engine; neither is a second Clifford implementation.

The checks include symbolic weights, nondecomposable calibrations, full fourteen-coordinate support, projector traces, the explicit excluded component, a genuine four-copy vector control and variable-weight derivatives. The artificial curvature coefficients in a baseline control are algebraic data, not a constructed observer metric. The contraction and representation arguments above supply universality beyond these finite controls.
