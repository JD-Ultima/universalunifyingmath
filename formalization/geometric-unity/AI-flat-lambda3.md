# AI. The Lambda3 endpoints over a flat base

Over flat-base Z1, the two single $W_A$ Lambda3 branches admit a complete endpoint classification. The contraction/exterior-Gamma branch retains one known constant odd, zero even solution at one scale. The opposite branch has no solution at positive scales. The result concerns the full reconstructed $\Upsilon=0$ equation on that fixed observer metric; [AK](AK-lambda3-classification.md) subsequently classifies both single branches over arbitrary smooth positive bases and all positive Z1 parameters. The present note supplies the flat-model converse and supporting controls.

Code: [AI_flat_lambda3.py](../../computation/geometric-unity/AI_flat_lambda3.py) and [output](../../computation/geometric-unity/AI_flat_lambda3.out). Record: metrics/geometric-unity/AI-flat-lambda3.yaml. [Independent internal Review33](../../reports/reviews/REVIEW33-GU-AI-2026-11-02.md) accepts AI1–AI3 and reproduces all twelve controls in 1.7 seconds. This is internal review, not external peer review.

## Assumptions and flat observer geometry

Use the true dualized equation, compact real $(14,0)$, gauge algebra $\mathfrak u(128)$, fixed parallel real $W_A$, real fields and constant $\kappa\ne0$. On a connected open observer domain, keep
$$T_o=i(t\Phi^1+u\Phi^1\Gamma),\qquad
T_e=\Sigma(dt,du)+k,\qquad \det M\ne0.$$
The field remains in one single Lambda3 kernel branch, away from all other non-vector singular factors.

The base is a flat Riemannian four-metric, with reconstructed Z1 at $\lambda=1/5$ and $c_H,c_V>0$. [F7](F-zorro-reduction.md) gives the Einstein constant on the entire flat observer model, not only along the section:
$$\operatorname{Ric}^Y=-\frac{5}{4c_V}g_Y,\qquad
G^Y=\frac{15}{2c_V}g_Y,\qquad
\operatorname{tr}_{g_Y}G^Y=\frac{105}{c_V}>0.$$
The Einstein-tensor coefficient follows in dimension fourteen. The flat base-coordinate rescaling in F7 allows every positive $c_H$; no relation $c_H=2c_V$ is imposed here.

Put $a=t/\kappa$, $w=u/\kappa$. AF and AG give the only possible odd endpoints:
$$\begin{array}{c|c}
\Lambda^3:\mathsf c,\mathsf w\Gamma & (a,w)=(-3/112,0)\text{ or }(-1/400,0)\\
\Lambda^3:\mathsf w,\mathsf c\Gamma & (a,w)=(-3/272,0)\text{ or }(-1/240,0).
\end{array}$$
The odd field is constant on the connected branch domain, so its slaved vector part is zero. The extra even field may still vary.

## AI1. Three endpoint exclusions

At any of the three endpoints
$$a=-3/112,\qquad -3/272,\qquad -1/240,$$
the full equation is impossible whenever $\operatorname{tr}_{g_Y}G^Y\ge0$.

*Proof.* T2's exact compact residual identity for a full solution is
$$\frac{\langle T_o,E_v\rangle}{\kappa}=\frac12\|k\|^2\ge0.$$
For a constant odd endpoint with $u=0$ and zero slaved vector field, its left side is
$$a\,\operatorname{tr}_{g_Y}G^Y
+14\kappa^2(1+312a)a^2.$$
Each displayed $a$ is strictly below $-1/312$. Hence its first term is nonpositive under the trace hypothesis and the second is strictly negative. This contradicts the positive norm. Derivatives of $k$ are included in T2's orthogonality argument, rather than set to zero. ∎

This applies to flat Z1 because its trace is positive. It removes the first endpoint of AF and both endpoints of AG without imposing an extra condition on the remaining kernel form.

## AI2. The remaining endpoint

On the first branch at $a=-1/400,w=0$, a full solution over the flat base has
$$k=0,\qquad c_V\kappa^2=\frac{150000}{11}.$$

*Proof.* The kernel equations force its chart to have $r=0,s\ne0$, so
$$k=s\mathsf w(\omega)\Gamma.$$
AA7 makes the real three-form $\omega$ zero or decomposable at each point. AA4 puts this copy's derivative in grades nine and thirteen, including all weight derivatives; it has no grade-one contribution. The vector part is identically zero, so its derivatives and mixed terms vanish.

For a unit simple form the full quadratic grade-one tensor has eigenvalues $220s^2/3$ on its three-plane and $180s^2$ on the perpendicular eleven-plane. By $SO(14)$ equivariance and quadratic homogeneity, putting $N=|\omega|^2$ gives
$$[Q(k)]_1/i=180s^2NI-\frac{320}{3}s^2NP,$$
where $P$ is the orthogonal three-plane projector and $NP=0$ at a zero of the form.

The constant odd baseline has grade-one coefficient
$$\frac{15}{2c_V}+\kappa t+312t^2
=\frac{15}{2c_V}-\frac{11\kappa^2}{20000}.$$
Thus the full grade-one equation is
$$\left[\frac{15}{2c_V}-\frac{11\kappa^2}{20000}+180s^2N\right]I
=\frac{320}{3}s^2NP.$$
The identity has rank fourteen unless its scalar is zero, whereas the right side has rank at most three. It follows that both sides vanish. Since $s\ne0$ and the form norm is positive, $N=0$ and hence $k=0$. The remaining scalar equation gives the stated scale. No integrability or global three-frame choice is required. ∎

## AI3. Exact flat-base classification

Under the assumptions above:

| Single Lambda3 branch | Full solutions on a connected branch domain |
|---|---|
| $\mathsf c,\mathsf w\Gamma$ | At $c_V\kappa^2=150000/11$, exactly $T_o=-i\kappa\Phi^1/400$, $T_e=0$. At every other positive scale, none. |
| $\mathsf w,\mathsf c\Gamma$ | None at positive scales. |

*Proof of existence in the first row.* The constant odd field is parallel. With $T_e=0$, its even equation and all terms involving that even field vanish. The remaining odd equation is
$$i\left[\frac{15}{2c_V}+\kappa t+312t^2\right]\Phi^1=0,$$
which holds at $t=-\kappa/400$ and the stated scale. This is the existing T/R10 constant-background solution, now shown to be the only solution in this single branch over flat Z1. The two vector determinant factors are nonzero there, as are all other non-vector factors. The preceding endpoint arguments give necessity and exclude the second row. ∎

The classification is local on the reconstructed observer domain and works for either sign of $\kappa$. [AM](AM-global-flat-models.md) globally realizes its known solution on complete flat-space and torus observer models. Other field sectors, coincident kernels, indefinite pairings, other operators and Lambda5 fields are not classified here. AK supplies the broader positive-Z1 Lambda3 result. It does not solve a separate weighted base-metric variation.

## Reproduction

Run from the repository root:

    python -X utf8 -I computation/geometric-unity/AI_flat_lambda3.py

Twelve exact controls pass in 1.7 seconds with Python 3.12 and SymPy 1.14.0. They check the three endpoint residuals with arbitrary algebraic curvature diagonal and symbolic $\kappa$, coordinate and non-coordinate simple-form tensors, derivative grade controls, the complete remaining baseline, the positive scale, vector-block factors and regularity of every other non-vector block. Algebraic curvature controls are not asserted observer metrics. The flat Einstein formula is imported from F7, rather than recomputed by this helper.

The hand endpoint and rank arguments supply the classification. The helper shares T's accepted Clifford engine; no new physical evidence or second implementation is claimed.
