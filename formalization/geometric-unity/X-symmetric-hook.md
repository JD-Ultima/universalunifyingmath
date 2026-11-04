# X. A pointwise obstruction for a symmetric hook family with one common axis

[AX5](AX-hook-curvature.md) later excludes the entire pure $H_2\Gamma$ branch over every smooth positive four-base with Z1 at $\lambda=1/5$, using the derivative equation, a flat-model curvature certificate and a fibre-density argument. This earlier note excludes a declared family by a separate pointwise argument. The calculation is symbolic in all thirteen eigenvalues; the extension to nondiagonal tensors uses equivariance, not parameter sampling.

Code: `computation/geometric-unity/X_symmetric_hook.py` → `.out`. It uses T's compact-real-form helpers and the earlier independent Clifford engine. Record: `metrics/geometric-unity/X-symmetric-hook.yaml`. Independently reviewed in `reports/reviews/REVIEW18-GU-X-2026-11-02.md`; all executable checks were reproduced. This is internal review, not external peer review.

## The declared family

Keep V3's assumptions: the true equation, signature (14,0), compact gauge real form $\mathfrak u(128)$, real member $W_A$, $\kappa\ne0$, and
$$t=3\kappa/136,\qquad u=0.$$
The odd field is $T_o=it\Phi^1$. Its even kernel is $H_2\Gamma$, and $\det M=(1313/17)\kappa^4>0$.

At a point choose a unit vector $r$ and a real symmetric traceless endomorphism $S$ of the thirteen-dimensional space $W=r^\perp$. Set
$$B_a=\sum_{b\in W}S_{ab}\gamma_b\gamma_r\quad(a\in W),\qquad B_r=0,
\qquad k_a=iB_a\Gamma.\tag{X1}$$
The contraction is $\sum_a\gamma^aB_a=(\operatorname{tr}S)\gamma_r=0$, so $B\in H_2$. The matrices $B_a\Gamma$ have grade 12 and are Hermitian; multiplication by $i$ makes $k$ a real anti-Hermitian field in $H_2\Gamma$. With the normalized positive pairing of T,
$$\|k\|^2=\operatorname{tr}(S^2)=:N\ge0.$$

For each fixed axis, this is a 90-dimensional linear family, since $\dim\operatorname{Sym}^2_0(\mathbb R^{13})=13\cdot14/2-1=90$. Allowing the axis to vary gives a union of these families, not every hook tensor. The unit axis may be chosen at each point. Neither the pointwise calculation nor the obstruction assumes that a smooth field in the family is parallel.

## The quadratic term

**Proposition X1.** For (X1), $Q(k)$ has grade 1 only. Written as an imaginary symmetric tensor in an orthonormal frame adapted to $W\oplus\mathbb Rr$,
$$\frac{Q(k)}i\ \longleftrightarrow\
\begin{pmatrix}4S^2-2N I_W&0\\0&-\tfrac23N\end{pmatrix}.\tag{X2}$$
In particular the grade-five equation of V3 vanishes identically on this family.

**Proof.** The stabilizer $SO(13)$ of $r$ diagonalizes the real symmetric $S$; the diagonalizing basis can be chosen oriented. Clifford multiplication, the Shiab operator and its cubic gradient are equivariant under its spin lift. It therefore suffices to calculate with $S=\operatorname{diag}(\lambda_0,\ldots,\lambda_{12})$.

First allow independent eigenvalues without imposing trace zero. Put $L=\sum_b\lambda_b$ and $N=\sum_b\lambda_b^2$. The exact differentiated cubic gives a diagonal grade-one output with
$$\frac{Q_a}{i\gamma_a}=4\!\sum_{b<c,\ b,c\ne a}\lambda_b\lambda_c
=2\big[(L-\lambda_a)^2-(N-\lambda_a^2)\big]\quad(a\in W),$$
$$\frac{Q_r}{i\gamma_r}=\frac43\sum_{b<c}\lambda_b\lambda_c
=\frac23(L^2-N).\tag{X3}$$
These are polynomial identities in thirteen indeterminates. The coefficient calculation has two types. A single diagonal basis field has zero quadratic term. For a pair of distinct diagonal basis fields $f_b,f_c$, the two terms in $3Q(f_b+f_c)$ have the following coefficients of $i\gamma_a$:

| output index $a$ | first gradient term | adjoint-gradient term |
|---|---|---|
| $b$ or $c$ | 0 | 0 |
| $W\setminus\{b,c\}$ | 4 | 8 |
| $r$ | 4 | 0 |

For example $[i\gamma_b\gamma_r\Gamma,i\gamma_c\gamma_r\Gamma]=-2\gamma_b\gamma_c$ gives the first column through the grade-3 anticommutator. The adjoint term uses the second commutator identity in V3's grade proof. Dividing the table by 3 gives the mixed quadratic coefficient $2\mathcal B(f_b,f_c)$, where $\mathcal B$ is Q's symmetric polarization. Permutation symmetry covers every pair, and quadratic polarization gives (X3). The script separately checks this coefficient table and all fourteen component polynomials, including the absence of other grades, with the signed Clifford product and differentiated cubic of the independent engine.

For trace-zero $S$, $L=0$, which gives (X2). Equivariance then covers every symmetric traceless $S$, not just its diagonal representatives. ∎

## The grade-one obstruction

**Proposition X2.** The tensor in (X2) is a scalar multiple of the fourteen-dimensional metric only when $S=0$, in which case it is zero.

**Proof.** Its trace on $W$ is $-22N$, whereas its $rr$ component is $-2N/3$. Isotropy would require
$$-22N=13(-2N/3),$$
so $N=0$. Positivity gives $S=0$. Equivalently,
$$\operatorname{tr}_{W}(Q/i)-13(Q/i)_{rr}=-\tfrac{40}3N$$
is a nonzero anisotropy trace whenever $S\ne0$. ∎

**Corollary X3 (flat-base Z1).** Let Z1 have a flat Riemannian base, $\lambda=1/5$ and $c_H,c_V>0$. No real full solution on the pure branch $t=3\kappa/136$, $u=0$ can have its even field in the family (X1), including its zero member.

**Proof.** F7 gives $G^Y=15g_Y/(2c_V)$. V3's grade-one equation requires the tensor (X2) to equal the scalar tensor
$$-\left(\frac{15}{2c_V}+\kappa t+312t^2\right)g_Y
=-\left(\frac{15}{2c_V}+\frac{201}{1156}\kappa^2\right)g_Y.$$
The scalar in parentheses is strictly positive for either sign of $\kappa$. X2 forces $S=0$, whose quadratic tensor is zero and cannot equal it. The derivative term ${*}^{-1}\circledast(d_0k)$ has grade 11 by V3; it cannot change the grade-one equation. ∎

This is a pointwise obstruction, so allowing $S$ or the axis $r$ to vary does not repair it. The separate grade-eleven differential equation has not been discarded or solved; it is simply unable to cancel this grade-one discrepancy.

## Evidence and limits

The script checks (X3) with thirteen independent real symbolic eigenvalues, then imposes their trace-zero relation and verifies the even-kernel condition and norm. A rational nondiagonal rotation is a separate convention check. It also checks all fourteen derivative symbols on that control, the positive flat-model scalar and the nonzero vector-block determinant.

The result does not classify arbitrary $H_2\Gamma$ fields, the other hook branches or coincident kernels. It does not address nonconstant motion on other singular curves, $\det M=0$, or indefinite real forms. AX5 supplies the later full-branch exclusion over arbitrary smooth positive bases at $\lambda=1/5$. The identities here remain useful algebraic information about the declared family; other parameter values and hook cases still require their full tensor and differential equations.
