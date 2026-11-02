# AK. The two single-Lambda3 branches on positive Z1

Both single $W_A$ Lambda3 branches admit a complete classification over smooth positive four-bases. On the contraction/exterior-Gamma branch, a full solution requires a flat base, $\lambda=1/5$, $c_V\kappa^2=150000/11$ and the known constant odd, zero even field. The opposite branch has no solution. Flatness and the value of $\lambda$ are conclusions, rather than initial restrictions.

The result concerns the specified reconstructed field equation and two single kernels. It does not classify arbitrary fields or establish a physical theory.

## Assumptions and endpoint tensors

Use a smooth Riemannian four-base, reconstructed Z1, constant $c_H,c_V>0$ and $\lambda<1/4$. The observer metric is positive of dimension fourteen. Keep the true dualized equation, compact gauge real form $\mathfrak u(128)$, fixed parallel real $W_A$, real fields, constant $\kappa\ne0$ and
$$T_o=i(t\Phi^1+u\Phi^1\Gamma),\qquad T_e=\Sigma(dt,du)+k,\qquad \det M\ne0.$$
Let a connected open observer domain remain in one single Lambda3 kernel, away from every other non-vector singular factor. The bundle, Clifford and parallel-operator data are assumed defined there; this is a local classification, with no separate global extension or topology claim.

[AF](AF-symmetric-gamma.md) and [AG](AG-other-lambda3.md) force $u=0$ and one constant endpoint $a=t/\kappa$. These algebraic reductions require a positive fourteen-dimensional observer metric; they do not impose a special $\lambda$ or base curvature. The slaved vector field therefore vanishes identically. Absorb the nonzero local kernel chart weight into a real three-form $\omega$. AA6–AA7 for the first branch, and AG2 for the opposite branch, make this form zero or simple at every point.

Write $G$ for the raised observer Einstein tensor, $z=|\omega|^2$ and $P$ for the orthogonal support three-plane at nonzero $\omega$. Use $zP=0$ at its zeros. The necessary symmetric grade-one equation is
$$H:=G+bI=-CzI-(A-C)zP,\qquad b=\kappa^2a(1+312a).\tag{1}$$
The four compact endpoint copies give:

| Branch | $a$ | Actual even copy | $b/\kappa^2$ | $A$ | $C$ |
|---|---:|---|---:|---:|---:|
| $\mathsf c,\mathsf w\Gamma$ | $-3/112$ | $\mathsf c(\omega)$ | $309/1568$ | $4$ | $4$ |
| $\mathsf c,\mathsf w\Gamma$ | $-1/400$ | $\mathsf w(\omega)\Gamma$ | $-11/20000$ | $220/3$ | $180$ |
| $\mathsf w,\mathsf c\Gamma$ | $-3/272$ | $i\mathsf c(\omega)\Gamma$ | $249/9248$ | $4/3$ | $12$ |
| $\mathsf w,\mathsf c\Gamma$ | $-1/240$ | $i\mathsf w(\omega)$ | $1/800$ | $220$ | $60$ |

Here $[Q(k)]_1/i=CzI+(A-C)zP$. Complete unit-simple Clifford calibrations fix the table. $SO(14)$ equivariance and quadratic homogeneity give the same tensor for every real simple form; zero forms require no choice of plane. These are pointwise identities, independent of base curvature and $\lambda$.

All covariant derivatives are retained. For the contraction copy, [AA4](AA-conic-traces.md) gives
$$\mathcal D(\mathsf c(\omega))=i\big[6\mathsf w(d\omega)-2\mathsf c(\delta_\nabla\omega)\big].\tag{2}$$
Its grade-one term is the antisymmetric two-tensor of $\delta_\nabla\omega$, so its symmetric part is zero. Its other term has grade five. The exterior-Gamma copy has derivative grades nine and thirteen. AG2 gives grades eleven for $i\mathsf c(\omega)\Gamma$ and three/seven for $i\mathsf w(\omega)$. Thus no derivative changes (1), including derivatives of the absorbed weight. The remaining derivative equations have not been assumed solved.

## AK1. The radial multiplicity criterion

Let $A\ne C$ be real constants, $D=3A+11C\ne0$, and $b$ a real constant. If $H=G+bI$ satisfies
$$\big(DH-A\tau I\big)\big(DH-C\tau I\big)=0,
\qquad \tau=\operatorname{tr}H,\tag{3}$$
on a nonempty open observer set, then
$$\lambda=\frac15,\qquad b=-\frac{15}{2c_V}.$$

*Proof.* Fix a base point in the open set's projection. By [S4(a)](S-z1-off-section.md), every coordinate entry of (3) is rational in the positive fibre matrix $y$, with denominators powers of $\det y$. Clearing denominators extends this metric identity from the nonempty fibre-open patch to the entire positive fibre. No field, inequality or projector is continued.

Use fixed $g$-normal coordinates and take $y=qg(x)$, $q\to\infty$. The parallel metric $qg$ has the same connection as $g$, so [F4/K](K-z1-ricci-general-n.md) applies along this ray. The base Ricci correction is proportional to $1/q$ and tends to zero. The limiting raised observer Ricci tensor has eigenvalues
$$\alpha_0=-\frac1c\quad\text{with multiplicity }5,
\qquad \beta_0=-\frac5{4c_V}\quad\text{with multiplicity }9,
\qquad c=4c_V(1-4\lambda)>0.$$
The first space is horizontal four-space plus the Euler line; the second is the vertical trace-free space. This fixed splitting follows directly from F4's vertical $G_{1/5}$ block, or N6(a).

The Einstein limits are
$$p_0=-\frac32\alpha_0-\frac92\beta_0,
\qquad q_0=-\frac52\alpha_0-\frac72\beta_0.$$
Consequently the limit of $H$ has eigenvalues $p=p_0+b$ with multiplicity five and $q=q_0+b$ with multiplicity nine. Its trace is $T=5p+9q$.

If $p=q=\delta$, (3) becomes $-33(A-C)^2\delta^2I=0$, so $\delta=0$. If $p\ne q$ and $T=0$, (3) instead forces $D^2p^2=D^2q^2=0$, a contradiction. If $p\ne q$ and $T\ne0$, they must be the two distinct roots $AT/D,CT/D$. Their two assignments to the five/nine spaces would require respectively
$$D=5A+9C\quad\text{or}\quad D=9A+5C.$$
These differ from $D=3A+11C$ by $2(A-C)$ and $6(A-C)$, both nonzero. Thus only the zero limit is possible.

It follows that $p_0=q_0=-b$. Their difference is $\alpha_0-\beta_0$, which vanishes exactly at $\lambda=1/5$. At that value both Einstein eigenvalues equal $15/(2c_V)$, giving the claimed $b$. This concludes a parameter restriction; it does not infer that $H$ vanishes away from the radial limit. ∎

## AK2. Three endpoint exclusions

The endpoints $a=-3/112,-3/272,-1/240$ admit no full solution under the assumptions above, for any $\lambda<1/4$.

*Proof.* At $a=-3/112$, (1) gives $G=-(b+4z)I$ with $b=309\kappa^2/1568>0$. In dimension fourteen this implies $\widehat{\operatorname{Ric}^Y}=(b+4z)I/6$, a positive scalar endomorphism. But [N1(d)](N-z1-neighbourhood.md) supplies a nonzero Euler field with $\widehat{\operatorname{Ric}^Y}E=-E/c$, where $c>0$. This is impossible. The skew derivative in (2) has not been discarded; only its symmetric part was used.

For the two opposite-branch endpoints, $A\ne C$. Equation (1) has eigenvalues $-Az$ with multiplicity three and $-Cz$ with multiplicity eleven. Its trace is $-Dz$, so (3) holds, including at $z=0$. AK1 forces $b=-15/(2c_V)<0$. Both rows instead have $b>0$, giving a contradiction. ∎

## AK3. Complete classification of the two branches

Under the stated assumptions:

| Single Lambda3 branch | Full solutions |
|---|---|
| $\mathsf c,\mathsf w\Gamma$ | Exactly $T_o=-i\kappa\Phi^1/400$, $T_e=0$, over a flat projected base, at $\lambda=1/5$ and $c_V\kappa^2=150000/11$. |
| $\mathsf w,\mathsf c\Gamma$ | None. |

*Proof.* AF/AG reduce each branch to the endpoint table. AK2 removes three rows. For the remaining row, (1) again gives (3); AK1 forces $\lambda=1/5$ and $b=-15/(2c_V)$, hence the exact scale. [AJ2](AJ-curved-endpoint.md) then forces a flat base on the domain's projection and $k=0$. [AI3](AI-flat-lambda3.md) supplies the complete-equation converse and regularity of the vector and other non-vector blocks. Every positive $c_H$ and either sign of $\kappa$ are allowed. ∎

## Evidence and scope

[AK_lambda3_classification.py](../../computation/geometric-unity/AK_lambda3_classification.py) and its [output](../../computation/geometric-unity/AK_lambda3_classification.out) give 28 exact controls: compact phases and endpoint kernels, full grade-one tensors, fourteen derivative directions for each copy, spectral and multiplicity identities, radial coefficients, and the complete surviving odd residual. The helper shares T's Clifford engine and uses the established geometric formulas as inputs. Finite calibrations support the stated equivariance arguments; they are not a substitute for the proof or a second metric implementation. [Internal Review38](../../reports/reviews/REVIEW38-GU-AK-2026-11-02.md) records the independent assessment.

Run from the repository root:

    python -X utf8 -I computation/geometric-unity/AK_lambda3_classification.py

The classification closes these two single kernels within the reconstruction. AJ's flatness input uses S5(b), including N10's computer-assisted curvature-pairing coefficient. Other kernels, coincident singular factors, Lambda5 and hook fields, other operators, $\det M=0$, indefinite metrics, non-invariant odd fields and a separate weighted base-metric variation remain open. The argument imposes no analyticity or curvature-trace sign. Literature novelty, external peer review and physical confirmation remain unestablished.
