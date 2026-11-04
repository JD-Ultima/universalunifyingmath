# BA. Covariant closure on the complementary pure hooks

The derivative equation on two remaining pure hook branches reduces to covariant exterior closure. This is a necessary condition on the full field, with every tensor index differentiated. It does not by itself prove existence or exclude these branches.

Code: [exact Clifford controls](../../computation/geometric-unity/BA_hook_closure.py) and [output](../../computation/geometric-unity/BA_hook_closure.out). The earlier first-hook case is [AX1](AX-hook-curvature.md); [AZ](AZ-first-hook-classification.md) supplies its subsequent positive-Z1 exclusion.

## Setting

Use the true dualized bosonic equation, the fixed parallel real member $W_A$ with V's normalization, observer signature $(14,0)$, compact gauge algebra $\mathfrak u(128)$, real smooth fields and constant $\kappa\ne0$. The Clifford generators are Hermitian and $\Gamma^2=-1$. Work locally in an oriented orthonormal frame; all displayed identities are tensorial.

For $q\in\{2,6,10\}$, let $B\in T^*Y\otimes\Lambda^qT^*Y$ be a real hook tensor:
$$
\operatorname{Alt}B=\sum_a e^a\wedge B_a=0,
\qquad \operatorname{Ctr}B=\sum_a\iota_aB_a=0.
$$
The real even field is $K=iB\Gamma$. Its actual Clifford degree is $14-q$, respectively $12,8,4$, so this phase is anti-Hermitian. Define
$$
D=\sum_b\nabla_bB_b,\qquad
(d_\nabla B)_a=\sum_b e^b\wedge\nabla_bB_a.
$$
Here $\nabla$ acts on the coefficient covector and every form index. $D$ is a $q$-form, not the form-slot codifferential of $B$.

## BA1. The complete derivative identity

Writing $\mathcal D K={*}^{-1}\circledast(d_0K)$ for the $W_A$ derivative term,
$$
\boxed{(\mathcal D K)_a\Gamma
=2\big((d_\nabla B)_a-e^a\wedge D\big).}\tag{BA1}
$$
In particular, $\mathcal D K=0$ is equivalent to $d_\nabla B=0$ for these three degrees. There is no additional lower-grade equation left out of (BA1).

*Proof.* At a normal-frame point put $H_{bc}=\nabla_bB_c-\nabla_cB_b$ and let $A_{abc}=e^a\wedge e^b\wedge e^c$ be the alternating Clifford trivector, zero for repeated indices. The real $W_A$ phase and right multiplication by $\Gamma$ give
$$
(\mathcal D K)_a\Gamma=\sum_{b<c}[A_{abc},H_{bc}].
$$
For any real even $q$-form $H$, Clifford multiplication gives
$$
\frac12[A_{abc},H]
=e^a\wedge e^b\wedge\iota_cH
-e^a\wedge e^c\wedge\iota_bH
+e^b\wedge e^c\wedge\iota_aH
+\iota_a\iota_b\iota_cH.
\tag{BA1.1}
$$
Interior operators act from right to left. The last term has degree $q-3$ and can be nonzero before imposing the differentiated hook conditions. Its full sum is
$$
\sum_{b<c}2\iota_a\iota_b\iota_cH_{bc}
=2\iota_a\sum_b\iota_b\nabla_b
  \left(\sum_c\iota_cB_c\right)=0.
$$
For the first two wedge terms use
$$
\sum_{b,c}e^b\wedge\iota_cH_{bc}=-D.
$$
Indeed, one summand differentiates $\operatorname{Ctr}B=0$; contracting the differentiated identity $\operatorname{Alt}B=0$ evaluates the other. The remaining wedge term satisfies
$$
\sum_{b<c}e^b\wedge e^c\wedge\iota_aH_{bc}
=\sum_b e^b\wedge\nabla_bB_a.
$$
This follows by applying $\iota_a$ to $\sum_c e^c\wedge\nabla_bB_c=0$. Substitution proves (BA1), including the cancellation of every lower-grade term.

The same differentiated hook identities give
$$
\sum_a\iota_a(d_\nabla B)_a=D,
\qquad
\sum_a\iota_a(e^a\wedge D)=(14-q)D.
$$
Thus contraction of (BA1) is $-2(13-q)D$. Its coefficients are $-22,-14,-6$, all nonzero. The forward implication follows; conversely $d_\nabla B=0$ already forces $D=0$. ∎

## BA2. Necessary conditions on the pure field branches

Restrict the invariant odd field to $T_o=i(t\Phi^1+u\Phi^1\Gamma)$ and the following constant pairs. By [R2](R-even-operator.md) and [V1](V-hook-intersections.md), the vector block is regular, its source is zero and only the indicated even block is singular:

| Odd pair | Even field | Hook variable | Derivative grade |
|---|---|---|---:|
| $t=3\kappa/136$, $u=0$ | $iB_2\Gamma$ | $B_2\in H_2$ | 11 |
| $t=3\kappa/104$, $u=0$ | $iB_6\Gamma$ | $B_6\in H_6$ | 7 |
| $t=\kappa/24$, $u=0$ | $iC_4$ | $B_{10}=-C_4\Gamma\in H_{10}$ | 3 |

Right $\Gamma$ interchanges the alternating and contraction traces, so the last identification is an isomorphism of the hook spaces. In each row an arbitrary real even field satisfying the even equation has exactly the displayed form.

The curvature term and the constant invariant baseline have grade one. For an actual even field $iA_r$, $r=12,8,4$, the real $W_A$ quadratic term $Q(iA_r)/i$ has real Clifford coefficients and is Hermitian. Consequently its odd degrees are $1\pmod4$. These are disjoint from the derivative degrees $11,7,3$ in the table. Each full solution therefore satisfies the entire derivative equation separately, hence
$$
\boxed{d_\nabla B_q=0,\qquad D=0,\qquad
F^{T^*Y}\wedge B_q=0.}\tag{BA2}
$$
The last condition follows from $d_\nabla^2B_q=F^{T^*Y}\wedge B_q$; this curvature acts on the coefficient covector. The form-slot connection is already included in the covariant exterior derivative.

These are necessary conditions, alongside every remaining quadratic and metric equation. AZ has closed the first row throughout positive Z1. BA2 does not extend that exclusion to the second or third row. For the third row, parallel metric Hodge duality also gives the equivalent condition $\delta_\nabla C_4=0$ on the form slots.

An injectivity argument of AZ's type cannot work on $H_{10}$ alone: its dimension is $11\,648$, whereas $T^*Y\otimes\Lambda^{12}T^*Y$ has dimension $1\,274$. Thus the curvature-wedge kernel there has dimension at least $10\,374$ for every metric. No rank for the $H_6$ curvature map is asserted here.

## Exact certificate and limits

The standalone helper checks the full Clifford commutator, the lower-grade cancellation, every output slot, the right-$\Gamma$ degree and the contracted coefficient. It uses two generator types, with distinct axes $p,s$ and a $(q-1)$-blade $E_J$ avoiding both:
$$
\begin{array}{lll}
\mathrm D:&B_p=e^s\wedge E_J,&B_s=e^p\wedge E_J,\\
\mathrm R:&B_p=e^p\wedge E_J,&B_s=-e^s\wedge E_J.
\end{array}
$$
All other external components vanish. For coefficients whose external index is absent from the form blade, fixing their union leaves one alternating constraint; the D generators span its kernel. For coefficients with the external index present, fixing the remaining blade leaves one contraction constraint; the R generators span its kernel. Their disjoint dimensions are
$$
q\binom{14}{q+1}+(14-q)\binom{14}{q-1}
=14\binom{14}{q}-\binom{14}{q+1}-\binom{14}{q-1}.
$$
They give $896,36\,608,11\,648$. Signed coordinate permutations carry every generator to its representative. A derivative covector adds at most one axis, so the input uses at most $q+2\le12$ axes. Reflection of an unused axis repairs the determinant to $SO(14)$ without changing the input. Linearity and equivariance therefore make the two representatives times fourteen derivative directions times three degrees a complete 84-case cover of the linear identity, rather than a sample of fields. Unconstrained blades separately show that the lower-degree terms in (BA1.1) exist before cancellation.

The hand proof supplies the universal identity and orbit coverage; exact integer Clifford arithmetic checks the conventions. No floating-point arithmetic, previous helper or recorded output is imported. Neither these controls nor the dimension bound solve the remaining nonlinear equations. Coincident kernels, the other pure hook copies, other real forms or operators and physical identification remain separate questions.
