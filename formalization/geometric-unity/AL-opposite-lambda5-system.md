# AL. The opposite Lambda5 constant endpoint system

The two remaining constant endpoints admit an exact reduction retaining the complete covariant jet, including its Cartan component. Each has one algebraic higher-grade constraint and one determined first-order equation. Exact first-jet controls permit nonzero norm currents. Smooth solutions on one fixed background, if they exist, are determined by one point value; existence remains unresolved.

Use [AH6](AH-lambda5-gram.md)'s hypotheses: a smooth positive four-base, reconstructed Z1, positive constant $c_H,c_V$, $\lambda<1/4$, positive fourteen-dimensional observer metric, compact $\mathfrak u(128)$, the true dualized equation, fixed parallel real $W_A$, real fields, constant $\kappa\ne0$ and $\det M\ne0$. On a connected open observer domain, fix $u\equiv0$ and one constant $a=t/\kappa$ throughout:

| $a$ | even field $k$ | $(q,\ell)$ | $b/\kappa^2$ |
|---|---|---|---|
| $-3/704$ | $\mathsf w(\omega)$ | $(112/3,320/3)$ | $87/61952$ |
| $-1/160$ | $\mathsf c(\omega)\Gamma$ | $(40,-32)$ | $19/3200$ |

Thus $T_o=i\kappa a\Phi^1$ and $T_e=k$. Remain in that single kernel, away from every other non-vector singular factor. Both compact phases are $1$. Absorb the nonzero chart weight into a smooth real five-form $\omega$, including all its derivatives. Here $\mathsf c(\eta)_a=\iota_a\eta$, $\mathsf w(\eta)_a=e^a\wedge\eta$ and $\Gamma^2=-1$. The slaved vector field vanishes identically, so $v=D(v)=0$. Isolated endpoints of varying odd fields are outside these assumptions.

## AL1. Complete derivative and inverse

For $J_a=\nabla_a\omega$, put $\operatorname{Alt}J=\sum_a e^a\wedge J_a=d\omega$, $\operatorname{Ctr}J=\sum_a\iota_aJ_a=-\delta_\nabla\omega$ and
$$K=I-\mathsf c\operatorname{Alt}-\mathsf w\operatorname{Ctr}.$$
The complete dualized derivatives are
$$\mathcal D(\mathsf w\omega)=-16iK(J),\qquad
\mathcal D(\mathsf c\omega\Gamma)=-8iK(J)\Gamma.\tag{AL1}$$
Their only grades are respectively five and nine; the previously allowed grades nine and thirteen vanish identically.

*Proof.* Every coordinate jet $e^j\otimes e^I$ has $j\in I$ or $j\notin I$. The helper compares every output slot for both representatives. Signed permutations cover every basis jet; an unused axis repairs determinant without changing the input. $SO(14)$ equivariance and linearity therefore prove (AL1) on the entire jet space. Parallel lifts include arbitrary covariant and weight derivatives.

On $T^*\otimes\Lambda^5=\Lambda^6\oplus\Lambda^4\oplus H_5$, the trace identities give $K$ eigenvalues $-5,-9,1$. Thus
$$K^{-1}=I-\tfrac15\mathsf c\operatorname{Alt}-\tfrac19\mathsf w\operatorname{Ctr}.$$
In particular, Cartan jets survive unchanged under $K$. ∎

## AL2. Equivalent full bosonic equation within the ansatz

Let $\Omega$ denote the Clifford image of $\omega$, $F_j=P_j(\Omega^2)$ and $Y_{d,a}=P_d(\Omega\gamma_a\Omega)$. Write $Q_d=[Q(k)]_d$. The universal formulas are

| copy | $Q_{5,a}/i$ | $Q_{9,a}/i$ |
|---|---|---|
| $\mathsf w\omega$ | $52Y_{5,a}+48P_5(\gamma_a\Omega^2)$ | $\frac{152}{3}Y_{9,a}+\frac{16}{3}P_9(\gamma_a\Omega^2)$ |
| $\mathsf c\omega\Gamma$ | $-\frac{44}{3}Y_{5,a}+\frac83P_5(\gamma_a\Omega^2)$ | $-\frac{40}{3}Y_{9,a}-\frac{56}{3}P_9(\gamma_a\Omega^2)$ |

Their gamma contractions are respectively $(272F_4,-512F_8/3)$ and $(256F_4/3,-176F_8/3)$; there are no other quadratic grades except one.

*Proof.* As in AH9, symbolic independent coefficients on two coordinate five-blades cover intersections zero through five, including diagonal terms. Complete output comparison, the unused-axis determinant repair and quadratic polarization prove the identities for arbitrary forms. Contracting uses $\sum_a\gamma_a\Omega\gamma_a=-4\Omega$ and square-term factors ten and six. ∎

Let $G$ be the raised observer Einstein tensor, $N=|\omega|^2$, $S_{ab}=\langle\iota_a\omega,\iota_b\omega\rangle$ and $b=\kappa^2a(1+312a)>0$. Within the standing ansatz, the full reconstructed bosonic equation $\Upsilon=0$ is equivalent to
$$G=-(b+qN)I-\ell S,\tag{AL2a}$$
and the appropriate row below:
$$\begin{array}{c|c|c}
\text{copy}&\text{algebraic constraint}&\text{determined covariant jet}\\\hline
\mathsf w\omega&Q_9=0&\nabla\omega=\frac1{16}K^{-1}(Q_5/i)\\
\mathsf c\omega\Gamma&Q_5=0&\nabla\omega=-\frac18K^{-1}(Q_9\Gamma/i).
\end{array}\tag{AL2b}$$
Indeed the even mass equation is the endpoint kernel identity. AH6 supplies grade one; (AL1) and the table exhaust all other odd grades. The conversion in the second row is **right multiplication by $\Gamma$**, giving $\mathcal D_9\Gamma/i=8K(J)$. No derivative is discarded. In particular, the algebraic constraints imply respectively $F_8=Y_9=0$ or $F_4=Y_5=0$, not the full AH9 simplicity criterion.

## AL3. Norm currents and nonconstant first jets

The determined jets simplify to
$$J_a=\tfrac{13}{4}Y_{5,a}+\tfrac{10}{9}e^a\wedge F_4,
\qquad
J_a=\tfrac53Y_{9,a}\Gamma+\tfrac{13}{15}\iota_a(F_8\Gamma),$$
respectively. Cyclic Clifford scalar trace, the positive five-form pairing and reversion give
$$dN=\tfrac{157}{18}[\Omega^3]_1,
\qquad dN=\tfrac85[\Omega^3\Gamma]_1.\tag{AL3}$$
For the first current the two pairings with $\Omega$ both give $[\Omega^3]_1$, so their coefficients sum to $13/2+20/9$. For the second they give $[\Omega^3\Gamma]_1$ and its negative, giving $10/3-26/15$. The axial current cannot be replaced by the ordinary vector current.

Exact real-form countercontrols are
$$\begin{aligned}
\mathsf w:\quad&\omega=e^{01234}+e^{12567}+e^{34567},&dN&=-\tfrac{157}{9}e^0,\\
\mathsf c\Gamma:\quad&\omega=e^{01234}+e^{05678}+e^{0,9,10,11,12},&dN&=\tfrac{48}{5}e^{13}.
\end{aligned}$$
Their forced jets satisfy both higher-grade equations and the first derivative of the unused algebraic constraint in all fourteen directions. They do not construct Z1 metric jets or smooth solutions. Synthetic $G$ controls check the complete odd residual, without asserting that such $G$ is an observer Einstein tensor. Curvature integrability, higher prolongations and actual metric realization remain unresolved.

## AL4. Conditional point-value uniqueness

Fix the smooth observer background, connection, parameters and endpoint copy. Two smooth solutions on a connected common domain agreeing at one point agree everywhere there.

*Proof.* Equation (AL2b) has the form $\nabla\omega=F(\omega)$ with smooth fibrewise quadratic $F$. Along any prescribed piecewise smooth path, parallel transport reduces it to an ordinary differential equation, locally Lipschitz in $\omega$. ODE uniqueness propagates equality along that path. Connected open domains are path connected. ∎

This supplies no existence theorem, admissibility of arbitrary point values, or uniqueness across different unknown metrics. The algebraic and metric constraints may prohibit the data.

**Flat-base boundary.** Neither remaining endpoint admits a solution domain over a locally flat base patch. Every positive fibre metric then extends locally as a parallel metric, so the entire observer lies in AH8's actual parallel set. The base Ricci rank is zero, and AH8 excludes both endpoints there. This requires a flat neighbourhood; curvature vanishing at just one base point is insufficient.

## Evidence and boundary

Run `python -X utf8 -I computation/geometric-unity/AL_opposite_lambda5_system.py`. Its exact controls share T's accepted Clifford engine and AH's Gram calculation. The exhaustive linear/quadratic basis certificates and hand trace proofs supply universality; finite controls do not replace AL4's ODE proof. This reduction covers only the two constant endpoints. Nonzero-$u$ fields, moving-odd isolated endpoints, coincident kernels, indefinite real forms and actual positive-Z1 solution existence remain open or outside its scope.
