# AV. Cubic Gram transport on the opposite conic

On every regular solution patch of [AS's opposite Lambda5 system](AS-opposite-lambda5-system.md), two symmetric cubic Clifford tensors vanish on the hyperplane perpendicular to the odd gradient. The proof reduces the derivative of the actual Gram using full AS compatibility. [AW](AW-lambda5-conic-classification.md) now excludes those actual domains by a separate normal-moment and conservation argument. Whether the cubic condition already follows algebraically from the zero-order equations remains unresolved.

The standalone [helper](../../computation/geometric-unity/AV_opposite_gram_transport.py) and [output](../../computation/geometric-unity/AV_opposite_gram_transport.out) record 13 exact controls. Its finite forms check conventions; universality follows from the hand identities below. [Metadata](../../metrics/geometric-unity/AV-opposite-gram-transport.yaml) records the scope. [Independent internal Review52](../../reports/reviews/REVIEW52-GU-AV-2026-11-03.md) accepts AV1–AV3 and reproduces the helper. External peer review remains outstanding.

## Setting and the canonical jet

Retain AS's smooth positive four-base, positive Z1, constant $c_H,c_V>0$, $\lambda<1/4$, compact real $(14,0)$ and $\mathfrak u(128)$, true dualized equation, fixed parallel real $W_A$, real fields, constant $\kappa\ne0$, $\det M\ne0$, and one opposite $\Lambda^5:\mathsf w,\mathsf c\Gamma$ kernel away from coincident factors. Work on a regular patch $da\ne0$ with $w\ne0$.

Absorb the nonzero weight $r$ into $\eta=r\omega$, and put $\tau=s/r$. Then
$$k=\mathsf w(\eta)+\tau\mathsf c(\eta)\Gamma,\qquad
\tau=\frac{3(1+704a/3)}{80w}=\frac{48w}{1+160a}.$$
This change retains the weight derivative: $\nabla\eta=r\nabla\omega+dr\otimes\omega$. Set $x=\tau^2$ and use AS's $X,Y,Z$ patterns with $\eta$:
$$X_i=C_i\eta,\qquad Y_i=\iota_i(C\wedge\eta),\qquad
Z_i=e^i\wedge\iota_C\eta,\qquad C=24da/\delta.$$
Here $x$ is a scalar chart variable, not a base point. Write $\Omega$ for the Clifford image of $\eta$, and retain AS's $F_j,Y_j,E_j$ notation. All Gamma products are right products and $\Gamma^2=-1$.

**Lemma AV1 (canonical jet).** AS2 gives $\nabla_i\eta=\mathcal F_i$, where
$$\begin{aligned}
\mathcal F_i={}&\left(\frac{13}{4}-\frac{11x}{12}\right)Y_{5,i}
+\left(\frac{10}{9}-\frac{23x}{54}\right)E_{5,i}
+\frac{7\tau}{2}Y_{9,i}\Gamma+\frac{23\tau}{10}E_{9,i}\Gamma\\
&+\gamma(a)X_i-\frac{44}{15}Y_i-\frac73Z_i,
\qquad \gamma(a)=\frac{246-9152a}{15}.
\end{aligned}\tag{AV1}$$

*Proof.* AS's parallel operator satisfies
$$K^{-1}X=X-Y/5-Z/9,\quad K^{-1}Y=-Y/5,\quad K^{-1}Z=-Z/9.$$
The $(\operatorname{Alt},\operatorname{Ctr})$ traces of the four quadratic patterns are respectively
$(0,-4F_4)$, $(0,10F_4)$, $(-4F_8\Gamma,0)$ and $(6F_8\Gamma,0)$. These follow by separating grades in $\sum_i\gamma_i\Omega\gamma_i=-4\Omega$. Hence
$$\begin{aligned}
K^{-1}Y_5&=Y_5+4E_5/9,&K^{-1}E_5&=-E_5/9,\\
K^{-1}(Y_9\Gamma)&=Y_9\Gamma+4E_9\Gamma/5,
&K^{-1}(E_9\Gamma)&=-E_9\Gamma/5.
\end{aligned}$$
Apply these identities to AS1a–AS1b and divide by sixteen. The $P_{14}$ mixed term adds $-208w\tau X/3$; since $w\tau=3(1+704a/3)/80$, its sum with $19X$ is the stated $\gamma X$. No derivative is frozen. ∎

## Pointwise reduction of the Gram variation

Use an oriented orthonormal frame with $C=\sqrt c\,e^0$, $c=|C|^2>0$; indices $i,j,k$ below are perpendicular to $C$. For a vector-valued five-form $U$, define
$$\mathcal B(U)_{ijk}
=\langle\iota_jU_i,\iota_k\eta\rangle+
\langle\iota_j\eta,\iota_kU_i\rangle.$$
This is the raw Gram variation. Taking its tracefree part in $j,k$ gives the tracefree transport tensor. The following reduction first uses the raw tensor, before removing traces.

**Lemma AV2 (zero-order coefficient identity).** Assume the actual zero-order Gram equation and full AS compatibility at the point. Put $H_5=\mathcal B(Y_5)$, $H_9=\mathcal B(Y_9\Gamma)$ and $D=\mathcal B(Y)$. Then $H_5,H_9$ are fully symmetric, $D$ has cyclic sum zero, and
$$A H_5+B H_9=0,\qquad dD=0,\tag{AV2a}$$
$$\mathcal B(\mathcal F)=\Lambda(x)H_5+fD,\qquad
\Lambda(x)=\frac{741+436x+55x^2}{6(38+11x)}>0.\tag{AV2b}$$
Here, with AS3's canonical $\alpha,\zeta,\vartheta$,
$$\begin{aligned}
A&=-\tfrac43\tau(41+11x),&B&=\tfrac83(38+11x),\\
d&=\tfrac83(4+x)(\zeta-\tau\alpha)-48\tau,&
f&=-\tfrac35+\tfrac{23}{10}\tau\zeta-
\left(\tfrac{10}{9}-\tfrac{23x}{54}\right)\alpha.
\end{aligned}$$
Both $A$ and $B$ are nonzero in a real nonzero chart.

*Universal trace proof.* Let $\operatorname{sc}$ denote scalar Clifford coefficient, equivalently normalized cyclic matrix trace. Define
$$T_{ijk}=\operatorname{sc}(\Omega\gamma_i\Omega\gamma_j\Omega\gamma_k),\qquad
T^\Gamma_{ijk}=\operatorname{sc}(\Omega\gamma_i\Omega\gamma_j\Omega\gamma_k\Gamma).$$
Both tensors are fully symmetric. Cyclicity cyclically permutes their indices. Reversion fixes $\Omega$ and the generators and reverses their order. In the Gamma case, reversion changes the sign of $\Gamma$; cyclically moving the last $\Omega$ to the front and moving $\Gamma$ through the resulting five odd factors supplies a compensating sign. Thus reversal of the index order also preserves both tensors.

The five-form identity
$$e^j\wedge\iota_k\eta+e^k\wedge\iota_j\eta
=\delta_{jk}\Omega+\tfrac12(\gamma_j\Omega\gamma_k+\gamma_k\Omega\gamma_j)$$
gives
$$\mathcal B(Y_5)_{ijk}=\delta_{jk}\xi_i+T_{ijk},\qquad
\mathcal B(Y_9\Gamma)_{ijk}=\delta_{jk}\upsilon_i-T^\Gamma_{ijk}.$$
AS3 and the actual Gram alignment imply, on tangent indices,
$$\xi_i=\langle F_4,\iota_i\eta\rangle=\alpha S(C,e_i)=0,\qquad
\upsilon_i=-\langle\iota_i(F_8\Gamma),\eta\rangle=0.$$
The second vanishes because it is a multiple of $C_i|\eta|^2-S(C,e_i)$. Therefore $H_5=T$ and $H_9=-T^\Gamma$. Full symmetry is proved without assuming tracefreeness.

For the remaining patterns write $\eta=e^0\wedge\phi+\psi$. Actual Gram alignment gives $\langle\phi,\iota_i\psi\rangle=0$. Direct contraction yields
$$\mathcal B(X)=0,\qquad \mathcal B(Z)=-D,$$
$$D_{ijk}=\sqrt c\,\{\langle\iota_i\iota_j\psi,\iota_k\phi\rangle+
\langle\iota_i\iota_k\psi,\iota_j\phi\rangle\},\qquad
D_{ijk}+D_{jki}+D_{kij}=0.$$

Apply $\mathcal B$ to the **full** AS3 compatibility, using $E_5=\alpha Z$ and $E_9\Gamma=\zeta Y$. Its result is $AH_5+BH_9+dD=0$. Fully symmetrizing, and then subtracting, proves (AV2a). Apply the same map to (AV1) and eliminate $H_9$ using $B\ne0$; ordinary rational coefficient simplification gives (AV2b). ∎

## Necessary cubic transport and its unresolved dependency

**Theorem AV3.** Every actual regular solution patch satisfies
$$\boxed{T(X,Y,Z)=T^\Gamma(X,Y,Z)=0\qquad(X,Y,Z\perp C).}\tag{AV3}$$

*Proof.* AQ's regular local identities and the actual Gram equation give
$$S=\lambda_\perp(a)g-\beta(a)C^2,\qquad |C|^2=c(a).$$
For example $\beta=13(1+160a)/19$ in the canonical chart; tracing AH2 makes $|\eta|^2$, and then $\lambda_\perp$, functions of $a$. Differentiating this tensor identity and evaluating on $X,Y,Z\perp C$ gives $\mathcal B(\mathcal F)=0$. The scalar derivative contains $da(X)$; every derivative of $C^2$ contains $C(Y)$ or $C(Z)$. All vanish. No derivative of a moving projector has been omitted.

It remains to retain the exceptional coefficient levels. The conic gives
$$a=\frac{5x-9}{32(66-25x)},\quad
w=\frac{7\tau}{16(66-25x)},\quad
\frac{dx}{da}=\frac{672}{5(1+160a)^2}>0.$$
Substitution in $d$ gives
$$d=\frac{16P(x)}{315\tau(16x+61)(19x+64)(25x-66)},$$
$$P(x)=1568750x^5+565335x^4-50705892x^3-81149677x^2
+138437232x+46051632.$$
The denominators are nonzero in a finite nonzero-$w$ chart. Since $P$ is nonzero, its exceptional real levels are finite. On a regular patch $a$ is a local coordinate and $x(a)$ is strictly increasing, so $d\ne0$ is dense. Equation (AV2a) gives $D=0$ there and smoothness extends this across the exceptional levels. Equation (AV2b) and strict positivity of $\Lambda$ now force $H_5=0$; (AV2a) forces $H_9=0$. ∎

At a fixed parameter with $d\ne0$, the zero-order equations already give $D=0$. The first tracefree tangential Gram-preservation question then reduces exactly to the pair-tracefree part of $T$ vanishing. Whether that part, or all of $T$, is already forced on the allowed real zero-order Gram and full AS compatibility data remains unresolved. No datum satisfying all those equations with nonzero residual $T$ is supplied. A polynomial-ideal identity could prove implication; ordinary ideal nonmembership alone would not supply real admissible counterdata. The theorem is a necessary transport condition, not a proof of its independence, five-form simplicity or field existence. Critical points are outside its regular-patch statement.

## Evidence and limits

Run from the repository root:

    python -X utf8 -I computation/geometric-unity/AV_opposite_gram_transport.py

The recorded run passes all 13 controls in 0.507 seconds; Review52 reproduces them in 0.514 seconds. The helper uses standalone exact positive Clifford multiplication and SymPy, with no repository-engine or private-file import. Nine finite convention checks compare the two trace identities and their symmetries in all 2744 slots on three declared forms. Those forms are not asserted to satisfy AS's zero-order equations and are not a universal basis certificate. Four scalar controls check the positive coefficient, nonzero compatibility factors, exact polynomial numerator and strictly positive chart derivative. Cyclic trace, reversion, contraction anticommutation and the covariant argument supply the universal statements.

The necessary condition retains AS/AQ's field, metric and single-kernel hypotheses. AW closes actual-domain existence under those hypotheses, while the separate algebraic implication and real zero-order feasibility questions remain open. Coincident kernels, $\det M=0$, other operators, non-invariant odd fields and indefinite real forms are outside this result. No physical or literature-first claim is made.