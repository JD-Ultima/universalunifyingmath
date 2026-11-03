# AH. Lambda5 Gram and curvature constraints

At a nonzero-$u$ point of either single $W_A$ Lambda5 conic, the symmetric Gamma equation forces contraction-Gram rank thirteen or fourteen. Left arcs require full support; on right arcs, rank thirteen is possible only when the slaved gradient is the annihilator direction. Forms supported in dimension at most twelve are excluded. AH5 gives an earlier cyclic-orbit obstruction. [AP](AP-first-lambda5-obstruction.md) excludes the entire first $c/w\Gamma$ conic. The remaining opposite $w/c\Gamma$ conic is restricted geometrically by [AR](AR-opposite-lambda5-geometry.md) and expressed as a full system in [AS](AS-opposite-lambda5-system.md).

On positive Z1, neither constant-odd endpoint of the first $c/w\Gamma$ pair admits a full solution (AH9). Its higher quadratic equations force a nonzero simple five-form, whose fixed normalized Gram spectrum contradicts AH7. For the opposite pair, AH6–AH8 supply necessary norm, spectrum and Ricci-index constraints. The later [AL reduction](AL-opposite-lambda5-system.md) and [AN Bianchi obstruction](AN-opposite-endpoint-obstruction.md) exclude those two constant endpoints as well. These four endpoint exclusions concern entire constant-odd domains; AP now closes the first moving conic, while the moving opposite-conic problem remains open.

Code: `computation/geometric-unity/AH_lambda5_gram.py` → `.out`, with the cyclic-orbit certificate in `AH_cyclic_obstruction.py` → `.out`. Record: `metrics/geometric-unity/AH-lambda5-gram.yaml`. Both helpers use T's accepted Clifford implementation and exact rational/polynomial controls, with no off-repository files. [Independent internal Review31](../../reports/reviews/REVIEW31-GU-AH-2026-11-02.md) accepts AH1–AH4 and reproduces all 25 original checks. [Review34](../../reports/reviews/REVIEW34-GU-AH-CYCLIC-2026-11-02.md) separately accepts AH5 and reproduces its 14 checks. [Review39](../../reports/reviews/REVIEW39-GU-AH-ENDPOINTS-2026-11-02.md) accepts AH6–AH7 and reproduces the 74-control helper through AH7. [Review40](../../reports/reviews/REVIEW40-GU-AH8-2026-11-02.md) accepts AH8 and reproduces the 94-control helper through AH8. [Review41](../../reports/reviews/REVIEW41-GU-AH9-2026-11-02.md) records the independent AH9 proof and coefficient audit. These are internal reviews, not external peer review.

## Setting and the two real kernels

Keep the true dualized equation, real compact $(14,0)$, gauge algebra $\mathfrak u(128)$, fixed parallel real $W_A$, constant $\kappa\ne0$, invariant odd field and $\det M\ne0$. Put
$$T_o=i\kappa(a\Phi^1+w\Phi^1\Gamma),\qquad
a=t/\kappa,\quad w=u/\kappa.$$
Stay in one single Lambda5 kernel, away from every other non-vector singular factor. The conic data are

| pair | $d$ | $e$ | $p$ | $q$ | compact kernel field |
|---|---:|---:|---:|---:|---|
| $c/w\Gamma$ | $96$ | $832/3$ | $48$ | $80/3$ | $k=i[r\mathsf c(\omega)+s\mathsf w(\omega)\Gamma]$ |
| $w/c\Gamma$ | $704/3$ | $160$ | $80/3$ | $48$ | $k=r\mathsf w(\omega)+s\mathsf c(\omega)\Gamma$ |

Here $\omega$ is a real five-form, $r,s$ are real smooth chart weights not both zero, and
$$w^2=f(a)=\frac{(1+da)(1+ea)}{pq},\qquad pq=1280.$$
The phase is $i$ for grades four/eight and $1$ for grades six/ten. In either pair,
$$(1+da)r=pws,\qquad (1+ea)s=qwr.$$
The charts $(r,s)=(pw,1+da)$ and $(1+ea,qw)$ cover the conic, including zero-$w$ endpoints. At those endpoints no smooth dependence $w=w(a)$ is presumed. Both weights are nonzero where $w\ne0$.

The full even field is $T_e=v+k$, with
$$v=P_2(C)+P_{14}(B_0),\quad
C=\frac{24}{\delta}da,\quad B_0=208wC,\quad
\delta=1+224a-3328f(a).$$
Both vector-block factors are nonzero by $\det M\ne0$.

## AH1. A universal five-form stress tensor

Define the form norm and contraction Gram tensor
$$N=|\omega|^2,\qquad
S_{ab}=\langle\iota_{e_a}\omega,\iota_{e_b}\omega\rangle.$$
These are the ordinary positive form inner products; $N$ is not the norm of $k$. For either compact pair above,
$$\boxed{\operatorname{Sym}([Q(k)]_{13}\Gamma/i)
=\frac{160}{3}rsN I-\frac{304}{3}rsS.}\tag{AH1}$$
The conversion is **right multiplication by $\Gamma$**, as in AF and AG.

*Proof.* After the Gamma products are simplified, a symmetric two-tensor built quadratically from two alternating five-tensors has only two metric-contraction patterns. Pairing the two output indices gives $NI$. Pairing one output index to each input and contracting the remaining four input pairs gives $S$. A contraction within one alternating input vanishes. There are twelve total input/output tensor slots, below dimension fourteen, so an epsilon contraction supplies no additional pattern.

For the unit simple form $e^{01234}$, the exact full symbolic-weight output has eigenvalues $-48rs$ on its five-plane and $160rs/3$ on its perpendicular nine-plane, for both real phases. On that form $N=1$ and $S$ is the five-plane projector, so the complete output fixes both coefficients in (AH1). This determines the formula for arbitrary five-forms, rather than only for the simple-form orbit. Nondecomposable, off-diagonal-Gram and full-support controls separately check it. ∎

The kernel derivative and mixed vector/kernel term have no symmetric converted grade-thirteen component. Their input type is $V\otimes\Lambda^5=\Lambda^6\oplus\Lambda^4\oplus H_5$, with no symmetric two-tensor summand. This includes every covariant derivative of the form and the chart weights; no derivative factorization from AA4 is imposed on the other pair. The invariant, curvature and vector terms are those already established in AF1. Every full solution therefore satisfies
$$\boxed{\left[-\kappa^2w(1+208a)+3328w|C|^2+\frac{160}{3}rsN\right]I
=3328wC\otimes C+\frac{304}{3}rsS.}\tag{AH2}$$
No equation for the slaved baseline or curvature-trace sign has been assumed. The other odd components, including the remaining derivative tensor, are retained.

## AH2. Rank and support at a nonzero-$w$ point

Suppose a full solution has $w\ne0$. Then $r,s\ne0$ and $N>0$.

*Proof that $N$ cannot vanish.* If $N=0$, positivity gives $\omega=0$ and $S=0$. Equation (AH2) has a scalar identity on the left and rank at most one on the right, so its scalar coefficient and $C$ vanish. It forces $a=-1/208$. But
$$f_{c/w\Gamma}(-1/208)=-7/49920<0,\qquad
f_{w/c\Gamma}(-1/208)=-1/43264<0,$$
contradicting real $w^2=f(a)$. ∎

Divide (AH2) by its nonzero $rs$ coefficient and write
$$S=\alpha I-bC\otimes C,\qquad b=\frac{624}{19}\frac{w}{rs}.\tag{AH3}$$
The Gram matrix is positive semidefinite, has trace $5N$, and has rank at least five. Indeed it is $L^*L$ for $L(X)=\iota_X\omega$, so its kernel is exactly the annihilator of $\omega$. The form is supported on that kernel's perpendicular complement; a nonzero five-form needs support dimension at least five.

If $C\ne0$, (AH3) has eigenvalue $\alpha$ on $C^\perp$ and $\alpha-b|C|^2$ on its line. Positivity gives $\alpha\ge0$; equality would give rank at most one, impossible. Thus $\alpha>0$, and the rank is thirteen or fourteen. If $C=0$, positivity and $N>0$ likewise give rank fourteen. In all cases,
$$\boxed{\operatorname{rank}S\in\{13,14\}.}\tag{AH4}$$
Equivalently, a full solution at nonzero $w$ cannot have its five-form supported in dimension at most twelve. No classification of the forms having the permitted Gram spectrum is asserted.

## AH3. Branch signs and gradient bounds

Put $A=1+da$. On the real conic with $w\ne0$, $A$ and $1+ea$ have the same nonzero sign. The kernel equations give
$$\frac{rs}{w}=\frac{p s^2}{A}.$$
Thus $b<0$ on a left branch and $b>0$ on a right branch. The branch endpoints are

| pair | left: $a<\min(-1/d,-1/e)$ | right: $a>\max(-1/d,-1/e)$ |
|---|---:|---:|
| $c/w\Gamma$ | $a<-1/96$ | $a>-3/832$ |
| $w/c\Gamma$ | $a<-1/160$ | $a>-3/704$ |

Tracing (AH3) yields
$$\alpha=\frac{5N+b|C|^2}{14},\qquad
\alpha-b|C|^2=\frac{5N-13b|C|^2}{14}.$$
On a left branch the second eigenvalue is larger than $\alpha>0$, so $S$ has full rank fourteen and
$$|b|\,|C|^2<5N.$$
In particular it has no nonzero annihilator direction. On a right branch,
$$b|C|^2\le\frac{5N}{13}.$$
Rank thirteen occurs exactly at equality, when $C\ne0$ and $\ker S=\mathbb RC^\sharp$, with $S$ identified as an endomorphism using $g_Y$; equivalently $\iota_{C^\sharp}\omega=0$. Otherwise the rank is fourteen. These are necessary algebraic restrictions, not sufficient conditions for the other equations.

## AH4. A full-support limitation of this channel

Let $B=\{0,1,3,6,10\}\subset\mathbb Z/14$, and take the sum of the fourteen basis five-forms whose supports are the cyclic shifts of $B$, each written in increasing order with coefficient $+1$. The supports are distinct and any two intersect in at most two indices. Hence two contractions cannot share a four-index support, so the Gram matrix has no off-diagonal entries. Every coordinate lies in five supports. Exactly,
$$N=14,\qquad S=5I.$$
The isotropic case permitted by (AH3) is therefore nonempty.

There are regular first-conic data satisfying this selected channel with $C=0$:
$$a=1/96,\quad w=\sqrt{14}/48,\quad r=\sqrt{14},\quad s=2,
\quad\kappa=96\sqrt{15/19}.$$
They obey $f(a)=7/1152$, both kernel relations, $\delta=-152/9$ and the other vector factor $-122/9$. All other non-vector factors are regular, as checked by the exact block certificates. For these data,
$$\kappa^2w(1+208a)=480\sqrt{14}=240rs,$$
which is precisely (AH2) when $N=14,S=5I,C=0$.

This is an algebraic channel control. No observer metric, field jet satisfying all components, or full PDE solution has been constructed. The Gram restriction admits this form, while AH5 below excludes its declared orbit from full solutions on the first branch at nonzero $w$.

## AH5. The cyclic form fails the full equation on the first branch

Let $\Omega$ be AH4's exact fourteen-term form. Keep the standing hypotheses, and restrict to the first $\Lambda^5:c/w\Gamma$ branch at a point with $w\ne0$.

**Proposition AH5.** No full solution has $\omega=hR\Omega$ at that point, where $h\ne0$ is real and $R\in SO(14)$. This is a pointwise statement about one scaling/rotation orbit, not arbitrary five-forms.

*Proof.* First take $R=I$. Its Gram tensor is $5h^2I$. Moving the isotropic terms in (AH2) to one side gives a scalar identity equal to $3328wC\otimes C$. Its two sides have incompatible ranks unless $C=0$, since $w\ne0$. Thus $v=0$ at the point. No derivative of $v$ is assumed zero.

The kernel inputs have Clifford grades four and eight. In the $W_A$ derivative, a triple anticommutator with an even monomial retains even intersection sizes. The first input therefore yields grades three/seven, and the second grades seven/eleven. Connection and chart-weight derivatives preserve this exclusion of grades five/nine. Vector derivatives have input type $V\otimes V$, with no component in $V\otimes\Lambda^5$ or its right-Gamma dual. The invariant baseline has only grades one/thirteen, the curvature term grade one, and the quadratic vector and mixed terms vanish at $v=0$. The full grade-five/nine equations are consequently just those components of $Q(k)$.

Direct integer Clifford polarization of the fixed form gives the following coefficients of $Q(k)/i$ for $h=1$:

| external index | Clifford indices | coefficient |
|---|---|---|
| $0$ | $(2,3,4,5,6)$ | $-\frac83(7r^2-27s^2)$ |
| $0$ | $(1,2,4,7,9)$ | $-160rs$ |
| $0$ | $(0,1,3,4,6,8,9,10,11)$ | $-\frac{16}{3}(4r^2+17s^2)$ |

The grade-nine coefficient is nonzero for every real nonzero chart. Scaling multiplies it by $h^2$. Equivariance makes the entire nonzero tensor remain nonzero under $SO(14)$ rotations. This proves the proposition. No general five-form contraction classification is used: the grade-five target has sixteen total input/output slots, so AA6's twelve-slot argument would not apply. ∎

At AH4's regular weights $(r,s)=(\sqrt{14},2)$, the displayed coefficients are $80/3$, $-320\sqrt{14}$ and $-1984/3$. A direct phased call to the accepted T engine checks the polarization in every grade-five/nine slot and confirms that the same control still balances AH4's symmetric Gamma channel. The three coefficient polynomials also generate $(r^2,rs,s^2)$ exactly.

The portable certificate `computation/geometric-unity/AH_cyclic_obstruction.py` includes the compact phase, real conic control and nonzero derivative controls. Arbitrary five-forms, the other Lambda-five pair, zero-$w$ endpoints and coincident kernels remain outside AH5.

## AH6. Constant-odd endpoint domains: exact Euler/Gram bounds

For this section also assume the reconstructed Z1 metric over a smooth positive four-base, with constant $c_H,c_V>0$ and $\lambda<1/4$. Put
$$c=4c_V(1-4\lambda)>0.$$
On an entire nonempty open observer domain, let $a$ equal one of the four endpoint constants below and let $u\equiv0$. Retain the single-kernel and regularity hypotheses above. Then $T_o=i\kappa a\Phi^1$ is parallel, the slaved vector field $v$ vanishes identically, and so does $D(v)$. A nonzero chart weight can be absorbed into a smooth real five-form $\omega$; the form itself may vanish. No simplicity assumption is made.

This **constant-odd domain** assumption is essential. An isolated zero-$u$ point on a varying conic need not have $D(v)=0$, even if $v=0$ there, and is outside AH6.

Write $G$ for the raised observer Einstein tensor and
$$b=\kappa^2a(1+312a)>0.$$
The complete quadratic grade-one tensor of each pure kernel copy is
$$[Q(k)]_1/i=q_1NI+\ell S.$$
Here $S$ and $N$ use the positive form norm defined in AH1. The constants and retained derivative grades are:

| $a$ | compact kernel $k$ | $b/\kappa^2$ | $q_1$ | $\ell$ | grades of $D(k)$ |
|---|---|---:|---:|---:|---|
| $-1/96$ | $i\mathsf c(\omega)$ | $3/128$ | $40/3$ | $32/3$ | $3,7$ |
| $-3/832$ | $i\mathsf w(\omega)\Gamma$ | $3/6656$ | $112$ | $-64$ | $7,11$ |
| $-3/704$ | $\mathsf w(\omega)$ | $87/61952$ | $112/3$ | $320/3$ | $5,9$ |
| $-1/160$ | $\mathsf c(\omega)\Gamma$ | $19/3200$ | $40$ | $-32$ | $9,13$ |

*Universal tensor and derivative proof.* The twelve-slot metric-contraction count used in AH1 applies to the grade-one output as well. The only two-tensor contractions of two alternating five-tensors are $NI$ and $S$, both symmetric; an epsilon term cannot occur. The exact output on a unit simple form fixes the two constants in each row, and hence the full tensor for arbitrary forms. Separate nondecomposable, off-diagonal-Gram and full-support controls check its conventions.

For the first pair, AA4 at form degree five gives the complete covariant identities, with $\delta_\nabla$ the codifferential:
$$D(i\mathsf c\omega)=-10\mathsf w(d\omega)+6\mathsf c(\delta_\nabla\omega),$$
$$D(i\mathsf w\omega\Gamma)=[-14\mathsf w(d\omega)+18\mathsf c(\delta_\nabla\omega)]\Gamma.$$
They include derivatives of every absorbed weight. For the opposite pair, an even Clifford input of grade six or ten has only even intersection sizes in the triple anticommutator, giving the last two grade sets in the table. This does not impose AA4's first-pair factorization on the opposite pair: nonzero Cartan jets can survive. All four derivatives have no grade-one component, while their other equations remain required. Since the invariant baseline is $bI$ in grade one and $v\equiv0$, every full solution obeys
$$\boxed{G=-(b+q_1N)I-\ell S.}\tag{AH6a}$$

In dimension fourteen, $\mathrm{Ric}^{\sharp}=G-\operatorname{tr}(G)I/12$ and $\operatorname{tr}S=5N$. Thus
$$\mathrm{Ric}^{\sharp}=\left[\frac b6+\frac{2q_1+5\ell}{12}N\right]I-\ell S.$$
[N1(d)](N-z1-neighbourhood.md) gives a nonzero Euler field $E$ with $|E|^2=c$ and $\mathrm{Ric}^{\sharp}E=-E/c$. Set $e=E/\sqrt c$ and $\mu=S(e,e)=|\iota_e\omega|^2$. Since every $\ell\ne0$, the full vector identity requires the alignment
$$Se=\mu e.$$
The orthogonal decomposition $\omega=e^\flat\wedge\alpha+\beta$ on $e^\perp$ gives $0\le\mu\le N$. Writing $L=b+6/c>0$, the exact scalar consequence is
$$\boxed{2L+(2q_1+5\ell)(N-\mu)+(2q_1-7\ell)\mu=0.}\tag{AH6b}$$
It yields precisely the following bounds for this channel:

| $a$ | required $\mu$ | lower bound on $N$ | possible $\mu/N$ |
|---|---|---|---|
| $-1/96$ | $(40N+L)/64$ | $N\ge L/24$ | $(5/8,1]$ |
| $-3/832$ | $(48N-L)/384$ | $N\ge L/48$ | $[0,1/8)$ |
| $-3/704$ | $(304N+L)/640$ | $N\ge L/336$ | $(19/40,1]$ |
| $-1/160$ | $(40N-L)/192$ | $N\ge L/40$ | $[0,5/24)$ |

In particular $N>0$. The coefficient pairs in (AH6b) are $(80,-48)$, $(-96,672)$, $(608,-672)$ and $(-80,304)$. Each pair has opposite signs, so these necessary Euler/Gram equations do not exclude the endpoint domains.

The bounds are sharp **only for the selected algebraic grade-one/Euler channel**. For any allowed pair $(N,\mu)$, the real form
$$\omega=\sqrt\mu\,e^{01234}+\sqrt{N-\mu}\,e^{56789},\qquad e=e_0,$$
has the prescribed norm and contraction, with $e$ a Gram eigenvector. Defining the synthetic tensor $G$ by (AH6a) then satisfies the full selected Euler vector equation whenever (AH6b) holds. At each lower bound use a simple form with the Euler direction inside its support when $\ell>0$, or in its annihilator when $\ell<0$. These examples are not curvatures of constructed Z1 metrics and do not solve the other nonlinear or differential equations.

As a separate consequence, every nonzero isotropic Gram $S=(5N/14)I$ is incompatible with these constant endpoint domains. In the four rows $q_1+5\ell/14$ equals $120/7$, $624/7$, $528/7$ and $200/7$, respectively. Equation (AH6a) would make $\mathrm{Ric}^{\sharp}$ a positive scalar, contradicting its negative Euler eigenvalue. This includes AH4's cyclic orbit at the four constant endpoints, without extending AH5 to general five-forms or isolated endpoints of moving fields.

## AH7. No constant normalized Gram spectrum along a fibre patch

Keep AH6's constant-odd endpoint domain $U$. Every full solution would have $N>0$ there. Fix any $x\in\pi(U)$ and any nonempty relative-open subset $\Omega\subset U\cap F_x$ of its positive fibre.

**Proposition AH7.** The unordered eigenvalue multiset of $S/N$ cannot be constant on $\Omega$.

*Proof.* Suppose this normalized spectrum is constant on $\Omega$. By AH6a, put
$$H=G+bI=-NA,\qquad A=q_1I+\ell S/N,\qquad
D=\operatorname{tr}A=14q_1+5\ell>0,\qquad \tau=\operatorname{tr}H=-DN.$$
For a real five-form, $0\le S\le NI$. Consequently the eigenvalues $\alpha_1,\ldots,\alpha_{14}$ of $A$ are positive: the four uniform lower bounds are $40/3$, $48$, $112/3$ and $8$. They are constant on $\Omega$, with repetitions allowed, and $D$ is respectively $240$, $1248$, $1056$ or $400$.

The characteristic polynomial of $A$ and Cayley–Hamilton give the homogeneous metric identity
$$\boxed{\prod_{j=1}^{14}(DH-\alpha_j\tau I)=0.}\tag{AH7a}$$
Indeed each factor on $\Omega$ is $-DN(A-\alpha_jI)$. All factors commute. No eigenvector frame or distinct-eigenvalue assumption is required.

Fixing $x$ makes the $\alpha_j$ real constants. [S4(a)](S-z1-off-section.md) makes every coordinate entry of (AH7a) rational in the fibre matrix $y$, with denominators powers of $\det y$. Clear those denominators. Since the resulting polynomial vanishes on the nonempty real open set $\Omega$, it vanishes identically, so (AH7a) holds on the entire positive fibre. Only this metric identity is extended; no field, form, eigenvector, norm or inequality is continued.

If $\tau=0$ anywhere on that fibre, (AH7a) gives $D^{14}H^{14}=0$. The positive observer metric makes $H$ self-adjoint, so all its eigenvalues vanish and $H=0$. Then $G=-bI$ and $\mathrm{Ric}^{\sharp}=bI/6$, contradicting N1's negative Euler eigenvalue $-1/c$. Thus the continuous scalar $\tau$ is nowhere zero on the whole positive fibre. That fibre is connected, so its sign is constant. On $\Omega$ its sign is negative, since $\tau=-DN<0$.

There is a contradiction along the radial ray $y=\rho g(x)$, $\rho>0$. The metric $\rho g$ has the same connection as $g$, so [F4/K](K-z1-ricci-general-n.md) applies at these section points. Taking the scalar trace gives exactly
$$\tau(\rho)=14b+\frac{30}{c}+\frac{135}{2c_V}
-\frac{6R_g(x)}{\rho c_H}.$$
Its limit as $\rho\to\infty$ is strictly positive. This contradicts the negative sign on the connected fibre. The sign is carried by continuity and the zero-trace barrier, not by analytic continuation of an inequality. ∎

The normalized spectrum may depend on the base point $x$; constancy on just one nonempty fibre-open patch is already impossible. In particular, a fixed nonzero five-form's $SO(14)$ and nonzero-scaling orbit has constant normalized Gram spectrum and cannot describe such a patch. This includes zero/simple-form domains (zeros already violate AH6) and isotropic-Gram families. AH6 separately excludes isotropic values pointwise inside a constant endpoint domain. AH7 does not exclude isolated simple-form values in a field with generally varying normalized spectrum, and does not assert a nonzero spectral derivative at every point.

These conclusions do not affect AH6's synthetic Gram/Euler examples: those establish only selected algebraic feasibility, not metric realizability on a fibre-open domain. AH9 below excludes the two first-pair constant endpoints; AN subsequently excludes the opposite pair using Bianchi rather than a fixed-spectrum assumption.

## AH8. Ricci inertia and the local-parallel set

Keep all of AH6's constant-odd endpoint-domain hypotheses: positive Z1 with constant $c_H,c_V>0$, $\lambda<1/4$ and $c=4c_V(1-4\lambda)>0$, the fixed real $W_A$ and compact real form, $\kappa\ne0$, $\det M\ne0$ and the declared single kernel. No simplicity or constant normalized spectrum is assumed. Let $n_{\le0}$ count all zero or negative eigenvalues of the self-adjoint raised observer Ricci tensor, with multiplicity.

**Proposition AH8.** Every full solution on such a domain satisfies:

| constant endpoint $a$ | strict trace-budget limit | upper bound on $n_{\le0}$ |
|---|---:|---:|
| $-1/96$ | $8$ | $7$ |
| $-3/832$ | $72/7$ | $10$ |
| $-3/704$ | $200/19$ | $10$ |
| $-1/160$ | $216/19$ | $11$ |

These also bound the strictly negative index. They are necessary bounds; no attainment by five-form Grams or observer metrics is asserted.

*Proof.* AH6 gives $N>0$, $0\le S\le NI$, $\operatorname{tr}S=5N$ and
$$\mathrm{Ric}^{\sharp}=\xi I-\ell S,\qquad
\xi=\frac b6+\frac{2q_1+5\ell}{12}N,\qquad b>0.$$
For $\ell>0$, $\xi>0$. Each nonpositive Ricci eigenvalue corresponds to an $S$ eigenvalue at least $\xi/\ell$. Simultaneous orthogonal diagonalization and the trace budget therefore give
$$n_{\le0}\le\frac{5\ell N}{\xi}
<\frac{60\ell}{2q_1+5\ell}.$$
For $\ell<0$, instead put
$$S'=NI-S\ge0,\quad\operatorname{tr}S'=9N,\quad
\zeta=\frac b6+\frac{2q_1-7\ell}{12}N>0.$$
Then $\mathrm{Ric}^{\sharp}=\zeta I-|\ell|S'$, giving
$$n_{\le0}\le\frac{9|\ell|N}{\zeta}
<\frac{108|\ell|}{2q_1-7\ell}.$$
The final inequalities are strict because $b>0$, including when some Ricci eigenvalues are exactly zero. Substitution of AH6's four coefficient pairs gives the table. In particular the first bound permits at most seven, rather than eight, nonpositive directions. ∎

Now let $\Pi$ be [S2's parallel set](S-z1-off-section.md): $(x,y)\in\Pi$ when $y$ extends to a positive $\aleph$-parallel metric near $x$, equivalently it is fixed by the relevant local holonomy. Pointwise annihilation by the curvature endomorphisms alone is not this definition. [S3](S-z1-off-section.md) transfers [F4/K](K-z1-ricci-general-n.md) through the parallel metric, whose connection and covariant base Ricci tensor agree with those of $g$. At a point of $\Pi$,
$$\mathrm{Ric}^Y_{HH}=\mathrm{Ric}_g-\frac{c_H}{c}y,\qquad
\mathrm{Ric}^Y_{HV}=0,\qquad
\mathrm{Ric}^Y_{VV}=-\frac54G_{1/5}|_y.$$
The vertical block is negative definite in all ten directions. Relative to the positive observer metric its trace-free nine-space has eigenvalue $-5/(4c_V)$ and its Euler line eigenvalue $-1/c$. These signs hold for every $\lambda<1/4$. The orthogonal block decomposition adds the Ricci indices, so at a possible domain point in $\Pi$:

- The $a=-1/96$ row is impossible: ten negative vertical directions already exceed seven. Its entire domain must avoid $\Pi$.
- The $a=-3/832,-3/704$ rows require $\mathrm{Ric}_g-(c_H/c)y$ to be positive definite; there is no remaining allowance for even a zero horizontal eigenvalue.
- The $a=-1/160$ row allows at most one nonpositive horizontal direction.

In particular, all four domains avoid $\Pi$ above every base point with $\operatorname{rank}\mathrm{Ric}_g\le2$. The kernel of the base Ricci tensor then has dimension at least two, and on it $\mathrm{Ric}_g-(c_H/c)y$ is strictly negative because $c_H/c>0$ and $y$ is positive. Together with the ten vertical directions this gives at least twelve negative directions. This argument needs no simultaneous diagonalization of $y$ and $\mathrm{Ric}_g$.

AH8 uses pointwise trace budgets and the geometric blocks only where the domain meets $\Pi$. It continues no inertia inequality or field to other fibre points. AH9 and AN supply stronger obstructions excluding all four constant-endpoint domains. Isolated zero-$u$ points of varying odd fields remain outside the constant-domain hypothesis. [Review40](../../reports/reviews/REVIEW40-GU-AH8-2026-11-02.md) records the original internal derivation and adversarial checks.

## AH9. The first-pair constant endpoint domains are empty

Keep AH6's standing constant-odd endpoint-domain hypotheses, for either $a=-1/96$ with $k=i\mathsf c(\omega)$ or $a=-3/832$ with $k=i\mathsf w(\omega)\Gamma$. The chart weight is absorbed into an arbitrary smooth real five-form. Then $v=D(v)=0$. AA4 retains every form and weight derivative in $D(k)$, whose grades are respectively $3,7$ or $7,11$. The curvature and constant invariant baseline have grade one, and the even equation gives $L_{T_o}k=0$. Hence every full solution requires
$$[Q(k)]_5=[Q(k)]_9=0.$$

Write $\Omega$ for the Clifford image of $\omega$, $P_j$ for grade projection and $F_j=P_j(\Omega^2)$. The following identities hold for every real five-form and every orthonormal generator $\gamma_a$:

| compact copy | $[Q(k)_a]_5/i$ | $[Q(k)_a]_9/i$ |
|---|---|---|
| $i\mathsf c(\omega)$ | $\frac{28}{3}P_5(\Omega\gamma_a\Omega)-\frac83P_5(\gamma_a\Omega^2)$ | $\frac{40}{3}P_9(\Omega\gamma_a\Omega)-24P_9(\gamma_a\Omega^2)$ |
| $i\mathsf w(\omega)\Gamma$ | $-36P_5(\Omega\gamma_a\Omega)+\frac{112}{3}P_5(\gamma_a\Omega^2)$ | $-40P_9(\Omega\gamma_a\Omega)-\frac{16}{3}P_9(\gamma_a\Omega^2)$ |

**Quadratic basis proof of these identities.** Both sides are $SO(14)$-equivariant quadratic maps. For $\omega=\alpha e^{01234}+\beta e^J$ with symbolic independent real $\alpha,\beta$, the helper checks every output slot for intersections $|\{0,1,2,3,4\}\cap J|=0,\ldots,5$. Any pair of coordinate five-blades occupies at most ten axes. A signed coordinate permutation maps it to the representative of its intersection type; flip an unused axis if needed to repair its determinant without changing either input blade. Input reordering signs are absorbed by $\alpha,\beta$. Thus the representative checks cover every basis pair under $SO(14)$, including the diagonal case. Polarization then proves the complete quadratic identities on $\Lambda^5\mathbb R^{14}$. This is an exhaustive basis certificate, with no choice of eigenvectors or reflection convention required.

For two distinct coordinate five-blades, the exact nonzero-slot counts for either copy are:

| overlap | grade five | grade nine |
|---|---:|---:|
| $0$ | $0$ | $10$ |
| $1$ | $0$ | $6$ |
| $2$ | $6$ | $0$ |
| $3$ | $10$ | $0$ |
| $4$ | $0$ | $0$ |

All nonzero coefficients are rational multiples of $\alpha\beta$. The universal zero-locus argument follows from the identities, rather than from this pair catalogue alone.

**Scalar square and vector sandwich.** Clifford reversion fixes a real grade-five $\Omega$. Consequently $\Omega^2$ has only grades $0,4,8$, with scalar part $N=|\omega|^2$, and $\Omega\gamma_a\Omega$ has only grades $1,5,9$. The standard contraction $\sum_a\gamma_a\Omega\gamma_a=-4\Omega$ and the square-term contractions $10F_4$, $6F_8$ give

| compact copy | $\sum_a\gamma_a[Q(k)_a]_5/i$ | $\sum_a\gamma_a[Q(k)_a]_9/i$ |
|---|---|---|
| $i\mathsf c(\omega)$ | $-64F_4$ | $-\frac{592}{3}F_8$ |
| $i\mathsf w(\omega)\Gamma$ | $\frac{1552}{3}F_4$ | $128F_8$ |

All four constants are nonzero. If both higher-grade equations vanish, then $F_4=F_8=0$ and $\Omega^2=N$. For $\omega\ne0$, positivity gives $N>0$ and $\Omega^{-1}=\Omega/N$. The nonzero sandwich coefficients in the first table then force its grades five and nine to vanish. Reversion leaves only grade one, so $T(v)=\Omega v\Omega^{-1}$ maps the vector space into itself.

**Simplicity lemma.** This map $T$ preserves Clifford squares, hence is orthogonal; $\Omega^2=N$ also gives $T^2=I$. Choose an orthonormal eigenbasis with eigenvalues $\epsilon_a\in\{1,-1\}$. The relation $\Omega\gamma_a=\epsilon_a\gamma_a\Omega$ holds coefficientwise in the grade-five monomial expansion. Each such monomial commutes with $\gamma_a$ exactly when its five-set contains $a$, and anticommutes otherwise. Every nonzero coefficient must therefore have support equal to the entire $+1$ eigenset. This set has size five, since $\Omega$ is nonzero and homogeneous of grade five; only its one five-blade can occur. Thus $\omega$ is simple. Conversely a simple five-form has scalar square and vector sandwich, so both displayed higher-grade equations vanish. Their common real zero locus is precisely zero or simple five-forms.

The full sandwich condition matters. For example,
$$\omega=e^{01}\wedge(e^{234}+e^{567}+e^{8,9,10}+e^{11,12,13})$$
has $N=4$ and $S=\operatorname{diag}(4,4,1,\ldots,1)$ of rank fourteen. It has scalar square $\Omega^2=4$ and $\omega\wedge\iota_X\omega=0$ for every $X$, yet is not simple. Both copies have $Q_9=0$ but a nonzero $Q_5$ coefficient: at external index $2$ and Clifford indices $(3,4,5,6,7)$, $Q_5/i$ is $-56/3$ and $72$, respectively. Thus neither the weaker wedge condition nor scalar square alone establishes simplicity.

**Endpoint exclusion.** A full solution has $N>0$ by AH6. The preceding proof would make $\omega$ nonzero simple everywhere, with constant normalized Gram spectrum (five ones and nine zeros). Every nonempty observer domain contains a nonempty relative-open fibre patch. AH7 excludes that spectrum there. Hence neither first-pair constant endpoint admits a full solution, for any $\lambda<1/4$ or either sign of $\kappa$ under the standing hypotheses.

This proof uses necessary higher-grade equations and then AH7. For the opposite pair, $D(k)$ can have grades five or nine, so this separation does not apply; AN instead uses its determined derivative and contracted Bianchi. Moving-odd endpoints, nonzero-$u$ conic fields, coincident kernels, $\det M=0$ and indefinite real forms remain outside these constant-domain exclusions. [Review41](../../reports/reviews/REVIEW41-GU-AH9-2026-11-02.md) records AH9's independent proof audit and shared-engine verification boundary.

## Reproduction and remaining boundary

Run `python -X utf8 -I computation/geometric-unity/AH_lambda5_gram.py` from the repository root. The original helper checks both compact phases, full symbolic-weight formulas on general forms, off-diagonal and full-support Grams, mixed and derivative controls, both exact kernel charts and the negative conic values. Its cyclic channel control uses the exact Gram calculation. The regular channel parameters and remaining block regularity are checked exactly.

Run `python -X utf8 -I computation/geometric-unity/AH_cyclic_obstruction.py` for AH5. It computes the three raw integer polarizations of the explicit cyclic form, then checks every grade-five/nine control slot with the actual compact phase. The grade-nine positive quadratic supplies the real nonzero-weight obstruction; the displayed polynomial ideal is an additional exact check.

AH1 and AH6's all-form stress identities follow from the contraction count, and the rank conclusions from positive Gram linear algebra. AH9's higher-channel identities use the complete quadratic basis-orbit certificate proved above. The additional general-form controls do not classify arbitrary five-form orbits. Both helpers share T's accepted Clifford engine and supply no independent physical evidence.

The Gram helper preserves all 94 controls through AH8 and adds 50 AH9 controls. The symbolic all-slot orbit cover checks both copies at overlaps zero through five, including the diagonal terms; general forms, exact gamma contractions, separate derivative grades and a full-support weak-condition counterexample check the conventions. The elementary simplicity proof supplies the scalar-square/vector-sandwich implication. The recorded run passes 144 controls in 11.5 seconds with Python 3.12 and SymPy 1.14.0. It uses F4, S3–S4 and N1 as geometric inputs; it does not reconstruct the metric engine. Finite controls do not establish metric or field existence. The separately recorded cyclic-obstruction run passes all 14 checks in 5.5 seconds.

AH1–AH5 retain their stated scopes. AH6–AH9 and AN exclude all four constant-odd endpoints; [AL](AL-opposite-lambda5-system.md) preserves the opposite-endpoint jet equations. AP excludes the entire first single conic, and AR excludes all opposite zero-$w$ contacts. The opposite nonzero-$w$ system, its critical points, coincident kernels, $\det M=0$, other operators, non-invariant odd fields and indefinite real forms remain open or outside these conclusions. AS retains the actual form/Gram equations; its traced constraints alone do not solve the system.
