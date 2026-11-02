# AH. Lambda5 Gram constraints and a cyclic-form obstruction

The symmetric Gamma equation restricts both single $W_A$ Lambda5 kernel branches without classifying their five-forms. At a nonzero-$u$ point, the contraction Gram matrix must have rank thirteen or fourteen. It has rank fourteen on a left conic branch; on a right branch, rank thirteen is possible only when the slaved gradient is its annihilator direction. These conditions exclude forms supported in dimension at most twelve.

An explicit form with full fourteen-dimensional support has isotropic Gram matrix and satisfies the selected symmetric channel at regular real conic parameters. AH5 examines its remaining equations and excludes its scaling/rotation orbit on the first branch at nonzero $u$. General five-forms remain unclassified.

Code: `computation/geometric-unity/AH_lambda5_gram.py` → `.out`, with the cyclic-orbit certificate in `AH_cyclic_obstruction.py` → `.out`. Record: `metrics/geometric-unity/AH-lambda5-gram.yaml`. Both helpers use T's accepted Clifford implementation and exact rational/polynomial controls, with no off-repository files. [Independent internal Review31](../../reports/reviews/REVIEW31-GU-AH-2026-11-02.md) accepts AH1–AH4 and reproduces all 25 original checks. [Review34](../../reports/reviews/REVIEW34-GU-AH-CYCLIC-2026-11-02.md) separately accepts AH5 and reproduces its 14 checks. These are internal reviews, not external peer review.

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

## Reproduction and remaining boundary

Run `python -X utf8 -I computation/geometric-unity/AH_lambda5_gram.py` from the repository root. The original helper checks both compact phases, full symbolic-weight formulas on general forms, off-diagonal and full-support Grams, mixed and derivative controls, both exact kernel charts and the negative conic values. Its cyclic channel control uses the exact Gram calculation. The regular channel parameters and remaining block regularity are checked exactly.

Run `python -X utf8 -I computation/geometric-unity/AH_cyclic_obstruction.py` for AH5. It computes the three raw integer polarizations of the explicit cyclic form, then checks every grade-five/nine control slot with the actual compact phase. The grade-nine positive quadratic supplies the real nonzero-weight obstruction; the displayed polynomial ideal is an additional exact check.

The all-form tensor identity follows from the contraction count, and the rank conclusions from positive Gram linear algebra. Finite controls check those conventions rather than classify five-forms. Both helpers share T's accepted Clifford engine and supply no independent physical evidence.

With Python 3.12 and SymPy 1.14.0, the recorded Gram run passes all 25 checks in 4.5 seconds. The cyclic-obstruction run passes all 14 checks in 5.5 seconds.

Zero-$w$ points, coincident kernels, $\det M=0$, other Shiab members, non-invariant odd fields and indefinite real forms remain outside these nonzero-$w$ conclusions. The full nonlinear Lambda5 system and its PDE sufficiency are unresolved.
