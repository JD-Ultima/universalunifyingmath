# AF. Endpoint reduction on the first exterior conic

The symmetric grade-thirteen equation rules out every nonzero-u point on the single $\Lambda^3:\mathsf c,\mathsf w\Gamma$ kernel branch. Together with AA7's pointwise three-form classification, it reduces this branch to two constant odd endpoints. The remaining even-field equations at those endpoints are not solved here.

Code: [AF_symmetric_gamma.py](../../computation/geometric-unity/AF_symmetric_gamma.py) and [output](../../computation/geometric-unity/AF_symmetric_gamma.out). Record: metrics/geometric-unity/AF-symmetric-gamma.yaml. [Independent internal Review28](../../reports/reviews/REVIEW28-GU-AF-2026-11-02.md) accepts AF1–AF3 and reproduces the ten controls in 2.6 seconds. This is internal review, not external peer review. The proof uses [AA](AA-conic-traces.md), rather than a classification of arbitrary singular fields.

## Setting and conversion

Keep the true dualized equation, compact real $(14,0)$, gauge algebra $\mathfrak u(128)$, fixed parallel real $W_A$, constant $\kappa\ne0$ and $\det M\ne0$. The invariant odd field and slaved vector part are
$$T_o=i\kappa(a\Phi^1+w\Phi^1\Gamma),\qquad
v=P_2(C)+P_{14}(B_0),\qquad
C=\frac{24}{\delta}\,da,\quad B_0=208wC,$$
where
$$w^2=f(a)=\frac{(1+\frac{112}{3}a)(1+400a)}{11264/3},\qquad
\delta=1+224a-3328f(a).$$
All other non-vector singular factors are excluded. Thus the only extra even field is
$$k=r\mathsf c(\omega)+s\mathsf w(\omega)\Gamma$$
for a real three-form and a smooth nonzero real kernel chart $(r,s)$. The weights obey
$$(1+\tfrac{112}{3}a)r=\tfrac{352}{3}ws,\qquad
(1+400a)s=32wr.$$
In particular both weights are nonzero wherever $w\ne0$.

AA6 and AA7 imply that every full solution has $\omega=0$ or a decomposable three-form at each point. Where it is nonzero, let $P$ be the orthogonal projector onto its three-dimensional covector support, and put
$$\nu=rs|\omega|^2.$$
At a zero of $\omega$, define $\nu P=0$. No smooth choice of a global three-frame is needed.

Convert the Clifford grade-thirteen output by **right multiplication by $\Gamma$**, then divide its grade-one coefficients by $i$. Since $\Gamma^2=-1$, this conversion has the opposite sign from its inverse. Identify its symmetric part with a symmetric covariant two-tensor using $g_Y$.

## The symmetric equation

**Proposition AF1.** Every full solution in the stated branch satisfies
$$\boxed{\left[-\kappa^2w(1+208a)+3328w|C|^2+40\nu\right]I
=3328w\,C\otimes C+\frac{208}{3}\nu P.}\tag{AF1}$$

*Proof.* The invariant contribution $\kappa T_o+Q(T_o)$ becomes
$-\kappa^2w(1+208a)I$. The $W_A$ curvature contribution has Clifford grade one, so supplies nothing in this converted channel.

For aligned vector covectors $C=ce^0$, $B_0=be^0$, the complete differentiated cubic, with symbolic real $b,c$, gives
$$[Q(v)]_{13}\Gamma/i=16bc(I-e^0\otimes e^0).$$
Every aligned pair is an $SO(14)$ rotation and rescaling of this pair. Substitution of $B_0=208wC$ therefore gives
$$[Q(v)]_{13}\Gamma/i=3328w(|C|^2I-C\otimes C).$$
This also holds at $C=0$.

For the unit simple form $\omega=e^{012}$, the complete symbolic-weight cubic output has eigenvalue $-88rs/3$ on its three-plane and $40rs$ on the perpendicular eleven-plane, with all other entries zero. The group $SO(14)$ acts transitively on unit simple three-forms; quadratic homogeneity includes their norms. Hence
$$[Q(k)]_{13}\Gamma/i=40\nu I-\frac{208}{3}\nu P.$$
The symbolic calibration covers every $r,s$, and the orbit argument covers all simple forms. It is not extrapolation from sampled nondecomposable forms.

The symmetric converted derivative of either vector lift is zero. Indeed an equivariant linear map from a vector derivative matrix $H$ to its symmetric two-tensor output has form $\alpha\,\operatorname{Sym}H+\beta\,\operatorname{tr}(H)I$. Complete coordinate-jet outputs at $H_{00}=1$ and $H_{10}=1$ give zero symmetric matrices, forcing $\alpha=\beta=0$. The coefficients are constant because the lifts and operator are parallel. Derivatives of $B_0$, including its variable weight, are included.

The mixed vector/kernel term and the kernel derivative also have zero symmetric part in this channel. Their input type is $V^*\otimes\Lambda^3V^*$. A metric contraction to two free output indices can only contract the vector index against one three-form index, leaving an antisymmetric pair. Every other pairing traces two indices of the alternating input and vanishes. An epsilon contribution cannot survive: there are six input/output tensor slots, below dimension fourteen, and contracting unused epsilon indices with the symmetric metric is zero. Equivalently the input splits as $\Lambda^4\oplus\Lambda^2\oplus H_3$, with no symmetric two-tensor summand. This includes all derivatives of $r,s,\omega$, consistently with AA4.

These are all terms in the untraced odd equation. Summing them proves (AF1). No curvature-trace sign or equation for the slaved baseline was assumed. ∎

## Rank obstruction and endpoints

**Theorem AF2.** A full solution in this single kernel branch has no point with $w\ne0$.

*Proof.* At any point, the right side of (AF1) has rank at most four. A nonzero scalar multiple of the fourteen-dimensional identity has rank fourteen. Thus its scalar coefficient must vanish, and
$$3328w\,C\otimes C=-\frac{208}{3}\nu P.$$
If $\nu\ne0$, then $\omega\ne0$ and the right side has rank three, whereas the left side has rank at most one. Consequently $\nu=0$. At a hypothetical point with $w\ne0$, both chart weights are nonzero; positivity of the form norm gives $\omega=0$. The last tensor equality then gives $C=0$.

The scalar coefficient in (AF1) now reduces to
$$-\kappa^2w(1+208a)=0.$$
Since $\kappa,w\ne0$, it forces $a=-1/208$. But the exact conic value is
$$f(-1/208)=-3/14872<0,$$
contradicting real $w$ and $w^2=f(a)$. ∎

**Corollary AF3.** On a connected solution domain remaining in this single branch, the odd field is constant and
$$w=0,\qquad a=-3/112\quad\text{or}\quad a=-1/400.$$
These are necessary endpoints, not asserted solutions.

*Proof.* AF2 gives $w=0$ throughout. The conic has exactly the two displayed zero-w endpoints. A continuous $a$ on a connected domain cannot switch between them. ∎

The endpoint kernel and its derivative remain in the other equations. Coincident kernels, the other Lambda3 conic, the Lambda5 conics, $\det M=0$, non-invariant odd fields and indefinite real forms are outside this result. In particular this does not classify all constant singular fields.

## Reproduction

From the repository root:

    python -X utf8 -I computation/geometric-unity/AF_symmetric_gamma.py

Ten exact controls pass in 2.7 seconds with Python 3.12 and SymPy 1.14.0. They include complete symbolic-weight and aligned-vector outputs, a non-coordinate simple three-plane, derivative and mixed controls, the right-Gamma baseline sign and the exact conic/endpoints. An arbitrary algebraic curvature diagonal in one control is test data, not a constructed observer metric.

The equivariance, contraction and rank arguments supply the universal statements. The helper shares T's accepted Clifford engine; it is not a second implementation. No new physical evidence or PDE sufficiency result is claimed.
