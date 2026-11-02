# AG. The other Lambda3 exterior conic

The second $W_A$ Lambda3 kernel branch also has no full solution at a nonzero-$u$ point. Its compact phase and nonlinear coefficients differ from AA6, while its simple-form symmetric Gamma tensor agrees with AF. A connected solution domain in this branch can therefore reach only two constant odd endpoints. Their remaining even-field equations are unresolved.

Code: `computation/geometric-unity/AG_other_lambda3.py` → `.out`. It uses T's accepted Clifford implementation, with no J/R evaluator or off-repository dependencies. Record: `metrics/geometric-unity/AG-other-lambda3.yaml`. [Independent internal Review30](../../reports/reviews/REVIEW30-GU-AG-2026-11-02.md) accepts AG1–AG4 and reproduces all 30 checks. This is internal review, not external peer review. The universal statements below use contraction, equivariance and rank arguments, rather than a sampled parameter range.

## Setting, phase and kernel charts

Keep the true dualized equation, real compact $(14,0)$, gauge algebra $\mathfrak u(128)$, fixed parallel real $W_A$, constant $\kappa\ne0$ and $\det M\ne0$. Write
$$T_o=i\kappa(a\Phi^1+w\Phi^1\Gamma),\qquad a=t/\kappa,\quad w=u/\kappa.$$
Stay on the single $\Lambda^3:w,c\Gamma$ factor, away from every other non-vector singular factor. Its conic is
$$w^2=f(a)=\frac{(1+240a)(1+\frac{272}{3}a)}{11264/3}.\tag{AG1}$$
The slaved vector field is unchanged from AA and AF:
$$v=P_2(C)+P_{14}(B_0),\qquad
C=\frac{24}{\delta}\,da,\quad B_0=208wC,\quad
\delta=1+224a-3328f(a).$$
Both vector-block factors are nonzero by $\det M\ne0$.

The two kernel copies have Clifford grades four and twelve and are Hermitian. Their real compact field therefore has the phase
$$\boxed{k=i[r\,\mathsf w(\omega)+s\,\mathsf c(\omega)\Gamma],}\tag{AG2}$$
with $\omega$ a real three-form and real chart weights $r,s$, not both zero. The lifts $\mathsf c,\mathsf w$ are those of AA4; the scalar $w$ remains distinct from the wedge lift. The kernel equations are
$$Ar=32ws,\qquad Bs=\frac{352}{3}wr,\qquad
A=1+240a,\quad B=1+\frac{272}{3}a.$$
Two smooth charts are
$$ (r,s)=(32w,A),\qquad (r,s)=(B,\tfrac{352}{3}w).\tag{AG3}$$
They cover the conic because $A,B$ cannot both vanish. In particular both weights are nonzero where $w\ne0$. At a zero-$w$ endpoint use a nonzero smooth chart in $(a,w)$; $w$ is not presumed to be a smooth function of $a$. The compact-phase kernel identities are checked symbolically modulo $w^2-f(a)$.

## AG1. Nonlinear hook equations

Use AA6's notation
$$Z_a=(\iota_{e_a}\omega)\wedge\omega,\qquad
F_4=\sum_a(\iota_{e_a}\omega)\wedge(\iota_{e_a}\omega),$$
and its grade-five hook projector
$$\mathcal H_5U=U-\frac16\mathsf c\left(\sum_a e^a\wedge U_a\right)
-\frac1{10}\mathsf w\left(\sum_a\iota_{e_a}U_a\right).$$
The conversion of grade nine is **right multiplication by $\Gamma$**, not its inverse. The full quadratic identities for (AG2) are
$$\boxed{[Q(k)]_5=\frac{8i}{3}(3s^2-59r^2)Z-8ir^2\mathsf w(F_4),}\tag{AG4}$$
$$\boxed{[Q(k)]_9\Gamma=\frac{224i}{3}rsZ.}\tag{AG5}$$
Thus their hook projections have coefficients $\frac{8i}{3}(3s^2-59r^2)$ and $\frac{224i}{3}rs$ multiplying the same $\mathcal H_5Z$.

*Proof of the universal identities.* The metric-contraction count in AA6 applies to these copies after their phase and Gamma products are included. Two alternating three-tensors and an output in $V^*\otimes\Lambda^5V^*$ have twelve tensor slots, below dimension fourteen. No epsilon pattern can add a map. A contraction within one alternating input, or between two alternating output slots, vanishes. With no cross-input contraction, the free vector output pairs with an input slot and gives $Z$. With one cross-input contraction it must pair with an alternating output slot, giving $\mathsf w(F_4)$. More cross-input contractions force a forbidden output pairing. These are the only patterns.

The exact differentiated cubic on $e^{012}+e^{345}$ and $e^{012}+2e^{034}$ fixes their coefficients polynomially in the real weights. The first has $F_4=0$ and nonzero $Z$; the second distinguishes the wedge-lift pattern through nonzero $F_4$ and hook part. With the phase $i$ in (AG2), the complete outputs give (AG4–5). Additional nondecomposable and full-support controls check the conventions. ∎

## AG2. The full equation forces a zero or simple form

Every full solution in this single branch satisfies $\mathcal H_5Z=0$, and hence has $\omega=0$ or a decomposable real three-form at each point by AA7.

*Proof.* The unsolved baseline and the vector/kernel mixed term have no $H_5$ component in these two grades. As in AA6, their relevant input types are $V\otimes V$ and $V\otimes\Lambda^3$, neither containing $H_5$.

The derivative here does not use AA4's formula for the other copies. Direct Clifford grades suffice: the $W_A$ derivative of the grade-four input is contained in grades three and seven; the grade-twelve input gives grade eleven. There is no grade-five or grade-nine output, including derivatives of the variable chart weights and of $\omega$. The lifts and operator are parallel, so this statement holds for arbitrary covariant jets.

Consequently both hook coefficients must annihilate $\mathcal H_5Z$. They have no common nonzero real weight zero: if $rs=0$, the other coefficient reduces to $8is^2$ or $-472ir^2/3$. The chart weights are not both zero, proving $\mathcal H_5Z=0$. AA7 applies in real positive dimension fourteen. ∎

This derivative assertion does not make the derivative tensor vanish. Explicit nonzero Cartan jets survive on each copy in other odd grades. Those equations remain part of the full system.

## AG3. The symmetric Gamma obstruction

For a unit simple form $\omega=e^{012}$, the symmetric part of the converted tensor $[Q(k)]_{13}\Gamma/i$ has eigenvalues $-88rs/3$ on its three-plane and $40rs$ on its perpendicular eleven-plane. SO(14) transitivity on unit simple forms and quadratic homogeneity therefore give, for every simple $\omega$,
$$\operatorname{Sym}([Q(k)]_{13}\Gamma/i)=40\nu I-\frac{208}{3}\nu P,
\qquad \nu=rs|\omega|^2,$$
where $P$ is the orthogonal projector onto its three-plane. At zero use $\nu P=0$. No smooth global factorization or frame is required.

The invariant and vector contributions are exactly those in AF1, since neither the odd field nor its slaved vector part changed. The new kernel derivative has no grade thirteen by the preceding grade proof. The mixed term has no symmetric converted component because $V\otimes\Lambda^3=\Lambda^4\oplus\Lambda^2\oplus H_3$ contains no symmetric two-tensor. Thus every full solution satisfies AF1's same equation
$$\left[-\kappa^2w(1+208a)+3328w|C|^2+40\nu\right]I
=3328wC\otimes C+\frac{208}{3}\nu P.\tag{AG6}$$
All derivatives are retained in the other odd components; they cannot contribute to this symmetric channel.

Suppose a full solution had a point with $w\ne0$. The right side of (AG6) has rank at most four, so its scalar coefficient is zero. The resulting rank-at-most-one versus rank-three equality forces $\nu=0$. Both chart weights are then nonzero, and positivity gives $\omega=0$; the tensor equality gives $C=0$. Its scalar part forces $a=-1/208$. But on the present conic
$$f(-1/208)=-1/43264<0,$$
contradicting real $w^2=f(a)$. There is therefore no full solution at a nonzero-$w$ point on this single branch. No curvature-trace sign is needed.

## AG4. The remaining necessary endpoints

On a connected solution domain remaining in this branch, the odd field is constant and
$$w=0,\qquad a=-3/272\quad\text{or}\quad a=-1/240.$$
Indeed AG3 gives $w=0$ throughout, these are the two roots of $f$, and continuity prevents switching between them on a connected domain. This does not settle the remaining even-field or differential equations at either endpoint.

## Evidence and boundary

Run `python -X utf8 -I computation/geometric-unity/AG_other_lambda3.py` from the repository root. The portable default uses exact rational and symbolic Clifford computations through T, with no off-repository files. It checks full symbolic-weight hook formulas, compact phases, both kernel charts, the coordinate and non-coordinate simple-form Gamma tensors, variable-weight derivative controls and the exact conic/endpoints. Nonzero Cartan controls on both copies confirm that the other derivative equations have not been removed. Artificial diagonal curvature coefficients in one baseline control are algebraic test data, not a constructed observer metric.

The recorded fresh run passes all 30 checks in 8.6 seconds with Python 3.12 and SymPy 1.14.0. The helper shares the accepted Clifford engine; the universal proof is not a second implementation of that engine.

The result concerns only this single real Lambda3 wedge/contraction-Gamma factor. Coincident kernels, Lambda5 branches, other Shiab members, non-invariant odd fields, $\det M=0$ and indefinite real forms are outside its scope. No endpoint solution, PDE sufficiency or classification of all singular strata is asserted.

[AI](AI-flat-lambda3.md) completes the endpoint analysis over flat-base Z1 and gives a separate nonnegative-curvature-trace exclusion. AG's sign-free endpoint statement retains its original scope.
