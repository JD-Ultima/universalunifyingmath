# AM. Complete observer models over flat bases

[F7](F-zorro-reduction.md) already identifies flat Z1 locally with a symmetric space. This note gives global coordinates, complete torus quotients and descent of the known [AI3](AI-flat-lambda3.md) bosonic field. [AO](AO-global-invariant-family.md) extends the construction to Q15's whole constant invariant family. The geometry is classical: Lovrić–Min-Oo–Ruh construct this symmetric-space metric on Gaussian parameters, distinct from the Fisher metric ([primary paper](https://math.mcmaster.ca/minoo/mypapers/lmr.pdf), §§2–3). Globke–Quiroga-Barranco explicitly discuss that distinction ([2019 paper](https://arxiv.org/abs/1904.12921)). No new symmetric-space or information-geometric metric is claimed.

Use the standard flat connection on $X=\mathbb R^4$, positive fibre matrices $y$, constant $c_H,c_V>0$, and $\lambda=1/5$. Then
$$
g_Y=c_H\,dx^T y\,dx+c_V\left[\operatorname{tr}((y^{-1}dy)^2)
-\frac15(\operatorname{tr}(y^{-1}dy))^2\right].
$$

Code: [AM_global_flat_models.py](../../computation/geometric-unity/AM_global_flat_models.py), [output](../../computation/geometric-unity/AM_global_flat_models.out). Record: [metric entry](../../metrics/geometric-unity/AM-global-flat-models.yaml). [Independent internal Review42](../../reports/reviews/REVIEW42-GU-AM-2026-11-03.md) records the proof assessment; external peer review remains outstanding.

## AM1. A global complete realization

Put $s=\sqrt{c_H/(2c_V)}$, $d=\det y$, and
$$
N(x)=\begin{pmatrix}I_4&-sx\\0&1\end{pmatrix},\qquad
Q(x,y)=d^{-1/5}N(x)^T\operatorname{diag}(y,1)N(x).
$$
This is a global isometry onto
$$
\mathcal P_5=\{Q=Q^T>0:\det Q=1\}\simeq SL(5,\mathbb R)/SO(5)
$$
with metric $c_V\operatorname{tr}(Q^{-1}dQ\,Q^{-1}dQ)$, and is complete.

*Proof.* Positivity and determinant one follow directly. For arbitrary
$Q=\begin{psmallmatrix}B&b\\b^T&h\end{psmallmatrix}\in\mathcal P_5$, its Schur complement is $h-b^TB^{-1}b=1/\det B$. The smooth inverse is
$$
x=-B^{-1}b/s,\qquad y=(\det B)B.
$$
Indeed $\det y=(\det B)^5$; substitution recovers every block. At $x=0$, the metric pullback has horizontal term $2c_Vs^2dx^Ty\,dx$, the displayed vertical term and zero cross term. Translation is the congruence $Q(x+t,y)=N(t)^TQ(x,y)N(t)$, so the identity holds everywhere.

For a tangent $V$ at $Q_0$, set $A=Q_0^{-1/2}VQ_0^{-1/2}$. It is symmetric and trace free. Every initial condition has the all-time geodesic
$$
Q(t)=Q_0^{1/2}\exp(tA)Q_0^{1/2}.
$$
It stays positive with determinant one. The trace-metric connection is
$\nabla_UV=D_UV-(UQ^{-1}V+VQ^{-1}U)/2$; the curve satisfies $Q''=Q'Q^{-1}Q'$. Thus geodesic completeness, and hence completeness by Hopf–Rinow, holds. ∎

## AM2. Complete torus quotients and spin descent

Let $L\subset\mathbb R^4$ be a full translation lattice. Its free, properly discontinuous isometric action $(x,y)\mapsto(x+n,y)$ gives exactly the flat torus observer $Y(T^4)$. In the symmetric-space coordinates it is the unipotent congruence group $N(L)$. This is a lattice in the four-dimensional translation subgroup; no finite-covolume lattice in $SL(5,\mathbb R)$ is asserted. The quotient of the complete cover is complete.

All sixteen torus spin structures admit the descended Levi-Civita spin connection and Clifford-soldering form $\Phi^1$.

*Proof.* Choose constant symmetric matrices $E_A$ orthonormal for $G_{1/5}|_I$. A global oriented orthonormal frame on the cover consists of horizontal vectors $c_H^{-1/2}y^{-1/2}\partial_x$ and vertical vectors $c_V^{-1/2}y^{1/2}E_Ay^{1/2}$. It depends only on $y$, so the deck action has identity $SO(14)$ cocycle.

Every character $\chi:L\to\{\pm1\}$ defines the spin lift
$$
(x,y,r)\longmapsto(x+n,y,\chi(n)r).
$$
The cover is contractible; these characters give all $2^4$ spin choices relative to the frame. The spin connection is translation invariant and its central, constant gluing leaves its connection form unchanged. Central signs also act trivially by conjugation on $\Phi^1=\sum_a\theta^a\otimes\gamma_a$. Thus the Hermitian spinor metric, compact adjoint bundle and soldering form descend for every choice. No canonical spin choice is claimed. ∎

## AM3. Global realization of the known bosonic field

Retain AI's true dualized equation, positive $(14,0)$, compact $\mathfrak u(128)$, fixed parallel real $W_A$, and constant real $\kappa\ne0$. On the cover and every quotient above,
$$
T_o=-i\kappa\Phi^1/400,\qquad T_e=0,\qquad
c_V\kappa^2=150000/11
$$
is a smooth global solution of the full reconstructed bosonic equation $\Upsilon=0$, for either sign of $\kappa$.

*Proof.* AM2 makes the anti-Hermitian field global for every spin structure. Adding it to the Levi-Civita spin connection gives a unitary connection. The identity gauge transformation is globally available. AI3 applies in every chart: $\Phi^1$ is parallel and the entire remaining residual is
$$
i\left[\frac{15}{2c_V}-\frac{11\kappa^2}{20000}\right]\Phi^1=0.
$$
The four-by-four vector block at $t=-\kappa/400,u=0$ is
$$
M=\operatorname{diag}(\kappa,\kappa+224t,\kappa+544t,\kappa),\qquad
\det M=-99\kappa^4/625\ne0.
$$
Local covariance and central descent establish the equation globally. ∎

No fermionic solution, joint base-metric variation, finite action, physical vacuum or empirical prediction follows. The field is the existing AI3 solution, globally realized here; it is not a nonzero Lambda5 field.

## Reproduction

Run `python -X utf8 -I computation/geometric-unity/AM_global_flat_models.py` from the repository root with Python 3.12 and SymPy 1.14.0. Twelve exact controls check the full 196-entry pullback at nonzero $x$, inverse, translations, field normalization, full determinant with nonunit $\kappa$ of both signs, and a boost geodesic with constant squared speed $2c_V$. They import no metric or Clifford engine. Completeness and spin descent are the hand proofs above, not conclusions from the finite controls.
