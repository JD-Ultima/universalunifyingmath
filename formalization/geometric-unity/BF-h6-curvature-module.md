# BF. A 121-dimensional submodule in the H6 curvature kernel

On the flat positive Z1 model at $\lambda=1/5$, the coefficient-curvature map on $H_6$ has a real 121-dimensional submodule in its kernel. It is the $SO(5)$ orbit span of BC's explicit closed tensor. This identifies a substantial space left by curvature compatibility; it does not classify the whole kernel or solve the nonlinear field equation.

Evidence: [Casimir helper](../../computation/geometric-unity/BF_h6_casimir.py), [output](../../computation/geometric-unity/BF_h6_casimir.out) and [record](../../metrics/geometric-unity/BF-h6-curvature-module.yaml). The helper imports BC's tensor actions as a library and checks every component of the spectral identities. The representation argument below is a hand proof. Internal [Review68](../../reports/reviews/REVIEW68-GU-BF-2026-11-05.md) and [Review69](../../reports/reviews/REVIEW69-GU-BF-2026-11-05.md) accept BF1–BF3 and disclose prior contributions; external peer review remains outstanding.

## Setting

Use a flat positive four-base and reconstructed positive Z1 with $\lambda=1/5$, $c_H,c_V>0$. At normalized scales $c_H=2,c_V=1$, [AM](AM-global-flat-models.md) and [AZ1](AZ-first-hook-classification.md) identify the observer with $SL(5,\mathbb R)/SO(5)$ and its tangent with
$$V=\operatorname{Sym}^2_0(\mathbb R^5),\qquad \dim V=14.$$
Set $K=SO(5)$ and let $L=SO(4)$ fix the last unit coordinate vector in $\mathbb R^5$. $L$ is also the isotropy of the affine transitive group used in [BC](BC-affine-hooks.md).

All hook slots are covariant:
$$H_6=\{B\in V^*\otimes\Lambda^6V^*: \operatorname{Alt}B=\operatorname{Ctr}B=0\},\qquad \dim H_6=36608.$$
The coefficient-curvature map is
$$\mathcal R:H_6\longrightarrow V^*\otimes\Lambda^8V^*,\qquad B\longmapsto F^{T^*Y}\wedge B.$$
This is the curvature in $d_\nabla^2$ on the coefficient covector, rather than the full tensor-curvature action on a parallel field. At the symmetric model it is $K$-equivariant.

BF is a statement about this positive metric and real tensor representation. In the pure $H_6\Gamma$ GU branch, [BA](BA-hook-closure.md) makes $\mathcal RB=0$ necessary. The geometric kernel bound alone requires no choice of $\kappa$ or claim that a tensor solves GU.

## BF1. Exact Casimir decomposition of the explicit closed value

Let $B_*$ be BC's nonzero real $L$-fixed hook tensor, given by the [1,358-entry certificate](../../computation/geometric-unity/BC_h6_closed_generator.json). Its affine-invariant extension has $d_\nabla B_*=0$. Therefore
$$\mathcal RB_*=d_\nabla^2B_*=0.$$
Since $\ker\mathcal R$ is $K$-invariant, it contains the entire $K$ orbit span of $B_*$.

For $J_{ab}=E_{ab}-E_{ba}$ in $\mathfrak{so}(5)$, define the positive Casimir on the covariant hook representation by
$$\mathcal C=-\sum_{a<b}\rho(J_{ab})^2.$$
The helper calibrates the standard generators against AZ's metric isometry $\Phi:V\to\operatorname{Sym}^2_0(\mathbb R^5)$. All 140 tangent columns verify the six stabilizer actions and
$$\Phi(\Lambda_{H_i}X)=\tfrac12[J_{i4},\Phi(X)],\qquad i=0,\ldots,3.$$
BC's induced action $A_i$ is minus the Nomizu map on every covariant slot. On an $L$-fixed tensor the stabilizer terms vanish, so
$$\mathcal C=-4\sum_{i=0}^3 A_i^2.$$
Each computed Casimir power is checked to remain $L$-fixed and hook-valued. The sign is independently checked by
$$\langle\mathcal CB_*,B_*\rangle=4\sum_i\|A_iB_*\|^2.$$

Exact arithmetic gives rank two for the Gram matrices on $(B_*,\mathcal CB_*)$ and $(B_*,\mathcal CB_*,\mathcal C^2B_*)$, with positive first determinant. Every tensor component satisfies
$$\mathcal C^2B_*-58\mathcal CB_*+720B_*=0.$$
Thus the minimal polynomial on this Casimir cyclic space is exactly $(c-18)(c-40)$. In particular,
$$v_{18}=\frac{40B_*-\mathcal CB_*}{22},\qquad v_{40}=\frac{\mathcal CB_*-18B_*}{22}$$
are nonzero, orthogonal, real $L$-fixed vectors in $\ker\mathcal R$. The helper checks their positive norms and the full identities $\mathcal Cv_{18}=18v_{18}$ and $\mathcal Cv_{40}=40v_{40}$.

## BF2. From two fixed vectors to an orbit span of dimension 121

For either $v=v_\eta$, $\eta\in\{18,40\}$, put $W_v=\operatorname{span}_{\mathbb R}(K\cdot v)$. The map
$$T:W_v\longrightarrow C^\infty(K/L),\qquad T(w)(kL)=\langle w,k\cdot v\rangle$$
is well-defined because $v$ is $L$-fixed. It intertwines the $K$ actions and is injective: if $T(w)=0$, then $w$ is perpendicular to its whole generating orbit, hence to $W_v$ itself. Centrality of the Casimir gives $\mathcal Cw=\eta w$ on $W_v$.

Here $K/L$ is the unit sphere $S^4$. With the displayed standard rotation generators, the function Casimir is $-\Delta_{S^4}$. For example it acts by four on the coordinate functions. The polar-coordinate formula for the five-dimensional Laplacian shows that restrictions of harmonic homogeneous polynomials of degree $\ell$ have eigenvalue $\ell(\ell+3)$.

These are the whole smooth eigenspaces: polynomial restrictions are dense by the real Stone–Weierstrass theorem, since they form a unital algebra separating sphere points. Decomposing polynomials into harmonic parts times powers of the radius, and using integration by parts on the compact sphere, makes an eigenfunction orthogonal to every different harmonic degree. Subtracting its projection to the matching degree leaves a function orthogonal to all polynomial restrictions and hence zero.

The real harmonic space $\mathcal H_\ell(S^4)$ is irreducible under $K$. A short proof fixes the north pole and writes $z\in\mathbb R^4$, $t\in\mathbb R$. An $L$-invariant homogeneous polynomial has the form
$$\sum_j a_jt^{\ell-2j}|z|^{2j}.$$
Harmonicity recursively determines all coefficients from $a_0$:
$$a_{j+1}=-\frac{(\ell-2j)(\ell-2j-1)}{4(j+1)(j+2)}a_j.$$
Thus the fixed harmonic subspace is one-dimensional, with nonzero evaluation at the pole. If $P$ is the orthogonal projection onto a $K$-invariant subspace of $\mathcal H_\ell$, it commutes with $K$. Its action on the pole's evaluation representer must be scalar. Transitivity gives the same scalar on every evaluation representer; these span $\mathcal H_\ell$. Consequently $P$ is a scalar identity, and $P^2=P$ forces $P=0$ or $I$.

The nonzero invariant image of $T$ is therefore the entire matching harmonic space. Eigenvalues 18 and 40 select degrees 3 and 5, giving
$$\dim W_{v_{18}}=\binom74-\binom54=30,\qquad \dim W_{v_{40}}=\binom94-\binom74=91.$$
The dimension formula follows by subtracting homogeneous polynomials of degree $\ell-2$ from those of degree $\ell$ in five variables, using their harmonic decomposition.

The $K$ orbit span of $B_*$ contains both spectral projections, since it is stable under the Lie algebra and Casimir. Conversely $B_*=v_{18}+v_{40}$. The two distinct Casimir eigenspaces have zero intersection, so
$$\boxed{\dim\operatorname{span}_{\mathbb R}(SO(5)\cdot B_*)=30+91=121,\qquad \dim\ker\mathcal R\ge121.}$$
Equivalently, $\operatorname{rank}\mathcal R\le36487$. This is an upper bound, not a computed rank of the full map.

## BF3. Scope and the remaining field problem

Homogeneity transports the contained module to every point of the flat model. BC's affine-conjugating horizontal dilation and positive homothety preserve the hook and curvature-kernel conditions, so the same dimension bound holds for all $c_H,c_V>0$ at this fixed $\lambda=1/5$.

The Casimir cyclic space above has dimension two; its $SO(5)$ orbit span has dimension 121. Neither is asserted to be the whole curvature kernel. BC's nonlinear anisotropy excludes amplitudes and rotations of the original $B_*$, but does not exclude arbitrary linear combinations in this larger module.

The calculation is pointwise. Applying the differential Casimir of Killing fields to a closed field preserves closure, and at their common fixed origin it gives the representation operations above. Away from that origin the Killing fields also differentiate spatial dependence. The displayed minimal polynomial is therefore not automatically a global field identity, nor does it construct global Casimir eigenfields or solve the remaining nonlinear equations.

The exact helper checks the geometric normalization and spectral identities; the matrix-coefficient and harmonic arguments supply the module dimension. No full ambient rank calculation, curved-base or other-parameter classification, physical conclusion or literature priority is claimed.
