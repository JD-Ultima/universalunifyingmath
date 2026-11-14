# CL. Constant-norm closed H6 fields of arbitrarily high frequency

On the positive observer \(Y=SL(5,\mathbb R)/SO(5)\), there are global real analytic \(H_6\) hooks of positive constant norm, covariantly closed and compatible with the full coefficient curvature, with arbitrarily large pointwise covariant gradients. They are constructed from a 52-dimensional subspace of the BR triangular closure kernel. This construction shows that a gradient upper bound for the full equation must use additional equations beyond the hooks, linear closure, and constant norm.

The [BV](BV-h6-triangular-positive-obstruction.md) receiver and
the [BW](BW-h6-triangular-orbit-averages.md) diagonal limit exclude
these examples and their bounded rotated superpositions from the
full nonlinear equation. The result complements the
[CE](CE-h6-oscillating-triangular-fields.md) small-frequency
obstruction by treating a narrower space of amplitudes at
every real frequency.

## Character and linear carriers

Let \(S=AN\) be the real positive upper triangular subgroup and write \(sK\) for its simply transitive parametrization of \(Y\). For a tensor \(b\) in the physical BC coframe at the identity, let \(B_b(sK)=s\cdot b\). Let \(W_0\) be the BR hook kernel of the full AN exterior covariant derivative. For \(b\in W_0\), \(B_b\) is globally analytic and \(d^\nabla B_b=0\).

Define the real additive character
\[
\chi(s)=-\tfrac14\sum_{j=0}^3\log s_{jj}+\log s_{44}.
\]
In the symmetric positive matrix model \(sK\leftrightarrow ss^T\) the infinitesimal upper triangular lift \(A\) of the tangent \(\Phi\) is \(A+A^T=\Phi\). The diagonal of \(\Phi_{13}\) is \((-2/5,-2/5,-2/5,-2/5,8/5)\). The exact character differential is therefore \(e^{13}\) at the identity and the corresponding left-invariant covector everywhere. In the normalization used here \(\|d\chi\|^2=5/16\). The verifier computes this normalization from every \(\Phi_a\) and checks that \(\chi\) annihilates the infinitesimal upper triangular commutator.

Let
\[
U_\chi=\{b\in W_0:e^{13}\wedge b=0\}.
\]
A tensor lies in the wedge kernel precisely when every form blade with nonzero coefficient contains axis 13. Retain only those columns in the BR hook/closure matrix. This operation commutes with the 16 determinant-one diagonal signs in the defining representation. The exact dimensions are 2 for the trivial character, 3 in each of ten characters of weight two, and 4 in each of five characters of weight four. The verifier computes each sparse nullspace over \(\mathbb Q\), verifies every resulting basis vector against the exact row equations, and recomputes the rank independently modulo 65537. Hence \(\dim_\mathbb R U_\chi=52\).

Two selected independent real carriers have 222 and 225 covariant coefficients. They pass every hook, closure, wedge, and curvature constraint exactly. Their Gram matrix is diagonal with entries \(1026505/524288\) and \(1710915/32768\). Dividing each by the positive square root of its squared norm gives orthonormal carriers \(u,v\).

## Analytic constant-norm profiles

For every smooth real function \(f\) on \(\mathbb R\) and every \(b\in U_\chi\),
\[
d^\nabla(f(\chi) B_b)=f'(\chi)d\chi\wedge B_b+f(\chi)d^\nabla B_b=0.
\]
All hook identities are pointwise linear. In particular, for every real \(\omega\),
\[
B_\omega(sK)=\cos(\omega\chi(s))B_u(sK)+\sin(\omega\chi(s))B_v(sK)
\]
is a global real analytic unit-norm solution of the full linear closure and hook system. At the identity its covariant derivative is \(\nabla B_u+\omega\,d\chi\otimes v\), with \(\|d\chi\otimes v\|=\sqrt5/4\). In fact the covariant derivative has a uniform lower bound on all of \(Y\),
\[
\|\nabla B_\omega\|\ge \tfrac{\sqrt5}{4}|\omega|-\|\nabla B_u\|-\|\nabla B_v\|.
\]
The invariant derivatives of \(B_u\) and \(B_v\) have constant norm. For every fixed frequency the field and its gradient are bounded, and there is no frequency-independent gradient upper bound coming from the linear system and norm. The pointwise values in the triangular frame remain in \(W_0\) for every frequency.

## Rotated bounded profiles and the full equation

Let \(F_y\) be BV's exact homogeneous quadratic receiver, positive definite on \(W_0\). It vanishes at full nonlinear solutions. The following argument allows finite real sums of rotated carrier profiles, any bounded smooth profiles if the full equations are imposed globally, and finite-variation superpositions under the stated boundedness condition.

Consider a real field represented by an integral of
\(k\cdot(f(\chi)B_b)\), with \(k\in SO(5)\), \(b\in U_\chi\), and
jointly measurable bounded profiles \(f\). Assume the integral of
\(\|b\|\|f\|_\infty\) against the variation of the representing
measure is finite. This includes finite sums of bounded profiles
and finite-frequency trigonometric profiles. The argument below
only uses the global representation and its integrable bound.

Choose a generic BW diagonal ray by Haar Fubini using the integrable variation weights on the rotations. In the physical moving frame the integrand has tensor factor \(Q(t,k)\cdot b\), where the QR compact factor tends almost everywhere to \(uD(k)\) for a single \(u\in SO(5)\) and an even diagonal sign \(D(k)\). The scalar factor at the inverse triangular QR component stays bounded by \(\|f\|_\infty\); no phase limit is needed. After transport by \(u^{-1}\), replace \(u^{-1}Q(t,k)\cdot b\) by \(D(k)\cdot b\). The replacement integral takes values in \(W_0\), since the even sign subgroup preserves \(W_0\). Its norm difference from the true field tends to zero by dominated convergence, using the bound \(2\|b\|\|f\|_\infty\).

A full solution would have positive constant norm \(r\) and \(F_y=0\). The replacement sequence has norm tending to \(r\), so \(F_y\) is bounded below there by a positive constant. The norm difference just proved contradicts \(F_y=0\). The same argument applies when a finite sum of analytic rotating profiles satisfies the full equation on a nonempty open set, since its analytic residual then vanishes globally. Other superpositions with known analyticity have that same local-to-global implication.

This proof covers arbitrary frequencies because the amplitudes already lie in \(W_0\). It does not cover arbitrary amplitudes in the larger complex kernels \(W(\zeta)\) of CE, which can contain vectors outside \(W_0\). Nor does it assert analytic continuation for merely smooth profiles which solve an equation only on an open set.

## Verification

Run the [verifier](../../computation/geometric-unity/CL_h6_oscillating_carriers.py)
with the pinned dependencies. The default checks the exact Cartan
differential, both carrier tensors and their Gram, ambient
curvature, and one rank block. The full command checks all 16.

```text
python -X utf8 computation/geometric-unity/CL_h6_oscillating_carriers.py
python -X utf8 computation/geometric-unity/CL_h6_oscillating_carriers.py --full
```

[Review101](../../reports/reviews/REVIEW101-GU-CL-2026-11-14.md)
records the exact checks. External peer review remains outstanding.
