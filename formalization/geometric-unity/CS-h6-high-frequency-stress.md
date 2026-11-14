# CS. Exact high-frequency stress obstruction for H6 triangular profiles

On the positive \(H_6\) observer, the complete grade-one stress equation excludes nonzero full solutions assembled from purely high-frequency triangular characters on one Cartan line. This holds for arbitrary finite sums or finite-variation superpositions of their compact rotations with axial frequency \(|\omega|\ge21\), under the global or analytic scope below. A fixed narrow cone about that line is covered by a related estimate. No assumption on how close a twisted amplitude is to the zero-frequency AN kernel is needed.

## The fixed axial stress exceeds zero on high-frequency amplitudes

Let \(N=\|b\|^2\), \(E_{\xi\xi}=\|b_\xi\|^2\) and \(I_{\xi\xi}=\sum_a\|\iota_\xi b_a\|^2\), where \(\xi\) is a unit physical tangent vector. For real \(H_6\) values the BS grade-one identity gives scalarity \(2E+I=(4/7)Ng\). At the fixed unit vector dual to \(e^{13}/\|e^{13}\|\), the corresponding scalar receiver is
\[
F_{13}(b)=2E_{\xi\xi}+I_{\xi\xi}-\tfrac47N.
\]
It is zero for every full solution. Extend it Hermitianly to complex tensors by summing its values on the real and imaginary parts. This maintains the formula with Hermitian norms. The contraction and wedge norm identity gives
\[
F_{13}(b)=2E_{\xi\xi}+\tfrac37N-\|\xi^\flat\wedge b\|^2
\ \ge\ \tfrac37N-\|\xi^\flat\wedge b\|^2.
\]

As in CL, the real left-invariant frame differential of \(\chi(s)=-(1/4)\sum_{j<4}\log s_{jj}+\log s_{44}\) is \(e^{13}\), of squared norm \(5/16\). In the AN frame, an actual complex twisted amplitude of \(e^{i\omega\chi(s)}\) satisfies
\[
D_S b+i\omega e^{13}\wedge b=0
\]
in addition to the hooks. All form blades in the second term contain axis 13. Let \(P_{13}\) be the orthogonal projection on that set of output blades, independently of the external slot. The exact operator bound verified here is
\[
\|P_{13}D_S\|_{\mathrm{op}}^2\le
\frac{9233556899517}{160000000000}<58.
\]
It holds on every ambient covariant degree-six tensor. Orthogonal projection preserves the twist equation and gives
\[
F_{13}(b)\ \ge\
\left(\tfrac37-\frac{\|P_{13}D_S\|_{\mathrm{op}}^2}{\omega^2(5/16)}\right)\|b\|^2.
\]
For \(|\omega|\ge21\) this supplies one fixed receiver and a uniform positive margin exceeding \(0.009\). There is no need to assume that \(b\) belongs to the untwisted AN kernel.

## Exact metric bound and the cone

The BR closure column is the derivative on all form and coefficient slots before projection. A covariant entry mapping input \((a,m)\) to output \((b,k)\) is scaled by the square root of
\[
\frac{g_a\prod_{i\in m}g_i}{g_b\prod_{i\in k}g_i}
\]
to obtain the coefficient in unit metrics. The verifier computes conservative rational upper bounds for these square roots and then the absolute row and column sums in each sign block. The standard matrix norm inequality gives \(\|D\|_{\mathrm{op}}^2\le\max_r\sum_c|D_{rc}|\ \max_c\sum_r|D_{rc}|\). The sixteen sign blocks are orthogonal for the input and output metrics. Their maximum gives the bound above. Using all output rows also yields \(\|D_S\|_{\mathrm{op}}^2<64\).

In particular let \(\alpha\) be the real physical covector of any Cartan character and let \(\xi^\flat=e^{13}/\sqrt{5/16}\). For a nonzero \(\alpha\), set \(\widehat\alpha=\alpha/\|\alpha\|\). The twist equation implies \(\|\widehat\alpha\wedge b\|\le8\|b\|/\|\alpha\|\). If \(\widehat\alpha\) is within \(1/20\) of either \(\xi^\flat\) or \(-\xi^\flat\), and \(\|\alpha\|\ge14\), then
\[
\frac{\|\xi^\flat\wedge b\|}{\|b\|}\le\frac1{20}+\frac8{14},
\]
whose square is strictly below \(3/7\). The same \(F_{13}\) therefore has a fixed positive margin throughout these two high-frequency cones.

## Frequency averaging and scope

Consider a real field represented by a finite-variation integral over pairs of frequencies and compact rotations, with actual twisted amplitudes \(b(\zeta,k)\) and integrable norm \(\int\|b\|\,d\mu<\infty\). Assume all frequencies are in the stated axial high-pass set, or all in the two stated cones. As in CE, use a generic BW ray and a common compact transport \(u^{-1}\). The ray values differ by \(o(1)\) in norm from the Fourier integral with asymptotic measures on frequencies. This uses only dominated convergence and unit modulus of real character phases. For each frequency atom its transformed coefficient remains in the same twisted kernel, because the diagonal sign group preserves the Cartan covector and the AN closure stack.

For each bounded Hermitian form, the Cesàro energy over the ray equals the sum of its values on the frequency atoms, by the pair-measure argument of CE. This does not discard a possible diffuse frequency measure by assumption. For a full solution, constant squared norm is a positive constant \(R^2\), and the receiver \(F_{13}\) has average zero. The norm average therefore forces \(\sum_\zeta\|\nu(\{\zeta\})\|^2=R^2\). The receiver average is at least the uniform positive margin times this sum, a contradiction. The receiver is fixed in the common compact ray frame, and full scalarity eliminates it in every physical frame.

The argument excludes an everywhere-full field under the stated integrability and frequency support. It also excludes a field satisfying the full equations on a nonempty open set when the integral represents a real analytic field. Finite sums of characters and finite-variation superpositions of bounded frequency support have this analyticity by the CE kernel argument. Unbounded frequency support is allowed for the everywhere-full conclusion; boundedness of the amplitudes and phases at real points is sufficient there.

A simpler application is a periodic single axial AN-frame profile. Write \(B(sK)=s\cdot b(\chi(s))\) with a differentiable periodic \(b\). If every Fourier coefficient that does not vanish is supported at actual frequencies \(|\omega|\ge21\), then differentiating coefficients gives the twisted kernel equation. Parseval computes the averages of squared norm and \(F_{13}\), giving the same positive margin contradiction. This does not require that the field be analytic.

Finally, for any differentiable axial AN-frame profile satisfying closure, complete stress scalarity, and fixed norm \(R>0\), put \(c=\xi^\flat\wedge b\). The equation gives \(\|c'\|\le\sqrt{58}R/\sqrt{5/16}\), while stress requires \(\|c\|\ge\sqrt{3/7}R\). Consequently the normalized direction \(c/\|c\|\) has speed at most \(\sqrt{58\cdot16\cdot7/(5\cdot3)}\) per unit \(\chi\). This does not bound derivatives in the kernel of the wedge map, or cover arbitrary mixed low and high frequencies, wider frequency sets, or bounded profiles without the declared spectral assumptions.

## Arbitrary bounded profiles with a distributional Fourier gap

The axial result also holds for a bounded differentiable profile on the whole real Cartan line, without requiring periodicity or a Fourier–Stieltjes representation. Assume its tempered distributional Fourier support is contained in \(\{|\omega|\ge21\}\). The same spectral statement holds if only the wedge projection \(c=\xi^\flat\wedge b\) has this gap. Choose a nonzero Schwartz function \(\phi\) whose Fourier transform has compact support in \([-M,M]\), and set \(\phi_L(t)=\phi(t/L)\). Write \(c_L=\phi_Lc\). The Fourier transform of \(c_L\) is a convolution with support broadened by at most \(M/L\), so Plancherel gives
\[
 \|c_L'\|_2\ge(21-M/L)\|c_L\|_2
 \ge(21-M/L)\sqrt{3/7}\,R\|\phi_L\|_2.
\]
The bounded profile and the closure equation imply
\[
 \|c_L'\|_2\le \frac{\|P_{13}D_S\|}{\sqrt{5/16}}\,R\|\phi_L\|_2
   +R\|\phi_L'\|_2 .
\]
The ratio \(\|\phi_L'\|_2/\|\phi_L\|_2\) is \(O(L^{-1})\). Letting \(L\) grow yields a contradiction by the certified margin at frequency 21. Multiplication by \(\phi_L\) and differentiation are legitimate in tempered distributions. The closure bound puts \(c_L\) and its derivative in \(L^2\).

The same argument covers an arbitrary bounded differentiable profile on the entire four-dimensional real Cartan, assuming its distributional Fourier support is contained in the two stated high-frequency cones in the physical covector norm. Write \(B(sK)=s\cdot b(a)\), where \(a\) is the four-dimensional log-diagonal part of \(s\). The linear closure equation is the fixed AN-frame operator \(D_S b\) plus exterior multiplication by the Cartan derivatives. Multiply by a Schwartz cutoff \(\phi_L(a)=\phi(a/L)\) of compact Fourier support and set \(b_L=\phi_Lb\). The Fourier support broadens by \(O(L^{-1})\) in the physical frequency coordinate \(\alpha\). The Fourier closure equation becomes
\[
 (D_S+i\,\alpha\wedge)\widehat b_L=\widehat e_L,\qquad
 \|e_L\|_2=O(L^{-1})R\|\phi_L\|_2 .
\]
Here \(e_L\) is the product-rule term from \(d\phi_L\). On the broadened support the physical radial minimum is \(14-O(L^{-1})\) and the angular distance to the two fixed axis directions is at most \(1/20+O(L^{-1})\), uniformly even if the support is unbounded. Dividing by \(\|\alpha\|\) and applying the exterior norm inequality gives
\[
 \|\xi^\flat\wedge\widehat b_L\|_2\le
 \bigl(1/20+8/14+O(L^{-1})\bigr)\|\widehat b_L\|_2
 +O(1)\|\widehat e_L\|_2 .
\]
The left side is at least \(\sqrt{3/7}\|\widehat b_L\|_2\) by Plancherel and pointwise scalarity for the original field, and \(\|b_L\|_2=R\|\phi_L\|_2\). The strict cone margin gives a contradiction as \(L\) grows. No measure representation or analyticity is needed in these global single-profile conclusions.

## Verification

Run from the root of the repository. Python, NumPy, and SymPy suffice.

```text
python -B computation/geometric-unity/CS_h6_highpass.py
```

The recomputation visits every one of the 42,042 physical columns for both the ambient and projected estimates, forms the exact rational margins, and writes stdout.

[Review106](../../reports/reviews/REVIEW106-GU-CR-CS-CV-2026-11-14.md) records the independent checks. External peer review remains outstanding.
