# CV. An all-frequency H6 obstruction on the sign-trivial axial kernel

There is one fixed Hermitian receiver which is strictly positive, with a uniform norm-relative margin, on the determinant-one diagonal-sign invariant amplitudes of every imaginary axial AN character. The frequency is arbitrary real and the axis is the CL Cartan covector \(e^{13}\). This excludes full nonlinear solutions assembled by finite-variation mixtures of compact translates of these particular amplitudes, even with mixed low, intermediate, and high frequencies.

Write a complex amplitude as \(b=u+iv\) in the BR coframe. The parameter \(z\) in the polynomial pencil below represents the differential \(z\,e^{13}\), with \(z=i\omega\) for the real frequency \(\omega\) multiplying the CL character. The theorem requires the trivial character of the sixteen even coordinate signs in the real defining five-dimensional representation.

## Completeness of the polynomial frame

Let \(K_0\) be the sign-invariant part of the BO pointwise H6 curvature kernel. The weight characters of the 15 known irreducible types, with their multiplicities, give \(\dim K_0=389\). The verifier reconstructs this count via the \(B_2\) Kostant weight multiplicity formula. Average the characters over the identity, the ten sign elements with two negatives, and the five with four negatives. The tangent and exterior slots use the same compact group action.

Let \(A\) denote covariant AN closure and \(D\) exterior multiplication by \(e^{13}\), with both in the complete rational physical frame. The original raw pencil also stacks the two hooks. Every hooked amplitude in the twisted kernel has coefficient curvature zero. Indeed the character is closed, and the square of \(d^\nabla+z\,d\chi\wedge\) is still the coefficient-curvature wedge. The converse after restriction to \(K_0\) follows from the definition. Thus the exact pencil on \(K_0^\mathbb C\) captures all of the declared sign-trivial amplitudes.

The data give 43 polynomial columns \(p(z)\) in the physical covariant frame. Their degrees occur with multiplicities \(2,23,13,5\) in degrees \(0,1,2,3\), respectively. Every coefficient of the hooks and \((A+zD)p(z)\), including the leading endpoint, vanishes exactly in the verifier. For a degree-\(d\) chain, \(A p_0=0\), \(A p_j=-D p_{j-1}\) for \(1\le j\le d\), and \(D p_d=0\). All coefficients lie in \(K_0\) by the preceding curvature identity applied to the polynomial, or by coefficient comparison after constructing the twisted fields.

Write \(R\subset K_0\) for the span of the 107 coefficient vectors and \(T\) for the span of the 64 vectors \(D p_j\) with \(j<d\) across all chains. Both dimensions follow by exact modular minors. In these bases the restriction \(A+zD:R\to T\) is the direct sum of the maps with equations \(z x_j-x_{j+1}\). Each such chain has rank \(d\) for every \(z\) and a one-dimensional kernel represented by \(p(z)\). The map \(K_0/R\to (\Lambda^7V^*\otimes V^*)/T\) induced by \(A+zD\) is injective for every complex \(z\). To see this, the verifier proves
\[
 \dim\bigl(T+A(K_0)+D(K_0)\bigr)\ge628=64+2(389-107).
\]
The two quotient images have dimension at most 282 each. Thus this bound shows both quotient maps injective and their images disjoint. It gives injectivity of their pencil for every \(z\), including no finite exceptional parameter.

The rank witnesses are modular reductions of rational data. For the curvature kernel itself, rational and modular dimension both equal 389. An invertible pivot minor modulo 65537 gives a rational kernel basis with denominators prime to 65537 by retaining the same free coordinates. Its reduction spans the entire modular kernel. Thus the rank lower bound using that modular kernel is also a characteristic-zero bound on the true quotients. Consequently the 43 polynomial columns span exactly the whole kernel at every finite \(z\).

## Exact Hermitian positivity

Let \(F_y\) have the same physical homogeneous receiver normalization as BV. The data give a new dyadic receiver \(y\), evaluated on the Clifford quadratic plus \((26/21)\|b\|^2\gamma\), in physical unit coordinates. For complex \(b=u+iv\), its Hermitian form is \(H_y(b,b)=F_y(u)+F_y(v)\). The associated real quadratic operator has norm at most \(16\sum|y_{a,m}|<60000\). One way to see the bound is to use normalized spinor trace in the six-product expansion of a physical Clifford coefficient. Each product of two input Clifford symbols and unit blades has absolute normalized trace at most the product of the two normalized Hilbert–Schmidt norms. Summing over distinct \(b,c\) external slots gives \(\sum_{b<c}\|b_b\|\|b_c\|\le(13/2)\|b\|^2\). With the six products and their \(1/3\) factor this yields at most \(13\|b\|^2\) per coefficient. The homogeneous scalar correction costs less than \(2\|b\|^2\). The declared factor 16 is therefore conservative.

To reproduce raw Gram matrices exactly, approximate each rational covariant chain coefficient in its physical unit frame by a 57-bit dyadic. The verifier bounds every inverse square root of the physical squared metric weight between adjacent 90-bit dyadics using integer square roots. From these bounds it obtains a rigorous sum of squared proxy errors across all vectors and degrees. The Frobenius norm of each of the four proxy coefficient matrices is less than seven. With the 60000 receiver bound, these estimates control the operator difference between the true and proxy receiver Grams.

The graded Pauli factorization of BV evaluates the dyadic basis polarizations modulo eleven split Gaussian primes. The 32-dimensional compressed spinor calculation uses only signed phases \(1,-1,i,-i\) and modular arithmetic. Exact Chinese remaindering reconstructs every Gram entry. For a dyadic numerator vector with largest coefficient \(\ell^1\) norm \(L\) and receiver numerator \(\ell^1\) norm \(Y\), the magnitude of an integer contracted core entry is at most \(39L^2Y\). This is the preceding trace estimate times three, followed by the operator norm bound for the real polarized form. The verifier uses the larger bound \(120L^2Y\). The prime product exceeds twice this bound, so integer recovery is unique. The homogeneous scalar correction and all complex powers are added exactly afterwards.

Evaluate each polynomial at \(z=i\omega\) to obtain the Hermitian Gram in chain coordinates. Order each chain basis column by its own degree. For \(0\le\omega\le21\), partition into 42 half-unit intervals and on each rescale an entire degree-\(d\) column by the positive constant \((1+\mathrm{upper\ endpoint})^{-d}\). This preserves Gram definiteness and gives a uniform bound for the effect of proxy errors. For the tail reverse each individual column in its degree and use the reciprocal \(t=1/\omega\in[0,1/21]\). Multiplication of a column by \((i\omega)^d\), and changing \(t\) to \(-t\), leave the definite signature unchanged. The chart at \(t=0\) includes the leading coefficient of every chain, with no common-degree omission. Negative frequencies follow by complex conjugation.

In these charts the Gram is a degree-at-most-six Hermitian matrix polynomial. The verified symmetric realifications of all 301 Bernstein coefficients are strictly positive after the exact rational correction for the proxy error. The verifier may use floating arithmetic to propose a dyadic Cholesky factor for a Bernstein matrix after subtracting a positive shift. It checks that factor by exact rational multiplication, and the exact maximum absolute row residual is strictly below the shift minus the rigorous proxy error. The factor itself enters only as a positive semidefinite rational matrix. Thus signs do not rely on floating eigenvalues. By positivity and the Bernstein partition of unity the true smallest coordinate Gram value exceeds \(10^{-6}\) throughout the two charts. Chart frame maps have operator norm less than 29 because each coefficient proxy has Frobenius norm less than seven and the interval weights on each power are at most one. Hence the common Hermitian receiver exceeds \(10^{-9}\|b\|^2\) at every nonzero sign-trivial axial amplitude.

## Field consequence and limits

Take an integrable finite-variation measure of axial imaginary frequencies and compact rotations, with amplitudes in the stipulated sign-trivial kernel. On the generic BW ray the QR compact factors converge to a common rotation followed by even diagonal signs. Those signs act as the identity on these amplitudes. The CE asymptotic Fourier measure therefore has every frequency atom in the same sign-trivial twisted kernel. Its Cesàro squared norm is the sum of the squared norms of these atoms, and the quadratic Cesàro identity holds for this one common receiver on all frequencies. The full nonlinear equation has zero receiver and positive constant squared norm \(R^2\), contradicting the uniform bound \(10^{-9}R^2\).

The averaging excludes global full solutions even when the measure has unbounded frequency support, subject to the finite variation with integrable amplitude norms. A solution on an open set is excluded when the represented field is real analytic, as for finite sums and bounded frequency support. The proof does not assert analytic continuation for arbitrary unbounded-frequency smooth superpositions. It also does not cover amplitudes in other sign characters or nonaxial frequencies.

## Arbitrary bounded single axial profiles

The same receiver excludes bounded differentiable sign-trivial profiles \(B(sK)=s\cdot b(\chi(s))\) on all of \(Y\), without any spectral representation. Use fixed orthogonal complements \(K_0=R\oplus R^\perp\) and \(W_{\mathrm{out}}=T\oplus T^\perp\). The exact quotient-rank argument shows that both quotient maps \(\bar A,\bar D:R^\perp\to T^\perp\) have disjoint images. In particular the linear map \((x,y)\mapsto\bar A x+\bar D y\) is injective and has a positive smallest singular value in any fixed physical metrics. It follows that
\[
 \|(\bar A+i\omega\bar D)v\|\ge c\sqrt{1+\omega^2}\|v\|.
\]
The inverse bound is uniform on the real frequency line. For the component in \(R\), each degree-\(d\) chain incidence matrix has entries \(i\omega x_j-x_{j+1}\). On the complement of its kernel, each squared singular value is at least \(1+\omega^2-2|\omega|\cos(\pi/(d+1))\), from the eigenvalues of its tridiagonal row Gram. Its minimum over all \(\omega\) is positive for \(1\le d\le3\). The degree-zero chains are already in the kernel. All coordinate maps from the chain basis and its image basis are fixed finite isomorphisms, so the lower bound remains uniform in the physical metric.

For \(x=r+v\) in the orthogonal decomposition, the quotient estimate bounds \(\|v\|\) by \(C(1+\omega^2)^{-1/2}\|(A+i\omega D)x\|\). Subtracting \((A+i\omega D)v\) from the residual costs at most a further constant factor, because the pencil has operator norm at most \(C'(1+|\omega|)\). The chain incidence estimate then bounds the distance of \(r\) to its kernel. This proves
\[
 \operatorname{dist}\bigl(x,\ker(A+i\omega D)\bigr)\le C_0\|(A+i\omega D)x\|
\]
uniformly for every real \(\omega\), with a finite constant \(C_0\) independent of \(x\).

Now let any bounded differentiable \(b:\mathbb R\to K_0\) solve the full equations globally in this axial parametrization. Its constant positive norm is \(R_0\). Choose a nonzero Schwartz cutoff \(\phi_L(t)=\phi(t/L)\) and write \(b_L=\phi_Lb\). Its Fourier closure equation has product-rule residual \(\widehat e_L\) with
\[
 \|e_L\|_2\le \|D\|_{\rm op}R_0\|\phi_L'\|_2=o(R_0\|\phi_L\|_2).
\]
The equation holds distributionally, and both sides used after multiplication lie in \(L^2\). The uniform inverse estimate gives a measurable pointwise orthogonal projection \(v_L(\omega)\) of \(\widehat b_L(\omega)\) into the twisted kernel with
\[
 \|\widehat b_L-v_L\|_2=o(\|b_L\|_2).
\]
For example, the explicit polynomial frame gives a continuous projection at every finite frequency. By Plancherel, boundedness of the receiver operator shows that replacing \(\widehat b_L\) by \(v_L\) changes its integrated receiver by \(o(\|b_L\|_2^2)\). Pointwise \(F_y(b)=0\) implies the integrated receiver of \(b_L\) is zero, while the uniform strictly positive lower bound on \(v_L\) gives a positive constant times \(\|b_L\|_2^2\) for large \(L\). This contradiction requires no analyticity or boundary behavior of \(b\).

## Local linear parameterization and a finite interval bound

The exact quotient splitting also parameterizes every local sign-trivial axial solution to the linear equation. If \(b(t)=r(t)+v(t)\) in \(R\oplus R^\perp\), the quotient equation \(\bar A v+\bar D v'=0\) and the independence of the two ranges imply \(v=0\). This works distributionally on any interval. In each chain, the coordinate \(x_{\ell,j+1}\) is the distributional derivative of \(x_{\ell,j}\). Thus
\[
 b(t)=\sum_{\ell=1}^{43}\sum_{j=0}^{d_\ell}
          f_\ell^{(j)}(t)\,p_{\ell,j}
\]
with independent free functions \(f_\ell\), interpreted as distributions when needed. The converse is coefficientwise for distributions and ordinary derivatives whenever the expressions are differentiable. In particular every local closed profile has all point values in the fixed 107-dimensional \(R\), irrespective of its frequency representation.

The uniform projection estimate also bounds the length of any local interval on which such a profile could satisfy the full nonlinear equations. Fix \(m=10^{-9}\), an operator bound \(M=60000\) for the receiver, and a finite uniform pencil inverse constant \(C_0\) from the preceding proof. On any interval choose a smooth nonzero scalar cutoff \(\phi\) of compact support, and extend \(\phi b\) by zero to the whole real line. Write \(u=\widehat{\phi b}\), let \(v\) be its pointwise orthogonal projection to the twisted kernel, and put
\[
\varepsilon=C_0\|D\|_{\rm op}\frac{\|\phi'\|_2}{\|\phi\|_2}.
\]
The product rule and positive constant norm give \(\|u-v\|_2\le\varepsilon\|u\|_2\). The integrated receiver of \(u\) is zero. By orthogonal projection \(\|v\|_2\le\|u\|_2\), and the integrated receiver difference is at most \(2M\varepsilon\|u\|_2^2\). Strict positivity on \(v\) gives \(m(1-\varepsilon)^2\|u\|_2^2\) as a lower bound if \(\varepsilon<1\). Hence \(\varepsilon\ge m/(8M)\): otherwise \(\varepsilon<1/2\), the lower bound exceeds \(m\|u\|_2^2/4\) and the difference is smaller than it. Rescale one fixed nonzero test function of compact support on \((0,1)\) into an interval of length \(L\). It has derivative ratio \(L^{-1}\|\phi_1'\|_2/\|\phi_1\|_2\). Thus every interval of a full local profile has length bounded by
\[
 L\le \frac{8M C_0\|D\|_{\rm op}}{m}\,\frac{\|\phi_1'\|_2}{\|\phi_1\|_2}.
\]
The constant is finite and independent of the profile. [CZ](CZ-h6-axial-local.md) subsequently gives a different receiver positive on the entire local value space, which excludes arbitrarily short intervals as well.

## Reproduction

Run from the root of the repository. Python, NumPy, and SymPy suffice.

```text
python -B computation/geometric-unity/CV_h6_axial.py
python -B computation/geometric-unity/CV_h6_axial.py --full
```

The default checks completeness, all geometry, proxy coordinate error, and all rational positivity matrix checks against the archived Gram. The full option also reconstructs every archived Gaussian-integer Gram entry via the eleven primes. Use the full option for independent verification of the whole certificate.

[Review106](../../reports/reviews/REVIEW106-GU-CR-CS-CV-2026-11-14.md) records the independent checks. External peer review remains outstanding.
