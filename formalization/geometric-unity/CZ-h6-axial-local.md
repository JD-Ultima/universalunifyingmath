# CZ. A local H6 obstruction for sign-trivial axial triangular profiles

On the positive \(H_6\) observer, no nonzero full nonlinear solution of the form
\[
 B(sK)=s\cdot b(\chi(s))
\]
exists on any open axial interval when \(b(t)\) is fixed by all sixteen determinant-one diagonal signs and \(\chi\) is the Cartan additive character of CL. The proof also excludes global finite-variation mixtures of compact rotations of bounded such linearly closed profiles. It requires no restriction on frequencies and no analyticity for the global conclusion.

## Complete local value space

Use the notation of CV. At each fixed \(t\) in an axial interval, \(b(t)\) is a real \(H_6\) hook satisfying the coefficient-curvature wedge, and lies in the exact sign-trivial curvature kernel \(K_0\) of dimension 389. Covariant closure has the constant-coefficient axial equation \(A b+D b'=0\), where \(A\) is the complete BR exterior covariant derivative and \(D=e^{13}\wedge\) on the form slots. Curve derivatives may be interpreted distributionally.

Let \(R\) and \(T\) be the exact coefficient and image spans of the CV polynomial chains, of dimensions 107 and 64. CV proves that the two maps from \(K_0/R\) to the output quotient by \(T\) induced by \(A\) and \(D\) are injective with disjoint images. Applying the quotient to closure separately forces \(\bar A[b]=0\) and \(\bar D[b]'=0\), hence \(b(t)\in R\). On \(R\) the direct-sum chain relations are \(x'_{\ell,j}=x_{\ell,j+1}\), so the complete local linear profile has the 43 free functions and their derivatives described in CV. Conversely any point value of \(R\) can be attained by prescribing finitely many scalar Taylor coefficients independently at one point. Therefore \(R\) is exactly the full local axial value space, regardless of frequency representation.

The verifier checks CV's rank and closure data, and independently checks a second rational family of 43 polynomial chains whose 107 coefficient vectors are independent. Those coefficient vectors are in \(R\). At every complex parameter the value of each second chain is in the exact CV kernel, hence is in \(R\), and polynomial coefficient extraction by finite interpolation gives the same for all its coefficients. Their count and independence make them another basis of \(R\). The certificate expresses each of the 107 rational test vectors in this second coefficient basis. The verifier reconstructs them and checks the hash of the canonical ordered rational coefficients without loading a duplicate basis archive.

## One physical pointwise receiver

Let \(F_y\) be the linear functional of the four grade outputs of the homogeneous full quadratic \(q+(26/21)\|B\|^2\gamma\), with the same physical sign, Clifford volume, and unit coframe as BV. For a full H6 solution, tracing the grade-one equation fixes the norm and gives \(F_y(B)=0\) for every \(y\). The certificate uses a particular dyadic vector \(y\) supported in sign-trivial receivers. It checks strict positivity on the whole real linear space \(R\), a stronger condition than positivity separately on the frequency kernels. The \(R\)-receiver is also positive on the Hermitian extension of \(R\), by splitting real and imaginary parts.

Verify the Gram as follows. A rational covariant tensor coefficient at \((a,m)\) converts to the unit coframe by dividing by \(\sqrt{g_a\prod_{i\in m}g_i}\). Integer square-root bounds give rigorous intervals for these reciprocal square roots. The verifier rounds to 48-bit dyadics, and bounds the Frobenius norm of the total error of the 107-column matrix by \(10^{-11}\). The Frobenius norm of its proxy matrix is less than 11. The exact rational upper bound \(16\sum |y_{a,m}|<30000\) for the quadratic operator is established by the normalized spinor trace estimate in CV. Therefore \(30000\bigl(22\cdot10^{-11}+10^{-22}\bigr)\) bounds the possible effect of the proxy perturbation on the Gram, with the exact smaller operator bound used in the verifier.

The dyadic proxy Gram is reconstructed using BV's 32-dimensional graded Pauli representation. Each signed Gaussian product is evaluated over nine split primes, and the exact coefficient bound \(120Y L^2\) makes signed CRT reconstruction unique. Here \(Y\) is the receiver numerator \(\ell^1\) norm and \(L\) is the largest input proxy numerator \(\ell^1\) norm. The \(26/21\) scalar adjustment is added exactly. This is the same compressed polarization normalization as CV.

A stored exact rational factor \(P\) verifies positivity. The maximum absolute row sum of the difference between the proxy Gram, \(2I/5\), and \(PP^T\) is below \(3.64\cdot10^{-10}\). The exact metric perturbation is below \(8.37\cdot10^{-7}\). Thus the actual rational-source Gram is bounded below by \(0.39I\). Its column map has squared operator norm below 121, so the receiver has the uniform physical bound
\[
 F_y(b)>\tfrac3{1000}\|b\|^2
 \qquad(0\ne b\in R).
\]
This shows the 107 verified rational source vectors are independent and completes the pointwise obstruction. The Clifford Gram and matrix positivity use no numerical eigenvalue acceptance. The optimizer served only to propose a receiver, rational source combinations, and the rational factor.

## Local and rotated consequences

For a nonzero linearly closed axial profile the sign-trivial field value is in \(R\) at every point. The strict receiver immediately rules out full nonlinear equality anywhere its norm is positive. This is a local exclusion on arbitrarily short intervals, with no extension assumption.

For the global conclusion, take a real finite-variation integral of compact translates of these axial profiles. The integrand may carry an arbitrary bounded globally defined axial profile, and its pointwise values in the AN frame belong to \(R\) by the linear argument. Assume the bound of its supremum norm is integrable against the variation measure. On the generic BW axis ray, the compact QR factor of each translated integrand tends almost everywhere to one common compact rotation followed by a diagonal sign. The sign acts as the identity on \(R\). Replacing the compact QR factor with its limit thus gives an \(R\)-valued integral in the common frame. The norm difference to the true value is bounded by the integral of the QR tensor-action difference times the profile supremum norm, so it tends to zero by dominated convergence. The profile itself can continue oscillating without a limit in that frame. The strict receiver bound on \(R\) contradicts the positive constant norm and vanishing receiver of any global full solution. For a globally analytic representation the same conclusion follows from full equality on a nonempty open set by analytic continuation. A merely smooth representation equal to the full equations only locally does not have that continuation property.

The scope is the sign-trivial axial profile and its compact rotations and bounded mixtures as stated. The certificate concerns the separate pure \(H_6\) branch. It leaves other sign characters, nonaxial triangular profiles, unrestricted local closed fields, and general mixtures of distinct pure hooks open.

## Reproduction

Run from the root of the repository. Python, NumPy, and SymPy suffice.

```text
python -B computation/geometric-unity/CZ_h6_axial_local.py
python -B computation/geometric-unity/CZ_h6_axial_local.py --full
```

The default recomputes CV geometry and the quotient ranks, verifies all additional rational span memberships, checks the proxy Gram modulo one prime, and verifies the exact positive shift in the archived rational matrix. Full mode checks all nine primes and thereby rebuilds the complete proxy Gram.

[Review106](../../reports/reviews/REVIEW106-GU-CR-CS-CV-2026-11-14.md) records the independent checks. External peer review remains outstanding.
