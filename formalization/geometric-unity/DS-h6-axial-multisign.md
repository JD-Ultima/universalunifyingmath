# DS. One local H6 axial obstruction on three sign characters

No nonzero pure \(H_6\) nonlinear solution on a positive observer open set can have the axial triangular form
\[
B(sK)=s\cdot b(\chi(s))
\]
with values in the sum of the diagonal-sign characters \(0,15,27\). Here \(s\in AN\), \(\chi\) is the additive axial Cartan coordinate of CL, and a character \(c\) evaluates on a determinant-one diagonal sign mask \(g\) by \((-1)^{|c\cap g|}\). The claim applies on an arbitrary open interval, with just \(C^1\) regularity for the profile. The same receiver excludes a global full solution represented by a finite-variation mixture of compact rotations of bounded such profiles.

## Exact local axial value spaces

The two hooks and coefficient-curvature wedge imply that each profile component takes values in its sign-character curvature kernel \(K_c\). The \(SO(5)\) representation and multiplicities of BO give \(\dim(K_0,K_{15},K_{27})=(389,443,443)\). The verifier recomputes the sign projector dimensions by B2 weight enumeration and character averaging. The triangular connection commutes with determinant-one diagonal signs. Therefore each component separately satisfies the constant-coefficient closure equation
\[
A b+D b'=0,\qquad D=e^{13}\wedge
\]
with \(A\) given by the BR exterior covariant derivative and \(D\) acting on the form slots. In the computation the two hook equations are included in this stack, with identically zero rows for their derivative term. Both maps are multiplied by the same exact denominator.

For a rational polynomial chain \(\sum_{j=0}^d r_j z^j\) in the kernel of \(A+zD\), the coefficient directions obey
\[
A r_0=0,\quad A r_j=-D r_{j-1}\quad (1\le j\le d),
\qquad D r_d=0.
\]
Every coefficient lies in \(K_c\). For real \(z\), its exponential closed field proves curvature compatibility of its evaluation. Interpolation proves the same for each coefficient. Write \(R\) for the span of all chain coefficients and \(T=A(R)+D(R)\). The certificates give the following dimensions and modular lower bounds.

| Sign mask | \(\dim K_c\) | Chains of degrees 0, 1, 2, 3 | \(\dim R\) | \(\dim T\) | \(\operatorname{rank}[\bar A,\bar D]\) lower bound | Additional modes |
|---|---:|---|---:|---:|---:|---|
| 0 | 389 | 2, 23, 13, 5 | 107 | 64 | 564 | |
| 15 | 443 | 4, 35, 26, 6 | 176 | 105 | 523 | 2 at \(z=0\), 7 at \(z=-4\), 2 at \(z=-32/5\) |
| 27 | 443 | 4, 31, 29, 7 | 181 | 110 | 524 | |

Here \(\bar A,\bar D\) are maps from \(K_c/R\) to the output quotient by \(T\). The ranks lift from \(\mathbf F_{65537}\) to characteristic zero. The curvature rows have denominators that are units at this prime, and their modular rank agrees with the characteristic-zero rank. A maximal pivot minor which is a unit at the prime therefore constructs a localized kernel basis whose reduction spans the full modular kernel. The chain and image sizes are exact by independent modular columns and the displayed relations.

To see why these ranks control all local profiles, take closure modulo \(T\). For \(v(t)=[b(t)]\) one has \(\bar A v+\bar D v'=0\). Thus \(v(t)\) is in the kernel of \(\bar A\) followed by reduction modulo \(\operatorname{im}\bar D\). Its dimension is at most
\[
\dim(K_c/R)-\bigl(\operatorname{rank}[\bar A,\bar D]-\operatorname{rank}\bar D\bigr).
\]
The bound is zero for \(c=0,27\) and eleven for \(c=15\). The eleven rational vectors in the last block are verified in the full original equations \((A+z D)b=0\) at their stated rates, with independent classes modulo \(R\). Consequently all profile values are in
\[
R_0^+\oplus R_{15}^+\oplus R_{27}^+,\qquad
(\dim R_0^+,\dim R_{15}^+,\dim R_{27}^+)=(107,187,181),
\]
where \(R_c^+\) adds the tabulated vectors to \(R\).

Conversely each full chain has the profile \(\sum_j f^{(j)}(t)r_j\) for an arbitrary sufficiently smooth \(f\), with arbitrary values for those derivatives at a prescribed point. The eleven extra modes use \(e^{zt}b\). Thus these are exactly the local point-value spaces. In character 15 the additional modular check proves \(\bar D\) injective. The quotient equation then gives the eleven scalar exponential ODEs in the displayed eigenbasis, with no additional polynomial factor. This local result is independent of a Fourier representation or restrictions on global boundedness.

## A common exact receiver

Use the physical homogeneous quadratic
\[
q^0(B)=q(B)+\tfrac{26}{21}\|B\|^2\gamma
\]
and the Clifford calibration of BV and CV. A full \(H_6\) solution has \(q^0(B)=0\) and positive constant norm. The archived dyadic functional \(y\) is supported on sign-trivial receivers in grades 1, 5, 9 and 13. It is a *single* functional across the three blocks. It gives \(F_y(B)=\langle y,q^0(B)\rangle\). The polarization of \(F_y\) between two different characters vanishes, since some diagonal sign multiplies the two inputs by opposite signs and fixes the receiver. Thus \(F_y\) adds across the blocks, as does squared norm.

The verifier reconstructs rational test columns as exact combinations of the complete chain and regular-mode columns. BC covariant coefficients at \((a,m)\) are converted into the unit coframe by dividing by \(\sqrt{g_a\prod_{i\in m}g_i}\). Exact enclosure of each reciprocal square root yields a 46-bit dyadic proxy and a total Frobenius error below \(10^{-10}\) for each block. The proxy squared Frobenius norm is below \(n+1\) for an \(n\)-column block. The exact operator estimate \(16\|y\|_1\) from CV bounds the change in the Gram. A stored rational factor and the maximum absolute row sum of its defect certify a lower bound at each block relative to the squared norm of its physical column map. The estimates are, in mask order,
\[
0.0014666792,\quad 0.0008453454,\quad 0.0008731017 .
\]
Every exact estimate exceeds \(1/1250\). In particular the rational test columns form bases, and the combination satisfies
\[
F_y(b)>\tfrac1{1250}\|b\|^2\qquad(0\ne b\in R_0^+\oplus R_{15}^+\oplus R_{27}^+).
\]
The linear block geometry and positivity inputs are rational. In exhaustive mode, the proxy Gram is rebuilt using BV's exact finite-field Clifford representation across nine split primes per block. With \(Y\) the receiver numerator \(\ell^1\) norm and \(L\) the largest proxy column numerator \(\ell^1\) norm, the signed coefficient bound \(120YL^2\) gives unique reconstruction. The scalar adjustment is added exactly.

## Consequences and scope

For a local profile taking values in the three allowed sign characters, any nonzero point value violates the full equations. For a global finite-variation integral of compact translates, assume the supremum norms of the bounded axial profiles are integrable against the variation measure. On the BW generic QR axis ray, the compact factors converge almost everywhere to a common rotation followed by an element of the diagonal sign group. Each limiting sign preserves the sum \(R_0^+\oplus R_{15}^+\oplus R_{27}^+\). Replacing the compact factors by their limits therefore gives a value in one common rotated value space, even if the profiles themselves continue to oscillate. The difference from the true field tends to zero by dominated convergence, and the common receiver lower bound then contradicts the full equations and positive constant norm. A globally analytic representation extends this consequence to full equality on any open set by analytic continuation.

The result is for one axial Cartan coordinate and the stated three characters. It does not decide the remaining sign characters, more general multi-Cartan profiles, or solutions outside this class.

## Reproduction

Run from the root of the repository. The verifier reuses CV's public
geometry and one CZ chain archive. Python, NumPy and SymPy suffice.

```text
python -B computation/geometric-unity/DS_h6_axial_multisign.py --character 15 --geometry-only
python -B computation/geometric-unity/DS_h6_axial_multisign.py --full
```

The selected regression checks the exact closure geometry and exceptional
modes of character 15. The default command checks every chain and mode,
all modular quotient ranks, each proxy Gram modulo one prime, and all
exact positive shifts. Exhaustive mode checks every Gram at all nine
primes, giving the full independent reconstruction.

[Review111](../../reports/reviews/REVIEW111-GU-DS-DH-DI-2026-11-14.md)
records the checks. External peer review remains outstanding.
