# DK. One real characteristic orbit for the three mixed gradients

Let \(E_{55}\) be the real \((5,5)\) fiber and let \(L_\xi b\) be the direct sum of the three [DK](DK-h6-mixed-gradients.md) generalized-gradient symbols \(G_{55},G_{54},G_{53}\). The complete real characteristic set away from zero is the cone of conjugates of \(\xi_0=\operatorname{diag}(1,1,-1,-1,0)\). Its kernel has real dimension four at each of those covectors.

Use \(E_{55}^{\mathbb C}=\operatorname{Sym}^{10}\mathbb C^4\). Normalize the five \(4\times4\) spin matrices by \(\gamma_i\gamma_j+\gamma_j\gamma_i=2\delta_{ij}\). The compact generator \(\rho(J_{ij})\) is the degree-ten polynomial action of \(\gamma_i\gamma_j/2\). On \(V_{14}\otimes E_{55}^{\mathbb C}\) the cross Casimir is \(X=-2\sum_{i<j}\rho_{14}(J_{ij})\otimes\rho(J_{ij})\). The compact Casimir normalization is \(C_{\lambda}=l(l+3)+m(m+1)\). DK's six summands have cross eigenvalues
\[
\begin{array}{c|rrrrrr}
\lambda&(7,5)&(6,5)&(6,4)&(5,5)&(5,4)&(5,3)\\
X&20&4&-6&-10&-20&-28.
\end{array}
\]
Thus \(R=(X-20)(X-4)(X+6)\) vanishes precisely on the three DK-allowed modules. For every \(\xi\), \(R_\xi b=R(\xi\otimes b)\) has the same kernel as \(L_\xi\).

In the polynomial spin model the verifier constructs all four \(R_\xi\) for a rational trace-orthogonal Cartan basis. On \(\xi_0\) it checks directly that the four coordinate degree-ten powers all vanish. The other 282 columns are independent over \(\mathbb Q(i)\) by the rank test below. These give the complete complex kernel. For a real \(\xi_0\), the complex kernel is stable under the compact real structure, so its real form has dimension four. The coordinate degree-ten powers describe the complex span.

For local isolation let \(K=\ker R_{\xi_0}\), let \(E\) be the coordinate complement of the four powers, and let \(F\) be the image of \(E\) under \(R_{\xi_0}\). The first, second and fourth standard Cartan axes form a complement to \(\mathbb C\xi_0\) in this Cartan basis. The verifier joins the 282 columns of \(R_{\xi_0}|_E\) and all twelve columns of the three transverse \(R_\eta|_K\). The matrix has rank 294 modulo 65537 with \(i\mapsto256\). Any nonzero maximal minor modulo that prime lifts to a nonzero Gaussian-integer minor. A common integer scale has already cleared the rational metric denominators. In particular, the projected linear map
\[
(\mathrm{Cartan}/\mathbb C\xi_0)\otimes K\longrightarrow
(V_{14}\otimes E_{55}^{\mathbb C})/F
\]
is injective.

Fix any positive norms. By finite-dimensional injectivity this map has a positive least singular value. On \(E\), the image component in \(F\) remains invertible onto \(F\) for \(\xi_0+\eta\) sufficiently close. A vector \(b=k+e\) in the kernel then satisfies \(\|e\|\le C\|\eta\|\|k\|\) by the \(F\) equation. The quotient equation gives \(c\|\eta\|\|k\|\le C'\|\eta\|\|e\|\). It forces \(k=e=0\) for all sufficiently small nonzero \(\eta\). Real compact diagonalization relates the Cartan result to the orbit near \(\xi_0\). The neighborhoods also extend to the three-dimensional complex diagonal projective chart.

## Exact cover of the real Weyl chamber

Write \(P\) for the orthogonal projection onto the three forbidden eigenspaces of \(X\). In the monomial spin basis \(z^\alpha\), the positive Hermitian form is proportional to \(\operatorname{diag}(\prod_j\alpha_j!)\). In the Cartan basis used by the verifier the trace metric is \(\operatorname{diag}(2,6,12,20)\). These are the source and target metrics in the calculations below. Let \(F_\xi=P(\xi\otimes\cdot)\), with the same kernel as \(R_\xi\) and \(L_\xi\). Spectral interpolation with rational coefficients gives the exact \(P\) from the six tabulated eigenvalues. Its Cartan-to-Cartan block supplies all relevant Gram products because \(P^*P=P\). The verifier checks this block's exact metric symmetry and proves the restriction of \(P\) to Cartan \(\otimes E_{55}^{\mathbb C}\) has squared norm below \(5/8\). Consequently
\[
\|F_\xi-F_\eta\|\le\sqrt{5/8}\,\|\xi-\eta\|
\]
for any Cartan covectors.

For the neighborhood of the known ray take a trace-unit representative \(x=\operatorname{diag}(1,1,-1,-1,0)/2\). Set \(K=\ker F_x\) and let \(E=K^\perp\), \(W=F_x(E)\). The exact Gram proof gives
\[
s_{\min}(F_x|_E)>1/4,\qquad
s_{\min}\bigl((I-\Pi_W)P|_{x^\perp\otimes K}\bigr)>3/5,
\]
where \(x^\perp\) is taken in the four-dimensional Cartan. The second check is the Schur complement of the unshifted Gram of \(F_x|_E\) in the joint Gram matrix of \(F_x|_E\) and the three transverse maps on \(K\), with \(9/25\) subtracted in the transverse metric. The first check is the axial complement Gram shifted by \(1/16\). All the factors are rational. In the exact program the vector \((0,2,1,0)\) of standard Cartan coordinates is \(6x\). Its scale in the base Gram and cross terms is retained.

For \(\xi=x+h\) with \(\|h\|<1/25\), write \(h=a x+\eta\) with \(\eta\perp x\). If \(F_\xi(e+k)=0\) with \(e\in E\) and \(k\in K\), projection to \(W\) yields
\[
\|e\|\le
\frac{(4/5)\|\eta\|}{(1-1/25)/4-(4/5)/25}\,\|k\|.
\]
Here \(|a|\le1/25\) and \(\sqrt{5/8}<4/5\). Projection to \(W^\perp\) yields \((3/5)\|\eta\|\|k\|\le(4/5)\|\eta\|\|e\|\). The two inequalities contradict each other for \(k\ne0\) and \(\eta\ne0\), since the first ratio is at most \(2/13\). If \(k=0\), the \(W\) equation yields \(e=0\). Thus throughout this ball only points on the known ray are characteristic.

For the remaining chamber use sorted diagonal entries dotting to zero, with sum of their successive gaps set to one. Their ordered gap variables \(u_1,\dots,u_4\) are nonnegative and sum to one. This normalizes every nonzero ordered ray within the simplex. The vertex \(u_k=1\) has its first \(k\) diagonal entries \((5-k)/5\) and its other entries \(-k/5\). The characteristic point corresponds to \(p=(0,1/2,1/2,0)\), the midpoint of the second and third vertices. Inserting this midpoint triangulates the simplex into two tetrahedra. Each recorded subdivision chooses an edge \(AB\) and replaces one endpoint by its midpoint in each child, which partitions that tetrahedron. The archived list records 700 visited nodes, of which 349 are subdivisions and 351 are leaves. The program reconstructs all coordinates by that recursion.

For 27 leaves each vertex is within Euclidean trace distance \(1/25\) of the characteristic point, so the entire leaf lies in the ball by convexity. These use the local proof after a fixed permutation. The other 324 leaves use their exact rational barycenter \(z\). Let \(r\) be the largest trace distance from \(z\) to a vertex, which bounds its distance to any point of the leaf. The verifier proves the strict matrix inequality \(F_z^*F_z>(5/8)r^2 I\). The Lipschitz estimate then shows \(F_\xi\) is injective throughout that leaf. This covers every real nonzero trace-free symmetric covector up to positive scaling and compact diagonalization. Any permutation of a diagonal is realizable by an \(SO(5)\) conjugation because a compensating diagonal sign commutes with it.

## Arithmetic and reproduction

All positive-definiteness checks use a dyadic factor proposal and exact integer acceptance. After a positive diagonal congruence and division by the largest matrix entry in absolute value, let \(A\) denote the tested matrix. The verifier forms dyadic \(T\) and \(J\) by rounding the proposed numerical factor and its inverse. It proves \(\|I-TJ\|_2<1/2\) with the square root of the product of the exact maximum absolute row and column sums. Thus \(T\) is invertible, \(\|T^{-1}\|_2\le2\|J\|_2\), and \(TT^T\ge I/(4\|J\|_2^2)\). It checks that the maximum absolute row sum of \(A-TT^T\) is below that last lower bound, using \(\|J\|_2^2\le\|J\|_1\|J\|_\infty\). These tests consume integers and reject any approximate factor that fails the bound. Python arbitrary integers handle the large Gram and factor numerators. Fixed-width modular products have operands below \(65537\) with at most two products per additive update. The projector's sparse Gaussian-integer multiplication bounds every real and imaginary component partial sum below \(2^{53}\). For large intermediate columns it splits coefficients into low and high 24-bit pieces, checks each multiplication independently, and recombines the result with Python integers before the next step.

The short verifier checks the exact characteristic kernel and transverse rank. The full verifier also checks the global Weyl chamber cover.

```text
python -X utf8 computation/geometric-unity/DK_h6_characteristic.py
python -X utf8 computation/geometric-unity/DK_h6_characteristic_cover.py
```

The second command rebuilds the rational projector, verifies the operator and local gaps, and replays the entire characteristic cover. The theorem concerns the projected \((5,5)\) principal symbol forced by DK for mixed closed fields.
