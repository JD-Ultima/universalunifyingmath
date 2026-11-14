# DL. A global mixed H6 field through a balanced characteristic

The real characteristic orbit of [DK](DK-h6-characteristic-orbit.md)
admits a nonzero real analytic closed H6 hook on the whole positive
observer. The new field takes values in the sum of types \((7,3)\)
and \((5,5)\). It permits arbitrary entire multipliers in six
globally analytic complex functions before taking real parts.

## DL1. An eight-dimensional leaf algebra

Use the [BO](BO-h6-full-curvature-kernel.md) split coframe
\(e^0,\ldots,e^{13}\) and tangent matrices \(p_0,\ldots,p_{13}\).
Set \(\xi=e^{12}-e^{13}\) and
\[
 \beta=e^0\wedge e^1\wedge e^3\wedge e^8\wedge e^9\wedge\xi,
 \qquad v_0=(e^1)^\sharp=(E_{12}+E_{30})/4.
\]
The two nonzero covariant entries of \(v_0\otimes\beta\) are
\((1,4875)=1\) and \((1,8971)=-1\).
The verifiers check both hooks, full coefficient curvature and the
exterior symbol \(\xi\wedge(v_0\otimes\beta)=0\). The physical direction
of \(\xi\) is proportional to \(\operatorname{diag}(1,1,-1,-1,0)\).
With BO highest tensors \(b_{55},b_{73}\) and the lowering operator
\(F_1\) for the first simple root, this value equals
\((F_1^2 b_{73}-2b_{55})/6\).

Take the complex tangent subspace and the compact centralizer
\[
 D=\operatorname{span}(p_2,p_4,p_5,p_6,p_7,p_{10},p_{11},p_{12}+p_{13}),
 \qquad \mathfrak k_{v_0}=\{A\in\mathfrak k_{\mathbb C}:[A,v_0]=0\}.
\]
Exact brackets give \([D,D]=\mathfrak k_{v_0}\) of dimension four
and \([\mathfrak k_{v_0},D]\subseteq D\). This makes
\(\mathfrak h=D\oplus\mathfrak k_{v_0}\) a 12-dimensional complex
Lie subalgebra. Compact conjugation does not preserve it.
The complementary tangent basis is
\((p_0,p_1,p_3,p_8,p_9,(p_{12}-p_{13})/2)\). Exponentiating this
complement and then \(D\) gives a holomorphic chart at the origin.

The same local leaf argument applies to any representative \(g\)
whose Maurer–Cartan form along a transverse level set takes values
in \(\mathfrak h\). Its compact part centralizes \(v_0\).
Thus \([g,v_0]\) is covariantly parallel along the leaves.
Wedging with all six transverse differentials gives a closed
coefficient-valued six-form. In the moving frame the differentials
annihilate \(D\), so this form is a scalar multiple of
\(v_0\otimes\beta\). The hooks and the two-type membership are
therefore preserved.

## DL2. Global formula

Use [CU](CU-h6-global-analytic-closure.md)'s real analytic spectral
graph data in the reordered split axes \((1,3,4,0,2)\).
In that basis write \(y=L_0^{-1}YL_0^{-T}\) and
\[
 y=\begin{pmatrix}A&b&C\\b^T&d&e\\C^T&e^T&\bar A\end{pmatrix},
 \quad
 n=\begin{pmatrix}I&0&0\\v^T&1&0\\S&0&I\end{pmatrix},
 \quad
 n^{-1} y n^{-T}=
 \begin{pmatrix}A'&b'&R\\b'^T&d'&0\\R^T&0&0\end{pmatrix}.
\]
Here \(S\) is the contractive symmetric graph, \(R=C-AS\),
and \(v\) is the vector determined by the CU block equations.
The positive eigenvalues of \(R\) give \(\operatorname{tr}R>0\).
Let \(N_W\) be the first two columns of \(n\) and put
\[
 (t_1,\ldots,t_5)=(S_{11},S_{12},v_1/2,S_{22},v_2/2),\qquad
 w=R_{11}/\operatorname{tr}R,
 \qquad
 H=\begin{pmatrix}R_{12}w&R_{11}(1-w)\\
 R_{11}(1-w)&R_{21}(1-w)\end{pmatrix}.
\]
The globally real analytic complex tangent coefficient in split
covariance coordinates is
\[
 B=(N_W H N_W^T)\,dt_1\wedge\cdots\wedge dt_5\wedge dw .
\]
Physical covariance components use
\(L_0N_W H N_W^T L_0^T\). The tangent metric turns these into
the covariant coefficient used in the hook equations.

To prove closure, work first near the identity. Set
\(z=R_{11}/R_{22}\), \(u=\tfrac12\log z\) and
\(J=\operatorname{diag}(z^{1/4},z^{-1/4})\) on each outer block.
Then \(\widetilde R=J^{-1}RJ^{-1}\) has equal diagonal entries.
Put \(Q=\left(\begin{smallmatrix}0&1\\1&0\end{smallmatrix}\right)\).
The matrix \(\widetilde RQ\) is symmetric, so choose \(L\) near
the identity with \(LQL^T=2\widetilde RQ\). Such a holomorphic
choice exists by the local submersion of the congruence map.
Set the last diagonal block to \(QLQ\), the middle scalar to
\((\det L)^{-2}\), and use the two upper blocks to fit \(A',b'\)
as in CU. Its determinant-one representative lies in the
upper block subgroup with last block \(QLQ\) and zero middle-right
block. That subgroup has Lie algebra \(\mathfrak h\), as verified
in the split basis. The representative at a point has the form
\(g=n\operatorname{diag}(J,1,J)h\). On a level set of
\((t_1,\ldots,t_5,u)\), both \(n\) and \(J\) are constant,
providing the Maurer–Cartan hypothesis of DL1.

In the split covariance model \(v_0\) is \(Q/8\) in the first
two axes, with all other entries zero. The \(h\) congruence
gives \(\widetilde RQ/4\). After the \(J\) congruence, the
resulting block in the \(n\) frame is
\[
 \frac14\begin{pmatrix}
 R_{12}\sqrt z&R_{11}/\sqrt z\\
 R_{11}/\sqrt z&R_{21}/\sqrt z
 \end{pmatrix}.
\]
Multiply the local transverse volume with last differential \(du\)
by \(8\exp(3u)/(1+\exp(2u))^3\). Using \(du=dz/(2z)\) and
\(w=z/(1+z)\), the square roots cancel and the formula is exactly
the globally defined \(B\) above. Multiplication by this transverse
function preserves closure.

The formula for \(B\) is globally analytic on the connected real
observer. The closure, hook and complementary-type residues are
analytic and vanish near the identity, so they vanish everywhere.
The local logarithm and factorization are only used on that
neighborhood. At the identity \(R=I/2\), \(H=Q/4\) and
\(dw=\xi/2\), so \(B(o)=v_0\otimes\beta\). Its compact real
part has squared norm \(1/131072\).

Every entire function of \((t_1,\ldots,t_5,w)\) is a global
multiplier because its differential wedges to zero with the
six transverse differentials. Taking real parts preserves every
linear equation. [BQ](BQ-h6-invariant-moments.md) has strictly
negative grade-nine coefficients on both constituent types.
Thus the family cannot solve the full nonlinear equations wherever
its value is nonzero.

    python -X utf8 computation/geometric-unity/DK_h6_characteristic_field.py
    python -X utf8 computation/geometric-unity/DK_h6_characteristic_global.py

The verifiers check the Lie brackets, exact transversality at the
origin, coefficient congruences, square-root cancellation and the
base normalization.
