# CU and CT. Complete global analytic lifts for H6 closure

The space of values of smooth closed H6 hooks on a neighborhood
of any point of the flat positive observer is exactly the
6,820-dimensional submodule of the curvature kernel with type
\((2,0)\) removed. Every value in this submodule has a globally
defined real analytic linear lift. The new \((5,5)\) fields
represent nonzero local classes in the curvature-kernel complex.

## CU1. A leaf with a parallel coefficient

Work on the complexification \(G_{\mathbb C}/K_{\mathbb C}\) of
\(SL(5,\mathbb R)/SO(5)\). Use the split coordinates of [BO](BO-h6-full-curvature-kernel.md),
with metric \(G_{01}=G_{10}=G_{23}=G_{32}=2\), \(G_{44}=1\) and
all remaining entries zero. Write \(\mathfrak{sl}_5=\mathfrak k\oplus
\mathfrak p\) for the \(G\)-skew and \(G\)-selfadjoint trace-free
parts. A tangent \(P\in\mathfrak p\) is normalized by
\(t\mapsto\exp(tP/2)K_{\mathbb C}\). In the homogeneous bundle
with coefficient \(\mathfrak k\) or \(\mathfrak p\), put
$$
 \theta=\frac12\sum_i e^i\wedge \operatorname{ad}_{P_i}.
$$
This alternates the bracket with the soldering form and
interchanges the coefficient bundles. The soldering form is
covariantly closed by torsion-freeness. Equivariance of the
bracket and the exterior product rule give
\(d^\nabla\theta=-\theta d^\nabla\).

Let \(W=\operatorname{span}(e_1,e_3)\) in the defining
five-dimensional representation, a maximal isotropic plane.
Set
$$
 D=\{P\in\mathfrak p:PW\subset W\},\qquad
 \mathfrak k_c=\{K\in\mathfrak k:KW\subset W,\
                   \operatorname{tr}(K|_W)=0\},\qquad
 c=(E_{12}-E_{30})/4.
$$
Here \(E_{ij}\) is the defining matrix unit. Under the equivariant
identification \(\Lambda^2\mathbb C^5\simeq\mathfrak k\), \(c\)
is a nonzero multiple of the volume of \(W\). Every element of
\(\mathfrak k_c\) acts trivially on it.

For \(P,Q\in D\), their commutator lies in \(\mathfrak k\),
preserves \(W\), and has trace zero there since its restriction
is itself a commutator. Similarly
\([\mathfrak k_c,D]\subset D\) and
\([\mathfrak k_c,\mathfrak k_c]\subset\mathfrak k_c\).
Thus \(\mathfrak h=D\oplus\mathfrak k_c\) is a Lie subalgebra.
The exact dimensions are 9 and 6. The BO tangent indices of a
complement \(T\) to \(D\) are \(0,1,3,8,9\).

## CU2. The actual analytic field

For \(t\in\mathbb C^5\), \(z\in\mathbb C^9\), let \(T(t)\)
and \(Z(z)\) denote the coordinate combinations of those bases
of \(T\) and \(D\). The map
$$
 g(t,z)=\exp(T(t)/2)\exp(Z(z)/2),\qquad
 (t,z)\longmapsto g(t,z)K_{\mathbb C}
$$
has invertible differential at zero. It gives a local
holomorphic chart by the inverse function theorem. In the gauge
given by \(g\) define the \(\mathfrak k\)-valued five-form
$$
 A=[g,c]\;dt^1\wedge\cdots\wedge dt^5.
$$
For derivatives in \(z\), the Maurer–Cartan form is that of
\(\exp(Z(z)/2)\). Because \(\mathfrak h\) is a Lie algebra,
its skew part lies in \(\mathfrak k_c\). The canonical
connection therefore differentiates \([g,c]\) to zero along
the \(z\)-leaves. The other derivative terms carry one of the
\(dt^i\), so wedge multiplication by the transverse volume
kills them. Hence \(d^\nabla A=0\) and the
\(\mathfrak p\)-valued form \(B=\theta A\) is closed.

The symmetric part of the same Maurer–Cartan form maps the
nine leaf directions isomorphically to \(D\) near zero. In the
\(g\)-frame the transverse volume is consequently a nonzero
scalar times
\(\alpha=e^0\wedge e^1\wedge e^3\wedge e^8\wedge e^9\).
Thus \(B\) is a scalar multiple in that frame of the fixed
tensor \(\theta(c\alpha)\). Its BO covariant coefficients
are \(1,-1,1,1\) in the respective slots
\((0,811),(1,4875),(1,8971),(8,783)\). The
[verifier](../../computation/geometric-unity/CU_h6_local_field.py)
checks the brackets, both hook traces, all coefficient
curvature slots, the weight \((5,5)\) and both highest-weight
equations. The chart argument therefore places the entire
holomorphic field in the complex \((5,5)\) subbundle.

For the real observer, restrict to the real analytic slice
and take compact real parts. The operators, hook conditions
and isotypic subbundle are real. The real part has squared
norm \(3/131072\) at the origin, by the exact compact conjugation
check. This gives the claimed nonzero real analytic closed
field on an actual neighborhood. An arbitrary holomorphic
factor \(f(t)\) may also multiply \(A\) in the construction,
since \(df\wedge dt^1\wedge\cdots\wedge dt^5=0\).

## CU3. Extension to the whole positive observer

Reorder the split defining axes as \((W,e_4,U)\) with
\(U=\operatorname{span}(e_0,e_2)\). Write the symmetric
matrix of a point in the model \(y=gG^{-1}g^T\) as
$$
 y=\begin{pmatrix}A_0&b&C\\ b^T&d&e\\ C^T&e^T&F\end{pmatrix}.
$$
Here \(A_0,F,C\) are two-by-two and \(b\) is a column. The
real observer consists of the transforms \(y=L_0^{-1}Y L_0^{-T}\)
of real positive definite \(Y\) with determinant one, where
the split basis matrix \(L_0\) satisfies \(L_0^T L_0=G\). Complex
conjugation interchanges \(W,U\) while fixing the middle axis.
Consequently \(F=\overline{A_0}\), \(C^T=\overline C\),
and the Hermitian matrix
$$
 M=\begin{pmatrix}C&-A_0\\-F&C^T\end{pmatrix}
$$
is positive definite. Indeed \(y\) multiplied by the conjugation
swap is the Hermitian congruence \(L_0^{-1}Y(L_0^{-1})^\dagger\).
Take its \(W,U\) principal block and change the sign on the
second block to obtain \(M\).

Let \(J=\operatorname{diag}(I_2,-I_2)\) and \(H=JM\).
This matrix is similar to the Hermitian \(M^{1/2}JM^{1/2}\).
Its positive and negative eigenspaces have dimension two each.
The positive eigenspace is \(J\)-positive, since eigenspaces
of different eigenvalues are \(J\)-orthogonal and an
eigenvector with eigenvalue \(\lambda>0\) satisfies
\(u^\dagger Mu=\lambda u^\dagger Ju\). It is therefore the
graph \(\{(x,Sx)\}\) of a two-by-two matrix \(S\). The matrix
\(H\) is Hamiltonian for the usual complex symplectic form
because \(A_0,F\) are symmetric. Pairing its positive
eigenvectors by that form gives zero, since a pairing is
annihilated by the sum of their eigenvalues. Hence the graph
is Lagrangian and \(S=S^T\). Its invariance gives
$$
 F-SC-C^TS+SA_0S=0,\qquad R=C-A_0S .
$$
The restriction of \(H\) to the graph is conjugate to \(R\).
In particular \(\det R>0\) and the eigenvalues of \(R\) are
positive. Spectral projectors separated from zero are locally
analytic functions of \(H\), for instance by the resolvent
contour formula. The unique positive graph therefore defines
one real analytic \(S(y)\) and \(R(y)\) on the whole observer,
with a holomorphic extension near each real point.

Define
$$
 v=(R^T)^{-1}(e^T-Sb),\qquad
 n(S,v)=\begin{pmatrix}I&0&0\\v^T&1&0\\S&0&I\end{pmatrix}.
$$
Then
$$
 y'=n^{-1}y n^{-T}=
 \begin{pmatrix}A'&b'&R\\b'^T&d'&0\\R^T&0&0\end{pmatrix},\qquad
 d'=\frac{1}{16(\det R)^2}.
$$
The last equality follows from \(\det y'=1/16\).
Use the positive root \(q=2\sqrt{\det R}\), and put
$$
 D_0=\operatorname{diag}(q,1),\quad L=2R D_0^{-T},\quad
 a=q^{-2},\quad \ell=b'/a,\quad M_0=(A'-\ell\ell^T)L^{-T}.
$$
The matrix
$$
 h=\begin{pmatrix}L&\ell&M_0\\0&a&0\\0&0&D_0\end{pmatrix}
$$
satisfies \(hG^{-1}h^T=y'\) and
\(\det L=\det D_0=q\), so \(\det h=1\).
This constructs \(g=nh\) analytically on the entire real observer.
The matrix \(h\) preserves \(W\subset W^\perp\) and its
determinants on \(W\) and on the last quotient agree.
The Lie algebra of this determinant constraint in the
flag-preserving subgroup is \(D\oplus\mathfrak k_c\).
Indeed the \(G\)-adjoint preserves the flag constraint and
interchanges the outer block traces. The selfadjoint part lies
in \(D\), and the skew part has trace zero on \(W\).

The five global transverse functions are
\(t=(S_{11},S_{12},v_1/2,S_{22},v_2/2)\).
Locally their holomorphic differentials have rank five.
The Riccati derivative in a symmetric variable \(Z\) is
\(-ZR-R^TZ\), an invertible Sylvester map because the eigenvalue
sums of \(R\) are positive. The other transverse equation
has invertible coefficient \(R^T\) in \(v\). To see surjectivity
of \(dt\), keep \(h\) fixed at the given point and vary
\((S,v)\) in \(nhG^{-1}h^T n^T\). The transverse equations
hold at these varied values and the invertible Jacobian selects
the same local holomorphic branch. Thus this slice is a local
right inverse for the five functions.

Use \(A=[g,c]\,dt^1\wedge\cdots\wedge dt^5\) again.
Along the kernel of \(dt\), \(dn=0\), so
\(g^{-1}dg=h^{-1}dh\) takes values in \(D\oplus\mathfrak k_c\).
In each local holomorphic extension the corresponding nine
directions map isomorphically to \(D\) in the moving frame.
The proof of CU2 now applies on all local extensions and gives
one globally real analytic closed complex field \(B=\theta A\).
At the identity \(g=I\), \(S=v=0\) and \(C=I_2/2\).
The coordinate derivative gives \(dS=dF\) and \(dv=2\,de^T\).
The BO tangent matrices therefore give exactly the five
transverse coframe axes of CU2. The real part has the same
nonzero base squared norm \(3/131072\).
Entire functions of \(t\) give global multipliers of \(A\)
and \(B\) by the same wedge argument.
The [symbolic block check](../../computation/geometric-unity/CU_h6_global_field.py)
reconstructs the positive Gram congruence, the Riccati derivative,
the determinant-one factorization and the five base covectors.

## CU4. All values and the nonlinear scope

The compact orbit of the nonzero real value spans the single
286-dimensional irreducible \((5,5)\) module. Finitely many
rotated global fields and their linear combinations realize
every value in it. [BR](BR-h6-closed-field-lifts.md) supplies
the complementary 6,534-dimensional module by global analytic
fields.
[BU](BU-h6-second-jet-obstruction.md) removes the final
14-dimensional type from all smooth closed hook values. This
proves the equality of the local analytic value space with the
6,820-dimensional module.

The single-type fields satisfy linear closure. The exact
[BQ](BQ-h6-invariant-moments.md) grade-five invariant on
that type is \(L_5=-65\|B\|^2/12\), which excludes any
nonzero pure \((5,5)\) value from the full equation. The
construction adds this type to the available space for mixed
nonlinear fields. [DB](DB-h6-bounded-analytic-lifts.md) subsequently
obtains bounded global lifts by compact averaging. General
nonlinear feasibility remains open.

## CT1. The local cohomology class

For each degree define the parallel bundle
\(E^q=\ker(R^{T^*}\wedge:T^*Y\otimes\Lambda^q T^*Y\to
T^*Y\otimes\Lambda^{q+2}T^*Y)\).
Parallel curvature and the exterior product rule imply that
\(d^\nabla\) preserves these bundles. On them
\((d^\nabla)^2=R^{T^*}\wedge=0\), so they form a differential
complex.

The [CT verifier](../../computation/geometric-unity/CT_h6_potential_obstruction.py)
shows that \(V^*\otimes E^5\) contains no \((5,5)\) component.
There are 88 raw columns of weight \((5,5)\) in
\(V^*\otimes V^*\otimes\Lambda^5 V^*\). Impose curvature
on each derivative slot and the two highest-weight operators
on all slots. A selected maximal minor is \(61{,}217\) modulo
the prime \(65{,}537\). Its denominators are invertible.
Thus the simultaneous highest-vector kernel is zero over
\(\mathbb Q\) and \(\mathbb C\). No hook traces are imposed
on the degree-five input.

For any section \(H\) of \(E^5\), parallelness puts its
derivative in \(V^*\otimes E^5\). Equivariance of exterior
alternation then gives zero \((5,5)\) component of
\(d^\nabla H\). Any potential of a closed form already
satisfies \(R^{T^*}\wedge H=(d^\nabla)^2H=0\).
The real closed form in CU2 therefore defines a nonzero germ
of cohomology in degree six. The auxiliary \(A\) has
\(\mathfrak k\) coefficients and enters through the soldering
bracket, which is consistent with this obstruction.

    python -X utf8 computation/geometric-unity/CU_h6_local_field.py
    python -X utf8 computation/geometric-unity/CU_h6_global_field.py
    python -X utf8 computation/geometric-unity/CT_h6_potential_obstruction.py

[Review105](../../reports/reviews/REVIEW105-GU-CU-CT-2026-11-14.md)
records the independent computations. External peer review
remains outstanding.
