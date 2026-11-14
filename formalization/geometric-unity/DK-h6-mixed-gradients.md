# DK. First-order equations for the type \((5,5)\) part of a mixed H6 field

Let \(B\) be a smooth closed H6 hook on the positive trace-metric
observer \(Y=SL(5,\mathbb R)/SO(5)\). By
[BU](BU-h6-second-jet-obstruction.md) its values take place in the
6,820-dimensional admissible bundle. Put \(b=\Pi_{55}B\), with the
parallel orthogonal projection to the type \((5,5)\) subbundle.
For a highest weight \(lm\), let \(G_{lm}b\) be the corresponding
projection of \(\nabla b\).

The universal first-order conditions are
\[
 G_{55}b=G_{54}b=G_{53}b=0,qquad
 \|d^\nabla b\|^2=\tfrac{26}{27}\|G_{64}b\|^2,qquad
 \|\delta^\nabla b\|^2=\tfrac{13}{27}\|G_{64}b\|^2.
\]
The nonzero real characteristic directions of the three combined
gradient equations form one projective compact orbit. Its diagonal
representative is \(\operatorname{diag}(1,1,-1,-1,0)\), with a
four-dimensional real kernel. The [characteristic proof](DK-h6-characteristic-orbit.md)
classifies the whole trace-free symmetric cotangent space.
[DL](DL-h6-balanced-characteristic-lift.md) realizes a characteristic
value by a global analytic mixed closed hook.

## DK1. The exact four cancellation tests

[DE](DE-h6-elliptic-periodic.md) gives
\[
 V_{14}\otimes E_{55}
 =E_{75}\oplus E_{65}\oplus E_{64}\oplus E_{55}\oplus E_{54}\oplus E_{53}.
\]
Exterior multiplication and contraction vanish on \(E_{75}\) and
\(E_{65}\). For the remaining four derivative types, compare the
exterior image ranks with and without the \(E_{55}\) fiber column.

| Derivative type | Other source fibers | All columns | Other columns |
| --- | --- | ---: | ---: |
| \((6,4)\) | \((7,3),(6,3),(6,2)\) | 1 | 1 |
| \((5,5)\) | None | 1 | 0 |
| \((5,4)\) | \((6,3),(5,2)\) twice | 1 | 0 |
| \((5,3)\) | \((7,3),(6,3),(6,2),(5,2)\) twice, \((5,1)\) | 5 | 4 |

The verifier constructs the B2 Racah–Speiser tensor inventories,
checks the raw hooks and coefficient curvature on every derivative
slice, and checks the weight, raising relations and fiber Casimir.
It verifies independence of all displayed highest copies, including
the two copies of \((5,2)\). The inventories identify all contributing
copies. The ranks are over the rationals in the BO split coordinates.

For the last three rows, the exterior image of the \((5,5)\) fiber
column is independent of the other columns. Closure therefore
forces those three gradient components to vanish. For \((6,4)\)
there is one common exterior image line in the highest-weight
space. The verifier uses the positive compact Hermitian form to
compute the source and target norms. The squared norm multipliers
of exterior multiplication and contraction on the pure \((5,5)\)
source copy are \(26/27\) and \(13/27\), respectively. Equivariance
gives the pointwise identities on the entire isotypic channel.

The squared exterior multipliers of the three lower source
fibers \((7,3),(6,3),(6,2)\) are \(221/198,35/18,980/3069\).
Their sum is \(943/279\). For a closed mixed field, the
Cauchy–Schwarz bound on the three contributions gives
\[
 \|d^\nabla b\|^2\le\tfrac{943}{279}
 \left(\|G_{64}B_{73}\|^2+\|G_{64}B_{63}\|^2+\|G_{64}B_{62}\|^2\right).
\]

## DK2. The isolated \((5,5)\) gradient

Write \(J_{ij}=E_{ij}-E_{ji}\), and let \(\rho\) be their compact
action on \(E_{55}\). For real trace-free symmetric \(\xi\) set
\[
 Q(\xi)=\sum_{i,k,j}\xi_{ik}\rho(J_{ij})\rho(J_{kj}),
 \qquad J_{ii}=0 .
\]
It is real, selfadjoint and equivariant. On a diagonal covector it
equals \(\sum_{i<j}(\xi_i+\xi_j)\rho(J_{ij})^2\).
The associated map \(T\colon V_{14}\otimes E_{55}\to E_{55}\) is
nonzero by the calibration below. Since this fiber appears once in
the source, \(T\) is a scalar times the orthogonal projection after
an isometric identification. Thus the universal equation
\(G_{55}b=0\) is equivalently
\[
 \sum_a Q(P_a)\nabla_{P_a} b=0
\]
in a trace-orthonormal tangent frame. Its parallel selfadjoint
symbol makes this first-order operator formally skew-adjoint.

For \(\xi_*=\operatorname{diag}(1,1,1,1,-4)\), use the
complex spinor-polynomial model \(E_{55}^{\mathbb C}=\operatorname{Sym}^{10}\mathbb C^4\).
Under \(\operatorname{Spin}(4)\) the four spinor coordinates split
into two two-dimensional half-spin representations. Splitting
degree ten into \(p\) and \(10-p\) gives components of dimension
\((p+1)(11-p)\). The compact \(SO(4)\) Casimir on each is
\((p-5)^2+35\), and the \(SO(5)\) Casimir is 70.
Consequently \(Q(\xi_*)=3C_{SO(5)}-5C_{SO(4)}\) has eigenvalue
\(35-5(p-5)^2\). Its positive and negative multiplicities are
170 and 116. Both axial directions have invertible symbol.

The trace of \(Q(\xi_*)^2\) is 386,100 and \(\|\xi_*\|^2=20\).
Equivariance gives \(TT^*=945I\), so
\(\|G_{55}(\xi\otimes v)\|^2=\|Q(\xi)v\|^2/945\).
In a unit axial direction its least squared singular value is
\(1/189\).

## DK3. A quantitative finite-energy condition

The [DE](DE-h6-elliptic-periodic.md) Weitzenböck identity applies to
arbitrary \(b\) in this parallel type. Its curvature term is
\(-7b\). If \(b\) is globally \(H^1\) and square-integrable,
integration by parts on the complete space gives
\[
 \int\|\nabla b\|^2=7\int\|b\|^2+
 \tfrac{13}{9}\int\|G_{64}b\|^2 .
\]
Combining the [DH](DH-h6-pure-type-spectral-gap.md) bundle bound
\(\nabla^*\nabla\ge29/4\) with DK1 gives
\[
 \int\|G_{64}b\|^2\ge\tfrac9{52}\int\|b\|^2,quad
 \int\|d^\nabla b\|^2\ge\tfrac16\int\|b\|^2,quad
 \int\|\delta^\nabla b\|^2\ge\tfrac1{12}\int\|b\|^2 .
\]
The combined energy in the three lower \((6,4)\) channels is
at least \((93/1886)\int\|b\|^2\) when those channels are in
\(L^2\). A full solution of constant norm need not meet these
integrability hypotheses.

## DK4. Compactly supported projected fields vanish

Suppose the support of \(b\) is compactly contained in the open
set on which \(B\) is smooth and closed. Fix a nonzero real vector
\(a\) and use \(\chi(Y)=\log(a^TYa)\). In a frame \(Y=gg^T\)
its trace-metric gradient is \(uu^T-I/5\) with
\(u=g^Ta/\|g^Ta\|\). The gradient is nonzero and has constant
axial eigenvalue pattern \((4,-1,-1,-1,-1)/5\). Thus the analytic
system of DK2 is noncharacteristic on every level of \(\chi\).

If the support were nonempty, \(\chi\) would attain its maximum
there at some \(x\). The section is zero on the side above that
level and has zero Cauchy data on the level through \(x\).
Holmgren uniqueness for square linear analytic systems gives
vanishing near \(x\), contradicting the definition of support.
The argument places no support condition on the other fiber types.
[Klainerman, Lecture Notes 2010, equation (216) and Theorem 4.1](https://web.math.princeton.edu/~seri/courses/theNotes.pdf).

    python -X utf8 computation/geometric-unity/DK_h6_mixed_gradient.py

The verifier independently checks the raw highest tensors and
reconstructs the ranks, norm multipliers and axial calibration.
