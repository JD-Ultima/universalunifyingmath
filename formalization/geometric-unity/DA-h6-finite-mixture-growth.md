# DA. Degeneration of finite isotropic H6 mixtures

Theorem. Let B_j be compact SO5 rotates of the normalized global CU
complex field, and let t_j be their five global complex transverse
coordinates. Let f_j be arbitrary entire functions on C5. Let G be
any globally bounded field. If
G + sum_j Re(f_j(t_j) B_j)
is globally bounded, then the finite sum of the isotropic fields vanishes identically.
More precisely, after grouping terms by their unoriented center line
a_j in R5, every such group is identically zero. The theorem does not
assert linear independence of fields within one group.

For a bounded globally real analytic [BW](BW-h6-triangular-orbit-averages.md) field G, a finite mixed
ansatz cannot satisfy the full nonlinear equations even on a nonempty
open set: its fixed norm identity would extend to all of connected Y,
DA would reduce the field to G, and the existing BW obstruction applies.
This is a statement about this verified global ansatz, not an exclusion
of arbitrary local or globally bounded (5,5) fields. Infinite CU measures
are not covered; estimates are not uniform near coincident centers.
## A global pointwise CU norm formula
Fix the CU center a, write the real positive metric as
Y=[[A4,b4],[b4^T,d4]] relative to a^perp plus a, and use the CU
split coordinates and positive spectral graph S. Set
R=C-AS, q=2sqrt(detR), Delta=I-S* S.
Here * is conjugate transpose; Delta is positive and Delta<=I.
The two eigenvalues lambda1,lambda2 of R are positive real. One has
q^4=det A4. For the normalized CU field,

 ||Re(f(t)B)(Y)|| / ||Re b0||
   = |f(t(Y))| q^(-5) [q/tr R] det(Delta)^(5/2).

In particular q/trR=2sqrt(lambda1 lambda2)/(lambda1+lambda2)<=1,
and det Delta<=1, so
 ||Re(f(t)B)(Y)|| <= ||Re b0|| |f(t(Y))| (det A4)^(-5/4).
The exact base norm is ||Re b0||^2=3/131072.

Here is a derivation that does not rely on the special axial example.
In the CU complex frame g=n h, the coefficient is c and the transverse
volume is J0 alpha. The three S derivatives are obtained by inverting
Z -> R^T Z+ZR on Sym^2(C2), whose determinant is 4detR trR=q^2 trR.
Their numerator is congruence by the lower diagonal block D of h,
whose determinant on Sym^2 is (detD)^3=q^3. The two v/2 derivatives,
modulo the S rows, have matrix (a/2)R^(-T)D with a=q^(-2), so their
determinant is q^(-5). Thus
J0=(q/trR)q^(-5)=q^(-4)/(trR).
The upper off-diagonal blocks of h only alter off-diagonal blocks in
this five-by-five derivative matrix, not its determinant.

It remains to convert the highest tensor norm to a real frame.
Put U=[I;S], J=diag(I,-I), and M=JH>0 as in CU. The graph equation
is MU=JUR, so U* J M^(-1) J U=Delta R^(-1).
If E4 consists of the physical split columns W,U, then
M=J E4* A4 E4 J/4. The Gram matrix of
Z=A4^(-1/2) E4 U
is exactly Delta R^(-1). In an arbitrary real orthonormal frame, Z
spans the adapted isotropic plane. Normalizing its two columns to the
standard compact isotropic basis gives determinant modulus
sqrt(det(Delta)/detR)/2. The factor q from det(h|W) makes the
coefficient c scale by sqrt(det Delta).
The BO highest b0 has determinant-line weight5 under the stabilizer
of W. Thus its norm scales by det(Delta)^(5/2). This character is
checked directly by all seven generators of the W-stabilizer in k:
K.b0=5 trace(K|W)b0. Combining with J0 proves the norm formula.
The complex highest tensor has bilinear norm zero. Hence real parts
have norm exactly 1/sqrt2 of the Hermitian norm, independently of the
phase of f; this proves the displayed formula for arbitrary complex f.
Finally detM=det A4/16 and the Hamiltonian eigenvalues are
lambda1,lambda2,-lambda1,-lambda2, proving q^4=det A4.
## Simultaneous growth of an entire same-center group
Fix any positive real base metric Y0 and a center a. Choose a real
determinant-one frame of block Cholesky form
F0=[[L,0],[w^T,delta]], delta>0, Y0=F0 F0^T.
For rho>0 set
F_rho=F0 diag(rho I4,rho^(-4)), Y_rho=F_rho F_rho^T.
The A4 and b4 blocks are multiplied by rho^2; their Schur complement
is multiplied by rho^(-8). For every CU choice W with this center,
S and v are unchanged, R is multiplied by rho^2, and q by rho^2.
Therefore every entire multiplier f(t) is unchanged on this path.

This gives tensor equality, not only equality of norms. To see this,
the isotropic vectors nW have center coordinate
v^T=b4^T A4^(-1) E4 U.
Consequently their coordinates in F_rho are rho^(-1) times their
coordinates in F0, with zero last component. Multiplication by q
cancels the two factors rho^(-1), so the associated coefficient c
is the same in these common real frames. Its highest determinant
line is also the same, with the same phase. Meanwhile J0 scales by
rho^(-10). Thus for every field in this center group,

 F_rho^* Re(f(t)B)(Y_rho) = rho^(-10) F0^* Re(f(t)B)(Y0).

The equality holds for both choices of orientation/chirality and
commutes with real conjugation. It therefore holds for their finite
sum, even when individual f(0) vanish or base values cancel.
## Other-center fields decay on the same path
The same path has the form
Y_rho=rho^2 Y0+(rho^(-8)-rho^2)delta^2 aa^T.
For another center a' not parallel to a, let Q4 be the positive
compression of Y0 to (a')^perp and let z be the nonzero projection
of delta a to (a')^perp. This also covers perpendicular centers.
The determinant lemma gives
det(Y_rho|(a')^perp)=rho^8 det(Q4)
                         (1+(rho^(-10)-1)z^T Q4^(-1)z)
which is asymptotic to a positive constant times rho^(-2).
The pointwise norm bound therefore gives O(rho^(5/2)) before the
multiplier, uniformly for a fixed base and fixed other center.

The multiplier stays bounded on this path. S remains contractive.
Write the cross block as rho^2 b0+(rho^(-8)-rho^2)z beta.
Sherman-Morrison shows
(Q4+(rho^(-10)-1)zz^T)^(-1)(b0+(rho^(-10)-1)z beta)
is bounded as rho approaches zero. Since v^T=b4^T A4^(-1)E4[I;S],
v is bounded as well. Thus all five transverse coordinates stay in
a compact subset of C5. Every entire f is bounded on that compact
set, so the full other-center term is O(rho^(5/2)).
## Elimination of finite groups
If a same-center group is not identically zero, choose a base Y0
where its value is nonzero.
Its own-center contribution grows exactly as rho^(-10), while every
other-center CU contribution tends to zero and G stays bounded.
This contradicts boundedness of the total sum. Therefore every
same-center group vanishes identically. This argument does not need
to classify internal cancellations or establish any spectral
independence of CU fields.

The entire-multiplier hypothesis can be weakened to holomorphic
multipliers bounded on the compact coordinate sets approached by
these degenerations. No such weakening is used in the theorem.

## Reproduction

    python -X utf8 computation/geometric-unity/DA_h6_mixed_growth.py

[Review107](../../reports/reviews/REVIEW107-GU-DA-DB-2026-11-14.md)
records the independent checks. External peer review remains outstanding.
