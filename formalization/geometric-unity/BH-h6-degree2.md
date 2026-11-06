# BH. Degree-two H6 values have no scalar grade-one quadratic

On the flat positive Z1 observer at $\lambda=1/5$, the new 14-dimensional degree-two summand of [BG's pointwise curvature kernel](BG-h6-fixed-values.md) contains no nonzero real value with scalar grade-one quadratic. This treats every value in that summand, rather than only its $SO(4)$-fixed line. It does not treat arbitrary sums mixing this summand with the harmonic-degree-three and degree-five summands.

Evidence: [exact coefficient helper](../../computation/geometric-unity/BH_h6_degree2.py), [output](../../computation/geometric-unity/BH_h6_degree2.out) and [BG's tensor certificate](../../computation/geometric-unity/BG_h6_curvature_kernel.json). The hand proof supplies diagonalization, the complete diagonal section and the spectrum lemma. Internal [Review70](../../reports/reviews/REVIEW70-GU-BG-BH-2026-11-05.md) and [Review71](../../reports/reviews/REVIEW71-GU-BG-BH-2026-11-05.md) accept BH1–BH3 and disclose shared contributions and libraries. External peer review remains outstanding.

## Setting

Use a flat positive four-base, reconstructed Z1 at $\lambda=1/5$ and constant $c_H,c_V>0$. At normalized scales $c_H=2,c_V=1$, [AM](AM-global-flat-models.md) identifies its tangent with
$$V=\operatorname{Sym}^2_0(\mathbb R^5),\qquad K=SO(5).$$
All slots of the real hook $B_6\in H_6\subset V^*\otimes\Lambda^6V^*$ are covariant. Let $v_{10}=K_2$ be BG's explicit nonzero, $SO(4)$-fixed Casimir-ten value. Its $K$ orbit span is the real 14-dimensional degree-two module $W_{10}$. This is a pointwise submodule of the coefficient-curvature kernel; no covariantly closed field lift is assumed.

For the field consequence retain [BA's hypotheses](BA-hook-closure.md): true dualized bosonic equation, fixed parallel real $W_A$, signature $(14,0)$, compact $\mathfrak u(128)$, real smooth fields, constant $\kappa\ne0$, and the pure pair
$$T_o=i(3\kappa/104)\Phi^1,\qquad u=0.$$
The even equation then gives exactly $T_e=iB_6\Gamma$. No coincident even kernel occurs at this pair.

## BH1. A complete diagonal section of the module

BF/BG's sphere argument identifies $W_{10}$ equivariantly with $\operatorname{Sym}^2_0(\mathbb R^5)$. Normalize the intertwiner $\iota$ by
$$\iota(S_4)=v_{10},\qquad S_j=e_je_j^T-I/5.$$
Let $A_i$ be BC's induced covariant Nomizu action along $H_i$. The calibrated action is $A_i=\rho(J_{i4})/2$, with $J_{i4}=E_{i4}-E_{4i}$. Directly in the matrix representation,
$$\rho(J_{i4})^2S_4=2(S_i-S_4).$$
Therefore the actual quarter rotations of the zonal seed are
$$v_i=v_{10}+2A_i^2v_{10}=\iota(S_i),\qquad i=0,\ldots,3.$$
The helper checks the full-component identity $A_i^3v_{10}=-A_iv_{10}$, the five-term sum $\sum_{i=0}^4v_i=0$, all hook conditions and the full ten-term Casimir eigen-identities. The first four values are independent: their Gram matrix has diagonal $18225/16$ and off-diagonal $-18225/64$.

For $z\in\mathbb R^4$, the tensor $\sum_{i<4}z_iv_i$ corresponds to the diagonal matrix with
$$s_i=z_i-\bar z\ (i<4),\qquad s_4=-\bar z,\qquad \bar z=\tfrac15\sum_{i<4}z_i.$$
Conversely every trace-free diagonal matrix occurs, with $z_i=s_i-s_4$. This is the whole diagonal section of $W_{10}$.

## BH2. Ten exact moments and an elementary spectrum lemma

For a diagonal $S=\operatorname{diag}(s_0,\ldots,s_4)$, let $q_{ij}$ be the metric-relative value of $Q_1(i\iota(S)\Gamma)/i$ on the tangent direction $E_{ij}+E_{ji}$. These are diagonal stress values on ten off-diagonal **tangent** directions. They are necessary scalarity tests, not a claim that ten entries determine the entire stress tensor.

The helper computes all four pure-amplitude and all six mixed coefficients for every one of these ten moments. Exact polynomial identities give
$$q_{ij}=a(s_i^2+s_j^2)+bs_is_j+c\operatorname{tr}(S^2),$$
$$\boxed{a=-2375/32,\qquad b=1125/16,\qquad c=-3625/2.}$$
The normalization uses the physical field $iB_6\Gamma$: evaluate on $B_6E$, where $\Gamma=\sqrt{\det g}\,E$, and multiply the quadratic by $\det g$. No extra chirality sign is inserted. The output records all ten four-amplitude polynomials. A separate sorted-generator, nested-gradient implementation checks ten targets on a generic exact mixture; these are convention checks for the same cubic formula.

**Spectrum lemma.** If five real $s_i$ sum to zero, $b\ne0$, $b-3a\ne0$, and all the displayed $q_{ij}$ coincide, then all $s_i=0$.

*Proof.* If the $s_i$ all coincide, trace zero gives the conclusion. Otherwise choose $s_j\ne s_k$. For every $i$ outside this pair,
$$0=q_{ij}-q_{ik}=(s_j-s_k)\big[a(s_j+s_k)+bs_i\big].$$
Since $b\ne0$, the other three values equal some $r$. If $s_j\ne r$, apply the same argument to $j$ and one of those three indices: the remaining values include $s_k$ and two copies of $r$, so $s_k=r$. The symmetric alternative gives the same result. Thus every nonconstant candidate has one exceptional entry $s$ and four copies of $r$.

Trace zero gives $s=-4r$. Comparing an exceptional/ordinary pair with an ordinary/ordinary pair yields
$$0=q_{sr}-q_{rr}=(s-r)[a(s+r)+br]=-5(b-3a)r^2.$$
Hence $r=s=0$. This proves the lemma. Here $b-3a=9375/32\ne0$. $\square$

Every real symmetric trace-free matrix is diagonalizable by $SO(5)$: changing one eigenvector sign repairs the orientation of an orthogonal eigenbasis. Both $\iota$ and $Q_1$ are $SO(5)$-equivariant. If $Q_1$ were scalar at a nonzero module value, it would remain scalar after diagonalization, making all ten necessary moments equal. The lemma contradicts nonzero $S$. Therefore
$$\boxed{B_6\in W_{10},\quad Q_1(iB_6\Gamma)/i\text{ scalar}\ \Longrightarrow\ B_6=0.}$$

## BH3. Full-field consequence and boundary

BA1 includes derivatives of every coefficient and form slot for arbitrary smooth hook fields. On this pure branch its derivative has only grade seven; variable amplitudes cannot create a grade-one term. The full grade-one equation is
$$Lg_Y+Q_1/i=0,\qquad L=\frac{15}{2c_V}+\frac{15}{52}\kappa^2>0.$$
Every nonzero value in $W_{10}$ fails scalarity by BH2, while its zero value fails $L>0$. Consequently a full solution on this branch cannot take a value in $W_{10}$ even at one point. No closure lift or constancy of the field amplitudes is used.

BC's affine-conjugating horizontal dilation reduces arbitrary positive scales to a positive homothety of the normalized model. It preserves the module and scalarity. For a fixed seven-covariant tensor, all metric-relative quadratic moments acquire the common positive factor $c_V^{-7}$. The positive actual baseline retains its displayed scale dependence. Thus the result holds for all $c_H,c_V>0$ and either sign of $\kappa$ at this fixed $\lambda$.

Arbitrary mixtures with the 30- and 91-dimensional summands, the whole ambient curvature kernel, closed-field reconstruction, curved bases, other parameters, other hooks and coincident factors remain separate questions. This is a conditional mathematical obstruction, not a physical test or a literature-priority claim.

## Reproduction

Run `python -X utf8 -I computation/geometric-unity/BH_h6_degree2.py` from the repository root in the pinned environment. It reads the adjacent BG certificate, imports BC/BF and the separate Clifford library without their mains, and writes stdout only. No private receipt, earlier output or previous helper run is required. Exact polarization supplies the polynomial certificate; the spectrum lemma and equivariant diagonalization supply its universal consequence.
