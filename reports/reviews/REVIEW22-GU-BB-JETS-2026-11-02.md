# Internal review of BB's higher-jet extension

November 2, 2026. Review of BB3–BB6 in `formalization/geometric-unity/BB-weyl-direction.md` and `computation/geometric-unity/BB_cubic_extension.py`, following `eb22da1`. The extension is valid as stated. Review21 remains the review of the auxiliary quadratic certificates; it has not been broadened retrospectively. This is internal review, not external peer review.

## Higher-jet invariance

The support arguments of Z do not require zero Weyl curvature. Fixing any quadratic metric jet fixes the linear base-connection term. Its quadratic term is linear in the cubic metric coefficients, and its cubic term is independent of those coefficients and affine in the quartics. Hence the observer first jet is fixed, the second is affine in the cubics, and the third is affine jointly in cubics and quartics. No product of two variable cubic terms enters these orders.

Only mixed metric second derivatives with two horizontal derivative indices vary with the cubics. For a diagonal vertical slot, the leading degree-two metric variation is a base gradient. Commuting its third derivatives eliminates the associated three-horizontal Riemann components. The remaining cubic curvature variation has three horizontal slots and one off-diagonal vertical slot; Ricci variation is horizontal/off-diagonal-vertical. The inverse metric preserves this support at the two diagonal witnesses.

Thus the 722 coefficient rows and the three compatibility rows remain fixed. Their selected RHS rows have even base-reflection parity, so their affine cubic variations vanish. The leading degree-three quartic variation is also a gradient in the diagonal-vertical slot; its principal Ricci-derivative terms cancel by commuting fourth derivatives. The companion vertical derivative vanishes at the origin. These are identities of the relevant Taylor jets for the fixed quadratic curvature. A smooth zero-four-jet remainder affects none of them.

## Kernel identities and the generic rank step

At each of the two witnesses and each sign, the four sparse columns

$$K_i=e_{28+i}+\frac{24d_i}{5}e_{14+i}$$

annihilate all 722 invariant rows at three distinct amplitudes. BB1's already proved degree-two coefficient bound makes this an identity for every real amplitude. The columns remain independent, so the invariant rank is at most 38. This is not an assertion that sampled ranks determine a universal rank.

Outside the finite exceptional sets from the quadratic calculation, choose a witness with both its compatibility polynomial and an existing coefficient 42-minor nonzero. Its first 38 rows are independent; their kernel is therefore exactly the span of $K$.

In an appended row $(i,O,O)$, $O=5$ or $6$, a cubic change can affect only the coefficient of $\xi_O$. All columns of $K$ have $\xi_O=0$, so this row change annihilates $K$. The annihilator is exactly the span of those 38 invariant rows. Subtracting that row combination restores the baseline appended row and preserves the same determinant for every cubic tensor. Quartics change no coefficient entries. This proves the needed rank transfer at the selected witness without assuming a global determinant identity or rank 38 at both witnesses for every amplitude.

## The twelve exceptional amplitudes

At every member of the two six-element exceptional sets, the full 722-row subsystem at the first witness has exact rank 38. Together with the polynomial kernel identity, this proves that its kernel is $K$, even if the original 38-row list becomes dependent. A fresh independent subset can be chosen from all 722 rows.

The completion restrictions are exactly

$$J_5K=-(3w+\sigma)^2\Delta_{j,5},\qquad
J_6K=-(3w-\sigma)^2\Delta_{j,6},$$

with the recorded strictly positive diagonal matrices. Their signs and constants match the output. Since $\sigma=\pm1$, the two scalar factors cannot vanish simultaneously. At least one completion is invertible on $K$, and cubic row changes do not alter that restriction. This gives coefficient rank 42 for every higher jet at all twelve remaining amplitudes. Their first-witness compatibility polynomials are already nonzero, so augmented rank is 43 there as well.

The generic selected-minor argument and these finite checks cover all amplitudes. No untested rank limit at an exceptional value is used.

## Smooth metrics and intrinsic scope

Normalization sends the Weyl amplitude to $\omega/|\rho|$, cubics to $Q_3/|\rho|^{3/2}$ and quartics to $Q_4/\rho^2$. Their arbitrary ranges and the zero-four-jet remainder are preserved.

Normal coordinates at a point with curvature $\rho R_0+\omega\mathcal W$ have precisely the prescribed quadratic term. Their actual cubic and quartic Taylor coefficients are now allowed, so BB6 genuinely removes the polynomial-metric restriction. It requires **one orthonormal frame simultaneously** giving the stated rank-one Ricci axis and Weyl tensor. It is not a condition on independent diagonalizations of the two tensors, nor a classification of arbitrary Weyl curvature.

For each fixed actual metric, full coefficient rank and inconsistency give a nonzero augmented 43-minor. Continuing that fixed minor at its positive coordinate fibre matrix, and retaining nonzero base curvature nearby, supplies the U5 argument. The original compatibility weights need not persist. The full-field conclusion retains S.0's compact $(14,0)$ setting, Z1 at $\lambda=1/5$, positive $c_H,c_V$, invariant odd field, slaved vector even field, $\kappa\ne0$ and $\det M\ne0$. No arbitrary even-field classification or parameter-uniform neighbourhood is claimed.

## Reproduction and boundary

`python -X utf8 -I computation/geometric-unity/BB_cubic_extension.py` passed all **41 checks** in 33.2 seconds with Python 3.12 and SymPy 1.14.0. Deterministic results matched the recorded output: polynomial kernel identities at both witnesses/signs, both completion formulas, all twelve exceptional invariant ranks, and the mixed cubic/quartic controls. The controls also verified reflection, affine dependence and preservation of the selected generic minor. The script reuses the accepted curvature engine and the reviewed quadratic certificates; it is not a second curvature implementation. The support and rank-nullity reasoning above were checked separately.

The extension places no additional higher-jet restriction inside the declared joint curvature orbit. Weyl tensors outside that orbit, Ricci-flat nonflat points, non-vector even branches and indefinite real forms remain outside this conclusion.
