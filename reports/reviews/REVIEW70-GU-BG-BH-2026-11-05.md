# Internal review of the H6 fixed-value and degree-two obstructions

2026-11-05. **BG1–BG3 and BH1–BH3 are accepted under their stated hypotheses.** I read both frozen proofs, metadata, the three helpers/outputs and the explicit tensor certificate. Three fresh public-helper runs passed; every non-timing line matches its author output:

| Helper | Fresh runtime | Matching non-timing lines |
|---|---:|---:|
| `BG_h6_curvature.py` | 7.648 seconds | 32 |
| `BG_h6_moments.py` | 2.876 seconds | 34 |
| `BH_h6_degree2.py` | 6.442 seconds | 27 |

No mathematical correction or theorem narrowing is needed. This is internal AI-assisted review, not external peer review. The reviewer authored the earlier private moment calculation and its public BG moment helper, contributed earlier BC tensor work and the prior spherical-module argument, and reviewed their prior certificates. The fresh runs share those public libraries. This is a disclosed collaborative proof/code assessment, not blind independent discovery or a second independent full geometry engine.

## Complete invariant curvature restriction

BG reconstructs the complete 30-dimensional raw invariant space and 23-dimensional hook space, rather than restricting the domain to BC's closed line. Hodge duality identifies the complete invariant degree-eight target; its nonzero 30-coordinate determinant supplies an isomorphism.

The integer matrices $Q_a=5\Phi_a$ have trace metric $25g$. Their full ambient curvature coefficients correctly divide the fourfold commutator pairing by $2500g_d$. All six stabilizer actions intertwine, and metric skewness, Bianchi and Ricci normalization are checked. Neither mixed/horizontal channels nor the external coefficient slot is omitted.

On covectors the curvature is $-R^T$. Both implementations of its wedge use the correct exterior signs. The selected-coordinate matrix has rank 20 and the recorded nonzero minor. Three independent nullvectors are reconstructed and checked against the **entire sparse curvature output**, establishing exact kernel dimension three independently of a deficient sampled rank. The public JSON is reconstructed and compared, not used to select a desired kernel.

The full-component Casimir action has eigenvalues 10,18,40. BF's matrix-coefficient/sphere argument identifies the three cyclic modules with harmonic degrees two, three and five. Their distinct Casimir eigenvalues make their sum direct, of dimension $14+30+91=135$. This is a contained module and a lower bound on the full kernel, not an exact ambient-kernel classification or a closed-field lift of the new direction.

## Quadratic obstruction and the any-point consequence

The physical field is $iB_6\Gamma$, not real $B_6$ alone. Right multiplication by the coordinate volume followed by the factor $\det g$ correctly implements the unit volume in the quadratic. Every amplitude cross term is retained. Direct norm, weighted trace, recovery of BC's moments, positive principal minors and the explicit positive-square identity agree. Twelve separate nested-Clifford targets test conventions; they implement the same mathematical cubic identity and do not constitute a second derivation of it.

BG's sum of squares makes $q_I-q_H$ strictly positive on every nonzero real vector in the complete fixed three-space. Its full-field consequence requires only a fixed **value at one point**, not an invariant field on a neighborhood: BA differentiates every slot universally, gives closure and coefficient-curvature compatibility, and leaves derivative grade seven only. Thus a varying field still has no grade-one derivative that could repair the scalar-stress obstruction. Zero value fails the positive baseline $\ell$. The final BG3 notation correctly distinguishes this scalar from the stabilizer $L=SO(4)$. Equivariance extends the conclusion to conjugate stabilizers. No averaging of a nonlinear solution is used.

## BH's whole degree-two module

The degree-two identification normalizes an equivariant map from $\operatorname{Sym}^2_0(\mathbb R^5)$ by its nonzero zonal seed. The four quarter rotations are correct: the calibrated infinitesimal action is half the standard rotation, and $A_i^3v=-A_iv$ gives $v+2A_i^2v$ at angle $\pi/2$. Their Gram determinant, five-term relation, hook conditions and **all ten** Casimir terms are checked. The resulting four tensors span the entire trace-free diagonal section.

Exact polarization retains all ten amplitude monomials for each of ten off-diagonal tangent moments. Their full coefficient identities verify
$$q_{ij}=a(s_i^2+s_j^2)+bs_is_j+c\operatorname{tr}(S^2),\quad
a=-2375/32,\ b=1125/16,\ c=-3625/2.$$
These ten necessary scalarity tests are not claimed to reconstruct the whole stress. The spectrum lemma is complete: an unequal pair forces the other three eigenvalues equal; a second unequal pair leaves one exception and four equal entries. Trace zero and comparison of the two pair types give $-5(b-3a)r^2=0$. Both required coefficients are nonzero. Orientation-preserving diagonalization covers every real symmetric trace-free matrix, hence every value in the 14-dimensional module. No closed lift is required for this negative pointwise test.

The positive-scale transfer in both notes is valid: the affine dilation preserves the model and constant homothety preserves the coefficient-curvature kernel under the natural tangent identification. Metric-relative quadratic moments of a fixed seven-covariant tensor have the common positive weight $c_V^{-7}$; scalarity is unchanged. The actual baseline remains strictly positive for either sign of nonzero real $\kappa$.

## Reproduction and trust boundary

Each invocation used `python -X utf8 -I computation/geometric-unity/<helper>.py` with pinned Python 3.12/SymPy 1.14. Transcripts were preserved separately; author outputs were untouched. Only the three new helpers ran, once each in this review. No older main, suite, private reconstruction or full ambient rank calculation was run.

Canonical LF SHA-256 receipts, without BOM; all files below are in `computation/geometric-unity/`:

| Artifact | SHA-256 |
|---|---|
| BG curvature source | 6d8c25a290757afd4576896bb44155ef152bc9957ed7cdd26d35415828380736 |
| BG curvature output | 945e737feb3e6c10e2291a5bdba1a6d74cd5b094b3828e4cd08a1705a5bd1533 |
| BG moment source | 6fa3e73dc5bc4eadf33de10467360ff6ee0133ba5597ca7a5fae92dacc858399 |
| BG moment output | f208ef43ad6c1e5549ddd9fedf3f40311eca76443daf0b6de41f9148a99fa2a5 |
| BH source | 8fb3dc6da9d37377145016e291fb8da073a9bc0e8119131ea19f7ee2abaae098 |
| BH output | c2c1dc889e2fc31d91d5314bfdcc2ace03ef6b351969562254d9e1a310138c55 |
| BG tensor certificate | eac97d439360d4e63468e16ad356f2e92d1e198e46a9f309d7a99157f44482e6 |
| BC library | 9c8871a93e3da243057248a3757f51a39e1dd4710f51108694e50eceb017caba |
| BF library | 72f509389a2292e55a7bb0b7769991ea7b6abbccf2bcf01818b747486a5514c5 |
| Separate Clifford library | 8b2849350623974b2a5a3900e3a7491372dd38ba4c93939c06f9b297cf4c5f11 |

These results retain the flat positive $\lambda=1/5$ setting and the stated pure compact-real $W_A$ branch assumptions for field consequences. They do not exclude arbitrary mixtures of the 14-, 30- and 91-dimensional modules, other ambient curvature components, curved bases, other parameters or coincident kernels. No global differential-Casimir identity, field construction, physical conclusion or literature priority is established. External peer review remains outstanding.
