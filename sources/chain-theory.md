# Sources — Chain Theory (A. D. Wissner-Gross)

All primary sources were obtained from the author's public repository.

**Repository:** `https://github.com/alexwg/chain-theory`, cloned 2026-10-30, commit
`a28db8998bce5bd5ef11195c05f13c86cd9dac5c` (2026-10-24 00:53 UTC). Licenses: Apache-2.0 (code), CC BY 4.0 (text).

| ID | Item | Path in repo | Version / date | sha256 | Status |
|---|---|---|---|---|---|
| CT-S1 | Letter: *A Parameter-Free Candidate for the Fine-Structure Constant from Entangled Qudits* | `papers/alpha/alpha.tex`, `alpha.pdf` | v1.0.0, dated 2026-10-23 | tex `2ed2833b…6fa7c6`; pdf `4bd42b02…22009` | **Primary.** Formal claims: Theorem 1 (Uniqueness), Theorem 2 (Value), Conjecture 1 (Photon). |
| CT-S2 | Supplemental Material to CT-S1 | `papers/alpha/supplement.tex` | same | `8b88c9f2…d10ed9336b722945c331da` | **Primary.** Proofs S1–S12, Table S-inputs, Table S-alternatives. |
| CT-S3 | Lean 4 check of CT-S1 | `papers/alpha/lean/ChainAlpha.lean` | Lean 4.35.0 (pinned in `lean-toolchain`) | `522844e9…c331da` | **Primary (code).** Core library only, no imports. |
| CT-S4 | Record table and enumeration script | `papers/alpha/anc/records.csv`, `enumerate_records.py` | same | csv `f4bfd5c9…abbc9`; py `12292270…dec` | **Primary (data/code).** |
| CT-S5 | Compendium *Chain Theory: Toward Fundamental Physics from Distilling Networks of Entangled Qudits* | `compendium/` (`pdf/compendium.pdf`, `index/`, `tex/`, `lean/`) | release 0.1.0 | pdf `cdd99bf0…bc9` | **Primary**, self-described "working draft". Audited in `reports/chain-theory-compendium.md`. |
| CT-S6 | Slides *Chain Theory: An Introduction* (71 slides) | `talks/introduction.pdf` | — | `a49239e7…540c` | **Secondary summary** by the author. The seed URL `alexwg.org/publications/ChainTheoryTalk.pdf` was not reachable from this sandbox (proxy 403). We have not checked that the two files are identical. |

Full hashes are reproducible with `sha256sum` on the files at the commit above.

## External results the Letter imports (to be checked as lemmas, not assumed)

| Lemma | Reference | Used for |
|---|---|---|
| Every positive map $M_2\to M_2$ is decomposable | Størmer, Acta Math. 110, 233 (1963) | read cap |
| Every positive map between $M_2$ and $M_3$ is decomposable; a non-decomposable positive map $M_2\to M_4$ exists | Woronowicz, Rep. Math. Phys. 10, 165 (1976) | read cap |
| A non-decomposable positive map $M_3\to M_3$ exists | Choi, Linear Algebra Appl. 12, 95 (1975) | read cap |
| Twirl commutant of $U\otimes\bar U$ is $\mathrm{span}\{\mathbb 1,P_d\}$ | M. & P. Horodecki, PRA 59, 4206 (1999) | unique readout |
| $S(U(3)\times U(2))\cong (SU(3)\times SU(2)\times U(1))/\mathbb Z_6$ | Baez & Huerta, Bull. AMS 47, 483 (2010) | gauge group |
| CODATA 2022: $\alpha^{-1}=137.035\,999\,177(21)$ | Mohr et al., RMP 97, 025002 (2025) | comparison |
| $\alpha^{-1}=137.035\,999\,166(15)$ (electron $g-2$); $137.035\,999\,206(11)$ (Rb); $137.035\,999\,046(27)$ (Cs) | Fan et al. PRL 130, 071801 (2023); Morel et al. Nature 588, 61 (2020); Parker et al. Science 360, 191 (2018) | comparison |

These are standard published theorems and measurements. We use them as stated and do not re-prove them.
