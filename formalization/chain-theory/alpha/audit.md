# Audit of the Chain Theory α(0) candidate

**Object of study.** Wissner-Gross, *A Parameter-Free Candidate for the Fine-Structure Constant from Entangled Qudits*, v1.0.0 (2026-10-23), with Supplement and `ChainAlpha.lean` (sources CT-S1…S4, commit `a28db89`).
**Claim IDs** refer to `claims/chain-theory-alpha.md`, and gap IDs to `gaps/chain-theory-alpha.md`.
**Code:** `computation/chain-theory/alpha/` (Python) and `proofs/lean/chain-theory/` (Lean). Every number below is printed by one of these scripts, and the outputs are committed beside them.

## Summary of findings

| ID | Finding | Status |
|---|---|---|
| F-CT-1 | The 172 records, sector sizes, closed form and inverse are reproduced (all 51 printed digits of the Supplement; 60 digits of the ancillary `enumerate_records.py`) by an independent operator-level implementation. | Verified (exact) |
| F-CT-2 | On the 169 product records, the two "commuting exchanges" $X$ and $Y$ are the same involution $E$. The weight law is a class function of (E-parity, slot class). | Proven |
| F-CT-3 | Selecting $(2,3)$ by an "acting central charge" depends on a conjugation convention and on a direct-sum trace. Of four natural readings, only the Letter's selects $(2,3)$ alone. Lean confirms that uniqueness fails under the conjugate convention. | Proven (Python + Lean) |
| F-CT-4 | Freedom analysis. The 12 weight-law alternatives of Table S5 were reproduced; its 3 other rows (counting unit, normalisation by 172, other configuration) lie outside the declared space. In a declared 46 080-rule space, the expected number of chance hits within 1σ is $4\times10^{-6}$. Among $1.6\times10^9$ small-coefficient elements of $\tfrac{1}{31104}\mathbb Z[\sqrt2,\sqrt3]$, the author's is the only one within 1σ (0.02 expected by chance). | Computed |
| F-CT-5 | `ChainAlpha.lean` builds under Lean 4.35.0 with standard axioms only. It certifies arithmetic of the definitions. The physically load-bearing choices (CT-G01–G05) enter as definitions. | Verified |
| F-CT-6 | Deviations from the measurements are reproduced. The candidate is algebraic of degree 4; its minimal polynomial is given. | Verified |
| F-CT-9 | The imported read-cap lemma is cross-checked numerically: Choi's map is positive and non-decomposable (SDP), also after embedding. | Computed |
| F-CT-10 | The Letter's exchange $E$ does not commute with the symmetry group $G_\Omega$ of its own joint state, and its record span is not $G_\Omega$-invariant. The covariant version $E'$ gives identical counts and value. Symmetry leaves a 10-parameter family of readouts. | Proven (numerical, exact structure) |
| F-CT-8 | Exploration: the record weights are not the vacuum amplitudes or norms of the record operators. | Computed (negative) |
| F-CT-7 | The construction's formula $f(d_Q,d_q)$ at other dimension pairs gives values far from $\alpha^{-1}$. $(3,2)$ is the only pair within 1% of $\alpha^{-1}$; the next closest, $(2,3)$, is 2.5% away. | Computed |

---

## F-CT-1 Independent reproduction

**Setting.** $H=\mathbb C^3_{Q_1}\otimes\mathbb C^2_{q_1}\otimes\mathbb C^3_{Q_2}\otimes\mathbb C^2_{q_2}$, $\dim H=36$. Letters are the Gell-Mann matrices $\lambda_a$ at $Q_i$ and the Pauli matrices $\sigma_\mu$ at $q_i$. $E(A)=SAS^{\top}$, where $S$ is the permutation unitary $(Q_1,q_1,Q_2,q_2)\mapsto(Q_2,q_2,Q_1,q_1)$.

**Proposition 1.**
(i) There are 169 two-slot products. Their slot-pair sizes are $24,64,24,24,9,24$.
(ii) $E$ permutes them, fixes 11, and pairs 158.
(iii) The three operators $D_{\mu\nu}=[\sigma_\mu,\sigma_\nu]^{(q_1)}-[\sigma_\mu,\sigma_\nu]^{(q_2)}$ satisfy $E(D)=-D$.
(iv) The 172 operators $\{p\pm E(p)\}\cup\{\text{fixed}\}\cup\{D_{\mu\nu}\}$ are pairwise Hilbert–Schmidt orthogonal.
(v) Under the Letter's survival rule the sector sizes are $(+,+)\!:\!79$, $(+,-)\!:\!48$, $(-,+)\!:\!42$, $(-,-)\!:\!3$, and
$$\frac1N\sum_r\sqrt{k_r}=\frac{79+48\sqrt2+42\sqrt3+3\sqrt6}{31104},\qquad \alpha_\text{cand}^{-1}=137.035\,999\,189\,566\,155\,029\,183\,744\,020\,866\,799\,566\,408\,684\,673\,854\,265\,541\ldots$$

*Proof.* Explicit $36\times36$ construction with numerical assertions on every step, then exact symbolic summation (`records_independent.py`, output in `records_independent.out`). The per-class counts also match the author's `records.csv`: 48+31 odd, 48 even mixed, 31+11 even same-end, 3 commutators. ∎

This confirms the counting and arithmetic. It does not test the rules that assign the weights.

## F-CT-2 The two exchanges are one involution

**Proposition 2.** Let $P$ be the span of the 169 product records. The "exchange of reads" $X$ and the "exchange of links" $Y$ act on $P$ through the same linear involution $E$. Write $\mathrm{cls}(r)\in\{\text{same-end},\text{mixed},\text{comm}\}$ and $e_r\in\{\pm1\}$ for the $E$-eigenvalue. The Letter's squared weight is the class function
$$k_r=\begin{cases}1 & e_r=-1,\ \mathrm{cls}\ne\text{comm}\\ 3 & e_r=+1,\ \text{same-end}\\ 2 & e_r=+1,\ \text{mixed}\\ 6 & \text{comm}.\end{cases}$$
In particular, the statement that $X$ and $Y$ "commute" (S7) holds trivially.

Machine-checked against the author's definitions: `ConventionCheck.weight_is_class_function` (Lean 4.35.0, `proofs/lean/chain-theory/`).

*Proof.* On slots, the read exchange is $Q_1\leftrightarrow Q_2,\ q_1\leftrightarrow q_2$. The link exchange carries each end of one read link to the same end of the other, which is the same permutation. The Supplement's own trace count writes $X=-E$ on same-end products and $Y=-E$ on mixed ones. In `ChainAlpha.lean`, `sigX` (via `Rec.reads`) and `sigY` both apply `exchAssign`. Substituting into $k_r=3^{[X_r=-1]}2^{[Y_r=-1]}$ gives the table. ∎

**Consequence.** For a center of dimension $d_Q$ and far ends of dimension $d_q$, with $a=d_Q^2-1$ and $b=d_q^2-1$, the numerator is
$$\Big[\tbinom a2+\tbinom b2+2ab\Big]+2ab\sqrt{d_q}+\Big[\tfrac{a(a+1)}2+\tfrac{b(b+1)}2\Big]\sqrt{d_Q}+\tbinom b2\sqrt{d_Qd_q}\,[d_q\text{ even}].$$
The content of the weight law is therefore two assignments:

- which dimension multiplies even same-end records (3, by the end rule);
- which multiplies even mixed records (2, by "copies are one qudit").

These are gap CT-G02.

## F-CT-3 Convention dependence of the selection

Take the centre $(ix\mathbb 1_{d_1}, iy\mathbb 1_{d_2})$ of $u(d_1)\oplus u(d_2)$. On $\mathbb C^{d_1}\otimes\mathbb C^{d_2}$ it acts as $i(x+y)$ if the far block acts by $\xi_2$, and as $i(x-y)$ if it acts by $\bar\xi_2$.

**Proposition 3.** Among readable types $\{(2,2),(2,3)\}$, the types whose centre acts nontrivially on the link are:

| Reading | Constraint on centre | Acting types | Hypercharge of (2,3) link |
|---|---|---|---|
| A (Letter) | $d_1x+d_2y=0$, far end by $\xi_2$ | $(2,3)$ only | $Y=1/6$ ($Q_L$) |
| A* | $d_1x+d_2y=0$, far end by $\bar\xi_2$ | $(2,2)$ and $(2,3)$ | $Y=5/6$ ($X$/$Y$ leptoquark block of the SU(5) adjoint) |
| B | $x+y=0$ (no phase on the matched directions) | none | — |
| C | centre mod kernel of the action | $(2,2)$ and $(2,3)$ | — |

*Proof.*
- A: $y=-d_1x/d_2$ gives $x+y=x(d_2-d_1)/d_2$, which vanishes iff $d_1=d_2$.
- A*: $x-y=x(d_1+d_2)/d_2\neq0$ for $x\neq0$.
- B: the action is identically zero.
- C: $(1,0)$ acts as $i$ on every link.

Symbolic check: `charge_conventions.py`. ∎

Both conventions give a 4-dimensional stabilizer of the pair ket at $(2,3)$, so CT-010 does not distinguish them. At $(2,2)$ the stabilizer has dim 4 under A and dim 3 under A*.

**Lean.** `proofs/lean/chain-theory/ConventionCheck.lean` is appended unchanged to the pinned upstream file. It proves `uniqueness_fails_under_conjugate_convention`: with `actingCharge` replaced by the A* value, the configuration "shared qubit, two read (2,2) links, one unread (2,2) link" satisfies every hypothesis of the upstream `configuration_unique` and differs from the object. Under reading A it is excluded (`upstream_excludes_qubitConfig`). A referee addition (`third_config_under_conjugate_convention`) shows that, because the upstream Lean predicate `readable` is $d\le3$ rather than the read cap, the configuration 3R 3R 3U at a shared qutrit also passes the Lean hypotheses under A*. The read cap excludes it (`config_enumeration.py`). Build log: `proofs/lean/chain-theory/run.out`.

**Full enumeration** (`config_enumeration.py`). The Supplement counts are reproduced independently: 76 configurations, 42 admitted by the read cap, 23 readable. Applying the end rule, the two-point definition and the acting-charge test, the selected configurations are:

| Reading | Selected configurations |
|---|---|
| A | the object only |
| A* | the object and centre 2: 2R 2R 2U |
| B | none |

Under A*, the competing configuration's readout under the Letter's own weight law is $(30+30\sqrt2)/2304$, so $1/\text{value}=31.81$ (`algebraic.out`, row (2,2)). The numerator agrees with the Supplement's own value for 2R 2R at a shared qubit. Theorem 1 would then deliver two two-point functions, and the Photon Conjecture would need a further rule to choose between them.

**Assessment.** Theorem 1 (Uniqueness) is correct under reading A. Reading A is an interpretive choice, not a consequence of the axiom. Its justification (S4, "sees no global phase of a pair" imposed over all $d_1+d_2$ directions of a direct sum) uses a space that the axiom does not supply. This is gap CT-G03/G04. Under A* the competitor is 2R 2R 2U at a shared qubit, $1/\text{value}=31.81$, shown below. (The Supplement's 2R 2R value, 3.535, has no unread link and fails the two-point definition under every convention.)

## F-CT-4 Freedom (look-elsewhere) analysis

**A. Published alternatives.** Table S5 has the computed row plus 15 closed-form alternatives. The computed row and the 12 weight-law alternatives were recomputed from the record classes $(79\text{ odd}, 48\text{ even mixed}, 36+6\text{ even same-end}, 3\text{ comm})$ and match the printed inverses to $<6\times10^{-4}$ (`freedom.py`, `freedom.out`). The other three rows are outside our declared space and far from 137: read-end counting (70.190), normalisation by the 172 records (0.758; an alternative to $N$) and the smallest readable configuration (3.535). The table omits these further alternatives to $N$:

| $N$ | 1/value |
|---|---|
| $24\cdot36^2$ (Letter) | 137.036 |
| $36^3$ (whole space at all three links) | 205.55 |
| $24^3$ (traceless at all three) | 60.90 |
| $24\cdot35^2$ | 129.53 |
| $35^3$ | 188.90 |

**B. A declared rule space.** We vary the following independently:

- the squared factor $k\in\{1,2,3,6\}$ for each of the four product classes (odd; even mixed; even $QQ$; even $qq$);
- the commutator factor ($\{1,2,3,4,6,9\}$) and count ($\{0,3,6\}$);
- amplitude vs Born weighting;
- the five values of $N$ above.

This gives $4^4\cdot6\cdot3\cdot2\cdot5=46\,080$ rules and 15 639 distinct values.

- Closest to CODATA 2022: the Letter's rule, at $0.598\sigma$. The next closest is $1.3\times10^5\sigma$ away.
- Distinct values within relative $10^{-2},10^{-3},10^{-4},10^{-6}$ of $\alpha^{-1}$: 288, 28, 4, 1.
- Density of distinct values with $1/\text{value}\in(120,155)$: 103.9 per unit. The expected number of chance values in the $\pm1\sigma$ window ($4.2\times10^{-8}$ wide) is therefore $4.4\times10^{-6}$.
- An $O(1)$ chance hit at this density needs about $10^{10}$ rules.

**C. Lattice density.** Consider $\{(a+b\sqrt2+c\sqrt3+d\sqrt6)/31104 : 0\le a,b,c,d\le M\}$.

- $M=100$: $1.04\times10^8$ elements, 1 hit within 1σ, expected 0.018.
- $M=200$: $1.63\times10^9$ elements, 1 hit within 1σ, expected 0.024.

The single hit is $(79,48,42,3)$ in both cases. At this denominator, 10-digit agreement is not a generic property of small-coefficient elements of the field.

**D. Lean-side cross-check.** We ran the author's `mutations.py` with Lean 4.35.0 on the pinned file. It applies 19 single-definition rewordings and has Lean evaluate each readout (`upstream_mutations_lean.out`; values in `mutations_values.out`).

- 13 of the 19 rewordings lie in the declared space B, and Lean reproduces our numerator exactly for each. 10 are Table S5 rows reproduced in Part A. Three are in B only: normalisation $36^3$ (205.55), reads at the qubit (150.79) and commutator records carrying no reads (138.94).
- Four rewordings lie outside B: keeping both commutator combinations (133.97), a shifted slot exchange (254.54), dropping the exchange-fixed records (149.59) and dropping the even records (295.10).
- (Corrected after review 2026-10-30, CT-2; see `reports/CORRECTIONS.md`.)
- Two rewordings leave the value unchanged: `fixesFar_always_true` (which breaks the end-rule statement but not the readout) and an equivalent definition of `wX`.
- No other rewording lands within $10^{-3}$ of $\alpha^{-1}$.

**Interpretation.** These are statements about trial factors, not about physics. Even allowing selection of the best rule within the declared space B, which contains the 12 weight-law alternatives the author lists, a chance agreement this close has probability $\approx4\times10^{-6}$; a single rule fixed in advance would have about $10^{-10}$. If the construction was selected from a search of $\gtrsim10^{9}$–$10^{10}$ implicit variants (wording clauses × configurations × normalizations × conventions such as F-CT-3), the agreement is unremarkable. The paper does not record the order in which rules and value were fixed (CT-G10), so the data in hand cannot distinguish these cases. Further digits can (F-CT-6).

**Limits.** The rule space B is ours. A different "natural" space would give a different density. The density estimate assumes local uniformity of rule outputs near 137, which holds approximately here (median gap 0.0070 between adjacent values).

## F-CT-5 Lean audit

- **Build.** `lean ChainAlpha.lean`, Lean 4.35.0 (commit `eb6cafb7`, release tarball sha256 `7289597a…268f6d`). It succeeds in 170 s on 4 cores. There is no `sorry`, `native_decide` or `implemented_by`.
- **Axioms reported.** `main`: `propext`. `records_orthogonal`: `propext, Quot.sound`. `selection_under_read_cap` and `configuration_unique`: `propext, Classical.choice, Quot.sound`. This matches S Table S-lean.
- **What `main` asserts.** A list equality, `readout objectWorld objectTwoPoint = ([(1,79),(2,48),(3,42),(6,3)], 31104)`: the radicand/coefficient list and the readout dimension. It is not a statement about real numbers.
- **What is a definition, not a theorem:**
  - `actingCharge d1 d2 := d1 ≠ d2`, so reading A is built in (F-CT-3);
  - `ReadCap` is a hypothesis;
  - `wX := readEnd` and `wY` encode the end assignment of CT-G02;
  - the $\sqrt{k}$ rule (CT-G08);
  - `readoutDim` (CT-G05).

The Lean file therefore certifies, as the paper states, that "the stated definitions produce Eq. (6)". It does not certify that the definitions follow from the axiom.

## F-CT-6 Numerical comparison and algebraic data

$\alpha_\text{cand}^{-1}-137.035\,999\,177=1.26\times10^{-8}$, which is $0.60\sigma$ (CODATA 2022, $\sigma=2.1\times10^{-8}$). The other deviations are $+1.57\sigma$ ($g-2$, Fan 2023), $-1.49\sigma$ (Rb, Morel 2020) and $+5.32\sigma$ (Cs, Parker 2018). The Rb and $g-2$ values differ from each other by $2.15\sigma$.

$\alpha_\text{cand}$ is a totally real algebraic number of degree 4. The minimal polynomial of $1/\alpha_\text{cand}$ is
$$62960159\,t^4-27464956416\,t^3-16967292715008\,t^2+9509021128065024\,t-935976560656121856 .$$
Its Galois conjugates are $465.71$, $406.51$ and $-573.03$ (`algebraic.py`).

**Falsification criterion (CT-020).** Suppose a determination of $\alpha(0)$ with standard uncertainty $u$, consistent with other determinations, finds $|\alpha^{-1}_\text{meas}-137.035\,999\,189\,566|>5u$. Then the Photon Conjecture is refuted as stated. At the current CODATA central value, this would take $u\lesssim2.5\times10^{-9}$, an 8× improvement. Settling the Rb/$g-2$ tension is the nearer-term test.

## F-CT-7 The same construction at other dimensions

$f(d_Q,d_q)$ is the Letter's construction with a center of dimension $d_Q$ and far ends of dimension $d_q$, with commutator records present only for even $d_q$. Values of $1/f$:

| $(d_Q,d_q)$ | (2,2) | (2,3) | (2,4) | (3,2) | (3,3) | (4,2) |
|---|---|---|---|---|---|---|
| $1/f$ | 31.81 | 140.40 | 216.04 | **137.036** | 791.66 | 314.67 |

(Full table in `algebraic.out`.) Most of these pairs are excluded by the read cap or the end rule. They are listed only to show that $1/f$ varies by orders of magnitude across neighbouring constructions, so the agreement is specific to $(3,2)$.

## F-CT-9 Cross-check of the imported read-cap lemmas

The read cap rests on published theorems (Størmer, Woronowicz, Choi), which are imported and not re-proved. As a numerical cross-check (`readcap_lemmas.py`) we take Choi's map, written as the generalized Choi map $\Phi[2,0,1]$.

- **Positivity.** Its Choi matrix is block-positive on $2\times10^5$ random product vectors.
- **Non-decomposability.** The SDP $\lambda^*=\max\{\lambda: C-\lambda\mathbb 1=P+Q^\Gamma,\ P,Q\succeq0\}$ gives $\lambda^*=-0.154701$; the value $2/\sqrt3-1=0.154701$ is observed numerically, not proven. So the map is not decomposable.
- **Embedding.** The same $\lambda^*$ holds after embedding into $M_3\to M_4$ and $M_4\to M_4$, which is the embedding step behind "$d_1d_2\ge 8$ ⇒ dark" for those types.
- **Control.** The transpose map gives $\lambda^*=0$ (decomposable), as it should.

The "every positive map is decomposable at $(2,2),(2,3)$" direction is a universal statement and is not checkable this way; it remains imported.

## F-CT-10 Covariance and the freedom left by symmetry

$G_\Omega$ is the symmetry group of the joint read state $\lvert\Omega\rangle=\lvert\omega_3\rangle_{Q_1Q_2}\lvert\omega_2\rangle_{q_1q_2}$: $U\otimes\bar U$ on $(Q_1,Q_2)$ and $V\otimes\bar V$ on $(q_1,q_2)$, with $U\in SU(3)$, $V\in SU(2)$. Results from `symmetry.py`:

1. **The plain exchange is not covariant.** The plain slot swap $E$ used by the Letter (and its Lean file) does not commute with $G_\Omega$; the relative defect is $1.4$. The exchange composed with transposition at the slots, $E'$, does commute (defect $5\times10^{-16}$). The reason is that $E$ identifies $Q_1$ with $Q_2$ (and $q_1$ with $q_2$) by the identity map, whereas the join and the qubit pair identify them by the transpose (CT-002 transpose identity).
2. **The Letter's records are not covariant.** The span of the 169 products is invariant. Its splitting into $E$-even and $E$-odd parts is not (leakage $0.89$), and neither is the span of the three commutator records (leakage $0.97$). So the Letter's record set depends on a choice of bases at the two read links, which the axiom says nothing lines up.
3. **The covariant version gives the same number.** Records built with $E'$ span a $G_\Omega\times\langle E'\rangle$-invariant 172-dimensional space. The commutator-type records become $D_x$, $D_z$ and $S_y=\sigma_y^{(q_1)}+\sigma_y^{(q_2)}$. The sector sizes are again $79/48/42/3$, so the value is unchanged. The construction can be repaired with no change to the number.
4. **Freedom left by symmetry.** Under $SU(3)\times SU(2)\times\langle E'\rangle$ the covariant record space decomposes into irreducible summands (SU(3), SU(2); $E'$-parity; multiplicity):

   | Parity | Summands |
   |---|---|
   | $+$ | $(1,1)^{\times2}$, $(1,5)$, $(8,1)$, $(8,3)^{\times2}$, $(27,1)$ |
   | $-$ | $(1,3)^{\times2}$, $(8,1)$, $(8,3)^{\times2}$, $(10,1)$, $(\overline{10},1)$ |

   That is 10 isotypic components over $\mathbb C$, and the commutant has real dimension 22. Any invariant dressing $M$ has $\operatorname{tr}M=\sum_{i=1}^{10}\dim\rho_i\,\operatorname{tr}A_i$. A dressing that maps Hermitian records to Hermitian records commutes with $A\mapsto A^\dagger$, which exchanges the $10$ and $\overline{10}$ components, so those two weights are equal. The result is a **9-parameter (real)** family of readouts; the count was corrected after review CT-5. The Letter's rule is one point in it, fixed by the wording clauses of CT-G02. In particular, within the single isotypic component $(1,3)_-$ the Letter gives weight 1 to the three antisymmetric qubit–qubit products and $\sqrt6$ to the three commutator records. Symmetry does not distinguish these two copies of the same representation; only the "one letter vs two letters" construction rule does.

## F-CT-8 Exploration: do natural state functionals give the weights? (negative)

Labelled exploration, toward CT-G01. The joint read state is exactly $\lvert\Omega\rangle=\lvert\omega_3\rangle_{Q_1Q_2}\otimes\lvert\omega_2\rangle_{q_1q_2}$ (verified: the projected vector has squared norm $1/6$, matching S10). For HS-normalized record operators $O_r$ we computed $|\langle\Omega|O_r|\Omega\rangle|$ and $36\lVert O_r\Omega\rVert^2$ (`weight_functionals.py`).

- Neither is constant on the weight classes. Odd records alone take 7 distinct values of $36\lVert O_r\Omega\rVert^2$, from 0 to 3.
- The class means are $0.848$ (odd), $1.083$ (even mixed), $1.190$ (even same-end) and $0.667$ (comm). They are not in the ratio $1:\sqrt2:\sqrt3:\sqrt6$.
- The vacuum expectation vanishes on every class except even same-end.

So the weight law is not the vacuum amplitude or norm of the record operators. This does not exclude other operator constructions; CT-G01 remains open.

## Graded ledger (α Letter)

| Claim | Mathematical validity | Physical validity |
|---|---|---|
| CT-014, 015, 017, 018 (counts, value) | Proven (independently) | n/a |
| CT-005 (read cap) | Proven from imported theorems + definition of "fix" | n/a |
| CT-007, 012 (acting charge, uniqueness) | Proven under reading A; fails under A*, B, C | n/a |
| CT-016 (weight law) | Definition; three wording clauses; alternatives give 112–147 | Untested |
| CT-019 (on-shell $Z_3$) | Analogy, not derivation | Open |
| CT-020 (Photon Conjecture) | Conjecture | Agrees at 0.6σ; decisive test needs $u\lesssim 2.5\times10^{-9}$; hypercharge vs EM gap (CT-G06) open |
