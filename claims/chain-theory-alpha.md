# Claim Register — Chain Theory, α Letter (CT-S1, CT-S2, CT-S3)

The compendium (CT-S5) has its own register: `claims/chain-theory-compendium.md`.
Location: L = Letter, S = Supplement. Type: Def / Ax / Thm / Pred / Ana / Num.
"Formalizable": whether the claim can be stated as precise mathematics *as written*.

| ID | Statement (minimal paraphrase) | Type | Formalizable? | Source |
|---|---|---|---|---|
| CT-001 | Axiom: nothing exists but finitely many maximally entangled pairs of qudits, some shared, known only by distilling them. | Ax | Partially. "Known", "distilling" need the five definitions; the axiom does not say whether a qudit's far end transforms by $U$ or $\bar U$ under basis changes (see CT-G03). | L §Model |
| CT-002 | Pair of type $(d_1,d_2)$ is $\lvert\omega_m\rangle=m^{-1/2}\sum_{i\le m}\lvert ii\rangle$, $m=\min(d_1,d_2)$. | Def | Yes | L Eq. 1 |
| CT-003 | Distilling = projection $P_m=\lvert\omega_m\rangle\langle\omega_m\rvert$, with no failure branch. | Def | Yes (as a definition). Physically it is postselection; S1 acknowledges this. | L, S Cor. S2 |
| CT-004 | A read at a prime-dimension end has $d+1$ MUB bases; their traceless projector parts span $su(d)$. | Thm (imported) | Yes. Standard (Wootters–Fields 1989). | L |
| CT-005 | **Read cap:** a type is readable (every positive map between its ends is decomposable) iff $d_1d_2\le 6$, so the readable types are $(2,2),(2,3)$. | Thm | Yes. It follows from imported theorems (Størmer, Woronowicz, Choi). The step "a read fixes two ends only when every positive map is decomposable" is a *definition of fix*. | L, S2 |
| CT-006 | The identification algebra of a $(d_1,d_2)$ link is $s(u(d_1)\oplus u(d_2))$, acting on the link as $\xi_1\otimes 1+1\otimes\xi_2$. | Def | Yes, but it is a choice (CT-G03, CT-G04). | L Eq. 3, S3 |
| CT-007 | Its centre acts on the link as $i(d_2-d_1)$; only $(2,3)$ among readable types has an acting central charge. | Thm | Yes. **Holds under the Letter's convention only.** See F-CT-3. | S Lemma S-center |
| CT-008 | $S(U(3)\times U(2))\cong(SU(3)\times SU(2)\times U(1))/\mathbb Z_6$, and $\mathcal Q=6iY$. | Thm (imported) | Yes. Standard. | L |
| CT-009 | The link carries $(3,2)_{1/6}$, the quantum numbers of $Q_L$. | Thm | Yes. Holds under the unconjugated convention. Under the conjugated one, $Y=5/6$, which is the SU(5) $X$/$Y$-boson block. See F-CT-3. | S3 |
| CT-010 | Stabilizer of the pair ket in the (2,3) algebra has dim 4, with centre $\mathcal U=i\,\mathrm{diag}(1,1;-1,-1,0)$. | Thm | Yes under A. dim 4 verified under both conventions; the centre $\mathcal U$ holds under A only (under A* it is $i\,\mathrm{diag}(1,1;1,1,-4)$ in labels; `conjugation.md` Prop. 6, 2026-10-30). | S Thm S-stabilizer |
| CT-011 | $\mathcal Q$ is conserved; $Z_1=Z_2$. | Thm | Partially. "Conserved count" is defined in-model; the analogy to the Ward identity is not a derivation of a QFT identity. | S4 |
| CT-012 | **Uniqueness (Thm 1):** exactly one configuration has readable read links, is a two-point function, is read at the qutrit end, and carries an acting charge on every link: a shared qutrit with three qubit links, two of them read. | Thm | Yes, given CT-005/007 and the definitions "two-point function" and "end rule". Lean checks it under these as hypotheses/definitions. | L Thm 1, S5 |
| CT-013 | Readout dimension $N=\dim R=24\cdot 36\cdot 36=31104$. | Def + count | Yes. Whole operator space at read links, traceless part at the unread link: a definitional choice. | L Eq. 4, S6 |
| CT-014 | 169 two-slot products; the read exchange fixes 11 and pairs 158; 3 commutator records; 172 records, mutually orthogonal. | Thm | Yes. **Reproduced independently** (F-CT-1). | L, S Cor. S-records |
| CT-015 | Sector sizes $(+,+)\!:\!79$, $(+,-)\!:\!48$, $(-,+)\!:\!42$, $(-,-)\!:\!3$ under exchanges $X$, $Y$. | Thm | Yes. **Reproduced.** $X$ and $Y$ are the same involution on product records (F-CT-2). | L Table I, S Lemma S-sectors |
| CT-016 | **Weight law:** $k_r=3^{[X_r=-1]}2^{[Y_r=-1]}$, weight $\sqrt{k_r}$. | Def (presented as derived) | Partially. The factors $\sqrt d$ come from "amplitude summed over matched directions". The assignment of $X$ to the qutrit and $Y$ to the qubit comes from three wording clauses (S Table S-inputs). | L Eq. 5, S Thm S-gain |
| CT-017 | **Value (Thm 2):** $G_\text{dressed}/G_\text{bare}=(79+48\sqrt2+42\sqrt3+3\sqrt6)/31104$. | Thm | Yes, as arithmetic of the definitions. **Reproduced exactly** (F-CT-1). | L Thm 2 |
| CT-018 | $\alpha_\text{cand}^{-1}=137.035\,999\,189\,566\,155\,029\,183\,744\ldots$ | Num | Yes. **Reproduced to 60 digits.** | L Eq. 6 |
| CT-019 | The normalized readout is the on-shell $Z_3$ with $\alpha_\text{bare}=1$. | Def / Ana | No (as physics). It is an analogy to on-shell QED; the model has no propagator, momentum, or spacetime (S8 says so). | L, S8 |
| CT-020 | **Photon Conjecture:** the (2,3) link is the elementary photon link and $\alpha_\text{cand}=\alpha(0)$. | Pred | Testable numerically. Physically underspecified: the charge is hypercharge, not electric charge (CT-G06). | L Conj. 1 |
| CT-021 | $\alpha_\text{cand}$ is $0.6\sigma$ from CODATA 2022, $+1.6\sigma$ from $g-2$, $-1.5\sigma$ from Rb, $+5.3\sigma$ from Cs. | Num | Yes. **Reproduced** (F-CT-6). | L Fig. 3 |
| CT-022 | Every alternative in Table S-alternatives changes the value; two of them (the ends swapped; commutator records surviving $X$ only) land within 1% of it, which the Supplement states. | Num | Yes. **All 15 closed-form alternatives reproduced** from their printed closed forms (`mutations_classified.out`; the 12 weight-law ones also from the record classes, F-CT-4). The smallest change is 0.83% in 1/value (0.84% in the value), for the ends swapped; the next is 0.96%; all others exceed 1%. (Corrected after the second-round review, N4: this entry said "≥0.9%" and "all 13".) | S9 |
| CT-023 | Lean's `main` theorem certifies that the stated definitions produce the closed form, using only `propext`. | Thm (meta) | Yes. Checked by building with Lean 4.35.0 (F-CT-5). | L §Formalization, S11 |
