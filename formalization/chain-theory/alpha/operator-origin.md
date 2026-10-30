# Is there an operator origin of the weight law? (gaps CT-G01, CT-G12)

**Object.** The dressed insertion $M$ of the α Letter (Wissner-Gross, v1.0.0, Eq. (5) and Supplement Thm S-gain). The Letter *defines* $M$ to carry $\sqrt{k_r}$ on record $r$, with $k_r=1$ (E′-odd products), $2$ (E′-even mixed), $3$ (E′-even same-end), $6$ (commutator records). The records here are the covariant ones (exchange $E'=E\circ T$, F-CT-10). The sector sizes and the value are those of the Letter.
**Question.** Is there an operator, built from the model's primitives and not from the record list, whose matrix, diagonal or per-isotypic traces on the record space $R$ reproduce $\sqrt{k_r}$?
**Code.** `computation/chain-theory/alpha/operator_origin.py` → `operator_origin.out`. It imports `records_independent.py` and `symmetry.py` and does not modify them. Every number below is printed there. "out §n" refers to section n of that output; "Section n" refers to this document.
**Date.** 2026-10-30. Exploratory analysis toward CT-G01; it does not change the audited value. Scope of Lemma 3 and Theorem 4, the frame of "local", and the description of $S\Pi_X$ revised after review REVIEW3-CT (R3-CT-1–3, R3-CT-7; `reports/CORRECTIONS.md`).

## Result in brief

1. **An explicit operator exists, but it does not derive the weights.** On $R$,
   $$M_L=\bigl(1+(a_3-1)\,S\,\Pi_X\bigr)\bigl(1+(a_2-1)\,S\,\Pi_Y\bigr),\qquad a_d=\langle\Phi_d|\omega_d\rangle=\sqrt d,$$
   where $S$, $\Pi_X$, $\Pi_Y$ are projectors built from partial traces and $E'$ alone (Prop. 2; numerical check $\lVert M_\text{nat}-M_L\rVert=4.6\times10^{-15}$). This meets the closure criterion that CT-G01 states. It does not close the gap. $S\Pi_X$ gives the qutrit amplitude $\sqrt3$ to 9 of the 12 records that have no qutrit letter: the 6 E′-even qubit–qubit records and the 3 commutator records. (The 3 E′-odd qubit–qubit records lie in $\Pi_X$, but $S=0$ there.) That assignment is the end-rule clause of CT-G02, written as a projector.
2. **No multiplicative construction reproduces the weights.** Family F1 has all 56 355 words of length ≤ 3 in a 38-letter alphabet of primitives, plus the identity. They give 11 973 distinct compressions to $R$. No word reproduces $\sqrt{k_r}$ or $k_r$ in any tested sense: as an operator, as its diagonal, as its per-isotypic traces, after a twirl, or under $|C|$ or $\sqrt C$. Two theorems extend parts of this negative to every word length.
   - *Galois* (Lemma 3): every primitive is rational over $\mathbb Q(i)$. This extends the $\sqrt{k_r}$, trace and value negatives to all lengths and all $\mathbb Q(i)$-linear combinations.
   - *Factorization* (Thm 4): every covariant primitive except $\Pi_\pm=(1\pm E')/2$ factorizes across qutrit|qubit slots. A product cannot give mixed records the ratio $\sqrt2:1$ (or $2:1$). This extends the $k_r$ negative to all products of factorizing letters and their twirls.
   - For words that contain one of the 8 non-factorizing letters ($\Pi_\pm$ and the six read-link letters), the $k_r$, $|C|$ and $\sqrt C$ negatives hold up to length 3 only. Functional calculus escapes Lemma 3 ($K_L\in\mathcal A$, Prop. 2), and Theorem 4 covers products only.
3. **No dressing local in the original frame reproduces the weights.** Family F2 is the operator form of the Letter's own mechanism: survival projector × end amplitude, for two exchanges. It has 321 489 ordered choices and 19 926 distinct dressings. It contains all nine survival/attribution rows of Supplement Table S5 (the rows that keep the 172 records and the amplitude level). $M_L$ is reached by 36 choices, and every one of them credits $\sqrt3$ to records with no qutrit letter. No choice that is *local in the original frame* reaches it. Such a choice credits an end's dimension only to records with a letter at that end, read where the letters sit. Prop. 5 proves this for every such product-of-dimensions dressing, not only for those in F2. Locality is frame-relative: after the qubit letters are transported to the qutrit copies through the read links, 6 frame-local choices reach $M_L$, with the frame assignment and the full-join credit as inputs (`transport.md`, Theorem 6, Prop. 7).

**Precise negative statement.** No operator in F1 (11 973 distinct compressions out of 56 355 words) is class-diagonal on $R$ with eigenvalue ratios $1:\sqrt2:\sqrt3:\sqrt6$, whether or not it is $G_\Omega$-invariant, and also after the $G_\Omega\times\langle E'\rangle$-twirl, $|C|$ or $\sqrt C$. None is so with ratios $1:2:3:6$ (the square-root route). None has per-isotypic traces proportional to the Letter's point. In F2, $M_L$ is attained only by attributions that are not local in the original frame (36 of 321 489 ordered choices).

**Consequence.** The weight law has an operator *form* (survival projector times $\langle\Phi_d|\omega_d\rangle$) but no operator *origin* in these families. CT-G01 reduces to CT-G02, the end at which each exchange acts, and to CT-G08, the amplitude level.

---

## 1. The primitive set

$H=\mathbb C^3_{Q_1}\otimes\mathbb C^2_{q_1}\otimes\mathbb C^3_{Q_2}\otimes\mathbb C^2_{q_2}$. Superoperators act on $B(H)$ (dim 1296); $R\subset B(H)$ is the 172-dim covariant record span. Each primitive is taken from the Letter:

| Primitive | Letter source | Superoperator(s) |
|---|---|---|
| Join $P_3=\lvert\omega_3\rangle\langle\omega_3\rvert_{Q_1Q_2}$ | S7 joint state | $L_{P_3},R_{P_3},\mathrm{Ad}_{P_3}$ |
| Qubit pair $P_2$ on $(q_1,q_2)$ | S10 joint state | $L_{P_2},R_{P_2},\mathrm{Ad}_{P_2}$ |
| Joint state $\lvert\Omega\rangle\langle\Omega\rvert=P_3P_2$ | S7 | $\mathrm{Ad}_{P_\Omega}$ |
| Read links $P_{L_i}=\lvert\omega_2\rangle\langle\omega_2\rvert_{Q_iq_i}$ (support levels 0,1 of $Q_i$) | Eq. (1), pair of type (2,3) | $L,R,\mathrm{Ad}$ of $P_{L_1},P_{L_2}$ |
| Exchange $E$ (plain swap), covariant $E'=T\circ E$ | S7; F-CT-10 | $E$, $E'$ |
| Exchange at one end type: $E'_Q$ (swap $Q_1\leftrightarrow Q_2$ + transpose there), $E'_q$; plain $F_Q,F_q$ | "exchange of reads/links" | $E'=E'_QE'_q$ |
| Transpose identity | Eq. (1)–(2) | $T$, partial $T_s$ |
| Partial trace at a slot | readout algebra | $\tau_s(A)=\tfrac{1}{d_s}\mathbb 1_s\otimes\mathrm{Tr}_sA$, $1-\tau_s$ |
| Exchange eigenprojectors | "sorts that space into its two eigenspaces" | $\Pi_\pm=\tfrac12(1\pm E')$ |
| Own sign | S7 "its own sign" | $\sigma=(-1)^{\hat N}$, $\hat N=\sum_s(1-\tau_s)$ |
| Twirls | "invariant under ... the isotropic twirl" | $\mathcal T_\Omega$ (onto $\mathrm{span}\{1,P_3\}\otimes\mathrm{span}\{1,P_2\}$), $\mathcal T_Q$, $\mathcal T_q$ |

**The "amplitude summed over matched directions".** Write $\Phi_d=\sum_{i=1}^d|ii\rangle$ for the flat vector and $\omega_d=\Phi_d/\sqrt d$. The Letter's $d\cdot(1/\sqrt d)$ is
$$a_d:=\langle\Phi_d|\omega_d\rangle=\operatorname{Tr}C_{\omega_d}=\operatorname{Tr}\sqrt{\rho_d}=\sqrt d,$$
with $C$ the coefficient matrix and $\rho_d=\mathbb 1/d$ the reduced state (out §1). The Supplement's "traced against the projector onto the orbit of the letter" gives the same number, $\operatorname{Tr}(\mathbb 1_d\cdot\mathbb 1_d/\sqrt d)$. As an operator, the amplitude is $\lvert\Phi_d\rangle\langle\omega_d\rvert=\sqrt d\,P_d$ (out §1). It is a scalar times a pair projector, and it carries no information about which record it attaches to.

**Slot classes are derivable.** The occupation-type projectors $[\hat N_Q=a,\hat N_q=b]=\sum_{|o\cap Q|=a,\,|o\cap q|=b}\prod_{s\in o}(1-\tau_s)\prod_{s\notin o}\tau_s$ are polynomials in the $\tau_s$. On $R$ there are four types: QQ (64 records), qq (9), mix (96) and comm (3) (out §0). The Letter's four classes are joint eigenspaces of $(E',\hat N_Q,\hat N_q)$ on $R$, so the class projectors $\Pi_c$ are natural.

**Properties (out §1, computed for each of the 38 letters).**
- 24 letters are $G_\Omega$-covariant (defect $<10^{-15}$); the other 14 are not (defect $\approx1$).
- 30 letters are Q|q *products*: operator-Schmidt rank 1 across $B(Q_1Q_2)\otimes B(q_1q_2)$. The exceptions are $\Pi_\pm$ (rank 2) and the read-link multiplications (rank 4 or 16).
- Every letter has entries in $2^{-a}3^{-b}\,\mathbb Z[i]$; the largest denominator is 144, at $\mathcal T_\Omega$.

## 2. Exact results

**Proposition 1 (symmetry structure, out §2).** $R$ splits into 10 isotypic components under $SU(3)\times SU(2)\times\langle E'\rangle$, with 10 and $\overline{10}$ separated, and the commutant has dimension 22, as in F-CT-10. $M_L$ lies in the commutant (defect $6.7\times10^{-14}$). Each component's class content and the exact per-isotypic traces of $M_L$ (the "Letter's point") are:

| component | mult | class content (odd, evmix, evsame, comm) | $\operatorname{tr}M_L$ |
|---|---|---|---|
| $(1,1)_+$ | 2 | (0,0,2,0) | $2\sqrt3$ |
| $(1,3)_-$ | 2 | (3,0,0,3) | $3+3\sqrt6$ |
| $(1,5)_+$ | 1 | (0,0,5,0) | $5\sqrt3$ |
| $(8,1)_\mp$ | 1, 1 | (8,0,0,0); (0,0,8,0) | $8$; $8\sqrt3$ |
| $(8,3)_\mp$ | 2, 2 | (48,0,0,0); (0,48,0,0) | $48$; $48\sqrt2$ |
| $(10,1)_-$, $(\overline{10},1)_-$ | 1, 1 | (10,0,0,0) each | $10$ each |
| $(27,1)_+$ | 1 | (0,0,27,0) | $27\sqrt3$ |

The two copies of $(1,3)_-$, which CT-G12 notes symmetry cannot separate (weights 1 and $\sqrt6$), are separated by $\hat N$: $\hat N=2$ on the q₁q₂-odd copy and $\hat N=1$ on the commutator copy. $\hat N$ is $G_\Omega\times\langle E'\rangle$-covariant and is built from partial traces. So the distinction CT-G12 refers to is available as an operator. The weights it receives are not.

**Proposition 2 (the Letter's $M$ as a polynomial in natural projectors).** Let
$$S=\tfrac12\bigl(1+\sigma E'\bigr),\qquad \Pi_X=[\hat N_Q\hat N_q=0],\qquad \Pi_Y=[\hat N_q=1].$$
These are commuting projectors on $B(H)$; each maps $R$ into $R$ and is diagonal on the records (out §3). Then on $R$
$$M_L=\bigl(1+(\sqrt3-1)S\Pi_X\bigr)\bigl(1+(\sqrt2-1)S\Pi_Y\bigr),\qquad K_L=M_L^2=(1+2S\Pi_X)(1+S\Pi_Y)=\tfrac12(2-\hat X)(3-\hat Y),$$
with the natural involutions $\hat X=1-2S\Pi_X$ and $\hat Y=1-2S\Pi_Y$. This is the Letter's Corollary S-gram, $k=\tfrac12(2-X)(3-Y)$, now with $X$ and $Y$ given as operators.

*Proof.* $S$ is the Letter's survival. $\sigma=+1$ on products (two occupied slots) and $-1$ on commutator records (one occupied slot per term), so $S$ projects onto "$E'$ returns $A$ times its own sign". The table:

| class | $S$ | $\Pi_X$ | $\Pi_Y$ | $M$ |
|---|---|---|---|---|
| E′-odd product | 0 | · | · | 1 |
| E′-even mixed ($\hat N_Q=\hat N_q=1$) | 1 | 0 | 1 | $\sqrt2$ |
| E′-even QQ or qq | 1 | 1 | 0 | $\sqrt3$ |
| commutator ($\hat N_Q=0,\hat N_q=1$) | 1 | 1 | 1 | $\sqrt6$ |

∎ (Numerical: leakage $2.4\times10^{-16}$, $\lVert M_\text{nat}-M_L\rVert=4.6\times10^{-15}$.)

*What this does and does not show.* The *form* of the weight law is operator-natural: the survival projector, the amplitude $a_d=\langle\Phi_d|\omega_d\rangle$, and Cor. S-gram. The *content* lies in $\Pi_X$ and $\Pi_Y$. $\Pi_X=[\hat N_Q\hat N_q=0]$ ("all letters at one end type") carries the qutrit amplitude $a_3$, and it contains the 12 records whose letters all sit at the qubits. A rule that credits each end only where the record has a letter there (Prop. 5) must instead give the qubit-only records factors of $a_2$. Table S5's "X where letters sit" is such a rule: $k=2$ on E′-even q₁q₂ records and $k=4$ on commutator records, value 139.030 (out §5). The choice between these rules is CT-G02.

**Lemma 3 (Galois grading).** Let $\mathbb F=\mathbb Q(i)$ and let $\mathcal A$ be the $\mathbb F$-algebra generated by the 38 letters (all $\mathbb F$-rational, out §1). Each covariant record is a positive multiple of an $\mathbb F$-rational matrix; the factors are 1, $1/\sqrt3$ or $1/3$, from $\lambda_8$. For every $\Phi\in\mathcal A$:
1. the Rayleigh quotients $\langle r,\Phi r\rangle/\langle r,r\rangle$ lie in $\mathbb F$;
2. the traces $\operatorname{tr}(\Pi_i\Pi_R\Phi\Pi_R)$ on the components $(8,1)_\pm$ lie in $\mathbb F$, because $\Pi_R$ and the spectral projectors of the rational Casimirs and of $E'$ are $\mathbb F$-rational.

Since $\{1,\sqrt2,\sqrt3,\sqrt6\}$ is a basis of $\mathbb Q(i,\sqrt2,\sqrt3)=\mathbb Q(\zeta_{24})$ over $\mathbb F$, no $\Phi\in\mathcal A$ has diagonal $\propto\sqrt{k_r}$ ($\sqrt2\notin\mathbb F$). None has per-isotypic traces $\propto$ the Letter's point ($8\sqrt3/8\notin\mathbb F$), and none has $\operatorname{tr}_R\Phi/N=\alpha_\text{cand}$. All quantities used are basis-free (Rayleigh quotients, traces on $\mathbb F$-rational invariant subspaces, $\operatorname{tr}_R$); in the orthonormal record basis the off-diagonal entries need not lie in $\mathbb F$. The lemma says nothing about the integer pattern $k_r$.

*Corollary (where irrationals can enter).* Adjoining the amplitudes $a_2,a_3$ (normalized versus flat pair vectors), every construction is uniquely $\Phi=\sum_{a,b\in\{0,1\}}a_2^a a_3^b\,\Phi_{ab}$ with $\Phi_{ab}\in\mathcal A$. Its compression equals $M_L$ iff
$$\Phi_{00}|_R=\Pi_\text{odd},\quad \Phi_{10}|_R=\Pi_\text{evmix},\quad \Phi_{01}|_R=\Pi_\text{evsame},\quad \Phi_{11}|_R=\Pi_\text{comm},$$
as identities of superoperators on $R\subset B(H)$, not entrywise in the orthonormal record basis. The only other source of $\sqrt{\cdot}$ is functional calculus, $M=\sqrt K$ with $K\in\mathcal A$, and then $K=K_L$ (eigenvalues 1, 2, 3, 6) is the same assignment problem over the integers. Lemma 3 does not exclude this route ($K_L\in\mathcal A$ by Prop. 2); Theorem 4 excludes it for products and their twirls. In short, an operator origin of the weights is equivalent to a derivation of *which class receives which end's amplitude*.

**Theorem 4 (factorization obstruction).** Let $W_Q\subset B(Q_1Q_2)$ be the 16-dim span of single-qutrit-letter operators and $W_q\subset B(q_1q_2)$ the 6-dim analogue. The 96 mixed records span $W_Q\otimes W_q$, on which $E'=E'_Q\otimes E'_q$. Let $\Phi$ be a Q|q-product superoperator, or its $G_\Omega$-twirl, or its $G_\Omega\times\langle E'\rangle$-twirl. If the compression of $\Phi$ to the mixed records is $\mu_+\Pi_+^\text{mix}+\mu_-\Pi_-^\text{mix}$, then $\mu_+=\pm\mu_-$; if $\Phi\ge0$ there, $\mu_+=\mu_-$. Hence no such $\Phi$, nor $\sqrt\Phi$, equals $cM_L$ (which needs $\sqrt2:1$) or $cK_L$ (which needs $2:1$) for any $c\in\mathbb C$.

*Proof.*
1. *Product case.* The compression of $\Phi_Q\otimes\Phi_q$ to $W_Q\otimes W_q$ is $A\otimes B$. Write $\mu_+\Pi_++\mu_-\Pi_-=\tfrac12(\mu_++\mu_-)\,1\otimes1+\tfrac12(\mu_+-\mu_-)\,E'_Q\otimes E'_q$. This has operator-Schmidt rank 2 unless one coefficient vanishes, since $E'_Q|_{W_Q}\not\propto1$.
2. *$G_\Omega$-twirl.* $G_\Omega=SU(3)\times SU(2)$ acts as a product, so its twirl preserves products.
3. *$G_\Omega\times\langle E'\rangle$-twirl.* Use $W_Q\cong 8\otimes\mathbb C^2$ (slot $Q_2$ identified with $Q_1$ by $\lambda\mapsto\lambda^T$) and $W_q\cong3\otimes\mathbb C^2$. After the $G_\Omega$-twirl, $A=1_8\otimes a$ and $B=1_3\otimes b$ with $a,b\in M_2$, and $E'_Q=1\otimes s$, where $s$ is the swap. The $E'$-average is
   $$\tfrac12(a\otimes b+sas\otimes sbs)=a_0\otimes b_0+a_1\otimes b_1,$$
   where $a_0=\tfrac12(a+sas)\in\mathrm{span}\{1,s\}$ and $a_1$ anticommutes with $s$. $M_2\otimes M_2$ is the direct sum of the four sectors (commuting/anticommuting with $s$, on each factor). The target $\alpha\,1\otimes1+\beta\,s\otimes s$ lies in the commuting⊗commuting sector, where $a_0\otimes b_0$ has rank 1. So $\alpha\beta=0$.

∎

*Machine check.* Both $M_L$ and $K_L$ have mixed-block Schmidt rank 2 (out §3). Among the 222 factorizing F1 words whose twirled mixed block is $E'$-scalar, none violates $|\mu_+|=|\mu_-|$ (out §4). The non-factorizing words with an $E'$-scalar twirled mixed block (they contain $\Pi_\pm$ or a read link) give $\mu_+/\mu_-\in\{\infty,0,\pm1\}$, except 72 words with ratios $\pm\tfrac12,\pm\tfrac35,\pm\tfrac53,\pm2$. The printed examples combine $\Pi_\pm$ with partial transposes and partial traces. None of the 72 is class-scalar after the twirl (out §4). Nine of them have the ratio $+2$ of $K_L$, which Theorem 4 forbids for products; they fail on the other classes, and nothing proved here excludes a longer non-product word with compression $cK_L$.

So the $\sqrt2$ that the Letter gives mixed records (its $Y$) cannot come from any product of the covariant model operators. It requires an affine combination such as $\Pi_+=(1+E')/2$, that is, a dressing $1+(a_2-1)P$, or the read-link projectors, which break $G_\Omega$.

**Proposition 5 (locality in the original frame).** Call a dressing *local* (in the original frame; below, "local" always means this) if its squared weight on each record is a product of end dimensions $d_e$, each credited only where the record has a letter at an $e$-slot. The 12 records with letters only at the qubits are the 9 q₁q₂ products and the 3 commutator records. On these a local dressing has $k\in\{1,d_q,d_q^2,\dots\}=\{1,2,4,\dots\}$. The Letter needs $k=3$ on the 6 E′-even q₁q₂ records and $k=6$ on the commutator records, so **no dressing local in the original frame reproduces $M_L$.** The content is that 3 and 6 are not powers of $d_q=2$. ∎ Locality is relative to where the letters are read: after transporting the qubit letters to the qutrit copies, a frame-local dressing does reproduce $M_L$, with the frame of each exchange and the full-join credit as inputs (`transport.md`, Theorem 6, Prop. 7). The non-local factor is the end-rule clause: "a read is an act at the end at which it is taken ... $X$ contributes 3, at the qubit–qubit records too" (Letter; Supplement Thm S-gain (ii), and the strand picture's own remark that "the picture ... does not supply" it).

## 3. Family F1: words in the primitives

**Declared family.** The 38 letters of Section 1 (24 covariant, 14 not), acting on $B(H)$. The words are all products $w=\pi_1\cdots\pi_\ell$ with $\ell\le3$, plus the identity (the bare insertion): $1+38+38^2+38^3=56\,355$ words. Words whose images of $R$ agree up to a scalar are merged; the shortest is kept, and its left extensions cover the rest. This leaves 11 973 distinct nonzero images of $R$ (out §4). Each is compressed to $C=\Pi_R w\Pi_R$ in the orthonormal covariant record basis.

**Tests.**
- (a) $G_\Omega$- and $E'$-invariance of $C$.
- (b) Is $C$ class-scalar, $C=\sum_c v_c\Pi_c$, with $v\propto(1,\sqrt2,\sqrt3,\sqrt6)$ or $\propto(1,2,3,6)$? (b′) Is the diagonal constant on classes? The fit residual of $\operatorname{diag}C$ against $c\sqrt{k_r}$ and $ck_r$ is also recorded.
- (c) Are the per-isotypic traces (10-vector) proportional to the Letter's point?
- Each test is repeated for the $G_\Omega\times\langle E'\rangle$-twirl of $C$, for $|C|=(C^\dagger C)^{1/2}$, and for $\sqrt C$ when $C\ge0$.

**Results (out §4).**
- **Size.** 56 355 words. 29 679 are formed explicitly; the rest are left extensions of images already seen. They give 12 021 distinct images of $R$; 48 compress to zero on $R$, leaving **11 973 analysed compressions**.
- **(a) Invariance.** 1629 compressions are $G_\Omega$-invariant (1575 from all-covariant words) and 295 are also $E'$-invariant. 950 words map $R$ into $R$, and 148 compressions are diagonal in the record basis.
- **(b) Class-scalar operators.**

  | form | class-scalar words | distinct patterns |
  |---|---|---|
  | raw $C$ | 19 | 8 |
  | $G_\Omega\times\langle E'\rangle$-twirled | 155 | 9 |
  | $\lvert C\rvert$ | 67 | 2 |
  | $\sqrt C$, $C\ge0$ (636 such words) | 1 | 1, namely $(1,1,1,1)$ |

  **Every one of these 242 patterns has entries in $\{0,\pm1\}$ up to scale.** They are indicator and sign operators: $1$, $E'$, $\sigma$, $\sigma E'$, $\Pi_\pm$, $\sigma\Pi_-$, and projectors isolating the commutator class or the even-mixed class. **None is $\propto(1,\sqrt2,\sqrt3,\sqrt6)$, as Lemma 3 requires, and none is $\propto(1,2,3,6)$, so the square-root route fails too.**
- **(b′) Class-constant diagonals.** 138 words have a diagonal constant on each class, with 10 patterns, all in $\{0,\pm1\}$.
- **Diagonal fits.** No word fits $\operatorname{diag}C\approx c\sqrt{k_r}$ better than the identity, the bare insertion (residual 0.2471 for both). The same holds for $|C|$ and $\sqrt C$. For $c\,k_r$ the best is $L_{P_{L1}}\Pi_+L_{P_{L1}}$, at 0.4572 against 0.4674 for the identity.
- **(c) Per-isotypic traces.** The projective residual against the Letter's point is 0.2070 for the identity.

  | form | best word | residual |
  |---|---|---|
  | raw $C$ | $L_{P_{L2}}\,\Pi_+\,L_{P_{L2}}$ (and its three read-link images) | 0.1008 |
  | $\lvert C\rvert$ | $\Pi_+\,R_{P_{L2}}\,T_{q_1}$ | 0.0912 |
  | $\sqrt C$ | $L_{P_{L2}}R_{P_{L1}}$ | 0.1007 |

  **Among $G_\Omega$-invariant words nothing beats the identity** (0.2070). So every closer miss has a non-invariant compression and contains a non-covariant letter. The listed ones use read-link projectors or partial transposes.
- **What distinguishes the best miss.** $L_{P_{L2}}\Pi_+L_{P_{L2}}$ has class means $(1,1.40,1.37,0.97)$ relative to the odd class. The even-mixed and even-same-end classes come out about equal, and the commutator class falls below the odd class; the Letter has $1<\sqrt2<\sqrt3<\sqrt6$. Its rescaled trace vector (out §4, "best (c) miss, detail") falls short on $(1,3)_-$ (6.21 against $3+3\sqrt6=10.35$), on $(27,1)_+$ (39.8 against 46.8) and on $(1,1)_+$, and overshoots the odd components.
- **Theorem 4 check.** 222 factorizing words have an $E'$-scalar twirled mixed block, and none violates $|\mu_+|=|\mu_-|$. Of the 1616 non-factorizing words with an $E'$-scalar twirled mixed block, 72 reach other ratios: $\pm\tfrac12,\pm\tfrac35,\pm\tfrac53,\pm2$, 9 words each. The printed examples combine $\Pi_\pm$ with a partial transpose and a partial trace (e.g. $T_{q_1}\Pi_-\tau_{q_1}$, $T_{q_2}\Pi_+T_{q_1}$). **None of the 72 is class-scalar after the twirl** (residuals 0.22–0.69 in the examples). Even when the mixed sector alone has the right ratio, the other classes do not follow.

## 4. Family F2: survival × attribution dressings

**Declared family.** This is the operator form of the Letter's mechanism, with every ingredient a natural operator on $B(H)$:
$$K=\prod_{x\in\{X,Y\}}\bigl(1+(d_Q-1)S_xA_x^Q+(d_q-1)S_xA_x^q\bigr),\qquad M=\sqrt K=\prod_x\bigl(1+(a_3-1)S_xA_x^Q+(a_2-1)S_xA_x^q\bigr).$$
For each exchange $x$:
- **survival operator** $S_x$, one of 7: $1$, $\Pi_+(E')$, $\Pi_+(\sigma E')$, $\Pi_+(E'_Q)$, $\Pi_+(\sigma E'_Q)$, $\Pi_+(E'_q)$, $\Pi_+(\sigma E'_q)$;
- **attribution**, one of 81: each occupation type (QQ, qq, mix, comm) is credited to Q, to q, or to neither.

That gives $(7\cdot81)^2=321\,489$ ordered choices. All ingredients commute and preserve $R$. They have 9 joint eigenspaces on $R$, of dimensions 28, 36, 3, 24, 24, 24, 24, 3, 6 (out §5). This makes the family exactly computable, with exact numerators in $\mathbb Z[\sqrt2,\sqrt3]$.

**Results (out §5).**
- 291 distinct single-exchange factors and **19 926 distinct dressings** on $R$, taking 11 570 distinct values.
- **$M_L$ is reached by exactly one distinct dressing**, realized by 36 ordered choices. Every one credits $d_Q=3$ to the E′-even q₁q₂ records and to the commutator records (checked: True), so none is local in the original frame.
- **Local dressings:** 1992 distinct, with 1054 distinct values, $1/\text{value}\in[61.71,180.84]$. **None equals $M_L$** (Prop. 5).
- **All nine survival/attribution rows of Table S5 lie in F2**, plus "no weights". The other weight rows (Born weights, CT-G08) and the record-set rows (commutator sums kept, no commutator records) are outside F2 by construction.

  | Table S5 row | numerator | 1/value | local? |
  |---|---|---|---|
  | computed (Letter) | $79+48\sqrt2+42\sqrt3+3\sqrt6$ | 137.036 | no |
  | two ends swapped | $79+42\sqrt2+48\sqrt3+3\sqrt6$ | 135.894 | no |
  | both at qubit end | $85+90\sqrt2$ | 146.524 | no |
  | both at qutrit end | $88+90\sqrt3$ | 127.536 | no |
  | X where letters sit | $85+54\sqrt2+36\sqrt3$ | 139.030 | **yes** |
  | qq moved by Y | $79+48\sqrt2+36\sqrt3+9\sqrt6$ | 134.485 | no |
  | copies separate | $79+42\sqrt3+51\sqrt6$ | 112.423 | no |
  | no own sign | $82+48\sqrt2+42\sqrt3$ | 139.713 | no |
  | comm survive X only | $79+48\sqrt2+45\sqrt3$ | 138.348 | no |
  | no weights | $172$ | 180.837 | yes |

  At the operator level, the Letter's choice and its listed alternatives stand on the same footing. Among them, only "X where the letters sit" credits each end's dimension to records that have a letter there.
- **Closest local misses to $M_L$** (HS distance $\lVert M-M_L\rVert/\lVert M_L\rVert$):

  | distance | numerator | 1/value | X | Y |
  |---|---|---|---|---|
  | 0.0570 | $97+48\sqrt2+36\sqrt3$ | 136.880 | $\Pi_+(\sigma E')$; qq, comm → q | $\Pi_+(\sigma E')$; QQ → Q, qq, mix, comm → q |
  | 0.0616 | $85+54\sqrt2+36\sqrt3$ | 139.030 | Table S5 "X where letters sit" | |
  | 0.0697 | $94+51\sqrt2+36\sqrt3$ | 136.135 | | |

  They agree with the Letter on the 64 QQ and 96 mixed records. They differ only on the 12 qubit-only records, where the local rule gives $\sqrt2$ or $2$ instead of $1$, $\sqrt3$ and $\sqrt6$: all three on the 6 E′-even qq and the 3 commutator records, and the third also gives $\sqrt2$ to the 3 E′-odd qq records, whose Letter weight is 1 ($51=48+3$).
- **Values near $1/\alpha$.** 49 F2 dressings (the Letter's among them) lie within $10^{-3}$ relative of $1/\alpha$; none of them is local. The closest local value is 137.186 ($1.1\times10^{-3}$ relative). Multiplying the local density 6.89 per unit by the window width gives the reference count 1.89 under an illustrative local-uniformity assumption. This is not a calibrated probability or significance test, and the reference count is not used as evidence.

## 5. Limits

- "Natural" here means declared. F1's alphabet and length bound and F2's survival and attribution sets are our choices. Lemma 3 extends the $\sqrt{k_r}$, per-isotypic-trace and value negatives of F1 to words of any length and to all $\mathbb Q(i)$-linear combinations. Theorem 4 extends the $k_r$ negative to all products of factorizing letters and their twirls. For words that contain $\Pi_\pm$ or a read link, the $k_r$, $|C|$ and $\sqrt C$ negatives are covered by neither theorem and hold up to length 3 only (F1). Prop. 5 extends the F2 negative to every product-of-dimensions rule local in the original frame.
- A sum of words with amplitude tags can reproduce $M_L$. Prop. 2 is one, and by Lemma 3 every such construction must contain the four class projectors with the Letter's tags. That is the content of the weight law and is not derived.
- The transport direction (carry qubit letters to the qutrit copies through the read links, then credit the join's three directions; it breaks $G_\Omega$ to $SU(2)\times U(1)$) is tested in `transport.md`. A frame-local dressing in the transported frame reproduces $M_L$, but only with the frame of each exchange and the full-join credit as inputs (Theorems 5, 6 there). No compression $T^\dagger DT$ and no Rayleigh quotient of any $D$ through a transport does (Theorem 4 there), and no word of length ≤ 3 with a transport letter does (Prop. 8 there).
- The value $1/\alpha_\text{cand}$ is not affected. Nothing here changes the counts, the value, or F-CT-1–F-CT-10.

## 6. Consequences for the gaps

- **CT-G01** (status of the weight law). The candidate resolution as stated ("construct $M$ as an explicit operator on $R$ ... show that its diagonal on the record basis is $\sqrt{k_r}$") is met by Prop. 2 without closing the gap. The gap is now precise: *an operator construction in which the end attribution is not an input.* For example, one local in the original frame, which Prop. 5 shows cannot give the Letter's numbers. The $G_\Omega$-breaking transport construction gives them only with the frame assignment and the full-join credit as inputs (`transport.md`). CT-G01 therefore reduces to CT-G02 (plus CT-G08 for the square root).
- **CT-G12** (freedom left by symmetry). The two copies of $(1,3)_-$ are separated by the natural operator $\hat N$. Symmetry plus partial traces and the exchanges leave the weights on the 9 joint eigenspaces of Section 4 free. The Letter's point is one choice among them and, by Prop. 5, one that is not local in the original frame (it is frame-local after transport, `transport.md` Prop. 7).
