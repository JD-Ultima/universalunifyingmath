# Transporting the qubit letters to the qutrit copies (gaps CT-G01, CT-G02)

**Object.** The weight law of the α Letter (Wissner-Gross, v1.0.0, Eq. (5), Supplement Thm S-gain): record $r$ has weight $\sqrt{k_r}$, with $k_r=1$ on E′-odd products, $2$ on E′-even mixed, $3$ on E′-even same-end (the 6 E′-even qubit–qubit records included), $6$ on the 3 commutator records. As in `operator-origin.md`, the records are the 172 covariant ones (exchange $E'=E\circ T$, F-CT-10), and $M_L=\operatorname{diag}\sqrt{k_r}$, $K_L=M_L^2$.
**Question.** `operator-origin.md` §5 left one direction untested. First carry the qubit letters to the qutrit copies through the read links $|\omega_2\rangle_{Q_iq_i}$ (support: levels 0, 1 of $Q_i$), by the transpose identity $(A\otimes1)|\omega\rangle=(1\otimes A^{\mathsf T})|\omega\rangle$. Then credit the three directions of the join $|\omega_3\rangle_{Q_1Q_2}$. Can a construction of this kind give $3$ and $6$ on the 12 qubit-only records and $1,2,3$ on the other records? It is tested as an operator, as its diagonal, and as per-isotypic traces.
**Code.** `computation/chain-theory/alpha/transport.py` → `transport.out`. It loads the definitions and tests of `operator_origin.py`, which itself imports `records_independent.py` and `symmetry.py`; none of the three is modified. "out §n" means section n of the output.
**Date.** 2026-10-30. Exploratory analysis toward CT-G01/CT-G02. The audited value is not changed.

## Result in brief

1. **The transport exists, is unique, and has less symmetry than the join** (Lemma 1, Prop. 2).
   - The only traceless Hermitian operator at $Q_i$ that acts on the read link as a qubit letter $B$ at $q_i$ is $t(B)=\iota B^{\mathsf T}\iota^\dagger$. So $t(\sigma_x)=\lambda_1$, $t(\sigma_y)=-\lambda_2$, $t(\sigma_z)=\lambda_3$, all supported on the 2-dim support $S$.
   - The subgroup of $G_\Omega$ that fixes both read links is $H\cong SU(2)\times U(1)$ (dim 4). Up to scale, $t$ is the only $H$-equivariant map from qubit letters to qutrit letters. The only $G_\Omega$-equivariant one is $0$.
2. **Transporting every record is many-to-one** (Theorem 4). Under each of the six declared transports, cross-mixed records land on QQ records, and commutator records on images of same-link mixed records. Colliding records then receive tied weights, for every operator $D$ applied in the transported frame. This rules out $M_L$ and $K_L$ both as $\operatorname{diag}(T^\dagger DT)$ and as Rayleigh quotients of $D$.
3. **Transported qubit letters reach two of the join's three directions** (Theorem 5).
   - Under any symmetry that the transport respects, the smallest invariant support of a transported letter is $S$, of rank 2. The credit is then $\operatorname{Tr}(\Pi_SC_{\omega_3})=2/\sqrt3$, or $\sqrt2$ after renormalizing the compressed pair, not $\sqrt3$.
   - Rank 3 is forced only by groups that move $S$, such as $SU(3)$, and under those the transport is not equivariant.
   - 6 of the 12 transported qubit-only records annihilate the join altogether.
4. **Frame-local dressings** (family F3, 338 688 ordered choices; Theorem 6, Prop. 7).
   - $M_L$ is reached by exactly 6 ordered choices. Every one assesses one exchange in the transported frame with the full-join credit $\sqrt3$ (J1), and the other in the original frame at the qubit end.
   - In these, the attribution of each record class is derived from the frame and locality; the Letter does not have to state it class by class.
   - No member with the transport-compatible credits J2 ($2/\sqrt3$) or J3 ($\sqrt2$) reaches $M_L$, or even its per-isotypic traces. The Letter-shaped members then give $1/\text{value}=140.702$ and $139.030$.
5. **Words with transport letters** (family F3a: 22 635 words containing a transport letter, 5953 distinct compressions analysed). No word is class-scalar or has a class-constant diagonal in any test. The best per-isotypic residual is 0.0947, against 0.2070 for the identity and 0.1008 for the best word of F1.

**Precise negative statement.**
- No linear construction through the transport reproduces $M_L$ or $K_L$. Uniform transport fails by Theorem 4, for any $D$. Selective transport needs the selector $[\hat N_Q=0]$ of the 12 qubit-only records, which is the content of $\Pi_X$ again (Prop. 3(5)).
- Among the record-wise frame rules F3, $M_L$ requires all three of the following:
  - (a) the exchange of reads assessed in the transported frame;
  - (b) the exchange of links assessed in the original frame;
  - (c) the rank-3 credit J1, given to letters that are supported on a rank-2 subspace.
  
  Each is a choice. (a) and (b) are forced by $M_L$ itself (Theorem 6), not derived. (c) is the $G_\Omega$-invariant credit, while the transport is only $H$-covariant.
- No word in F3a reproduces $\sqrt{k_r}$ or $k_r$ as operator, diagonal or per-isotypic traces. This holds also after the $G_\Omega\times\langle E'\rangle$-twirl, and for $|C|$ and $\sqrt C$.

**Consequence.** The transport gives the end rule an operator-level reading: "the exchange of reads is assessed where every letter has been carried to the end at which its read is taken." In that frame, crediting the qutrit to the qubit–qubit records is local: their transported images *are* QQ records of the same E′-parity (Prop. 3(1)), and the pull-back of $M_L$ agrees with $M_L$ on them. The cost is two separate choices. One is which exchange is assessed in which frame. The other is whether the join is credited with the symmetry of the join ($\sqrt3$) or with that of the transport ($2/\sqrt3$ or $\sqrt2$). The Letter's value needs the first option in both. CT-G01 is not closed. It reduces to CT-G02, now as these two choices, and to CT-G08.

---

## 1. The transport at the letter level

Notation: $S=\operatorname{span}\{|0\rangle,|1\rangle\}\subset\mathbb C^3$ is the common support of a read link, $\iota:\mathbb C^2\to\mathbb C^3$ the inclusion, $\Pi_S=\iota\iota^\dagger$, $C_{\omega_d}=\mathbb 1_d/\sqrt d$ the coefficient matrix of $|\omega_d\rangle$. The Gell-Mann matrices are numbered in the standard way.

**Lemma 1 (uniqueness of the transport; exact, out §1).** For $B\in M_2$, the solutions $A\in M_3$ of $(A\otimes1)|\omega_2\rangle=(1\otimes B)|\omega_2\rangle$ on $\mathbb C^3_{Q_i}\otimes\mathbb C^2_{q_i}$ are exactly
$$A=\begin{pmatrix}B^{\mathsf T}&a\\0&b\end{pmatrix},\qquad a\in\mathbb C^2,\ b\in\mathbb C .$$
1. If $B$ is traceless Hermitian, the unique traceless Hermitian solution is $t(B)=\iota B^{\mathsf T}\iota^\dagger$.
2. For $B=1$, the Hermitian solutions are $\Pi_S+b\,|2\rangle\langle2|$ with $b\in\mathbb R$ free. The read link does not fix what the transported identity does at the third direction.
3. $t(\sigma_x)=\lambda_1$, $t(\sigma_y)=-\lambda_2$, $t(\sigma_z)=\lambda_3$. Every $t(B)$ satisfies $t(B)=\Pi_St(B)\Pi_S$ and $t(B)|2\rangle=0=\langle2|t(B)$. For $B\neq0$ traceless Hermitian, $t(B)$ has rank 2 and support exactly $S$.

*Proof.* Compare components. $(A\otimes1)|\Phi_2\rangle$ has component $A_{ji}$ at $|j\rangle|i\rangle$ ($i<2$), and $(1\otimes B)|\Phi_2\rangle$ has $B_{ij}$ there for $j<2$ and $0$ for $j=2$. So $A_{ji}=B_{ij}$ on $S$, $A_{2i}=0$, and the column of $|2\rangle$ is free.
- (1) Hermiticity forces $a=\overline{(A_{20},A_{21})}=0$ and $b\in\mathbb R$; tracelessness forces $b=-\operatorname{tr}B=0$.
- (2) Set $B=1$.
- (3) A nonzero traceless Hermitian $2\times2$ matrix has eigenvalues $\pm\lambda\neq0$.

`sympy` solves both systems in out §1. ∎

**Proposition 2 (the symmetry of the transport; out §1).** Let $H\subset G_\Omega$ be the subgroup fixing both read-link projectors $P_{L_1},P_{L_2}$.
1. $\dim\operatorname{Lie}H=4$. Explicitly, $H=\{(U,V):U=c\bar V\oplus c^{-2}\}\cong SU(2)\times U(1)$. This is the diagonal $su(2)$ on $S$ plus $\lambda_8$; compare the Supplement's stabilizer of the pair ket, "$\zeta\oplus\bar\zeta$ … on the two-dimensional common support" (§S3, Thm S-stabilizer). $E'$ exchanges $P_{L_1}$ and $P_{L_2}$.
2. On $\mathbb C^3$, $H$ has commutant $\operatorname{span}\{\Pi_S,|2\rangle\langle2|\}$, so its invariant subspaces are $0$, $S$, $\mathbb C|2\rangle$ and $\mathbb C^3$. $SU(3)$ has commutant $\mathbb C\,1$, with invariant subspaces $0$ and $\mathbb C^3$.
3. Consider the linear maps from qubit letters at $q_i$ to qutrit letters at $Q_i$ that intertwine a group. For $H$ they form a 1-dim space, spanned by $t$. For $G_\Omega$ there are none (the space is 0).
4. Let $G'\le G_\Omega$ be any subgroup for which $t$ is equivariant. Then $G'\le H$ and $S$ is $G'$-invariant.

*Proof.*
1. $(U\otimes V)|\omega_2\rangle\propto|\omega_2\rangle$ iff $U\iota V^{\mathsf T}\propto\iota$, that is, $US=S$ and $U|_S\propto\bar V$. The second link gives the same condition, since there $G_\Omega$ acts by $\bar U\otimes\bar V$. Then use $\det U=1$.
2. $\mathbb C^3=S\oplus\mathbb C|2\rangle$ under $H$, a sum of two inequivalent irreducibles.
3. Under $H$, the qubit letters form the adjoint $\mathbf 3$ of $SU(2)$. The qutrit letters decompose as $\mathbf 8=\mathbf 3_0\oplus\mathbf 2_{+}\oplus\mathbf 2_{-}\oplus\mathbf 1_0$, where $\mathbf 3_0$ is the block on $S$ and the subscripts are $U(1)$ charges. So the intertwiners form a 1-dim space. Under $G_\Omega$ the two sides are $(\mathbf 1,\mathbf 3)$ and $(\mathbf 8,\mathbf 1)$, which are inequivalent (Schur).
4. Equivariance means $Ut(B)U^\dagger=t(VBV^\dagger)$. The left side has support $US$ and the right side $S$, so $US=S$. Then $\iota^\dagger U\iota$ intertwines $B^{\mathsf T}\mapsto(VBV^\dagger)^{\mathsf T}$, so it is $\propto\bar V$ by Schur on $su(2)$.

Numerical checks in out §1: $\dim\mathfrak h=4$, relation residual $1.6\times10^{-16}$, kernel of the projection to $su(2)_q$ equal to $\lambda_8$ (defect $1.9\times10^{-16}$), commutants 2 and 1, intertwiner spaces 1 and 0. ∎

## 2. The transport on the record space

**Declared transports.** At read link $i$, with $b\in\{0,1\}$, let $t_b(B)=\iota B^{\mathsf T}\iota^\dagger+b\,\tfrac{\operatorname{tr}B}{2}|2\rangle\langle2|$. On $B(\mathbb C^3_{Q_i}\otimes\mathbb C^2_{q_i})$ define
$$A\otimes B\ \mapsto\ A\,t_b(B)\ (\text{R}),\quad t_b(B)\,A\ (\text{L}),\quad \tfrac12\{A,t_b(B)\}\ (\text{S}),$$
and leave the identity at $q_i$.
- With $b=0$, these maps are literally the read link plus a partial trace: R$=4\,\tau_{q_i}\circ R_{P_{L_i}}$, L$=4\,\tau_{q_i}\circ L_{P_{L_i}}$ (defect 0, out §2).
- With $b=1$, the transported identity is $1_3$, and native qutrit letters are left in place.
- The full transports apply both links. The E′-covariant ones are the orderings (R,L), (L,R) and (S,S); (R,R) is not, since $T(At(B))=t(B)^{\mathsf T}A^{\mathsf T}$ turns R into L.

That gives **six transports** $T_b(o_1,o_2)$. The reference is $T_*=T_1(\text{R},\text{L})$.

**Proposition 3 (out §2, §5).** Each of the six transports leaves no letter at a qubit slot (max q-occupation 0) and keeps E′-parity. It is $H$-equivariant (defect $\approx3\times10^{-15}$) and not $G_\Omega$-equivariant (defect 1.2–1.4). Under $T_*$:
1. the 9 qubit–qubit records go, up to scale, onto the 9 QQ records whose letters lie in $\{\lambda_1,\lambda_2,\lambda_3\}$ (6 E′-even, 3 E′-odd; parity kept);
2. the 48 cross-mixed records (letters at $Q_i$ and $q_j$, $i\ne j$) go onto QQ records of the same parity, except 3 E′-odd ones, which vanish;
3. the 48 same-link mixed records and the 3 commutator records go to single-slot qutrit operators outside $R$;
4. the QQ records are fixed, and $T_*$ has rank 76 on $R$ (36 for $T_0(\text{R,L})$, where $A\mapsto A\Pi_S$ truncates native letters);
5. *selective transport.* $T'=T_*[\hat N_Q=0]$ has rank 12 on $R$, with support exactly the qq and commutator records. On $R$, $\Pi_X=\Pi_X[\hat N_Q\ge1]+\operatorname{supp}(T'^\dagger T')$. Prop. 2's $\Pi_X$ thus becomes "qutrit letters present, before or after transport", but the selector $[\hat N_Q=0]$ is the list of the 12 qubit-only records.
6. The linear pull-back $T_*^\dagger[\hat N_Q\ge1]T_*$ is not diagonal in the record basis (off-diagonal part 0.757 of its norm).

So the transported frame can be read record by record (Section 4), but not as one linear operator diagonal on $R$.

**Theorem 4 (fiber obstruction).** Let $r,r'$ be unit records with $Tr=c\,Tr'\neq0$. For every superoperator $D$ on $B(H)$:
$$\langle r,T^\dagger DT\,r\rangle=|c|^2\langle r',T^\dagger DT\,r'\rangle,\qquad \frac{\langle Tr,DTr\rangle}{\langle Tr,Tr\rangle}=\frac{\langle Tr',DTr'\rangle}{\langle Tr',Tr'\rangle}.$$
Under $T_*$ there are 168 colliding pairs, and they rule out every reproduction:

| $r$ | $r'$ | $\kappa=\lvert c\rvert^2$ | $M_L$ ratio | $K_L$ ratio |
|---|---|---|---|---|
| odd cross-mixed (21) | odd QQ | $2/3$ | 1 | 1 |
| even cross-mixed (24) | even QQ | $2/3$ (21), $4/3$ (3) | $\sqrt{2/3}$ | $2/3$ |
| even qq (6) | even QQ | $4/9$ | 1 | 1 |
| odd same-link mixed (9) | commutator | $1/2$, $3/2$ | $1/\sqrt6$ | $1/6$ |

- The diagonal form fails on rows 1 and 3, where $\kappa\neq1$ but both $M_L$ and $K_L$ need 1.
- The Rayleigh form fails on rows 2 and 4, where the Letter's weights differ: $\sqrt2$ vs $\sqrt3$, and $1$ vs $\sqrt6$.

Hence no $D$ gives $M_L$ or $K_L$ in either form. The same holds for each of the six transports (out §3).

*Proof.* Linearity. Random-$D$ check: $7.7\times10^{-15}$. ∎

*Remark.* On 21 of the 24 even pairs, the norm loss $\kappa=2/3=d_q/d_Q$ is exactly the $K_L$ ratio $2:3$ between even mixed and even same-end records. A qubit letter carried onto a qutrit trades $\operatorname{tr}1_3=3$ for $\operatorname{tr}1_2=2$. But the same $2/3$ is forced on the odd pairs, where $K_L$ has $1:1$.

*Pull-back of $M_L$ (out §3).* $\rho(r)=\langle Tr,M_LTr\rangle/\langle Tr,Tr\rangle$ is defined on the records whose image lies in $R$. It equals $M_L$ on the 64 QQ, the 9 qq and the 21 odd cross-mixed records. On the 24 even cross-mixed records it gives $\sqrt3$ instead of $\sqrt2$. In the transported frame, the qubit–qubit records inherit the QQ weights. The weight $\sqrt2$ of the exchange of links does not survive transport.

## 3. What the transported records see of the join

**Theorem 5 (two of three directions).** Let $A=t(B)\ne0$ be a transported qubit letter, and $G'\le G_\Omega$ any subgroup for which the transport is equivariant.
1. The smallest $G'$-invariant projector containing $\operatorname{supp}A$ is $\Pi_S$, of rank 2.
2. The amplitude summed over the matched directions that this projector covers is $a=\operatorname{Tr}(\Pi_SC_{\omega_3})=2/\sqrt3$, with square $4/3$ (credit **J2**).
3. Compressing the join to the support and renormalizing gives $(\Pi_S\otimes\Pi_S)|\omega_3\rangle=\sqrt{2/3}\,|\omega_2\rangle_{S\otimes S}$, with amplitude $\sqrt2$ and square 2 (credit **J3**). This is also the amplitude of the qubit pair that the joint state carries from $q_1$ to $q_2$ through link, join and link (out §4: the $q_1q_2$ marginal of $|\Omega\rangle$ is pure, $\operatorname{Tr}C=\sqrt2$).
4. The Letter's $\sqrt3=\operatorname{Tr}(1_3C_{\omega_3})$ (credit **J1**) needs the projector $1_3$. It is $G'$-invariant but not minimal: it adds the line $\mathbb C|2\rangle$, on which $A$ vanishes from both sides.
5. $1_3$ is the smallest invariant projector containing $S$ only for groups that do not preserve $S$, such as $SU(3)$. By Prop. 2(4), the transport is not equivariant under any such group.

*Proof.* (1) By Prop. 2(4), $S$ is $G'$-invariant, and by Lemma 1(3), $\operatorname{supp}A=S$. (2) and (3) are computations (out §4: $2/\sqrt3=1.154700538379$, $\sqrt2=1.414213562373$, compression defect $1.6\times10^{-16}$). (4) is Lemma 1(3) and Prop. 2(2). (5): a group preserving $S$ has $\Pi_S$ as an invariant projector containing $S$; then use Prop. 2(4). ∎

*Further facts (out §4).*
- *Transport variants.* For all six transports, every transported qubit-only record vanishes against $|22\rangle\langle22|_{Q_1Q_2}$ from both sides. With $b=1$, the letter supports of all 12 records lie in $S$. With $b=0$, the identity left at the other copy becomes $\Pi_S$. Its traceless part, which is proportional to $\lambda_8$, reaches $|2\rangle$ on the 3 commutator records. That is an artefact of the truncation, not a letter.
- *Annihilation of the join.* 6 of the 12 transported records annihilate the join, $T_*r\,|\omega_3\rangle=0$, by the transpose identity. These are the 3 commutator records (E′-odd single-slot sums $X^{Q_1}-X^{\mathsf T Q_2}$) and the 3 E′-even qq records with two different, anticommuting letters. The Letter credits exactly these with $\sqrt6$ and $\sqrt3$. 8 native QQ records also annihilate the join, and the Letter gives all E′-even QQ records $\sqrt3$ regardless. So "acting on the join" is not what the Letter's credit counts.
- *The $H$ credit applied to every record.* Assess X in the transported frame with squared credit $(\text{rank of the smallest }H\text{-invariant projector containing the record's letter supports})^2/3$, and Y as the Letter does. The X-survivors then split into 30 QQ records of rank 3 and, at rank 2, 6 QQ records with letters in $su(2)_S$, the 6 qq records and the 3 commutator records. This gives numerator $79+48\sqrt2+38\sqrt3+2\sqrt6$, $1/\text{value}=142.942$, distance $0.137$ from $M_L$, $G_\Omega\times\langle E'\rangle$-twirl defect $0.075$, and per-isotypic residual $0.046$. Under $H$, native QQ records whose letters lie on $S$ lose the third direction too.

## 4. Frame-local dressings (family F3)

**Definition (frame-local dressing).** Each exchange $x\in\{X,Y\}$ gets three things:
- a survival operator $S_x$, one of the 7 of `operator-origin.md` F2: $1$, $\Pi_+(E')$, $\Pi_+(\sigma E')$, $\Pi_+(E'_Q)$, $\Pi_+(\sigma E'_Q)$, $\Pi_+(E'_q)$, $\Pi_+(\sigma E'_q)$;
- a domain $D_x\subseteq\{\text{QQ},\text{qq},\text{mix},\text{comm}\}$, one of 16;
- a frame $\varphi_x$, one of 3:
  - **orig-Q**: original frame, crediting the end where the record's letters sit (Q when both ends carry letters);
  - **orig-q**: the same, but q when both ends carry letters;
  - **tr**: after transport, where every letter sits at a qutrit copy.

A record $r\in D_x$ that survives $S_x$ gets a squared credit:
- 2 at a qubit end, 3 at a qutrit end with a native letter;
- in the tr frame, for records whose qutrit letters are all transported (qq and comm), the squared credit $j$ of a global choice J1 ($j=3$), J2 ($4/3$) or J3 ($2$).

Otherwise the credit is 1. Then $K=\prod_x(\text{credit})$ and $M=\sqrt K$. Every member is local in its own frame.

**Size.** $7\cdot16\cdot3=336$ choices per exchange, $336^2=112\,896$ ordered pairs, times 3 credits: **338 688 ordered choices**. They give 7504 distinct dressings (J1 and J2: 4748 each; J3: 1992) and 4186 distinct values. All ingredients are functions on the 9 joint eigenspaces of $(E',E'_Q,E'_q,\sigma,\text{type})$ (dimensions 28, 36, 3, 24, 24, 24, 24, 3, 6), so the enumeration is exact. Per-isotypic traces are obtained from the class content of these eigenspaces.

**Theorem 6 (both frames, and the full-join credit, are forced).** If a frame-local dressing equals $M_L$ on the 172 records, then:
- (a) some exchange is assessed in the original frame with qubit priority (orig-q) and credits the qubit end to the 48 E′-even mixed records;
- (b) the *other* exchange is assessed in the transported frame with credit J1, and credits the 6 E′-even qq records;
- (c) the 3 commutator records get 2 from (a)'s exchange and 3 from (b)'s.

So the two exchanges are assessed in different frames, and the credit is J1.

*Proof.*
- E′-even mixed records have a native qutrit letter, so each credit they can receive is 2 (orig-q, qubit end) or 3 (orig-Q, or tr). $k=2$ needs exactly one credit 2 and no credit 3.
- E′-even qq records can receive 2 (original frame) or $j$ (tr). Products of at most one factor from each exchange lie in $\{1,2,4,j,2j,j^2\}$, which contains 3 only if $j=3$ and exactly one factor is $j$.
- One exchange has one frame. Hence the exchange that gives evmix its 2 (original) differs from the one that gives qq-even its 3 (tr).
- Commutator records carry qubit letters only, so they can receive 2 (original) or $j=3$ (tr). $k=6$ needs one of each.

∎

**Proposition 7 (F3 results, out §5).**
1. $M_L$ is reached by **6 ordered choices** (3 unordered), all with J1 and frames {tr, orig-q}. Example: $\Pi_+(\sigma E')$ on {QQ, qq, comm} in tr, with $\Pi_+(\sigma E')$ on {mix, comm} in orig-q (the Letter's survivals and domains). Two original frames never give $M_L$ (Prop. 5 of `operator-origin.md`), and neither do two transported frames.
2. **Per-isotypic traces.** Only $M_L$ itself, under J1, has traces proportional to the Letter's point. No member does under J2 or J3.
3. **Letter-shaped members** (the Letter's survivals and domains, frames and credit varied):

   | credit | X frame | Y frame | numerator / 31104 | 1/value | $\lVert M-M_L\rVert/\lVert M_L\rVert$ |
   |---|---|---|---|---|---|
   | J1 | tr | orig-q | $79+48\sqrt2+42\sqrt3+3\sqrt6$ | **137.036** | 0 |
   | J2 | tr | orig-q | $79+48\sqrt2+40\sqrt3+2\sqrt6$ | **140.702** | 0.112 |
   | J3 | tr | orig-q | $85+54\sqrt2+36\sqrt3$ | **139.030** | 0.062 |
   | any | orig | orig-q | $85+54\sqrt2+36\sqrt3$ | 139.030 | 0.062 |
   | J1 | tr | tr | $88+90\sqrt3$ | 127.536 | 0.134 |
   | J2 | tr | tr | $83+88\sqrt3$ | 132.121 | 0.182 |
   | J3 | tr | tr | $85+6\sqrt2+84\sqrt3$ | 130.154 | 0.138 |
   | J1 | orig | tr | $79+6\sqrt2+84\sqrt3+3\sqrt6$ | 129.424 | 0.131 |
   | J1 | tr | orig-Q | $79+90\sqrt3+3\sqrt6$ | 128.405 | 0.123 |

   Credit J3 in the transported frame reproduces Table S5's "X where the letters sit" (139.030) exactly. Carrying the qubit letters to the qutrit and crediting the two directions they reach gives the same number as not carrying them. J1 with both exchanges transported is Table S5's "both at the qutrit end" (127.536).
4. **Closest non-$M_L$ members.**
   - J1: $79+54\sqrt2+36\sqrt3+3\sqrt6$ (138.197, distance 0.044). It differs from $M_L$ only on the 6 E′-even qq records, which get $\sqrt2$.
   - J2 and J3: the original-frame local misses of `operator-origin.md`, $97+48\sqrt2+36\sqrt3$ (136.880, 0.057), $85+54\sqrt2+36\sqrt3$ (139.030, 0.062) and $94+51\sqrt2+36\sqrt3$ (136.135, 0.070).
5. **Table S5 rows inside F3.**
   - Under J1: the Letter, "both at qutrit end", "X where letters sit", "qq moved by Y", "copies separate", "no own sign", "comm survive X only", "no weights".
   - Under J2/J3: only "X where letters sit" and "no weights".
   - Not in F3 under any credit: "two ends swapped" and "both at qubit end". They credit the qubit end to QQ records, which is local in no frame.
6. **Values near $1/\alpha$.** 5 of the 4186 distinct values lie within $10^{-3}$ relative of $1/\alpha$: the Letter's and four others (137.018, 137.010, 137.062, 137.066), all under J1. This count is not used as evidence.

## 5. Words with transport letters (family F3a)

**Declared family.** The alphabet has 33 letters:
- the 24 $G_\Omega$-covariant primitives of `operator-origin.md` §1;
- 9 transport letters: $T_0(\text{R,L})$, $T_0(\text{S,S})$, $T_1(\text{R,L})$, $T_1(\text{S,S})$ and their Hilbert–Schmidt adjoints, and the support compression $C_S=\mathrm{Ad}_{\Pi_S\otimes\Pi_S}$.

All 9 are E′-covariant and $H$-equivariant, and $\mathbb Q$-rational (denominators 1, 4 or 16). The eight transports are not Q|q products (Schmidt rank 16); $C_S$ is one (rank 1). The words are all words of length ≤ 3: $37\,060$ in total, $22\,635$ of them containing a transport letter. They are enumerated as in F1: left extensions of distinct images, with duplicates up to scale merged. This gives 7781 distinct images; 6194 are first reached by a word with a transport letter, and 5953 of these are nonzero on $R$ and analysed with `operator_origin.py`'s `analyse` (the same tests as F1).

**Proposition 8 (out §6).**
1. No word is class-scalar on $R$ (raw, $G_\Omega\times\langle E'\rangle$-twirled or $|C|$), none has a class-scalar $\sqrt C$ (352 PSD words), and none has a class-constant diagonal. Hence none matches $(1,\sqrt2,\sqrt3,\sqrt6)$ or $(1,2,3,6)$.
2. Minimum residuals, with the identity in brackets:

   | test | residual |
   |---|---|
   | diagonal vs $\sqrt k$ | 0.3760 (0.2471) |
   | diagonal vs $k$ | 0.5113 (0.4674) |
   | per-isotypic traces | 0.0947 (0.2070) |
   | $\lvert C\rvert$ traces | 0.0940 (0.2070) |
   | $\sqrt C$ traces | 0.1007 (0.2070) |

   The F1 minima are 0.2471, 0.4572, 0.1008, 0.0912 and 0.1007. Transport letters improve only the raw per-isotypic residual, from 0.1008 to 0.0947.

   The best trace miss, $T_0(\text{R,L})^\dagger(1-\tau_{Q_2})T_0(\text{R,L})$, has class means $(1,1.21,1.33,0.88)$. Normalized so that odd records weigh 1, it gives $1/\text{value}=159.258$. It is non-invariant and has the wrong class ordering (the commutator class is below the odd class).
3. By Lemma 3 of `operator-origin.md`, which applies since all letters are rational, no word of any length matches $\sqrt{k_r}$ on the diagonal or in the per-isotypic traces. Theorem 4 there does not apply to words containing a transport, because the eight transports are not Q|q products. The $K$-route negative of item 1 is therefore new evidence, and it holds only up to length 3.

## 6. What is derived, what is input, what symmetry survives

- **Derived (record-wise).** In the 6 F3 choices that give $M_L$, the end credited to each class follows from the frame and locality. In particular, the qutrit credit of the 9 qq records follows from the operator identity $T_*(\text{qq})=$ QQ records (Prop. 3(1)).
- **Input.**
  - (i) The frame of each exchange: X transported, Y original. By Theorem 6 this is forced by $M_L$, but nothing in the construction selects it. It restates the end-rule clause ("a read is an act at the end at which it is taken") and the copies clause ("the exchange of links moves only what sits at the far ends").
  - (ii) The credit at the join, J1. It is the $G_\Omega$-invariant credit, while the transport is only $H$-covariant. Under $H$, the minimal invariant support of a transported letter has rank 2 (Theorem 5).
  - (iii) At the level of linear operators, the selector $[\hat N_Q=0]$ (Prop. 3(5)), or else the fiber obstruction (Theorem 4).
- **Symmetry.**
  - Transport: $H\times\langle E'\rangle$, with $H\cong SU(2)\times U(1)$ of dimension 4.
  - Credit J1: $G_\Omega$-invariant.
  - Resulting $M$: in the $G_\Omega\times\langle E'\rangle$ commutant, since it is a function of the joint eigenspaces.
  - The $H$-consistent alternative (the $H$ credit for every record, §3): only $H\times\langle E'\rangle$-invariant (twirl defect 0.075).

## 7. Limits

- **Declared families.** The families are declared choices:
  - six transports ($b\in\{0,1\}$, three E′-covariant orderings);
  - three frames (original with Q- or q-priority, fully transported);
  - three credits;
  - F2's seven survivals and 16 domains;
  - F3a's 33 letters up to length 3.

  Theorems 4, 5 and 6 are general within their stated hypotheses. Theorem 4 covers every $D$, Theorem 5 every subgroup for which the transport is equivariant, and Theorem 6 every frame-local product-of-credits dressing. The word negative (Prop. 8) holds only up to length 3.
- **Not tested.**
  - Half transports (one link at a time, averaged). They superpose occupation types, so "where the letters sit" is not sharp.
  - Reverse transport (qutrit letters compressed to the qubits). It annihilates $\lambda_4,\dots,\lambda_7$.
  - Transport composed through the join to the other read link, other than via the joint state's $q_1q_2$ marginal (J4 $=\sqrt2$, which equals J3).
  - The non-covariant F1 letters in the F3a alphabet.
- **Bookkeeping and reference transport.** The tr frame in F3 is read record by record from the transported letters. Prop. 3(6) shows that no single linear pull-back does this. The letter supports are taken with $b=1$; with $b=0$, the transported identity $\Pi_S$ adds a $\lambda_8$ component.
- **The value.** $1/\alpha_\text{cand}$ is not affected. Nothing here changes the counts, the value or F-CT-1–F-CT-11.

## 8. Consequences for the gaps

- **CT-G01.** The direction left open in `operator-origin.md` §5 is tested.
  - A transport construction reproduces $M_L$ exactly. In it, the class-by-class attribution is derived from a frame and locality, and the qubit–qubit records inherit the QQ weight because their transported images are QQ records.
  - It is not an operator origin in CT-G01's sharpened sense ("an operator construction in which the end attribution is not an input"):
    - as a linear operator it either fails for every $D$ (Theorem 4) or needs the qubit-only selector (Prop. 3(5));
    - as a record-wise rule it needs the frame assignment and the credit J1 as inputs (Theorem 6, Theorem 5).
  - With the transport-compatible credits the value moves to 140.702 (J2) or 139.030 (J3). The gap stays open.
- **CT-G02.** "The end at which an exchange acts" now has an operator-level reading: the frame obtained by transporting along the read links. It no longer refers to the wording. Two things still have to be chosen: which exchange is read in which frame (of the 9 frame pairs, only {tr, orig-q} works), and whether the join is credited with the symmetry of the join ($\sqrt3$) or of the transport ($2/\sqrt3$, $\sqrt2$). The Letter's value corresponds to the first option in each.
