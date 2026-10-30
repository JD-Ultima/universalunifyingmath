# A uniform definition of the records? (gap CT-G09)

**Object.** The record space of the α Letter (Wissner-Gross, v1.0.0, "The records"; Supplement Cor. S-known, Cor. S-records, Thm S-further; `ChainAlpha.lean` `recs`, `commutatorRecs`): the 169 two-slot products and the three commutator records (F-CT-1), and the covariant version of F-CT-10.
**Question.** CT-G09: is there a uniform definition of "record", with no exception clause, that yields exactly these 172 records with the Letter's classes? Why is the difference of the two commutators kept, and not the sum?
**Code.** `computation/chain-theory/alpha/records_uniform.py` → `records_uniform.out`. Exact arithmetic: rational linear algebra (sympy `DomainMatrix`), and sympy for the Gell-Mann and Weyl checks. 1964 asserted checks with `--upstream`, 1914 without. It imports `symmetry.py`, which loads `records_independent.py`; neither is modified, and they are used only for a numerical cross-check of the two record spans. "out §n" is section n of the output.
**Date.** 2026-10-30. Analysis toward CT-G09 and CT-G11. The audited value is not changed.

**Notation.**
- Slots $Q_1,q_1,Q_2,q_2$. $V_s$ is the space of traceless operators at slot $s$ (the letters). $V_q=V_{q_1}\oplus V_{q_2}$ (dim 6), $V_Q=V_{Q_1}\oplus V_{Q_2}$ (dim 16).
- $P_2$: the span of the 169 two-slot products. $P_{\le2}$: operators with at most two letters, $1\oplus V_Q\oplus V_q\oplus P_2$ (dim 192).
- Own sign $\sigma=(-1)^{\hat N}$, with $\hat N$ the number of letters.
- Exchanges. $E$ is the plain slot swap: the Letter's, and `exchAssign` in Lean. $E'=E\circ T$ (F-CT-10). In general $f_{\varepsilon_Q,\varepsilon_q}=\varepsilon_Q^{\hat N_Q}\varepsilon_q^{\hat N_q}E'$; in particular $F_m:=f_{+,-}=(-1)^{\hat N_q}E'$.
- For an exchange $f$: $R_f:=P_2\oplus V_q^{f=-1}$. This is the Letter's construction carried out with $f$: products of both parities, plus the single-qubit combinations at the order sign.
- $C_L:=V_q^{E=-1}=\operatorname{span}\{\sigma_\lambda^{(q_1)}-\sigma_\lambda^{(q_2)}\}$ (the Letter's commutator records). $C:=V_q^{E'=-1}$ (F-CT-10's $D_x,D_z,S_y$). $C^+:=V_q^{E'=+1}$.
- $I:=P_{\le1}$: the image of the two read links' identification algebras, acting by $\xi_1\otimes1+1\otimes\xi_2$.
- $\mathrm{Lie}(G_\Omega)=\operatorname{span}\{h^{(1)}-(h^{\mathsf T})^{(2)}\}$ at the join $(Q_1,Q_2)$ and at the qubit pair $(q_1,q_2)$, with dim 11.
- Weights. $M_L(f)$ is `operator-origin.md` Prop. 2 with $f$ in place of $E'$. The "word rule" is the Letter's $X$, $Y$, extended to single-letter composites (`w2` in the script): $X$ reorders two letters at one slot, so a single letter of order sign $-1$ survives it; $Y$ moves a single qubit letter between the read qubits and does not act on the qutrit copies. The two laws agree on every $R_f$.

## Result in brief

1. **Source** (Sec. 1). A record is two letters at two distinct slots, and a single letter is excluded. The exception is the commutator of two anticommuting reads at a qubit, combined across the two reads as the difference only.
   - The Lean file builds these records as a separate list: at the far end, under `weylIsSign`, with a parity filter that the product list does not have.
   - The qutrit is excluded by an argument for invertible (Weyl) reads. The Letter's letters, the Gell-Mann matrices, contain five anticommuting pairs, and so does Lean's basis (Prop. 3).
2. **The qubit/qutrit half of the clause is uniform** (Thm 2). The order of the two insertions has a well-defined sign on exactly $D=1\oplus P_2\oplus V_q$.
   - At a qubit the anticommutator of two letters is a scalar.
   - At a qutrit every traceless single letter is also a symmetric composite.
   - Requiring a definite order sign (OD) therefore admits the qubit single letters and no qutrit ones, with no basis and no determinant.
3. **The parity half is not a sign rule** (Thm 4). On $D_0=D\cap1^\perp$, $R_f$ is "own sign $+1$, or survives $f$", three of the four $(\sigma,f)$ classes. No single sign test gives it, and no sum of isotypic components of $P_{\le2}$ has the Letter's classes (Prop. 5).
4. **One identification rule gives it** (Thms 6, 7). The rule (ID): a record is orthogonal to every identification whose action the exchange preserves. Then $R=D\cap\mathrm{Comp}_f^\perp$.
   - **Plain $E$:** $R=R_E$, the Letter's 172 records exactly, with no exception clause.
   - **Every covariant exchange:** $\mathrm{Comp}_f=1\oplus\mathrm{Lie}(G_\Omega)$, the infinitesimal identifications that fix the joint state. So $R=P_2\oplus C^+$, canonically.
   - This equals $R_f$, with the Letter's sectors, value and isotypic table, iff $\varepsilon_q=-1$. Then the exchange acts on the read qubits by the unitary $(\sigma_y\otimes\sigma_y)\,\mathrm{SWAP}$, the matching that undoes the conjugation in the Letter's readers remark (alpha.tex:79).
   - For $E'$ ($\varepsilon_q=+1$, F-CT-10) the Letter's commutator records are the matchings $C$, which (ID) removes. The rule's records sorted by $E'$ give $1/\text{value}=139.713$ ($M_L$) or $138.348$ (word rule).
5. **Difference vs sum** (Sec. 6). The suggested reason, "the sum is the action of the identification algebra and is not a record, as the global phase is not", is exactly (ID).
   - It selects the Letter's difference in the plain frame and in the covariant frames with $\varepsilon_q=-1$.
   - In the $E'$ frame the kept combination *is* the matching algebra of the qubit pair, and the reason selects the other three operators.
   - It is consistent with the Letter's definitions but is not stated there.
6. **The transpose-identity collapse points the other way** (Prop. 8). Under $X\mapsto X|\Omega\rangle$ the antisymmetric qubit–qubit products collapse onto the sums $C^+$. The commutator records $C$ collapse to 0.
7. **Declared family** (Props. 9, 10). The family has 864 members: 12 composite spaces (8, four of them with and without the identity) × 24 conditions × 3 frames. Every member with the Letter's classes is its frame's $R_f$. In the $E'$ frame the only such member is $D\cap(I^{E'})^\perp$, the Letter's rule restated, and $D\cap I^{E'}$ is not a Lie algebra. The images of $P_2$ under 18 natural superoperators never have $C$, $C_L$ or $C^+$ as their single-letter part.
8. **What the clause adds** (Cor. 11). It adds the 3-dimensional space $V_q^{f=-1}$, in operator form $R_f=D_0\cap\ker\bigl((1+f)[\hat N=1]\bigr)$. Its qubit/qutrit half is OD. Its parity half is (ID) when the exchange is automorphic on the qubits, and "the order-definite matchings are records" for $E'$.

**Consequence.** CT-G09's closure criterion is met conditionally, on two readings the Letter does not state and, in a covariant frame, one choice of exchange.
- (OD) replaces the Weyl/determinant argument.
- (ID) extends the dropped "multiple of the identity" to the exchange-compatible identifications.
- (R): the read qubits are identified by $\sigma_y$, as in the Letter's readers remark, not by the transpose.

Under these, records are defined by one rule and the value is unchanged. Under $E'$ alone, the Letter's covariant records need a second clause.

---

## 1. What the source says

Commit `a28db89`. The passages below are checked at these lines by out §0, which also checks the sha256 of `ChainAlpha.lean` against the pin.

**Record and exclusion.**
- `alpha.tex:109`: "A *measurement record* carries the two insertions as two letters. Two letters at the same slot act on the same end, where their product is a multiple of the identity plus one traceless observable, a single insertion. So the two letters of a record sit at different slots [Fig. 1(b)], except in the commutator records below."
- `supplement.tex:300` (Cor. S-records): "A product record is a product of one letter at each of two of the four slots."
- `supplement.tex:303` (proof): "A single letter is one insertion, and three or four letters are more than two. Two letters at one slot multiply to a multiple of the identity plus one traceless observable, a single insertion, which is not a record of two insertions, except where anticommuting letters make their commutator one."
- `supplement.tex:35` (Table S-inputs): "a record: two letters, the two insertions", entered as a "step of the model". `supplement.tex:573` (glossary): "product of letters at two slots, or a commutator record".

**The exception.**
- `alpha.tex:62` (the known): "Two reads fix their product jointly, and exchanging them returns the product times the phase by which they fail to commute. An exchange done twice is no exchange, so that phase must be a sign, and two anticommuting reads thus fix one thing more than either alone, their product, with their order as its sign."
- `alpha.tex:109`: "At a qubit two reads in different bases anticommute (the Weyl phase there is $-1$) … At the qutrit no two reads anticommute, since $AB=-BA$ would give $\det(AB)=-\det(BA)$ in odd dimension, and no such record arises. The exchange carries the commutator at one read's qubit to the same commutator at the other's, and of the two combinations only their difference reads the order sign $-1$ under it, so there are three commutator records".
- `supplement.tex:78`, `:81` (Cor. S-known): "Two reads that anticommute therefore fix a further record, their product, with their order as its sign". In the proof: "$AB=\tfrac12\{A,B\}+\tfrac12[A,B]$, the first part fixed whatever their order and the second by it … For two invertible reads with $AB=\lambda BA$ … a sign is possible only at even $d$". The reads there are the clock and shift letters of the mutually unbiased bases (`:70`).
- `supplement.tex:310`, `:313` (Thm S-further): "The commutator records are the three differences $[s_\mu,s_\nu]_{q_1}-[s_\mu,s_\nu]_{q_2}$ … on their span it reads $+1$ at the sum and $-1$ at the difference. The definition of the known makes the order sign the sign of the record, so the record is the difference, and the sum, which reads $+1$, is not a commutator record."
- `supplement.tex:364`: the three span "a copy of $\Lambda^2V_3$".

**Both eigenvalues kept for products.** `alpha.tex:109` says the exchange sorts the space "into its two eigenspaces, keeping both, and does not quotient it". `supplement.tex:306` adds that "both eigenvalues are kept, the eigenvalue labelling the sector".

**Own sign.** `alpha.tex:113`: "$+1$ for a product … and the order sign $-1$ for a commutator record". `supplement.tex:575` says the same.

**Lean.**
- `ChainAlpha.lean:216–220`: the letters are "the real traceless basis of an end (the paper's Gell-Mann and Pauli matrices, up to factors of i)"; "a commutator record is a single letter at the far end of one read minus its image at the far end of the other".
- `:300–302`: `commutatorCandidates` takes the commutators at slots 1 and 3 (the far ends) and returns their sum and difference.
- `:324–329`: `commutatorRecs` is `if weylIsSign dB then …` and keeps a candidate `if eqZ (parity fg.1) (orderSign dB)`.
- `:280`: `weylIsSign` is the clock–shift phase. `:295`: `orderSign`.
- `:350`: `recs` is `odd ++ even ++ fixedP ++ commutatorRecs dA dB`, so the products are kept at both parities.
- `:638`, `:665`: `weyl_sign` and `record_count`.

So the Lean file encodes the exception as a separate list with three restrictions the products do not have: the far end only, the Weyl sign of the clock and shift letters, and parity equal to the order sign. Its own qutrit basis (`tracelessBasis 3`) is never tested for anticommutation.

**Readers remark** (used in Thm 7).
- `alpha.tex:79`: "the far end can take that basis as its own only if a matching undoes the conjugation, as $\sigma_y$ does at a qubit and no unitary does at a qutrit".
- `supplement.tex:170`: "on traceless letters $A^{\mathsf T}=-\sigma_yA\sigma_y$", "so this route is a second one and not the one the Letter takes".

## 2. The record spaces

**Proposition 1 (Proven; out §1, §2).**
1. $R_E$ is the span of the Letter's 172 records, and $R_{E'}$ is the span of F-CT-10's covariant records. Both identities are checked numerically against `records_independent.py` and `symmetry.py`.
2. $C_L$, $C$ and $C^+$ are three 3-dimensional subspaces of $V_q$, with $\dim(C_L\cap C)=2$, $\dim(C_L\cap C^+)=1$ and $C\cap C^+=0$. The Letter's commutator records $[b_\mu,b_\nu]^{(q_1)}-[b_\mu,b_\nu]^{(q_2)}$ span $C_L$. Each is one letter at $q_1$ minus the same letter at $q_2$ (up to $2i$).
3. $\mathrm{Lie}(G_\Omega)$ is the $E'$-odd part of $V_Q\oplus V_q$, and $C$ is its qubit part. These are the infinitesimal matchings $V\otimes\bar V$ of the qubit pair. In the plain frame, $C_L\cap\mathrm{Lie}(G_\Omega)$ is spanned by the $x$ and $z$ differences, and $C_L^+\cap\mathrm{Lie}(G_\Omega)$ by the $y$ sum.
4. For each $f\in\{E,E',F_m,f_{-,+},\sigma E'\}$, $R_f$ has the classes $(79,48,42,3)$ and the numerator $79+48\sqrt2+42\sqrt3+3\sqrt6$ under both weight laws.
5. $R_f$ depends on $f$ only through $\varepsilon_q$: $R_{E'}=R_{f_{-,+}}=P_2\oplus C$ and $R_{F_m}=R_{\sigma E'}=P_2\oplus C^+$. With $R_E=P_2\oplus C_L$ these are three distinct spaces with one value.

*Proof.*
- (3): $\tfrac{d}{d\theta}\bigl(e^{i\theta h}\otimes\overline{e^{i\theta h}}\bigr)=i(h\otimes1-1\otimes\bar h)$ and $\bar h=h^{\mathsf T}$ for Hermitian $h$. Also $E'(h^{(q_1)})=h^{\mathsf T(q_2)}$. In the plain frame, $\sigma_x^{\mathsf T}=\sigma_x$, $\sigma_z^{\mathsf T}=\sigma_z$ and $\sigma_y^{\mathsf T}=-\sigma_y$.
- (4), (5): the classes are cell dimensions (out §1), and the cell dimensions are the same in every frame.

∎

## 3. Order-definite composites

Let $V=\bigoplus_sV_s$ (22 letters), $m:V\otimes V\to B(H)$ the product, $\tau$ the flip, and $m_\pm=m(\mathrm{Sym}^2V)$, $m(\Lambda^2V)$.

**Theorem 2 (order-definiteness; Proven, out §3).**
1. At an end of dimension $d$, the traceless part of $\{A,B\}$ vanishes for all letters $A,B$ iff $d=2$. Always $m(\Lambda^2V_d)=V_d$. Hence $m(\mathrm{Sym}^2V_d)\cap m(\Lambda^2V_d)$ is $0$ at $d=2$ and $V_d$ at $d\ge3$ (computed at $d=2,3,4$).
2. $m(V\otimes V)=P_{\le2}$; $m_+=1\oplus P_2\oplus V_Q$ (dim 186); $m_-=V_Q\oplus V_q$ (dim 22). So $A:=m_+\cap m_-=V_Q$.
3. $\tau$ induces a well-defined map on $m(V\otimes V)/B$ (that is, $m(x)\in B\Rightarrow m(\tau x)\in B$) iff $B\supseteq A$. On $m(V\otimes V)/A$ it is the involution with eigenspaces $m_\pm/A$.
4. With the trace pairing, $D:=A^\perp\cap P_{\le2}=1\oplus P_2\oplus V_q$ (dim 176). The induced order sign is $+1$ on $1\oplus P_2$ and $-1$ on $V_q$; it equals the own sign $\sigma$ there.

*Proof.*
1. For $d=2$, $\{a\cdot\sigma,b\cdot\sigma\}=2(a\cdot b)\,1$. For $d\ge3$, take $A=B=E_{01}+E_{10}$: $A^2=E_{00}+E_{11}$ has a nonzero traceless part. The map $\mathrm{Sym}^2V_d\to V_d$ is $\mathfrak{su}(d)$-equivariant and $V_d$ is irreducible, so the map is onto. $[V_d,V_d]=V_d$ because $\mathfrak{sl}(d)$ is perfect.
2. For $s\ne t$, $a^{(s)}b^{(t)}=b^{(t)}a^{(s)}$. So $\Lambda^2$ maps to $0$ there and $\mathrm{Sym}^2$ maps onto $V_s\otimes V_t$. Then add the same-slot terms from (1).
3. Necessity: for $y\in A$, write $y=m(x_+)=m(x_-)$. Then $x=x_+-x_-$ has $m(x)=0$ and $m(\tau x)=2y$. Sufficiency: if $m(x)\in A$, then $m(x_+)=m(x)-m(x_-)\in m_-\cap m_+=A$, so $m(\tau x)=2m(x_+)-m(x)\in A$.
4. This follows from (2), and from the fact that $A$ is a sum of slot blocks orthogonal to the rest.

∎

**Reading (OD).** A record has a definite order sign: it lies in $D$, and, since it is traceless, in $D_0=P_2\oplus V_q$ (dim 175). This is the Supplement's "the first part fixed whatever their order and the second by it" applied to spans instead of to a pair of letters. It needs no basis and no invertibility. Under (OD) the qutrit copies carry no single-letter record because every traceless single letter there is also a symmetric composite. The read qubits carry six, all with order sign $-1$. At $d=4$, (OD) admits none, as Lean's `weylIsSign 4 = false` does.

**Proposition 3 (anticommuting letters; Proven, out §3).**
1. Among the eight Gell-Mann matrices exactly five pairs anticommute: $(1,2)$, $(1,3)$, $(2,3)$, $(4,5)$, $(6,7)$. Their commutators span $\{\lambda_1,\lambda_2,\lambda_3,\lambda_8\}$. The same holds for Lean's real basis.
2. $\lambda_1,\dots,\lambda_7$ are singular. No two of the Weyl letters $X^aZ^b$ anticommute at $d=3$, and neither do any two letters of a generic orthogonal Hermitian basis (one exhibited exactly, by a rational Cayley transform).
3. At $d=2$ every orthogonal traceless Hermitian basis anticommutes pairwise.

*Consequences.*
- The Letter's sentence "at the qutrit no two reads anticommute" is true of invertible reads, and the determinant argument covers only those. Its records, however, are built from Gell-Mann letters.
- The rule taken literally, with those letters, adds four $E$-odd qutrit records, giving $79+48\sqrt2+46\sqrt3+3\sqrt6$ and $1/\text{value}=132.977$ (out §9).
- At the qutrit the count is basis-dependent (0 to 5 pairs). At the qubit it is 3 in every basis.
- (OD) gives the Letter's exclusion without a basis.

## 4. Sign rules and isotypic sums

**Theorem 4 (sign rules; Proven, out §4).** Let $f\in\{E,E',F_m\}$.
1. On $D_0$ the commuting involutions $\sigma$ and $f$ give four classes, $(\sigma,f)=(+,+),(+,-),(-,+),(-,-)$, of dimensions $90, 79, 3, 3$.
2. The seven rules with at most one sign test (no test, or $g=\pm1$ for one $g\in\{\sigma,f,\sigma f\}$) give dimensions $\{175,169,6,93,82,93,82\}$.
3. Exactly two unions of classes have dimension 172: $\{(+,+),(+,-),(-,-)\}=R_f$ and $\{(+,+),(+,-),(-,+)\}$.
4. Neither union is cut out by one sign test. $R_f=\{\sigma=+1\}\cup\{\sigma f=+1\}$, that is, "products, and the survivors among the order-sign composites".
5. On $P_{\le2,0}$ (no OD) the classes are $(90,79,11,11)$, and no union has dimension 172. With the identity included ($D$ or $P_{\le2}$), no sign rule has dimension 172.

*Proof.* Cell dimensions, then arithmetic. The sets cut out by at most one sign test are the cosets of subgroups of $\mathbb Z_2^2$ of index $\le2$, so they have 2 or 4 elements. $R_f$ has 3. ∎

The Letter's own wording matches the disjunction: products keep "both eigenvalues" (`supplement.tex:306`), while the commutator combinations are kept only at "the order sign" (`:313`).

**Proposition 5 (isotypic sums; Proven, out §4).**
1. Under $SU(3)\times SU(2)\times\langle E'\rangle$, $P_{\le2}$ decomposes as follows (quadratic Casimirs, exact; $10$ and $\overline{10}$ are separated by $\Lambda^2\mathbf8=\mathbf8+\mathbf{10}+\overline{\mathbf{10}}$):

   | $+$ | $(1,1)$ 3, $(1,3)$ 3, $(1,5)$ 5, $(8,1)$ 16, $(8,3)$ 48, $(27,1)$ 27 |
   |---|---|
   | $-$ | $(1,3)$ 6, $(8,1)$ 16, $(8,3)$ 48, $(10,1)$ 10, $(\overline{10},1)$ 10 |

2. $(1,3)_-$ is the odd qubit–qubit products together with $C$. $(8,1)_\pm$ is QQ together with the qutrit single letters.
3. Of the $2^{11}$ sums of isotypic components, exactly one has dimension 172: all but $(10,1)_-$ and $(\overline{10},1)_-$. Its classes are wrong: it keeps the identity and all 22 single letters.
4. $R_{E'}$ is not a sum of isotypic components: it meets $(8,1)_-$ in 8 of 16 dimensions.
5. $R_{F_m}$ under $\langle F_m\rangle$ has the same table as $R_{E'}$ under $\langle E'\rangle$ (F-CT-10 item 4).

## 5. Identification rules and the covariant exchanges

Let $\mathrm{ad}_X=[X,\cdot\,]$ on $P_{\le2}$. Define:
- $I^f$: the $f$-even part of $I$;
- $\mathrm{Comp}_f=\{X\in I: f\,\mathrm{ad}_X f=\mathrm{ad}_X\}$, the identifications whose action the exchange preserves;
- $\mathrm{stab}_\Omega=\{X\in I: X|\Omega\rangle\in\mathbb C|\Omega\rangle\}$.

**Theorem 6 (identification rules; Proven, out §5).**
1. $\mathrm{stab}_\Omega=1\oplus\mathrm{Lie}(G_\Omega)$ (dim 12).
2. $\mathrm{Comp}_E=I^E$, the $E$-even part: the diagonal action $\xi\mapsto\xi^{(L_1)}+\xi^{(L_2)}$ of the read links' identification algebra, with the global phase. It is not $\mathrm{stab}_\Omega$. For each of the four covariant exchanges, $\mathrm{Comp}_f=\mathrm{stab}_\Omega$.
3. (OD)+(ID): $R^{\rm id}_f:=D\cap\mathrm{Comp}_f^\perp$.
   - $R^{\rm id}_E=R_E$: **the Letter's 172 records, exactly.**
   - For every covariant $f$: $R^{\rm id}_f=D\cap\mathrm{stab}_\Omega^\perp=P_2\oplus C^+$.
   - Hence $R^{\rm id}_f=R_f$ iff $f=E$ or $\varepsilon_q=-1$.
4. In every frame, $R_f=D\cap(I^f)^\perp$. The part that $D$ sees, $D\cap I^f=1\oplus V_q^{f=+1}$, is a Lie algebra iff $f=E$ or $\varepsilon_q=-1$, and then it equals $D\cap\mathrm{Comp}_f$. For $\varepsilon_q=+1$ ($E'$, $f_{-,+}$) it is $1\oplus C^+$, and $[C^+,C^+]=C$.
5. Without (OD): $P_{\le2}\cap\mathrm{stab}_\Omega^\perp=P_2\oplus$ ($E'$-even single letters at both ends), of dimension 180.

*Proof.*
1. Let $X=\sum_sx_s^{(s)}$. By the transpose identity, $X|\Omega\rangle=(x_{Q_1}+x_{Q_2}^{\mathsf T})^{(Q_1)}|\Omega\rangle+(x_{q_1}+x_{q_2}^{\mathsf T})^{(q_1)}|\Omega\rangle$. The two terms lie in independent tensor factors, so $X|\Omega\rangle\in\mathbb C|\Omega\rangle$ forces each bracket to be a multiple of $1$. Since the letters are traceless, it is $0$.
2. Automorphism $E$: $E\,\mathrm{ad}_XE=\mathrm{ad}_{E(X)}$, so the condition is $E(X)-X\in\mathbb C1$. Anti-automorphism $E'$: $E'([X,Y])=-[E'X,E'Y]$, so $E'\mathrm{ad}_XE'=-\mathrm{ad}_{E'X}$ and the condition is $E'X=-X$ modulo $1$. In both cases this is a condition on each slot block separately. The other three exchanges act slot by slot (Thm 7), and $\mathrm{Comp}_f$ is computed for all five (out §5).
3. Follows from (2) and Thm 2(4).
4. Bracket computation.

∎

**Theorem 7 (the covariant exchanges; Proven, out §5).** Call an exchange *letter-wise* if it is $\mathrm{swap}\circ\bigotimes_s\varphi_s$, with unital linear maps $\varphi_s$ from slot $s$ to its partner slot that carry letters to letters and preserve the trace pairing.
1. The letter-wise involutive exchanges that commute with $G_\Omega$ are exactly the four $f_{\varepsilon_Q,\varepsilon_q}$. By Schur, the $G_\Omega$-intertwiners $V_{q_1}\to V_{q_2}$ and $V_{Q_1}\to V_{Q_2}$ are lines, spanned by $b\mapsto b^{\mathsf T}$. So $\varphi_s=c_s\,T$. Preserving the pairing gives $c_s=\pm1$, and being an involution gives $c_s=c_{s'}$ on partner slots. $E$ is an automorphism and is not covariant.
2. $E'$ is an anti-automorphism of $B(H)$. The other three are neither automorphisms nor anti-automorphisms (out §5 tests all products of two letters).
3. All four give the Letter's sectors $79/48/42/3$, its value, and the same isotypic table (Prop. 1(4), Prop. 5(5)).
4. For $\varepsilon_q=-1$ the exchange acts on the read qubits as $\mathrm{Ad}_W$ with $W=(\sigma_y\otimes\sigma_y)\,\mathrm{SWAP}_{q}$, a unitary. This is checked on all 16 basis operators of $P_{\le2}$ without qutrit letters, using $-b^{\mathsf T}=\sigma_yb\sigma_y$. At the qutrit copies no unitary intertwines $U$ with $\bar U$, so the transpose is forced there (conjugation.md Remark (b); `supplement.tex:170`).

So covariance fixes the exchange only up to the sign $\varepsilon_q$ (and $\varepsilon_Q$, which no record or value sees). The commutator-type span is $C$ for $\varepsilon_q=+1$ and $C^+$ for $\varepsilon_q=-1$.

## 6. Why the difference, not the sum

- **What the Letter gives as its reason.** It is `supplement.tex:313`: the record's exchange eigenvalue must equal its order sign. That is survival used as a condition of membership. The Letter does not use it for products: their non-survivors stay as records of weight 1, and "both eigenvalues are kept" (`:306`).
- **Plain frame.** The sums $\sigma_\lambda^{(q_1)}+\sigma_\lambda^{(q_2)}$ are $V_q\cap\mathrm{Comp}_E$. They are the diagonal action of the read links' identification algebra, whose $\mathfrak u(1)$ is the global phase. Removing the exchange-compatible identifications, as the global phase is removed, removes exactly the sums (Thm 6(3)). So the suggested operator reason holds verbatim for the Letter's own exchange. Two caveats:
  - $E$ is not $G_\Omega$-covariant (F-CT-10).
  - In $G_\Omega$ terms the Letter's kept span $C_L$ contains two matchings (the $x$ and $z$ differences) and the Letter's discarded span contains one (the $y$ sum) (Prop. 1(3)).
- **Covariant frames.** The identifications that every covariant exchange preserves are the matchings and the phase.
  - With $\varepsilon_q=-1$: the kept combination is $F_m$-odd, $C^+=V_q\cap\mathrm{stab}_\Omega^\perp$, and the discarded "sum" is $C$, the matching algebra $V\otimes\bar V$ of the qubit pair. The reason holds covariantly.
  - With $E'$: the kept combination $C$ is that matching algebra, and the reason would remove it.
- **Does it follow from the Letter's other definitions?** No.
  - The Letter defines the identification algebra and its action (`alpha.tex:81`). It drops the identity part of a same-slot product (`:109`). It states the $\sigma_y$ identification at the qubit only in a remark that it does not use (`supplement.tex:170`).
  - It states neither (ID) nor which covariant exchange it means. Its own exchange is $E$, and "nothing lines up the bases of a shared qudit's copies" (`alpha.tex:60`).
  - The reason is consistent with the definitions, not implied by them.

**Proposition 8 (collapse under the transpose identity; Proven, out §6).** Let $\kappa(X)=X|\Omega\rangle$ (rank 36 on $P_{\le2}$).
1. $\ker\kappa\cap I=\mathrm{Lie}(G_\Omega)$.
2. $\kappa(C)=0$. $\kappa(C^+)=\kappa\bigl((qq)^{E'=-1}\bigr)$, of dimension 3, and it lies inside $\kappa(P_2)$.
3. In the plain frame, $\dim\kappa(C_L)=1$ and $\dim\kappa(C_L^+)=2$.
4. $P_2+(D_0\cap\ker\kappa)=D_0$.

So a condition read through $X|\Omega\rangle$ cannot keep $P_2$ and drop $C^+$. The relation "the commutator records are the antisymmetric part of a two-slot insertion that collapses under the transpose identity" fails: $(a\otimes b^{\mathsf T}-b\otimes a^{\mathsf T})|\omega_2\rangle=([a,b]\otimes1)|\omega_2\rangle$ identifies the antisymmetric qubit–qubit records with the sums.

## 7. The declared family

**Proposition 9 (family F; Computed, out §7).**

*Declared family.* A member is $W=B\cap K$, in a frame $f\in\{E,E',F_m\}$.
- The composite space $B$ is one of 8: $P_2$; $m(V\otimes V)$; $m(\mathrm{Sym}^2V)$; $m(\Lambda^2V)$; $D$; $D_+$; $D_-$; one letter per read ($QQ+qq+$cross-mixed, 121). Where $B$ contains the identity it is taken with and without it.
- The condition $K$ is one of 24:
  - none; traceless;
  - $f$-even; $f$-odd; $\sigma=\pm1$; survivor; non-survivor;
  - $\perp I$, $\perp I^f$, $\perp\mathrm{Comp}_f$, $\perp\mathrm{stab}_\Omega$, $\perp\mathrm{Lie}(G_\Omega)$;
  - $\perp$ $G_\Omega$-invariants (the kernel of the twirl); $G_\Omega$-invariant (its image);
  - $\perp$ read-link matchings $t(B)^{(Q_i)}-B^{(q_i)}$; $\perp$ all pair matchings;
  - $\ker\kappa$; $\langle\Omega|X|\Omega\rangle=0$; $X\Omega\in\mathbb C\Omega$; $[X,P_\Omega]=0$;
  - $X|L_1L_2\rangle=0$; $[X,P_{L_i}]=0$; $[X,P_3]=[X,P_2]=0$.

The occupation projectors $[\hat N=k]$ and the slot-class projectors are not conditions, because selecting by letter count or end type is the exception clause itself.

*Results.* 864 members.
- 21 members have dimension 172. Every member with the Letter's classes is its own frame's $R_f$.
- $E$: $D$ or $D_0$, with $\perp I^E$ or $\perp\mathrm{Comp}_E$.
- $E'$: only $D$ or $D_0$ with $\perp I^{E'}$, where $D\cap I^{E'}$ is not a Lie algebra (Thm 6(4)).
- $F_m$: $D$ or $D_0$ with $\perp I^{F_m}$, $\perp\mathrm{Comp}_{F_m}$, $\perp\mathrm{stab}_\Omega$; also $D_0\cap\mathrm{Lie}(G_\Omega)^\perp$.
- The other dimension-172 members are $P_2\oplus C^+$ seen from $E$ (classes $(79,48,42,1)+2$; $1/\text{value}$ 138.809 or 137.908) or from $E'$ (classes $(79,48,42,0)+3$; 139.713 or 138.348).

**Proposition 10 (one-step images; Computed, out §8).** Take $\Phi\in\{L_P,R_P,\mathrm{Ad}_P\}$ for $P\in\{P_3,P_2,P_\Omega,P_{L_1},P_{L_2},P_{L_1}P_{L_2}\}$, 18 maps in all. None gives a single-letter part of $\Phi(P_2)$ equal to $C$, $C_L$ or $C^+$.
- Multiplication by the qubit-pair projector gives $1\oplus C^+$, the sums.
- Multiplication by the join gives $1\oplus(8,1)_+$.
- The read links give subspaces that are not graded by $E'$.

## 8. Variants and their values

Values are given under $M_L$ of the stated frame and under the word rule (out §9).

| record space | frame | dim | $M_L$ numerator | 1/value | word numerator | 1/value |
|---|---|---|---|---|---|---|
| Letter: $P_2\oplus C_L$ | $E$ | 172 | $79+48\sqrt2+42\sqrt3+3\sqrt6$ | 137.036 | same | 137.036 |
| F-CT-10: $P_2\oplus C$ | $E'$ | 172 | same | 137.036 | same | 137.036 |
| $D\cap\mathrm{stab}_\Omega^\perp=P_2\oplus C^+$ | $F_m$ | 172 | same | 137.036 | same | 137.036 |
| $P_2\oplus C^+$ | $E'$ | 172 | $82+48\sqrt2+42\sqrt3$ | 139.713 | $79+48\sqrt2+45\sqrt3$ | 138.348 |
| $P_2$ (no commutator records) | any | 169 | $79+48\sqrt2+42\sqrt3$ | 141.621 | same | 141.621 |
| $D_0$ (sums and differences) | $E'$ | 175 | $82+48\sqrt2+42\sqrt3+3\sqrt6$ | 135.248 | $79+48\sqrt2+45\sqrt3+3\sqrt6$ | 133.969 |
| $P_2\oplus\mathrm{Lie}(G_\Omega)$ | $E'$ | 180 | $79+48\sqrt2+50\sqrt3+3\sqrt6$ | 129.152 | same | 129.152 |
| $P_{\le2,0}$ | $E'$ | 191 | $90+48\sqrt2+50\sqrt3+3\sqrt6$ | 123.510 | $79+48\sqrt2+61\sqrt3+3\sqrt6$ | 119.683 |
| Gell-Mann literal: $R_E$ + 4 qutrit differences | $E$ | 176 | $79+48\sqrt2+46\sqrt3+3\sqrt6$ | 132.977 | same | 132.977 |

- Table S5's "sums and differences kept" (132.739) gives the sums $\sqrt6$. The Letter's survival rule gives them $\sqrt3$ at own sign $-1$ (133.969, as the word rule here and the Lean rewording of F-CT-4 D do), or $\sqrt2$ at own sign $+1$ ($79+51\sqrt2+42\sqrt3+3\sqrt6$, 134.522). It gives $\sqrt6$ at neither (out §9).
- The two values for $P_2\oplus C^+$ under $E'$ coincide numerically with Table S5's "no own sign" (139.713) and "comm survive X only" (138.348).

## 9. What the exception clause adds

**Corollary 11 (Proven).**
1. In each frame the clause adds exactly the 3-dimensional space $V_q^{f=-1}$: $C_L$ for $E$, $C$ for $E'$, $C^+$ for $\varepsilon_q=-1$. In operator form $R_f=D_0\cap\ker\bigl((1+f)[\hat N=1]\bigr)$ (out §5). The occupation projector $[\hat N=1]$, which restricts the condition to single letters, is the clause.
2. It has two halves.
   - (a) Single-letter composites are admitted at the qubit and not at the qutrit. This is (OD), a uniform rule (Thm 2).
   - (b) Only the exchange-odd combination is kept. This is (ID) when the exchange acts on the read qubits as an automorphism ($E$, or $\varepsilon_q=-1$) (Thm 6).
3. For $E'$, (b) is "the order-definite matchings are records". The covariant two-clause form is $R_{E'}=P_2\oplus(D_0\cap\mathrm{stab}_\Omega)=D_0\cap(P_2+\mathrm{Lie}(G_\Omega))$.
4. Without (b), the record space is $D_0$ (175 records; 135.248 or 133.969). Without (a), it is $P_2\oplus\mathrm{Lie}(G_\Omega)$ (180; 129.152) or $P_{\le2,0}$ (191).

## 10. Status

- **Proven.** Props. 1, 3, 5, 8; Theorems 2, 4, 6, 7; Corollary 11. These are finite linear algebra and Schur's lemma. Every finite statement is checked exactly in `records_uniform.out`, except Prop. 1(1), which is a numerical cross-check.
- **Computed.** Props. 9 and 10, and the variants table. They hold within the declared families.
- **Not stated in the Letter (inputs of the resolution).**
  - (OD): a record has a definite order sign. It reads `supplement.tex:81` on spans and replaces the determinant argument, which covers invertible reads only.
  - (ID): an identification whose action the exchange preserves is not a record. The Letter applies this only to the identity part of a same-slot product.
  - (R), covariant frames only: the exchange identifies the read qubits by $\sigma_y$. This is the readers remark's matching (`alpha.tex:79`), which the Letter says it does not use.
- **Open.**
  - An axiom-level reason for (ID) and (R).
  - Or, for $E'$, a reason why the order-definite matchings are records.

## 11. Consequences for the gaps

- **CT-G09. Resolved conditionally.**
  - In the Letter's frame, (OD)+(ID) define the 172 records with no exception clause: $R_E=D\cap\mathrm{Comp}_E^\perp$.
  - Covariantly, the same rule gives the canonical $G_\Omega$-invariant space $P_2\oplus C^+$. It carries the Letter's sectors and value iff the exchange identifies the read qubits by $\sigma_y$ ($\varepsilon_q=-1$).
  - Under $E'$ (CT-G11 as written) the Letter's records need the second clause of Cor. 11(3). The uniform rule then gives 139.713 or 138.348.
- **CT-G11.**
  - The covariant letter-wise exchanges are the four $f_{\varepsilon_Q,\varepsilon_q}$, not $E'$ alone. All four give the same sectors, value and isotypic table, so CT-G11 is resolved without change of value under each.
  - The commutator-type records are $D_x,D_z,S_y$ (F-CT-10) for $\varepsilon_q=+1$, and $S_x,D_y,S_z$ for $\varepsilon_q=-1$.
  - CT-G09 prefers $\varepsilon_q=-1$. $E'$ is the only one that is an (anti-)automorphism of $B(H)$.
- **CT-G12.** The commutant structure (22, with 9 real trace parameters) is the same for $R_{F_m}$ (Prop. 5(5)).
- **Committed statements this sharpens.**
  - `gaps/chain-theory-alpha.md` (CT-G09 row): "Only the difference is kept, not the sum" is true for the Letter's $E$. Covariantly the invariant statement is "only the exchange-odd combination". Under $E'$ that combination contains the $y$ sum.
  - `formalization/chain-theory/alpha/audit.md:52` and `:174`, and `computation/chain-theory/alpha/algebraic.py:23`, count $\binom b2$ commutator records at even $d_q$. At $d_q=4$ that is 105 records. But the commutators of the 15 letters at one slot span only 15 dimensions, so 105 orthogonal records do not exist. Lean's `commutatorRecs` gives none at $d_q=4$ (`weylIsSign 4 = false`), and so does (OD).
  - F-CT-7's $(2,4)$ entry, 216.04, becomes 331.40 with no commutator records, or 307.91 with 15. Its conclusion, that $(3,2)$ is the only pair within 1%, is unaffected.
