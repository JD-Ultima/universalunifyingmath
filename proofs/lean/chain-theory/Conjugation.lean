/-! # Far-end conjugation and the direct-sum trace (alpha audit, gaps CT-G03, CT-G04)

Standalone file. Lean 4.35.0, core library only (no Mathlib, no imports, no upstream file). Run with
`run_conj.sh`; the output is `Conjugation.out`. Companion of `formalization/chain-theory/alpha/conjugation.md`
and `computation/chain-theory/alpha/conjugation.py`.

What is checked (section numbers below):

* §1 (the acting scalar of a reading). A reading of a link with ends of dimensions c (the shared, read end) and
  e (the far end) is a trace condition α tr X_c + β tr X_e = 0 and a variance sign ε of the far factor of the
  link (ε = 1: X_e acts; ε = -1: conj X_e acts). On the centre (i x 1_c, i y 1_e) the condition is
  α c x + β e y = 0, whose solutions are spanned by (β e, -α c); that element acts on the link as
  i (x + ε y) = i (β e - ε α c). Readings A, A*, B of the audit (F-CT-3) are instances. The scalar under the
  relabellings ε = (ε1, ε2) of the Letter's trace (α = β = 1) depends on ε only through ε1 ε2.
* §2 (matching closure). The matching phase generator has traces (i m, -i m), m = min(c, e) ≥ 1, at the two
  ends. It satisfies the trace condition iff α = β; in ε-labels with the Letter's trace, iff ε1 = ε2.
* §3 (Theorem 1's selection). An independent enumeration of the 76 configurations of Supplement Sec. S5
  (42 admitted by the read cap, 23 readable), as in `config_enumeration.py`. For every acting predicate P,
  the end rule, the two-point definition and "P on every link" select exactly the qubit configuration
  (2R 2R 2U at a shared qubit) if P 2 2 and the object (2R 2R 2U at a shared qutrit) if P 3 2. Hence the
  selection for every fixed (α, β, ε), for readings A, A*, B, C, for the four relabellings, and for the
  natural trace spaces F of `conjugation.py` Sec. 2 (their coefficients α_F, β_F entered as the polynomials
  that the Python file computes and checks against representation matrices).

What is NOT checked here: ℂ, Lie algebras, representations, the trace calculus of the constructions
(Λ², Sym², End, ⊗, ⊕) and the identification of the readings with subalgebras; these are in
`conjugation.py` (exact sympy) and in the proofs of `conjugation.md`. The coefficients α_F, β_F are inputs. -/

namespace Conjugation

/-! ## §1. The acting scalar -/

/-- The scalar (divided by i) by which the spanning element (β e, -α c) of the centre of the reading
{α tr X_c + β tr X_e = 0} acts on the link, the far factor having variance ε. -/
def actScalar (α β ε c e : Int) : Int := β * e - ε * (α * c)

/-- Acting central charge of a reading: with no condition (α = β = 0) the whole centre u(1) ⊕ u(1) acts;
otherwise the centre is a line and acts iff its scalar is nonzero. -/
def acting (α β ε : Int) (c e : Nat) : Bool :=
  if α = 0 ∧ β = 0 then true else decide (actScalar α β ε c e ≠ 0)

/-- Reading A: the Letter's trace on H_c ⊕ H_e, the far end acting by X_e. -/
theorem actScalar_A (c e : Int) : actScalar 1 1 1 c e = e - c := by unfold actScalar; omega

/-- Reading A*: the same trace, the far end acting by conj X_e. -/
theorem actScalar_Astar (c e : Int) : actScalar 1 1 (-1) c e = e + c := by unfold actScalar; omega

theorem actScalar_Astar_ne (c e : Int) (hc : 0 < c) (he : 0 < e) : actScalar 1 1 (-1) c e ≠ 0 := by
  rw [actScalar_Astar]; omega

/-- Reading B: the trace on the link H_c ⊗ H_e, i.e. e tr X_c + c tr X_e = 0. -/
theorem actScalar_B (c e : Int) : actScalar e c 1 c e = 0 := by
  unfold actScalar; rw [Int.mul_comm e c]; omega

/-- Reading B with the conjugate far end: the trace on H_c ⊗ conj H_e is e tr X_c - c tr X_e, and the
centre again acts trivially. -/
theorem actScalar_B_conj (c e : Int) : actScalar e (-c) (-1) c e = 0 := by
  unfold actScalar; rw [Int.neg_mul, Int.mul_comm e c]; omega

/-- The relabelling ε = (ε1, ε2) of the Letter's h_W: its centre (e, -c) acts on the link as ε1 e - ε2 c. -/
def actEps (ε1 ε2 c e : Int) : Int := ε1 * e - ε2 * c

theorem actEps_zero_iff (ε1 ε2 c e : Int) (h1 : ε1 = 1 ∨ ε1 = -1) (h2 : ε2 = 1 ∨ ε2 = -1)
    (hc : 0 < c) (he : 0 < e) : actEps ε1 ε2 c e = 0 ↔ ε1 = ε2 ∧ c = e := by
  unfold actEps
  constructor
  · intro h
    rcases h1 with rfl | rfl <;> rcases h2 with rfl | rfl <;> exact ⟨by omega, by omega⟩
  · intro ⟨h3, h4⟩
    subst h3; subst h4; omega

/-- The acting predicate of the relabelled readings depends on ε only through the product ε1 ε2. -/
theorem actEps_product (ε1 ε2 c e : Int) (h1 : ε1 = 1 ∨ ε1 = -1) (h2 : ε2 = 1 ∨ ε2 = -1)
    (hc : 0 < c) (he : 0 < e) :
    (actEps ε1 ε2 c e = 0 ↔ actEps 1 (ε1 * ε2) c e = 0) := by
  unfold actEps
  rcases h1 with rfl | rfl <;> rcases h2 with rfl | rfl <;>
    exact ⟨fun h => by omega, fun h => by omega⟩

/-! ## §2. Matching closure -/

/-- The matching phase generator (i 1_m at the end of dimension m = min(c,e), -i on the support of the
other end, 0 off it) has traces (m, -m) (divided by i). It satisfies α tr X1 + β tr X2 = 0 iff α = β. -/
theorem mc_iff (α β m : Int) (hm : 0 < m) : α * m + β * (-m) = 0 ↔ α = β := by
  constructor
  · intro h
    rw [Int.mul_neg] at h
    have h' : (α - β) * m = 0 := by rw [Int.sub_mul]; omega
    rcases Int.mul_eq_zero.mp h' with h1 | h1 <;> omega
  · intro h
    subst h
    rw [Int.mul_neg]; omega

/-- In ε-labels the matching phase generator has label traces (ε1 m, -ε2 m) (tr conj X = -tr X), so it lies
in the Letter's h_W (α = β = 1) iff ε1 = ε2: exactly the relabellings with ε1 ε2 = 1. -/
theorem mc_eps (ε1 ε2 m : Int) (hm : 0 < m) : ε1 * m + ε2 * (-m) = 0 ↔ ε1 = ε2 := mc_iff ε1 ε2 m hm

/-! ## §3. Theorem 1's selection over the 76 configurations -/

/-- An arm: far dimension and read flag. A configuration: centre dimension and its arms. -/
abbrev Arm := Nat × Bool
abbrev Cfg := Nat × List Arm

/-- Count vectors (a, b, c, d) with a + b + c + d = L, in lexicographic order. -/
def countVecs (L : Nat) : List (Nat × Nat × Nat × Nat) :=
  (List.range (L + 1)).flatMap fun a =>
    (List.range (L + 1 - a)).flatMap fun b =>
      (List.range (L + 1 - a - b)).map fun c => (a, b, c, L - a - b - c)

/-- The multiset with those counts of 2R, 2U, 3R, 3U, as a sorted list. -/
def armsOf (v : Nat × Nat × Nat × Nat) : List Arm :=
  List.replicate v.1 (2, true) ++ List.replicate v.2.1 (2, false) ++
  List.replicate v.2.2.1 (3, true) ++ List.replicate v.2.2.2 (3, false)

/-- Supplement Sec. S5: centre 2 or 3, between 2 and c + 1 links, at least one read. -/
def configs : List Cfg :=
  [2, 3].flatMap fun c => (List.range' 2 c).flatMap fun L =>
    (countVecs L).filterMap fun v => if (armsOf v).any (·.2) then some (c, armsOf v) else none

/-- The read cap: a read arm needs c e ≤ 6. -/
def admitted : List Cfg := configs.filter fun k => k.2.all fun a => !a.2 || decide (k.1 * a.1 ≤ 6)
/-- At least two reads. -/
def readableCfgs : List Cfg := admitted.filter fun k => decide (2 ≤ (k.2.filter (·.2)).length)

theorem config_counts : configs.length = 76 ∧ admitted.length = 42 ∧ readableCfgs.length = 23 := by decide

/-- (a) the end rule: a read arm is read at the centre, which is the larger end. -/
def endRule (k : Cfg) : Bool := k.2.all fun a => !a.2 || decide (a.1 ≤ k.1)
/-- (b) the two-point definition: three arms of one type, two of them read. -/
def twoPoint (k : Cfg) : Bool :=
  k.2.length == 3 && (k.2.filter (·.2)).length == 2 && k.2.all (fun a => a.1 == (k.2.headD (0, false)).1)
/-- (c) an acting central charge on every link, for an acting predicate P (centre dim, far dim). -/
def actingAll (P : Nat → Nat → Bool) (k : Cfg) : Bool := k.2.all fun a => P k.1 a.1

def select (P : Nat → Nat → Bool) : List Cfg :=
  readableCfgs.filter fun k => actingAll P k && (endRule k && twoPoint k)

def qubitCfg : Cfg := (2, [(2, true), (2, true), (2, false)])
def objectCfg : Cfg := (3, [(2, true), (2, true), (2, false)])

/-- (a) and (b) alone leave exactly two configurations. -/
theorem endRule_twoPoint : readableCfgs.filter (fun k => endRule k && twoPoint k) = [qubitCfg, objectCfg] := by
  decide

/-- The selection for every acting predicate. -/
theorem select_general (P : Nat → Nat → Bool) :
    select P = (if P 2 2 then [qubitCfg] else []) ++ (if P 3 2 then [objectCfg] else []) := by
  have h : select P = [qubitCfg, objectCfg].filter (actingAll P) := by
    unfold select; rw [← endRule_twoPoint, List.filter_filter]
  rw [h]
  cases h22 : P 2 2 <;> cases h32 : P 3 2 <;>
    simp [List.filter, actingAll, qubitCfg, objectCfg, h22, h32]

/-- Theorem 1 selects the object alone iff the acting test fails at (2,2) and holds at (3,2). -/
theorem select_object_iff (P : Nat → Nat → Bool) : select P = [objectCfg] ↔ (P 2 2 = false ∧ P 3 2 = true) := by
  rw [select_general]
  cases P 2 2 <;> cases P 3 2 <;> decide

/-- For a fixed reading (α, β) with the canonical link (ε = 1): the selection, in closed form. -/
theorem select_fixed (α β : Int) (h : ¬ (α = 0 ∧ β = 0)) :
    select (acting α β 1) = (if α = β then [] else [qubitCfg]) ++ (if 2 * β = 3 * α then [] else [objectCfg]) := by
  rw [select_general]
  have e22 : acting α β 1 2 2 = !decide (α = β) := by
    unfold acting actScalar; simp only [h, ↓reduceIte]
    by_cases hab : α = β <;> simp [hab] <;> omega
  have e32 : acting α β 1 3 2 = !decide (2 * β = 3 * α) := by
    unfold acting actScalar; simp only [h, ↓reduceIte]
    by_cases hab : 2 * β = 3 * α <;> simp [hab] <;> omega
  rw [e22, e32]
  by_cases h1 : α = β <;> by_cases h2 : 2 * β = 3 * α <;> simp [h1, h2]

/-- Matching closure (α = β) with a nonzero condition selects the object alone. -/
theorem select_mc (α : Int) (h : α ≠ 0) : select (acting α α 1) = [objectCfg] := by
  rw [select_fixed α α (by omega)]
  have : ¬ (2 * α = 3 * α) := by omega
  simp [this]

/-- Conversely, for a fixed nonzero (α, β), the object alone forces α = β. -/
theorem object_forces_mc (α β : Int) (h : ¬ (α = 0 ∧ β = 0)) (hs : select (acting α β 1) = [objectCfg]) :
    α = β := by
  rw [select_fixed α β h] at hs
  by_cases h1 : α = β
  · exact h1
  · simp [h1] at hs
    exact absurd hs.1 (by decide)

/-- The four readings of F-CT-3. B is type-dependent: its coefficients are (e, c). -/
def actA (c e : Nat) : Bool := acting 1 1 1 c e
def actAstar (c e : Nat) : Bool := acting 1 1 (-1) c e
def actB (c e : Nat) : Bool := acting e c 1 c e
def actC (c e : Nat) : Bool := acting 0 0 1 c e

theorem readings :
    select actA = [objectCfg] ∧ select actAstar = [qubitCfg, objectCfg] ∧
    select actB = [] ∧ select actC = [qubitCfg, objectCfg] := by decide

/-- The four relabellings of the Letter's trace. -/
def actRel (ε1 ε2 : Int) (c e : Nat) : Bool := decide (actEps ε1 ε2 c e ≠ 0)

theorem relabellings :
    select (actRel 1 1) = [objectCfg] ∧ select (actRel (-1) (-1)) = [objectCfg] ∧
    select (actRel 1 (-1)) = [qubitCfg, objectCfg] ∧ select (actRel (-1) 1) = [qubitCfg, objectCfg] := by decide

/-- Natural trace spaces F (centre end H1 = H_c, far end H2 = H_e): coefficients (α_F, β_F) as computed and
checked against representation matrices in `conjugation.py` Sec. 2. -/
def coefF : List (String × (Int → Int → Int × Int)) :=
  [("W = H1+H2", fun _ _ => (1, 1)),
   ("Lambda2 W", fun c e => (c + e - 1, c + e - 1)),
   ("Sym2 W", fun c e => (c + e + 1, c + e + 1)),
   ("W (x) W", fun c e => (2 * (c + e), 2 * (c + e))),
   ("W' = H1+conj H2", fun _ _ => (1, -1)),
   ("Lambda2 W'", fun c e => (c + e - 1, -(c + e - 1))),
   ("End W", fun _ _ => (0, 0)),
   ("End H1 + End H2", fun _ _ => (0, 0)),
   ("L = H1 (x) H2", fun c e => (e, c)),
   ("Lambda2 L", fun c e => (e * (c * e - 1), c * (c * e - 1))),
   ("Sym2 L", fun c e => (e * (c * e + 1), c * (c * e + 1))),
   ("End L", fun _ _ => (0, 0)),
   ("L' = H1 (x) conj H2", fun c e => (e, -c)),
   ("H1", fun _ _ => (1, 0)),
   ("H2", fun _ _ => (0, 1)),
   ("H1+H1+H2", fun _ _ => (2, 1))]

def actF (f : Int → Int → Int × Int) (c e : Nat) : Bool := acting (f c e).1 (f c e).2 1 c e

/-- Number of configurations selected for each natural trace space: 1 = the object alone, 2 = both, 0 = none.
Matching closure at (2,2) and (3,2) (α_F = β_F there) holds for the first four rows, for the End rows
(vacuously) and, at (2,2) only, for the L rows. -/
theorem natural_spaces :
    coefF.map (fun p => (p.1, (select (actF p.2)).length)) =
      [("W = H1+H2", 1), ("Lambda2 W", 1), ("Sym2 W", 1), ("W (x) W", 1),
       ("W' = H1+conj H2", 2), ("Lambda2 W'", 2), ("End W", 2), ("End H1 + End H2", 2),
       ("L = H1 (x) H2", 0), ("Lambda2 L", 0), ("Sym2 L", 0), ("End L", 2),
       ("L' = H1 (x) conj H2", 2), ("H1", 2), ("H2", 2), ("H1+H1+H2", 2)] := by decide

/-- Matching closure per natural space at (2,2) and at (3,2). -/
def mcF (f : Int → Int → Int × Int) (c e : Int) : Bool := decide ((f c e).1 = (f c e).2)

theorem natural_spaces_mc :
    coefF.map (fun p => (p.1, mcF p.2 2 2, mcF p.2 3 2)) =
      [("W = H1+H2", true, true), ("Lambda2 W", true, true), ("Sym2 W", true, true), ("W (x) W", true, true),
       ("W' = H1+conj H2", false, false), ("Lambda2 W'", false, false), ("End W", true, true),
       ("End H1 + End H2", true, true), ("L = H1 (x) H2", true, false), ("Lambda2 L", true, false),
       ("Sym2 L", true, false), ("End L", true, true), ("L' = H1 (x) conj H2", false, false),
       ("H1", false, false), ("H2", false, false), ("H1+H1+H2", false, false)] := by decide

end Conjugation

#print axioms Conjugation.actScalar_A
#print axioms Conjugation.actScalar_Astar_ne
#print axioms Conjugation.actScalar_B
#print axioms Conjugation.actScalar_B_conj
#print axioms Conjugation.actEps_zero_iff
#print axioms Conjugation.actEps_product
#print axioms Conjugation.mc_iff
#print axioms Conjugation.mc_eps
#print axioms Conjugation.config_counts
#print axioms Conjugation.endRule_twoPoint
#print axioms Conjugation.select_general
#print axioms Conjugation.select_object_iff
#print axioms Conjugation.select_fixed
#print axioms Conjugation.select_mc
#print axioms Conjugation.object_forces_mc
#print axioms Conjugation.readings
#print axioms Conjugation.relabellings
#print axioms Conjugation.natural_spaces
#print axioms Conjugation.natural_spaces_mc
