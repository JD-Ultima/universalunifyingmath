/-
  HowardArith.lean — core lemmas for the "1 ⊗ 1 = 2" reconstructions (claims TH-001, TH-002).
  Lean 4 core only (no Mathlib, no imports). Checked with Lean 4.35.0 (exit 0, output = HowardArith.lean.out):
      lean proofs/lean/howard/HowardArith.lean
  The ℤ part (everything except the sections marked "added 2026-10-30") was also checked with Lean 4.24.0.
  The added sections (generic one-sided no-go, ℚ instances, A6) were checked with 4.35.0 only.
  Write-up: formalization/howard/arithmetic.md (theorem numbers A0–A7 refer to it).

  Axioms (full audit at the end of the file): every ℤ / Fin 4 / generic theorem uses at most
  `propext` and `Quot.sound`. The ℚ instances also use `Classical.choice`, inherited from core's own
  `Rat` order and cast lemmas (e.g. `Rat.le_trans`, `Rat.natCast_add` depend on it).

  Notation. `op` is a candidate operation ⊗. H1 : op 1 1 = 2.  H2 : op 2 2 = 4.
  LD (left distributivity):  op a (b + c) = op a b + op a c.
  RD (right distributivity): op (a + b) c = op a c + op b c.
-/

namespace Howard

/-! ### A0. Standard arithmetic: 1 is the multiplicative identity (definitional in Peano). -/

theorem std_one_mul_one : (1 : Nat) * 1 = 1 := by decide
theorem std_peano (a : Nat) : a * 1 = a := Nat.mul_one a

/-- Any identity element is idempotent ("has no effect on itself"). -/
theorem identity_idempotent {α : Type} (op : α → α → α) (e : α)
    (he : ∀ b, op e b = b) : op e e = e := he e

/-- H1 forbids 1 from being a left or a right identity of ⊗. -/
theorem h1_not_left_identity (op : Int → Int → Int) (h1 : op 1 1 = 2) :
    ¬ (∀ b, op 1 b = b) := by
  intro h; have := h 1; omega

theorem h1_not_right_identity (op : Int → Int → Int) (h1 : op 1 1 = 2) :
    ¬ (∀ a, op a 1 = a) := by
  intro h; have := h 1; omega

/-! ### A1. Generic no-go: two-sided distributivity turns H1 into 2̄ ⊗ 2̄ = 8̄.
   Works for an arbitrary carrier with an arbitrary binary `add`; numerals are iterated sums of `one`
   (`one` is NOT assumed to be a ⊗-identity). -/

theorem two_two_eq_eight {α : Type} (add op : α → α → α) (one : α)
    (hl : ∀ a b c, op a (add b c) = add (op a b) (op a c))
    (hr : ∀ a b c, op (add a b) c = add (op a c) (op b c))
    (h1 : op one one = add one one) :
    op (add one one) (add one one)
      = add (add (add one one) (add one one)) (add (add one one) (add one one)) := by
  rw [hr, hl, h1]

/-- A1 (generic): H1 + H2 + LD + RD force 4̄ + 4̄ = 4̄ in the additive structure. -/
theorem nogo_generic {α : Type} (add op : α → α → α) (one : α)
    (hl : ∀ a b c, op a (add b c) = add (op a b) (op a c))
    (hr : ∀ a b c, op (add a b) c = add (op a c) (op b c))
    (h1 : op one one = add one one)
    (h2 : op (add one one) (add one one) = add (add one one) (add one one)) :
    add (add (add one one) (add one one)) (add (add one one) (add one one))
      = add (add one one) (add one one) := by
  rw [← two_two_eq_eight add op one hl hr h1, h2]

/-- A1 (ℤ): no operation on ℤ satisfies LD, RD, H1, H2. -/
theorem nogo_int (op : Int → Int → Int)
    (hl : ∀ a b c, op a (b + c) = op a b + op a c)
    (hr : ∀ a b c, op (a + b) c = op a c + op b c)
    (h1 : op 1 1 = 2) (h2 : op 2 2 = 4) : False := by
  have e1 : op 2 2 = op 1 2 + op 1 2 := by have := hr 1 1 2; simpa using this
  have e2 : op 1 2 = op 1 1 + op 1 1 := by have := hl 1 1 1; simpa using this
  omega

/-- A1 (ℤ/4ℤ loophole): in characteristic 4 the operation a⊗b = 2ab satisfies LD, RD, H1, H2. -/
def op4 (a b : Fin 4) : Fin 4 := 2 * a * b
theorem char4_h1 : op4 1 1 = 2 := by decide
theorem char4_h2 : op4 2 2 = 4 := by decide   -- (4 : Fin 4) = 0 = 8
theorem char4_ld : ∀ a b c : Fin 4, op4 a (b + c) = op4 a b + op4 a c := by decide
theorem char4_rd : ∀ a b c : Fin 4, op4 (a + b) c = op4 a c + op4 b c := by decide

/-! ### A2. One-sided distributivity + one instance of associativity (or commutativity). -/

theorem nogo_int_ld_assoc (op : Int → Int → Int)
    (hl : ∀ a b c, op a (b + c) = op a b + op a c)
    (hassoc : op (op 1 1) 1 = op 1 (op 1 1))
    (h1 : op 1 1 = 2) (h2 : op 2 2 = 4) : False := by
  have e1 : op 2 2 = op 2 1 + op 2 1 := by have := hl 2 1 1; simpa using this
  have e2 : op 1 2 = op 1 1 + op 1 1 := by have := hl 1 1 1; simpa using this
  rw [h1] at hassoc
  omega

theorem nogo_int_ld_comm (op : Int → Int → Int)
    (hl : ∀ a b c, op a (b + c) = op a b + op a c)
    (hcomm : op 2 1 = op 1 2)
    (h1 : op 1 1 = 2) (h2 : op 2 2 = 4) : False := by
  have e1 : op 2 2 = op 2 1 + op 2 1 := by have := hl 2 1 1; simpa using this
  have e2 : op 1 2 = op 1 1 + op 1 1 := by have := hl 1 1 1; simpa using this
  omega

theorem nogo_int_rd_assoc (op : Int → Int → Int)
    (hr : ∀ a b c, op (a + b) c = op a c + op b c)
    (hassoc : op (op 1 1) 1 = op 1 (op 1 1))
    (h1 : op 1 1 = 2) (h2 : op 2 2 = 4) : False := by
  have e1 : op 2 2 = op 1 2 + op 1 2 := by have := hr 1 1 2; simpa using this
  have e2 : op 2 1 = op 1 1 + op 1 1 := by have := hr 1 1 1; simpa using this
  rw [h1] at hassoc
  omega

theorem nogo_int_rd_comm (op : Int → Int → Int)
    (hr : ∀ a b c, op (a + b) c = op a c + op b c)
    (hcomm : op 2 1 = op 1 2)
    (h1 : op 1 1 = 2) (h2 : op 2 2 = 4) : False := by
  have e1 : op 2 2 = op 1 2 + op 1 2 := by have := hr 1 1 2; simpa using this
  have e2 : op 2 1 = op 1 1 + op 1 1 := by have := hr 1 1 1; simpa using this
  omega

/-! ### A1/A2, generic and over ℚ (added 2026-10-30).
   The one-sided cores hold in any magma `(α, add)`: each forces 8̄ = 4̄, where 2̄ = 1+1, 4̄ = 2̄+2̄,
   8̄ = 4̄+4̄. Hence they are inconsistent wherever 4̄ + 4̄ ≠ 4̄ (e.g. any cancellative monoid with
   4̄ ≠ 0: ℕ, ℤ, ℚ, ℝ, rings of characteristic not dividing 4). ℤ/4 (`char4_*`) shows this is sharp. -/

theorem nogo_generic_ld_assoc {α : Type} (add op : α → α → α) (one : α)
    (hl : ∀ a b c, op a (add b c) = add (op a b) (op a c))
    (hassoc : op (op one one) one = op one (op one one))
    (h1 : op one one = add one one)
    (h2 : op (add one one) (add one one) = add (add one one) (add one one)) :
    add (add (add one one) (add one one)) (add (add one one) (add one one))
      = add (add one one) (add one one) := by
  have e1 : op one (add one one) = add (add one one) (add one one) := by rw [hl, h1]
  have e2 : op (add one one) one = add (add one one) (add one one) := by
    rw [h1] at hassoc; rw [hassoc, e1]
  calc _ = add (op (add one one) one) (op (add one one) one) := by rw [e2]
    _ = op (add one one) (add one one) := (hl _ _ _).symm
    _ = _ := h2

theorem nogo_generic_ld_comm {α : Type} (add op : α → α → α) (one : α)
    (hl : ∀ a b c, op a (add b c) = add (op a b) (op a c))
    (hcomm : op (add one one) one = op one (add one one))
    (h1 : op one one = add one one)
    (h2 : op (add one one) (add one one) = add (add one one) (add one one)) :
    add (add (add one one) (add one one)) (add (add one one) (add one one))
      = add (add one one) (add one one) := by
  have e1 : op one (add one one) = add (add one one) (add one one) := by rw [hl, h1]
  have e2 : op (add one one) one = add (add one one) (add one one) := by rw [hcomm, e1]
  calc _ = add (op (add one one) one) (op (add one one) one) := by rw [e2]
    _ = op (add one one) (add one one) := (hl _ _ _).symm
    _ = _ := h2

theorem nogo_generic_rd_assoc {α : Type} (add op : α → α → α) (one : α)
    (hr : ∀ a b c, op (add a b) c = add (op a c) (op b c))
    (hassoc : op (op one one) one = op one (op one one))
    (h1 : op one one = add one one)
    (h2 : op (add one one) (add one one) = add (add one one) (add one one)) :
    add (add (add one one) (add one one)) (add (add one one) (add one one))
      = add (add one one) (add one one) := by
  have e1 : op (add one one) one = add (add one one) (add one one) := by rw [hr, h1]
  have e2 : op one (add one one) = add (add one one) (add one one) := by
    rw [h1] at hassoc; rw [← hassoc, e1]
  calc _ = add (op one (add one one)) (op one (add one one)) := by rw [e2]
    _ = op (add one one) (add one one) := (hr _ _ _).symm
    _ = _ := h2

theorem nogo_generic_rd_comm {α : Type} (add op : α → α → α) (one : α)
    (hr : ∀ a b c, op (add a b) c = add (op a c) (op b c))
    (hcomm : op (add one one) one = op one (add one one))
    (h1 : op one one = add one one)
    (h2 : op (add one one) (add one one) = add (add one one) (add one one)) :
    add (add (add one one) (add one one)) (add (add one one) (add one one))
      = add (add one one) (add one one) := by
  have e1 : op (add one one) one = add (add one one) (add one one) := by rw [hr, h1]
  have e2 : op one (add one one) = add (add one one) (add one one) := by rw [← hcomm, e1]
  calc _ = add (op one (add one one)) (op one (add one one)) := by rw [e2]
    _ = op (add one one) (add one one) := (hr _ _ _).symm
    _ = _ := h2

/-- Numerals in core `Rat`. -/
theorem rat_one_add_one : (1 + 1 : Rat) = 2 := by
  rw [show (1 : Rat) = ((1 : Nat) : Rat) from rfl, ← Rat.natCast_add]; rfl
theorem rat_two_add_two : (2 + 2 : Rat) = 4 := by
  rw [show (2 : Rat) = ((2 : Nat) : Rat) from rfl, ← Rat.natCast_add]; rfl
theorem rat_eight_ne_four : (4 + 4 : Rat) ≠ 4 := by
  rw [show (4 : Rat) = ((4 : Nat) : Rat) from rfl, ← Rat.natCast_add]
  intro h; exact absurd (Rat.natCast_inj.mp h) (by decide)
theorem rat_eight_ne_four' : ((1 + 1) + (1 + 1)) + ((1 + 1) + (1 + 1)) ≠ ((1 + 1) + (1 + 1) : Rat) := by
  rw [rat_one_add_one, rat_two_add_two]; exact rat_eight_ne_four

/-- A1 over ℚ. -/
theorem nogo_rat (op : Rat → Rat → Rat)
    (hl : ∀ a b c, op a (b + c) = op a b + op a c)
    (hr : ∀ a b c, op (a + b) c = op a c + op b c)
    (h1 : op 1 1 = 2) (h2 : op 2 2 = 4) : False :=
  rat_eight_ne_four' (nogo_generic (· + ·) op 1 hl hr
    (by show op 1 1 = 1 + 1; rw [rat_one_add_one]; exact h1)
    (by show op (1 + 1) (1 + 1) = (1 + 1) + (1 + 1); rw [rat_one_add_one, rat_two_add_two]; exact h2))

/-- A2 over ℚ (the four one-sided cores). -/
theorem nogo_rat_ld_assoc (op : Rat → Rat → Rat)
    (hl : ∀ a b c, op a (b + c) = op a b + op a c)
    (hassoc : op (op 1 1) 1 = op 1 (op 1 1))
    (h1 : op 1 1 = 2) (h2 : op 2 2 = 4) : False :=
  rat_eight_ne_four' (nogo_generic_ld_assoc (· + ·) op 1 hl hassoc
    (by show op 1 1 = 1 + 1; rw [rat_one_add_one]; exact h1)
    (by show op (1 + 1) (1 + 1) = (1 + 1) + (1 + 1); rw [rat_one_add_one, rat_two_add_two]; exact h2))

theorem nogo_rat_ld_comm (op : Rat → Rat → Rat)
    (hl : ∀ a b c, op a (b + c) = op a b + op a c)
    (hcomm : op 2 1 = op 1 2)
    (h1 : op 1 1 = 2) (h2 : op 2 2 = 4) : False :=
  rat_eight_ne_four' (nogo_generic_ld_comm (· + ·) op 1 hl
    (by show op (1 + 1) 1 = op 1 (1 + 1); rw [rat_one_add_one]; exact hcomm)
    (by show op 1 1 = 1 + 1; rw [rat_one_add_one]; exact h1)
    (by show op (1 + 1) (1 + 1) = (1 + 1) + (1 + 1); rw [rat_one_add_one, rat_two_add_two]; exact h2))

theorem nogo_rat_rd_assoc (op : Rat → Rat → Rat)
    (hr : ∀ a b c, op (a + b) c = op a c + op b c)
    (hassoc : op (op 1 1) 1 = op 1 (op 1 1))
    (h1 : op 1 1 = 2) (h2 : op 2 2 = 4) : False :=
  rat_eight_ne_four' (nogo_generic_rd_assoc (· + ·) op 1 hr hassoc
    (by show op 1 1 = 1 + 1; rw [rat_one_add_one]; exact h1)
    (by show op (1 + 1) (1 + 1) = (1 + 1) + (1 + 1); rw [rat_one_add_one, rat_two_add_two]; exact h2))

theorem nogo_rat_rd_comm (op : Rat → Rat → Rat)
    (hr : ∀ a b c, op (a + b) c = op a c + op b c)
    (hcomm : op 2 1 = op 1 2)
    (h1 : op 1 1 = 2) (h2 : op 2 2 = 4) : False :=
  rat_eight_ne_four' (nogo_generic_rd_comm (· + ·) op 1 hr
    (by show op (1 + 1) 1 = op 1 (1 + 1); rw [rat_one_add_one]; exact hcomm)
    (by show op 1 1 = 1 + 1; rw [rat_one_add_one]; exact h1)
    (by show op (1 + 1) (1 + 1) = (1 + 1) + (1 + 1); rw [rat_one_add_one, rat_two_add_two]; exact h2))

/-! ### A3. Classification on ℤ: every bi-additive ⊗ is a ⊗ b = a·b·(1 ⊗ 1). -/

theorem right_additive_linear (op : Int → Int → Int)
    (hr : ∀ a b c, op (a + b) c = op a c + op b c) :
    ∀ a c, op a c = a * op 1 c := by
  intro a c
  have h0 : op 0 c = 0 := by
    have := hr 0 0 c; simp at this; omega
  have hnat : ∀ n : Nat, op (n : Int) c = (n : Int) * op 1 c := by
    intro n
    induction n with
    | zero => simpa using h0
    | succ n ih =>
      have h := hr (n : Int) 1 c
      rw [show ((n + 1 : Nat) : Int) = (n : Int) + 1 by omega, h, ih, Int.add_mul, Int.one_mul]
  have hneg : ∀ x : Int, op (-x) c = - op x c := by
    intro x
    have := hr x (-x) c
    rw [Int.add_right_neg, h0] at this
    omega
  obtain ⟨n, rfl | rfl⟩ := Int.eq_nat_or_neg a
  · exact hnat n
  · rw [hneg, hnat, Int.neg_mul]

theorem left_additive_linear (op : Int → Int → Int)
    (hl : ∀ a b c, op a (b + c) = op a b + op a c) :
    ∀ a b, op a b = b * op a 1 := by
  intro a b
  exact right_additive_linear (fun x y => op y x) (fun x y z => hl z x y) b a

/-- A3: bi-additive ⊗ on ℤ is determined by k = 1 ⊗ 1:  a ⊗ b = a·(b·k). -/
theorem biadditive_classification (op : Int → Int → Int)
    (hl : ∀ a b c, op a (b + c) = op a b + op a c)
    (hr : ∀ a b c, op (a + b) c = op a c + op b c) :
    ∀ a b, op a b = a * (b * op 1 1) := by
  intro a b
  rw [right_additive_linear op hr a b, left_additive_linear op hl 1 b]

/-- A3 + H1: the only bi-additive ⊗ on ℤ with 1 ⊗ 1 = 2 is a ⊗ b = 2ab. -/
theorem forced_two_ab (op : Int → Int → Int)
    (hl : ∀ a b c, op a (b + c) = op a b + op a c)
    (hr : ∀ a b c, op (a + b) c = op a c + op b c)
    (h1 : op 1 1 = 2) : ∀ a b, op a b = a * (b * 2) := by
  intro a b; rw [biadditive_classification op hl hr, h1]

/-! ### A4. The forced operation ⊗₂ : a ⊗ b = 2ab on ℤ. -/

def op2 (a b : Int) : Int := a * (b * 2)

theorem op2_h1 : op2 1 1 = 2 := by decide
theorem op2_not_h2 : op2 2 2 = 8 := by decide
theorem op2_ld (a b c : Int) : op2 a (b + c) = op2 a b + op2 a c := by
  unfold op2; rw [Int.add_mul, Int.mul_add]
theorem op2_rd (a b c : Int) : op2 (a + b) c = op2 a c + op2 b c := by
  unfold op2; rw [Int.add_mul]
theorem op2_comm (a b : Int) : op2 a b = op2 b a := by
  unfold op2; simp only [Int.mul_left_comm]
theorem op2_assoc (a b c : Int) : op2 (op2 a b) c = op2 a (op2 b c) := by
  unfold op2; simp only [Int.mul_comm, Int.mul_left_comm]

/-- φ(x) = 2x is an injective ring homomorphism (ℤ, +, ⊗₂) → (ℤ, +, ·) with image 2ℤ. -/
theorem phi_add (a b : Int) : 2 * (a + b) = 2 * a + 2 * b := by omega
theorem phi_mul (a b : Int) : 2 * op2 a b = (2 * a) * (2 * b) := by
  unfold op2; simp only [Int.mul_comm, Int.mul_left_comm]
theorem phi_inj (a b : Int) (h : 2 * a = 2 * b) : a = b := by omega
theorem phi_image (y : Int) : (∃ x, 2 * x = y) ↔ y % 2 = 0 := by
  constructor
  · rintro ⟨x, rfl⟩; omega
  · intro h; exact ⟨y / 2, by omega⟩

/-- (ℤ, ⊗₂) has no identity element ... -/
theorem op2_no_identity : ¬ ∃ e : Int, ∀ b, op2 e b = b := by
  rintro ⟨e, he⟩
  have := he 1
  unfold op2 at this
  rw [Int.one_mul] at this
  omega

/-- ... and no idempotent other than 0. -/
theorem op2_idempotent (e : Int) (h : op2 e e = e) : e = 0 := by
  unfold op2 at h
  have h3 : e * (e * 2 - 1) = 0 := by rw [Int.mul_sub, Int.mul_one, h]; omega
  rcases Int.mul_eq_zero.mp h3 with h4 | h4
  · exact h4
  · omega

/-- ⊗₂ strictly increases positive integers: a, b ≥ 1 ⇒ a ⊗₂ b > a (and > b by symmetry). -/
theorem op2_increasing (a b : Int) (ha : 1 ≤ a) (hb : 1 ≤ b) : a < op2 a b := by
  unfold op2
  have := Int.mul_pos (show 0 < a by omega) (show 0 < b * 2 - 1 by omega)
  rw [Int.mul_sub, Int.mul_one] at this
  omega

/-- ⊗₂ dominates addition on positive integers ("balance": both inputs accounted for). -/
theorem op2_ge_add (a b : Int) (ha : 1 ≤ a) (hb : 1 ≤ b) : a + b ≤ op2 a b := by
  unfold op2
  have e : a * (b * 2) = a * b + a * b := by
    rw [Int.mul_comm b 2, Int.mul_left_comm, Int.two_mul]
  -- a + b ≤ 2ab  ⟸  (a-1)(b-1) ≥ 0  i.e.  ab ≥ a + b - 1
  have h2 := Int.mul_nonneg (show 0 ≤ a - 1 by omega) (show 0 ≤ b - 1 by omega)
  rw [Int.sub_mul, Int.mul_sub, Int.mul_sub, Int.one_mul, Int.mul_one, Int.one_mul] at h2
  omega

/-! ### Witnesses for the consistent axiom subsets (each fails exactly where stated). -/

/-- ⊗ = + : H1, H2, associative, commutative; not left-distributive. -/
def opAdd (a b : Int) : Int := a + b
theorem opAdd_h1 : opAdd 1 1 = 2 := by decide
theorem opAdd_h2 : opAdd 2 2 = 4 := by decide
theorem opAdd_assoc (a b c : Int) : opAdd (opAdd a b) c = opAdd a (opAdd b c) := by
  unfold opAdd; omega
theorem opAdd_comm (a b : Int) : opAdd a b = opAdd b a := by unfold opAdd; omega
theorem opAdd_not_ld : opAdd 1 (0 + 0) ≠ opAdd 1 0 + opAdd 1 0 := by decide
theorem opAdd_ld_fails_everywhere (a b c : Int) (h : opAdd a (b + c) = opAdd a b + opAdd a c) :
    a = 0 := by unfold opAdd at h; omega

/-- a ⊗ b = 2b : H1, H2, LD; not RD, not commutative, not associative. -/
def opLeft (_a b : Int) : Int := 2 * b
theorem opLeft_h1 : opLeft 1 1 = 2 := by decide
theorem opLeft_h2 : opLeft 2 2 = 4 := by decide
theorem opLeft_ld (a b c : Int) : opLeft a (b + c) = opLeft a b + opLeft a c := by
  unfold opLeft; omega
theorem opLeft_not_rd : opLeft (1 + 1) 1 ≠ opLeft 1 1 + opLeft 1 1 := by decide
theorem opLeft_not_comm : opLeft 2 1 ≠ opLeft 1 2 := by decide
theorem opLeft_not_assoc : opLeft (opLeft 1 1) 1 ≠ opLeft 1 (opLeft 1 1) := by decide

/-- Tropical reading (⊕ = max, ⊗ = +): ⊗ distributes over ⊕; H1 and H2 hold; but 1 ⊕ 1 = 1. -/
theorem tropical_ld (a b c : Int) : a + max b c = max (a + b) (a + c) := by omega
theorem tropical_idempotent_add : max (1 : Int) 1 = 1 := by decide

/-! ### A6. Archimedean no-go (added 2026-10-30): one distributive law + "multiplication makes more".
   Generic form: a carrier `α` with `add`, a relation `le` and a predicate `pos` such that positives are
   closed under `add`, `add` is monotone for `le`, every positive element has positive (n+1)-th parts,
   and `le` is Archimedean on positives. Then no `op` satisfies LD (needed only at positive right
   arguments) together with `le a (op a b)` for all positive a, b; mirror for RD and `le b (op a b)`.
   The generic theorems use no axioms. ℚ (core `Rat`) is instantiated below; ℝ satisfies the same four
   hypotheses (ordered field + Archimedean property), but ℝ is not in Lean core, so that instance is
   on paper. The positive-only hypotheses also cover the carriers ℚ_{>0}, ℝ_{>0}. -/

/-- `nsm add n x` is the (n+1)-fold sum `x + x + ⋯ + x` (left-nested). -/
def nsm {α : Type} (add : α → α → α) : Nat → α → α
  | 0, x => x
  | n + 1, x => add (nsm add n x) x

theorem nsm_pos {α : Type} (add : α → α → α) (pos : α → Prop)
    (hpa : ∀ x y, pos x → pos y → pos (add x y)) (x : α) (hx : pos x) :
    ∀ n, pos (nsm add n x)
  | 0 => hx
  | n + 1 => hpa _ _ (nsm_pos add pos hpa x hx n) hx

theorem nsm_mono {α : Type} (add : α → α → α) (le : α → α → Prop)
    (hmono : ∀ x y x' y', le x y → le x' y' → le (add x x') (add y y'))
    (x y : α) (h : le x y) : ∀ n, le (nsm add n x) (nsm add n y)
  | 0 => h
  | n + 1 => hmono _ _ _ _ (nsm_mono add le hmono x y h n) h

/-- LD (needed only at positive arguments) commutes with (n+1)-fold sums in the right slot. -/
theorem ld_nsm {α : Type} (add : α → α → α) (pos : α → Prop) (op : α → α → α)
    (hpa : ∀ x y, pos x → pos y → pos (add x y))
    (hl : ∀ a b c, pos b → pos c → op a (add b c) = add (op a b) (op a c))
    (a t : α) (ht : pos t) : ∀ n, op a (nsm add n t) = nsm add n (op a t)
  | 0 => rfl
  | n + 1 => by
      show op a (add (nsm add n t) t) = add (nsm add n (op a t)) (op a t)
      rw [hl a _ _ (nsm_pos add pos hpa t ht n) ht, ld_nsm add pos op hpa hl a t ht n]

theorem archimedean_ld_generic {α : Type} (add : α → α → α) (le : α → α → Prop) (pos : α → Prop)
    (hpa : ∀ x y, pos x → pos y → pos (add x y))
    (hmono : ∀ x y x' y', le x y → le x' y' → le (add x x') (add y y'))
    (hdiv : ∀ n x, pos x → ∃ y, pos y ∧ nsm add n y = x)
    (harch : ∀ x z, pos x → ∃ n, ¬ le (nsm add n x) z)
    (p : α) (hp : pos p) (op : α → α → α)
    (hl : ∀ a b c, pos b → pos c → op a (add b c) = add (op a b) (op a c))
    (hge : ∀ a b, pos a → pos b → le a (op a b)) : False := by
  obtain ⟨n, hn⟩ := harch p (op p p) hp
  obtain ⟨y, hy, hny⟩ := hdiv n p hp
  apply hn
  have e : op p p = nsm add n (op p y) := by
    have h := ld_nsm add pos op hpa hl p y hy n
    rw [hny] at h
    exact h
  rw [e]
  exact nsm_mono add le hmono _ _ (hge p y hp hy) n

theorem archimedean_rd_generic {α : Type} (add : α → α → α) (le : α → α → Prop) (pos : α → Prop)
    (hpa : ∀ x y, pos x → pos y → pos (add x y))
    (hmono : ∀ x y x' y', le x y → le x' y' → le (add x x') (add y y'))
    (hdiv : ∀ n x, pos x → ∃ y, pos y ∧ nsm add n y = x)
    (harch : ∀ x z, pos x → ∃ n, ¬ le (nsm add n x) z)
    (p : α) (hp : pos p) (op : α → α → α)
    (hr : ∀ a b c, pos a → pos b → op (add a b) c = add (op a c) (op b c))
    (hge : ∀ a b, pos a → pos b → le b (op a b)) : False :=
  archimedean_ld_generic add le pos hpa hmono hdiv harch p hp (fun x y => op y x)
    (fun a b c hb hc => hr b c a hb hc) (fun a b ha hb => hge b a hb ha)

/-! Instance ℚ (core `Rat`). -/

theorem rat_lt_of_lt_of_le {a b c : Rat} (h1 : a < b) (h2 : b ≤ c) : a < c :=
  Rat.not_le.mp fun h3 => Rat.not_le.mpr h1 (Rat.le_trans h2 h3)

theorem rat_natCast_one : ((1 : Nat) : Rat) = 1 := rfl

theorem rat_nsm (x : Rat) : ∀ n, nsm (· + ·) n x = ((n + 1 : Nat) : Rat) * x
  | 0 => by
      show x = ((0 + 1 : Nat) : Rat) * x
      rw [Nat.zero_add, rat_natCast_one, Rat.one_mul]
  | n + 1 => by
      show nsm (· + ·) n x + x = ((n + 1 + 1 : Nat) : Rat) * x
      rw [rat_nsm x n, Rat.natCast_add (n + 1) 1, Rat.add_mul, rat_natCast_one, Rat.one_mul]

theorem rat_pos_add (x y : Rat) (hx : 0 < x) (hy : 0 < y) : 0 < x + y := by
  have h : 0 + y < x + y := Rat.add_lt_add_right.mpr hx
  rw [Rat.zero_add] at h
  exact rat_lt_of_lt_of_le hy (Rat.le_of_lt h)

theorem rat_mono (x y x' y' : Rat) (h : x ≤ y) (h' : x' ≤ y') : x + x' ≤ y + y' :=
  Rat.le_trans (Rat.add_le_add_right.mpr h) (Rat.add_le_add_left.mpr h')

theorem rat_succ_pos (n : Nat) : (0 : Rat) < ((n + 1 : Nat) : Rat) :=
  Rat.natCast_pos.mpr (Nat.succ_pos n)

theorem rat_div (n : Nat) (x : Rat) (hx : 0 < x) :
    ∃ y : Rat, 0 < y ∧ nsm (· + ·) n y = x := by
  have hc := rat_succ_pos n
  have hc0 : ((n + 1 : Nat) : Rat) ≠ 0 := fun h => Rat.lt_irrefl (h ▸ hc)
  refine ⟨x / ((n + 1 : Nat) : Rat), ?_, ?_⟩
  · apply (Rat.lt_div_iff hc).mpr
    rw [Rat.zero_mul]; exact hx
  · rw [rat_nsm, Rat.mul_comm, Rat.div_mul_cancel hc0]

theorem rat_arch (x z : Rat) (hx : 0 < x) :
    ∃ n, ¬ nsm (· + ·) n x ≤ z := by
  refine ⟨(z / x).floor.toNat, ?_⟩
  rw [rat_nsm]
  apply Rat.not_le.mpr
  apply (Rat.div_lt_iff hx).mp
  apply rat_lt_of_lt_of_le (Rat.lt_floor_add_one (z / x))
  rw [← Rat.intCast_natCast]
  apply Rat.intCast_le_intCast.mpr
  omega

/-- A6 over ℚ: no ⊗ satisfies LD (at positive arguments) and a ⊗ b ≥ a for all positive a, b. -/
theorem archimedean_ld_rat (op : Rat → Rat → Rat)
    (hl : ∀ a b c, 0 < b → 0 < c → op a (b + c) = op a b + op a c)
    (hge : ∀ a b, 0 < a → 0 < b → a ≤ op a b) : False :=
  archimedean_ld_generic (· + ·) (· ≤ ·) (0 < ·) rat_pos_add rat_mono rat_div rat_arch
    1 (Rat.natCast_pos.mpr Nat.one_pos) op hl hge

theorem archimedean_rd_rat (op : Rat → Rat → Rat)
    (hr : ∀ a b c, 0 < a → 0 < b → op (a + b) c = op a c + op b c)
    (hge : ∀ a b, 0 < a → 0 < b → b ≤ op a b) : False :=
  archimedean_rd_generic (· + ·) (· ≤ ·) (0 < ·) rat_pos_add rat_mono rat_div rat_arch
    1 (Rat.natCast_pos.mpr Nat.one_pos) op hr hge

/-- Lattice core {LD, INC} over ℚ, with INC : a ⊗ b > max(a, b) for positive a, b. -/
theorem core_ld_inc_rat (op : Rat → Rat → Rat)
    (hl : ∀ a b c, op a (b + c) = op a b + op a c)
    (hinc : ∀ a b, 0 < a → 0 < b → a < op a b ∧ b < op a b) : False :=
  archimedean_ld_rat op (fun a b c _ _ => hl a b c)
    (fun a b ha hb => Rat.le_of_lt (hinc a b ha hb).1)

theorem core_rd_inc_rat (op : Rat → Rat → Rat)
    (hr : ∀ a b c, op (a + b) c = op a c + op b c)
    (hinc : ∀ a b, 0 < a → 0 < b → a < op a b ∧ b < op a b) : False :=
  archimedean_rd_rat op (fun a b c _ _ => hr a b c)
    (fun a b ha hb => Rat.le_of_lt (hinc a b ha hb).2)

/-! ### G1. Euler bound used for the Lynchpin (claim TH-020):
   a polyhedral 2-sphere whose faces are all pentagons and whose vertices all have degree ≥ 3
   has at least 12 faces; in particular six pentagonal faces are impossible. -/

theorem pentagon_faces_ge_12 (V E F : Int)
    (hface : 5 * F = 2 * E) (hdeg : 3 * V ≤ 2 * E) (heuler : V - E + F = 2) : 12 ≤ F := by
  omega

theorem no_six_pentagon_sphere (V E : Int)
    (hface : 5 * 6 = 2 * E) (hdeg : 3 * V ≤ 2 * E) (heuler : V - E + 6 = 2) : False := by
  omega

end Howard

-- Full axiom audit (added 2026-10-30): one `#print axioms` line for each of the file's
-- 80 theorems, in file order. Expected: at most propext / Quot.sound, plus Classical.choice
-- for the ℚ instances (inherited from core `Rat`); never sorryAx.
#print axioms Howard.std_one_mul_one
#print axioms Howard.std_peano
#print axioms Howard.identity_idempotent
#print axioms Howard.h1_not_left_identity
#print axioms Howard.h1_not_right_identity
#print axioms Howard.two_two_eq_eight
#print axioms Howard.nogo_generic
#print axioms Howard.nogo_int
#print axioms Howard.char4_h1
#print axioms Howard.char4_h2
#print axioms Howard.char4_ld
#print axioms Howard.char4_rd
#print axioms Howard.nogo_int_ld_assoc
#print axioms Howard.nogo_int_ld_comm
#print axioms Howard.nogo_int_rd_assoc
#print axioms Howard.nogo_int_rd_comm
#print axioms Howard.nogo_generic_ld_assoc
#print axioms Howard.nogo_generic_ld_comm
#print axioms Howard.nogo_generic_rd_assoc
#print axioms Howard.nogo_generic_rd_comm
#print axioms Howard.rat_one_add_one
#print axioms Howard.rat_two_add_two
#print axioms Howard.rat_eight_ne_four
#print axioms Howard.rat_eight_ne_four'
#print axioms Howard.nogo_rat
#print axioms Howard.nogo_rat_ld_assoc
#print axioms Howard.nogo_rat_ld_comm
#print axioms Howard.nogo_rat_rd_assoc
#print axioms Howard.nogo_rat_rd_comm
#print axioms Howard.right_additive_linear
#print axioms Howard.left_additive_linear
#print axioms Howard.biadditive_classification
#print axioms Howard.forced_two_ab
#print axioms Howard.op2_h1
#print axioms Howard.op2_not_h2
#print axioms Howard.op2_ld
#print axioms Howard.op2_rd
#print axioms Howard.op2_comm
#print axioms Howard.op2_assoc
#print axioms Howard.phi_add
#print axioms Howard.phi_mul
#print axioms Howard.phi_inj
#print axioms Howard.phi_image
#print axioms Howard.op2_no_identity
#print axioms Howard.op2_idempotent
#print axioms Howard.op2_increasing
#print axioms Howard.op2_ge_add
#print axioms Howard.opAdd_h1
#print axioms Howard.opAdd_h2
#print axioms Howard.opAdd_assoc
#print axioms Howard.opAdd_comm
#print axioms Howard.opAdd_not_ld
#print axioms Howard.opAdd_ld_fails_everywhere
#print axioms Howard.opLeft_h1
#print axioms Howard.opLeft_h2
#print axioms Howard.opLeft_ld
#print axioms Howard.opLeft_not_rd
#print axioms Howard.opLeft_not_comm
#print axioms Howard.opLeft_not_assoc
#print axioms Howard.tropical_ld
#print axioms Howard.tropical_idempotent_add
#print axioms Howard.nsm_pos
#print axioms Howard.nsm_mono
#print axioms Howard.ld_nsm
#print axioms Howard.archimedean_ld_generic
#print axioms Howard.archimedean_rd_generic
#print axioms Howard.rat_lt_of_lt_of_le
#print axioms Howard.rat_natCast_one
#print axioms Howard.rat_nsm
#print axioms Howard.rat_pos_add
#print axioms Howard.rat_mono
#print axioms Howard.rat_succ_pos
#print axioms Howard.rat_div
#print axioms Howard.rat_arch
#print axioms Howard.archimedean_ld_rat
#print axioms Howard.archimedean_rd_rat
#print axioms Howard.core_ld_inc_rat
#print axioms Howard.core_rd_inc_rat
#print axioms Howard.pentagon_faces_ge_12
#print axioms Howard.no_six_pentagon_sphere
