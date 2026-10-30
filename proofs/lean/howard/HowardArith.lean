/-
  HowardArith.lean — core lemmas for the "1 ⊗ 1 = 2" reconstructions (claims TH-001, TH-002).
  Lean 4 core only (no Mathlib, no imports). Checked with Lean 4.24.0 and 4.35.0 (identical output):
      lean proofs/lean/howard/HowardArith.lean
  Write-up: formalization/howard/arithmetic.md (theorem numbers A1–A9 refer to it).

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

/-! ### A5. Witnesses for the consistent axiom subsets (each fails exactly where stated). -/

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

-- Axiom audit: every theorem should depend on at most propext / Quot.sound / Classical.choice.
#print axioms Howard.nogo_int
#print axioms Howard.nogo_generic
#print axioms Howard.nogo_int_ld_assoc
#print axioms Howard.nogo_int_rd_comm
#print axioms Howard.biadditive_classification
#print axioms Howard.forced_two_ab
#print axioms Howard.op2_assoc
#print axioms Howard.op2_idempotent
#print axioms Howard.op2_ge_add
#print axioms Howard.char4_ld
#print axioms Howard.no_six_pentagon_sphere
