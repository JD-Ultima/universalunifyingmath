/-
  GStueckelberg.lean — algebraic core of Theorem G1 (formalization/geometric-unity/G-full-el-fermions.md, G.2):
  the repaired tilted gauge group τ'(H) acts on 𝒢 = H ⋉ N freely with orbit space N, the quotient map being
  π'(ε,ϖ) = ε(ϖ − ε⁻¹d_{A₀}ε)ε⁻¹; every τ'-invariant function on 𝒢 is a function of π' alone, namely its value
  in the "unitary gauge" ε = 1.  Same abstract model as GUTilt.lean (definitions copied verbatim so that this
  file is self-contained).  Lean 4 core only (no Mathlib, no imports).
  Check:  lean proofs/lean/geometric-unity/GStueckelberg.lean
-/
namespace GStueckelberg

set_option linter.unusedSectionVars false

class Grp (H : Type) where
  mul : H → H → H
  one : H
  inv : H → H
  mul_assoc : ∀ a b d : H, mul (mul a b) d = mul a (mul b d)
  one_mul : ∀ a : H, mul one a = a
  mul_one : ∀ a : H, mul a one = a
  inv_mul : ∀ a : H, mul (inv a) a = one
  mul_inv : ∀ a : H, mul a (inv a) = one

class AbG (N : Type) where
  add : N → N → N
  zero : N
  neg : N → N
  add_assoc : ∀ x y z : N, add (add x y) z = add x (add y z)
  add_comm : ∀ x y : N, add x y = add y x
  zero_add : ∀ x : N, add zero x = x
  neg_add : ∀ x : N, add (neg x) x = zero

class Cocycle (H : Type) (N : outParam Type) [Grp H] [AbG N] where
  act : N → H → N
  c : H → N
  act_add : ∀ (x y : N) (h : H), act (AbG.add x y) h = AbG.add (act x h) (act y h)
  act_mul : ∀ (x : N) (h₁ h₂ : H), act x (Grp.mul h₁ h₂) = act (act x h₁) h₂
  act_one : ∀ x : N, act x Grp.one = x
  cocycle : ∀ h₁ h₂ : H, c (Grp.mul h₁ h₂) = AbG.add (act (c h₁) h₂) (c h₂)

open Grp AbG Cocycle

instance {N : Type} [AbG N] : Std.Associative (AbG.add (N := N)) := ⟨AbG.add_assoc⟩
instance {N : Type} [AbG N] : Std.Commutative (AbG.add (N := N)) := ⟨AbG.add_comm⟩

section General
variable {H N : Type} [Grp H] [AbG N] [Cocycle H N]

/-! ### Elementary consequences of the axioms -/

theorem add_zero' (x : N) : add x zero = x := by rw [add_comm, zero_add]
theorem add_neg (x : N) : add x (neg x) = zero := by rw [add_comm, neg_add]

theorem neg_unique {x y : N} (h : add x y = zero) : y = neg x := by
  have e : add (neg x) (add x y) = y := by rw [← add_assoc, neg_add, zero_add]
  rw [h, add_zero'] at e
  exact e.symm

theorem add_left_cancel' {a x y : N} (h : add a x = add a y) : x = y := by
  have hx : add (neg a) (add a x) = x := by rw [← add_assoc, neg_add, zero_add]
  have hy : add (neg a) (add a y) = y := by rw [← add_assoc, neg_add, zero_add]
  rw [← hx, ← hy, h]

theorem neg_neg' (x : N) : neg (neg x) = x := (neg_unique (neg_add x)).symm

theorem neg_add_distrib (x y : N) : neg (add x y) = add (neg x) (neg y) := by
  refine (neg_unique ?_).symm
  rw [show add (add x y) (add (neg x) (neg y)) = add (add x (neg x)) (add y (neg y)) by ac_rfl,
    add_neg, add_neg, zero_add]

theorem act_zero (h : H) : act (zero : N) h = zero := by
  have e : add (act (zero : N) h) (act zero h) = add (act zero h) zero := by
    rw [← act_add, zero_add, add_zero']
  exact add_left_cancel' e

theorem act_neg (x : N) (h : H) : act (neg x) h = neg (act x h) := by
  apply neg_unique
  rw [← act_add, add_neg, act_zero]

theorem act_inv_cancel (x : N) (h : H) : act (act x h) (inv h) = x := by
  rw [← act_mul, mul_inv, act_one]

theorem act_inv_cancel' (x : N) (h : H) : act (act x (inv h)) h = x := by
  rw [← act_mul, inv_mul, act_one]

theorem c_one : c (one : H) = (zero : N) := by
  have e := cocycle (one : H) one
  rw [one_mul, act_one] at e
  have e' : add (c (one : H)) zero = add (c (one : H)) (c (one : H)) := by rw [add_zero']; exact e
  exact (add_left_cancel' e').symm

theorem c_inv (h : H) : c (inv h) = neg (act (c h) (inv h)) := by
  apply neg_unique
  rw [← cocycle, mul_inv, c_one]

theorem inv_unique {a b : H} (h : mul a b = one) : b = inv a := by
  have e : mul (inv a) (mul a b) = b := by rw [← mul_assoc, inv_mul, one_mul]
  rw [h, mul_one] at e
  exact e.symm

theorem inv_mul_rev (a b : H) : inv (mul a b) = mul (inv b) (inv a) := by
  refine (inv_unique ?_).symm
  rw [mul_assoc, ← mul_assoc b, mul_inv, one_mul, mul_inv]

def gmul (g₁ g₂ : H × N) : H × N := (mul g₁.1 g₂.1, add (act g₁.2 g₂.1) g₂.2)

theorem gmul_assoc (a b d : H × N) : gmul (gmul a b) d = gmul a (gmul b d) := by
  simp only [gmul, mul_assoc, act_add, act_mul, add_assoc]

/-- repaired tilt (Prop. C3′): τ'(h) = (h, +h⁻¹ d_{A₀} h) -/
def tau' (h : H) : H × N := (h, c h)

theorem tau'_hom (h₁ h₂ : H) : tau' (mul h₁ h₂) = gmul (tau' h₁) (tau' h₂) := by
  simp only [tau', gmul, cocycle]

/-- π'(ε, ϖ) = ε (ϖ − c(ε)) ε⁻¹ = εϖε⁻¹ − (d_{A₀}ε)ε⁻¹  (Prop. C3′ item 4) -/
def pi' (g : H × N) : N := act (add g.2 (neg (c g.1))) (inv g.1)

theorem pi'_tau' (g : H × N) (h : H) : pi' (gmul g (tau' h)) = pi' g := by
  obtain ⟨e, w⟩ := g
  simp only [pi', gmul, tau', cocycle, inv_mul_rev, neg_add_distrib]
  have k : add (add (act w h) (c h)) (add (neg (act (c e) h)) (neg (c h))) = act (add w (neg (c e))) h := by
    rw [act_add, act_neg, show add (add (act w h) (c h)) (add (neg (act (c e) h)) (neg (c h)))
        = add (add (act w h) (neg (act (c e) h))) (add (c h) (neg (c h))) by ac_rfl, add_neg, add_zero']
  rw [k, act_mul, act_inv_cancel]

/-! ### Theorem G1 (algebraic core): unitary-gauge decomposition and factorisation -/

theorem inv_one' : inv (one : H) = one := (inv_unique (one_mul (one : H))).symm

/-- every g = (ε, ϖ) is the unitary-gauge point (1, π'g) moved by τ'(ε) -/
theorem decomp (g : H × N) : gmul (one, pi' g) (tau' g.1) = g := by
  obtain ⟨e, w⟩ := g
  simp only [gmul, tau', pi', one_mul, act_inv_cancel']
  rw [add_assoc, neg_add, add_zero']

/-- π' restricted to the unitary gauge ε = 1 is the identity of N (so π' is onto) -/
theorem pi'_unit (S : N) : pi' ((one : H), S) = S := by
  simp only [pi', c_one]
  rw [show neg (zero : N) = zero from (neg_unique (add_zero' zero)).symm, add_zero', inv_one', act_one]

/-- the fibres of π' are exactly the τ'(H)-orbits; the action is free (h is determined by ε₂ = ε₁h) -/
theorem fiber (g₁ g₂ : H × N) : pi' g₁ = pi' g₂ ↔ ∃ h : H, g₂ = gmul g₁ (tau' h) := by
  constructor
  · intro k
    refine ⟨mul (inv g₁.1) g₂.1, ?_⟩
    have e1 : g₁ = gmul (one, pi' g₁) (tau' g₁.1) := (decomp g₁).symm
    calc g₂ = gmul (one, pi' g₂) (tau' g₂.1) := (decomp g₂).symm
      _ = gmul (one, pi' g₁) (tau' (mul g₁.1 (mul (inv g₁.1) g₂.1))) := by
          rw [k, ← mul_assoc, mul_inv, one_mul]
      _ = gmul (gmul (one, pi' g₁) (tau' g₁.1)) (tau' (mul (inv g₁.1) g₂.1)) := by
          rw [tau'_hom, gmul_assoc]
      _ = gmul g₁ (tau' (mul (inv g₁.1) g₂.1)) := by rw [← e1]
  · rintro ⟨h, rfl⟩
    exact (pi'_tau' g₁ h).symm

theorem free (g : H × N) (h₁ h₂ : H) (k : gmul g (tau' h₁) = gmul g (tau' h₂)) : h₁ = h₂ := by
  have e := congrArg Prod.fst k
  simp only [gmul, tau'] at e
  have e2 := congrArg (mul (inv g.1)) e
  simp only [← mul_assoc, inv_mul, one_mul] at e2
  exact e2

/-- Theorem G1: any τ'(H)-invariant function I on 𝒢 (e.g. the action (9.4) under R-C1′, Prop. C4.4) satisfies
    I = J ∘ π' with J(S) := I(1, S), the same functional evaluated in the unitary gauge ε = 1. -/
theorem factor {X : Type} (I : H × N → X) (hI : ∀ (g : H × N) (h : H), I (gmul g (tau' h)) = I g)
    (g : H × N) : I g = I ((one : H), pi' g) := by
  calc I g = I (gmul (one, pi' g) (tau' g.1)) := by rw [decomp]
    _ = I (one, pi' g) := hI _ _

/-- and conversely every function of π' is τ'-invariant -/
theorem factor_converse {X : Type} (J : N → X) (g : H × N) (h : H) :
    J (pi' (gmul g (tau' h))) = J (pi' g) := by rw [pi'_tau']

end General
end GStueckelberg

#print axioms GStueckelberg.decomp
#print axioms GStueckelberg.pi'_unit
#print axioms GStueckelberg.fiber
#print axioms GStueckelberg.free
#print axioms GStueckelberg.factor
#print axioms GStueckelberg.factor_converse
