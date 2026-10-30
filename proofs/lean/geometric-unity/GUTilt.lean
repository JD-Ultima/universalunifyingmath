/-
  GUTilt.lean — abstract algebraic core of Props. C1, C2, C3, C3′
  (formalization/geometric-unity/C-gauge-shiab.md; scope and limits: E-lean-tilt.md).
  Lean 4 core only (no Mathlib, no imports).  Check:  lean proofs/lean/geometric-unity/GUTilt.lean

  Model.  H a group; N an abelian group; `act x h` a right action of H on N by additive automorphisms
  (draft: ϖ ↦ h⁻¹ ϖ h on Ω¹(Y, ad)); `c : H → N` a 1-cocycle, c (h₁h₂) = act (c h₁) h₂ + c h₂
  (draft: c(h) = h⁻¹ d_{A₀} h, lemma L-C1).  A connection A₀ + a is represented by its offset a ∈ N, so
  the gauge action A ↦ h⁻¹Ah + h⁻¹dh becomes  gauge a h = act a h + c h.
  Draft formulas (round-1 reading R-C1):
    (5.9)  (ε₁,ϖ₁)(ε₂,ϖ₂) = (ε₁ε₂, ε₂⁻¹ϖ₁ε₂ + ϖ₂)          gmul
    (5.10) A·(ε,ϖ) = A·ε + ϖ                                  actR
    (5.11) (ε,ϖ)·A = (A + ϖ)·ε⁻¹                              actL
    (6.4)  τ(h) = (h, −h⁻¹d_{A₀}h)        repair: τ'(h) = (h, +h⁻¹d_{A₀}h)
    (6.11) π(ε,ϖ) = εϖε⁻¹ + (d_{A₀}ε)ε⁻¹  (6.15) g·γ = ε(ϖ+γ)ε⁻¹ + (d_{A₀}ε)ε⁻¹
    (6.18) μ(ε,ϖ) = (A₀ + ϖ, A₀·ε)        (7.3)  T(ε,ϖ) = ϖ − ε⁻¹d_{A₀}ε
-/
namespace GUTilt

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

/-! ### C1. The inhomogeneous gauge group 𝒢 = H ⋉ N and its actions -/

def gmul (g₁ g₂ : H × N) : H × N := (mul g₁.1 g₂.1, add (act g₁.2 g₂.1) g₂.2)
def gone : H × N := (one, zero)
def ginv (g : H × N) : H × N := (inv g.1, neg (act g.2 (inv g.1)))

theorem gmul_assoc (a b d : H × N) : gmul (gmul a b) d = gmul a (gmul b d) := by
  simp only [gmul, mul_assoc, act_add, act_mul, add_assoc]

theorem gmul_gone (g : H × N) : gmul g gone = g := by
  simp only [gmul, gone, mul_one, act_one, add_zero']

theorem gone_gmul (g : H × N) : gmul gone g = g := by
  simp only [gmul, gone, one_mul, act_zero, zero_add]

theorem gmul_ginv (g : H × N) : gmul g (ginv g) = gone := by
  simp only [gmul, ginv, gone, mul_inv, add_neg]

theorem ginv_gmul (g : H × N) : gmul (ginv g) g = gone := by
  simp only [gmul, ginv, gone, inv_mul, act_neg, act_inv_cancel', neg_add]

/-- gauge action on connections A₀ + a:  (A₀ + a)·h = A₀ + h⁻¹ a h + h⁻¹ d_{A₀} h -/
def gauge (a : N) (h : H) : N := add (act a h) (c h)
def actR (a : N) (g : H × N) : N := add (gauge a g.1) g.2            -- (5.10)
def actL (g : H × N) (a : N) : N := gauge (add a g.2) (inv g.1)      -- (5.11)

theorem gauge_mul (a : N) (h₁ h₂ : H) : gauge (gauge a h₁) h₂ = gauge a (mul h₁ h₂) := by
  simp only [gauge, act_add, act_mul, cocycle, add_assoc]

theorem actR_mul (a : N) (g₁ g₂ : H × N) : actR (actR a g₁) g₂ = actR a (gmul g₁ g₂) := by
  simp only [actR, gauge, gmul, act_add, act_mul, cocycle]
  ac_rfl

theorem actL_mul (g₁ g₂ : H × N) (a : N) : actL (gmul g₁ g₂) a = actL g₁ (actL g₂ a) := by
  simp only [actL, gauge, gmul, inv_mul_rev, act_mul, act_add, cocycle, act_inv_cancel]
  ac_rfl

/-! ### C2. The printed tilt τ (6.4): homomorphism, right stabiliser of A₀, Lemma 6.3 -/

def tau (h : H) : H × N := (h, neg (c h))
def tau' (h : H) : H × N := (h, c h)

theorem tau_hom (h₁ h₂ : H) : tau (mul h₁ h₂) = gmul (tau h₁) (tau h₂) := by
  simp only [tau, gmul, cocycle, neg_add_distrib, act_neg]

theorem stab_right (g : H × N) : actR (zero : N) g = zero ↔ g = tau g.1 := by
  obtain ⟨e, w⟩ := g
  simp only [actR, gauge, act_zero, zero_add, tau]
  constructor
  · intro h; rw [neg_unique h]
  · intro h; cases h; exact add_neg _

def pi (g : H × N) : N := act (add g.2 (c g.1)) (inv g.1)                 -- (6.11)
def actN (g : H × N) (γ : N) : N := act (add (add g.2 γ) (c g.1)) (inv g.1)   -- (6.15)

theorem pi_tau (g : H × N) (h : H) : pi (gmul g (tau h)) = pi g := by
  obtain ⟨e, w⟩ := g
  simp only [pi, gmul, tau, cocycle, inv_mul_rev]
  have k : add (add (act w h) (neg (c h))) (add (act (c e) h) (c h)) = act (add w (c e)) h := by
    rw [act_add, show add (add (act w h) (neg (c h))) (add (act (c e) h) (c h))
        = add (add (act w h) (act (c e) h)) (add (neg (c h)) (c h)) by ac_rfl, neg_add, add_zero']
  rw [k, act_mul, act_inv_cancel]

theorem pi_mul (g₁ g₂ : H × N) : pi (gmul g₁ g₂) = actN g₁ (pi g₂) := by
  simp only [pi, actN, gmul, cocycle, inv_mul_rev, act_mul, act_add, act_inv_cancel]
  ac_rfl

/-! ### C3. Lemmas 6.5 and 7.2 for the printed τ: exact residual −2c(h) -/

def mu (g : H × N) : N × N := (g.2, c g.1)                                 -- (6.18), offsets from A₀
def gauge2 (m : N × N) (h : H) : N × N := (gauge m.1 h, gauge m.2 h)
def T (g : H × N) : N := add g.2 (neg (c g.1))                             -- (7.3)

theorem mu_tau_residual (g : H × N) (h : H) :
    mu (gmul g (tau h)) = (add (gauge (mu g).1 h) (neg (add (c h) (c h))), gauge (mu g).2 h) := by
  obtain ⟨e, w⟩ := g
  simp only [mu, gmul, tau, gauge, cocycle, neg_add_distrib]
  rw [show add (add (act w h) (c h)) (add (neg (c h)) (neg (c h)))
      = add (act w h) (add (add (c h) (neg (c h))) (neg (c h))) by ac_rfl, add_neg, zero_add]

theorem T_tau_residual (g : H × N) (h : H) :
    T (gmul g (tau h)) = add (act (T g) h) (neg (add (c h) (c h))) := by
  obtain ⟨e, w⟩ := g
  simp only [T, gmul, tau, cocycle, neg_add_distrib, act_add, act_neg]
  ac_rfl

/-- Printed Lemma 7.2 holds for a given h iff 2 c(h) = 0. -/
theorem T_tau_iff (g : H × N) (h : H) :
    T (gmul g (tau h)) = act (T g) h ↔ add (c h) (c h) = (zero : N) := by
  rw [T_tau_residual]
  constructor
  · intro k
    have k' : add (act (T g) h) (neg (add (c h) (c h))) = add (act (T g) h) zero := by
      rw [add_zero']; exact k
    have z := add_left_cancel' k'
    rw [← neg_neg' (add (c h) (c h)), z]
    exact (neg_unique (add_zero' zero)).symm
  · intro k
    rw [k, show neg (zero : N) = zero from (neg_unique (add_zero' zero)).symm, add_zero']

/-! ### C3′. The repair τ'(h) = (h, +c(h)) -/

theorem tau'_hom (h₁ h₂ : H) : tau' (mul h₁ h₂) = gmul (tau' h₁) (tau' h₂) := by
  simp only [tau', gmul, cocycle]

theorem stab_left (g : H × N) : actL g (zero : N) = zero ↔ g = tau' g.1 := by
  obtain ⟨e, w⟩ := g
  simp only [actL, gauge, zero_add, c_inv, tau']
  rw [← act_neg, ← act_add]
  constructor
  · intro k
    have k2 := congrArg (fun x => act x e) k
    simp only [act_inv_cancel', act_zero] at k2
    rw [add_comm] at k2
    rw [neg_unique k2, neg_neg']
  · intro k
    cases k
    rw [add_neg, act_zero]

theorem B_equivariant (e h : H) : gauge (c e : N) h = c (mul e h) := by   -- (A₀·ε)·h = A₀·(εh), any tilt
  simp only [gauge, cocycle]

theorem mu_tau' (g : H × N) (h : H) : mu (gmul g (tau' h)) = gauge2 (mu g) h := by
  simp only [mu, gmul, tau', gauge2, gauge, cocycle]

theorem T_tau' (g : H × N) (h : H) : T (gmul g (tau' h)) = act (T g) h := by
  obtain ⟨e, w⟩ := g
  simp only [T, gmul, tau', cocycle, neg_add_distrib, act_add, act_neg]
  rw [show add (add (act w h) (c h)) (add (neg (act (c e) h)) (neg (c h)))
      = add (add (act w h) (neg (act (c e) h))) (add (c h) (neg (c h))) by ac_rfl, add_neg, add_zero']

def pi' (g : H × N) : N := act (add g.2 (neg (c g.1))) (inv g.1)

theorem pi'_tau' (g : H × N) (h : H) : pi' (gmul g (tau' h)) = pi' g := by
  obtain ⟨e, w⟩ := g
  simp only [pi', gmul, tau', cocycle, inv_mul_rev, neg_add_distrib]
  have k : add (add (act w h) (c h)) (add (neg (act (c e) h)) (neg (c h))) = act (add w (neg (c e))) h := by
    rw [act_add, act_neg, show add (add (act w h) (c h)) (add (neg (act (c e) h)) (neg (c h)))
        = add (add (act w h) (neg (act (c e) h))) (add (c h) (neg (c h))) by ac_rfl, add_neg, add_zero']
  rw [k, act_mul, act_inv_cancel]

/-- (ε, ϖ) ↦ (ε, −ϖ) is an automorphism of 𝒢 exchanging τ and τ' (P9, P10). -/
def phi (g : H × N) : H × N := (g.1, neg g.2)

theorem phi_hom (a b : H × N) : phi (gmul a b) = gmul (phi a) (phi b) := by
  simp only [phi, gmul, neg_add_distrib, act_neg]

theorem phi_tau (h : H) : phi (tau h) = (tau' h : H × N) := by
  simp only [phi, tau, tau', neg_neg']

theorem phi_invol (g : H × N) : phi (phi g) = g := by
  simp only [phi, neg_neg']

end General

/-! ### Counter-model: the printed Lemmas 6.5, 7.2 and the τ'-invariance of π do not follow.
    H = (ℤ, +) acting trivially on N = (ℤ, +), c = id (a cocycle: c(h₁+h₂) = c h₁ + c h₂). -/

instance : Grp Int where
  mul := (· + ·)
  one := 0
  inv := (- ·)
  mul_assoc := by intros; omega
  one_mul := by intros; omega
  mul_one := by intros; omega
  inv_mul := by intros; omega
  mul_inv := by intros; omega

instance : AbG Int where
  add := (· + ·)
  zero := 0
  neg := (- ·)
  add_assoc := by intros; omega
  add_comm := by intros; omega
  zero_add := by intros; omega
  neg_add := by intros; omega

instance : Cocycle Int Int where
  act := fun x _ => x
  c := fun h => h
  act_add := by intros; rfl
  act_mul := by intros; rfl
  act_one := by intros; rfl
  cocycle := by intros; rfl

theorem lemma72_as_printed_false :
    ¬ (∀ (g : Int × Int) (h : Int), T (gmul g (tau h)) = act (T g) h) := by
  intro k
  have := k (0, 0) 1
  revert this
  decide

theorem lemma65_as_printed_false :
    ¬ (∀ (g : Int × Int) (h : Int), mu (gmul g (tau h)) = gauge2 (mu g) h) := by
  intro k
  have := k (0, 0) 1
  revert this
  decide

theorem tau'_not_right_stabiliser : ¬ (∀ h : Int, actR (zero : Int) (tau' h) = zero) := by
  intro k
  have := k 1
  revert this
  decide

theorem pi_not_tau'_invariant : ¬ (∀ (g : Int × Int) (h : Int), pi (gmul g (tau' h)) = pi g) := by
  intro k
  have := k (0, 0) 1
  revert this
  decide

end GUTilt

-- Axiom audit (no `sorry`; only core axioms may appear):
#print axioms GUTilt.gmul_assoc
#print axioms GUTilt.actL_mul
#print axioms GUTilt.stab_right
#print axioms GUTilt.pi_mul
#print axioms GUTilt.T_tau_residual
#print axioms GUTilt.T_tau_iff
#print axioms GUTilt.stab_left
#print axioms GUTilt.T_tau'
#print axioms GUTilt.pi'_tau'
#print axioms GUTilt.phi_hom
#print axioms GUTilt.lemma72_as_printed_false
#print axioms GUTilt.pi_not_tau'_invariant
