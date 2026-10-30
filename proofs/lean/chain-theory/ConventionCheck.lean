/-! # Convention check appended to the author's ChainAlpha.lean (sha256 522844e9…c331da)

The upstream file defines `actingCharge d1 d2 := d1 ≠ d2`. This encodes ONE convention: the
centre (x·1, y·1) of s(u(d1) ⊕ u(d2)), traceless on the direct sum (d1 x + d2 y = 0, solution
x = d2, y = -d1, primitive when gcd(d1,d2) = 1), acting on the link as x + y (far end by ξ₂).
Only whether that value vanishes is used, which does not depend on the normalisation.
If the far end acts by conj ξ₂ (the far end transforms by Ū under a matching), it acts as x - y.
We show that under that convention the upstream uniqueness theorem's hypotheses are satisfied by a
second configuration, a shared qubit with two read (2,2) links and one unread (2,2) link.
Nothing here modifies an upstream definition. -/

namespace ConventionCheck
open ChainAlpha

/-- Value by which the central element (d2·1, -d1·1) acts on a (d1,d2) link, upstream convention
(primitive when gcd(d1,d2) = 1; only its vanishing is used). -/
def actA (d1 d2 : Nat) : Int := (d2 : Int) + (-(d1 : Int))
/-- Same element, far end acting by the complex conjugate. -/
def actConj (d1 d2 : Nat) : Int := (d2 : Int) - (-(d1 : Int))

theorem actA_ne_zero_iff (d1 d2 : Nat) : actA d1 d2 ≠ 0 ↔ actingCharge d1 d2 := by
  unfold actA actingCharge; omega

theorem actConj_ne_zero (d1 d2 : Nat) (h : 2 ≤ d1) : actConj d1 d2 ≠ 0 := by
  unfold actConj; omega

def actingChargeConj (d1 d2 : Nat) : Prop := actConj d1 d2 ≠ 0

/-- The shared-qubit configuration 2R 2R 2U. -/
def qubitConfig : Config := ⟨2, [⟨2, true⟩, ⟨2, true⟩, ⟨2, false⟩]⟩

theorem qubitConfig_ne_object : qubitConfig ≠ objectConfig := by decide

theorem qubitConfig_twoPoint : qubitConfig.twoPoint := by
  unfold Config.twoPoint qubitConfig; decide

theorem qubitConfig_endRule : qubitConfig.endRule := by
  unfold Config.endRule qubitConfig
  intro a ha _
  simp at ha
  rcases ha with rfl | rfl | rfl <;> decide +kernel

theorem qubitConfig_readable : ∀ a ∈ qubitConfig.arms, 2 ≤ a.far ∧ readable qubitConfig.center a.far := by
  intro a ha
  simp [qubitConfig] at ha
  rcases ha with rfl | rfl | rfl <;>
    exact ⟨by decide, by unfold readable readers spansEnd; decide⟩

theorem qubitConfig_actingConj : ∀ a ∈ qubitConfig.arms, actingChargeConj qubitConfig.center a.far := by
  intro a _
  exact actConj_ne_zero _ _ (by decide)

/-- Under the conjugate convention, the hypotheses of `configuration_unique` (with `actingCharge`
replaced by `actingChargeConj`) hold for a configuration other than the object, so the
uniqueness statement does not survive the change of convention. -/
theorem uniqueness_fails_under_conjugate_convention :
    ∃ c : Config, c ≠ objectConfig ∧ 2 ≤ c.center ∧ c.twoPoint ∧ c.endRule ∧
      (∀ a ∈ c.arms, 2 ≤ a.far ∧ readable c.center a.far) ∧
      (∀ a ∈ c.arms, actingChargeConj c.center a.far) :=
  ⟨qubitConfig, qubitConfig_ne_object, by decide, qubitConfig_twoPoint, qubitConfig_endRule,
   qubitConfig_readable, qubitConfig_actingConj⟩

/-- And the upstream convention excludes it, as the paper states. -/
theorem upstream_excludes_qubitConfig : ¬ (∀ a ∈ qubitConfig.arms, actingCharge qubitConfig.center a.far) := by
  intro h
  have := h ⟨2, true⟩ (by simp [qubitConfig])
  exact this rfl

/-! ## Proposition 2 (audit F-CT-2): the weight law is a class function

On the object, the squared weight `kOf` that the upstream `readout` uses equals a function of
three data only: whether the record is a commutator record (own sign -1), its parity under the
single slot exchange `exchAssign`, and whether its two letters sit at one end. -/

def classK (r : Rec) : Nat :=
  if r.sign == -1 then 6
  else if parity r.form == 1 then (if sameEnd r.form then 3 else 2) else 1

theorem weight_is_class_function :
    (recs 3 2).all (fun r => kOf (wX 3 2) objectWorld objectTwoPoint r == classK r) = true := by
  decide +kernel

end ConventionCheck

#print axioms ConventionCheck.uniqueness_fails_under_conjugate_convention
#print axioms ConventionCheck.actA_ne_zero_iff
#print axioms ConventionCheck.weight_is_class_function

/-! ## Referee addition (review 2026-10-30, CT-4): the upstream `readable` predicate is d ≤ 3, not the
read cap d1 d2 ≤ 6, so under the conjugate convention a third configuration (3R 3R 3U at a shared qutrit)
also satisfies the hypotheses of `configuration_unique`. The #evals confirm the readouts of the (2,2) and (3,3) worlds. -/
namespace ReviewCheck
open ChainAlpha ConventionCheck

/-- Shared qutrit, two read (3,3) links and one unread (3,3) link. -/
def cfg33 : Config := ⟨3, [⟨3, true⟩, ⟨3, true⟩, ⟨3, false⟩]⟩

theorem cfg33_ne_object : cfg33 ≠ objectConfig := by decide
theorem cfg33_twoPoint : cfg33.twoPoint := by unfold Config.twoPoint cfg33; decide
theorem cfg33_endRule : cfg33.endRule := by
  unfold Config.endRule cfg33
  intro a ha _
  simp at ha
  rcases ha with rfl | rfl | rfl <;> decide +kernel
theorem cfg33_readable : ∀ a ∈ cfg33.arms, 2 ≤ a.far ∧ readable cfg33.center a.far := by
  intro a ha
  simp [cfg33] at ha
  rcases ha with rfl | rfl | rfl <;>
    exact ⟨by decide, by unfold readable readers spansEnd; decide⟩
theorem cfg33_actingConj : ∀ a ∈ cfg33.arms, actingChargeConj cfg33.center a.far := by
  intro a _
  exact actConj_ne_zero _ _ (by decide)

/-- A third configuration satisfying all hypotheses of `configuration_unique` under A*. -/
theorem third_config_under_conjugate_convention :
    ∃ c : Config, c ≠ objectConfig ∧ c ≠ qubitConfig ∧ 2 ≤ c.center ∧ c.twoPoint ∧ c.endRule ∧
      (∀ a ∈ c.arms, 2 ≤ a.far ∧ readable c.center a.far) ∧
      (∀ a ∈ c.arms, actingChargeConj c.center a.far) :=
  ⟨cfg33, cfg33_ne_object, by decide, by decide, cfg33_twoPoint, cfg33_endRule,
   cfg33_readable, cfg33_actingConj⟩

def qubitWorld : World where
  qudits := 4
  dim := fun _ => 2
  links := [(0, 1), (0, 2), (0, 3)]
def qubitTP : TwoPoint := { center := 0, read1 := (0, 1), read2 := (0, 2), state := (0, 3) }
def quditWorld : World where
  qudits := 4
  dim := fun _ => 3
  links := [(0, 1), (0, 2), (0, 3)]

#eval readout qubitWorld qubitTP
#eval (recs 2 2).length
#eval readout quditWorld qubitTP
end ReviewCheck
#print axioms ReviewCheck.third_config_under_conjugate_convention
