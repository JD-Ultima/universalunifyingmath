/-! # Convention check appended to the author's ChainAlpha.lean (sha256 522844e9…c331da)

The upstream file defines `actingCharge d1 d2 := d1 ≠ d2`. This encodes ONE convention: the
centre (x·1, y·1) of s(u(d1) ⊕ u(d2)), traceless on the direct sum (d1 x + d2 y = 0, primitive
solution x = d2, y = -d1), acting on the link as x + y (far end by ξ₂).
If the far end acts by conj ξ₂ (the far end transforms by Ū under a matching), it acts as x - y.
We show that under that convention the upstream uniqueness theorem's hypotheses are satisfied by a
second configuration, a shared qubit with two read (2,2) links and one unread (2,2) link.
Nothing here modifies an upstream definition. -/

namespace ConventionCheck
open ChainAlpha

/-- Value by which the primitive central element acts on a (d1,d2) link, upstream convention. -/
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

end ConventionCheck

#print axioms ConventionCheck.uniqueness_fails_under_conjugate_convention
#print axioms ConventionCheck.actA_ne_zero_iff
