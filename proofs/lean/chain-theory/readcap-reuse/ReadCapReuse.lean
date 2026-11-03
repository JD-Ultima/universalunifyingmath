import MathLab.ChainTheory.Woronowicz
import ChainAlpha

/-!
Adapters between the pinned operator-theoretic read cap and the author's
unchanged abstract binder. These reuse existing proofs; they do not formalize
single-copy certainty, physical distillation, the weight law or the photon claim.
See verification.json for exact source, dependency and compiler revisions.
-/

namespace ReadCapReuse

theorem author_read_cap : ChainAlpha.ReadCap ChoiDark.Readable :=
  ChainWoronowicz.readcap

theorem operator_iff_readableByReads (d1 d2 : Nat) (h1 : 2 ≤ d1) (h2 : 2 ≤ d2) :
    ChoiDark.Readable d1 d2 ↔ ChainAlpha.readableByReads d1 d2 :=
  (ChainWoronowicz.readcap d1 d2 h1 h2).trans
    (ChainAlpha.readers_reproduce_the_cap d1 d2 h1 h2).symm

theorem selection (d1 d2 : Nat) (h1 : 2 ≤ d1) (h2 : 2 ≤ d2)
    (hr : ChoiDark.Readable d1 d2) (hq : ChainAlpha.actingCharge d1 d2) :
    (d1 = 2 ∧ d2 = 3) ∨ (d1 = 3 ∧ d2 = 2) :=
  ChainAlpha.selection_under_read_cap ChoiDark.Readable
    author_read_cap d1 d2 h1 h2 hr hq

theorem weaker_readable_three_three : ChainAlpha.readable 3 3 :=
  (ChainAlpha.readable_iff 3 3).2 ⟨by decide, by decide⟩

theorem not_operator_readable_three_three : ¬ ChoiDark.Readable 3 3 := by
  intro hr
  have h := (ChainWoronowicz.readcap 3 3 (by decide) (by decide)).1 hr
  exact (by decide : ¬ ((3 : Nat) * 3 ≤ 6)) h

end ReadCapReuse

#print ChainReadCap.ReadCap
#print ChainAlpha.ReadCap
#check ChainWoronowicz.woronowicz
#check ChainWoronowicz.readable_two_three
#check ChainWoronowicz.readcap
#check ReadCapReuse.author_read_cap
#check ReadCapReuse.operator_iff_readableByReads
#check ReadCapReuse.selection
#check ReadCapReuse.weaker_readable_three_three
#check ReadCapReuse.not_operator_readable_three_three

#print axioms ChainWoronowicz.sep_of_extreme
#print axioms ChainWoronowicz.woronowicz
#print axioms ChainWoronowicz.readable_two_three
#print axioms ChainWoronowicz.readcap
#print axioms ReadCapReuse.author_read_cap
#print axioms ReadCapReuse.operator_iff_readableByReads
#print axioms ReadCapReuse.selection
#print axioms ReadCapReuse.weaker_readable_three_three
#print axioms ReadCapReuse.not_operator_readable_three_three
