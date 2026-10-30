/-! # Gauge-group selection (compendium 64.2, 65.1; audit note 03-gauge-group.md, G1-G3)

Standalone file. Lean 4.35.0, core library only (no Mathlib, no imports). Run with `run_gauge.sh`;
the output is `GaugeSelection.out`.

What is checked (section numbers below):

* §1 (G3, selection by read cap and acting centre). The six genres of the alphabet {2,3,4}; the readable
  types under the cap d1 d2 ≤ 6; the selection of (3,2) by "readable and acting centre"; each input
  dropped or varied in turn (cap, d1 ≠ d2, alphabet upper end, alphabet lower end, the conjugate
  convention of the alpha audit, the upstream "both ends ≤ 3" readability predicate); the exact set of
  caps for which (3,2) is the unique selection, over {2,3,4} and over the alphabet of all integers ≥ 2
  (hence over every alphabet of integers ≥ 2 that contains 2, 3 and 4); alphabets without 4 have
  other cap sets ({2,3}: every cap ≥ 6; {2,3,5}: caps 6-9).
* §2 (G2/G3, the genre-matrix formulas). dim f = d1² + d2² - 1, derived dim = dim f - 1,
  rank = d1 + d2 - 1, central scalar d2 - d1: the printed genre matrix; the "dim 12 and rank 4"
  uniqueness cell, and each of its two conditions alone, over the six genres and over all of ℕ.
* §3 (G1/G2, matrix level, any commutative ring). Matrices are functions `Fin n → Fin n → α` over a
  type with a core `Lean.Grind.CommRing` structure. Kernel of the face map (X, Y) ↦ X⊗1 + 1⊗Y; action
  of the central generator on the link; kernel of (A, B) ↦ A⊗B; triviality of that kernel on the
  determinant slice for (2,3); the rescaling (cA, c⁻¹B) that puts U(2) × U(3) onto the slice with the
  same image; the covering (A, B, z) ↦ (z³A, z⁻²B) maps into the slice and its kernel is {z : z⁶ = 1};
  for every genre, the kernel of the scalar slice map and of the covering in terms of z.
* §4 (finite exponent models of the roots of unity). In μ_M ≅ ℤ/M: the six kernel triples of the
  (2,3) covering and their generator (3,2,1) (the 65.1 wording); kernel orders lcm(d1, d2) for the six
  genres with the primitive covering, d1 d2 with the unreduced one; tensor-kernel orders |d1 - d2|.

What is NOT checked here (stated again in README.md): ℂ, unitarity, the Lie groups U(n), SU(n) and
S(U(2) × U(3)) as Lie groups, the isomorphism S(U(2) × U(3)) ≅ (SU(2) × SU(3) × U(1))/ℤ₆ as an
isomorphism of Lie groups, the mixed-product (homomorphism) property of ⊗, rank-nullity, the real
dimension of u(n), the derived algebra and the centre of the face algebra (Lie theory), and the fact
that μ_N ⊂ U(1) has N elements. The finite results of §1-§2 are about the selection formulas; that the
formulas are the Lie-algebra data rests on §3 plus those paper steps. -/

namespace GaugeSelection

/-! ## §1. Selection of the (3,2) genre (G3, second bullet) -/

/-- Genres over an alphabet: unordered pairs written (d1, d2) with d1 ≥ d2, listed by d2 then d1
(the order of `gauge_audit.py`). -/
def genres (A : List Nat) : List (Nat × Nat) :=
  A.flatMap (fun b => (A.filter (fun a => b ≤ a)).map (fun a => (a, b)))

/-- The admitted alphabet of 9.1: "the dimensions two, three and four". (Not 2.7 "The act alphabet",
an OPEN box about the two move letters attach and join.) -/
def alphabet : List Nat := [2, 3, 4]

/-- The read cap of 9.1 ("readable at the complete-read bound"), with the cap value as a parameter. -/
def readable (cap : Nat) (g : Nat × Nat) : Bool := decide (g.1 * g.2 ≤ cap)

/-- The acting centre of 9.2 in the convention of the source and of the upstream alpha Lean
(`ChainAlpha.actingCharge d1 d2 := d1 ≠ d2`): the central generator acts on the link as i(d2 - d1). -/
def acting (g : Nat × Nat) : Bool := g.1 != g.2

/-- The acting centre when the far end acts by the complex conjugate (alpha audit F-CT-3,
`ConventionCheck.actConj d1 d2 = d2 - (-d1)`): the generator acts as i(d2 + d1). Written in ℕ, where
d2 + d1 is the same number. -/
def actingConj (g : Nat × Nat) : Bool := g.2 + g.1 != 0

/-- The readability predicate of the upstream alpha Lean (`ChainAlpha.readable_iff`: both ends ≤ 3). -/
def readableEnds (g : Nat × Nat) : Bool := decide (g.1 ≤ 3) && decide (g.2 ≤ 3)

/-- "Readable and acting centre". -/
def select (cap : Nat) (A : List Nat) : List (Nat × Nat) :=
  (genres A).filter (fun g => readable cap g && acting g)

theorem six_genres : genres alphabet = [(2,2), (3,2), (4,2), (3,3), (4,3), (4,4)] := by decide

/-- The readable types under the cap 6. This is also the selection with the acting centre dropped. -/
theorem readable_types : (genres alphabet).filter (readable 6) = [(2,2), (3,2)] := by decide

/-- G3: readable and acting centre select (3,2) only. -/
theorem selection : select 6 alphabet = [(3,2)] := by decide

/-- Drop the read cap: the acting centre alone leaves three genres. -/
theorem drop_cap : (genres alphabet).filter acting = [(3,2), (4,2), (4,3)] := by decide

/-- The selection as a function of the cap value, for every cap c. -/
theorem cap_profile (c : Nat) : select c alphabet =
    (if 6 ≤ c then [(3,2)] else []) ++ (if 8 ≤ c then [(4,2)] else []) ++
      (if 12 ≤ c then [(4,3)] else []) := by
  simp only [select, genres, alphabet, readable, acting]
  by_cases h6 : 6 ≤ c <;> by_cases h8 : 8 ≤ c <;> by_cases h12 : 12 ≤ c <;>
    simp [List.filter, h6, h8, h12]

/-- The cap values for which (3,2) is the unique selection are exactly 6 and 7. -/
theorem unique_iff_cap (c : Nat) : select c alphabet = [(3,2)] ↔ 6 ≤ c ∧ c < 8 := by
  rw [cap_profile]
  by_cases h6 : 6 ≤ c <;> by_cases h8 : 8 ≤ c <;> by_cases h12 : 12 ≤ c <;>
    simp [h6, h8, h12] <;> omega

/-- The sensitivity values of 03-gauge-group.md G3. -/
theorem cap_4 : select 4 alphabet = [] := by decide
theorem cap_8 : select 8 alphabet = [(3,2), (4,2)] := by decide
theorem cap_9 : select 9 alphabet = [(3,2), (4,2)] := by decide
theorem cap_12 : select 12 alphabet = [(3,2), (4,2), (4,3)] := by decide

/-- Alphabet {2,...,8} at cap 6 (28 genres): the selection is unchanged. -/
theorem alphabet_2_to_8 : (genres (List.range' 2 7)).length = 28 ∧
    select 6 (List.range' 2 7) = [(3,2)] := by decide

/-- Core arithmetic of the selection, with no upper end on the alphabet: for d1 > d2 ≥ 2,
d1 d2 ≤ 6 iff (d1, d2) = (3, 2). -/
theorem cap6_acting_unbounded (d1 d2 : Nat) (h2 : 2 ≤ d2) (hlt : d2 < d1) :
    d1 * d2 ≤ 6 ↔ d1 = 3 ∧ d2 = 2 := by
  constructor
  · intro h
    have a : d1 * 2 ≤ d1 * d2 := Nat.mul_le_mul_left _ h2
    have b : 3 * d2 ≤ d1 * d2 := Nat.mul_le_mul_right _ (by omega)
    exact ⟨by omega, by omega⟩
  · rintro ⟨rfl, rfl⟩; decide

/-- The readable types with no upper end on the alphabet: (2,2) and (3,2). -/
theorem readable_unbounded (d1 d2 : Nat) (h2 : 2 ≤ d2) (hle : d2 ≤ d1) :
    d1 * d2 ≤ 6 ↔ (d1 = 2 ∧ d2 = 2) ∨ (d1 = 3 ∧ d2 = 2) := by
  constructor
  · intro h
    have a : d1 * 2 ≤ d1 * d2 := Nat.mul_le_mul_left _ h2
    have b : 2 * d2 ≤ d1 * d2 := Nat.mul_le_mul_right _ (by omega)
    have hd1 : d1 ≤ 3 := by omega
    have hd2 : d2 = 2 := by
      have : d2 * d2 ≤ d1 * d2 := Nat.mul_le_mul_right _ hle
      have : ¬ 3 ≤ d2 := fun h3 => by
        have : 3 * 3 ≤ d2 * d2 := Nat.mul_le_mul h3 h3
        omega
      omega
    subst hd2
    by_cases e : d1 = 2
    · exact Or.inl ⟨e, rfl⟩
    · exact Or.inr ⟨by omega, rfl⟩
  · rintro (⟨rfl, rfl⟩ | ⟨rfl, rfl⟩) <;> decide

/-- Any alphabet of integers ≥ 2 that contains 2 and 3: at cap 6 the selection is (3,2) and nothing
else (membership form; duplicates in the alphabet would only repeat the entry). -/
theorem selection_any_alphabet (A : List Nat) (hA : ∀ a ∈ A, 2 ≤ a) (h2 : 2 ∈ A) (h3 : 3 ∈ A)
    (g : Nat × Nat) : g ∈ select 6 A ↔ g = (3,2) := by
  obtain ⟨a, b⟩ := g
  simp only [select, genres, readable, acting, List.mem_filter, List.mem_flatMap, List.mem_map]
  constructor
  · rintro ⟨⟨b', hb', a', ⟨_, hle⟩, he⟩, hc⟩
    simp only [Prod.mk.injEq] at he
    obtain ⟨rfl, rfl⟩ := he
    have hb2 := hA _ hb'
    simp at hle hc
    have := (cap6_acting_unbounded a' b' hb2 (by omega)).mp hc.1
    simp [this]
  · intro he
    simp only [Prod.mk.injEq] at he
    obtain ⟨rfl, rfl⟩ := he
    exact ⟨⟨2, h2, 3, ⟨h3, by decide⟩, rfl⟩, by decide⟩

/-- Over the alphabet of all integers ≥ 2 (every pair d1 > d2 ≥ 2), (3,2) is the only readable acting
pair exactly for caps 6 and 7. Hence the same holds over every alphabet of integers ≥ 2 that contains
2, 3 and 4; it is not a statement about alphabets without 4 (`unique_iff_cap_23`, `unique_iff_cap_235`). -/
theorem unique_iff_cap_unbounded (c : Nat) :
    (∀ d1 d2, 2 ≤ d2 → d2 < d1 → (d1 * d2 ≤ c ↔ d1 = 3 ∧ d2 = 2)) ↔ 6 ≤ c ∧ c < 8 := by
  constructor
  · intro h
    have h32 := (h 3 2 (by decide) (by decide)).mpr ⟨rfl, rfl⟩
    have h42 := h 4 2 (by decide) (by decide)
    have : ¬ (4 * 2 ≤ c) := fun hc => by have := h42.mp hc; omega
    exact ⟨by omega, by omega⟩
  · rintro ⟨h6, h8⟩ d1 d2 h2 hlt
    constructor
    · intro h
      have a : d1 * 2 ≤ d1 * d2 := Nat.mul_le_mul_left _ h2
      have b : 3 * d2 ≤ d1 * d2 := Nat.mul_le_mul_right _ (by omega)
      have : ¬ (d1 ≥ 4) := fun h4 => by
        have : 4 * 2 ≤ d1 * d2 := Nat.mul_le_mul h4 h2
        omega
      exact ⟨by omega, by omega⟩
    · rintro ⟨rfl, rfl⟩; omega

/-- Alphabets of integers ≥ 2 without 4 have other cap sets. Over {2,3}: every cap ≥ 6. -/
theorem unique_iff_cap_23 (c : Nat) : select c [2, 3] = [(3,2)] ↔ 6 ≤ c := by
  simp only [select, genres, readable, acting]
  by_cases h6 : 6 ≤ c <;> simp [List.filter, h6]

/-- Over {2,3,5}: caps 6 to 9 ((5,2) enters at 10, (5,3) at 15). -/
theorem unique_iff_cap_235 (c : Nat) : select c [2, 3, 5] = [(3,2)] ↔ 6 ≤ c ∧ c < 10 := by
  simp only [select, genres, readable, acting]
  by_cases h6 : 6 ≤ c <;> by_cases h10 : 10 ≤ c <;> by_cases h15 : 15 ≤ c <;>
    simp [List.filter, h6, h10, h15] <;> omega

/-- Two instances (review R3-CT-10), and the interval 6 ≤ c < 8 failing at {2,3,5}. -/
theorem alphabet_23_cap100 : select 100 [2, 3] = [(3,2)] := by decide
theorem alphabet_235_cap9 : select 9 [2, 3, 5] = [(3,2)] := by decide
theorem alphabet_235_not_interval : ¬ (select 9 [2, 3, 5] = [(3,2)] ↔ 6 ≤ 9 ∧ 9 < 8) := by decide

/-- Drop the lower end of the alphabet (admit d = 1): at cap 6 four genres over {1,2,3,4}. -/
theorem alphabet_with_one : select 6 [1, 2, 3, 4] = [(2,1), (3,1), (4,1), (3,2)] := by decide

/-- ... and, with no upper end, (d, 1) for d = 2..6 besides (3,2). -/
theorem cap6_acting_with_one (d1 d2 : Nat) (h1 : 1 ≤ d2) (hlt : d2 < d1) :
    d1 * d2 ≤ 6 ↔ (d2 = 1 ∧ d1 ≤ 6) ∨ (d1 = 3 ∧ d2 = 2) := by
  by_cases e : d2 = 1
  · subst e; constructor
    · intro h; exact Or.inl ⟨rfl, by omega⟩
    · rintro (⟨_, h⟩ | ⟨_, h⟩) <;> omega
  · have h2 : 2 ≤ d2 := by omega
    rw [cap6_acting_unbounded d1 d2 h2 hlt]
    constructor
    · exact Or.inr
    · rintro (⟨h, _⟩ | h)
      · exact absurd h e
      · exact h

/-- Drop the cap with no upper end on the alphabet: every (d, 2), d ≥ 3, has an acting centre. -/
theorem drop_cap_unbounded (d : Nat) (h : 3 ≤ d) : acting (d, 2) = true := by
  simp only [acting, bne_iff_ne, ne_eq]; omega

/-- Conjugate convention: the acting criterion holds at every genre (d ≥ 1) ... -/
theorem actingConj_always (d1 d2 : Nat) (h : 1 ≤ d1) : actingConj (d1, d2) = true := by
  obtain ⟨k, rfl⟩ : ∃ k, d1 = k + 1 := ⟨d1 - 1, by omega⟩
  rfl

/-- ... so "readable and acting centre" selects what the cap alone selects: (2,2) and (3,2). -/
theorem conj_selection :
    (genres alphabet).filter (fun g => readable 6 g && actingConj g) = [(2,2), (3,2)] := by decide

/-- Upstream readability (both ends ≤ 3): three readable types; with the acting centre (3,2) is still
unique; under the conjugate convention all three survive (cf. `ReviewCheck`, review CT-4). -/
theorem upstream_readable : (genres alphabet).filter readableEnds = [(2,2), (3,2), (3,3)] := by decide
theorem upstream_selection :
    (genres alphabet).filter (fun g => readableEnds g && acting g) = [(3,2)] := by decide
theorem upstream_conj_selection :
    (genres alphabet).filter (fun g => readableEnds g && actingConj g) = [(2,2), (3,2), (3,3)] := by
  decide

/-! ## §2. Genre-matrix formulas and the "dim 12, rank 4" cell (G2; G3, first bullet) -/

/-- dim f(d1, d2) = d1² + d2² - 1 (the kernel of §3 `face_kernel` is one-dimensional). -/
def dimF (d1 d2 : Nat) : Nat := d1 * d1 + d2 * d2 - 1
/-- Derived algebra su(d1) ⊕ su(d2). -/
def derivedDim (d1 d2 : Nat) : Nat := (d1 * d1 - 1) + (d2 * d2 - 1)
/-- Rank: diagonal part, d1 + d2 directions with the same one-dimensional overlap. -/
def rankF (d1 d2 : Nat) : Nat := d1 + d2 - 1
/-- The central generator (i d2 1, -i d1 1) acts on the link as i(d2 - d1) (§3 `centre_acts`). -/
def centralScalar (d1 d2 : Nat) : Int := (d2 : Int) - (d1 : Int)

def sixGenres : List (Nat × Nat) := [(2,2), (3,2), (4,2), (3,3), (4,3), (4,4)]

theorem sixGenres_eq : sixGenres = genres alphabet := by decide

/-- The printed genre matrix (01-framework.tex:52-56), (dimension; central scalar):
(7;0), (12;1), (19;2), (17;0), (24;1), (31;0). -/
theorem genre_matrix :
    sixGenres.map (fun g => (dimF g.1 g.2, (centralScalar g.1 g.2).natAbs)) =
      [(7,0), (12,1), (19,2), (17,0), (24,1), (31,0)] := by decide

theorem derived_and_rank :
    sixGenres.map (fun g => (derivedDim g.1 g.2, rankF g.1 g.2)) =
      [(6,3), (11,4), (18,5), (16,5), (23,6), (30,7)] := by decide

/-- The printed scalar is |d2 - d1|, i.e. the normalisation (d2, -d1) of the generator. The primitive
integer solution of d1 x + d2 y = 0, (d2/g, -d1/g) with g = gcd, acts as (d2 - d1)/g: it differs at
(4,2). Only zero versus nonzero enters the selection. -/
theorem primitive_scalar :
    sixGenres.map (fun g => ((g.2 : Int) - g.1).natAbs / Nat.gcd g.1 g.2) = [0, 1, 1, 0, 1, 0] := by
  decide

/-- The acting criterion: nonzero central scalar iff d1 ≠ d2 (cf. upstream `acting_charge`). -/
theorem acting_iff (d1 d2 : Nat) : centralScalar d1 d2 ≠ 0 ↔ d1 ≠ d2 := by
  unfold centralScalar
  constructor
  · intro h e; subst e; omega
  · intro h e; omega

/-- dim f = dim su(d1) + dim su(d2) + 1 and rank = rank su(d1) + rank su(d2) + 1: the targets 12 and
4 are the dimension and rank of su(3) ⊕ su(2) ⊕ u(1). -/
theorem dimF_split (d1 d2 : Nat) (h1 : 1 ≤ d1) (h2 : 1 ≤ d2) :
    dimF d1 d2 = (d1 * d1 - 1) + (d2 * d2 - 1) + 1 ∧ derivedDim d1 d2 + 1 = dimF d1 d2 ∧
      rankF d1 d2 = (d1 - 1) + (d2 - 1) + 1 := by
  have a : 1 * 1 ≤ d1 * d1 := Nat.mul_le_mul h1 h1
  have b : 1 * 1 ≤ d2 * d2 := Nat.mul_le_mul h2 h2
  unfold dimF derivedDim rankF
  exact ⟨by omega, by omega, by omega⟩

theorem targets_are_SM : dimF 3 2 = 12 ∧ rankF 3 2 = 4 ∧ (3 * 3 - 1) + (2 * 2 - 1) + 1 = 12 ∧
    (3 - 1) + (2 - 1) + 1 = 4 := by decide

/-- G3 first bullet: dim 12 and rank 4 select (3,2) only. -/
theorem dim12_rank4 : sixGenres.filter (fun g => dimF g.1 g.2 == 12 && rankF g.1 g.2 == 4) = [(3,2)] := by
  decide
/-- Each condition alone already selects (3,2) among the six. -/
theorem dim12_alone : sixGenres.filter (fun g => dimF g.1 g.2 == 12) = [(3,2)] := by decide
theorem rank4_alone : sixGenres.filter (fun g => rankF g.1 g.2 == 4) = [(3,2)] := by decide

/-- dim 12 alone selects (3,2) among all pairs of natural numbers (13 = 3² + 2² only). -/
theorem dim12_unbounded (d1 d2 : Nat) (h : d2 ≤ d1) : dimF d1 d2 = 12 ↔ d1 = 3 ∧ d2 = 2 := by
  unfold dimF
  constructor
  · intro e
    have h13 : d1 * d1 + d2 * d2 = 13 := by omega
    have b1 : d1 < 4 := Nat.lt_of_not_le fun h4 => by
      have : 4 * 4 ≤ d1 * d1 := Nat.mul_le_mul h4 h4
      omega
    have b2 : d2 < 4 := by omega
    have key : ∀ a, a < 4 → ∀ b, b < 4 → b ≤ a → a * a + b * b = 13 → a = 3 ∧ b = 2 := by decide
    exact key d1 b1 d2 b2 h h13
  · rintro ⟨rfl, rfl⟩; decide

/-- rank 4 alone selects (3,2) among all pairs with d2 ≥ 2; with d2 = 1 admitted, (4,1) also. -/
theorem rank4_unbounded (d1 d2 : Nat) (h2 : 2 ≤ d2) (h : d2 ≤ d1) : rankF d1 d2 = 4 ↔ d1 = 3 ∧ d2 = 2 := by
  unfold rankF
  constructor
  · intro e; exact ⟨by omega, by omega⟩
  · rintro ⟨rfl, rfl⟩; decide
theorem rank4_with_one : rankF 4 1 = 4 := by decide

/-! ## §3. Matrix-level statements over an arbitrary commutative ring (G1, G2) -/

section Algebra
open Lean.Grind

variable {α : Type} [CommRing α]

abbrev Mat (α : Type) (n : Nat) := Fin n → Fin n → α
/-- Operators on the link C^m ⊗ C^n, indexed by pairs. -/
abbrev PMat (α : Type) (m n : Nat) := Fin m × Fin n → Fin m × Fin n → α

def scal (c : α) : Mat α n := fun i j => if i = j then c else 0
def scalP (c : α) : PMat α m n := fun x y => if x = y then c else 0
def smul (c : α) (A : Mat α n) : Mat α n := fun i j => c * A i j
/-- A ⊗ B. -/
def kron (A : Mat α m) (B : Mat α n) : PMat α m n := fun x y => A x.1 y.1 * B x.2 y.2
/-- The face map (X, Y) ↦ X ⊗ 1 + 1 ⊗ Y. -/
def face (X : Mat α m) (Y : Mat α n) : PMat α m n :=
  fun x y => (if x.2 = y.2 then X x.1 y.1 else 0) + (if x.1 = y.1 then Y x.2 y.2 else 0)
/-- Leibniz determinants in dimensions 2 and 3. -/
def det2 (A : Mat α 2) : α := A 0 0 * A 1 1 - A 0 1 * A 1 0
def det3 (A : Mat α 3) : α :=
  A 0 0 * (A 1 1 * A 2 2 - A 1 2 * A 2 1) - A 0 1 * (A 1 0 * A 2 2 - A 1 2 * A 2 0)
    + A 0 2 * (A 1 0 * A 2 1 - A 1 1 * A 2 0)

theorem eq_neg_of_add_eq_zero {a b : α} (h : a + b = 0) : a = -b := by
  have : a + b + -b = -b := by rw [h, AddCommMonoid.zero_add]
  rwa [Semiring.add_assoc, AddCommGroup.add_neg_cancel, Semiring.add_zero] at this

theorem mul_sub' (a b c : α) : a * (b - c) = a * b - a * c := by
  rw [Ring.sub_eq_add_neg, Ring.sub_eq_add_neg, Semiring.left_distrib, Ring.mul_neg]

theorem pow_mul' (z : α) (e d : Nat) : (z ^ e) ^ d = z ^ (e * d) := by
  induction d with
  | zero => rw [Semiring.pow_zero, Nat.mul_zero, Semiring.pow_zero]
  | succ d ih => rw [Semiring.pow_succ, ih, Nat.mul_succ, Semiring.pow_add]

theorem pow3 (c : α) : c ^ 3 = c * (c * c) := by
  rw [Semiring.pow_succ, Semiring.pow_two, CommSemiring.mul_comm]

/-- G2 / dimension: the kernel of the face map is the line {(c 1, -c 1)}, for every genre. With
rank-nullity and dim_ℝ u(d) = d² (not checked here) this gives dim f = d1² + d2² - 1; applied to
diagonal X, Y it gives rank d1 + d2 - 1. -/
theorem face_kernel (p0 : Fin m) (q0 : Fin n) (X : Mat α m) (Y : Mat α n) :
    face X Y = (fun _ _ => 0) ↔ ∃ c : α, X = scal c ∧ Y = scal (-c) := by
  constructor
  · intro h
    have e : ∀ x y, face X Y x y = 0 := fun x y => by rw [h]
    refine ⟨X p0 p0, ?_, ?_⟩
    · funext p p'
      by_cases hp : p = p'
      · subst hp
        have h1 := e (p, q0) (p, q0)
        have h2 := e (p0, q0) (p0, q0)
        simp only [face, ite_true] at h1 h2
        simp only [scal, ite_true]
        rw [eq_neg_of_add_eq_zero h1, eq_neg_of_add_eq_zero h2]
      · have h1 := e (p, q0) (p', q0)
        simp only [face, ite_true, hp, ite_false, Semiring.add_zero] at h1
        simp only [scal, hp, ite_false, h1]
    · funext q q'
      by_cases hq : q = q'
      · subst hq
        have h2 := e (p0, q) (p0, q)
        simp only [face, ite_true] at h2
        simp only [scal, ite_true]
        rw [Semiring.add_comm] at h2
        exact eq_neg_of_add_eq_zero h2
      · have h1 := e (p0, q) (p0, q')
        simp only [face, ite_true, hq, ite_false, AddCommMonoid.zero_add] at h1
        simp only [scal, hq, ite_false, h1]
  · rintro ⟨c, rfl, rfl⟩
    funext x y
    simp only [face, scal]
    by_cases h1 : x.1 = y.1 <;> by_cases h2 : x.2 = y.2 <;>
      simp only [h1, h2, ite_true, ite_false, Semiring.add_zero, AddCommGroup.add_neg_cancel]

/-- G2: the central element (a 1, -b 1) acts on the link as the scalar a - b. With a = i d2,
b = i d1 (the generator (i d2 1, -i d1 1)) this is i(d2 - d1). -/
theorem centre_acts (a b : α) :
    face (m := m) (n := n) (scal a) (scal (-b)) = scalP (a - b) := by
  funext x y
  simp only [face, scal, scalP]
  by_cases h1 : x.1 = y.1 <;> by_cases h2 : x.2 = y.2
  · have : x = y := Prod.ext h1 h2
    simp only [this, ite_true, Ring.sub_eq_add_neg]
  · have : x ≠ y := fun e => h2 (congrArg Prod.snd e)
    simp only [h1, h2, this, ite_true, ite_false, AddCommMonoid.zero_add]
  · have : x ≠ y := fun e => h1 (congrArg Prod.fst e)
    simp only [h1, h2, this, ite_true, ite_false, Semiring.add_zero]
  · have : x ≠ y := fun e => h1 (congrArg Prod.fst e)
    simp only [h1, h2, this, ite_false, Semiring.add_zero]

/-- The generator is traceless on C^{d1} ⊕ C^{d2}: d1 (t d2) + d2 (-(t d1)) = 0, written for
arbitrary ring elements a = d1, b = d2. -/
theorem centre_traceless (t a b : α) : a * (t * b) + b * (-(t * a)) = 0 := by
  rw [Ring.mul_neg, CommSemiring.mul_left_comm, CommSemiring.mul_comm a b,
    CommSemiring.mul_left_comm b, AddCommGroup.add_neg_cancel]

/-- G1: A ⊗ B = 1 iff A = λ1, B = μ1 with λμ = 1 (any commutative ring, any m, n ≥ 1). -/
theorem kron_kernel (p0 : Fin m) (q0 : Fin n) (A : Mat α m) (B : Mat α n) :
    kron A B = scalP 1 ↔ ∃ l u : α, l * u = 1 ∧ A = scal l ∧ B = scal u := by
  constructor
  · intro h
    have e : ∀ x y, kron A B x y = scalP 1 x y := fun x y => by rw [h]
    have hlu : A p0 p0 * B q0 q0 = 1 := by have := e (p0, q0) (p0, q0); simpa [kron, scalP] using this
    refine ⟨A p0 p0, B q0 q0, hlu, ?_, ?_⟩
    · funext p p'
      have h1 := e (p, q0) (p', q0)
      simp only [kron, scalP, Prod.mk.injEq, and_true] at h1
      have key : A p p' = (A p p' * B q0 q0) * A p0 p0 := by
        rw [Semiring.mul_assoc, CommSemiring.mul_comm (B q0 q0), hlu, Semiring.mul_one]
      rw [key, h1]
      by_cases hp : p = p'
      · simp [scal, hp, Semiring.one_mul]
      · simp [scal, hp, Semiring.zero_mul]
    · funext q q'
      have h1 := e (p0, q) (p0, q')
      simp only [kron, scalP, Prod.mk.injEq, true_and] at h1
      have key : B q q' = A p0 p0 * B q q' * B q0 q0 := by
        rw [Semiring.mul_assoc, CommSemiring.mul_comm (B q q'), ← Semiring.mul_assoc, hlu,
          Semiring.one_mul]
      rw [key, h1]
      by_cases hq : q = q'
      · simp [scal, hq, Semiring.one_mul]
      · simp [scal, hq, Semiring.zero_mul]
  · rintro ⟨l, u, hlu, rfl, rfl⟩
    funext x y
    simp only [kron, scal, scalP]
    by_cases h1 : x.1 = y.1 <;> by_cases h2 : x.2 = y.2
    · simp [Prod.ext h1 h2, hlu]
    · have : x ≠ y := fun e => h2 (congrArg Prod.snd e)
      simp [h1, h2, this, Semiring.mul_zero]
    · have : x ≠ y := fun e => h1 (congrArg Prod.fst e)
      simp [h1, h2, this, Semiring.zero_mul]
    · have : x ≠ y := fun e => h1 (congrArg Prod.fst e)
      simp [h1, h2, this, Semiring.zero_mul]

theorem det2_scal (l : α) : det2 (scal l) = l ^ 2 := by
  simp [det2, scal, Semiring.mul_zero, Ring.sub_eq_add_neg, AddCommGroup.neg_zero, Semiring.add_zero,
    Semiring.pow_two]

theorem det3_scal (l : α) : det3 (scal l) = l ^ 3 := by
  simp [det3, scal, Semiring.mul_zero, Semiring.zero_mul, Ring.sub_eq_add_neg, AddCommGroup.neg_zero,
    Semiring.add_zero, pow3]

theorem det2_smul (c : α) (A : Mat α 2) : det2 (smul c A) = c ^ 2 * det2 A := by
  simp only [det2, smul, mul_sub', Semiring.pow_two, Semiring.mul_assoc, CommSemiring.mul_left_comm,
    CommSemiring.mul_comm]

theorem det3_smul (c : α) (A : Mat α 3) : det3 (smul c A) = c ^ 3 * det3 A := by
  simp only [det3, smul, mul_sub', Semiring.left_distrib, pow3, Semiring.mul_assoc,
    CommSemiring.mul_left_comm, CommSemiring.mul_comm]

/-- Scalar core of G1: λμ = 1 and λ²μ³ = 1 force λ = μ = 1. -/
theorem slice_scalar_23 (l u : α) (h1 : l * u = 1) (h2 : l ^ 2 * u ^ 3 = 1) : l = 1 ∧ u = 1 := by
  have hu : u = 1 := by
    have : l ^ 2 * u ^ 3 = u := by
      rw [Semiring.pow_succ u 2, ← Semiring.mul_assoc, ← CommSemiring.mul_pow, h1, Semiring.one_pow,
        Semiring.one_mul]
    rw [← this, h2]
  refine ⟨?_, hu⟩
  rw [hu, Semiring.mul_one] at h1
  exact h1

/-- G1, injectivity: the kernel of S(U(2) × U(3)) → U(6), (A, B) ↦ A ⊗ B, is trivial (over any
commutative ring; unitarity is not used). -/
theorem tensor_kernel_trivial_23 (A : Mat α 2) (B : Mat α 3) (hdet : det2 A * det3 B = 1)
    (h : kron A B = scalP 1) : A = scal 1 ∧ B = scal 1 := by
  obtain ⟨l, u, hlu, rfl, rfl⟩ := (kron_kernel 0 0 A B).mp h
  rw [det2_scal, det3_scal] at hdet
  obtain ⟨rfl, rfl⟩ := slice_scalar_23 l u hlu hdet
  exact ⟨rfl, rfl⟩

/-- G1, image: (cA, c⁻¹B) with c = det A det B lies in the slice and has the image of (A, B). -/
theorem rescale_into_slice (A : Mat α 2) (B : Mat α 3) (c ci : α) (hc : c * ci = 1)
    (hdet : det2 A * det3 B = c) :
    det2 (smul c A) * det3 (smul ci B) = 1 ∧ kron (smul c A) (smul ci B) = kron A B := by
  constructor
  · rw [det2_smul, det3_smul]
    have : c ^ 2 * det2 A * (ci ^ 3 * det3 B) = (c * ci) ^ 2 * (ci * (det2 A * det3 B)) := by
      rw [CommSemiring.mul_pow, pow3, Semiring.pow_two, Semiring.pow_two]
      simp only [CommSemiring.mul_left_comm, CommSemiring.mul_comm]
    rw [this, hc, hdet, Semiring.one_pow, Semiring.one_mul, CommSemiring.mul_comm, hc]
  · funext x y
    simp only [kron, smul]
    have : c * A x.1 y.1 * (ci * B x.2 y.2) = (c * ci) * (A x.1 y.1 * B x.2 y.2) := by
      simp only [Semiring.mul_assoc, CommSemiring.mul_left_comm]
    rw [this, hc, Semiring.one_mul]

/-- The covering (A, B, z) ↦ (z³A, w²B), w = z⁻¹, maps SU(2) × SU(3) × U(1) into the slice. -/
theorem cover_into_slice (A : Mat α 2) (B : Mat α 3) (z w : α) (hzw : z * w = 1)
    (hA : det2 A = 1) (hB : det3 B = 1) : det2 (smul (z ^ 3) A) * det3 (smul (w ^ 2) B) = 1 := by
  rw [det2_smul, det3_smul, hA, hB, Semiring.mul_one, Semiring.mul_one, pow_mul', pow_mul',
    ← CommSemiring.mul_pow, hzw, Semiring.one_pow]

/-- If c A = 1 then A is the scalar matrix of its (0,0) entry, and c A₀₀ = 1. -/
theorem smul_eq_one (c : α) (A : Mat α (k + 1)) (h : smul c A = scal 1) :
    A = scal (A 0 0) ∧ c * A 0 0 = 1 := by
  have e : ∀ i j, c * A i j = scal 1 i j := fun i j => by
    have := congrFun (congrFun h i) j; simpa [smul] using this
  have h00 : c * A 0 0 = 1 := by have := e 0 0; simpa [scal] using this
  refine ⟨?_, h00⟩
  funext i j
  have key : A i j = (c * A i j) * A 0 0 := by
    rw [CommSemiring.mul_comm c, Semiring.mul_assoc, h00, Semiring.mul_one]
  rw [key, e i j]
  by_cases hij : i = j
  · simp [scal, hij, Semiring.one_mul]
  · simp [scal, hij, Semiring.zero_mul]

/-- G1, the Z6 kernel: (A, B, z) with det A = det B = 1 maps to (1, 1) iff z⁶ = 1 and
(A, B) = (z³ 1, z² 1). So the kernel is {z : z⁶ = 1}, carried bijectively by z. -/
theorem cover_kernel_23 (A : Mat α 2) (B : Mat α 3) (z w : α) (hzw : z * w = 1) :
    (det2 A = 1 ∧ det3 B = 1 ∧ smul (z ^ 3) A = scal 1 ∧ smul (w ^ 2) B = scal 1) ↔
      (z ^ 6 = 1 ∧ A = scal (z ^ 3) ∧ B = scal (z ^ 2)) := by
  have hzw2 : z ^ 2 * w ^ 2 = 1 := by rw [← CommSemiring.mul_pow, hzw, Semiring.one_pow]
  constructor
  · rintro ⟨hA, hB, h1, h2⟩
    obtain ⟨eA, ha⟩ := smul_eq_one _ A h1
    obtain ⟨eB, hb⟩ := smul_eq_one _ B h2
    -- B = scal b with w² b = 1, hence b = z²
    have hbz : B 0 0 = z ^ 2 := by
      have : B 0 0 = (z ^ 2 * w ^ 2) * B 0 0 := by rw [hzw2, Semiring.one_mul]
      rw [this, Semiring.mul_assoc, hb, Semiring.mul_one]
    have h6 : z ^ 6 = 1 := by
      rw [eB, hbz, det3_scal, pow_mul'] at hB; exact hB
    have haz : A 0 0 = z ^ 3 := by
      have : A 0 0 = A 0 0 * z ^ 6 := by rw [h6, Semiring.mul_one]
      rw [this, show (6 : Nat) = 3 + 3 from rfl, Semiring.pow_add, ← Semiring.mul_assoc,
        CommSemiring.mul_comm (A 0 0), ha, Semiring.one_mul]
    exact ⟨h6, by rw [eA, haz], by rw [eB, hbz]⟩
  · rintro ⟨h6, rfl, rfl⟩
    refine ⟨?_, ?_, ?_, ?_⟩
    · rw [det2_scal, pow_mul']; exact h6
    · rw [det3_scal, pow_mul']; exact h6
    · funext i j
      simp only [smul, scal]
      by_cases hij : i = j
      · simp only [hij, ite_true]
        rw [← Semiring.pow_add]; exact h6
      · simp only [hij, ite_false, Semiring.mul_zero]
    · funext i j
      simp only [smul, scal]
      by_cases hij : i = j
      · simp only [hij, ite_true]
        rw [CommSemiring.mul_comm]; exact hzw2
      · simp only [hij, ite_false, Semiring.mul_zero]

/-- G2, every genre (scalar level): on the line λμ = 1 the slice condition λ^{d1} μ^{d2} = 1 is
λ^{d1 - d2} = 1 (d1 ≥ d2). So the kernel of S(U(d1) × U(d2)) → U(d1 d2) is μ_{d1-d2}, and the whole
circle when d1 = d2. (Uses `kron_kernel` and det(λ1_d) = λ^d, checked above for d = 2, 3 only.) -/
theorem tensor_kernel_scalar (d1 d2 : Nat) (h : d2 ≤ d1) (l u : α) (hlu : l * u = 1) :
    l ^ d1 * u ^ d2 = 1 ↔ l ^ (d1 - d2) = 1 := by
  have e : l ^ d1 * u ^ d2 = l ^ (d1 - d2) := by
    rw [show d1 = (d1 - d2) + d2 from (Nat.sub_add_cancel h).symm, Nat.add_sub_cancel,
      Semiring.pow_add, Semiring.mul_assoc, ← CommSemiring.mul_pow, hlu, Semiring.one_pow,
      Semiring.mul_one]
  rw [e]

/-- G2, every genre (scalar level). For the covering (A, B, z) ↦ (z^{e2} A, z^{-e1} B) with
e2 d1 = e1 d2 = L, a kernel element has A = a 1, B = b 1 with z^{e2} a = 1 and b = z^{e1}; then
A ∈ SU(d1), B ∈ SU(d2) (a ∈ μ_{d1}, b ∈ μ_{d2}) holds iff z^L = 1. With e1 = d1/g, e2 = d2/g
(g = gcd, the primitive covering) L = lcm(d1, d2); with e1 = d1, e2 = d2, L = d1 d2. -/
theorem cover_kernel_scalar (d1 d2 e1 e2 L : Nat) (h1 : e2 * d1 = L) (h2 : e1 * d2 = L) (z a b : α) :
    (a ^ d1 = 1 ∧ b ^ d2 = 1 ∧ z ^ e2 * a = 1 ∧ b = z ^ e1) ↔
      (z ^ L = 1 ∧ z ^ e2 * a = 1 ∧ b = z ^ e1) := by
  constructor
  · rintro ⟨_, hb, ha, rfl⟩
    refine ⟨?_, ha, rfl⟩
    rw [← h2, ← pow_mul']; exact hb
  · rintro ⟨hz, ha, rfl⟩
    refine ⟨?_, ?_, ha, rfl⟩
    · have : (z ^ e2 * a) ^ d1 = 1 := by rw [ha, Semiring.one_pow]
      rw [CommSemiring.mul_pow, pow_mul', h1, hz, Semiring.one_mul] at this
      exact this
    · rw [pow_mul', h2]; exact hz

/-- The kernel element is determined by z (a is the inverse of z^{e2}). -/
theorem cover_kernel_unique (z a a' : α) (e : Nat) (ha : z ^ e * a = 1) (ha' : z ^ e * a' = 1) :
    a = a' := by
  have : a = (z ^ e * a') * a := by rw [ha', Semiring.one_mul]
  rw [this, CommSemiring.mul_comm (z ^ e), Semiring.mul_assoc, ha, Semiring.mul_one]

end Algebra

/-- The core class is satisfiable by a concrete ring (sanity: the hypotheses are not vacuous). -/
example : Lean.Grind.CommRing Int := inferInstance
theorem kron_kernel_int_instance :
    kron (scal (-1 : Int) : Mat Int 2) (scal (-1 : Int) : Mat Int 3) = scalP 1 := by
  funext x y
  simp only [kron, scal, scalP]
  by_cases h1 : x.1 = y.1 <;> by_cases h2 : x.2 = y.2
  · simp [Prod.ext h1 h2]
  · have : x ≠ y := fun e => h2 (congrArg Prod.snd e)
    simp [h1, h2, this]
  · have : x ≠ y := fun e => h1 (congrArg Prod.fst e)
    simp [h1, h2, this]
  · have : x ≠ y := fun e => h1 (congrArg Prod.fst e)
    simp [h1, h2, this]
/-- (-1, -1) is in the kernel of ⊗ but not in the slice: det = (-1)² (-1)³ = -1. -/
theorem minus_one_not_in_slice : det2 (scal (-1 : Int)) * det3 (scal (-1 : Int)) = -1 := by
  decide

/-- The 65.1 slice representative (3, 3, -2, -2, -2) is the centre generator (d2 1_{d1}, -d1 1_{d2})
at (d1, d2) = (2, 3) (factor i dropped), and it is traceless. -/
def centreDiag (d1 d2 : Nat) : List Int := List.replicate d1 (d2 : Int) ++ List.replicate d2 (-(d1 : Int))
theorem slice_representative : centreDiag 2 3 = [3, 3, -2, -2, -2] ∧ (centreDiag 2 3).sum = 0 := by
  decide
theorem centreDiag_traceless_six : sixGenres.map (fun g => (centreDiag g.1 g.2).sum) = [0, 0, 0, 0, 0, 0] := by
  decide

/-! ## §4. Finite exponent models of the roots of unity -/

/-- Exponent model of the kernel of the covering inside μ_M ≅ ℤ/M (ω = exp(2πi/M)):
triples (α, β, γ) with a = ω^α ∈ μ_{d1}, b = ω^β ∈ μ_{d2}, z = ω^γ, z^{e2} a = 1, b = z^{e1}. -/
def kerCover (M d1 d2 e1 e2 : Nat) : List (Nat × Nat × Nat) :=
  (List.range M).flatMap fun a => (List.range M).flatMap fun b => ((List.range M).filter fun c =>
    d1 * a % M == 0 && d2 * b % M == 0 && (e2 * c + a) % M == 0 && b == e1 * c % M).map
      fun c => (a, b, c)

/-- 65.1: "the kernel of the covering is re-met from the three register centres at the six triples
... with the generator (3,2,1)". In μ_6 the kernel of (A, B, z) ↦ (z³A, z⁻²B) is k·(3, 2, 1), k = 0..5. -/
theorem six_triples :
    kerCover 6 2 3 2 3 = [(0,0,0), (0,2,4), (0,4,2), (3,0,3), (3,2,1), (3,4,5)] ∧
    ((kerCover 6 2 3 2 3).all fun t => t == (3 * t.2.2 % 6, 2 * t.2.2 % 6, t.2.2)) = true ∧
    ((List.range 6).all fun k => (kerCover 6 2 3 2 3).contains (3 * k % 6, 2 * k % 6, k)) = true := by
  decide +kernel

/-- In the larger ambient μ_12 the same kernel still has six elements, all with z ∈ μ_6
(even exponent), generated by (6, 4, 2) in units of 2π/12, i.e. (3, 2, 1) in units of 2π/6. -/
theorem six_triples_in_mu12 : (kerCover 12 2 3 2 3).length = 6 ∧
    ((kerCover 12 2 3 2 3).all fun t =>
      t.2.2 % 2 == 0 && t == (6 * (t.2.2 / 2) % 12, 4 * (t.2.2 / 2) % 12, t.2.2)) = true := by
  decide +kernel

/-- G2: kernel orders of the primitive coverings (e1, e2) = (d1/g, d2/g) for the six genres, inside
μ_12 (every lcm divides 12): lcm(d1, d2) = 2, 6, 4, 3, 12, 4. -/
theorem cover_kernel_orders :
    sixGenres.map (fun g => (kerCover 12 g.1 g.2 (g.1 / Nat.gcd g.1 g.2) (g.2 / Nat.gcd g.1 g.2)).length)
      = sixGenres.map (fun g => Nat.lcm g.1 g.2) ∧
    sixGenres.map (fun g => Nat.lcm g.1 g.2) = [2, 6, 4, 3, 12, 4] := by decide +kernel

/-- The lcm needs the primitive normalisation: at (4,2), inside μ_24, the covering with
(e1, e2) = (2, 1) has kernel order 4 = lcm, the unreduced (4, 2) has 8 = d1 d2. -/
theorem cover_kernel_normalisation :
    (kerCover 24 4 2 2 1).length = 4 ∧ (kerCover 24 4 2 4 2).length = 8 := by decide +kernel

/-- Exponent model of the kernel of S(U(d1) × U(d2)) → U(d1 d2) on scalars inside μ_M:
pairs (α, β) with λμ = 1 and λ^{d1} μ^{d2} = 1. -/
def kerTensor (M d1 d2 : Nat) : List (Nat × Nat) :=
  (List.range M).flatMap fun a => ((List.range M).filter fun b =>
    (a + b) % M == 0 && (d1 * a + d2 * b) % M == 0).map fun b => (a, b)

/-- G2: inside μ_12 the tensor kernel has |d1 - d2| elements for d1 ≠ d2 (1, 2, 1) and all 12 for
d1 = d2 (the whole circle in U(1)). -/
theorem tensor_kernel_orders :
    sixGenres.map (fun g => (kerTensor 12 g.1 g.2).length) = [12, 1, 2, 12, 1, 12] := by decide +kernel

end GaugeSelection

/-! ## Readable tables (not certificates; the theorems above are) -/
#eval (List.range 17).map (fun c => (c, GaugeSelection.select c GaugeSelection.alphabet))
#eval GaugeSelection.kerCover 6 2 3 2 3

/-! ## Axiom audit -/
#print axioms GaugeSelection.six_genres
#print axioms GaugeSelection.readable_types
#print axioms GaugeSelection.selection
#print axioms GaugeSelection.drop_cap
#print axioms GaugeSelection.cap_profile
#print axioms GaugeSelection.unique_iff_cap
#print axioms GaugeSelection.cap_4
#print axioms GaugeSelection.cap_8
#print axioms GaugeSelection.cap_9
#print axioms GaugeSelection.cap_12
#print axioms GaugeSelection.alphabet_2_to_8
#print axioms GaugeSelection.cap6_acting_unbounded
#print axioms GaugeSelection.readable_unbounded
#print axioms GaugeSelection.selection_any_alphabet
#print axioms GaugeSelection.unique_iff_cap_unbounded
#print axioms GaugeSelection.unique_iff_cap_23
#print axioms GaugeSelection.unique_iff_cap_235
#print axioms GaugeSelection.alphabet_23_cap100
#print axioms GaugeSelection.alphabet_235_cap9
#print axioms GaugeSelection.alphabet_235_not_interval
#print axioms GaugeSelection.alphabet_with_one
#print axioms GaugeSelection.cap6_acting_with_one
#print axioms GaugeSelection.drop_cap_unbounded
#print axioms GaugeSelection.actingConj_always
#print axioms GaugeSelection.conj_selection
#print axioms GaugeSelection.upstream_readable
#print axioms GaugeSelection.upstream_selection
#print axioms GaugeSelection.upstream_conj_selection
#print axioms GaugeSelection.sixGenres_eq
#print axioms GaugeSelection.genre_matrix
#print axioms GaugeSelection.derived_and_rank
#print axioms GaugeSelection.primitive_scalar
#print axioms GaugeSelection.acting_iff
#print axioms GaugeSelection.dimF_split
#print axioms GaugeSelection.targets_are_SM
#print axioms GaugeSelection.dim12_rank4
#print axioms GaugeSelection.dim12_alone
#print axioms GaugeSelection.rank4_alone
#print axioms GaugeSelection.dim12_unbounded
#print axioms GaugeSelection.rank4_unbounded
#print axioms GaugeSelection.rank4_with_one
#print axioms GaugeSelection.eq_neg_of_add_eq_zero
#print axioms GaugeSelection.mul_sub'
#print axioms GaugeSelection.pow_mul'
#print axioms GaugeSelection.pow3
#print axioms GaugeSelection.face_kernel
#print axioms GaugeSelection.centre_acts
#print axioms GaugeSelection.centre_traceless
#print axioms GaugeSelection.kron_kernel
#print axioms GaugeSelection.det2_scal
#print axioms GaugeSelection.det3_scal
#print axioms GaugeSelection.det2_smul
#print axioms GaugeSelection.det3_smul
#print axioms GaugeSelection.slice_scalar_23
#print axioms GaugeSelection.tensor_kernel_trivial_23
#print axioms GaugeSelection.rescale_into_slice
#print axioms GaugeSelection.cover_into_slice
#print axioms GaugeSelection.smul_eq_one
#print axioms GaugeSelection.cover_kernel_23
#print axioms GaugeSelection.tensor_kernel_scalar
#print axioms GaugeSelection.cover_kernel_scalar
#print axioms GaugeSelection.cover_kernel_unique
#print axioms GaugeSelection.kron_kernel_int_instance
#print axioms GaugeSelection.minus_one_not_in_slice
#print axioms GaugeSelection.slice_representative
#print axioms GaugeSelection.centreDiag_traceless_six
#print axioms GaugeSelection.six_triples
#print axioms GaugeSelection.six_triples_in_mu12
#print axioms GaugeSelection.cover_kernel_orders
#print axioms GaugeSelection.cover_kernel_normalisation
#print axioms GaugeSelection.tensor_kernel_orders
