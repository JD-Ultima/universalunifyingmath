# Research Charter: Formal Investigation of Three Unconventional Unification Frameworks

**Status:** Initial specification (v0.2, seed sources added)
**Subjects:** (1) Chain Theory (A. D. Wissner-Gross), (2) Geometric Unity (E. Weinstein), (3) Terrence Howard's "Grand Unified Theory" / Terryology

## 1. Mission

For each framework, determine which components can be stated as precise, internally consistent mathematics, and which of those yield testable or peer-reviewable results. The working question is:

> *What is the strongest well-defined mathematical structure this framework could be, and what does it actually imply?*

The investigation is hypothesis-neutral. Verdicts come from derivation and verification, not from authority, consensus, or reputation.

## 2. Operating principles

1. **Charitable reconstruction.** Construct the strongest precise interpretation(s) of ambiguous claims before evaluating them. Document every interpretive choice and test the alternatives.
2. **No inherited conclusions.** Critiques and endorsements are leads (objections to test, claims to check), never evidence.
3. **First-principles derivation.** Standard results are used as explicitly imported lemmas, with references.
4. **Pattern-seeking is welcome; pattern-claiming is gated.** Numerology and analogies are allowed as conjectures, and graduate to results only after formalization and verification.
5. **Falsifiability by construction.** Every claim is paired with a condition under which it would fail.
6. **Honest accounting.** Negative results, inconsistencies and dead ends are reported with the same rigor as positive ones.
7. **Separate the layers:**
   - (a) what the source asserts;
   - (b) our formal reconstruction;
   - (c) what is proven;
   - (d) what is conjectured;
   - (e) what is empirically testable.

## 3. Phase 1 — Research and source reconstruction

- **Source acquisition.** Primary sources with exact citations, versions, dates, URLs/DOIs and hashes. Informal sources are marked as such and quoted verbatim with locations. A source that cannot be located is flagged, not guessed.
- **Claim Register.** Fields: ID (`GU-…`, `CT-…`, `TH-…`), Statement, Type (Definition / Axiom / Theorem-claim / Physical prediction / Analogy / Numerological observation), Formalizable? (with reason), Source.
- **Structural mapping.** Mathematical objects relied on, candidate bridges to established mathematics, and Open Specification Gaps.
- **Pattern catalogue.** Conjectures awaiting formalization.

Seed sources:

- Geometric Unity: geometricunity.org (2013 Oxford lecture; April 1, 2021 manuscript draft).
- Chain Theory: alexwg.org talk; the α Letter with Lean 4 file and audit script; the compendium (release 0.1).
- Howard: to be located.

Chain Theory seed result:

- The closed form $(79+48\sqrt2+42\sqrt3+3\sqrt6)/31104$ has inverse $137.035\,999\,189\,566\ldots$.
- The author conjectures that it is $\alpha(0)$.
- Phase 0 confirmed the arithmetic only.
- Phase 2 must:
  1. re-derive the records and weights independently;
  2. audit each definitional clause;
  3. check the Lean and Python and build an independent implementation;
  4. assess the construction's freedom against the digits matched;
  5. evaluate the compendium's conditional and open results by the author's own grades.

## 4. Phase 2 — Mathematical formalization and verification

**Protocol.** For each claim:

1. Define the objects.
2. State the claim as a theorem, conjecture or definition.
3. Check that it is well-posed.
4. Attempt a proof or disproof.
5. Characterize the outcome: holds / fails / holds under conditions C / undecidable from the given axioms.
6. Explore variants.

**Verification standards:**

- Exact where possible. Numerical results state precision.
- Reproducible: code, pinned dependencies, run instructions.
- Significant results are cross-checked independently.
- Proof assistants are encouraged.
- Mathematical validity and physical validity are reported separately.

**Metric record per result:**

- formalization status;
- internal consistency;
- proof status;
- novelty;
- reproducibility;
- empirical testability;
- dependencies;
- confidence, with justification.

## 5. Approach

Rigorous but unconstrained. Look for the idea inside each claim, and cross-pollinate between frameworks. Keep a strict wall between exploration (labeled as speculative) and reporting (precise, verified, graded).

## 6–8. Repository, conventions, success criteria

The layout is given in `README.md`.

- Mathematics is written in LaTeX/Markdown.
- Code is Python (SymPy/NumPy/mpmath); formal proofs are in Lean 4.
- Every result links to its Claim ID, metric record and code.
- Corrections and retractions are logged, not erased (`reports/CORRECTIONS.md`).

Success means, for each framework:

1. a reviewable formal reconstruction, or a documented account of why none exists;
2. a graded ledger;
3. any new verifiable results, written to submission standard;
4. a clear separation between what the mathematics establishes and what remains speculation.
