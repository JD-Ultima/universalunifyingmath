# Research Charter: Evidence-Led Mathematical Exploration

**Revision:** 2026-10-30. This charter updates the initial three-framework specification while retaining its evidence and correction requirements.

## Purpose

Investigate the explicit mathematical rules, computational mechanisms, and relationships that could explain physical regularities and yield testable predictions. The project's informal goal of finding the **“source code”** means making those rules precise enough to implement, analyze, and test. Their form and explanatory scope are questions to investigate.

The initial case studies are Chain Theory, Geometric Unity, and arithmetic and geometry associated with Terrence Howard. Additional approaches may be investigated when suitable sources and precise questions are available. No framework is assigned a preferred outcome before its statements are examined.

## Evidence standards

1. **Traceable sources.** Record the source, version, location, acquisition status, and a commit or hash where available. Distinguish primary material, transcripts, summaries, and commentary. Unavailable material does not supply premises.
2. **Explicit interpretation.** Identify the objects, operations, domains, units, conventions, and assumptions. Explore plausible alternative readings when the source is ambiguous. Attribute intentions only when the source explicitly supports them.
3. **Verifiable reasoning.** Support conclusions with a derivation, formal proof, exact counterexample, reproducible computation, or appropriately sourced measurement. State imported lemmas and precision limits.
4. **Scoped conclusions.** A result applies to the assumptions and family actually examined. A bounded search does not exclude constructions outside that family. A software check certifies the property it tests.
5. **Statistical accountability.** Numerical proximity is an observation. A probability or significance claim requires a defined random model, its justification, dependence and selection assumptions, and the comparison family. Local density and uniform-phase calculations must be identified as model-dependent estimates.
6. **Equal treatment of outcomes.** Record supporting results, contradictions, inconclusive checks, and failed approaches. Evaluate an assertion through its evidence rather than an author's reputation or membership in a particular school of thought.
7. **Separation of mathematics and observation.** Formal consistency, a derivation within a model, and empirical agreement are distinct results. A physical prediction needs an observable, an identification map, units, uncertainty, and a failure criterion.
8. **Visible corrections.** Correct active documents when an error is established and append the change to reports/CORRECTIONS.md. Preserve historical reports with their scope and any superseding notices.

Unsupported assertions are excluded from factual conclusions. Exploratory work is expressed as a question with the definitions or evidence needed to answer it.

## Investigation protocol

For each question:

1. Identify the source statement and assign a stable claim ID.
2. Define the mathematical or computational setting and record reconstruction choices.
3. State the question precisely, including hypotheses and scope.
4. Seek a proof, counterexample, computation, or measurement comparison.
5. Record the method, dependencies, precision, and reproducible command.
6. Classify the result as established, computed, conditional, contradicted under specified assumptions, or open.
7. Test relevant alternatives and state what would change the conclusion.
8. Update the claim register, result record, summary, and correction history when needed.

A research question can remain open. Missing definitions and unavailable sources are documented gaps, rather than reasons to infer either correctness or error.

## Exploratory observations

The patterns/ directory records observations and research questions. Each entry should identify its supporting calculation or source, its current limits, and a next test. An analogy, dimension match, or numerical coincidence supplies a question; an explanation requires a demonstrated relationship and, for physical claims, an empirical test.

Claims of novelty require a literature comparison. Separate implementations are useful checks, but internal AI-assisted reviews are not external peer review. Record which computations and formal proofs have actually been rerun in each validation pass.

## Reproducibility and reporting

- Use exact arithmetic where practical and disclose numerical tolerances elsewhere.
- Pin source snapshots and software versions; document platform requirements.
- Keep source statements, reconstruction choices, proofs, exploratory questions, and physical tests distinguishable.
- Connect significant results to claim IDs, evidence files, and result records.
- State whether a check uses committed logs or a newly executed build.
- Keep negative and inconclusive results available alongside supporting results.

The repository layout and selected verification commands are documented in README.md. Detailed commands and recorded outputs are maintained with each computation or proof.

## Criteria for a useful outcome

A useful investigation yields a reviewable reconstruction, a demonstrable result, a reproducible counterexample, a testable prediction, or a precise account of the evidence still required. Progress is measured by clearer questions and stronger evidence. Acceptance or rejection of a particular framework is an outcome of the investigation, never a prerequisite.
