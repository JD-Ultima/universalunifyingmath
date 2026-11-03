# Research charter

**Scope: Chain Theory and Geometric Unity. Updated November 3, 2026.**

Investigate explicit mathematical rules, computational mechanisms and relationships that could explain physical regularities. The informal goal of finding an underlying “source code” means making those rules precise enough to implement and test. Its existence, form and explanatory reach remain questions.

Neither framework has a preferred outcome. Assess supporting results, contradictions and inconclusive checks by the same standards.

## Philosophy of discovery

Combine mathematical abstraction, physical intuition and empirical testing. Use imagination to generate possibilities and disciplined reasoning to establish what follows.

- **Seek hidden structure.** Begin with the working belief that difficult problems may contain useful symmetries, invariants or simple underlying relations. Test that belief against alternatives; elegance alone supplies no premise or evidence.
- **Let intuition open the question.** Use pictures, analogies and thought experiments before demanding a finished formalism. Turn promising insights into explicit hypotheses that can be proved, tested or rejected.
- **Let nature judge physical models.** Reliable empirical conflict requires revising or rejecting the model, its assumptions or its physical identification. A valid conditional mathematical theorem retains its stated conclusion; whether nature satisfies its hypotheses is a separate question.

## Evidence standards

1. **Sources.** Record the version, location and acquisition status, with a commit or hash where available. Distinguish primary material from transcripts, summaries and commentary. Unavailable material supplies no premises.
2. **Definitions.** State objects, operations, domains, units, conventions, reality conditions and other assumptions. Examine plausible alternative readings when sources are ambiguous. Attribute intentions only when the source supports them.
3. **Reasoning.** Supply a derivation, formal proof, exact counterexample, reproducible computation or sourced measurement comparison. Identify imported lemmas, software versions and precision limits. Prefer exact arithmetic where practical; disclose numerical tolerances elsewhere.
4. **Scope.** A theorem applies to its hypotheses, a finite search to its declared family, and a computer check to its implemented statement. State whether a result is proved, computed, conditional, contradicted under specified assumptions, or open.
5. **Physical comparisons.** Separate mathematical consistency from empirical agreement. Specify the observable, identification map, units and uncertainty. A significance claim additionally requires a justified probability model, selection rule and account of dependence; numerical proximity alone is insufficient.
6. **Review and corrections.** Internal AI-assisted reviews are distinct from external peer review. Novelty claims require a literature comparison. Correct active documents when an error is established and append the change to the [correction log](reports/CORRECTIONS.md), preserving relevant evidence and historical scope.

## Investigation workflow

1. **Map the structure.** State the question, core variables and relationships. Distinguish definitions, established constraints and provisional assumptions. Look for invariants, conservation laws and symmetries. Check primary sources and reusable results; the [prior-work survey](reports/RELATED-WORK.md) records existing candidates.
2. **Run thought experiments.** Examine limiting cases, extreme scales and simplified worlds. Use visual models and analogies to generate competing explanations, then identify what would distinguish them.
3. **Build a useful approximation.** Use order-of-magnitude estimates, heuristics or bounded prototypes when they advance the question. State their domain, error estimates and unresolved uncertainty. A provisional construction can guide the next test before its full mechanism or proof is known.
4. **Test and revise.** Seek proofs and exact counterexamples for mathematical claims, and discriminating observations or experiments for physical claims. State what would count against a hypothesis. When reliable evidence conflicts with it, revise the assumptions or reject the proposal; preserve the result and its scope.

Move between these phases as needed. If one approach stalls, change the representation, scale or method. Keep process proportionate to the question, and choose elegance for its explanatory or practical value. First-principles reasoning, verified prior work and useful approximations all have a place; none overrides the evidence standards.

Assign stable claim IDs and record sources, assumptions, methods, dependencies, code, outputs and limitations. Distinguish recorded build logs from newly executed builds. Update the relevant claim, gap and result records, then the [current synthesis](reports/SYNTHESIS.md).

Keep useful negative and inconclusive results. Missing definitions are documented gaps rather than grounds for assuming either correctness or error.

## Exploratory observations

The [pattern catalogue](patterns/) records supported observations and research questions with their limits and next tests. An analogy, matching dimension or numerical coincidence motivates a question; a shared mechanism requires an explicit relationship and, for a physical claim, empirical evidence.

Progress means clearer questions and stronger evidence. An investigation may end with a scoped result, a reproducible failure or a precise unresolved gap.
