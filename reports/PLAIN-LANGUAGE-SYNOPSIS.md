# What we built and tested

Imagine someone gives us two sets of building instructions and says they might help explain how the universe works. We are checking the instructions, one step at a time.

Our goal is to find clear rules that other people can check. We do not yet have a demonstrated “source code” for the universe.

## Our math workshop

We built a workshop with five kinds of tools:

- **A source shelf:** records of where each idea came from and which version we read.
- **Careful instructions:** definitions, assumptions and written mathematical arguments.
- **Math calculators:** Python programs that recount objects, calculate formulas, try alternatives and look for contradictions.
- **A proof checker:** Lean files that let a computer check certain mathematical statements, step by step.
- **An inspection notebook:** reviews, corrections and a list of questions we have not answered.

The programs and their recorded answers are saved so another person can repeat the work.

## Two sets of instructions

### Chain Theory

One construction makes 172 mathematical pieces and uses them to calculate a particular number. We counted the pieces again and got the same number.

Then we changed some of the rules to see which choices matter. The answer depends on choices that still need an independent reason. Getting the number again establishes that we can follow the recipe; showing why nature uses it would require more evidence.

We also checked the larger compendium. We asked whether its steps depend on unfinished steps, and whether its computer proofs prove the claims written beside them. Some proofs establish useful mathematical facts. Many check narrower statements, such as counting entries in a supplied list.

### Geometric Unity

This proposal uses complicated geometry. We wrote down explicit versions of its objects and equations, checked how the parts fit, and repaired a sign problem. We also found that the printed equation does not follow from the proposed action under the stated assumptions.

For specified versions of the model, we worked out which arrangements can satisfy the equations and which cannot. The newer work settles two particular branches on flat starting geometry and excludes several families of curved starting geometry. One result covers points with zero Ricci curvature but nonzero total curvature. Those are two different ways of measuring curvature.

Every such result has a named set of assumptions. Other fields and geometric cases remain open. The detailed map is in the [synthesis](SYNTHESIS.md#geometric-unity).

## What “the tests passed” means

The [verification record](MORNING-HANDOFF-2026-11-02.md#validation) shows which Chain Theory and Geometric Unity checks were run, their results, and their limits. The [README](../README.md#reproduce-selected-checks) explains how to run the current selected suite.

Some tests confirm an answer. Others confirm that a proposed rule fails, or that a known error is detected. Passing means the program produced the expected, checked result.

This suite covers selected work. The recorded Lean checks and manuscript builds have their own scope; the selected suite does not rebuild them. Internal reviews include AI-assisted work, and external peer review remains outstanding.

We have a better map of the ideas: parts that fit, parts that need repair, and parts we still need to explore. The next discoveries must earn their place through the same checks.
