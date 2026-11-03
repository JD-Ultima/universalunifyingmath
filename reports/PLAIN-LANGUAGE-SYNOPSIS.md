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

We also rebuilt an existing computer-checked proof of the mathematical read cap and connected it to the author's dimension-selection theorem. This checks those precise definitions; explaining how the proposed physical reads implement them still needs work.

We also checked the larger compendium. We asked whether its steps depend on unfinished steps, and whether its computer proofs prove the claims written beside them. Some proofs establish useful mathematical facts. Many check narrower statements, such as counting entries in a supplied list.

### Geometric Unity

This proposal uses complicated geometry. We wrote down explicit versions of its objects and equations, checked how the parts fit, and repaired a sign problem. We also found that the printed equation does not follow from the proposed action under the stated assumptions.

For specified versions of the model, we worked out which fields can satisfy the equations. We completely settled two field families across the positive Z1 models we defined. One allows only a fixed field over flat four-dimensional starting geometry, at one exact parameter relation. The other allows no solution. The larger observer geometry is curved in this example.

For a harder family, called Lambda5, we proved that rotating and resizing one fixed kind of shape cannot give a solution on an open region. A further proof rules out two specified constant-field endpoint cases: their equations force exactly that forbidden kind of shape. Another counting argument limits how many curvature directions can be zero or negative. More general changing fields remain open.

We showed that one known working field extends across a complete global model, including versions with a periodic four-dimensional starting space. Its geodesics—the geometry's straightest paths—have no artificial stopping point. This establishes mathematical examples within our reconstruction.

For two remaining endpoint cases, we rewrote the equations so that a field's value fixes how it must change. If a solution exists on a chosen geometry, its value at one point determines the rest. Finding values that fit the geometry everywhere remains unfinished.

Separate results exclude several families of curved starting geometry, including points with zero Ricci curvature but nonzero total curvature. Those are different ways of measuring curvature. Every result has a named set of assumptions; the [synthesis](SYNTHESIS.md#geometric-unity) gives the precise map.

## Existing work we can use

Other people have already checked parts of these ideas. We found published programs and computer-checked mathematics that overlap ours. We will verify those pieces and their assumptions, give credit, and use them where they fit. The [prior-work survey](RELATED-WORK.md) lists the useful starting points. None of the work we assessed proves either entire physical theory.

## What “the tests passed” means

The [verification guide](../computation/README.md#completed-validation) shows which Chain Theory and Geometric Unity checks were run, their results, and their limits. The [README](../README.md#reproduce-selected-checks) explains how to run the current selected suite.

Some tests confirm an answer. Others confirm that a proposed rule fails, or that a known error is detected. Passing means the program produced the expected, checked result.

This suite covers selected work. The recorded Lean checks and manuscript builds have their own scope; the selected suite does not rebuild them. Internal reviews include AI-assisted work, and external peer review remains outstanding.
