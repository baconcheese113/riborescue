# ADR-0021 — The shuffle control is a permutation test, not a permutation

**Status:** accepted · **Date:** 2026-07-25 · **Deciders:** Joseph

## Context

**This record exists because a control failed, and it says so.** ADR-0020's shuffle control refused
the kinetics claim on GSE144140. Examining why exposed a defect in how the control was specified —
one that is a defect whichever way the result had fallen, and which would have gone unnoticed had
the control passed. The order matters and is stated plainly rather than smoothed over.

ADR-0020 fixes the control as a single permutation, with a bootstrap interval taken over the ten
rounds of that one shuffled mapping. Those ten training folds overlap heavily by construction — a
90% training fraction, ten times — so the interval describes how *stably* one arrangement helps. It
does not describe how often an arbitrary arrangement helps at all, which is the question a shuffle
control is asking.

The consequence is that a single draw decides the verdict. On this data, seed 721 gave:

| shuffle | interval | verdict under ADR-0020 |
|---|---|---|
| global | +0.00072 [+0.00029, +0.00113] | fails — excludes zero |
| within-gene | +0.00045 [−0.00017, +0.00103] | clears |
| context-matched | +0.00168 [−0.00054, +0.00358] | clears |

A forty-draw diagnostic put the null's spread at roughly 0.0006 with individual draws reaching
+0.0027 — wider than the largest observed gain. On that same seed the global shuffle made SRI look
like a leak and SJ6986 look clean, from the arrangement alone. Four codon tables built from numbers
that were never a measurement gave −0.00038, **+0.00131**, −0.00003 and −0.00008: one in four landed
near a real drug's real gain. The control is a lottery, and its ticket number is `SEED`.

A second defect sits beside it. Six drugs are predeclared and each was tested against its own null,
with no correction for testing six. They are not six independent experiments — one library, the same
genes, the same contexts — so a correction that assumed independence would be wrong in the other
direction.

## Decision

**The control becomes an ordinary permutation test, and this record does not disturb ADR-0020.**

**ADR-0020's outcome stands exactly as recorded.** Its rule was frozen before the numbers existed,
it was applied as written, and it refused the claim. That verdict is the project's primary,
preregistered result for Wave 2 and is not revisited, revised or reinterpreted here. What follows is
a **secondary, method-corrected analysis**, labelled as one wherever it appears. It is not a
retroactive preregistered pass, and no outcome of it converts ADR-0020's refusal into an
acceptance.

**Nine hundred and ninety-nine independent permutations per shuffle family.** All three families —
global, within-gene, context-matched — because deciding from the global shuffle alone would repeat
the mistake of resting a conclusion on one arrangement of one procedure. The empirical p-value is
`(1 + #{null ≥ observed}) / (n + 1)`, which includes the observation in its own null and therefore
cannot report zero; at 999 draws the smallest attainable value is 0.001, and that resolution is
reported beside every p-value so a reader can see what the test was capable of.

**The permutation is synchronised across drugs.** One shuffled mapping is applied to every drug
within a permutation, so the null carries the dependence between drugs that the data has. Six
measurements of one library are one family.

**The statistic is the maximum over drugs.** Each drug's observed gain is compared against the
distribution of the largest gain any drug reached under the same shuffled mapping — single-step
Westfall-Young. Chosen because it matches the worst-case logic ADR-0020's controls already use, and
fixed here rather than selected once the p-values existed. Where a per-drug procedure is also
wanted, Benjamini-Hochberg is applied across **all six** predeclared drugs and never across a subset
chosen by looking at the gains.

**Folds are identical between the observed and every permuted fit.** Generated once from the fixed
seed under the grouped-by-gene split, which is where ADR-0020 places the claim. A null built on
freshly drawn splits would carry the split variation as well as the shuffle's.

**The baseline is fitted once per drug and round.** It reads no kinetic column, so its held-out
score is identical under every permutation. This is what makes a thousand permutations a few hours
rather than a few days, and it also removes a source of drift: the observed and permuted gains are
differences against the same numbers.

**Four quantities are reported, and they are not interchangeable.**

- absolute ΔR², which is what was measured;
- ΔR² as a share of the reliability headroom that drug's baseline leaves under its own ceiling
  (ADR-0020), which is what makes six drugs comparable;
- the familywise empirical p-value and its resolution;
- the ten-round bootstrap, named **split stability** — it says how consistent a gain is across
  overlapping folds, and it is never described as independent scientific uncertainty.

**The vocabulary is restrained and fixed.** Where the corrected analysis does not clear, the finding
is "**no multiplicity-controlled evidence** that codon occupancy carries information beyond local
sequence context, at this permutation resolution". It is not "kinetics carries no information",
which is a claim about the world that a failure to reject does not license. Nor is the published
factorisation described as **adequate**: adequacy is an equivalence claim and would need a margin
declared in advance, which this project has not declared and does not assert. Where a single drug
clears, it is reported as a **small therapy-specific secondary signal**, never as a general
transferable kinetics result.

**A failing control is a scientific outcome, not broken software.** The control suite records what
each control returned and asserts the conjunction — `Criteria.supported` is the only thing entitled
to describe kinetics as carrying transferable information. An unsupported claim is a valid green
state. What the suite refuses is the claim being made while a criterion does not hold, and the
recorded outcome is pinned so that a change making it supported turns the suite red until somebody
decides that deliberately.

## Consequences

Wave 2 now has two results and they are labelled differently wherever they are written down: a
preregistered procedural refusal under ADR-0020, and a method-corrected secondary analysis under
this record. A reader who meets only one of them should still be able to tell which.

The correction generalises. Any future shuffle control in this project is a permutation test with a
declared number of draws, a synchronised permutation across whatever family it is correcting over,
and a reported resolution. The single-permutation form is not used again.

Nine hundred and ninety-nine permutations across three families and six drugs is roughly thirty
core-hours. It is sharded by permutation index, and it is the only heavy job that runs at a time.

The defect was invisible while the control passed and became visible only when it failed. That is
worth stating: a control which decides on one random draw will agree with a correct control most of
the time, and the times it does not are exactly the times anyone is paying attention.
