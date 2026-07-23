# ADR-0011 — Footprint length selection, and parity with the authors of GSE144140

**Status:** accepted · **Date:** 2026-07-23 · **Deciders:** Joseph

## Context

The calibration keeps footprints of 26 to 34 nucleotides. That range was never derived from anything
measured; it was a plausible default carried from the first library and never revisited.

GSE144140 makes it matter. Its footprints are shorter than the libraries the range was set on —
surveyed medians of 27 nt for the three DMSO libraries, 26 for the three G418 libraries, and 25 to 26
for the three SRI-37240 libraries. A floor at 26 discards roughly 40 to 45% of a control library and
about half of each treated one. That loss is correlated with treatment, in the same direction for
both compounds, inside the exact comparison the assay makes.

The authors of that dataset published their analysis code (`github.com/jrw24/SRI37240`, Python 2.7)
and deposited their count tables in the series. Reading it settles several questions that would
otherwise have been answered by invention.

**Their library handling matches ours.** All nine libraries use linker
`NNNNNNCACTCGGGCACCAAGGAC` with four nucleotides trimmed from the 5' end. Our samplesheet already
declares exactly this. (`trim5Plist` in their settings selects the linker; it is not a trim length.)

**Their estimator is ours.** The Ribosome ReadThrough Score is ribosome density between the normal
termination codon and the first in-frame termination codon, divided by coding-sequence density, per
transcript. That is `downstream_occupancy` as ADR-0010 defines it, arrived at independently.

**Four parameters differ, and each is a number, not a judgement.**

| | Theirs | Ours |
|---|---|---|
| Footprint lengths | 28–35 | 26–34 |
| Coding denominator | excludes first 18 nt and last 15 nt | the whole coding sequence |
| Downstream window opens | 6 nt past the stop | 1 nt past the stop |
| Normalisation | reads per million before ratios | pooled counts, ratio after |

They also define an empty-A-site population at 21 to 24 nucleotides and an empty-E-site population at
18 to 19, separately from the 28 to 35 full-length class.

Two consequences of our coding denominator are worth stating plainly. It contains the initiation
peak, which inflates it and deflates every ratio built on it. And it contains the termination peak,
which our `termination_occupancy` also counts in its numerator — so that quantity is a part measured
against a whole that includes it, and a real fall in termination density shrinks its own denominator.
Neither is fatal and neither invalidates ADR-0010, but both are differences from the published
analysis of the same libraries, and they were not chosen.

## Decision

**The coordinates, with their endpoints.** riboWaltz counts `psite_from_start` 0 at the first
nucleotide of the coding sequence and `psite_from_stop` 0 at its last, and `l_cds` includes the stop
codon. All three windows below are inclusive at both ends, and they are disjoint by construction:

```
coding denominator   psite_from_start >= 18   and   psite_from_stop <= -16
termination          psite_from_stop  in [-15, 0]
extension            psite_from_stop  in [1, extension - 3]
```

The coding denominator excludes sixteen nucleotides at its 3' end rather than the authors' fifteen,
because our termination window is the sixteen positions `[-15, 0]`. Fifteen would leave one position
counted in both, which is the kind of overlap that hides in a range written as `-15..0`.

**The coding denominator excludes the initiation and termination peaks.** This changes the primary
assay, not just its documentation. ADR-0010 divides by every coding P-site, so the denominator
carries the initiation peak, which inflates it and damps all three quantities, and carries the
termination peak, which `termination_occupancy` also counts in its numerator — a part measured
against a whole containing it, where a real fall shrinks its own denominator. Reading the authors'
code is what surfaced it.

It is changed now because GSE144140 has not been calibrated and no number from it has been seen. The
cost is stated plainly: the HEK293T exploratory figures were computed without insets and are **not
comparable** to anything produced under this record. They are not restated, and HEK293T is
recalibrated under this rule before any of it is repeated.

**The length window is selected per dataset, once, from periodicity alone.** One shared set of
lengths is chosen for all libraries of a dataset — not per library, which would make libraries
incomparable within a contrast. The first calibration pass surveys 18 to 40 nucleotides, wide enough
to contain the empty-E-site, empty-A-site and full-length populations the authors of GSE144140
resolve at 18–19, 21–24 and 28–35. A length joins the set when, in every library of the dataset, its
coding-frame-0 share over the coding denominator above exceeds each of the two off-frames and the
length carries at least 1% of that library's assigned P-sites. With the depth floor below, that 1%
is at least ten thousand P-sites for every retained length in every library.

**The selection is blind to the readthrough result.** It is computed from `frame_by_length.tsv` and
the metagene profiles, both of which exist before any contrast is run, and it is written into this
record before the assay runs. A window chosen after seeing which window produces a pass is not a
window, it is a result.

**Every library must pass on the selected set, or the dataset is inconclusive.** Predeclared: at
least 1,000,000 assigned P-sites per library over the selected lengths; an inferred 5' offset of 11
to 14 nucleotides on the dominant length; and coding-frame-0 share above 40% pooled over the selected
set. A library that fails is not dropped and the thresholds are not moved. The dataset is reported as
inconclusive, exactly as the misspecified frame control was.

A failed calibration is not repaired by borrowing the authors' fixed window. That substitution would
be a length set chosen because the predeclared one did not pass, which is the whole thing this record
exists to prevent. Inconclusive is the answer.

**A calibration manifest gates the assay.** Calibration writes the selected lengths, the inferred
offsets, the per-length periodicity, the per-library QC verdict and the checksum of the script that
produced them. The readthrough command requires that manifest and refuses to run when it is missing
or records a failure. Nothing about library quality currently reaches the assay, and a library with
no periodicity would flow into a frozen contrast without a word.

**The second arm is a published-parameter sensitivity analysis, and is not called parity.** GSE144140
is analysed again under the authors' numbers — lengths 28 to 35, coding denominator inset by 18 and
15, downstream window opening 6 nucleotides past the stop — and both are reported side by side.
That is our estimator run with their parameters. It is not a reproduction of their result, and the
word parity is reserved for something this does not do: run their code-defined estimator entire, with
their transcript-level calculation, their aggregation, their annotation and reference assumptions,
and a comparison against the count tables deposited with the series. Until that is done, agreement
here is agreement about parameters, not about implementations.

RiboRescue's assay is primary. The sensitivity arm cannot replace it, and cannot rescue it: a
primary result that fails is not repaired by a parameter set under which it would have passed. If
the two disagree, the disagreement is the finding and is reported as one.

**The adapter survey stops deciding whether a footprint length is acceptable.** It measures the
distance from a read's start to its adapter, so a read too short to show enough adapter is not
counted at all: with a 50 nt read, a 6 nt minimum overlap and this chemistry, nothing above about 34
nucleotides is observable. The survey reports that ceiling per library beside the distribution.
Its lower bound on plausible length is removed, because censoring is exactly what would trip it.

Its upper bound stays. Censoring can only pull a measured length down, never up, so a median far
above the plausible range is never an artifact of it — and that is the check that catches a series
whose declared adapter is present in the reads but has a profiling linker in front of it, which is
the case this survey exists for.

## Consequences

Calibration runs once more per dataset after the lengths are chosen, because the choice needs a first
pass to look at.

HEK293T was calibrated at 26 to 34 and its exploratory result stands under that window. It is
recalibrated under this rule before any of it is repeated, and the two are not mixed.

The reasoning that a treatment-correlated length shift reflects empty-A-site ribosomes is withdrawn.
That population is 21 to 24 nucleotides in these authors' hands, and a median moving from 27 to 26
does not reach it. The confound stands on the differential discard alone, which needs no mechanism.
