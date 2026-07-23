# ADR-0011 — Footprint length selection, and parity with the authors of GSE144140

**Status:** proposed · **Date:** 2026-07-22 · **Deciders:** Joseph, Mahan

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

**The length window is selected per dataset, once, from periodicity alone.** One shared set of
lengths is chosen for all libraries of a dataset — not per library, which would make libraries
incomparable within a contrast. A length joins the set when, in every library of the dataset, its
coding-frame-0 share over internal coding sequence exceeds the two off-frames and the length carries
at least 1% of that library's assigned P-sites. Internal means excluding the first 18 and last 15
nucleotides of the coding sequence, so the initiation and termination peaks cannot decide which
lengths look periodic.

**The selection is blind to the readthrough result.** It is computed from `frame_by_length.tsv` and
the metagene profiles, both of which exist before any contrast is run, and it is written into this
record before the assay runs. A window chosen after seeing which window produces a pass is not a
window, it is a result.

**Every library must pass on the selected set, or the dataset is inconclusive.** Predeclared: at
least 1,000,000 assigned P-sites per library over the selected lengths; an inferred 5' offset of 11
to 14 nucleotides on the dominant length; and coding-frame-0 share above 40% pooled over the selected
set. A library that fails is not dropped and the thresholds are not moved. The dataset is reported as
inconclusive, exactly as the misspecified frame control was.

**A calibration manifest gates the assay.** Calibration writes the selected lengths, the inferred
offsets, the per-length periodicity, the per-library QC verdict and the checksum of the script that
produced them. The readthrough command requires that manifest and refuses to run when it is missing
or records a failure. Nothing about library quality currently reaches the assay, and a library with
no periodicity would flow into a frozen contrast without a word.

**The parity arm is predeclared and reported whatever it shows.** GSE144140 is analysed a second
time under the authors' parameters — lengths 28 to 35, coding denominator inset by 18 and 15,
downstream window opening 6 nucleotides past the stop — and both are reported side by side. If the
verdict depends on which parameter set is used, that dependence is the finding and is stated as such.
Their deposited count tables are the third comparison, and needing no reprocessing they are the
cheapest.

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
