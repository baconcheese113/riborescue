# ADR-0019 — What a refused dataset may still say: the detectability arm

**Status:** accepted · **Date:** 2026-07-25 · **Deciders:** Joseph

## Context

GSE179274 is the only public dataset that profiles a suppressor tRNA, which is the modality this
project ultimately cares about, and its primary analysis has concluded: five of its ten libraries
failed the calibration ADR-0011 fixed before they were looked at, so `read_manifest` refuses the
dataset and no contrast runs on it. That verdict is correct and stands.

It is also unsatisfying in a specific way. A refusal says the dataset cannot answer the question. It
does not say *what would*. The libraries exist, the calibration ran, and the reason each library
failed is on record — which is enough to state how large an effect this design could have resolved,
and what a design that resolved a smaller one would have to look like. That statement is worth more
than the contrast would have been, because it is not a claim about readthrough at all. It is a claim
about measurement, and its uncertainty is the point rather than a caveat on the end.

**The two failure modes are not the same kind of thing, and only one of them tolerates this.**

| Library | P-sites | 5′ offset | Fails on |
|---|---|---|---|
| `fibroblast_egfp_rep1` | 1,519,474 | 12 | — |
| `fibroblast_egfp_rep2` | 659,427 | 12 | depth |
| `fibroblast_suptrna_tyr_rep1` | 999,855 | 12 | depth |
| `fibroblast_suptrna_tyr_rep2` | 815,668 | 13 | depth |
| `fibroblast_untreated_rep1` | 1,938,153 | 10 | offset |
| `fibroblast_untreated_rep2` | 2,351,378 | 9 | offset |

A depth failure is a **precision** failure. Fewer P-sites widen the interval around a quantity that
is still the quantity intended; nothing is being measured wrongly, only imprecisely, and an interval
is the honest and complete description of it.

An offset failure is a **validity** failure. The inferred P-site offset is how a read is assigned to
a codon, so an offset outside the window means the counts are attributed to the wrong position. That
does not widen an interval, it displaces the estimate, and a confidence interval around a displaced
estimate is a precise statement about the wrong quantity. No amount of sequencing repairs it and no
uncertainty measure represents it.

So the suppressor-tRNA contrast — whose four libraries all calibrate to an in-window offset and fail
only on depth — admits this treatment, and the G418 contrasts, which rest on untreated libraries at
9 and 10 nt, do not.

## Decision

**An exploratory detectability arm is run, and it estimates how large an effect this design could
have resolved. It is not a contrast and does not report one as a finding.** Everything below is
fixed before any GSE179274 effect is computed.

**Scope: the EGFP against tyrosine-suppressor-tRNA contrast only.** The two G418 contrasts of
ADR-0010 are reported as not analysable, naming the offset of each untreated library as the reason.
They are not run in this arm under any window, and the absence of a number for them is a result of
this record rather than of anything observed.

**A library whose calibration failure is validity is refused here as it is refused in production.**
The classification comes from the manifest's own fields rather than from its failure text: an offset
outside `OFFSET_FROM_5` or a coding frame-0 share below `MINIMUM_FRAME0_SHARE` is validity; a P-site
count below `MINIMUM_PSITES`, alone, is precision. The arm accepts precision failures and nothing
else, so it cannot be pointed at a dataset that failed for a reason an interval does not cover.

**Both length windows are run and both are reported.** The selected set the manifest records
(28, 29, 31 nt) and the published 28–35 nt window, exactly as the sensitivity arm of the passing
datasets does. Neither is chosen after the fact, and a disagreement between them is a result about
the window rather than an occasion to prefer one.

**The quantities are the ones ADR-0010 already fixes** — `downstream_occupancy`,
`termination_occupancy`, `frame_gap` — computed by the same pooled estimators over the same
qualifying universe. Nothing about the assay is redefined here.

**The reference effect is GSE144140's G418 contrast, per quantity.** A minimum detectable effect
means nothing without something to detect, and the only pre-declarable size is one measured
elsewhere: the HEK293T G418-against-DMSO differences from a dataset that passed its calibration and
completed its analysis (`downstream_occupancy` +0.00954, `termination_occupancy` −0.01147,
`frame_gap` +0.18301, selected set; the published-window figures are used for the published-window
arm). The question the arm answers is therefore fixed and concrete: *could this design have seen an
effect the size of the one G418 produces in HEK293T?* Using a suppressor tRNA's effect size, which
is unknown, would make the question circular; using the fibroblast effect would make it a
post-hoc rationalisation of whatever was found.

**Uncertainty is decomposed, because the two components buy different experiments.** Per library and
per quantity the arm reports:

- a **counting standard error**, computed analytically from the counts the ratio is built from — a
  Poisson error on an occupancy, a binomial error on a frame share. This is the component that falls
  as depth rises, and it is the only component that does;
- a **transcript-bootstrap interval**, resampling the library's qualifying transcripts with
  replacement. This is the honest within-library precision, wider than the counting error because it
  also carries transcript-to-transcript heterogeneity, and it is what a per-library value is quoted
  with.

Across the arm's libraries the observed variance is split into a counting part and a residual
**between-library** part, which no depth removes and only replication does. Where the residual comes
out negative — entirely possible at two libraries per arm, where the variance estimate has one
degree of freedom — it is reported as zero and flagged, not hidden.

**The minimum detectable effect is stated at 5% two-sided significance and 80% power**, on the
difference in means between unpaired arms of the sizes the dataset has. Its error combines the two
arms' variances without assuming they are equal, as ADR-0010's Welch interval does, but its
reference distribution takes the design's degrees of freedom — libraries in both arms, less two —
rather than Welch's approximation. That differs from how ADR-0010 reports an observed contrast, and
deliberately: Welch's degrees of freedom fall as the arms' variances diverge, so a hypothetical
depth would move the reference distribution as well as the error, and a design could come out
*less* able to resolve an effect at unlimited depth than at the depth it has. How many libraries
were sequenced is a property of the design and does not change with how deeply each was read.
ADR-0010's interval is unchanged; it reports a contrast, and this sizes a design.

Two further figures follow from it and are the arm's deliverable:

- the **depth multiplier** — the factor by which P-site depth would have to rise, at these arm
  sizes, for the reference effect to become detectable, or `unreachable` where the between-library
  component alone already exceeds it;
- the **replicates required** — the libraries per arm needed at the observed depth.

The depth multiplier is a lower bound and is labelled one. It scales the counting component as
1/√depth, which is right for that component, while the transcript heterogeneity inside a real
library's precision does not shrink with reads at all; the true requirement is therefore at least
this and probably more.

**At two libraries per arm the between-library variance is estimated from one degree of freedom, and
every figure resting on it is reported with that stated.** These numbers describe the design that
was run. They are not precise, and the arm's own headline is a range of what would be required
rather than a threshold to be met.

**The output is exploratory and structurally kept apart.** It writes under `results/exploratory/`,
the command refuses any destination outside it, and no production task reads that tree. Nothing from
this arm reaches the scored variant × therapy table, the coverage frontier, the web export or any
figure that is not itself labelled exploratory.

**The vocabulary is fixed.** This arm produces neither a pass, a confirmation, a validation nor a
negative, whichever way any number falls. GSE179274 remains what ADR-0010 calls it — independent
supporting evidence, which cannot complete a rule defined for three paired replicates — and what
ADR-0011 refused. Its verdict in the falsification ledger is `untestable`, with the depth and
replication that would make it testable recorded beside it.

**The primary manifest, the calibration thresholds and ADR-0013 are untouched.** No threshold is
lowered, no offset window is widened, and the unanimity amendment ADR-0013 anticipates is not
written here. This arm runs beside the refusal, not instead of it, and it exists because the refusal
holds.

## Consequences

A refused dataset acquires an output, which is the failure mode this record most has to guard
against: a table of numbers derived from GSE179274 now exists, and a reader who meets it out of
context could mistake a detectability estimate for a measurement of readthrough. The directory name,
the column names, the refusal to write elsewhere and the fixed vocabulary above are what stand
between those two readings, and they are load-bearing rather than decorative.

The reference effect ties this arm to GSE144140, so if that dataset's analysis is ever revised, the
detectability figures move with it and must be recomputed rather than quoted from here.

The same machinery applies to any future dataset that fails on depth alone, which is the more useful
consequence: what this project can say about a dataset it cannot use is now a computation rather
than a judgement.
