# ADR-0010 — The frame control, respecified

**Status:** accepted · **Date:** 2026-07-22 · **Deciders:** Joseph, Mahan

## Context

ADR-0009 fixed the readthrough assay before its answer was known, which is why the flaw described
here was found rather than absorbed. That decision stands as written and is not edited; this one
supersedes its frame control.

The rule required the in-frame ratio to rise by more than the out-of-frame ratio in every replicate.
In the implementation, in-frame counts one phase of the reading frame and out-of-frame pools the
other two. Under a background that has no frame preference, the pooled quantity is expected to be
about twice the single one, so the rule asked a single phase to gain more than two phases combined.
That is not a test of frame specificity; it is a much stricter and differently shaped question, and
the G418 positive control failed it.

Reanalysis of the same libraries — exploratory, on data already seen, and therefore not evidence for
anything — indicates the misspecification mattered. Pooled across the qualifying transcripts, the
downstream frame-0 share is 32.8% to 34.0% in the untreated libraries, which is the value expected
if downstream density had no frame preference at all, and 39.8% to 52.5% in the treated libraries,
within two to three points of each library's own coding-sequence frame composition. The shape of
that is what readthrough should look like. It is not a result, because the rule that would judge it
was chosen after seeing it.

## Decision

**The G418 positive control is inconclusive, not negative.** The frozen assay failed its predefined
rule and that is recorded as a failure of the rule. A misspecified control cannot support a
biological negative any more than it could support a positive.

**Frames are counted and reported separately.** Frame 0, +1 and +2 in the extension window, per
transcript per library, rather than one phase against two pooled.

**The null is the library's own coding sequence, not one third.** Coding-frame occupancy ranges from
42% to 64% across these libraries, so a fixed expectation of one third understates the bar in some
and overstates it in others. Readthrough continues translation, so its downstream frame composition
should resemble the composition upstream in the same library. The quantity is the downstream frame-0
share minus that library's coding frame-0 share, compared treated against untreated, paired within
replicate as before.

**Zero-heavy measurements are handled explicitly.** In the frozen run all three untreated medians
were exactly zero: the median qualifying transcript had no downstream P-sites at all, so each paired
difference was a treated value minus a floor. The respecified assay reports the share of transcripts
with any downstream signal, states the pooled counts beside the per-transcript medians, and does not
present a median at the floor as though it were a measurement.

**This respecification is not tested on the libraries that motivated it.** Everything above was
chosen while looking at those six libraries, so their reanalysis is exploratory and stays labelled
that way. The rule is confirmed on data that took no part in shaping it — the Calu-6 pair, and the
suppressor-tRNA libraries of GSE179274 — before the control is called passed.

## Consequences

The extension counts are recomputed with three frame columns, which is another pass over the
alignments.

Two claims are now separated and must not be merged: that the assay can see readthrough, and that
this dataset shows it. The first is what the control is for; the second needs the corrected rule
confirmed somewhere it was not designed.

Calu-6 carries the same drug and is Mycoplasma-contaminated in a way that tracks treatment, so it
tests the rule's mechanics rather than settling the biology.
