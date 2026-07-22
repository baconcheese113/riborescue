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

**The statistic is pooled within a library, not a median across transcripts.** A frame composition is
a proportion, and the frozen run showed the per-transcript median sitting at exactly zero in all
three untreated libraries. Pooling the qualifying transcripts of a library gives a share that is
defined even where most transcripts carry nothing. The per-transcript median and the share of
transcripts with any downstream signal are both reported beside it, so a pooled figure resting on
few transcripts is visible rather than hidden. This choice was made after seeing the frozen run,
which is one more reason the rule is confirmed elsewhere before it is believed.

**The pass rule, in full.** Writing `frame gap` for a library's downstream frame-0 share minus its
own coding frame-0 share, the control passes only when all three hold:

1. **Downstream occupancy rises.** The paired difference in downstream in-frame occupancy is
   positive in every replicate and positive on average.
2. **Termination occupancy falls.** The paired difference in termination occupancy is negative in
   every replicate and negative on average.
3. **Frame composition moves toward the coding frame.** The paired difference in `frame gap` is
   positive in every replicate, positive on average, **and its 95% interval excludes zero.**

The interval is required to exclude zero for the third condition only. It is the claim the control
exists to make, so it carries the stricter bar. The first two are directional requirements, where
agreement across three replicates is the evidence and a three-point interval is too wide to add
much. A result meeting one or two conditions is reported as what it is and does not pass.

**This respecification is not tested on the libraries that motivated it.** Everything above was
chosen while looking at those six libraries, so their reanalysis is exploratory and stays labelled
that way.

**Calu-6 is a mechanical check, not a confirmation.** One treated library against one untreated one
gives no replication, and its Mycoplasma burden tracks treatment, so it can show that the
computation runs and produces sane numbers on a second cell line. It cannot establish biological
specificity and is not counted toward the rule.

**The supporting contrasts are named here, and all of them are reported.** GSE179274 holds patient
fibroblasts as five separate libraries of two replicates each — untreated, G418 at 0.1 mg/mL, G418
at 0.5 mg/mL, an EGFP control construct, and a tyrosine suppressor tRNA — and a mouse liver arm of
three treated against three control. The contrasts are fixed as:

- untreated against G418 0.1 mg/mL;
- untreated against G418 0.5 mg/mL;
- EGFP construct against tyrosine suppressor tRNA;
- mouse liver control against treated.

Each is analysed and reported separately. Whichever looks best is not selected afterwards, and the
mouse arm needs its own annotation and extension windows before it can be run at all.

**No pairing is invented there.** Each treatment in GSE179274 is its own library, with the two
replicates sitting inside it, so treatment and library preparation are confounded by the design and
`rep1` of one arm has no relationship to `rep1` of another. Those contrasts are therefore unpaired
with two replicates per arm — weak, and able to corroborate a direction rather than carry the claim.
That the statistic is normalised inside each library against its own coding sequence is what makes
an unpaired comparison tolerable at all.

**The quantities, exactly.** Within one library, over the qualifying transcripts of the contrast's
shared universe, counts are summed before any ratio is taken:

```
downstream_occupancy = Σ extension_frame0 / Σ cds_frame0
termination_occupancy = Σ termination      / Σ cds_frame0
downstream_share0     = Σ extension_frame0 / Σ (extension_frame0 + extension_frame1 + extension_frame2)
cds_share0            = Σ cds_frame0       / Σ (cds_frame0 + cds_frame1 + cds_frame2)
frame_gap             = downstream_share0 - cds_share0
```

Paired difference for a replicate is the treated library's value minus its partner's. Unpaired
difference is the mean over treated libraries minus the mean over untreated ones.

**What can confirm, and what can only support.** The rule is defined for a paired design with three
replicates, so only a dataset with that structure can confirm it. That is **GSE144140**: HEK293T,
DMSO against G418 500 µg/mL against SRI-37240 10 µM, three replicates of each, and no part in
shaping this rule, and none of its runs appears in the series that did. It carries its own negative
control — SRI-37240 raises occupancy at stop codons without raising occupancy beyond them, which is
stalling rather than readthrough — so the contrast discriminates instead of merely responding.

The control is confirmed when G418 against DMSO meets all three conditions **and** SRI-37240 against
DMSO fails the signature in the direction its mechanism predicts: termination occupancy does not
fall, and downstream occupancy does not rise. Failure is required of the signature as a whole rather
than of the frame condition specifically, because a frame composition built from sparse downstream
counts is unstable, and a compound that produces almost no downstream signal would fail that
condition for arithmetic reasons that say nothing about stalling. A rule that fires on both
compounds has not demonstrated specificity.

**The comparison there is unpaired.** Its replicate letters run across the treatments (`1_dmso_A`,
`4_g418_A`, `7_sri37240_A`), which is consistent with matched batches, but the series metadata does
not say so. Absent an explicit statement of blocking, the primary analysis is the unpaired
three-against-three rule below. The paired analysis is reported as a pre-specified sensitivity
check, and it is reported whatever it shows: which of the two is primary is fixed here and is not
chosen later according to whether they agree.

For an unpaired contrast, conditions 1 and 2 require the difference in means to have the stated sign
with every treated library on the correct side of every untreated one, and condition 3 requires in
addition that the Welch interval on `frame_gap` excludes zero.

GSE179274 supports and cannot confirm. Its arms are separate libraries of two replicates each, so
its contrasts are unpaired at n = 2, where an interval is too wide to decide anything. They are
reported for direction and for the suppressor tRNA, which is the modality this project ultimately
cares about and which no other public dataset profiles.

**The unpaired rule, for the supporting contrasts.** The estimator is the difference in means, the
interval is a Welch t interval on the two groups, and consistency means every treated library
exceeds every untreated one for that quantity. A supporting contrast is described as agreeing with
the control, disagreeing with it, or indeterminate. It never converts into a pass.

**Results are kept apart by dataset.** Each dataset writes to its own directory and its own
combined tables, and the combiner reads only the directory it was given. Files from HEK293T,
Calu-6, the fibroblast arm and the mouse arm are never merged by a wildcard, which is how stale
per-library tables from an earlier run reached a combined table once already.

## Consequences

The extension counts are recomputed with three frame columns, which is another pass over the
alignments.

Two claims are now separated and must not be merged: that the assay can see readthrough, and that
this dataset shows it. The first is what the control is for; the second needs the corrected rule
confirmed somewhere it was not designed.

Calu-6 carries the same drug and is Mycoplasma-contaminated in a way that tracks treatment, so it
tests the rule's mechanics rather than settling the biology.
