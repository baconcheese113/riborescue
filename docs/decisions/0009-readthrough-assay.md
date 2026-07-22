# ADR-0009 — The readthrough assay, fixed before the answer is known

**Status:** accepted · **Date:** 2026-07-22 · **Deciders:** Joseph, Mahan

## Context

Gate 5 asks whether these libraries can detect stop-codon readthrough at all. It is a sensitivity
check on the assay, run against a positive control, and it has to pass before any null result for a
suppressor tRNA means anything: an assay that cannot see G418 working cannot be trusted to report
that something else is not.

The measurement is easy to talk oneself into. Ribosome density past a stop codon rises for reasons
that are not readthrough — a downstream gene, an unannotated exon, a short untranslated region
carrying spillover from the termination peak — and a metric tuned after seeing the treated samples
will find whichever of those makes the answer come out. So the analysis is fixed here first.

The plan was written after the periodicity figure was made, which necessarily shows both conditions
and in which the treated stop-codon peaks look smaller. No readthrough quantity has been computed.
Nothing below is changed once numbers exist; if the definition turns out to be wrong, that is
recorded as a failure and a new decision, not edited into a success.

## Decision

**Readthrough is a redistribution, and both halves are required.** G418 should lower ribosome
occupancy at the termination site *and* raise in-frame occupancy downstream of it. Either alone is a
different phenomenon: raised density at the stop without downstream signal is stalling, which is
what SRI-37240 does, and downstream signal without a fall at the stop is more likely an annotation
problem than a decoding one. A library passes only with both.

**The unit of inference is the biological replicate, not the read or the gene.** HEK293T carries
three untreated and three treated libraries, so n = 3 per condition. Reads within a library are not
independent of each other and genes within a library are not either; treating them as replicates
would manufacture significance from depth.

**The comparison is paired within replicate.** The three replicate experiments were prepared
differently — rep2 uses one linker chemistry, rep1 and rep3 another — and preparation drives
coding-frame occupancy from 42% to 64%, a spread far wider than any treatment effect is likely to
be. Comparing treated against untreated across the pooled set would confound preparation with
treatment. Each replicate is compared against its own partner, and the three differences are what
carry the claim.

**The quantity, per transcript per library.** In-frame P-sites in the 3' untranslated region,
divided by in-frame P-sites in the coding sequence of the same transcript. Normalising within a
transcript removes library depth and transcript abundance, neither of which is measured here.

**Out-of-frame downstream occupancy is the negative control.** Readthrough continues the reading
frame, so it must raise the in-frame fraction specifically. A treatment that raises in-frame and
out-of-frame occupancy equally has changed something else — degradation, contamination, mapping —
and the claim fails.

**Transcripts are excluded before the comparison, for reasons that do not involve treatment:**

- annotated programmed-readthrough genes, which read through natively;
- transcripts whose 3' untranslated region overlaps a downstream coding region on the same strand;
- untranslated regions shorter than 50 nt, where the termination peak spills into the window;
- transcripts without enough coding-sequence coverage to give a stable denominator.

**What is reported.** The per-replicate paired differences, their consistency, and a confidence
interval — not a pooled read count and not a p-value standing alone.

## Consequences

The HEK293T libraries have no matched RNA-seq, so a transcript's abundance is unknown and cannot be
controlled for. The result is therefore a statement about **where ribosomes sit on a transcript**,
relative to that transcript's own coding sequence, and not about absolute translation or about
transcript level. Calu-6 has matched RNA-seq but is Mycoplasma-contaminated in a way that tracks
treatment, so it cannot repair this.

Three replicates support a direction and a rough interval, not a precise effect size.

Passing Gate 5 establishes that the assay can see readthrough when it is present. It does not
establish the safety burden of any therapy, which needs the native-stop atlas across many
transcripts and is a later piece of work.
