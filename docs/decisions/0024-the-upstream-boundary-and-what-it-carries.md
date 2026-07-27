# ADR-0024 — The upstream boundary, and what it is allowed to carry

**Status:** accepted · **Date:** 2026-07-26 · **Deciders:** Joseph

## Context

ADR-0003 put `nf-core/riboseq` outside this repository and made `UpstreamHandoff` the only door
through which its outputs arrive. What happened instead is that a local `READS` workflow grew —
trimming, contaminant depletion, alignment, quality control — and every frozen result was produced
through it. The handoff has never fed an analysis.

The schema drifted the same way. It declares `psite_offsets`, `codon_coverage` and `cds_coverage`,
which are riboWaltz tables `nf-core/riboseq` genuinely publishes. Nothing reads them, and nothing
can:

- `codon_coverage()` and `cds_coverage()` take no read-length argument, so their output is summed
  over footprint lengths. ADR-0011 chooses which lengths a dataset keeps *from the periodicity the
  same pass measures*, and ADR-0013's sensitivity arm re-sums the same counts over the published
  window. A table that has already collapsed lengths cannot answer either question.
- The readthrough assay counts P-sites in the extension window between a transcript's own stop and
  the next in-frame stop, split by frame, and in the termination window before it. No upstream
  output carries those at all.

So `run_psite.R` is not a convenience that upstream could absorb. It is the one pass that produces
the length-stratified, extension-window counts every readthrough and occupancy result rests on.

## Decision

**`nf-core/riboseq` owns alignment and quality control for sequencing acquired from here on.**
Trimming, depletion, alignment and the QC report are a solved problem and are not this project's
work.

**`run_psite.R` stays local, and stays the single custom measurement pass.** It reads alignments and
writes the length-stratified counts, which is the part no external pipeline emits.

**The frozen corpus keeps the path it was measured on.** ADR-0012 closed the sequencing scope with
four series already aligned. Re-running them through a different preprocessor would reproduce inputs
already in hand and change no result, so `READS` remains until the corpus it produced is superseded.

**A handoff contract declares only outputs something reads.** `psite_offsets`, `codon_coverage` and
`cds_coverage` are removed from `UpstreamHandoff`, from the committed manifest fixture, and from the
fixture tree. A field is added back when a consumer exists, not in anticipation of one.

Removing the offsets deserves its own reason, because upstream's are per length and could in
principle be read. `run_psite.R` has to run regardless, and it infers offsets as part of the same
pass that produces the counts. Consuming upstream's as well would leave two offset estimates for the
same libraries, from different assignments of reads to transcripts, with nothing saying which
governs — the failure ADR-0022 has already cost this project once.

## Consequences

The handoff now names four outputs: the alignments, the merged RNA-seq counts and TPM, and the
MultiQC report. None has a consumer today either, because the local path produces its own; they stay
because they are what the boundary carries when it is used, and because each is readable as it
stands rather than needing a shape no upstream tool emits.

`READS` and the conventions that call the pipeline minimal both still describe a repository that
authors little of its own preprocessing, which is not what it does. That mismatch is recorded here
rather than repaired: repairing it means running the next dataset through `nf-core/riboseq`, and
that is the change this record authorises rather than performs.
