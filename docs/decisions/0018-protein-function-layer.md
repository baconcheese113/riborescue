# ADR-0018 — The protein-function layer: conservation and domain context, kept apart

**Status:** proposed · **Date:** 2026-07-24 · **Deciders:** Joseph

## Context

The NMD layer answers whether a premature stop leaves transcript to rescue; it says nothing about
whether the protein position that stop truncates *matters*. Two variants that both escape decay and both
read through are not equal if one sits in a disordered linker and the other in a catalytic domain at a
residue conserved across vertebrates. The PRD's next function layer is protein-function context —
conservation, domain and active-site membership, and, later, a learned missense-tolerance model
(ESM-1v). This records how the first two are built, before any of their numbers are looked at, and why
they stay separate fields rather than a single score.

Two decisions frame the work. First, **ADR-0005's restraint carries over**: conservation and domain
membership are evidence a curator weighs, not a priority number the tool computes. A residue being
conserved does not rank a variant above another; it is shown, with its source, beside the readthrough
and NMD layers. Second, **the layers must not be blended.** A conserved position, a domain hit, and a
learned tolerance score are three different measurements with three different failure modes; averaging
them into one "functional impact" number would hide which one is driving a call and manufacture a
confidence the inputs do not support.

## Decision

**Conservation — phyloP and phastCons at the stop, from a pinned UCSC track.** For each nonsense
variant, read the per-base phyloP (evolutionary rate) and phastCons (conserved-element probability) at
the variant's genomic position and across its stop codon, from the UCSC hg38 100-way vertebrate tracks.
The score carried is the codon-mean and the single-base value, both reported — a highly conserved stop
position is evidence the residue lost is constrained.

- **Source and pinning.** The tracks are large bigWigs (`hg38.phyloP100way.bw` ≈ 9.6 GB,
  `hg38.phastCons100way.bw` similar). Two ways to consume them, to be decided with the download budget
  in view: (a) fetch once, pin by checksum through the inputs system, and read locally with pyBigWig —
  reproducible and offline, at the cost of ~20 GB on disk; or (b) read the ~70k needed positions from
  the UCSC-hosted bigWig over HTTP range requests, no bulk download, at the cost of a runtime network
  dependency that the provenance model otherwise avoids. **(a) is preferred** for consistency with the
  fetch-with-checksum ethos; (b) is the fallback if the download budget forbids it. Either way the
  score is a lookup, CPU-light, no GPU.

**Domain and active-site context — from a pinned protein-feature source.** For each variant's protein
position on its MANE transcript, report whether it falls inside an annotated domain, and whether it is
at or adjacent to an active/binding site, from UniProt's per-protein feature table (pinned release).
The field is categorical — domain name(s) and site membership — not a score.

**Three fields, never one.** `conservation` (phyloP/phastCons), `domain` (membership + names), and, when
it arrives, `esm1v` (learned missense tolerance) are separate columns in the exports and separate slots
in the viewer, each with its own provenance and caveat. No combination into a single functional number.

**Frozen before results.** The tracks, the release versions, and the aggregation (codon-mean phyloP,
single-base phastCons, domain membership by interval overlap) are fixed here; none is tuned after a
distribution is seen. A smoke over a handful of variants validates the schema before the full pass.

## Consequences

The patient card and researcher view gain a protein-function slot that says, for a given stop, how
conserved the position is and whether it sits in a domain — beside, not blended with, the readthrough,
residue, safety and NMD layers. The curator sees three independent signals and weighs them; the tool
does not weigh them for anyone. ESM-1v follows under its own record once conservation and domain are in,
because a learned score is a heavier dependency and a different kind of claim than a lookup. Nothing here
becomes a rank, and nothing here shares a field with a measurement of a different thing.
