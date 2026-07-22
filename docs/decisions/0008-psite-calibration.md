# ADR-0008 — P-site calibration with riboWaltz

**Status:** accepted · **Date:** 2026-07-21 · **Deciders:** Joseph, Mahan

## Context

Gate 3 asks whether the footprints behave like ribosomes: read lengths near 30 nt, a stable P-site
offset per length, three-nucleotide periodicity and frame bias across the coding sequence. The
offset — the distance from a footprint's end to the codon in the ribosome's P-site — has to be
inferred per read length, and the frame has to be read in transcript coordinates, where an exon
junction does not interrupt the count.

The GLM is reimplemented in Python because it is served and must be a deployable artifact we control.
A P-site offset is not served; it is an upstream quantity the downstream code consumes. The PRD names
riboWaltz for exactly this step.

## Decision

Calibrate with **riboWaltz 2.0**, the published tool, rather than a reimplementation. It infers a
coherent offset for every footprint length from the pile-up of read ends at start codons — the
two-stage method its paper describes — and reports frame distributions and metaprofiles from the
same data.

riboWaltz pins **R 4.4**, one minor version behind the caret oracle's R 4.5, so it resolves in its
own locked Pixi environment (`psite`) rather than the default one. The calibration runs as a Pixi
task over the pipeline's footprint alignments, the same shape as the reproduction oracle.

Its input is a **transcriptome-coordinate alignment**, which STAR emits directly for the footprint
libraries with `--quantMode TranscriptomeSAM`. Placing the P-site in transcript space removes the
exon-junction arithmetic a genome coordinate would need, so the frame is read without a class of
silent error.

## Consequences

There are two locked R environments: R 4.5 for the caret oracle, R 4.4 for riboWaltz. Each is a
feature with its own solve, so neither constrains the other.

riboWaltz's annotation builder needs the Bioconductor data package `GenomeInfoDbData`, which ships
its payload in a conda post-link script rather than the package archive. Pixi disables post-link
scripts by default because they run arbitrary code, so the workspace enables them in a committed
`.pixi/config.toml` (`run-post-link-scripts = "insecure"`). This lets the install scripts of every
package in the workspace run outside a sandbox — a real widening of trust, bounded by installing
only from the committed lockfile and the pinned conda-forge and Bioconda channels, whose data
packages are checksum-verified before their scripts run.

The offsets come from the start-codon region and the frame bias is measured across the coding
sequence, so the signal that defines the offset is not the signal that tests it. A shift check
confirms the inferred offset maximises frame-0 enrichment against a one- or two-nucleotide neighbour.

Gate 3 is judged across all eight footprint libraries — stable offsets within a library preparation,
frame enrichment that holds — not from a single library. Low-support read lengths are reported rather
than forced to a verdict.
