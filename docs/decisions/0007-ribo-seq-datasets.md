# ADR-0007 — The Ribo-seq datasets

**Status:** accepted · **Date:** 2026-07-21 · **Deciders:** Joseph, Mahan

## Context

Two things need ribosome profiling, and they pull in different directions. The course requires raw
FASTQ from three samples carried through quality control, trimming and alignment. The safety layer
requires the readthrough signature — ribosome density falling at termination codons while 3'UTR
density rises — which only appears in a compound-treated library judged against an untreated one,
with enough replication to distinguish the effect from noise.

A dataset serving both must be human, deposit raw FASTQ openly rather than processed counts, and
carry matched RNA-seq so that absent ribosomes can be told apart from absent transcript.

## Decision

Both come from **PRJNA576648** (Wangen & Green, *eLife* 2020, `10.7554/eLife.52611`), in two arms:

- **Calu-6**, ribosome profiling and matched RNA-seq, untreated against G418. Calu-6 is homozygous
  for `TP53` p.Arg196Ter, so the course arm runs in a background carrying a pathogenic nonsense
  variant the scored table already contains.
- **HEK293T**, ribosome profiling only, three biological replicates each of untreated and G418.
  The Calu-6 arm is one library per condition, which cannot support the safety comparison; these
  can.

Twelve runs, 18 GB. `pipeline/assets/riboseq_samples.tsv` names them with the checksums ENA
publishes, and nothing is committed.

**The adapter is a declared property of each run, and the declaration is checked.** GEO records the
Ribo-seq linker for this project as `NNNNNNCACTCGGGCACCAAGGAC`. The Calu-6 reads do not contain it;
they carry `CTGTAGGCACCATCAAT` and no UMI on either end. Trimming against the recorded linker would
have removed nothing and passed 17 nt of adapter into alignment, where reads still map — softly
clipped — and periodicity quietly disappears. So the samplesheet carries the adapter per run and
`riborescue trim-summary` refuses any library whose declared adapter was found in under half its
reads.

`GSE179274` — ribosome profiling of `IDUA` p.Trp402Ter patient fibroblasts, with both a G418 and a
suppressor-tRNA arm — is the second dataset, once this one is through the pipeline. It is the only
public ribosome profiling of a suppressor tRNA, the modality the coverage ranking proposes, and
nothing substitutes for it. It has no RNA-seq and two replicates per arm, so it supplements rather
than replaces the choice above.

## Consequences

The index is built sampling every second suffix array entry. That halves what an alignment holds in
memory, so two samples align at once rather than one; each search is slower and the pair is not.

Two thirds of a footprint library is structural RNA and it aligns to many loci rather than none, so
it is depleted before the genome sees the reads. GENCODE annotates almost no rRNA — the rDNA arrays
are not in the primary assembly — so the repeating unit comes from its own accession, `U13369.1`.
Selecting GENCODE biotypes alone removes a tenth of the library; with the repeat unit, three fifths.

Read length is not uniform across the project — 50 nt and 100 nt single-end footprints, 2×100
paired transcriptome — so trimming is parameterised per run rather than set once.

`RPF_` and `RNA_` in the sample title are the only reliable way to tell the assays apart: SRA
labels every run in the project `RNA-Seq`, ribosome footprints included. The assay is recorded in
the samplesheet rather than derived from the archive.
