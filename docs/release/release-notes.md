# RiboRescue 0.1.0 — release notes

RiboRescue matches a patient's nonsense variant to candidate readthrough therapies. It is a Nextflow
pipeline that processes ribosome-profiling and RNA-seq libraries, a Python package that builds and
serves a scored variant × therapy table, and a static web app that presents it.

## What this release contains

**A scored table.** 422,256 variant × therapy predictions with confidence intervals, over 70,376 of
the 71,378 pathogenic nonsense variants in ClinVar release 20260715, for six readthrough compounds.
The model is a reproduction of Toledano *et al.* 2024, validated fold-by-fold against the authors'
own R implementation and cross-validated grouped by gene (median held-out R² 0.58–0.89 by drug).

**A confirmed frame control.** On GSE144140, G418 completes the pre-registered three-part ribosome
signature — downstream occupancy rises (+0.00962, CI 0.00505–0.01418), termination occupancy falls
(−0.01160), and downstream reads move into the coding reading frame (+0.18059, CI 0.11964–0.24154) —
while SRI-37240 fails it in the direction its mechanism predicts. 486 transcripts, three replicates
per arm, endpoint fixed in ADR-0010 before the contrast was computed.

**A native-stop safety atlas.** Downstream occupancy past *normal* stop codons under G418 in
HEK293T, as the cost side of readthrough.

**An NMD layer with four predictors** — a guideline rule, the fuller Lindeboom rule set, `aenmd`,
and NMDetective-AI — reported separately with their disagreements, never merged into a consensus.

**A suppressor-tRNA panel and coverage frontiers**, with per-disease coverage and a Pareto frontier
over experiment programmes.

**A read-processing pipeline** — FastQC, cutadapt, contaminant depletion, STAR, MultiQC,
featureCounts — over 31 archived libraries, with per-run checksum verification.

**A web viewer** with three audience-specific pages, deployed as a static export.

**A provenance manifest** (`results/provenance.json`): the commit, the lockfile hash, every input
with its checksum, pinned reference releases, container tags, the command behind each output, and an
explicit list of what could not be obtained.

## Results that did not go the way the hypothesis did

These are part of the release, not omissions from it.

**Ribosome kinetics adds no transferable information beyond local sequence.** The best held-out gain
over the reproduced baseline was +0.00356 R² (SJ6986). It beat a fully scrambled table (p = 0.020)
and beat gene composition alone (p = 0.005), but against the decisive control — permuting occupancy
within synonymous codon families, which preserves the amino acid and destroys only the tRNA-level
decoding demand the hypothesis names — the observed value sits inside the null at **p = 0.21**. The
gain is carried by amino-acid identity, not decoding speed. Not supported under the pre-registered
rule. Failing to reject that null is not equivalence; no equivalence margin was declared.

**The suppressor-tRNA dataset is inconclusive, not negative.** GSE179274 is the only public ribosome
profiling of a suppressor tRNA. Five of ten libraries fall below the pre-declared depth floor or
outside the P-site offset window, including both suppressor-tRNA libraries and their EGFP control.
The thresholds were not moved. The load-bearing frame-gap condition has a minimum detectable effect
of 0.393 against a G418-sized reference of 0.181, and 92–98% of the variance is between libraries
rather than in the counts — so the design could not have answered the question whatever it showed.
More replicates would reach it; more sequencing would not.

## Boundaries this release does not cross

- **Therapy scores are research rankings, not clinical recommendations.** Nothing here establishes
  that a compound helps a patient.
- **The G418 validation is within one laboratory and protocol family.** It is not independent
  replication, and it licenses no null about a therapy of different modality or data quality.
- **Native-stop occupancy is not protein production, toxicity, or clinical safety.** It is measured
  for G418 in HEK293T only.
- **Base-editing reachability is exploratory geometry.** Not editing efficiency, off-target
  activity, delivery, tissue access, splice consequence, or eligibility. Its negative class is named
  "not base-editable under the declared panel", never "no route exists" — prime editing is
  unevaluated. See `docs/release/adr-0023-scope-audit.md`.
- **Layers are never multiplied into one score.** They have no shared ground truth, so each is
  validated against its own evidence and reported as its own column.

## Known limitations

- **predNMD is unavailable** — no public release; requested from the authors and not obtained. The
  NMD ensemble ships without that tier.
- **Container images are pinned by tag, not digest.** Regeneration is reproducible to the tag, not
  to the byte.
- **The Calu-6 RNA-seq arm is Mycoplasma-contaminated in a way that tracks treatment** (35–42% of
  untreated reads, 13% of treated). G418 is an antibacterial. Every abundance comparison between
  those arms is confounded, and none is made. HEK293T libraries screen at 0.0% and the central
  result is untouched.
- **SRI-37240's replicate arms overlap.** The registered endpoint passes on means; a stricter
  replicate-separation diagnostic does not, and is reported beside the verdict rather than folded
  into it.
- **The dev server does not hydrate `/researcher` under WSL2.** The static export renders it in
  full, and the export is what ships. Use `pixi run app-preview`.

## Reproducing it

```
pixi install                     # resolve the toolchain from pixi.lock
pixi run check                   # lint + types + tests + controls
pixi run provenance              # the manifest every claim above is traceable through
```

Everything else is a named Pixi task, and the command behind each output is recorded in the
manifest.
