# Numbers and citations, checked against their sources

Every quantitative claim the release documentation makes, and where it comes from. A number that
could not be traced to an artifact does not appear in the release documentation, and a number quoted
anywhere else should be taken from here rather than re-derived.

## Numbers

| Claim | Value | Artifact |
|---|---|---|
| Pathogenic/LP variants in ClinVar 20260715 | 343,459 | `data/clinvar/clinvar_20260715.vcf.gz`, counted on `CLNSIG` |
| Of those, annotated nonsense | 81,359 (23.7%) | same, counted on `MC` |
| Qualifying pathogenic nonsense variants | 71,378 | `results/clinvar_nonsense.tsv` |
| Scoreable (placed on MANE Select) | 70,376 | `results/amenability_landscape.tsv` |
| Not placeable | 1,002 | difference of the two above |
| Scored rows | 422,256 | `results/variant_therapy_scores.tsv` (70,376 × 6) |
| Therapies scored | 6 | `results/variant_therapy_scores.tsv` |
| Median held-out R², grouped by gene | 0.891 / 0.873 / 0.772 / 0.749 / 0.707 / 0.575 | `results/evaluation.tsv` |
| Reaching 1% readthrough, point estimate | 59,797 | `results/landscape_summary.tsv` |
| Reaching 1% readthrough, lower bound | 15,488 | `results/landscape_summary.tsv` |
| Transcripts in the readthrough analysis | 486 | `results/readthrough/gse144140/*_by_library.tsv` |
| G418 − DMSO downstream occupancy | +0.00962 (0.00505–0.01418) | `results/readthrough/gse144140/g418_vs_dmso.unpaired.tsv` |
| G418 − DMSO termination occupancy | −0.01160 (−0.02619–0.00299) | same |
| G418 − DMSO frame gap | +0.18059 (0.11964–0.24154) | same |
| SRI-37240 − DMSO, all three | −0.00024 / +0.00092 / −0.06284 | `results/readthrough/gse144140/sri37240_vs_dmso.unpaired.tsv` |
| SRI-37240 termination, per replicate | 0.0704, 0.0574, 0.0610 vs DMSO 0.0667, 0.0628, 0.0566 | `sri37240_vs_dmso.unpaired_by_library.tsv` |
| GSE144140 libraries, calibration failures | 9, 0 failures | `results/psite/gse144140/calibration.json` |
| GSE179274 libraries, calibration failures | 10, 5 failures | `results/psite/gse179274_fibroblast/calibration.json` |
| Frame-gap minimum detectable effect | 0.39272 | `results/exploratory/detectability/gse179274_fibroblast/suptrna_tyr_vs_egfp.selected.tsv` |
| Frame-gap reference effect | 0.18059 | same |
| Between-library variance share | 92.5% / 98.1% / 97.1% | same |
| Best kinetic gain (SJ6986) | +0.00356 | `results/kinetics/permutation_null/familywise.tsv` |
| Synonymous-family null mean, p | +0.00299, p = 0.21 | same, `shuffle = context_matched` |
| Global shuffle p | 0.020 | same, `shuffle = global` |
| Within-gene shuffle p | 0.005 | same, `shuffle = within_gene` |
| Permutations per family | 199 | same |
| Escape map: total / scoreable / unscoreable | 71,378 / 70,660 / 718 | `results/base_editing.tsv` |
| Escape map: exact / alternative / not editable | 7,622 / 14,420 / 48,618 | same |
| Escape map: reachable | 22,042 (31.2% of scoreable) | same |
| Every reachable variant's editor | ABE7.10; BE4max contributes 0 guides | same, `editor` column |
| Staged sequencing runs | 31 | `results/staged_runs.tsv` |
| Calu-6 RNA-seq raw reads | 21.75 M / 13.64 M / 16.30 M / 13.76 M | `results/reads/qc/course_qc_summary.tsv` |
| Calu-6 RNA-seq duplication | 87.7 / 82.5 / 73.1 / 70.5% | same |
| Calu-6 RNA-seq GC | 46 / 46.5 / 51 / 51% | same |
| Retention after cleaning | 99.19–99.62% (RNA-seq), 85.5–99.9% (Ribo-seq) | `results/reads/qc/trim_summary.tsv` |
| Adapter rate, footprint libraries | 88.6–99.6% | same |
| GSE144140 dmso rep1 base loss | 13.26 Gbp → 3.53 Gbp (73%) | same |
| Calu-6 RNA-seq overall mapped | 21.7 / 19.1 / 67.7 / 68.3% | `results/reads/qc/alignment_summary.tsv` |
| Unmapped-too-short, untreated pair | 78.3 / 80.9% | same |
| Footprint mean mapped length | 26.7–30.2 nt | same |
| Footprint multi-mapping | 17.3–41.7% | same |
| Depletion, calu6 untreated rep1 riboseq | 25.34 M → 10.53 M (58%) | `trim_summary.tsv` and `alignment_summary.tsv` |
| featureCounts assigned | 2.87 / 1.85 / 7.31 / 5.67 M | `results/rnaseq/counts.tsv.summary` |
| Top 20 expressed genes | MT-RNR2 125,651 mean TPM … | `results/rnaseq/top_expressed.tsv` |
| Library composition | mito 56.1/13.0/9.5/17.9% | `results/rnaseq/composition.tsv` |
| NMD scoreable / guideline escape / full-rules escape | 70,376 / 7,120 / 21,351 | `frontend/public/riborescue_research.json` |

### Not a result — exploratory and confounded

The four Calu-6 RNA-seq libraries show chaperone genes higher and translation-machinery genes lower
in the G418 arm (HSP90AA1 472 / 461 against 10,527 / 9,870, in
`results/rnaseq/top_expressed_nuclear.tsv`). **This is a descriptive difference between two ranking
tables, not differential expression.** Two replicates per arm, no test, no pre-declared analysis, and
the arms differ in GC content, usable depth and Mycoplasma burden — a live bacterial infection is
itself a stressor whose contribution cannot be separated here. It supports no biological claim, and
would need a pre-declared differential-expression design with adequate replication in a
Mycoplasma-free line before it could.

### Corrected while checking

- **422,316 → 422,256.** `docs/map.md` and `docs/overview.html` carried a scored-row count that no
  longer matched the table. The correct value is 70,376 × 6.
- **"Roughly one in ten disease-causing variants are nonsense."** Not supported for ClinVar, where
  the counted figure is 23.7% of pathogenic and likely-pathogenic records. Use the counted figure.
- **"85–97% overall" for footprint alignment.** True of the two HEK293T datasets only; the
  fibroblast libraries run 73–86%. The range must be scoped per dataset.
- **Held-out R² values.** The artifact medians are SRI 0.772 and SJ6986 0.707, not 0.766 and 0.712.

## Citations

| Cited for | Source | Verified |
|---|---|---|
| Readthrough labels, ~5,800 PTCs × 8 drugs, RTDetective | Toledano, Supek & Lehner, *Nat. Genet.* 2024, doi:10.1038/s41588-024-01878-5 | Confirmed: ~5,800 pathogenic PTCs, eight drugs. The paper reports >140,000 individual measurements in total; the ~46,000 figure used here is variant–drug pairs (Extended Data Table 2), which is the level the model is fitted at |
| Published per-drug r² | same | 0.89 / 0.87 / 0.76 / 0.76 / 0.71 / 0.55, against this reproduction's 0.891 / 0.873 / 0.772 / 0.749 / 0.707 / 0.575 under a stricter split |
| SpCas9-NG recognises `NG` | Nishimasu *et al.*, *Science* 2018, doi:10.1126/science.aas9129 | Confirmed |
| SpG recognises `NGN` | Walton *et al.*, *Science* 2020, 368:290–296, doi:10.1126/science.aba8853 | Confirmed |
| SpRY is near-PAMless, `NRN` > `NYN` | same | Confirmed. The code implements `NRN` only; the `NYN` preference tier is not modelled, and the sensitivity arm is unrun in any case |
| BE4max window, protospacer 4–8 | base-editor literature | Confirmed as the standard BE4max window |
| ABE7.10 window, protospacer 4–7 | base-editor literature | Confirmed as the primary window. Sources differ on whether position 8 belongs to it; ADR-0023 takes the narrower 4–7, which is the conservative choice for a reachability count |
| Ingolia 3′ linker `CTGTAGGCACCATCAAT` | Derived from the data, not taken from the archive | The archive's declared linker was found in **zero** reads; this one in 99.4%, with a footprint length distribution peaking at 31 nt |
| Reporter context length | Resolved from the oligo design table | The preprint says 150 nt (75/75); the published Methods say 72 + 3 + 72 = 147. The design table decides |
| `nf-core/riboseq` v1.2.0 | Pinned in `pipeline/nextflow.config` | Consumed as an external pinned pipeline; no nf-core modules are vendored |
| predNMD | No public release | Requested from the authors, not obtained. Recorded in `results/provenance.json` under `unavailable` |
