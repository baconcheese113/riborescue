# Results

This directory tracks only the small outputs quoted elsewhere or needed to support the
main findings. Large and intermediate files stay out of Git.

## BIFS 619 RNA-seq analysis

The analysis uses four paired-end Calu-6 RNA-seq samples from PRJNA576648: two untreated and two
treated with G418. The [sample manifest](../pipeline/assets/riboseq_samples.tsv) lists each run.

| Requirement | Evidence |
|---|---|
| Quality control | [Report](reads/qc/multiqc_report.html), [figure](reads/qc/qc_overview.png), [table](reads/qc/course_qc_summary.tsv) |
| Read cleaning | [Figure](reads/qc/pre_post_cleaning.png), [table](reads/qc/trim_summary.tsv) |
| STAR alignment | [Figure](reads/qc/alignment_metrics.png), [table](reads/qc/alignment_summary.tsv) |
| Gene expression | [TPM](rnaseq/tpm.tsv), [top 20](rnaseq/top_expressed.tsv), [heatmap](rnaseq/top_expressed.png) |

Regenerate the figures with `pixi run course-figures` and `pixi run expression-figures`.

Cleaning retains more than 99% of read pairs. G418 samples map at about 68%; untreated samples map
at about 20%. Mitochondrial and small RNA differ between samples, so the
[nuclear-gene table](rnaseq/top_expressed_nuclear.tsv) and
[heatmap](rnaseq/top_expressed_nuclear.png) are also shown.

## RiboRescue conclusions

| Question | Result |
|---|---|
| Model reproduction | [Evaluation](evaluation.tsv) |
| Variant coverage | [Summary](landscape_summary.tsv) |
| G418 readthrough | [G418 result](readthrough/gse144140/g418_vs_dmso.paired.published.tsv) |
| Negative control | [SRI-37240 result](readthrough/gse144140/sri37240_vs_dmso.paired.published.tsv) |
| Codon occupancy test | [Permutation result](kinetics/permutation_null/familywise.tsv) |
| Follow-up work | [Experiment designs](experiments.tsv) |

G418 increases reads past stop codons; SRI-37240 does not. The codon occupancy test finds no clear
improvement over sequence alone.
