# The RiboRescue pipeline

Four workflows, selected with `--step`. `amenability` places ClinVar's pathogenic nonsense variants
on their MANE Select transcripts, scores them under each therapy and ranks suppressor designs;
`reads` takes the sequencing runs through quality control, adapter trimming and alignment;
`validate-labels` checks measured readthrough labels; `triage` classifies submitted variants. Every
process calls the installed `riborescue` command, so all scientific logic lives in Python and is
tested there.

`validate-labels` and `triage` require an upstream handoff. `amenability` and `reads` do not.

```bash
nextflow run pipeline --step amenability -profile local \
    --clinvar data/clinvar/clinvar_20260715.vcf.gz \
    --mane_annotation data/mane/MANE.GRCh38.v1.5.refseq_genomic.gff.gz \
    --mane_transcripts data/mane/MANE.GRCh38.v1.5.refseq_rna.fna.gz \
    --mane_proteins data/mane/MANE.GRCh38.v1.5.refseq_protein.faa.gz \
    --training 'tests/fixtures/oracle/features_*.tsv.gz' \
    --held_out 'tests/fixtures/oracle/predictions_*.tsv.gz' \
    --outdir results/pipeline
```

## Reads

`pipeline/assets/riboseq_samples.tsv` declares the sequencing runs with the checksums ENA publishes
for them. Staging fetches each one and writes a sheet naming it on disk; nothing is committed.

```bash
pixi run stage-runs      # fetch the FASTQ, verifying every digest
pixi run reads           # FastQC, cutadapt, FastQC again, MultiQC
```

The adapter is declared per run rather than detected, because the archive's own record of it is not
reliable, and `riborescue trim-summary` refuses a footprint library whose declared adapter was found
in under half its reads. A transcriptome library is not held to that floor: most of its fragments
are longer than the read, so the adapter is never reached.

Naming a reference adds alignment; leaving it out stops after trimming, which keeps the quick path
quick:

```bash
nextflow run pipeline --step reads -profile local \
    --samplesheet results/staged_runs.tsv \
    --genome data/gencode/GRCh38.primary_assembly.genome.fa.gz \
    --annotation data/gencode/gencode.v50.primary_assembly.annotation.gtf.gz \
    --outdir results/reads
```

The index samples every second suffix array entry, because a dense human index needs more memory
than a 30 GB machine has. It is kept in `data/star` rather than the work directory, so the hours it
takes are paid once however often the pipeline reruns.

## nf-core

The committed results use the `reads` workflow. `riborescue validate-handoff` can also check the
outputs of a pinned `nf-core/riboseq` run. `pipeline/tests/data/handoff.json` is an example.

## Running

```bash
nextflow run pipeline --step triage -profile test,local          # the committed fixtures
nextflow run pipeline --step triage -profile docker \
    --handoff handoff.json --results_root results/riboseq --variants variants.tsv
```

`-profile docker` runs every process in the `riborescue` image, built from `pixi.lock` so the
environment inside the container is the one the tests ran against; `-profile local` uses the
`riborescue` on `PATH`. Both produce byte-identical artifacts.

| Parameter | What it names |
|---|---|
| `--step` | `amenability`, `reads`, `validate-labels` or `triage` |
| `--samplesheet` | The staged runs to take through quality control and trimming (`reads`) |
| `--genome`, `--annotation` | GENCODE primary assembly and GTF; naming both adds alignment (`reads`) |
| `--clinvar` | The ClinVar VCF to draw the variant population from |
| `--mane_annotation`, `--mane_transcripts`, `--mane_proteins` | MANE Select annotation, sequences and reference proteins |
| `--training`, `--held_out` | The oracle's per-therapy feature tables and held-out rounds |
| `--handoff` | The upstream handoff manifest |
| `--results_root` | The `nf-core/riboseq` results tree the manifest describes |
| `--labels` | Measured readthrough efficiency per variant × therapy (`validate-labels`) |
| `--variants` | Variants to classify (`triage`) |
| `--outdir` | Where published results land |

## Tests

```bash
pixi run nf-test    # end-to-end over the fixtures in tests/data
pixi run nf-lint
```
