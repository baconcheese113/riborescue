# The RiboRescue pipeline

Three workflows, selected with `--step`. `amenability` places ClinVar's pathogenic nonsense variants
on their MANE Select transcripts, scores them under each therapy and ranks suppressor designs;
`train` fits the readthrough model from measured labels; `score` triages submitted variants. Every
process calls the installed `riborescue` command, so all scientific logic lives in Python and is
tested there.

Only the workflows consuming Ribo-seq output need the upstream handoff. The amenability path runs
from public annotation alone:

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

## Upstream

Ribo-seq processing is `nf-core/riboseq`, run separately at a pinned revision and never vendored:

```bash
nextflow run nf-core/riboseq -r 1.2.0 -profile docker \
    --input samplesheet.csv --outdir results/riboseq
```

RiboRescue then reads a handoff manifest naming the outputs it consumes — P-site offsets, codon and
CDS coverage, RNA-seq counts and TPM, alignments, MultiQC — and reads nothing else from that tree.
`pipeline/tests/data/handoff.json` is the worked example; `riborescue validate-handoff` checks it,
and the run stops before any analysis if a declared output is absent.

## Running

```bash
nextflow run pipeline --step score -profile test,local          # the committed fixtures
nextflow run pipeline --step score -profile docker \
    --handoff handoff.json --results_root results/riboseq --variants variants.tsv
```

`-profile docker` runs every process in the `riborescue` image, built from `pixi.lock` so the
environment inside the container is the one the tests ran against; `-profile local` uses the
`riborescue` on `PATH`. Both produce byte-identical artifacts.

| Parameter | What it names |
|---|---|
| `--step` | `amenability`, `train` or `score` |
| `--clinvar` | The ClinVar VCF to draw the variant population from |
| `--mane_annotation`, `--mane_transcripts`, `--mane_proteins` | MANE Select annotation, sequences and reference proteins |
| `--training`, `--held_out` | The oracle's per-therapy feature tables and held-out rounds |
| `--handoff` | The upstream handoff manifest |
| `--results_root` | The `nf-core/riboseq` results tree the manifest describes |
| `--labels` | Measured readthrough efficiency per variant × therapy (`train`) |
| `--variants` | Variants to triage and score (`score`) |
| `--outdir` | Where published results land |

## Tests

```bash
pixi run nf-test    # end-to-end over the fixtures in tests/data
pixi run nf-lint
```
