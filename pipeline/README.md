# The RiboRescue pipeline

Two workflows over one validated upstream handoff. `--step train` fits the readthrough model from
measured labels; `--step score` applies it to submitted variants. Every process calls the installed
`riborescue` command, so all scientific logic lives in Python and is tested there.

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

`-profile docker` runs every process in the pinned `riborescue` image; `-profile local` uses the
`riborescue` on `PATH`, which is how the Pixi environment and the tests run it.

| Parameter | What it names |
|---|---|
| `--step` | `train` or `score` |
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
