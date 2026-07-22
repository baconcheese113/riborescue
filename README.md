# RiboRescue

**Which readthrough therapy, for this nonsense variant, in this tissue?**

Nonsense mutations create a premature stop codon that truncates a protein and causes roughly 11% of
inherited genetic disease. Readthrough therapies — small molecules and engineered suppressor tRNAs —
can make the ribosome ignore that premature stop, but efficacy varies enormously between mutations
and between therapies, and no public tool matches a patient's variant to a candidate therapy.

RiboRescue integrates readthrough efficiency, transcript availability (nonsense-mediated decay),
inserted-residue compatibility, and native-stop safety into one scored **variant × therapy** table,
with uncertainty shown by default.

> **Research use only.** RiboRescue makes no clinical or diagnostic claims and produces no
> patient-facing recommendations.

## Repository layout

| Path | What lives here |
|---|---|
| `src/riborescue/` | The Python package: data contracts, variant triage, the evaluation harness, and the CLI |
| `tests/` | `pytest` suites — unit, property-based (`hypothesis`), and protected negative controls |
| `pipeline/` | The hand-authored Nextflow pipeline: Ribo-seq read processing and the scored variant × therapy table |
| `scripts/` | R analysis (reproduction oracle, P-site calibration) and data-fetch scripts |
| `frontend/` | The static web app (Next.js) that presents the table — *planned* |
| `data/` | Fetched inputs — from public sources, verified by checksum, never committed |
| `results/` | Pipeline and analysis outputs — regenerated from inputs, never committed |
| `docs/decisions/` | Architecture decision records — the reasoning behind each major technical choice |
| `Dockerfile` | The runtime image Nextflow processes use, carrying the `riborescue` command |
| `.github/workflows/` | The gate that runs on every push: lint, types, tests, reproduction parity, pipeline, image |

## Getting started

The toolchain lives in Pixi on `linux-64` (WSL2 or the dev container), so the host machine stays
clean. One lockfile pins Python, Nextflow, R and every tool; `pixi install` resolves it. It holds
three environments: `default` to develop in, `runtime` for what the container ships, and `psite`
for riboWaltz, which pins an older R than the reproduction oracle.

```bash
pixi install        # resolve the pinned toolchain
pixi run check      # lint, type-check, and run the full test suite
pixi run fetch      # fetch the public inputs (ClinVar, MANE, GENCODE), verified by checksum
```

**Reproduce the published model** — refit it and check parity to the R oracle:

```bash
pixi run oracle
pixi run test-slow
```

**Score the variants** — ClinVar's nonsense variants placed on transcripts and scored per therapy:

```bash
cd pipeline && nextflow run . --step amenability -profile local \
    --clinvar ../data/clinvar/*.vcf.gz \
    --mane_annotation ../data/mane/*.gff.gz \
    --mane_transcripts ../data/mane/*.rna.fna.gz \
    --mane_proteins ../data/mane/*.protein.faa.gz \
    --training '../tests/fixtures/oracle/features_*.tsv.gz' \
    --held_out '../tests/fixtures/oracle/predictions_*.tsv.gz'
```

**Process the Ribo-seq data** — QC, trimming, rRNA depletion, alignment, then P-site calibration:

```bash
pixi run stage-runs         # fetch the FASTQ named in pipeline/assets/, verified
pixi run reads              # FastQC → cutadapt → deplete → STAR → metrics
pixi run -e psite psite     # riboWaltz P-site offsets and periodicity
```

A single `riborescue` command backs every pipeline step; `pixi run triage --help` lists them.

## How it fits together

- **Development** happens through Pixi — one lockfile pins Python, Nextflow and R identically across
  machines.
- **The pipeline** runs offline on a workstation. It is not hosted; anyone with Docker and Nextflow
  reproduces its outputs from pinned inputs.
- **Results** publish to Zenodo with a DOI; a static app (planned) will read them and deploy to the
  edge with no backend.

The reasoning behind each major technical choice is recorded in
[`docs/decisions/`](docs/decisions/).

---

A BIFS 619 group project.
