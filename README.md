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

## Two arms, and which one you can run

Everything here is either **predicted** or **measured**, and the project never blends the two. Which
arm you can run depends entirely on the machine you have.

| | **Predicted** | **Measured** |
|---|---|---|
| What it is | A model scoring 70k ClinVar variants against six compounds | Ribosome profiling of real cells, reading where ribosomes actually sit |
| What you need | A laptop. **16 GB RAM is enough** | A workstation: STAR builds a ~32 GB index |
| Time | Minutes | Hours, plus a large download |
| What you get | The scored table and the whole web app | The readthrough control, the native-stop safety atlas, the evidence page |

**If you have a normal computer, run the predicted arm.** It is the majority of the project and it
needs no sequencing data at all. Start there.

## Getting started

The toolchain lives in Pixi on `linux-64` (Ubuntu, WSL2, or the dev container), so your machine
stays clean — Pixi is the only thing you install, and it brings Python, Node, Nextflow and R itself.

```bash
curl -fsSL https://pixi.sh/install.sh | bash   # then restart the shell
pixi install                                   # resolve the pinned toolchain
pixi run check                                 # lint, type-check, run the tests
```

### The predicted arm, on a laptop

```bash
pixi run fetch clinvar_grch38 mane_annotation mane_transcripts mane_proteins
pixi run clinvar && pixi run contexts       # ClinVar nonsense variants, placed on MANE transcripts
pixi run score && pixi run landscape        # the scored variant × therapy table
pixi run diseases                           # ClinVar conditions → MedGen/OMIM/Orphanet
pixi run site                               # build the viewer payloads and the static app
pixi run app-dev                            # http://localhost:3000
```

`site` writes `riborescue.json` (the per-variant example the explorer and patient views read) and
`riborescue_research.json` (the coverage aggregate behind the researcher dashboard). The app has
three views over them: **/patient** lays each therapy out as a card with its evidence slots,
**/researcher** draws the coverage frontiers with their denominators, and **/explorer** is the full
per-variant table.

To refit the published model and check it still matches the R oracle:

```bash
pixi run oracle && pixi run test-slow
```

### The measured arm, on a workstation

Sequencing is the half that needs real hardware. Reads are fetched and verified, aligned, then
calibrated with riboWaltz before any contrast runs:

```bash
pixi run stage-runs         # fetch the FASTQ named in pipeline/assets/, checksum-verified
pixi run reads              # FastQC → cutadapt → deplete → STAR → metrics
pixi run -e psite psite     # riboWaltz P-site offsets and periodicity
```

With those in place, the measured results and the panels that depend on them:

```bash
pixi run readthrough-gse144140   # does a drug really push ribosomes past a stop?
pixi run atlas && pixi run atlas-predict
pixi run export-web-safety       # adds the native-stop safety panel to the app
pixi run export-evidence         # the evidence payload: controls, calibration, codon signature
```

These need the sequencing results to exist and stop rather than invent them, which is why the
predicted arm never calls them.

### Running it as a pipeline

A single `riborescue` command backs every step — `pixi run riborescue --help` lists them all. The
Nextflow entry workflows live in [`pipeline/`](pipeline/README.md), which documents how to run
variant scoring and read processing as pipeline steps rather than one task at a time.

## For the course

The BIFS 619 requirements — quality control, read cleaning, alignment, and gene expression over the
Calu-6 RNA-seq samples — are documented with their outputs in
[`results/README.md`](results/README.md), which names the figure and table behind each requirement
and the command that regenerates them.

## How it fits together

Development happens through Pixi, so one lockfile pins every tool identically across machines. The
pipeline runs offline and is not hosted — anyone with Docker and Nextflow reproduces its outputs
from pinned inputs. Results publish to Zenodo with a DOI, and the static app reads a compact JSON
export with no backend.

| Path | What lives here |
|---|---|
| `src/riborescue/` | The Python package: data contracts, variant triage, the evaluation harness, and the CLI |
| `tests/` | `pytest` — unit, property-based (`hypothesis`), and protected negative controls |
| `pipeline/` | The hand-authored Nextflow pipeline |
| `scripts/` | R analysis (reproduction oracle, P-site calibration) and data-fetch scripts |
| `frontend/` | The static web app (Next.js) that presents the table |
| `data/` | Fetched inputs — public sources, verified by checksum, never committed |
| [`results/`](results/) | Small course results and project conclusions |
| [`docs/decisions/`](docs/decisions/) | Decision records — why each major choice was made, indexed with a reading order |
| `Dockerfile` | The runtime image Nextflow processes use |
| `.github/workflows/` | The gate on every push: lint, types, tests, parity, pipeline, image |

---

A BIFS 619 group project.
