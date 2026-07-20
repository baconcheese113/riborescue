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
| `pipeline/` | The Nextflow pipeline that consumes `nf-core/riboseq` output and produces the scored table |
| `frontend/` | The static web app (Next.js) that presents the table |
| `scripts/` | The R reproduction oracle and data-fetch scripts |
| `data/` | Working inputs and intermediates — fetched from public sources, never committed |
| `docs/decisions/` | Architecture decision records — the reasoning behind each major technical choice |

## Getting started

The toolchain lives in a single Pixi environment on `linux-64` (WSL2 or the dev container), so the
host machine stays clean.

```bash
pixi install        # resolve the pinned toolchain
pixi run check      # lint, type-check, and run the full test suite
pixi run oracle     # fetch the published data and regenerate the reproduction fixtures
pixi run triage --help
```

The Java toolchain the Nextflow editor tooling needs lives in the Pixi environment, so launch the
editor from an activated shell to inherit `JAVA_HOME`:

```bash
pixi shell
code .
```

## How it fits together

- **Development** happens through Pixi — one lockfile pins Python, Nextflow, R, Node, and every
  dependency identically across machines.
- **The pipeline** runs offline on a workstation. It is not hosted; anyone with Docker and Nextflow
  reproduces its outputs from pinned inputs.
- **Results** are published to Zenodo with a DOI. The static app reads them and deploys to the edge
  with no backend.

The reasoning behind each major technical choice is recorded in
[`docs/decisions/`](docs/decisions/).

---

A BIFS 619 group project.
