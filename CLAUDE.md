# RiboRescue — working conventions

RiboRescue matches a patient's nonsense variant to candidate readthrough therapies. It is a
Nextflow pipeline that produces a scored variant × therapy table, a Python package that builds and
serves that table, and a static web app that presents it. It is also a BIFS 619 course project.

## Writing style — read this first

**Everything in the repo reads as if one author wrote it all at once.** Code, comments, docs, and
commit-adjacent text are minimal, present-tense, and timeless.

- No diary or changelog voice. Never write "revert", "now using X", "previously", "post-approval",
  "as discussed", "see ADR-000x", "TODO from review", or anything describing how the code got here.
- A comment explains *what the code does or why it must be this way* — never its history.
- Prefer no comment to an obvious one. Delete a comment before letting it drift stale.
- Decision history has exactly one home: `docs/decisions/` (ADRs). Keep it out of code entirely.
- What broke and how it resolved has exactly one home: `docs/notebook.md`. Keep it out of code too.
- Don't leave stub files for deferred work. A directory or module appears when it holds real code.

## Notebook

`docs/notebook.md` is the running record of problems. Append an entry under **Challenges and
resolutions** whenever one resolves — a wrong result traced to its cause, a check that fired, a
benchmark that misled, an approach abandoned, an external claim that failed verification. Each
entry states the challenge, what went wrong, the resolution, what it cost, and the lesson.

Write it while the problem is fresh; a resolved problem nobody wrote down is gone. It is
chronological and informal by design — the single exception to the writing style above — and it is
gitignored, because it is working material rather than part of the repo.

## Layout

```
pyproject.toml pixi.lock the riborescue package, the Pixi toolchain, and one lockfile
pyproject.toml           the riborescue Python package (src-layout) + tool config
src/riborescue/          scientific logic: contracts, triage, evaluation, cli
tests/                   pytest — unit, property (hypothesis), controls (protected)
pipeline/                hand-authored minimal Nextflow; consumes nf-core/riboseq externally
frontend/                Next.js static export, served from Cloudflare Pages
scripts/                 R oracle and data-fetch scripts
data/                    gitignored; inputs are fetched, never committed
docs/decisions/          ADRs — the only place decision history lives
docs/notebook.md         gitignored; problems, resolutions, and report material
```

## Development

The host stays clean; everything runs through Pixi on `linux-64` inside WSL2.

```
pixi install             resolve the toolchain from pixi.lock
pixi run check           lint + types + tests + controls (the full gate)
pixi run test            pytest
pixi run lint            ruff
pixi run types           pyright
pixi run app-dev         Next.js dev server
```

Pixi runs tasks through a portable shell, not bash: loops and other reserved words are rejected. A
task that needs them calls a script in `scripts/` instead of inlining the shell.

Environments are `default` to develop in, `runtime` for what the container ships, and `psite` for
riboWaltz, which pins an older R than the reproduction oracle. Run a task elsewhere with
`pixi run -e <env> <task>`.

## Code conventions

- Python ≥ 3.12, native `X | Y` unions, type hints on every signature, checked with **pyright**.
- **ruff** for lint and format. Line length 100.
- **pytest** with `--import-mode=importlib` and `--strict-markers`; no `__init__.py` under `tests/`.
  Markers: `control`, `parity`, `integration`, `slow`.
- **`tables.py` is the validation boundary**: every table crossing a pipeline step is read and
  written through it, checked against **pandera** schemas that refuse bad state rather than coercing
  it. `contracts.py` is the vocabulary those schemas validate *to* — the closed sets a column may
  hold, and the reporter geometry every feature window is measured against. **Pydantic** models
  where a record is passed around rather than a table; plain dataclasses and **Click** elsewhere.
- Scientific logic is Python called from Nextflow; no analysis logic in Groovy.
- Prefer a well-maintained library over new code when it genuinely reduces complexity.

## Scaffolding

Generate files with the tool that owns them only when the output is lean; hand-write the rest.
`pixi init` / `pixi add` own the toolchain manifest. The Nextflow pipeline is **hand-authored and
minimal** — `main.nf` with the `HANDOFF`, `READS`, `LABELS`, `TRIAGE` and `AMENABILITY` entry
workflows, `nextflow.config`, a small `conf/`, and nf-test — because it consumes `nf-core/riboseq`
as an external pinned pipeline and authors almost no modules of its own. `nf-core pipelines create` is not used: it emits dozens of official-pipeline
files we do not need. Tooling runs natively in the WSL2 environment on `linux-64`.

## Provenance

Reproducibility rests on regeneration, not stored bytes: `pixi.lock` + pinned `nf-core/riboseq` +
container digests + fetch scripts with checksums + committed oracle fixtures. Large artifacts go to
Zenodo (DOI) and Cloudflare R2, never into git. See `docs/decisions/`.
