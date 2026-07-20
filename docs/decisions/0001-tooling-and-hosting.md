# ADR-0001 — Tooling, dev environment, and hosting

**Status:** accepted · **Date:** 2026-07-19 · **Deciders:** Joseph, Mahan

## Context

RiboRescue has two execution contexts that shape every tooling choice:

- **An offline pipeline** — Nextflow orchestrating STAR (~32 GB index), Ribo-seq QC, and
  ESM-1v on the GPU. It runs on the workstation to *produce* the scored table. It is never
  hosted. Its dependencies include non-PyPI tools (`nextflow`, `nf-test`, `R`).
- **A static app** — serves the *precomputed* table. Hosted free at the edge, no backend,
  zero compute at serve time.

A sibling project (GhostFrame, BIFS 617) uses `uv` + Makefile + a single Dockerfile. That is
correct *for GhostFrame*, which is pure-Python + REST with an all-PyPI dependency set. It does
not transfer here: `uv` is PyPI-only and cannot install `nextflow`, `nf-test`, or `R`.

## Decision

**Dependency + environment manager: Pixi.** One `pixi.lock` locks conda-forge/bioconda tools
(`nextflow`, `nf-test`, `nf-core`, `r-base`, `nodejs`, `python`) *and* every PyPI dependency
together, reproducibly, on both team machines. This is a single-lockfile provenance story over a
mixed conda+PyPI+JVM+R toolchain — precisely what `uv` + a Dockerfile would split into two
systems. Env is project-local, so the host stays clean.

**Heavy bio tools stay out of the dev env.** Nextflow pulls STAR/samtools/etc. as per-process
biocontainers at run time; the dev env carries only orchestration + analysis tooling.

**Quality tooling:** `pyright` (types), `ruff` (lint+format), `pytest`
(`--import-mode=importlib`, `--strict-markers`), `hypothesis` (property tests), `pre-commit`.

**Contracts use Pydantic v2 + pandera.** RiboRescue has no API, but the data contracts *are* the
validation boundary (refuse-to-construct `WindowSpec`, missing≠zero, separate rescue factors), so
Pydantic/pandera earn their place there. Plain dataclasses + Click elsewhere.

**Frontend: Next.js with static export (`output: 'export'`).** Chosen for **discoverability**:
`generateStaticParams` pre-renders one static, crawlable page per gene / common variant — the
largest organic-search lever — and Next's metadata/sitemap/JSON-LD tooling is the most mature.

**App hosting: Cloudflare Pages.** The only free tier with genuinely unlimited bandwidth and
requests. Large scored-table Parquet lives on **Cloudflare R2** (free, zero egress) or the
**Zenodo DOI**, fetched by the app at runtime, sidestepping the 25 MB/file Pages limit.

## Consequences

- Mental model: **pixi = how we develop · containers = how the pipeline reproducibly runs ·
  Cloudflare/Zenodo = how results are served.** Each layer, one job.
- `pixi` runs on `linux-64` (WSL2 or the dev container); bioconda has no `win-64` builds. Both
  members develop on WSL2, keeping Windows hosts clean.
- Provenance is layered: `pixi.lock` (toolchain) + `CONTRACTS_VERSION` + ADRs + committed oracle
  fixtures + controls manifest + model bundles with digests + pinned `nf-core/riboseq` version +
  container digests + Zenodo DOI.
- Not adopting `uv` here is deliberate, not an oversight; see Context.
