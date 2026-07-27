# ADR-0003 — Pipeline orchestration and the upstream handoff

**Status:** accepted, amended in part by [ADR-0024](0024-the-upstream-boundary-and-what-it-carries.md)
· **Date:** 2026-07-19 · **Deciders:** Joseph, Mahan

> The workflow inventory below is the one this record fixed and is no longer what the pipeline
> holds: `READS` was added for local preprocessing, `TRAIN` and `SCORE` are now `LABELS` and
> `TRIAGE`, and the handoff carries four outputs rather than seven. ADR-0024 states where the
> upstream boundary falls today. Everything else here — the external upstream, processes calling
> the CLI and nothing else, the validated manifest, nf-test over committed fixtures — still governs.

## Context

Ribo-seq processing is a solved problem that `nf-core/riboseq` already does well, and RiboRescue's
own orchestration is small: validate inputs, fit a model, score variants. The risks are the ones
that make scientific pipelines quietly wrong — analysis logic drifting into Groovy where it cannot
be tested, downstream steps reaching into arbitrary paths of an upstream results tree, and an
upstream version that moves under the results.

`nf-core pipelines create` was evaluated and rejected: it emits roughly sixty files of
official-pipeline machinery (governance, badges, module boilerplate, a schema builder) against three
processes of real work.

## Decision

**A hand-authored Nextflow pipeline in `pipeline/`, minimal and fully linted.**

- `nf-core/riboseq` stays an external dependency, run separately at a pinned revision. Nothing from
  it is vendored.
- **Two workflows, `TRAIN` and `SCORE`**, selected with `--step`. Nextflow's strict parser does not
  support `-entry`, so the entry workflow dispatches on the parameter.
- **Processes call the `riborescue` CLI and nothing else.** Groovy composes; Python decides. Every
  scientific step is therefore unit-tested in `tests/` rather than only end-to-end.
- **The upstream handoff is a validated manifest** (`riborescue.handoff.UpstreamHandoff`): it names
  each consumed output as a path relative to a results root, refuses paths that escape that root,
  refuses a moving revision such as `main` or `dev`, and refuses any pipeline other than
  `nf-core/riboseq`. A run stops before analysis if a declared output is absent.
- **Two profiles.** `docker` runs each process in the pinned `riborescue` image; `local` uses the
  `riborescue` on `PATH`, which is how the Pixi environment and CI run the tests.
- **`nf-test` over committed fixtures** exercises both workflows and the refusal path, with
  `pipefail` set so a failing command that feeds a pipe cannot report success.

## Consequences

- The pipeline is small enough to read in one sitting, and its correctness is mostly Python
  correctness.
- Upstream and downstream evolve independently; upgrading `nf-core/riboseq` is a revision bump in
  the manifest, and the handoff validation says immediately whether the outputs still line up.
- Handoff validation gives one place to answer "which upstream files does this project actually
  consume?" — the schema is the answer, and no process can quietly widen it.
- The pipeline runs offline on a workstation, so there is no scheduler or cloud executor config to
  maintain; a cluster profile can be added later without touching the workflows.
