# ADR-0012 — The sequencing scope is frozen

**Status:** accepted · **Date:** 2026-07-23 · **Deciders:** Joseph

## Context

The project has enough sequencing to make its point. Four series are aligned: the HEK293T and Calu-6
arms of PRJNA576648, the confirmatory GSE144140, and the GSE179274 fibroblast arm. Each further
dataset costs staging, alignment, calibration and a place in the pre-registration, and none of them
changes what the readthrough workflow can demonstrate. ADR-0010 named a mouse liver arm as future
work; this record closes it rather than leaving it open to pull the project along.

## Decision

**No mouse arm.** The GSE179274 mouse liver libraries are not staged, aligned, or analysed. They
need their own GRCm39 reference, index and extension windows, which is a second pipeline for an arm
that supports and cannot confirm.

**No newly acquired sequencing datasets.** The four aligned series are the whole sequencing corpus.
A dataset is not added because it would be interesting; the corpus is closed.

**GSE144140 is the final required sequencing validation.** It is the one dataset that can confirm the
readthrough control within its laboratory and protocol family, and it is the last piece of sequencing
the project needs to run.

**The fibroblast arm may be processed, but never blocks.** The GSE179274 fibroblast libraries are
already aligned, so running them through the frozen workflow is cheap and permitted. They are
independent supporting evidence: they cannot block project completion, and a result from them cannot
trigger a change to the method. If they are never processed, the project is still complete.

## Consequences

The mouse-arm contrast fixed in ADR-0010 is withdrawn from the plan of work. The other contrasts of
that record stand.

Sequencing is no longer on the critical path. What remains is the native-stop safety atlas, the
frontend and a reproducible release — neither of which needs another alignment run.
