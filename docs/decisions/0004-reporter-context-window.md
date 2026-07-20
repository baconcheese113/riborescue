# ADR-0004 — The reporter context window

**Status:** accepted · **Date:** 2026-07-19 · **Deciders:** Joseph, Mahan

## Context

The published sources disagree about the reporter context the readthrough library was measured in.
The preprint describes 150 nt as 75 upstream and 75 downstream of the premature stop; the journal
Methods describe 147 nt as 72 + 3 + 72. An off-by-three shifts every positional feature, and the
error is silent — features still compute, models still fit, and the numbers are simply wrong. The
contracts therefore refused to construct a feature window at all until the question was settled.

## Decision

**Settle it from the data rather than from either description.** In `treated_samples.rds` the
upstream triplet, stop codon and downstream triplet align at exactly one offset, for all 48,387
measured rows across every treatment:

| Quantity | Value |
|---|---|
| Upstream context | **72 nt**, for every variant |
| Downstream context | **72 nt or 75 nt**, depending on oligo design |
| Total | 147 nt or 150 nt |

Neither published description is right. The library holds two designs that differ **only
downstream**; upstream context is a constant, and 150 nt is 72 + 3 + 75 rather than 75 + 75.

`REPORTER_UPSTREAM_NT` and `REPORTER_DOWNSTREAM_NT` record this, and `WindowSpec` now accepts a
window exactly when it fits inside *every* variant's context — at most 72 nt either side. A window
reaching further is refused, because one that fits only the longer design silently changes which
variants a feature can be computed for.

`scripts/run_oracle.R` re-derives the offset on every run and writes it to `provenance.json`; a
parity test asserts the constants still match. The measurement is reproducible rather than a claim
made once.

## Consequences

- Feature work is unblocked, with a window whose limits are measured.
- Downstream features are capped at 72 nt even though 26% of variants carry 75, which is the price
  of a feature set defined identically for every variant. A feature that genuinely needs the extra
  three nucleotides must declare the shorter design as missing with a reason code, never as zero.
- The asymmetry is a property of the library, so any external validation set must be checked for its
  own context length rather than assumed to match.
- `CONTRACTS_VERSION` moves to 0.3.0.
