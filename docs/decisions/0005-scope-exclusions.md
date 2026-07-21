# ADR-0005 — Scope exclusions

**Status:** accepted · **Date:** 2026-07-21 · **Deciders:** Joseph, Mahan

## Context

The scored table invites summaries it cannot support. Once variants carry disease names, prevalence
figures and modality labels, several natural-looking aggregations become available — which disease
has the largest rescue gap, which indication is most neglected, how many patients a suppressor
design would reach. Each is a ranking whose axis measures curation effort rather than biology.

ClinVar is a submission archive: its variant counts track which genes sit on commercial panels.
Orphanet is a curation archive: what it omits records what has been curated, not what exists.
Searching either for evidence of a therapy and finding none establishes that the search found
nothing. A figure built on that quantity is a figure of its own methods.

The pull toward these summaries is strong because they are the ones an outside reader wants, so the
boundary is written down rather than left to judgement at figure-drafting time.

## Decision

**Eight exclusions, each with the supported statement that replaces it.**

| Excluded | Why | What ships instead |
|---|---|---|
| Therapeutic-attention or unmet-need rankings | The quantity is absence of evidence in searched sources, bounded by which sources are accessible | Treatment status as a categorical annotation, separately sourced |
| Any figure axis built on absence of evidence | The axis measures search effort | Coverage, ambiguity and missingness reported as results in their own right |
| Patient counts from prevalence multipliers | A point estimate manufactured from a range and a population figure | Prevalence in its native range or class, with geography and measure type attached |
| Disease rankings framed as commercial or clinical targets | §6.5 conflict-of-interest clearance applies, and the framing outruns the evidence | Research opportunity, with every contributing factor visible and sourced |
| Blending theoretical and observed variant universes | One is ascertainment-independent, the other is not; a combined number is neither | Two universes reported side by side, never summed or averaged |
| Novelty claims for selection against leaky native stops | Stop-codon context bias is long established | Whether *this* model recovers it unprompted, which is external validation |
| Unregularised saturated interaction models | Capacity that random folds reward and grouped splits expose | ADR-0006 |
| Public processed Ribo-seq treated as equivalent to raw | P-site calibration, QC thresholds and provenance are not consistent across sources | Processed data for exploration, labelled as such; raw data where a claim depends on it |

Two statements already load-bearing elsewhere are restated here because they are the ones most
easily lost in aggregation: **ClinVar gives variant coverage, never patient coverage**, and
**absence from ClinVar or Orphanet is not biological absence**.

Treatment status is categorical and its states never collapse into each other:

```
approved disease-modifying treatment
approved symptomatic management
investigational therapy
no therapy found in searched sources
treatment status unknown
```

`no therapy found in searched sources` is not `no therapy exists`, and the schema keeps them
distinct so that no downstream aggregation can merge them.

## Consequences

- The disease layer ships as annotation and enumeration. It produces no ordering of diseases.
- Every prevalence figure carries its range, geography and measure type, so a reader can see the
  width of what is being claimed.
- Join coverage and ambiguous condition mappings are headline numbers, not footnotes. A layer that
  maps 60% of variants to a disease reports that fraction wherever it reports anything else.
- The native-stop analysis is framed as validation, and its figure caption says so.
- Where an excluded quantity is genuinely wanted, the route is a new ADR, not an exception.
