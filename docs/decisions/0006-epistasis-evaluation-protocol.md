# ADR-0006 — The upstream–downstream interaction comparison

**Status:** accepted · **Date:** 2026-07-21 · **Deciders:** Joseph, Mahan

## Context

The published model is factorized. Its linear predictor is

```
RT_binomial ~ 0 + stop_type + down_123nt + up_123nt + stop_type:down_123nt
```

so the upstream triplet enters additively and there is no `up_123nt:down_123nt` term. The model
composes an unseen combination of triplets from marginals it has estimated separately, and its
prediction surface — 3 × 64 × 64 contexts per drug — rests on a few hundred design columns.

That is a structural assumption about the biology, not a fitted result: it asserts that upstream and
downstream context act independently. Deep mutational scanning of programmed readthrough
(`10.64898/2026.07.14.738498`) reports determinants that combine nonlinearly, with a limited number
of strong pairwise interactions. Whether upstream–downstream dependence carries held-out signal for
drug-induced readthrough is therefore an open and testable question against a baseline already
reproduced to parity.

The naive test fails. A saturated categorical `up_123nt:down_123nt` term is 4,096 levels against
roughly six thousand variants per drug, so most levels are absent or carry a single observation.
The fit aliases heavily, memorises repeated contexts, and is rewarded for it by random folds —
producing an apparent improvement that is variance. The comparison has to account for the
interaction model's far greater capacity before it means anything.

The published model already leaves support holes at the interaction it does fit: an empty
`stop_type × down_123nt` cell yields a dependent column, which `identifiable_columns` drops and the
prediction treats as zero effect, falling back to the marginal. How many scored variants land in
such a cell is unmeasured.

## Decision

**Four models, one fixed protocol, and a decision rule written before any number is seen.**

| # | Model | Isolates |
|---|---|---|
| 1 | Published four-term GLM | The reproduction baseline, unchanged |
| 2 | Regularised `up × down` interaction | Whether dependence exists at all, with capacity controlled |
| 3 | Sparse or low-rank interaction | Whether it is concentrated in few pairs rather than diffuse |
| 4 | 47-feature ElasticNet | The published comparator, kept distinct from the four-term model |

**Selection is nested and never touches the evaluation folds.** Outer folds are the authors' own
round assignments from the oracle fixtures. Regularisation strength and interaction rank are chosen
by validation *inside* each training fold. Model 1 is refit unchanged so that Toledano parity holds
throughout.

**The evaluation grid is fixed in advance**, and every model reports every row:

- published random 90/10 rounds
- grouped by gene
- grouped by sequence-context cluster
- stratified by interaction-cell support, so gains concentrated in well-observed cells are visible
  separately from gains claimed where the model has little to go on

**Controls.** A permuted-interaction control shuffles the upstream triplet within matched
`(stop_type, down_123nt)` cells, breaking the dependence while preserving both marginals and the
fitted interaction. Its ΔR² must include zero. Per-drug bootstrap intervals throughout, squared
Pearson correlation, normalised against per-drug replicate reliability.

**Decision rule.** Upstream–downstream dependence is reported as an improvement only when the gain's
bootstrap interval excludes zero under **grouped** splits *and* the permuted control's interval
includes zero. A gain appearing only under random splits is reported as capacity, not signal.

## Consequences

- The protocol is fixed before the numbers exist, so a marginal result cannot be rescued by choosing
  a split or a metric after the fact.
- The negative outcome is a result: it would show the published factorization is adequate for
  drug-induced readthrough and that the programmed-readthrough epistasis does not transfer.
- The support atlas is a prerequisite rather than a companion — the support strata and the
  aliased-cell audit are inputs to this comparison, and the aliased-cell count is reported whether or
  not the interaction models are pursued further.
- The extended sequence window remains a separate ablation under ADR-0004's measured limits. Window
  length and interaction structure are not varied together.
