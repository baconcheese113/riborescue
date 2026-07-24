# ADR-0014 — Suppressor-tRNA panel coverage, and its frozen predicate

**Status:** accepted · **Date:** 2026-07-24 · **Deciders:** Joseph

## Context

The per-design ranking (`trna-coverage`) answers "how many variants would this one suppressor tRNA
reach?" It does not answer the design question a wet lab actually faces: given a budget of *k*
engineered tRNAs, which set reaches the most variants, and what does each added tRNA buy? That is a
maximum-coverage problem, and it is the first researcher-dashboard result that needs no data the
repository does not already hold.

Two things have to be fixed before a frontier is computed, or the frontier is not interpretable.

**The coverage predicate.** A suppressor tRNA is a stop codon it decodes and a residue it inserts.
"Covered" could mean the design restores the exact native residue, or that it inserts any residue a
conservative substitution tolerates, or something weaker still. These give different frontiers. The
predicate has to be frozen *before* the frontier is examined, so the answer is not chosen by the
picture it produces.

**The objective.** Variants covered, genes covered, and patients affected are different numbers over
the same panel. Reporting one as "coverage" without saying which invites the inflated single count
this project exists to avoid.

## Decision

**Predicate, v1 — exact restoration.** A design covers a variant when it decodes the variant's stop
codon *and* inserts the residue the protein natively carries there. This is *model coverage*: the
design would put the correct residue back. It is not a claim the protein folds, functions, or that
any tRNA is safe or clinical, and it is labelled that way everywhere it surfaces.

**No safety axis.** The frontier optimises coverage alone. The native-stop atlas measures G418 in one
cell line (ADR-0009 assay, HEK293T) and does not transfer to engineered tRNAs; there is no comparable
suppressor-tRNA safety evidence in the repository, so there is no safety objective to trade against.
Inventing one would blend a measured layer into a design score — the merge this project refuses.

**Objectives kept separate.** `variants` and `genes` are computed as independent frontiers, each with
its own greedy order. Diseases and affected-patient counts are deferred to the disease-mapping work
(ADR to follow), where "a disease is covered" needs its own defined denominator — one variant, a
fraction of a gene's known nonsense variants, or an unmet-need threshold — and must not be produced by
swapping variant ids for disease labels.

**Engine.** Greedy maximum coverage: at each step the design adding the most uncovered elements, ties
broken by design id so the panel is deterministic. It emits the full auditable frontier — the panel
at every size *k*, cumulative and marginal coverage, and the count still uncovered — so any panel is
the first *k* rows. Greedy carries the standard `1 − 1/e` guarantee for overlapping sets, checked
against a brute-force optimum on small fixtures.

## Consequences

Under the frozen predicate the optimisation is degenerate *by construction*, and this is stated
rather than hidden: each variant has one stop and one native residue, so it is covered by exactly one
design. The covered sets are therefore disjoint — they partition the variants — and greedy maximum
coverage is provably optimal with no approximation gap, reducing to sorting designs by size. The
`1 − 1/e` machinery and the brute-force parity test are exercised on a synthetic *overlapping*
fixture, because that is the case a future predicate creates, not the case v1 is in.

The generic engine is the reason the eventual conservative-substitution predicate — where a variant
is covered many ways and the sets genuinely overlap — drops in without a rewrite. When it does, it
becomes a new record; until then v1's exact-restoration frontier is what `trna-panel` reports, and
its degeneracy is a property of the predicate, not a bug in the engine.
