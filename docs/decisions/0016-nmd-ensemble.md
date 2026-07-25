# ADR-0016 — The NMD layer: a rule ensemble now, model predictors specified

**Status:** accepted · **Date:** 2026-07-24 · **Deciders:** Joseph

## Context

Whether a premature stop escapes nonsense-mediated decay sets how much mutant transcript survives to
be read through, and it is Layer 2 of the PRD — Claim 2 asks *"do independent NMD models agree?"*, the
variant-curator persona's whole value is *"NMD efficiency with disagreement exposed"*, and §7.3 names
an NMD disagreement atlas. It is also the weakest load-bearing thing in the product: `escapes_decay`
has rested on a single rule — the 50-nt rule (escape if the PTC is in the last exon or within 55 nt of
the last exon-exon junction) — which is coarse: the fuller Lindeboom rule set classifies many more
stops as escaping, and where the two split the single rule gives a curator no signal that the call is
contested. ClinGen's PVS1 decision tree uses that same 50-nt rule to decide whether a PTC triggers NMD
and thus PVS1 strength, so the boundary has clinical reach.

The PRD names three predictors: **NMDetective-AI** (`Vejni/NMDetectiveAI`, an Orthrus-encoder model,
bioRxiv Mar 2026), **predNMD** (a random forest precomputed for all 13,968,776 GRCh38 stop-gain SNVs,
Docker + tables, bioRxiv Jun 2026), and **aenmd** (rule-based, transcript-dependent, an R/Bioconductor
package, Bioinformatics 2023). It also records a shared-provenance hazard: Toledano, Supek and Lehner
authored both the readthrough labels *and* NMDetective-AI, so training-set overlap must be checked
before that model is used as a feature beside those labels.

## Decision

**Build the rule tier now, from transcript geometry already in hand.** The context step already places
each PTC on its MANE Select transcript; extended to carry distance from the start codon and the length
of the PTC's exon, that geometry supports the full classic rule set with no download. Two named
predictors are computed, each deterministic:

- **`guideline`** — the 50-nt rule as ClinGen PVS1 uses it: escape if the PTC is in the last exon or
  within 50 nt (5') of the last exon-exon junction.
- **`full_rules`** — the Lindeboom–Supek–Lehner rule set (Nat. Genet. 2016): the guideline conditions
  *plus* start-proximal escape (PTC within 150 nt of the start codon) and long-exon escape (PTC in an
  exon longer than 407 nt).

The per-rule flags are exposed, not just the verdicts. Because `full_rules` is a strict superset of
`guideline`, the two disagree in one direction only, and the disagreement set is exactly the
start-proximal and long-exon escapes the guideline does not call. On the full ClinVar set (70,386
scoreable stops) the guideline escapes 10.1%, the full rule set 30.3%, and the two split on 20.2%
(14,231 stops — long-exon 10,593, start-proximal 2,933, both 705). The split is a change of
classification, not a measured error rate: which rule is right at a split needs an empirical NMD
benchmark, which the model tier below is meant to supply. Thresholds (50 nt, 150 nt, 407 nt) are the
published ones and are constants, cited in the module; the last-junction distance is the 50 nt
ClinGen PVS1 (Abou Tayoun et al. 2018) applies, the conservative end of Nagy & Maquat's "50-55 nt".

**Do not integrate the model predictors in this pass.** NMDetective-AI and predNMD stay specified, not
wired, until their licensing, versions and inputs are verified and — for NMDetective-AI — the overlap
with the readthrough labels is measured. predNMD's precomputed tables are the most tractable next
addition (a lookup, no inference), but their size and licence are verified before they are fetched;
aenmd is rule-based and would largely reproduce `full_rules`. No model is downloaded or stubbed here:
the ensemble is honestly the two rule predictors, and the ADR names what completes it.

**Transcript selection is MANE Select**, matching every other layer, so NMD status is read on the same
transcript the readthrough and residue layers use. A per-transcript sensitivity analysis is deferred
with the model predictors.

## Consequences

The curator column and the patient NMD slot stop resting on one rule: they show a verdict from each
predictor, which rules fired, and whether the two agree. The disagreement atlas is the aggregate of
that — how many variants the guideline and the full rule set split on, and which rule drives each
split — and it is honest that it is a *rule* disagreement, the boundary between the guideline and the
fuller rules, not yet the three-model disagreement the PRD's §7.3 ultimately wants. That larger atlas
arrives when NMDetective-AI and predNMD are integrated under a follow-up record, with the overlap
check done. Until then nothing here claims to be those models, and the single 50-nt rule is no longer
the only voice.
