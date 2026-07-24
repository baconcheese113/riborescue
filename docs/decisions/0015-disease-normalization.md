# ADR-0015 — Disease normalization from ClinVar's own cross-references

**Status:** accepted · **Date:** 2026-07-24 · **Deciders:** Joseph

## Context

Coverage has so far been counted over variants and genes. To count it over *diseases* — which
diseases have a model-covered variant, and what fraction of each disease's nonsense variants are
reached — the conditions attached to each variant need stable identifiers and a defensible source.
The scored table already carries ClinVar's condition *names* (CLNDN), but names are not keys: the same
disease appears under several spellings, and a name cannot be joined to a coverage denominator.

ClinVar's VCF carries the identifiers alongside the names. `CLNDISDB` is a list parallel to `CLNDN`,
one entry per asserted condition, each a comma-separated set of `Source:ID` cross-references —
`MedGen:C0035334,OMIM:268000,Orphanet:791,MONDO:MONDO:0019200,MeSH:D012174`. Placeholder assertions
("not provided", "not specified") carry a MedGen `CN*` concept rather than a real disease. So the
mapping from a scored variant to MedGen, OMIM and Orphanet is already inside the release we pinned.

This matters for licensing. OMIM's database is licensed and not freely redistributable; Orphadata is
CC BY 4.0. But taking the cross-reference *identifiers ClinVar publishes* is not importing either
resource — NCBI distributes ClinVar, including `CLNDISDB`, as public-domain data. Using it needs no
separate download, no license, and no second version to pin, and it is aligned by construction to the
exact variant set being scored.

## Decision

**Source.** Disease identifiers come from `CLNDISDB`/`CLNDN` in the same pinned ClinVar release
(`clinvar_20260715.vcf.gz`, checksum-verified by the fetch) that defines the variant set. No external
mapping file is downloaded; nothing licensed is redistributed. OMIM and Orphanet appear only as the
identifiers ClinVar itself publishes.

**Identifier precedence.** The MedGen CUI is the normalized primary key. It is NCBI's public-domain
concept hub, present on essentially every ClinVar condition, and the natural join key for a coverage
denominator. OMIM, Orphanet, MONDO and MeSH are retained as attached cross-references, never as the
key. A disease is one MedGen concept.

**One-to-many is preserved, never resolved silently.** A variant that asserts several conditions
becomes several rows. A condition that carries several OMIM or Orphanet ids keeps them all, joined,
not reduced to one. The normalization never picks a single disease for a multi-condition variant.

**Unmapped and ambiguous records are kept with an explicit reason.** Every condition row carries a
`mapping_status`: `placeholder` for the `CN*` "not provided/specified" concepts (a real MedGen id but
not a disease), `medgen_only` for a real concept with no OMIM or Orphanet cross-reference, and
`mapped` for one that has at least one. None is dropped; the completeness of the crosswalk is itself
reported data, not a filter applied out of sight.

**Terminology: condition entity, not disease.** A MedGen concept reached this way is a ClinVar
*condition*, which is not always a disease — the set includes findings, susceptibilities, and broad
umbrella labels ("Hereditary cancer-predisposing syndrome", "Cardiovascular phenotype"). Counts and
UI therefore say *condition entities*, and a headline never calls the 6,058 MedGen concepts
"diseases" unless each has been verified as one, which they have not.

## Consequences

Condition-entity aggregation keys on MedGen, and its denominators (ADR-0016-style coverage semantics)
are computed per MedGen concept over the non-placeholder rows. Provenance is exactly the ClinVar
release already pinned — no licensed input, no second source to version, and the crosswalk cannot
drift out of step with the variant set because both are read from one file.

The cost is that the crosswalk is only as complete as ClinVar's submitters made it: some conditions
resolve to MedGen alone, and the placeholder concepts carry no disease at all. That incompleteness is
surfaced per row rather than hidden, so a disease coverage number is always reported beside how much
of its denominator could be mapped. Should a later need require OMIM phenotype detail or Orphanet
prevalence, importing those resources is a separate decision with its own licensing record (Orphadata
CC BY 4.0 with attribution; OMIM under its own terms) — this record covers only the identifiers
ClinVar publishes.
