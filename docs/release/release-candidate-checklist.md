# Release-candidate checklist

Everything above the line is reviewable and reversible. Everything below it is not, and none of it
runs without Joseph's explicit approval.

## Gates

- [ ] `pixi run check` — lint, types, the test suite, control tests — exit 0
- [ ] `pixi run nf-lint` — Nextflow lint clean
- [ ] `pixi run nf-test` — pipeline tests pass
- [ ] `cd frontend && npx tsc --noEmit` — TypeScript clean
- [ ] `pixi run app-build` — static export builds
- [ ] `pixi run app-preview` — every page renders at 390 px and 1280 px with no horizontal overflow
- [ ] `pixi run provenance` — manifest regenerates with an empty `missing` list

## Outstanding before release

- [ ] **ADR-0023's own acceptance criteria.** The scope is accepted; four of the items the record
      requires do not exist yet — a hand-checked known-variant panel, a per-variant `reason` field, a
      reachable sensitivity arm, and a minus-strand test. Either build them or state in the release
      notes that the layer ships ahead of them. See `docs/release/adr-0023-scope-audit.md`.

## Content review

- [ ] `docs/release/release-notes.md` — the boundary list still matches what the surfaces say
- [ ] Every number in the release documentation traces to a file in `results/` — see
      `docs/release/numbers-verified.md`
- [ ] `results/rnaseq` observations are described as exploratory and treatment-confounded, never as
      differential expression

## Boundary audit — must hold on every surface

- [ ] No page or document calls a therapy score a clinical recommendation
- [ ] The G418 confirmation is stated as within one laboratory and protocol family
- [ ] Native-stop occupancy is never described as protein production, toxicity, or safety
- [ ] GSE179274 is described as inconclusive under the frozen gate, never as negative
- [ ] Wave 2 is described as a small upstream-residue interaction, never as transferable kinetic
      information
- [ ] predNMD's absence is stated wherever the NMD ensemble is described
- [ ] Base-editing reachability is labelled geometric and exploratory wherever it appears

## Packaging

- [ ] `results/provenance.json` regenerated at the release commit
- [ ] `pixi.lock` unchanged since the last full gate run
- [ ] Release notes final

---

## Requires explicit approval — do not run unprompted

- [ ] Create the `v0.1.0` git tag
- [ ] Deposit artifacts to Zenodo and mint a DOI
- [ ] Upload large artifacts to Cloudflare R2
- [ ] Publish the viewer to its public URL

A Zenodo deposition is permanent and a DOI cannot be withdrawn, only superseded. Nothing in this
section runs until every item above the line is checked and the open ADR-0023 decision is made.
