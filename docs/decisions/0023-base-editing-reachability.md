# ADR-0023 — Base-editing reachability, and the cross-modality escape map

**Status:** accepted · **Date:** 2026-07-26 · **Decider:** Joseph

## Context

The scored table answers one question — which of six readthrough compounds scores highest for a
nonsense variant — and Wave 2 showed those six are barely separable on efficacy. A variant that no
compound reaches looks, in the table as it stands, the same as one that three reach well. The table
has no way to say *readthrough is the wrong tool here*.

Base editing is a second route to a premature stop, and whether a stop is reachable by it is
**geometry, computable from sequence already fetched**. A cytosine or adenine base editor changes one
target base within a narrow window of a protospacer, positioned by a PAM a fixed distance away. A
premature stop is one substitution from the codon it replaced, so the reverting edit is defined; what
is not defined, and what discriminates, is whether an editor can be *placed* on it — a compatible PAM
at the right distance, the target base in the window, and what else in the window would also change.

This is not a therapeutic-eligibility claim and must never be read as one. Editing efficiency,
off-target activity, product purity, delivery, tissue access, splice consequence and toxicity are all
out of scope by construction. What is in scope is a single, honest geometric question and the
negative space it exposes.

The number that matters — how many pathogenic nonsense variants no modality reaches — is worthless if
it is computed first and defined afterward. **The endpoint is frozen here, before the analysis runs
on the full ClinVar set.**

## Decision

**A predeclared base-editing reachability layer, geometric only, and a cross-modality classification
built on it.** Prime editing is scoped but not built here (see Consequences).

**The primary editor panel names concrete enzymes**, fixed in advance. Their published windows are
not identical and are not rounded to a shared range:

| Editor | Edit | Protospacer | PAM | Window |
|---|---|---|---|---|
| BE4max (cytosine) | C→T on the edited strand | 20 nt | SpCas9 `NGG` | positions 4–8 |
| ABE7.10 (adenine) | A→G on the edited strand | 20 nt | SpCas9 `NGG` | positions 4–7 |

**Numbering, stated explicitly.** The protospacer is the 20 nt immediately 5′ of the PAM on the
edited (protospacer) strand. Position 1 is the PAM-**distal** end; position 20 is adjacent to the
PAM; the `NGG` sits at positions 21–23. Windows above are in this 1-indexed frame. Both strands are
searched — a stop reverted by a G→A or T→C change on the transcript strand is an ABE or CBE edit on
the antisense strand, where the target base reads as A or C.

Guides are **genomic** features — the protospacer and PAM are read from the reference genome around
the variant, and a guide is not designed across an exon junction into a coordinate that does not
exist as contiguous DNA. Coding consequence of every edited base, target and bystander alike, is
evaluated **transcript-aware** on MANE Select.

**Per variant, per candidate guide, the layer reports** — never a single reachable/not flag:

```
editing.reachable                 any panel editor can place the reverting edit
editing.editor                    BE4max | ABE7.10
editing.strand                    + | -
editing.guide                     candidate 20-nt protospacer
editing.pam                       the NGG placing the target in the window
editing.window_position           the target base's position in the protospacer
editing.restores                  exact_wildtype | alternative_sense
editing.bystanders[]              each in-window off-target base: position, ref/alt codon,
                                  aa_from → aa_to, and whether it is silent
editing.bystander_free            no in-window base other than the target changes an amino acid
editing.reason                    when unreachable: no_pam | target_out_of_window | not_a_be_edit
```

`restores` separates an edit that puts back the original residue from one that removes the stop by
installing a different sense codon — the base-editing analogue of near-cognate insertion of a
non-original residue, and reported the same way rather than collapsed. Bystanders are enumerated with
their amino-acid consequences, not counted: a silent bystander and a missense bystander are different
facts, and `bystander_free` is the honest headline a guide can carry.

**Expanded-PAM editors are a separately labelled sensitivity arm**, each a named enzyme with its own
PAM rule, window and pinned citation — never a generic "PAM relaxation":

| Editor | PAM rule | Source (pinned in the module's editor table, verified against source) |
|---|---|---|
| SpCas9-NG | `NG` | Nishimasu et al. 2018 |
| SpG | `NGN` | Walton et al. 2020 |
| SpRY | `NRN` > `NYN` (near-PAMless) | Walton et al. 2020 |

Each pairs with BE4max or ABE7.10 at the windows above. They relax the constraint that most
discriminates reachability, so they are reported as their own columns and never merged into the
primary count.

**The cross-modality map overlays a categorical axis on a continuous one.** Base-editing reachability
is genuine geometry and is categorical:

```
base-editable (exact wild-type restoration)
base-editable (alternative sense only)
not base-editable under the declared panel
```

Readthrough is **not** forced into a binary route. It has no evidence-based reachable/not threshold
declared, and inventing one would smuggle a decision into the map. Instead the continuous readthrough
amenability already produced — predicted readthrough per therapy with its interval — is **overlaid**
on the base-editing class, so a variant reads as *(not base-editable under the panel; best readthrough
1.4% [0.9–2.1])* rather than being sorted into a bucket a threshold invented.

The categorical negative class is therefore **"not base-editable under the declared panel,"** never
"therapeutic negative space." Prime editing is unevaluated, expanded-PAM editors sit in the
sensitivity arm, and readthrough is a continuous overlay — so no cell of this map is entitled to mean
"no route exists." When prime editing is added it becomes a second categorical axis, and only then
does a joint negative class become nameable, with its own predeclaration.

## Consequences

- **Prime editing is scoped, not built.** After base editing works, its integration is time-boxed:
  it pins an established designer (PrimeDesign) rather than a hand-rolled pegRNA heuristic, and if the
  integration proves disproportionate it is left explicitly partial. It enters the map as a second
  categorical axis when it lands, not before.
- **Every surface says geometric reachability.** No view, card or column implies eligibility. The
  excluded quantities are listed wherever the layer is presented.
- The layer reuses the reference genome and MANE annotation already fetched; it adds no dataset and no
  model. Its gate is the full lint/type/test/control suite, plus a hand-checked panel of known
  variants confirming PAM placement and strand on both editors.
- After this unit passes its gate, the scientific scope is closed. Remaining work is presentation,
  reports, the reproducible release, and the forked-repository submission requirement.
- This ADR is **accepted before aggregate coverage is computed**. The per-variant computation may be
  implemented and unit-tested against fixtures beforehand; the ClinVar-wide negative-space number is
  produced only once the endpoint above is accepted.
