# ADR-0023 — release-scope audit

A record of what the base-editing layer computes, what it already publishes, and what is and is not
behind it. It states facts and asks for a decision. It changes no status.

## 1 · What it computes

One question per variant: **can a base editor be physically placed on this premature stop?**

A base editor is a Cas protein that cannot cut, fused to an enzyme that rewrites a single letter.
It only works where three things line up at once:

1. A **PAM** — a short fixed motif the Cas protein must find in the DNA — sits nearby.
2. The letter to be rewritten falls inside a narrow **editing window**, a few positions along the
   20-letter guide.
3. Rewriting that letter actually removes the stop rather than producing another one.

The layer searches both DNA strands, refuses any guide that would have to span an exon junction
(there the spliced sequence is not contiguous DNA), and reads the consequence of every edited
letter in the transcript's own reading frame. Per variant it returns each valid placement, which
editor and strand it uses, where in the window the target sits, whether the result puts back the
**original amino acid** or merely a **different one that is not a stop**, and every other letter in
the window the editor would change along the way (a *bystander*).

`escape_map` then carries the continuous readthrough prediction alongside that class, so a variant
reads as *(not base-editable under the panel; best readthrough 1.4% [0.9–2.1])*. Readthrough is
never bucketed into reachable/not — no threshold for that exists.

**Out of scope by construction, in the module docstring, the CLI help, the caveat string, the
patient card and the researcher panel:** editing efficiency, off-target activity, product purity,
delivery, tissue access, splice consequence, toxicity, clinical eligibility.

## 2 · What is already exposed

`escape` is published in both `riborescue.json` and `riborescue_research.json`, rendered as a
three-bar panel on `/researcher` and as a per-variant "Possible rescue routes" card on `/patient`.

| | count | share of scoreable |
|---|---:|---:|
| ClinVar pathogenic nonsense variants | 71,378 | |
| placed on a stop-forming codon (scoreable) | 70,660 | |
| not placeable | 718 | |
| exact wild-type restoration | 7,622 | 10.8% |
| alternative sense codon only | 14,420 | 20.4% |
| not base-editable under the panel | 48,618 | 68.8% |
| reachable (exact + alternative) | 22,042 | 31.2% |
| reachable by a bystander-free guide | 18,249 | |
| exact restorations that are bystander-free | 6,384 | |

The negative class is named for the panel, not for the absence of any route, on every surface.

## 3 · What was known before it was built

- `docs/prd-updated.md` §7.6 already scoped this: geometric reachability by base and prime editing,
  computable from sequence in hand, labelled as geometry, with the negative space as the output that
  matters. The layer is planned scope, not a new direction.
- ADR-0023 was written before the module and before the ClinVar-wide numbers were produced. The
  endpoint freeze it claims for itself held in that order.
- That a cytosine editor cannot correct a nonsense variant is a fact of the genetic code, not a
  result: no stop codon contains a C, and turning either G of TGA/TAG into an A yields TAA, still a
  stop. This was knowable in advance and is asserted as a test.

## 4 · Evidence behind it, and its limits

**What is tested.** Eight unit tests on hand-built transcripts where every PAM is placed on purpose:
exact restoration of a tryptophan stop with the window position and PAM checked; a cytosine editor
producing no guide; a guide refused for crossing an exon junction; a coding bystander reported with
its His→Arg change; a substitution that makes no stop returning nothing; one row per variant with an
unmodelled gene marked not-scoreable rather than crashing; the panel's identity.

**Five limits, each of which narrows what the 31.2% may be called.**

1. **The two-editor panel is one editor in practice.** BE4max contributes zero guides across all
   70,660 scoreable variants — for the code-level reason above. Every reachable variant is reached
   by ABE7.10. "BE4max + ABE7.10" describes what was searched, not what contributed, and the number
   is an ABE7.10-with-NGG number.
2. **The hand-checked panel of known variants ADR-0023 requires as part of its gate does not
   exist.** Every test transcript is synthetic. Nothing has confirmed a PAM placement or a strand
   call against a real, independently designed guide.
3. **The PAM-relaxation sensitivity arm is defined but unreachable.** `SENSITIVITY_PANEL` exists;
   no CLI flag, export or test selects it, so SpCas9-NG, SpG and SpRY are declared and unrun.
4. **The per-variant `reason` field ADR-0023 specifies is absent.** A variant that fails is recorded
   as unreachable without saying whether no PAM was available, the target sat outside the window, or
   no single-letter edit removes the stop. 68.8% is therefore an undifferentiated denominator.
5. **No transcript on the minus strand appears in any test.** Real minus-strand transcripts are
   scored in production; nothing small and readable checks that path.

Two more facts that belong in a release audit:

- **Nothing was written to `docs/notebook.md`** — this layer left no record of what it cost or what
  went wrong while building it.
- **Two decision records briefly claimed the number 0022.** Base-editing reachability now sits at
  0023 and the stalling-endpoint reconciliation keeps 0022, which is what the falsification ledger
  cites. No scientific content changed in either.

## 5 · Status, and what still has to be built

**The layer is accepted scope, decided by Joseph.** Base editing was wanted, it is built, and it
publishes.

**Accepting the scope does not discharge the gate.** ADR-0023 names what its own acceptance
requires, and four of those items do not yet exist: the hand-checked known-variant panel, the
per-variant `reason` field, a reachable sensitivity arm, and a minus-strand test. Until they do, the
count is arithmetic on declared geometry that nothing outside this repository has checked.

**What the status never changes** is how the layer is described. ADR-0023 itself requires every
surface to say geometric reachability and to name the negative class for the panel rather than for
the absence of a route, and that holds whether the record reads proposed or accepted. Editing
efficiency, off-target activity, delivery, tissue access, splice consequence and eligibility stay
out of scope by construction.

## 6 · Resolution

Four of the five limits in §4 are closed; the fifth is a true description, not a defect.

- **Hand-checked real-variant panel** — added. `test_hand_checked_real_variants_reproduce_their_reachability_call` reads three real ClinVar variants on their real MANE Select transcripts (SAMD11 exact, ISG15 alternative, AGRN not-editable) and confirms the class, editor and strand. It runs wherever the fetched annotation is present and skips in CI, where it is not. Published, independently designed guide/PAM placements are still not among the checks, so the 200-variant run is described as *executed on real data*, not external validation.
- **Per-variant `reason`** — added. A scoreable variant that no guide reaches now carries `no_pam` or `target_outside_window`; an unscoreable one carries `unscoreable_context`. The 48,618 not-editable split into **17,931 no_pam** and **30,687 target_outside_window**. The `no_stop_removing_edit` reason is defined but structurally unreachable for a real premature stop: every stop begins with T, and a T→C edit (an adenine editor on the antisense strand) always yields a non-stop codon, so a stop-removing edit always exists — the same kind of code-level fact as BE4max never contributing.
- **Reachable sensitivity arm** — `riborescue base-editing --sensitivity` runs `PRIMARY_PANEL + SENSITIVITY_PANEL`; the row's `arm` names which editor reached the variant, keeping the primary count separate. `test_a_relaxed_pam_reaches_a_stop_the_ngg_panel_cannot` shows an NG PAM reaching a stop NGG cannot.
- **Minus-strand test** — added. `test_a_minus_strand_transcript_is_scored_on_its_own_reading_frame` scores the same stop through the minus-strand path, where offsets count down and the alleles are complemented.
- **The one-editor characterisation stands.** BE4max still contributes zero guides across the scoreable set; "BE4max + ABE7.10" names what was searched, and the number remains an ABE7.10-with-NGG number. This is reported, not fixed.

A notebook entry now records what building the layer cost. Both ADRs no longer collide: this layer's record is **ADR-0023**; `0022-the-stalling-endpoint-reconciled.md` keeps 0022. The ADR's `Status` line reads `accepted · Decider: Joseph`.

## 7 · The PAM-flexibility sensitivity analysis

Running the sensitivity arm returned a result worth stating carefully, because it is easy to
overstate. It is a **PAM-flexibility geometric sensitivity analysis: the primary BE4max/ABE7.10
editing windows are held fixed while PAM recognition is relaxed** to what SpCas9-NG (`NG`), SpG
(`NGN`) and SpRY (`NRN`) recognise. It is **not** a panel of individually validated editor
architectures — a validated NG-ABEmax construct has its own characterised window, which this
abstraction does not adopt; it varies only the PAM. This interpretation was clarified after the first
sensitivity numbers were seen.

Against the **70,660** scoreable denominator, base-editing placement rises from **31.2% (22,042)**
under canonical NGG to **89.0% (62,869)** under relaxed PAM. The honest incremental is **40,827**
variants that no canonical-NGG guide reaches but a relaxed-PAM guide does — read from the
`requires_relaxed_pam` flag, which is `reachable ∧ ¬primary_reachable`, **not** from the
representative-guide `arm` (that label counts 47,100, because a variant a primary guide reaches can
still pick a cleaner relaxed guide as its representative; a test pins this distinction).

Four things this number is not:

- **`31.2% → 89.0%` is geometric placement, not expected editing.** Activity, specificity, delivery
  and tissue compatibility are unmodelled. A placed guide is a candidate, not an outcome.
- **`NRN` is conservative only within this fixed-window abstraction.** SpRY also has lower `NYN`
  activity, not searched here, so 89.0% could rise under this same geometry — it is not an
  experimental lower bound.
- **The residual 7,791 is not entirely a window limit.** 6,667 are `target_outside_window` and 1,124
  lack any modelled PAM; both are stated.
- **It does not merge into the primary count.** The headline base-editing number stays 31.2% NGG;
  the 89.0% is the labelled sensitivity arm.

The finding that survives all of that: the canonical-NGG negative space is **largely a PAM-stringency
artifact of the geometry, not a fundamental limit** — which is a stronger reason to leave prime
editing deferred, since expanded-PAM base editing already spans most of the gap.
