# ADR-0020 — Codon occupancy, and what the kinetics head-to-head can actually test

**Status:** accepted · **Date:** 2026-07-25 · **Deciders:** Joseph

## Context

Claim 1 of the PRD asks whether Ribo-seq kinetic features add transferable information beyond local
sequence context, measured against the reproduced Toledano baseline. Everything upstream of it is
built: P-site calibration passes on GSE144140, the baseline reproduces to parity, and the three
shuffle controls exist as strict-xfail tests. What has not been fixed is what the comparison is, and
one structural fact has to be settled before any number exists, because it decides the experiment.

**The published design already spans every function of the features it has.** The baseline is

```
RT_binomial ~ 0 + stop_type + down_123nt + up_123nt + stop_type:down_123nt
```

with `up_123nt` and `down_123nt` as 64-level factors. That is all the sequence a scored variant
carries: three nucleotides upstream, the stop, three downstream. Nine bases. A per-codon kinetic
score `k(up_123nt)` is therefore a linear combination of the 64 indicator columns already in the
design, and adding it as a term adds a column the design matrix already contains.
`identifiable_columns` would drop it as rank-deficient, and the fitted model would be identical.
The improvement is not small, or noisy, or dependent on the split: it is exactly zero, by
construction, and it would be zero for any kinetic quantity whatever — measured well, measured
badly, or invented.

So the head-to-head as the build plan states it — "baseline vs baseline + kinetics" — cannot be run
as written. Discovering that after computing a codon table would look like moving the goalposts, so
it is settled here instead.

It also forces an honest restatement of Claim 1. Anything computable from nine nucleotides is a
function of local sequence context. There is no feature of the reporter that is *beyond* local
sequence, because local sequence is all the reporter records. What a kinetic feature can be is a
**structured hypothesis about which function of that sequence matters** — and the saturated
categorical is the comparator that has the information but not the structure. That is a smaller
claim than the PRD's wording, and it is the one the data can carry.

Two places remain where a kinetic quantity is not already spanned:

- **Across the up/down boundary.** The baseline interacts stop type with the downstream triplet and
  nothing else. `k(up) · k(down)` is one column that no combination of the baseline's columns
  reproduces. It is the low-rank version of the saturated 4,096-level `up × down` term that ADR-0006
  records as untestable naively — five degrees of freedom in place of four thousand, with a
  biological reason for that particular shape rather than a regularisation constant.
- **As a replacement rather than an addition.** A 64-level factor estimates every codon separately
  and knows nothing about a level it did not see; a decoding rate shares strength across codons and
  gives an unobserved triplet a value. Whether that helps is a question about generalisation, and it
  can only be asked under a split where levels actually go missing.

The measurement side has its own two decisions. GSE144140's calibrated length set is 21, 22 and
31 nt, and the short lengths are not the same ribosome conformation as the long one — 21–24 nt
footprints are the empty-A-site state, 28–35 nt the full-length state. Pooled over the three DMSO
libraries the selected set is 15.4% length 21, 18.4% length 22 and 66.2% length 31, so a third of it
is a ribosome waiting for a tRNA rather than holding one. That is not obviously wrong for measuring
decoding demand — a waiting ribosome is the paused state — but it is a different quantity from the
full-length one, and which is being reported has to be named. And the libraries that build the table
must not be the libraries whose response is being predicted: constructing a predictor of G418
response from G418-treated occupancy would be circular.

## Decision

**The quantity is codon occupancy, also called a pause score. It is never called dwell time,
decoding rate, translation rate or velocity.** Occupancy is the number of P-sites assigned to a
position, normalised within its transcript. A rate is inferred from occupancy under assumptions
about initiation and flux that this project does not make and cannot check. Existing prose that
calls it dwell time is wrong and is corrected as it is touched; `map.md` and the build plan are
working documents and are not evidence for the name.

### The export

**One site convention, one alternative, both named.** The codon scored is the one in the **A site**,
the position being decoded, because tRNA demand acts there. riboWaltz reports the P site, and the
A-site codon is three nucleotides 3′ of it. The P-site table is exported beside it as a pre-declared
sensitivity, not as a second result to choose between.

**The window is applied to whichever site is being scored**, not to the P site and then inherited.
A codon contributes when its own offset satisfies the ADR-0011 internal-CDS rule — at least 18 nt
from the start codon and at most −16 nt from the stop — so the initiation and termination peaks are
excluded from the quantity as they are from every other denominator in this project.

**Footprint lengths come from the dataset's calibration manifest and nowhere else.** For GSE144140
that is 21, 22, 31 nt, and the export refuses a dataset whose manifest does not pass. The published
28–35 nt window runs as the pre-declared sensitivity arm, exactly as the readthrough assay's does.
Both are reported. Because the two differ in conformation and not only in depth, a disagreement
between them is a finding about which ribosome state carries the signal, and is reported as one
rather than resolved by preference.

**Only untreated and vehicle libraries build the table.** For GSE144140 that is the three DMSO
libraries. No treated library contributes to any feature used to predict response to that treatment
or any other. This is not a preference about noise; it is what keeps the comparison from being
circular, and it is why the second HEK293T control set below is worth having.

**Normalisation is within transcript, before anything is pooled.** A position's occupancy is its
P-site count divided by the mean P-site count per codon over that transcript's own internal coding
sequence, so a transcript's abundance cancels and a highly expressed gene does not become the codon
table. A transcript contributes only when its internal coding sequence carries at least 100 P-sites
over the selected lengths — the same floor `MINIMUM_CDS_PSITES` sets for the readthrough assay,
because one convention is worth more here than a separately optimal one — and at least 64 scorable
codon positions, so a transcript cannot enter on a single window.

**One alignment per read, as everywhere else.** The primary placement only, matching `run_psite.R`.
Isoform identity is not inferred, and a footprint is counted once rather than once per isoform it
fits.

**Replicates are aggregated after the per-library table exists, never before.** The codon table is
computed per library and the per-library tables are written. The aggregate is the mean across
libraries of the per-library codon scores, so a library with more reads does not weight the table,
and the spread across libraries is reported per codon beside the score.

**Codon frequency is reported, not corrected for.** Each codon's score carries the number of
positions and the number of transcripts it was measured over. A rare codon has a noisier score and
that has to be visible, because the features below are lookups into this table and a codon measured
over few positions propagates its noise into every variant carrying it.

**A second control set is admitted only through the same gate.** The `hek293t` untreated libraries
are a separate HEK293T series with no calibration manifest. `select-lengths` is run on them under
ADR-0011 before any codon table is built from them; if they pass they become a pre-declared
robustness arm — the head-to-head is repeated with the codon table built from that series instead —
and if they fail they are reported as refused, as GSE179274 was. Their admission is decided by the
calibration, not by whether the arm agrees.

### The features

**Global-prior features only.** Each variant receives `kinetics.global_prior.a_site_up`, the codon
table's score for its `up_123nt` triplet, and `kinetics.global_prior.a_site_down`, the score for its
`down_123nt` triplet. The upstream triplet is the codon in the P site when the ribosome sits at the
premature stop; the downstream triplet is the first codon decoded if it reads through.

**Locus-observed occupancy is out of the primary comparison, with a reason.** The labels come from a
reporter construct, not from the endogenous gene, so occupancy at the endogenous locus is a
measurement of a different molecule that happens to share a gene name. Any gain it produced would be
gene-identity information, which is precisely what the within-gene shuffle exists to catch. It is
not computed for this comparison. If it is ever wanted, it enters as its own ADR.

**No feature is built from a treated library, and no feature is a function of the label.** Both are
asserted in code rather than remembered.

### The models

Four fits per drug, all named before any is run, all reported whatever they show.

| # | Model | What it isolates |
|---|---|---|
| B | The published four-term GLM, unchanged | The parity baseline |
| K1 | B + `k_up:stop_type` + `k_up:k_down` + `k_up:k_down:stop_type` | Whether decoding demand carries structure the factorised baseline cannot express |
| K2 | B with `up_123nt` replaced by `k_up` | Whether a rate generalises where a 64-level factor has nothing to say |
| S | B + the saturated `up_123nt:down_123nt` term | The capacity comparator: the information K1 uses, without its structure |

K1 is the primary claim and adds six design columns, **five** of which are identifiable beyond the
baseline: the three stop-type-specific `k_up` slopes sum to the marginal `k_up` effect, which the
baseline already spans, so they carry two new dimensions and not three. Only how decoding demand
varies with stop type is new, never its overall level — which is the span argument again, one level
down, and is why the term is written as an interaction rather than as a main effect. S is what ADR-0006 records as failing
naively at 4,096 levels; it is fitted here not to be believed but so that K1's gain can be read
against the gain available to something with all the information and no structure. Model B is refit
unchanged in every round so that Toledano parity holds throughout, and no selection touches an
evaluation fold.

### The evaluation

The grid is ADR-0006's and is fixed here: published random 90/10 rounds, grouped by gene, grouped by
sequence-context cluster, and stratified by support. Every model reports every row. The metric is
squared Pearson correlation on held-out data, per drug, over the authors' ten rounds, with
percentile bootstrap intervals and normalisation against that drug's own replicate reliability —
never against 1.0 and never against a single pan-drug ceiling. The assay's readthrough ceiling is
handled as the authors handle it in the primary, with the saturated observations excluded in a
pre-specified sensitivity.

**Grouped by gene is where the claim lives.** The published random split leaks near-identical
contexts across the fold boundary, and a gain that appears only there is reported as capacity.

### The controls

The three shuffles are respecified for a feature that is a lookup table rather than a positional
profile. Each permutes the codon-to-score assignment; none permutes the labels, and none touches the
sequence columns, so the baseline is identical under every one of them.

| Shuffle | What is permuted | What it catches |
|---|---|---|
| Global | The 64-codon score table, uniformly at random | Any dependence on the measured scores at all |
| Within-gene | Scores among the variants of a gene | Improvement carried by gene identity rather than by codon |
| Context-matched | The score table **within synonymous codon families** | Improvement carried by amino-acid identity rather than by decoding demand |

**Context-matched is the decisive one, and the synonymous restriction is what makes it decisive.**
Codons of the same amino acid encode the same residue and differ mainly in the tRNA that reads them.
Permuting scores within a family therefore preserves the amino acid, the chemistry and most of the
sequence context while destroying exactly the quantity the kinetic hypothesis names. If the
improvement survives it, the gain is amino-acid identity, and Claim 1 is not supported. Families are
taken from the standard genetic code; the three stop codons are not in the table, and the two
single-codon families (Met, Trp) are invariant under this shuffle and are reported as such.

Each shuffle is run over the same rounds as the models it controls, with the same bootstrap
procedure. A shuffle's interval that does not include zero is a leak, and the finding is that the
pipeline is wrong rather than that the shuffle is interesting.

### The passing rule, in full

Kinetics is reported as carrying transferable information only when **all** of the following hold,
per drug:

1. K1's gain over B has a bootstrap interval **excluding zero under the grouped-by-gene split**.
2. All three shuffle controls on K1 have intervals **including zero**.
3. The gain survives the support stratification — it is not confined to the cells where the
   baseline is aliased, which would make it an artefact of what B cannot fit rather than of what
   kinetics knows.

A gain that appears only under the published random split is reported as capacity, not signal. A
gain smaller than S's is reported as such, because a structured hypothesis that underperforms an
unstructured one with the same information has not demonstrated its structure. **A null result is
the expected outcome and is a result**: it would say the published factorisation is adequate for
drug-induced readthrough, which is worth reporting and is not a failure of the pipeline.

The four xfail controls in `tests/controls/` are removed only when their real analyses run and
assert against these frozen criteria. Until then they stay strict-xfail, and the guard that fails
when a graduated control keeps its marker stays as it is.

## Consequences

`run_psite.R` gains a codon-occupancy export and nothing else changes in it. The export is another
pass over the alignments it already reads, so it is added to the same per-sample run rather than
made a separate walk of every BAM, and one heavyweight process runs at a time.

The comparison is now small. Seven columns against roughly five thousand variants per drug is a
cheap fit, and the cost of Wave 2 has moved almost entirely into the export and the controls. That
is the right place for it.

ADR-0006 and this record now describe one experiment from two sides. Its Model 3 — "sparse or
low-rank interaction" — is K1 here, with the rank-1 structure supplied by biology instead of chosen
by validation. Its Model 2 and Model 4 remain W5.3's, unchanged.

Claim 1's wording in the PRD overstates what nine nucleotides can support. The claim this project
tests is the restated one above, and the report says so rather than quietly narrowing it.
