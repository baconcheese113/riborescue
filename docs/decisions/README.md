# Decision records

Every major technical and scientific choice in RiboRescue is recorded here, one file per decision,
in the order they were made. A record states what was decided and why; it is not updated when the
code moves on. Where a later record changes an earlier one, both say so.

## Start here

Five records carry most of what the project claims. Read them in this order and the rest are
context.

| | |
|---|---|
| [0010](0010-frame-control.md) | The control the whole Ribo-seq arm rests on: what a readthrough compound must do to endogenous ribosomes, and what its negative control must fail to do. |
| [0011](0011-length-selection.md) | Which footprint lengths a dataset keeps, chosen from its own periodicity before any contrast is run. The reason `run_psite.R` cannot be replaced by an upstream tool. |
| [0020](0020-codon-occupancy-and-the-kinetics-head-to-head.md) | What a kinetic feature could add over the published sequence model, and the four fits and three shuffles that were named before any of them ran. |
| [0022](0022-the-stalling-endpoint-reconciled.md) | Why a pre-registered endpoint governs over a stricter rule that only ever lived in code. Read this before changing any predicate that decides a verdict. |
| [0024](0024-the-upstream-boundary-and-what-it-carries.md) | Where `nf-core/riboseq` ends and this repository begins, and why the local read-processing workflow is still here. |

For how the pipeline is put together rather than what it measures, read
[0003](0003-pipeline-orchestration.md) alongside 0024 — its workflow inventory is superseded, but
Groovy composes and Python decides is still the rule.

## All records

| | Status | |
|---|---|---|
| [0001](0001-tooling-and-hosting.md) | accepted | Pixi on `linux-64` for the toolchain, and where the built things are served from. |
| [0002](0002-data-and-artifact-storage.md) | accepted | Nothing large in git: inputs are fetched and verified, artifacts go to Zenodo and R2. |
| [0003](0003-pipeline-orchestration.md) | amended by 0024 | A hand-authored Nextflow pipeline whose processes call the CLI and nothing else, with a validated upstream handoff. |
| [0004](0004-reporter-context-window.md) | accepted | How much sequence either side of a premature stop the measured library actually supports. |
| [0005](0005-scope-exclusions.md) | accepted | What this project will not claim, and what each excluded claim would need instead. |
| [0006](0006-epistasis-evaluation-protocol.md) | accepted | How an interaction between sequence positions is tested without rewarding a model for memorising. |
| [0007](0007-ribo-seq-datasets.md) | accepted | Which public Ribo-seq series are used, and what each one can and cannot support. |
| [0008](0008-psite-calibration.md) | accepted | P-site offsets come from riboWaltz, the published tool, in its own pinned R environment. |
| [0009](0009-readthrough-assay.md) | accepted | The readthrough quantities and their windows, fixed before the answer was known. |
| [0010](0010-frame-control.md) | accepted | The three-part signature a readthrough compound must complete, and the negative control that must fail it. |
| [0011](0011-length-selection.md) | accepted | Footprint lengths are chosen per dataset from periodicity, and the same counts serve the sensitivity arm. |
| [0012](0012-sequencing-scope-freeze.md) | accepted | The sequencing corpus is closed at four series; no dataset is added because it would be interesting. |
| [0013](0013-length-selection-unanimity.md) | proposed | Whether requiring every library to agree on a length is too strict, left unapplied to the data that raised it. |
| [0014](0014-suppressor-panel-coverage.md) | accepted | What it means for a suppressor tRNA design to cover a variant, frozen as exact restoration. |
| [0015](0015-disease-normalization.md) | accepted | ClinVar conditions are keyed on MedGen concepts, with placeholders kept and labelled rather than dropped. |
| [0016](0016-nmd-ensemble.md) | accepted | Nonsense-mediated decay as two named rule predictors that disagree in public, not one verdict. |
| [0017](0017-nmd-model-tier.md) | accepted | The published NMD tools added beside the rules, with their sources and integration order named. |
| [0018](0018-protein-function-layer.md) | proposed | Conservation and domain context as separate evidence, never blended into one impact score. |
| [0019](0019-detectability-arm.md) | accepted | What a dataset too thin to answer its question can still say: how large an effect it could have resolved. |
| [0020](0020-codon-occupancy-and-the-kinetics-head-to-head.md) | accepted | The kinetic hypothesis, its four fits, its three shuffles, and its passing rule, all named in advance. |
| [0021](0021-shuffle-control-as-a-permutation-test.md) | accepted | One shuffle is a draw, not a null: the control becomes a familywise permutation test. |
| [0022](0022-the-stalling-endpoint-reconciled.md) | accepted | A pre-registered endpoint governs over a stricter predicate that was never recorded. |
| [0023](0023-base-editing-reachability.md) | accepted | Whether a base editor can be placed on a premature stop, treated as geometry rather than efficacy. |
| [0024](0024-the-upstream-boundary-and-what-it-carries.md) | accepted | `nf-core/riboseq` owns alignment and QC; `run_psite.R` stays local; the handoff declares only what is read. |

## Writing one

Number it next, state the context, the decision and its consequences, and say what would have to be
true for the decision to be wrong. A record is written before the work it governs, not after — that
is what makes it a pre-registration rather than a summary. When a later record changes an earlier
one, name it in both.
