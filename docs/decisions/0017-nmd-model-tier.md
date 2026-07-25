# ADR-0017 — The NMD model tier: real tools, verified sources, integration order

**Status:** accepted · **Date:** 2026-07-24 · **Deciders:** Joseph

## Context

ADR-0016 built the rule tier (`guideline` and `full_rules`) and named three real predictors to
complete the §7.3 disagreement atlas — **aenmd**, **predNMD**, and **NMDetective-AI** — but left them
unwired until their sources, licences, inputs and any training-set overlap were verified from primary
sources rather than assumed. That verification is now done. Each fact below carries its primary source;
nothing here rests on a preprint summary a fetcher paraphrased.

### The tools, verified

- **aenmd** — a rule-based R package (Bioinformatics 2023) that annotates variants for NMD escape by
  the classic rules on a chosen transcript set. It is the license-clean, CPU-only, download-light tool
  and reproduces the rule tier's logic through an independent implementation — a genuine second opinion
  on the same rules, with its own edge-case handling. *(Repo, licence, data package and exact rule
  thresholds pinned in the integration notes; it is the first tool integrated.)*

- **predNMD** — Su & Brenner, UC Berkeley; bioRxiv 2026, DOI `10.64898/2026.06.20.733449`; repo
  `github.com/BrennerLab/predNMD`. A **random forest** trained on 5,304 nonsense variants from GTEx,
  TCGA, GEUVADIS and GREGoR, predicting a **continuous NMD-trigger probability (0–1)** plus a
  truncation-rescue class for variants it calls escape. **Two licences:** the *code* is MIT (SPDX:
  `MIT`, redistributable); the *manuscript* is CC BY-NC-ND (constrains the paper, not the software).
  GRCh38 (Ensembl release 104), variants keyed by forward-strand `chrom/pos/ref/alt` per transcript.
  Precomputed predictions for all **13,968,776** GRCh38 stop-gain SNVs are advertised at `predNMD.org`
  (→ `compbio.berkeley.edu`), but that host was unreachable to verification (unverifiable TLS
  certificate), the GitHub repo ships **no releases** and no table download, and no Zenodo/Figshare
  mirror was found — so the precomputed-table hosting, format and size are **unconfirmed**.

- **NMDetective-AI** — Veiner, Toledano, Palou-Márquez, Lehner & Supek; bioRxiv 2026, DOI
  `10.64898/2026.03.24.714003`; repo `github.com/Vejni/NMDetectiveAI` (handle `Vejni` = Marcell
  Veiner, Supek lab, IRB Barcelona — *not* Charles Vejnar). A deep model: a fine-tuned **Orthrus**
  RNA foundation encoder (Mamba-based, PyTorch), weights `models/NMDetectiveAI.pt`, licence **MIT**.
  Inputs are transcript sequence + 1-based PTC transcript position via GenomeKit on GRCh38. Trained on
  ~14,000 TCGA/GTEx somatic PTCs with allele-specific expression plus deep-mutational-scanning reporter
  libraries; tested on ~1,800 germline PTCs. Distinct from the **rule-based NMDetective-A/-B** decision
  scorers (Lindeboom, Supek & Lehner, *Nat. Genet.* 2019, whose genome-wide scores are on Figshare).

### The leakage question, resolved

ADR-0016 flagged a shared-provenance hazard: Toledano and Supek/Lehner authored both the readthrough
labels *and* NMDetective-AI, so training overlap had to be checked before using the model as an NMD
feature beside those labels. The verification result is favourable and specific:

- **NMDetective-A/-B and NMDetective-AI train on TCGA/GTEx/GEUVADIS allele-specific expression and DMS
  reporters — not ClinVar, not the readthrough panel.** Scoring NMD on our ClinVar nonsense variants
  therefore does **not** leak into any of these models' training. The NMD tier is clean.
- The real overlap is elsewhere and belongs to a different layer: the **Toledano readthrough panel is
  built from ~5,800 ClinVar pathogenic stop-gains**, so evaluating *readthrough* predictions on ClinVar
  variants shares population with that training set. That is an ADR-0012 / readthrough-evaluation
  concern — hold out or de-duplicate by variant there — and it is independent of the NMD tier.

## Decision

**Integrate the tools individually, cleanest first, and never substitute one for another.** A tool that
is genuinely blocked (licence, unavailable artifact, incompatible transcript definition, unverifiable
source, or an inference stack that does not run in this CPU-only environment) is recorded as an explicit
blocker with evidence, its column stays `unavailable` with a reason, and the remaining tools proceed —
the block is not spent retrying one model.

Order, by tractability:

1. **aenmd** — pinned, CPU-only, license-clean; a real independent implementation of the same rules.
   **Done.** Over the 48,148 stops aenmd scores, the hand-rolled `full_rules` agrees with it on
   99.36%; the 32% aenmd does not score is reported by cause (a build-version gap for most, its own
   splice filtering for the rest). See `variants/aenmd.py`, `scripts/aenmd_nmd.R`, `feature.aenmd`.

2. **predNMD** — **blocked on the tractable path.** The precomputed prediction tables (all
   13,968,776 GRCh38 stop-gain SNVs) are the ideal — a keyed lookup, no inference — but they are
   distributed "available upon request" from `compbio.berkeley.edu/proj/prednmd` (the host serves an
   untrusted TLS certificate; the page's own download tab lists both the GRCh37 and GRCh38 tables as
   *available upon request*, not as files). The GitHub code is MIT and runnable, but its random forest
   needs a heavy feature bundle first — Ensembl-104 genome/GTF/CDS + VEP, gnomAD LOEUF, a phyloP
   bigWig (~10 GB), plus m6A density and TranslationAI features — a pipeline disproportionate to this
   pass and with real feature-parity risk against the authors' exact computation. Recorded as blocked;
   the follow-up is an email for the tables or a scoped build of the feature pipeline.

3. **NMDetective-AI** — **integrated, running on the local RTX 4070 Ti.** Its CUDA stack lives in the
   `nmdetective` Pixi feature (PyTorch-GPU, GenomeKit, `mamba-ssm`/`causal-conv1d` from conda-forge —
   no local compilation), with the Orthrus HuggingFace pins and the Git-LFS weights installed from a
   pinned commit by `scripts/install_nmdetective.sh`. `scripts/nmdetective_predict.py` encodes each
   ClinVar variant on its MANE transcript into the model's six-track input and scores it; the full pass
   is a resumable background job (GenomeKit-encoding-bound, GPU near-idle). Unlike the rules and aenmd,
   NMDetective-AI emits a continuous NMD-*efficiency* score, not a verdict, so nothing thresholds it:
   the atlas reports the separation between the mean efficiency of rule-escape and rule-decay stops,
   which is positive (escaping stops score lower), so the deep model tracks the rules' direction
   without sharing their labels. See `variants/nmdetective.py`.

**Provenance and honesty.** Every fetched artifact is pinned by URL + checksum with a fetch script;
long inference is resumable and atomic per batch, verifies its output count, and refuses stale or mixed
output. Each tool's prediction rides its own column with an availability/reason field; the atlas reports
pairwise and three-way disagreement over the tools actually produced. **Rule-tier and model-tier results
stay visibly distinct** in every surface — the model tier never silently backfills a rule verdict, and a
tool that did not run is shown as not-run, never as agreement.

**Thresholds are frozen before results are seen.** The rule-tier constants (50 / 150 / 407 nt, ADR-0016)
and each model's decision threshold are fixed in code and cited; none is tuned after a distribution is
inspected.

## Consequences

The disagreement atlas grows from a two-rule split toward the three-to-five-way tool disagreement §7.3
describes, with the leakage question answered rather than deferred. Where a tool is blocked, the atlas is
smaller but honest about why, and the record here is enough for a later pass to finish it. The one thing
that does not happen is a manufactured ensemble: no averaged "NMD score", no tool standing in for
another, no rule verdict wearing a model's name.
