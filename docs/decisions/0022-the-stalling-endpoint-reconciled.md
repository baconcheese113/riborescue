# ADR-0022 — The stalling endpoint, reconciled

**Status:** accepted · **Date:** 2026-07-26 · **Deciders:** Joseph

## Context

ADR-0010 fixed what the negative control must show for the frame control to be confirmed:
SRI-37240 against DMSO fails the readthrough signature **in the direction its mechanism predicts —
termination occupancy does not fall, and downstream occupancy does not rise.** Two conditions, both
on means.

The implementation did not stay there. Read against the history, on 2026-07-22:

| time | event |
|---|---|
| 02:05 | ADR-0010 commits the two-condition endpoint |
| 02:14 | `stalling()` is written with exactly those two conditions, quoting the record |
| 02:23 | a third is added — the termination arms must separate completely across replicates |

The commit that added the third condition says why: *"a stalling verdict now needs the groups to
separate rather than only their means to differ."* It was deliberate, and it landed twenty-three
hours before the commit that made the SRI-37240 contrast computable at all, so it is not post-hoc.
It was also never written into any decision record — no ADR in this repository states it.

The divergence surfaced when the corrected 3' UTR exclusion prompted a rerun. On GSE144140 the two
definitions disagree: SRI-37240 gives termination +0.00092 and downstream −0.00024, meeting
ADR-0010's endpoint, while its arms overlap across replicates (0.0704, 0.0574, 0.0610 against
0.0667, 0.0628, 0.0566) and so fail the stricter rule.

Both definitions predate the data, so neither can be dismissed as fitted to the result. What decides
between them is authority, not chronology.

## Decision

**ADR-0010 governs this dataset.** The endpoint is the two mean-direction conditions it states, and
`stalling()` is restored to them.

A pre-registration is only worth what it forbids. This repository assigns decision authority to
`docs/decisions/` alone; if an implementation change that was never recorded could redefine a
declared endpoint, then no endpoint is declared and pre-registration protects nothing. That the
stricter rule was deliberate and pre-data makes it a defensible rule — it does not make it the
registered one. Chronology cannot settle authority, because the later definition is later only in
the artefact that has no authority to settle it.

**Replicate separation is retained as a named diagnostic, not a verdict.**
`termination_arms_separate` reports whether the *termination* arms separate completely — named for
the quantity it reads, so it cannot be mistaken for separation of downstream occupancy or of the
signature as a whole. It sits beside the stalling verdict and is never folded into it. A
directional result whose arms overlap is weaker evidence than one whose arms do not, and that
difference is worth reporting rather than either hiding or promoting to a threshold.

**The frame control is confirmed, with the overlap disclosed.** G418 completes the signature;
SRI-37240 fails it in the predicted direction. The web export states both facts and does not
describe the control as unresolved.

**Any stronger endpoint is adopted prospectively.** A three-condition rule requiring separation may
be the better rule. If it is to govern, it is recorded before the dataset it judges, not after.

## Consequences

The SRI37240-STALLING ledger row stays supported, and gains the qualification that the registered
endpoint passes while the separation diagnostic does not. The G418-READTHROUGH row is unaffected:
ADR-0010's confirmation is a conjunction, and both limbs hold.

Two libraries per arm more would settle the diagnostic on this design. Nothing about the registered
endpoint waits on it.

The general failure is worth naming, because it was invisible for four days: an endpoint written
in prose in one record and in code in another will drift, and ordinary implementation tests cannot
detect disagreement with an independently maintained prose endpoint. A registered acceptance example
can catch a mismatch once someone knows to look for it; prose and code still cannot validate each
other automatically. Where a decision record states a rule the code implements, the code should cite
the record — as `stalling()` now does.
