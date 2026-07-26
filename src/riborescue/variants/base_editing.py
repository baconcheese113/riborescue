"""Which nonsense variants a base editor can reach, as geometry, on the declared editor panel.

Reachability is whether an editor can be *placed* on the premature stop: a compatible PAM at the
right distance and the stop-removing base inside the editing window. It is not eligibility. Editing
efficiency, off-target activity, product purity, delivery, tissue access, splice consequence and
toxicity are out of scope, and no field here implies otherwise.

Guides are genomic. Within the exon that holds the stop the transcript-sense strand is one genomic
strand and its complement is the other, so a PAM found on either is a real PAM; a placement whose
protospacer or PAM would cross an exon junction is refused, because there the spliced sequence is
not contiguous DNA. Coding consequences of the edited base and every in-window bystander are read in
the transcript's own reading frame.
"""

from dataclasses import dataclass
from enum import StrEnum
from typing import cast

import pandas as pd
from Bio.Seq import Seq

from riborescue.variants.transcripts import STOP_CODONS, TranscriptModel

__all__ = [
    "PRIMARY_PANEL",
    "SENSITIVITY_PANEL",
    "Bystander",
    "Editor",
    "Guide",
    "Reach",
    "ReachClass",
    "escape_map",
    "reachability_for",
    "reachability_table",
]

PROTOSPACER_LEN = 20
_COMPLEMENT = str.maketrans("ACGT", "TGCA")
_MARGIN = PROTOSPACER_LEN + 6


def _revcomp(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


@dataclass(frozen=True)
class Editor:
    """A base editor: the transition it makes on its protospacer strand, its window, and its PAM.

    `window` is 1-indexed from the PAM-distal end of the 20-nt protospacer, inclusive. `pam` is the
    rule the three nucleotides immediately 3' of the protospacer must satisfy. `arm` names whether
    the editor is part of the primary panel or a labelled PAM-relaxation sensitivity arm.
    """

    name: str
    edit_from: str
    edit_to: str
    window: tuple[int, int]
    pam: str
    arm: str


# Named enzymes, canonical windows, SpCas9 NGG. The two windows are not identical and are not
# rounded to a shared range.
PRIMARY_PANEL = (
    Editor("BE4max", "C", "T", (4, 8), "NGG", "primary"),
    Editor("ABE7.10", "A", "G", (4, 7), "NGG", "primary"),
)

# PAM relaxation. Each is a named enzyme with its own rule and a citation pinned in the docstring;
# reported as its own arm and never merged into the primary count.
SENSITIVITY_PANEL = tuple(
    Editor(f"{base.name}+{fusion}", base.edit_from, base.edit_to, base.window, pam, "sensitivity")
    for base in PRIMARY_PANEL
    for fusion, pam in (("SpCas9-NG", "NG"), ("SpG", "NGN"), ("SpRY", "NRN"))
)
"""SpCas9-NG (Nishimasu 2018): NG · SpG (Walton 2020): NGN · SpRY (Walton 2020): NRN>NYN."""

_PAM_BASES = {
    "N": set("ACGT"),
    "G": {"G"},
    "R": {"A", "G"},
    "Y": {"C", "T"},
}


def _pam_ok(observed: str, rule: str) -> bool:
    """Whether the three nucleotides 3' of a protospacer satisfy a PAM rule.

    A two-letter rule (NG) is read against the first two nucleotides; the third is unconstrained.
    """

    return len(observed) >= len(rule) and all(
        observed[i] in _PAM_BASES[symbol] for i, symbol in enumerate(rule)
    )


class ReachClass(StrEnum):
    """The categorical base-editing reachability of a variant under the declared panel.

    Readthrough is not a category here — it is a continuous overlay. The negative class is named for
    the panel, not for the absence of any route.
    """

    exact = "base_editable_exact"
    alternative = "base_editable_alternative"
    none = "not_base_editable_under_panel"


@dataclass(frozen=True)
class Bystander:
    """An in-window base, other than the target, an editor would also change."""

    offset: int
    aa_from: str
    aa_to: str
    silent: bool
    coding: bool


@dataclass(frozen=True)
class Guide:
    """One placement that removes the stop: the protospacer, its PAM, and what it edits."""

    editor: str
    arm: str
    strand: str
    protospacer: str
    pam: str
    window_position: int
    restores: str
    new_codon: str
    bystanders: tuple[Bystander, ...]

    @property
    def bystander_free(self) -> bool:
        return not any(not b.silent for b in self.bystanders)


@dataclass(frozen=True)
class Reach:
    """Every stop-removing guide found for one variant, and the class they sum to."""

    transcript_id: str
    protein_position: int
    stop_type: str
    reference_codon: str
    guides: tuple[Guide, ...]
    reason: str = ""
    """Why no guide was found, when there is none: no_stop_removing_edit, target_outside_window, or
    no_pam. Empty when the variant is reachable."""

    @property
    def reach_class(self) -> ReachClass:
        if any(g.restores == "exact_wildtype" for g in self.guides):
            return ReachClass.exact
        if self.guides:
            return ReachClass.alternative
        return ReachClass.none

    @property
    def bystander_free(self) -> bool:
        """True when any stop-removing guide edits no amino acid but the target's."""

        return any(g.bystander_free for g in self.guides)


def _aa(codon: str) -> str:
    return str(Seq(codon).translate())


def _codon_consequence(
    sequence: str, index: int, new_base: str, coding: int, coding_end: int
) -> tuple[str, str, bool]:
    """The amino-acid change an edit at a transcript offset makes, in the frame from `coding`.

    A base past the natural stop is non-coding: no residue, and never counted as a missense
    bystander.
    """

    if index < coding or index > coding_end:
        return ("", "", False)
    start = coding + ((index - coding) // 3) * 3
    old = sequence[start : start + 3]
    if len(old) < 3:
        return ("", "", False)
    edited = old[: index - start] + new_base + old[index - start + 1 :]
    return (_aa(old), _aa(edited), True)


def _placements(
    seq: str, junctions: tuple[int, ...], target: int, editor: Editor
) -> list[tuple[int, str, int]]:
    """Protospacer starts placing `target` in the editing window on the forward strand of `seq`.

    Returns each start, its PAM, and the target's 1-indexed window position. A placement whose
    protospacer or PAM crosses an exon junction is dropped.
    """

    lo, hi = editor.window
    found = []
    for start in range(target - hi + 1, target - lo + 2):
        if start < 0 or start + PROTOSPACER_LEN + 3 > len(seq):
            continue
        pam = seq[start + PROTOSPACER_LEN : start + PROTOSPACER_LEN + 3]
        if not _pam_ok(pam, editor.pam):
            continue
        span_end = start + PROTOSPACER_LEN + 3
        if any(start < j < span_end for j in junctions):
            continue
        found.append((start, pam, target - start + 1))
    return found


def _bystanders(
    window: str,
    sense_lo: int,
    sense_hi: int,
    sense_from: str,
    sense_to: str,
    target: int,
    coding: int,
    coding_end: int,
) -> tuple[Bystander, ...]:
    """Every in-window sense base, other than the target, the editor would also change.

    Consequences are read in the transcript's own frame, whichever strand the protospacer sits on.
    """

    out = []
    for offset in range(sense_lo, sense_hi + 1):
        if offset == target or window[offset] != sense_from:
            continue
        aa_from, aa_to, coding_here = _codon_consequence(
            window, offset, sense_to, coding, coding_end
        )
        out.append(
            Bystander(
                offset,
                aa_from,
                aa_to,
                silent=aa_from == aa_to or not coding_here,
                coding=coding_here,
            )
        )
    return tuple(out)


def _strand_windows(
    strand: str,
    window: str,
    rc: str,
    length: int,
    junctions: tuple[int, ...],
    rc_junctions: tuple[int, ...],
    target: int,
    editor: Editor,
) -> list[tuple[str, str, int, int, int]]:
    """PAM placements for a target on one strand: protospacer, PAM, window position, and the sense
    index range the editing window covers.

    The PAM geometry is found on whichever strand carries the protospacer — the antisense strand via
    the reverse complement — but the returned window range is always in sense coordinates, because
    the edit's consequence is read there.
    """

    lo, hi = editor.window
    if strand == "sense":
        return [
            (window[s : s + PROTOSPACER_LEN], pam, wp, s + lo - 1, s + hi - 1)
            for s, pam, wp in _placements(window, junctions, target, editor)
        ]
    tr = length - 1 - target
    return [
        (rc[s : s + PROTOSPACER_LEN], pam, wp, length - s - hi, length - s - lo)
        for s, pam, wp in _placements(rc, rc_junctions, tr, editor)
    ]


def _pam_reaches_target(
    strand: str,
    window: str,
    rc: str,
    length: int,
    junctions: tuple[int, ...],
    rc_junctions: tuple[int, ...],
    target: int,
    editor: Editor,
) -> bool:
    """Whether any PAM places the target anywhere in a protospacer, editing window aside.

    A stop-removing edit with no in-window placement but a PAM covering the target somewhere in the
    guide failed on window position, not on PAM availability — the two failures are reported apart.
    """

    widened = Editor(
        editor.name, editor.edit_from, editor.edit_to, (1, PROTOSPACER_LEN), editor.pam, editor.arm
    )
    return bool(
        _strand_windows(strand, window, rc, length, junctions, rc_junctions, target, widened)
    )


def reachability_for(
    model: TranscriptModel,
    position: int,
    ref: str,
    alt: str,
    *,
    panel: tuple[Editor, ...] = PRIMARY_PANEL,
) -> Reach | None:
    """Base-editing reachability of the stop a substitution creates, or None if it makes none.

    The patient's mutant sequence is what an editor acts on, so the alternate allele is written in
    before the search. Each editor is tried on both strands; the strand fixes the sense transition
    and the PAM geometry, while every codon consequence is read in the transcript's reading frame.
    """

    offset = model.offset_of(position)
    coding = model.coding_offset
    coding_end = model.coding_end_offset
    if offset is None or coding is None or coding_end is None:
        return None
    expected, substituted = (
        (ref, alt)
        if model.strand == "+"
        else (ref.translate(_COMPLEMENT), alt.translate(_COMPLEMENT))
    )
    if model.sequence[offset] != expected:
        return None
    covering = model.codon_covering(position)
    if covering is None:
        return None
    stop_start, wildtype_codon = covering
    mutant = model.sequence[:offset] + substituted + model.sequence[offset + 1 :]
    if mutant[stop_start : stop_start + 3] not in STOP_CODONS or wildtype_codon in STOP_CODONS:
        return None

    low = max(0, stop_start - _MARGIN)
    high = min(len(mutant), stop_start + 3 + _MARGIN)
    window = mutant[low:high]
    length = len(window)
    rc = _revcomp(window)
    stop_local = stop_start - low
    coding_local, coding_end_local = coding - low, coding_end - low
    junctions = tuple(j - low for j in model.junction_offsets if low < j < high)
    rc_junctions = tuple(length - j for j in junctions)

    guides = []
    # (strand, editor, target) for each single-letter edit that would leave a sense codon; kept so a
    # failure can say whether it was the PAM or only the window that no guide cleared.
    targets: list[tuple[str, Editor, int]] = []
    for editor in panel:
        for strand in ("sense", "antisense"):
            sense_from, sense_to = (
                (editor.edit_from, editor.edit_to)
                if strand == "sense"
                else (
                    editor.edit_from.translate(_COMPLEMENT),
                    editor.edit_to.translate(_COMPLEMENT),
                )
            )
            for pos in range(3):
                target = stop_local + pos
                if window[target] != sense_from:
                    continue
                stop = window[stop_local : stop_local + 3]
                new_codon = stop[:pos] + sense_to + stop[pos + 1 :]
                if new_codon in STOP_CODONS:
                    continue
                restores = "exact_wildtype" if new_codon == wildtype_codon else "alternative_sense"
                targets.append((strand, editor, target))
                for proto, pam, wp, win_lo, win_hi in _strand_windows(
                    strand, window, rc, length, junctions, rc_junctions, target, editor
                ):
                    guides.append(
                        Guide(
                            editor=editor.name,
                            arm=editor.arm,
                            strand=strand,
                            protospacer=proto,
                            pam=pam,
                            window_position=wp,
                            restores=restores,
                            new_codon=new_codon,
                            bystanders=_bystanders(
                                window,
                                win_lo,
                                win_hi,
                                sense_from,
                                sense_to,
                                target,
                                coding_local,
                                coding_end_local,
                            ),
                        )
                    )

    reason = ""
    if not guides:
        if not targets:
            reason = "no_stop_removing_edit"
        elif any(
            _pam_reaches_target(strand, window, rc, length, junctions, rc_junctions, t, editor)
            for strand, editor, t in targets
        ):
            reason = "target_outside_window"
        else:
            reason = "no_pam"

    return Reach(
        transcript_id=model.transcript_id,
        protein_position=(stop_start - coding) // 3 + 1,
        stop_type=mutant[stop_start : stop_start + 3],
        reference_codon=wildtype_codon,
        guides=tuple(guides),
        reason=reason,
    )


def reachability_table(
    variants: pd.DataFrame,
    models: dict[int, TranscriptModel],
    *,
    panel: tuple[Editor, ...] = PRIMARY_PANEL,
) -> pd.DataFrame:
    """Base-editing reachability of every variant, one summary row each.

    A representative guide is chosen for the row — exact restoration first, then bystander-free,
    then fewest missense bystanders — the full guide set stays available through `reachability_for`.
    The expanded-PAM `SENSITIVITY_PANEL` is a separately requested arm, not the primary count; the
    row's `arm` names which editor reached the variant.
    """

    rows = []
    for variant in variants.itertuples():
        model = models.get(int(cast(int, variant.gene_id)))
        base: dict[str, object] = {
            "variant_id": variant.variant_id,
            "gene_id": int(cast(int, variant.gene_id)),
            "gene_symbol": getattr(variant, "gene_symbol", ""),
        }
        reach = (
            None
            if model is None
            else reachability_for(
                model, int(cast(int, variant.pos)), str(variant.ref), str(variant.alt), panel=panel
            )
        )
        if reach is None:
            rows.append(
                base
                | {
                    "scoreable": False,
                    "reach_class": "",
                    "reachable": False,
                    "reason": "unscoreable_context",
                }
            )
            continue
        best = min(
            reach.guides,
            key=lambda g: (
                g.restores != "exact_wildtype",
                not g.bystander_free,
                sum(not b.silent for b in g.bystanders),
            ),
            default=None,
        )
        rows.append(
            base
            | {
                "scoreable": True,
                "reachable": bool(reach.guides),
                "reach_class": reach.reach_class.value,
                "reason": reach.reason,
                "protein_position": reach.protein_position,
                "stop_type": reach.stop_type,
                "n_guides": len(reach.guides),
                "n_exact": sum(g.restores == "exact_wildtype" for g in reach.guides),
                "bystander_free": reach.bystander_free,
                "editor": best.editor if best else "",
                "arm": best.arm if best else "",
                "strand": best.strand if best else "",
                "restores": best.restores if best else "",
                "window_position": best.window_position if best else None,
            }
        )
    return pd.DataFrame(rows)


def escape_map(reachability: pd.DataFrame, scored: pd.DataFrame) -> pd.DataFrame:
    """Overlay continuous readthrough on the categorical base-editing class, variant by variant.

    Readthrough is not bucketed into a route: the best predicted readthrough and its interval are
    carried alongside the base-editing class, so a variant reads as its editability *and* its
    readthrough, never sorted by a threshold neither number declares. The negative class stays "not
    base-editable under the declared panel"; nothing here means "no route exists".
    """

    best = (
        scored.sort_values("readthrough_predicted", ascending=False)
        .drop_duplicates("variant_id")
        .loc[
            :,
            [
                "variant_id",
                "therapy_id",
                "readthrough_predicted",
                "readthrough_low",
                "readthrough_high",
            ],
        ]
        .rename(
            columns={
                "therapy_id": "best_therapy",
                "readthrough_predicted": "best_readthrough",
                "readthrough_low": "best_readthrough_low",
                "readthrough_high": "best_readthrough_high",
            }
        )
    )
    return reachability.merge(best, on="variant_id", how="left")
