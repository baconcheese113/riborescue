# ADR-0013 — The length-selection unanimity veto is too strict

**Status:** proposed · **Date:** 2026-07-23 · **Deciders:** Joseph

## Context

ADR-0011 keeps a footprint length for a dataset only when frame 0 strictly outnumbers both other
frames in *every* library. On GSE144140 that rule selected lengths 21, 22 and 31 and excluded 30 —
the most abundant footprint, whose frame-0 count is the plurality in all nine libraries and the
outright winner in seven. Two libraries, where frame 2 edged frame 0 by a fraction of a percent at
length 30, vetoed it for the whole dataset. The rule behaves as a unanimous veto over nine noisy
measurements, and it discards a length the periodicity plainly supports while keeping two
empty-A-site lengths that happen to separate cleanly.

The GSE144140 result did not depend on this: the 28–35 nt sensitivity analysis, which contains 30
and none of the short lengths, reproduced the downstream and frame effects with intervals excluding
zero. But the primary arm ran on a set no one would choose, and the next dataset deserves better.

## Decision

This record only names the limitation. It changes nothing, and in particular it does not re-analyse
GSE144140, whose result stands under the rule as it was frozen. Amending the rule now, with this
result in hand, would be choosing a length rule by the answer it produces — the exact move the
pre-registration exists to forbid.

A future amendment is expected to replace the unanimous requirement with a tolerant one — a majority
of libraries, or periodicity judged on the pooled counts rather than library by library — decided
and frozen before it is applied to any dataset that has not yet been analysed. Because the sequencing
scope is frozen (ADR-0012), no dataset is waiting on it, so it blocks nothing.

## Consequences

`select_lengths` is unchanged. When the amendment is written, it becomes a new accepted record and
this one is superseded; until then the unanimity behaviour is documented here so a reader of the
GSE144140 result knows why its primary length set looks the way it does.
