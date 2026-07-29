import { readFileSync } from "node:fs";
import { join } from "node:path";
import evidence from "../public/riborescue_evidence.json";
import research from "../public/riborescue_research.json";
import web from "../public/riborescue.json";
import type { EvidenceTable, LandscapeThreshold, ResearchAggregate, WebTable } from "../app/types";

/** The evidence payload is a build input, not a runtime fetch: the export is static, so a schema
 * change surfaces as a type error here rather than as a blank page in a browser. */
export const EVIDENCE = evidence as unknown as EvidenceTable;
export const RESEARCH = research as unknown as ResearchAggregate;
export const WEB = web as unknown as WebTable;

/** The threshold sweep stays in the scientific result table. The static build reads it directly so
 * the site does not acquire a second, hand-copied version of the same result. */
export const LANDSCAPE: LandscapeThreshold[] = readFileSync(
  join(process.cwd(), "..", "results", "landscape_summary.tsv"),
  "utf8",
)
  .trim()
  .split("\n")
  .slice(1)
  .map((line) => {
    const [
      readthrough_threshold,
      variants,
      escapes_decay,
      reaches_threshold,
      reaches_threshold_lower_bound,
      most_insertions_tolerable,
      all_conditions,
      all_conditions_lower_bound,
    ] = line.split("\t").map(Number);
    return {
      readthrough_threshold,
      variants,
      escapes_decay,
      reaches_threshold,
      reaches_threshold_lower_bound,
      most_insertions_tolerable,
      all_conditions,
      all_conditions_lower_bound,
    };
  });

/** A section is absent when its arm never ran, so every reader of one asks before it draws. */
export const has = <T,>(section: T | null | undefined): section is T =>
  section !== null && section !== undefined;

export const QUANTITY_LABELS: Record<string, string> = {
  downstream_occupancy: "Ribosomes past the stop",
  termination_occupancy: "Ribosomes sitting at the stop",
  frame_gap: "Reading frame kept",
};

export const CONTRAST_LABELS: Record<string, { drug: string; expectation: string }> = {
  g418_vs_dmso: { drug: "G418", expectation: "reads through" },
  sri37240_vs_dmso: { drug: "SRI-37240", expectation: "stalls" },
};

/** Treatment arms in the order they are drawn, control first. */
export const arms = (libraries: { treatment: string | null }[]): string[] => {
  const seen = [...new Set(libraries.map((l) => l.treatment ?? "unknown"))];
  return seen.sort((a, b) => (a === "dmso" ? -1 : b === "dmso" ? 1 : a.localeCompare(b)));
};

export const percent = (value: number | null | undefined, digits = 1) =>
  value === null || value === undefined ? "—" : `${(value * 100).toFixed(digits)}%`;

export const signed = (value: number | null | undefined, digits = 4) =>
  value === null || value === undefined
    ? "—"
    : `${value >= 0 ? "+" : "−"}${Math.abs(value).toFixed(digits)}`;

export const count = (value: number) => value.toLocaleString("en-US");
