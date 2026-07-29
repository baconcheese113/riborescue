/** The build history, written for a reader who has not seen the code. Each entry is one decision or
 * one failure, the number it moved, and the record that holds it. Authored rather than exported:
 * no table knows why a threshold stayed where it was when it cost us the answer we wanted. */
export type Milestone = {
  phase: string;
  kind: "built" | "broke" | "decided";
  title: string;
  body: string;
  figure?: string;
  figureLabel?: string;
  record?: string;
};

export const TIMELINE: Milestone[] = [
  {
    phase: "Setting up",
    kind: "decided",
    title: "One lockfile for four languages",
    figure: "4", figureLabel: "languages in one lockfile",
    body: "Four languages and a dozen tools in one lockfile, so a stranger reproduces a number rather than approximating it.",
    record: "ADR-0001",
  },
  {
    phase: "Setting up",
    kind: "decided",
    title: "Nothing large goes into version control",
    figure: "0 MB", figureLabel: "large files committed",
    body: "Large files are refetched from checksummed recipes, never committed. Reproducibility rests on rebuilding, not storing.",
    record: "ADR-0002",
  },
  {
    phase: "Reading the data",
    kind: "broke",
    title: "The published reporter window was wrong in both places",
    body: "The preprint said 150, the journal said 147, the data said neither. An off-by-three shifts every feature silently.",
    figure: "72 nt",
    figureLabel: "upstream context, measured not quoted",
    record: "ADR-0004",
  },
  {
    phase: "Reading the data",
    kind: "broke",
    title: "The adapter the archive recorded was not the adapter used",
    body: "Uncaught, every read carries seventeen stray bases into alignment, nothing errors, and the rhythm the project rests on vanishes.",
    figure: "0",
    figureLabel: "reads matching the recorded adapter",
  },
  {
    phase: "Reading the data",
    kind: "broke",
    title: "Most of a footprint library was not messenger RNA",
    body: "Annotation-based filtering caught a sixth of it. The obvious implementation looks reasonable and removes a fraction of its target.",
    figure: "62.6%",
    figureLabel: "structural RNA in one library",
  },
  {
    phase: "Measuring",
    kind: "built",
    title: "Ribosome positions, resolved to the codon",
    body: "Reading position is inferred from each fragment's length and ends, turning reads into a map of where ribosomes sat.",
    figure: "24.8M",
    figureLabel: "ribosome positions mapped",
    record: "ADR-0008",
  },
  {
    phase: "Measuring",
    kind: "decided",
    title: "The pass rule was written before the answer was known",
    figure: "3", figureLabel: "conditions fixed before the data",
    body: "Density past a stop rises for reasons that are not readthrough, so three conditions were fixed in writing beforehand.",
    record: "ADR-0009 · ADR-0010",
  },
  {
    phase: "Measuring",
    kind: "built",
    title: "The assay tells a readthrough drug from a stalling one",
    body: "SRI-37240 failed in the direction its stalling mechanism predicts. An assay that cannot separate the two cannot report a null.",
    figure: "3 of 3",
    figureLabel: "conditions met by G418; 0 of 3 by SRI-37240",
    record: "ADR-0010",
  },
  {
    phase: "Measuring",
    kind: "broke",
    title: "The drug was also treating a contamination",
    figure: "1", figureLabel: "confounded cell line",
    body: "G418 is an antibacterial, and one cell line was infected. The drug under study was also clearing the contamination.",
  },
  {
    phase: "Measuring",
    kind: "broke",
    title: "A safeguard that never ran, for the whole project",
    body: "It searched for a feature name the annotation does not use, so it excluded nothing. Correcting it changed no conclusion.",
    figure: "6,724",
    figureLabel: "transcripts the fixed filter removes",
  },
  {
    phase: "Predicting",
    kind: "built",
    title: "A published model, reproduced and checked against an independent implementation",
    body: "Checked against an independent R implementation over sixty rounds. Reproducing someone exactly is the cheapest way to catch a misreading.",
    figure: "<1e-8",
    figureLabel: "agreement with the independent oracle",
  },
  {
    phase: "Predicting",
    kind: "broke",
    title: "A quarter of the diseases were one non-disease",
    body: "A placeholder meaning \u201ccondition not provided\u201d looks like a real identifier. It slipped the filter and topped the list.",
    figure: "25,917",
    figureLabel: "variants under a placeholder condition",
    record: "ADR-0015",
  },
  {
    phase: "Predicting",
    kind: "decided",
    title: "Four separate conditions, never multiplied together",
    figure: "4", figureLabel: "columns, never one score",
    body: "Any one at zero makes the product zero, so they stay four columns. A single score would hide which one failed.",
    record: "ADR-0014 · ADR-0016",
  },
  {
    phase: "What we found",
    kind: "broke",
    title: "Honest uncertainty dissolved the ranking",
    body: "With intervals rather than point estimates, the ranking dissolves. Uncertainty was not decoration on the result — it was the result.",
    figure: "5 / 70,384",
    figureLabel: "variants where one drug is distinguishable",
  },
  {
    phase: "What we found",
    kind: "broke",
    title: "Our own hypothesis did not survive its decisive control",
    body: "Two shuffles supported it. The third, scrambling only within an amino acid, reproduced most of the gain — the signal is the residue.",
    figure: "p = 0.21",
    figureLabel: "the control that decided it",
    record: "ADR-0020 · ADR-0021",
  },
  {
    phase: "What we found",
    kind: "decided",
    title: "The pre-registration cost us the answer we wanted",
    figure: "30 s", figureLabel: "to widen the rule; it stayed",
    body: "Three positive differences and a published result agreeing, with one interval in the way. Widening the rule would have taken seconds.",
    record: "ADR-0022",
  },
  {
    phase: "What we found",
    kind: "broke",
    title: "One dataset could not answer the question it was collected for",
    body: "It missed our quality floor by 0.015%. The threshold did not move — and depth was never the constraint; replication was.",
    figure: "4",
    figureLabel: "libraries per arm actually needed",
    record: "ADR-0019",
  },
  {
    phase: "Where it leads",
    kind: "built",
    title: "A third route, decidable from sequence alone",
    body: "Can an editor be placed on this stop? Just under a third can, and relaxing the targeting rule reaches 89% — the limit is placement.",
    figure: "31.2%",
    figureLabel: "reachable under the declared editor panel",
    record: "ADR-0023",
  },
];

export const PHASES = [...new Set(TIMELINE.map((m) => m.phase))];
