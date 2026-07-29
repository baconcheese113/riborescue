import * as d3 from "d3";
import {
  C,
  axisX,
  axisY,
  draw,
  fmtFixed,
  fmtInt,
  fmtSigned,
  gridX,
  outcome,
  tip,
  fmtPct1,
  type Svg,
} from "../lib/d3";
import type {
  CodonOccupancy,
  EscapeSummary,
  EvidenceContrast,
  Experiment,
  FrameByLength,
  FrontierStep,
  KineticsNullRow,
  LandscapeThreshold,
  ModelParity,
  PeriodicityPoint,
  ResearchAggregate,
  SafetyConcordance,
} from "../app/types";

const ARM_LABEL: Record<string, string> = {
  dmso: "DMSO (control)",
  g418: "G418",
  sri37240: "SRI-37240",
};
const ARM_COLOUR: Record<string, string> = {
  dmso: C.absent,
  g418: C.measured,
  sri37240: C.alarm,
};

function xTitle(svg: Svg, text: string, x: number, y: number) {
  svg
    .append("text")
    .attr("class", "group-label")
    .attr("x", x)
    .attr("y", y)
    .attr("text-anchor", "middle")
    .text(text);
}

function yTitle(svg: Svg, text: string, x: number, y: number) {
  svg
    .append("text")
    .attr("class", "group-label")
    .attr("transform", `translate(${x},${y}) rotate(-90)`)
    .attr("text-anchor", "middle")
    .text(text);
}

/** The three-nucleotide beat, drawn per treatment arm. Pooling arms would average away the
 * stop-codon difference, which is the one place the arms are supposed to differ. */
export function periodicity(rows: PeriodicityPoint[]): string {
  const W = 1060;
  const H = 350;
  const PAD = { top: 42, bottom: 60, left: 72, right: 16, gap: 46 };
  const panelW = (W - PAD.left - PAD.right - PAD.gap) / 2;
  const regions: { key: "start" | "stop"; label: string }[] = [
    { key: "start", label: "around the start codon" },
    { key: "stop", label: "around the stop codon" },
  ];
  const yMax = d3.max(rows, (r) => r.scaled ?? 0) ?? 1;
  const y = d3.scaleLinear().domain([0, yMax]).nice().range([H - PAD.bottom, PAD.top]);
  const arms = [...new Set(rows.map((r) => r.treatment ?? ""))];

  return draw(
    W,
    H,
    `Ribosome footprint density around start and stop codons for three treatment arms, showing a repeating three-nucleotide pattern with peak density ${(yMax * 100).toFixed(1)} percent.`,
    (svg) => {
      yTitle(svg, "normalised ribosome footprints", 16, (PAD.top + H - PAD.bottom) / 2);

      regions.forEach((region, index) => {
        const left = PAD.left + index * (panelW + PAD.gap);
        const slice = rows.filter((r) => r.region === region.key);
        const x = d3
          .scaleLinear()
          .domain(d3.extent(slice, (r) => r.distance) as [number, number])
          .range([left, left + panelW]);

        svg
          .append("text")
          .attr("class", "group-label")
          .attr("x", left)
          .attr("y", PAD.top - 14)
          .text(region.label);

        gridX(svg, x, PAD.top, H - PAD.bottom, x.ticks(5));
        axisX(svg, x, H - PAD.bottom, x.ticks(5), (v) => String(v));
        if (index === 0) {
          svg
            .append("g")
            .attr("class", "axis")
            .attr("transform", `translate(${left},0)`)
            .call(d3.axisLeft(y).ticks(5).tickFormat((v) => `${((v as number) * 100).toFixed(0)}%`))
            .call((g) => g.select(".domain").remove());
        }

        svg
          .append("line")
          .attr("x1", x(0))
          .attr("x2", x(0))
          .attr("y1", PAD.top)
          .attr("y2", H - PAD.bottom)
          .attr("stroke", C.line)
          .attr("stroke-dasharray", "3,3");

        const line = d3
          .line<PeriodicityPoint>()
          .x((d) => x(d.distance))
          .y((d) => y(d.scaled ?? 0));

        arms.forEach((arm) => {
          const series = slice.filter((r) => (r.treatment ?? "") === arm).sort((a, b) => a.distance - b.distance);
          svg
            .append("path")
            .attr("d", line(series) ?? "")
            .attr("fill", "none")
            .attr("stroke", ARM_COLOUR[arm] ?? C.absent)
            .attr("stroke-width", 1.5)
            .attr("stroke-linejoin", "round")
            .append("title")
            .text(
              `${ARM_LABEL[arm] ?? arm}: normalised ribosome footprints around the ${region.key} codon`,
            );
        });

        xTitle(
          svg,
          `distance from ${region.key} codon (nucleotides)`,
          left + panelW / 2,
          H - 12,
        );
      });

      svg
        .append("text")
        .attr("class", "good-label")
        .attr("x", W - PAD.right)
        .attr("y", 18)
        .attr("text-anchor", "end")
        .text("EXPECTED QUALITY SIGNAL: PEAKS REPEAT EVERY 3 NT");
    },
  );
}

/** Which footprint lengths were kept, and the one-in-three a random fragment would produce. The
 * selection is shown being made rather than asserted. */
export function frameByLength(rows: FrameByLength[]): string {
  const W = 1060;
  const H = 320;
  const PAD = { top: 34, bottom: 58, left: 72, right: 18 };
  const kept = rows.filter((r) => r.kept).map((r) => r.length);
  const x = d3
    .scaleBand<number>()
    .domain(rows.map((r) => r.length))
    .range([PAD.left, W - PAD.right])
    .padding(0.22);
  const y = d3.scaleLinear().domain([0, 0.9]).range([H - PAD.bottom, PAD.top]);

  return draw(
    W,
    H,
    `Share of reads in the coding frame by footprint length. Lengths ${kept.join(", ")} were kept; the rest sit near the one-third a random fragment would give.`,
    (svg) => {
      yTitle(svg, "reads in the correct coding frame", 16, (PAD.top + H - PAD.bottom) / 2);

      svg
        .append("g")
        .attr("class", "axis")
        .attr("transform", `translate(${PAD.left},0)`)
        .call(d3.axisLeft(y).ticks(5).tickFormat((v) => `${((v as number) * 100).toFixed(0)}%`))
        .call((g) => g.select(".domain").remove());

      tip(
        svg
          .selectAll("rect.bar")
          .data(rows)
          .join("rect")
          .attr("class", "bar")
          .attr("x", (d) => x(d.length) ?? 0)
          .attr("y", (d) => y(d.frame0_share ?? 0))
          .attr("width", x.bandwidth())
          .attr("height", (d) => y(0) - y(d.frame0_share ?? 0))
          .attr("rx", 2)
          .attr("fill", (d) => (d.kept ? C.held : C.absent)),
        (d) =>
          `${d.length} nt — ${((d.frame0_share ?? 0) * 100).toFixed(1)}% in frame, ${(
            (d.library_share ?? 0) * 100
          ).toFixed(1)}% of reads${d.kept ? " · kept" : ""}`,
      );

      svg
        .append("line")
        .attr("x1", PAD.left)
        .attr("x2", W - PAD.right)
        .attr("y1", y(1 / 3))
        .attr("y2", y(1 / 3))
        .attr("stroke", C.muted)
        .attr("stroke-dasharray", "4,3");

      svg
        .append("text")
        .attr("class", "label")
        .attr("x", W - PAD.right)
        .attr("y", y(1 / 3) - 7)
        .attr("text-anchor", "end")
        .text("what a random fragment would give");

      svg
        .selectAll("text.keep")
        .data(rows.filter((r) => r.kept))
        .join("text")
        .attr("class", "value-label keep")
        .attr("x", (d) => (x(d.length) ?? 0) + x.bandwidth() / 2)
        .attr("y", (d) => y(d.frame0_share ?? 0) - 8)
        .attr("text-anchor", "middle")
        .attr("fill", C.held)
        .text("kept");

      axisX(svg, d3.scaleLinear().domain([0, 1]).range([0, 1]), H - PAD.bottom, [], () => "");
      svg
        .selectAll("text.tickl")
        .data(rows)
        .join("text")
        .attr("class", "axis tickl")
        .attr("x", (d) => (x(d.length) ?? 0) + x.bandwidth() / 2)
        .attr("y", H - PAD.bottom + 16)
        .attr("text-anchor", "middle")
        .attr("fill", C.muted)
        .style("font-size", "10px")
        .text((d) => d.length);

      xTitle(svg, "footprint length (nucleotides)", (PAD.left + W - PAD.right) / 2, H - 10);
    },
  );
}

const CONTRAST: { key: string; drug: string }[] = [
  { key: "g418_vs_dmso", drug: "G418" },
  { key: "sri37240_vs_dmso", drug: "SRI-37240" },
];
const QUANTITY: {
  key: "downstream_occupancy" | "termination_occupancy" | "frame_gap";
  label: string;
  desired: string;
  positive: boolean;
}[] = [
  {
    key: "downstream_occupancy",
    label: "Ribosomes continuing past the stop",
    desired: "DESIRED → more ribosomes pass",
    positive: true,
  },
  {
    key: "termination_occupancy",
    label: "Ribosomes sitting at the stop",
    desired: "fewer ribosomes stall ← DESIRED",
    positive: false,
  },
  {
    key: "frame_gap",
    label: "Reading frame kept past the stop",
    desired: "DESIRED → more stay in frame",
    positive: true,
  },
];

/** Each condition's effect with its interval, filled when the condition held and hollow when it
 * did not. Three panels because the three quantities are not on one scale. */
export function discriminationEffects(contrasts: Record<string, EvidenceContrast>): string {
  const W = 1060;
  const PANEL = 96;
  const H = QUANTITY.length * PANEL + 30;
  const PAD = { left: 250, right: 104 };

  return draw(
    W,
    H,
    "Mean difference and 95 percent interval for three readthrough conditions under G418 and under SRI-37240. G418 meets all three; SRI-37240 meets none.",
    (svg) => {
      QUANTITY.forEach((quantity, panelIndex) => {
        const top = panelIndex * PANEL + 26;
        const rows = CONTRAST.map(({ key, drug }) => ({
          drug,
          q: contrasts[key]?.quantities.find((x) => x.quantity === quantity.key),
        })).filter((r) => r.q);
        const lows = rows.map((r) => r.q!.ci_low ?? 0);
        const highs = rows.map((r) => r.q!.ci_high ?? 0);
        const span = d3.extent([...lows, ...highs, 0]) as [number, number];
        const pad = (span[1] - span[0]) * 0.12 || 0.01;
        const x = d3
          .scaleLinear()
          .domain([span[0] - pad, span[1] + pad])
          .range([PAD.left, W - PAD.right]);
        const desiredX = quantity.positive ? x(0) : PAD.left;
        const desiredWidth = quantity.positive ? W - PAD.right - x(0) : x(0) - PAD.left;

        svg
          .append("rect")
          .attr("x", desiredX)
          .attr("y", top - 2)
          .attr("width", Math.max(0, desiredWidth))
          .attr("height", 54)
          .attr("fill", C.heldSoft)
          .attr("opacity", 0.62);

        svg
          .append("text")
          .attr("class", "group-label")
          .attr("x", 8)
          .attr("y", top - 8)
          .text(quantity.label);

        svg
          .append("text")
          .attr("class", "good-label")
          .attr("x", quantity.positive ? W - PAD.right : PAD.left)
          .attr("y", top - 8)
          .attr("text-anchor", quantity.positive ? "end" : "start")
          .text(quantity.desired);

        svg
          .append("line")
          .attr("x1", x(0))
          .attr("x2", x(0))
          .attr("y1", top - 2)
          .attr("y2", top + 52)
          .attr("stroke", C.muted)
          .attr("stroke-dasharray", "3,3");

        svg
          .append("text")
          .attr("class", "group-label")
          .attr("x", x(0))
          .attr("y", top + 65)
          .attr("text-anchor", "middle")
          .text("no change");

        rows.forEach((row, index) => {
          const y = top + 14 + index * 28;
          const q = row.q!;
          const held = q.consistent;
          const colour = row.drug === "G418" ? C.measured : C.alarm;

          svg
            .append("text")
            .attr("class", "label")
            .attr("x", PAD.left - 12)
            .attr("y", y + 3)
            .attr("text-anchor", "end")
            .text(row.drug);

          svg
            .append("line")
            .attr("x1", x(q.ci_low ?? 0))
            .attr("x2", x(q.ci_high ?? 0))
            .attr("y1", y)
            .attr("y2", y)
            .attr("stroke", colour)
            .attr("stroke-width", 2)
            .attr("stroke-linecap", "round");

          const marker = outcome(held, colour);
          tip(
            svg
              .append("circle")
              .datum(q)
              .attr("cx", x(q.mean_difference ?? 0))
              .attr("cy", y)
              .attr("r", 6)
              .attr("fill", marker.fill)
              .attr("stroke", marker.stroke)
              .attr("stroke-width", marker.strokeWidth),
            (d) =>
              `${row.drug} · ${quantity.label}: ${fmtSigned(5)(d.mean_difference ?? 0)} [${fmtSigned(5)(
                d.ci_low ?? 0,
              )}, ${fmtSigned(5)(d.ci_high ?? 0)}] — ${d.consistent ? "condition held" : "condition did not hold"}`,
          );

          svg
            .append("text")
            .attr("class", "value-label")
            .attr("x", W - 8)
            .attr("y", y + 3)
            .attr("text-anchor", "end")
            .attr("fill", held ? C.ink : C.muted)
            .text(fmtSigned(4)(q.mean_difference ?? 0));

          if (!held) {
            svg
              .append("text")
              .attr("class", "annotation")
              .attr("x", x(q.mean_difference ?? 0))
              .attr("y", y - 12)
              .attr("text-anchor", "middle")
              .text("NOT MET");
          }
        });
      });
    },
  );
}

/** Every library as its own dot: separation and overlap are read off the page rather than inferred
 * from a summary statistic. */
export function discriminationLibraries(
  contrasts: Record<string, EvidenceContrast>,
  quantity: "downstream_occupancy" | "termination_occupancy" | "frame_gap",
): string {
  const W = 1060;
  const H = 198;
  const PAD = { top: 34, bottom: 60, left: 130, right: 30 };
  const points = CONTRAST.flatMap(({ key, drug }) =>
    (contrasts[key]?.libraries ?? []).map((lib) => ({
      drug,
      arm: lib.treatment === "dmso" ? "control" : "treated",
      sample: lib.sample,
      value: (lib[quantity] as number | null) ?? 0,
    })),
  );
  const y = d3
    .scaleBand<string>()
    .domain(CONTRAST.map((c) => c.drug))
    .range([PAD.top, H - PAD.bottom])
    .padding(0.42);
  const x = d3
    .scaleLinear()
    .domain(d3.extent(points, (p) => p.value) as [number, number])
    .nice()
    .range([PAD.left, W - PAD.right]);

  return draw(
    W,
    H,
    "Per-library values for both drug contrasts. The G418 arms separate completely; the SRI-37240 arms interleave.",
    (svg) => {
      gridX(svg, x, PAD.top, H - PAD.bottom, x.ticks(8));
      axisX(svg, x, H - PAD.bottom, x.ticks(8), (v) => fmtFixed(3)(v));
      axisY(svg, y, PAD.left);

      CONTRAST.forEach(({ drug }) => {
        const cy = (y(drug) ?? 0) + y.bandwidth() / 2;
        svg
          .append("line")
          .attr("x1", PAD.left)
          .attr("x2", W - PAD.right)
          .attr("y1", cy)
          .attr("y2", cy)
          .attr("stroke", C.grid)
          .attr("stroke-width", 14);

        svg
          .append("text")
          .attr("class", drug === "G418" ? "good-label" : "bad-label")
          .attr("x", W - PAD.right)
          .attr("y", cy - 14)
          .attr("text-anchor", "end")
          .text(drug === "G418" ? "SEPARATES · READTHROUGH SIGNAL" : "OVERLAPS · NO READTHROUGH SIGNAL");
      });

      tip(
        svg
          .selectAll("circle.lib")
          .data(points)
          .join("circle")
          .attr("class", "lib")
          .attr("cx", (d) => x(d.value))
          .attr("cy", (d) => (y(d.drug) ?? 0) + y.bandwidth() / 2)
          .attr("r", 6.5)
          .attr("fill", (d) => (d.arm === "control" ? C.absent : C.measured))
          .attr("stroke", C.panel)
          .attr("stroke-width", 1.5),
        (d) => `${d.sample} — ${fmtFixed(5)(d.value)} (${d.arm})`,
      );

      xTitle(
        svg,
        "ribosomes past the stop (fraction of transcript signal)",
        (PAD.left + W - PAD.right) / 2,
        H - 10,
      );
    },
  );
}

const SITE_LABEL: Record<string, string> = { a: "A site", p: "P site" };

/** Occupancy per amino acid in each ribosome site. Nothing was tuned to reproduce this ordering,
 * so recovering it checks the instrument rather than reporting a result. */
export function codonOccupancy(rows: CodonOccupancy[]): string {
  const W = 1060;
  const PAD = { top: 38, bottom: 58, left: 86, right: 66 };
  const mean = (values: number[]) => values.reduce((a, b) => a + b, 0) / (values.length || 1);
  const acids = [...new Set(rows.map((r) => r.amino_acid))]
    .map((acid) => ({
      acid,
      a: mean(rows.filter((r) => r.amino_acid === acid && r.site === "a").map((r) => r.occupancy ?? 0)),
      p: mean(rows.filter((r) => r.amino_acid === acid && r.site === "p").map((r) => r.occupancy ?? 0)),
      codons: rows.filter((r) => r.amino_acid === acid && r.site === "a").length,
    }))
    .sort((x, z) => z.a - x.a);
  const H = acids.length * 17 + PAD.top + PAD.bottom;
  const y = d3
    .scaleBand<string>()
    .domain(acids.map((a) => a.acid))
    .range([PAD.top, H - PAD.bottom])
    .padding(0.25);
  const x = d3
    .scaleLinear()
    .domain(d3.extent(acids.flatMap((a) => [a.a, a.p])) as [number, number])
    .nice()
    .range([PAD.left + 26, W - PAD.right]);

  return draw(
    W,
    H,
    `Relative ribosome occupancy per amino acid in the A and P sites, from ${acids[0].acid} at ${acids[0].a.toFixed(2)} down to ${acids[acids.length - 1].acid}.`,
    (svg) => {
      yTitle(svg, "amino acid (one-letter code)", 16, (PAD.top + H - PAD.bottom) / 2);
      gridX(svg, x, PAD.top, H - PAD.bottom, x.ticks(7));
      axisX(svg, x, H - PAD.bottom, x.ticks(7), (v) => fmtFixed(1)(v));
      axisY(svg, y, PAD.left);

      svg
        .append("line")
        .attr("x1", x(1))
        .attr("x2", x(1))
        .attr("y1", PAD.top)
        .attr("y2", H - PAD.bottom)
        .attr("stroke", C.muted)
        .attr("stroke-dasharray", "3,3");

      svg
        .append("text")
        .attr("class", "label")
        .attr("x", x(1) + 6)
        .attr("y", PAD.top + 11)
        .text("average codon");

      acids.forEach((acid) => {
        const cy = (y(acid.acid) ?? 0) + y.bandwidth() / 2;
        svg
          .append("line")
          .attr("x1", x(acid.a))
          .attr("x2", x(acid.p))
          .attr("y1", cy)
          .attr("y2", cy)
          .attr("stroke", C.line)
          .attr("stroke-width", 2);
        tip(
          svg
            .append("circle")
            .datum(acid)
            .attr("cx", x(acid.p))
            .attr("cy", cy)
            .attr("r", 4.5)
            .attr("fill", C.panel)
            .attr("stroke", C.predicted)
            .attr("stroke-width", 2),
          (d) => `${d.acid} · P site ${d.p.toFixed(2)} (${d.codons} codons)`,
        );
        tip(
          svg
            .append("circle")
            .datum(acid)
            .attr("cx", x(acid.a))
            .attr("cy", cy)
            .attr("r", 4.5)
            .attr("fill", C.measured),
          (d) => `${d.acid} · A site ${d.a.toFixed(2)} (${d.codons} codons)`,
        );
      });

      xTitle(
        svg,
        "relative occupancy (farther right = ribosome waits longer)",
        (PAD.left + W - PAD.right) / 2,
        H - 10,
      );
    },
  );
}

/** Held-out accuracy beside the ceiling the experiment's own repeats impose. The gap that matters
 * is to the tick, not to one. */
export function modelParity(rows: ModelParity[]): string {
  const W = 1060;
  const PAD = { top: 42, bottom: 58, left: 96, right: 120 };
  const sorted = [...rows].sort((a, b) => (b.r2_mean ?? 0) - (a.r2_mean ?? 0));
  const H = sorted.length * 34 + PAD.top + PAD.bottom;
  const y = d3
    .scaleBand<string>()
    .domain(sorted.map((r) => r.drug))
    .range([PAD.top, H - PAD.bottom])
    .padding(0.34);
  const x = d3.scaleLinear().domain([0, 1]).range([PAD.left, W - PAD.right]);

  return draw(
    W,
    H,
    `Held-out accuracy for six drugs against each drug's replicate-reliability ceiling, from ${sorted[0].drug} at ${(sorted[0].r2_mean ?? 0).toFixed(2)} of a ${(sorted[0].ceiling ?? 0).toFixed(2)} ceiling.`,
    (svg) => {
      svg
        .append("text")
        .attr("class", "good-label")
        .attr("x", W - PAD.right)
        .attr("y", 20)
        .attr("text-anchor", "end")
        .text("DESIRED: BAR ENDS CLOSE TO BLACK CEILING TICK");
      gridX(svg, x, PAD.top, H - PAD.bottom, x.ticks(6));
      axisX(svg, x, H - PAD.bottom, x.ticks(6), (v) => fmtFixed(1)(v));
      axisY(svg, y, PAD.left);

      sorted.forEach((row) => {
        const top = y(row.drug) ?? 0;
        svg
          .append("rect")
          .attr("x", PAD.left)
          .attr("y", top)
          .attr("width", x(row.ceiling ?? 0) - PAD.left)
          .attr("height", y.bandwidth())
          .attr("rx", 3)
          .attr("fill", C.caveatSoft);
        tip(
          svg
            .append("rect")
            .datum(row)
            .attr("x", PAD.left)
            .attr("y", top)
            .attr("width", x(row.r2_mean ?? 0) - PAD.left)
            .attr("height", y.bandwidth())
            .attr("rx", 3)
            .attr("fill", C.predicted),
          (d) =>
            `${d.drug}: held-out r² ${(d.r2_mean ?? 0).toFixed(3)} against a ceiling of ${(d.ceiling ?? 0).toFixed(
              3,
            )}, over ${d.rounds} rounds`,
        );
        svg
          .append("line")
          .attr("x1", x(row.ceiling ?? 0))
          .attr("x2", x(row.ceiling ?? 0))
          .attr("y1", top - 3)
          .attr("y2", top + y.bandwidth() + 3)
          .attr("stroke", C.ink)
          .attr("stroke-width", 2);
        svg
          .append("text")
          .attr("class", "value-label")
          .attr("x", W - PAD.right + 8)
          .attr("y", top + y.bandwidth() / 2 + 3)
          .text(`${(row.r2_mean ?? 0).toFixed(2)} of ${(row.ceiling ?? 0).toFixed(2)}`);
      });

      xTitle(
        svg,
        "prediction accuracy, r² (higher is better; black tick = assay ceiling)",
        (PAD.left + W - PAD.right) / 2,
        H - 10,
      );
    },
  );
}

const SHUFFLE: { key: string; label: string; preserves: string }[] = [
  { key: "within_gene", label: "Shuffled inside each gene", preserves: "gene composition" },
  { key: "global", label: "Shuffled everywhere", preserves: "the distribution only" },
  { key: "context_matched", label: "Shuffled between synonymous codons", preserves: "the amino acid" },
];

/** What each shuffle preserves against how much of the gain it reproduces. The last one preserves
 * amino-acid identity and reproduces most of it, which is what turns the claim down. */
export function permutationNull(rows: KineticsNullRow[]): string {
  const W = 1060;
  const H = 290;
  const PAD = { top: 38, bottom: 58, left: 300, right: 190 };
  const best = SHUFFLE.map((family) => {
    const group = rows.filter((r) => r.shuffle === family.key);
    const top = group.reduce((a, b) => ((b.gain ?? 0) > (a.gain ?? 0) ? b : a));
    return { ...family, ...top };
  });
  const y = d3
    .scaleBand<string>()
    .domain(best.map((b) => b.label))
    .range([PAD.top, H - PAD.bottom])
    .padding(0.45);
  const x = d3
    .scaleLinear()
    .domain([
      Math.min(0, d3.min(best, (b) => (b.null_mean ?? 0) - (b.null_sd ?? 0)) ?? 0),
      (d3.max(best, (b) => Math.max(b.gain ?? 0, b.null_max ?? 0)) ?? 0) * 1.06,
    ])
    .range([PAD.left, W - PAD.right]);

  return draw(
    W,
    H,
    `Observed improvement against the permutation null for three shuffle families, with familywise p values of ${best
      .map((b) => b.p_familywise)
      .join(", ")}.`,
    (svg) => {
      gridX(svg, x, PAD.top, H - PAD.bottom, x.ticks(6));
      axisX(svg, x, H - PAD.bottom, x.ticks(6), (v) => fmtFixed(3)(v));
      axisY(svg, y, PAD.left);

      best.forEach((row) => {
        const cy = (y(row.label) ?? 0) + y.bandwidth() / 2;
        const held = (row.p_familywise ?? 1) < 0.05;

        svg
          .append("text")
          .attr("class", "label")
          .attr("x", PAD.left - 12)
          .attr("y", cy + 21)
          .attr("text-anchor", "end")
          .attr("fill", C.faint)
          .text(`preserves ${row.preserves}`);

        tip(
          svg
            .append("rect")
            .datum(row)
            .attr("x", x((row.null_mean ?? 0) - (row.null_sd ?? 0)))
            .attr("y", cy - 9)
            .attr("width", x((row.null_mean ?? 0) + (row.null_sd ?? 0)) - x((row.null_mean ?? 0) - (row.null_sd ?? 0)))
            .attr("height", 18)
            .attr("rx", 2)
            .attr("fill", C.grid),
          (d) => `chance range: mean ${fmtFixed(5)(d.null_mean ?? 0)} ± ${fmtFixed(5)(d.null_sd ?? 0)}`,
        );

        svg
          .append("line")
          .attr("x1", x(row.null_max ?? 0))
          .attr("x2", x(row.null_max ?? 0))
          .attr("y1", cy - 11)
          .attr("y2", cy + 11)
          .attr("stroke", C.muted)
          .attr("stroke-dasharray", "2,2");

        const marker = outcome(held, held ? C.held : C.alarm);
        tip(
          svg
            .append("circle")
            .datum(row)
            .attr("cx", x(row.gain ?? 0))
            .attr("cy", cy)
            .attr("r", 7)
            .attr("fill", marker.fill)
            .attr("stroke", marker.stroke)
            .attr("stroke-width", marker.strokeWidth),
          (d) => `${d.drug}: observed gain ${fmtFixed(5)(d.gain ?? 0)}, familywise p = ${d.p_familywise}`,
        );

        svg
          .append("text")
          .attr("class", held ? "good-label" : "bad-label")
          .attr("x", W - 8)
          .attr("y", cy - 3)
          .attr("text-anchor", "end")
          .text(held ? "GAIN EXCEEDS SHUFFLE" : "SHUFFLE EXPLAINS GAIN");
        svg
          .append("text")
          .attr("class", "value-label")
          .attr("x", W - 8)
          .attr("y", cy + 13)
          .attr("text-anchor", "end")
          .attr("fill", held ? C.held : C.alarm)
          .text(`p = ${row.p_familywise}`);
      });

      svg
        .append("text")
        .attr("class", "group-label")
        .attr("x", PAD.left)
        .attr("y", 16)
        .text("grey band = what chance produces · dashed = the best chance ever reached");

      xTitle(
        svg,
        "improvement in r² over sequence alone (farther right = more gain)",
        (PAD.left + W - PAD.right) / 2,
        H - 10,
      );
    },
  );
}

/** Readthrough cutoffs are shown twice: first as prediction alone, then after every biological
 * gate. The separate scales keep the final hundreds visible beside the initial tens of thousands. */
export function amenabilityThresholds(rows: LandscapeThreshold[]): string {
  const W = 1100;
  const H = 420;
  const PAD = { top: 78, bottom: 72, left: 130, right: 26, gap: 68 };
  const panelW = (W - PAD.left - PAD.right - PAD.gap) / 2;
  const y = d3
    .scaleBand<string>()
    .domain(rows.map((row) => `≥ ${fmtPct1(row.readthrough_threshold)}`))
    .range([PAD.top, H - PAD.bottom])
    .padding(0.36);
  const xPrediction = d3
    .scaleLinear()
    .domain([0, d3.max(rows, (row) => row.reaches_threshold) ?? 1])
    .nice()
    .range([PAD.left, PAD.left + panelW]);
  const secondLeft = PAD.left + panelW + PAD.gap;
  const xAll = d3
    .scaleLinear()
    .domain([0, d3.max(rows, (row) => row.all_conditions) ?? 1])
    .nice()
    .range([secondLeft, secondLeft + panelW]);

  return draw(
    W,
    H,
    `Predicted small-molecule amenability at three readthrough thresholds. At 0.5 percent, 70,269 variants clear the point estimate and 43,419 clear its lower bound; 3,472 clear every gate and 1,509 clear every gate on the lower bound.`,
    (svg) => {
      yTitle(svg, "minimum predicted readthrough", 16, (PAD.top + H - PAD.bottom) / 2);
      const panels = [
        {
          left: PAD.left,
          scale: xPrediction,
          title: "Prediction clears the threshold",
          point: (row: LandscapeThreshold) => row.reaches_threshold,
          lower: (row: LandscapeThreshold) => row.reaches_threshold_lower_bound,
        },
        {
          left: secondLeft,
          scale: xAll,
          title: "Decay + insertion tolerance also clear",
          point: (row: LandscapeThreshold) => row.all_conditions,
          lower: (row: LandscapeThreshold) => row.all_conditions_lower_bound,
        },
      ];

      svg
        .append("g")
        .attr("class", "axis")
        .attr("transform", `translate(${PAD.left},0)`)
        .call(d3.axisLeft(y).tickSize(0))
        .call((g) => g.select(".domain").remove());

      panels.forEach((panel) => {
        svg
          .append("text")
          .attr("class", "group-label")
          .attr("x", panel.left)
          .attr("y", 28)
          .text(panel.title);
        svg
          .append("text")
          .attr("class", "label")
          .attr("x", panel.left)
          .attr("y", 48)
          .text("pale = model midpoint only · dark = interval clears");

        const ticks = panel.scale.ticks(4);
        gridX(svg, panel.scale, PAD.top, H - PAD.bottom, ticks);
        axisX(svg, panel.scale, H - PAD.bottom, ticks, (value) => fmtInt(value));

        rows.forEach((row) => {
          const label = `≥ ${fmtPct1(row.readthrough_threshold)}`;
          const top = y(label) ?? 0;
          const height = y.bandwidth();
          const point = panel.point(row);
          const lower = panel.lower(row);

          tip(
            svg
              .append("rect")
              .datum(row)
              .attr("x", panel.left)
              .attr("y", top)
              .attr("width", Math.max(0, panel.scale(point) - panel.left))
              .attr("height", height)
              .attr("rx", 4)
              .attr("fill", C.predictedSoft),
            () => `${fmtInt(point)} variants clear the point estimate`,
          );
          tip(
            svg
              .append("rect")
              .datum(row)
              .attr("x", panel.left)
              .attr("y", top + height * 0.27)
              .attr("width", Math.max(lower === 0 ? 2 : 0, panel.scale(lower) - panel.left))
              .attr("height", height * 0.46)
              .attr("rx", 2)
              .attr("fill", lower === 0 ? C.alarm : C.predicted),
            () => `${fmtInt(lower)} variants clear the interval lower bound`,
          );

          svg
            .append("text")
            .attr("class", "value-label")
            .attr("x", Math.min(panel.scale(point) + 7, panel.left + panelW - 4))
            .attr("y", top + 11)
            .attr("text-anchor", panel.scale(point) > panel.left + panelW - 80 ? "end" : "start")
            .text(fmtInt(point));
          svg
            .append("text")
            .attr("class", lower === 0 ? "annotation" : "value-label")
            .attr("x", lower === 0 ? panel.left + 8 : panel.scale(lower) + 7)
            .attr("y", top + height - 3)
            .text(lower === 0 ? "0 · NONE" : fmtInt(lower));
        });

        const axisTop = panel.scale.domain()[1];
        xTitle(
          svg,
          `candidate variants (farther right = more) · axis runs 0 to ${fmtInt(axisTop)}`,
          panel.left + panelW / 2,
          H - 12,
        );
      });
    },
  );
}

/** Exact restoration, alternative restoration and no placement partition every scoreable variant.
 * The second track exposes how much of geometric reach survives the bystander filter. */
export function editingReach(summary: EscapeSummary): string {
  const W = 1100;
  const H = 330;
  const PAD = { top: 84, bottom: 54, left: 212, right: 24 };
  const x = d3.scaleLinear().domain([0, summary.scoreable]).range([PAD.left, W - PAD.right]);
  const tracks = [
    {
      label: "editor can be placed",
      exact: summary.exact,
      alternative: summary.alternative,
      remainder: summary.not_editable,
    },
    {
      label: "and no extra nearby edit",
      exact: summary.exact_bystander_free,
      alternative: summary.reachable_bystander_free - summary.exact_bystander_free,
      remainder: summary.scoreable - summary.reachable_bystander_free,
    },
  ];
  const colours = [C.held, C.heldSoft, C.absent];

  return draw(
    W,
    H,
    `${fmtInt(summary.reachable)} of ${fmtInt(summary.scoreable)} scoreable variants are geometrically reachable by the declared base-editor panel; ${fmtInt(summary.reachable_bystander_free)} are also free of predicted bystander edits.`,
    (svg) => {
      svg
        .append("text")
        .attr("class", "group-label")
        .attr("x", PAD.left)
        .attr("y", 28)
        .text(summary.panel);
      svg
        .append("text")
        .attr("class", "good-label")
        .attr("x", PAD.left)
        .attr("y", 52)
        .text("BEST OUTCOME: EXACT REFERENCE AMINO ACID RESTORED");

      tracks.forEach((track, rowIndex) => {
        const y = 88 + rowIndex * 100;
        svg
          .append("text")
          .attr("class", "label")
          .attr("x", PAD.left - 12)
          .attr("y", y + 27)
          .attr("text-anchor", "end")
          .text(track.label);

        let cursor = 0;
        [track.exact, track.alternative, track.remainder].forEach((value, index) => {
          const labels = ["EXACT", "ALTERNATIVE", "NO PLACEMENT"];
          const start = cursor;
          cursor += value;
          tip(
            svg
              .append("rect")
              .datum({ value, index })
              .attr("x", x(start))
              .attr("y", y)
              .attr("width", Math.max(1, x(cursor) - x(start)))
              .attr("height", 52)
              .attr("fill", colours[index]),
            () =>
              `${["exact restoration", "alternative restoration", "not reachable"][index]}: ${fmtInt(value)} (${fmtPct1(value / summary.scoreable)})`,
          );
          if (value / summary.scoreable > 0.08) {
            svg
              .append("text")
              .attr("class", "segment-label")
              .attr("x", (x(start) + x(cursor)) / 2)
              .attr("y", y + 14)
              .attr("text-anchor", "middle")
              .text(labels[index]);
            svg
              .append("text")
              .attr("class", "value-label")
              .attr("x", (x(start) + x(cursor)) / 2)
              .attr("y", y + 31)
              .attr("text-anchor", "middle")
              .text(fmtInt(value));
            svg
              .append("text")
              .attr("class", "segment-pct")
              .attr("x", (x(start) + x(cursor)) / 2)
              .attr("y", y + 45)
              .attr("text-anchor", "middle")
              .text(fmtPct1(value / summary.scoreable));
          }
        });
      });

      svg
        .append("text")
        .attr("class", "group-label")
        .attr("x", PAD.left)
        .attr("y", H - 12)
        .text("green = candidate placement · grey = editor cannot be positioned");
    },
  );
}

/** The expanded NMD rule adds two declared biological exceptions. The added region is exactly the
 * disagreement, so it is shown as one visible band instead of as two overlapping percentages. */
export function nmdRuleExpansion(summary: ResearchAggregate["nmd"]): string {
  const W = 1100;
  const H = 290;
  const PAD = { top: 78, bottom: 78, left: 72, right: 24 };
  const segments = [
    { label: "last-exon guideline escape", value: summary.escape_guideline, colour: C.held },
    { label: "added by full rules", value: summary.disagree, colour: C.alarm },
    {
      label: "predicted decay",
      value: summary.scoreable - summary.escape_full_rules,
      colour: C.absent,
    },
  ];
  const x = d3.scaleLinear().domain([0, summary.scoreable]).range([PAD.left, W - PAD.right]);

  return draw(
    W,
    H,
    `The last-exon guideline predicts escape for ${fmtInt(summary.escape_guideline)} variants. Adding start-proximal and long-exon rules brings this to ${fmtInt(summary.escape_full_rules)}, a disagreement of ${fmtInt(summary.disagree)} variants.`,
    (svg) => {
      svg
        .append("text")
        .attr("class", "label")
        .attr("x", PAD.left)
        .attr("y", 28)
        .text("What happens to the RNA before a therapy can act?");
      svg
        .append("text")
        .attr("class", "label")
        .attr("x", PAD.left)
        .attr("y", 49)
        .text(`${fmtInt(summary.scoreable)} scoreable variants`);

      let cursor = 0;
      segments.forEach((segment, index) => {
        const start = cursor;
        cursor += segment.value;
        tip(
          svg
            .append("rect")
            .datum(segment)
            .attr("x", x(start))
            .attr("y", PAD.top)
            .attr("width", x(cursor) - x(start))
            .attr("height", 74)
            .attr("fill", segment.colour),
          (d) => `${d.label}: ${fmtInt(d.value)} (${fmtPct1(d.value / summary.scoreable)})`,
        );
        if (segment.value / summary.scoreable > 0.08) {
          svg
            .append("text")
            .attr("class", "segment-label")
            .attr("x", (x(start) + x(cursor)) / 2)
            .attr("y", PAD.top + 20)
            .attr("text-anchor", "middle")
            .text(["STRICT ESCAPE", "RULE-SENSITIVE", "PREDICTED DECAY"][index]);
          svg
            .append("text")
            .attr("class", "value-label")
            .attr("x", (x(start) + x(cursor)) / 2)
            .attr("y", PAD.top + 42)
            .attr("text-anchor", "middle")
            .text(fmtInt(segment.value));
          svg
            .append("text")
            .attr("class", "segment-pct")
            .attr("x", (x(start) + x(cursor)) / 2)
            .attr("y", PAD.top + 60)
            .attr("text-anchor", "middle")
            .text(fmtPct1(segment.value / summary.scoreable));
        }
      });

      svg
        .append("text")
        .attr("class", "good-label")
        .attr("x", x(summary.escape_guideline / 2))
        .attr("y", PAD.top + 102)
        .attr("text-anchor", "middle")
        .text("RNA REMAINS");
      svg
        .append("text")
        .attr("class", "bad-label")
        .attr("x", x(summary.escape_guideline + summary.disagree / 2))
        .attr("y", PAD.top + 102)
        .attr("text-anchor", "middle")
        .text("NEEDS MEASUREMENT");
      svg
        .append("text")
        .attr("class", "group-label")
        .attr("x", x(summary.escape_full_rules + (summary.scoreable - summary.escape_full_rules) / 2))
        .attr("y", PAD.top + 102)
        .attr("text-anchor", "middle")
        .text("RNA REMOVED · READTHROUGH UNLIKELY");
    },
  );
}

/** Each design adds only variants the earlier designs did not cover. Bars show marginal reach;
 * the line shows the exact-restoration universe closing as designs accumulate. */
export function suppressorFrontier(rows: FrontierStep[]): string {
  const W = 1100;
  const H = 390;
  const PAD = { top: 62, bottom: 64, left: 82, right: 62 };
  const x = d3
    .scaleBand<number>()
    .domain(rows.map((row) => row.rank))
    .range([PAD.left, W - PAD.right])
    .padding(0.26);
  const yCoverage = d3.scaleLinear().domain([0, 1]).range([H - PAD.bottom, PAD.top]);
  const total = rows.at(-1)?.cumulative ?? 1;
  const center = (rank: number) => (x(rank) ?? 0) + x.bandwidth() / 2;

  return draw(
    W,
    H,
    `Suppressor tRNA exact-restoration frontier. UAG-Q covers ${fmtInt(rows[0]?.marginal ?? 0)} variants first; 19 stop-residue designs cover all ${fmtInt(rows.at(-1)?.cumulative ?? 0)} scoreable variants.`,
    (svg) => {
      yTitle(svg, "share of variants with exact restoration", 16, (PAD.top + H - PAD.bottom) / 2);
      svg
        .append("text")
        .attr("class", "group-label")
        .attr("x", PAD.left)
        .attr("y", 24)
        .text("DESIGN ID = STOP CODON – RESTORED AMINO ACID");
      svg
        .append("text")
        .attr("class", "good-label")
        .attr("x", W - PAD.right)
        .attr("y", 24)
        .attr("text-anchor", "end")
        .text("HIGHER BLUE LINE = BROADER EXACT-RESTORATION REACH");
      [0, 0.25, 0.5, 0.75, 1].forEach((tick) => {
        svg
          .append("line")
          .attr("class", "frontier-grid")
          .attr("x1", PAD.left)
          .attr("x2", W - PAD.right)
          .attr("y1", yCoverage(tick))
          .attr("y2", yCoverage(tick))
          .attr("stroke", C.grid);
      });
      svg
        .append("g")
        .attr("class", "axis")
        .attr("transform", `translate(${PAD.left},0)`)
        .call(
          d3
            .axisLeft(yCoverage)
            .tickValues([0, 0.25, 0.5, 0.75, 1])
            .tickFormat((value) => fmtPct1(value as number)),
        )
        .call((g) => g.select(".domain").remove());
      svg
        .append("g")
        .attr("class", "axis")
        .attr("transform", `translate(0,${H - PAD.bottom})`)
        .call(
          d3
            .axisBottom(x)
            .tickValues([1, 5, 10, 15, 19])
            .tickFormat((rank) => String(rank))
            .tickSizeOuter(0),
        );

      tip(
        svg
          .selectAll("rect.frontier-bar")
          .data(rows)
          .join("rect")
          .attr("class", "frontier-bar")
          .attr("x", (row) => x(row.rank) ?? 0)
          .attr("y", (row) => yCoverage(row.marginal / total))
          .attr("width", x.bandwidth())
          .attr("height", (row) => H - PAD.bottom - yCoverage(row.marginal / total))
          .attr("fill", C.held),
        (row) =>
          `${row.design_id} adds ${fmtInt(row.marginal)} variants (${fmtPct1(row.marginal / total)})`,
      );

      const line = d3
        .line<FrontierStep>()
        .x((row) => center(row.rank))
        .y((row) => yCoverage(row.cumulative_fraction))
        .curve(d3.curveMonotoneX);
      svg
        .append("path")
        .datum(rows)
        .attr("class", "frontier-line")
        .attr("d", line)
        .attr("fill", "none")
        .attr("stroke", C.measured)
        .attr("stroke-width", 3);

      tip(
        svg
          .selectAll("circle.frontier-point")
          .data(rows)
          .join("circle")
          .attr("class", "frontier-point")
          .attr("cx", (row) => center(row.rank))
          .attr("cy", (row) => yCoverage(row.cumulative_fraction))
          .attr("r", (row) => ([1, 2, 19].includes(row.rank) ? 6 : 3))
          .attr("fill", C.measured)
          .attr("stroke", C.panel)
          .attr("stroke-width", 2),
        (row) =>
          `${row.rank}. ${row.design_id}: ${fmtInt(row.cumulative)} variants (${fmtPct1(row.cumulative_fraction)})`,
      );

      rows
        .filter((row) => [1, 2, 19].includes(row.rank))
        .forEach((row) => {
          svg
            .append("text")
            .attr("class", "value-label")
            .attr("x", center(row.rank) + (row.rank === 19 ? -8 : 8))
            .attr("y", yCoverage(row.cumulative_fraction) - 12)
            .attr("text-anchor", row.rank === 19 ? "end" : "start")
            .text(`${row.design_id} · ${fmtPct1(row.cumulative_fraction)}`);
        });

      xTitle(
        svg,
        "designs added in frontier order",
        (PAD.left + W - PAD.right) / 2,
        H - 12,
      );
    },
  );
}

const ATLAS_GROUP = {
  both: { label: "High in model + measurement", colour: C.alarm },
  "predicted only": { label: "High in model only", colour: C.predicted },
  "measured only": { label: "High in measurement only", colour: C.measured },
  neither: { label: "Lower in both", colour: C.absent },
} as const;

/** The atlas has a large prediction denominator and a much smaller matched measurement
 * denominator. The bar makes that loss visible before any correlation is interpreted. */
export function atlasScope(summary: SafetyConcordance): string {
  const W = 1100;
  const H = 260;
  const PAD = { top: 72, bottom: 72, left: 66, right: 26 };
  const unmatched = summary.canonical_stops_scored - summary.analysed;
  const segments = [
    { label: "prediction + measurement", value: summary.analysed, colour: C.measured },
    { label: "prediction only", value: unmatched, colour: C.absent },
  ];
  const x = d3
    .scaleLinear()
    .domain([0, summary.canonical_stops_scored])
    .range([PAD.left, W - PAD.right]);

  return draw(
    W,
    H,
    `${fmtInt(summary.canonical_stops_scored)} canonical stop codons were scored by the model; ${fmtInt(summary.analysed)}, or ${fmtPct1(summary.analysed / summary.canonical_stops_scored)}, also had a usable G418 measurement in HEK293T cells.`,
    (svg) => {
      svg
        .append("text")
        .attr("class", "label")
        .attr("x", PAD.left)
        .attr("y", 28)
        .text("One normal stop codon per MANE transcript");
      svg
        .append("text")
        .attr("class", "bad-label")
        .attr("x", W - PAD.right)
        .attr("y", 28)
        .attr("text-anchor", "end")
        .text("GREY = NO MATCHED MEASUREMENT");

      let cursor = 0;
      segments.forEach((segment, index) => {
        const start = cursor;
        cursor += segment.value;
        tip(
          svg
            .append("rect")
            .datum(segment)
            .attr("x", x(start))
            .attr("y", PAD.top)
            .attr("width", x(cursor) - x(start))
            .attr("height", 72)
            .attr("rx", index === 0 ? 4 : 0)
            .attr("fill", segment.colour),
          (d) =>
            `${d.label}: ${fmtInt(d.value)} canonical stops (${fmtPct1(
              d.value / summary.canonical_stops_scored,
            )})`,
        );

        const center = (x(start) + x(cursor)) / 2;
        svg
          .append("text")
          .attr("class", "segment-label")
          .attr("x", center)
          .attr("y", PAD.top + 25)
          .attr("text-anchor", "middle")
          .text(index === 0 ? "MATCHED" : "NO MATCHED MEASUREMENT");
        svg
          .append("text")
          .attr("class", "value-label")
          .attr("x", center)
          .attr("y", PAD.top + 48)
          .attr("text-anchor", "middle")
          .text(fmtInt(segment.value));
        svg
          .append("text")
          .attr("class", "segment-pct")
          .attr("x", center)
          .attr("y", PAD.top + 64)
          .attr("text-anchor", "middle")
          .text(fmtPct1(segment.value / summary.canonical_stops_scored));
      });

      svg
        .append("text")
        .attr("class", "good-label")
        .attr("x", x(summary.analysed / 2))
        .attr("y", PAD.top + 100)
        .attr("text-anchor", "middle")
        .text("CAN COMPARE MODEL WITH DATA");
      svg
        .append("text")
        .attr("class", "group-label")
        .attr("x", x(summary.analysed + unmatched / 2))
        .attr("y", PAD.top + 100)
        .attr("text-anchor", "middle")
        .text("PREDICTION EXISTS · EMPIRICAL CHECK DOES NOT");
    },
  );
}

/** Prediction and measurement use distinct endpoints, so this is a rank-concordance view rather
 * than a calibration curve. The descriptive cuts reproduce the exported four-group summary. */
export function atlasScatter(summary: SafetyConcordance): string {
  const W = 1100;
  const H = 590;
  const PAD = { top: 82, bottom: 88, left: 108, right: 34 };
  const rows = summary.points.filter(
    (row): row is typeof row & { predicted: number; measured: number } =>
      row.predicted !== null && row.measured !== null,
  );
  const predictedCut = d3.quantile(
    rows.map((row) => row.predicted).sort(d3.ascending),
    0.75,
  ) ?? 0;
  const measuredCut = 0.05;
  const measuredTop = 0.15;
  const x = d3
    .scaleLinear()
    .domain([0, Math.max(0.045, d3.max(rows, (row) => row.predicted) ?? 0)])
    .range([PAD.left, W - PAD.right]);
  const y = d3
    .scaleLinear()
    .domain([-0.02, measuredTop])
    .range([H - PAD.bottom, PAD.top]);
  const xCut = x(predictedCut);
  const yCut = y(measuredCut);
  const panels = [
    { x: PAD.left, y: PAD.top, width: xCut - PAD.left, height: yCut - PAD.top, fill: C.measuredSoft },
    { x: xCut, y: PAD.top, width: W - PAD.right - xCut, height: yCut - PAD.top, fill: C.alarmSoft },
    {
      x: PAD.left,
      y: yCut,
      width: xCut - PAD.left,
      height: H - PAD.bottom - yCut,
      fill: C.panel,
    },
    {
      x: xCut,
      y: yCut,
      width: W - PAD.right - xCut,
      height: H - PAD.bottom - yCut,
      fill: C.predictedSoft,
    },
  ];

  return draw(
    W,
    H,
    `Predicted G418 readthrough compared with measured change in ribosome occupancy after the normal stop for ${fmtInt(summary.analysed)} genes. Spearman rank correlation is ${summary.rho.toFixed(2)}, with a 95 percent confidence interval from ${summary.low.toFixed(2)} to ${summary.high.toFixed(2)}.`,
    (svg) => {
      svg
        .selectAll("rect.atlas-zone")
        .data(panels)
        .join("rect")
        .attr("class", "atlas-zone")
        .attr("x", (d) => d.x)
        .attr("y", (d) => d.y)
        .attr("width", (d) => d.width)
        .attr("height", (d) => d.height)
        .attr("fill", (d) => d.fill)
        .attr("opacity", 0.48);

      gridX(svg, x, PAD.top, H - PAD.bottom, x.ticks(5));
      axisX(svg, x, H - PAD.bottom, x.ticks(5), (value) => fmtPct1(value));
      svg
        .append("g")
        .attr("class", "axis")
        .attr("transform", `translate(${PAD.left},0)`)
        .call(
          d3
            .axisLeft(y)
            .tickValues([-0.02, 0, 0.05, 0.1, 0.15])
            .tickFormat((value) => {
              const numeric = value as number;
              const sign = numeric > 0 ? "+" : numeric < 0 ? "−" : "";
              return `${sign}${Math.abs(numeric * 100).toFixed(0)} pp`;
            }),
        )
        .call((g) => g.select(".domain").remove());

      svg
        .append("line")
        .attr("x1", xCut)
        .attr("x2", xCut)
        .attr("y1", PAD.top)
        .attr("y2", H - PAD.bottom)
        .attr("stroke", C.predicted)
        .attr("stroke-width", 1.5)
        .attr("stroke-dasharray", "5,4");
      svg
        .append("line")
        .attr("x1", PAD.left)
        .attr("x2", W - PAD.right)
        .attr("y1", yCut)
        .attr("y2", yCut)
        .attr("stroke", C.measured)
        .attr("stroke-width", 1.5)
        .attr("stroke-dasharray", "5,4");

      tip(
        svg
          .selectAll("circle.atlas-point")
          .data(rows)
          .join("circle")
          .attr("class", "atlas-point")
          .attr("cx", (row) => x(row.predicted))
          .attr("cy", (row) => y(Math.min(row.measured, measuredTop)))
          .attr("r", (row) => (row.measured > measuredTop ? 4.3 : 2.7))
          .attr("fill", (row) => ATLAS_GROUP[row.group as keyof typeof ATLAS_GROUP]?.colour ?? C.absent)
          .attr("stroke", (row) => (row.measured > measuredTop ? C.ink : "none"))
          .attr("stroke-width", 1)
        .attr("opacity", 0.68),
        (row) =>
          `${row.gene}: predicted G418 readthrough ${fmtPct1(row.predicted)}; measured change after the stop ${fmtSigned(2)(row.measured * 100)} percentage points; ${ATLAS_GROUP[row.group as keyof typeof ATLAS_GROUP]?.label ?? row.group}`,
      );

      const labels = [
        {
          x: PAD.left + 12,
          y: PAD.top + 21,
          anchor: "start",
          className: "label",
          text: `MEASUREMENT HIGH ONLY · ${fmtInt(summary.quadrants["measured only"] ?? 0)}`,
        },
        {
          x: W - PAD.right - 12,
          y: PAD.top + 21,
          anchor: "end",
          className: "bad-label",
          text: `HIGH IN BOTH · ${fmtInt(summary.quadrants.both ?? 0)} · MORE BURDEN`,
        },
        {
          x: PAD.left + 12,
          y: H - PAD.bottom - 14,
          anchor: "start",
          className: "group-label",
          text: `LOWER IN BOTH · ${fmtInt(summary.quadrants.neither ?? 0)} · NOT “SAFE”`,
        },
        {
          x: W - PAD.right - 12,
          y: H - PAD.bottom - 14,
          anchor: "end",
          className: "label",
          text: `MODEL HIGH ONLY · ${fmtInt(summary.quadrants["predicted only"] ?? 0)}`,
        },
      ];
      svg
        .selectAll("text.atlas-quadrant-label")
        .data(labels)
        .join("text")
        .attr("class", (d) => `atlas-quadrant-label ${d.className}`)
        .attr("x", (d) => d.x)
        .attr("y", (d) => d.y)
        .attr("text-anchor", (d) => d.anchor)
        .text((d) => d.text);

      svg
        .append("text")
        .attr("class", "group-label")
        .attr("x", xCut + 7)
        .attr("y", H - PAD.bottom + 34)
        .text(`MODEL “HIGH” CUT = TOP QUARTILE (${fmtPct1(predictedCut)})`);
      svg
        .append("text")
        .attr("class", "group-label")
        .attr("x", PAD.left + 8)
        .attr("y", yCut - 8)
        .text("MEASUREMENT “HIGH” CUT = +5 PERCENTAGE POINTS");
      svg
        .append("text")
        .attr("class", "annotation")
        .attr("x", W - PAD.right)
        .attr("y", PAD.top - 13)
        .attr("text-anchor", "end")
        .text("8 VALUES ABOVE +15 POINTS ARE PINNED TO THE TOP EDGE");

      xTitle(
        svg,
        "model-predicted G418 readthrough (higher → more continuation)",
        (PAD.left + W - PAD.right) / 2,
        H - 16,
      );
      yTitle(
        svg,
        "measured change after stop (percentage points)",
        19,
        (PAD.top + H - PAD.bottom) / 2,
      );
    },
  );
}

/** Four non-overlapping groups preserve the denominator and distinguish agreement from biological
 * direction: both-high is agreement, but it is also the greatest observed native-stop burden. */
export function atlasAgreement(summary: SafetyConcordance): string {
  const W = 1100;
  const H = 390;
  const PAD = { top: 70, bottom: 72, left: 250, right: 32 };
  const order = ["both", "predicted only", "measured only", "neither"] as const;
  const rows = order.map((key) => ({
    key,
    ...ATLAS_GROUP[key],
    value: summary.quadrants[key] ?? 0,
    agreement: key === "both" || key === "neither",
    short:
      key === "both"
        ? "Both high"
        : key === "predicted only"
          ? "Model high only"
          : key === "measured only"
            ? "Measurement high only"
            : "Lower in both",
  }));
  const x = d3.scaleLinear().domain([0, summary.analysed]).range([PAD.left, W - PAD.right]);
  const y = d3
    .scaleBand<string>()
    .domain(rows.map((row) => `${row.agreement ? "AGREES" : "DISAGREES"} · ${row.short}`))
    .range([PAD.top, H - PAD.bottom])
    .padding(0.28);
  const rowLabel = (row: (typeof rows)[number]) =>
    `${row.agreement ? "AGREES" : "DISAGREES"} · ${row.short}`;

  return draw(
    W,
    H,
    `Of ${fmtInt(summary.analysed)} canonical stops, ${fmtInt(summary.quadrants.both ?? 0)} were high in both prediction and measurement, ${fmtInt(summary.quadrants["predicted only"] ?? 0)} in prediction only, ${fmtInt(summary.quadrants["measured only"] ?? 0)} in measurement only, and ${fmtInt(summary.quadrants.neither ?? 0)} in neither.`,
    (svg) => {
      gridX(svg, x, PAD.top, H - PAD.bottom, [0, 425, 850, 1275, 1700]);
      axisX(svg, x, H - PAD.bottom, [0, 425, 850, 1275, 1700], fmtInt);
      axisY(svg, y, PAD.left);

      tip(
        svg
          .selectAll("rect.atlas-agreement")
          .data(rows)
          .join("rect")
          .attr("class", "atlas-agreement")
          .attr("x", x(0))
          .attr("y", (row) => y(rowLabel(row)) ?? 0)
          .attr("width", (row) => x(row.value) - x(0))
          .attr("height", y.bandwidth())
          .attr("rx", 3)
          .attr("fill", (row) => row.colour),
        (row) =>
          `${row.label}: ${fmtInt(row.value)} of ${fmtInt(summary.analysed)} (${fmtPct1(
            row.value / summary.analysed,
          )}); ${row.agreement ? "model and measurement agree on the descriptive side" : "model and measurement disagree"}`,
      );

      svg
        .selectAll("text.atlas-agreement-value")
        .data(rows)
        .join("text")
        .attr("class", "value-label atlas-agreement-value")
        .attr("x", (row) => x(row.value) + 9)
        .attr("y", (row) => (y(rowLabel(row)) ?? 0) + y.bandwidth() / 2 + 4)
        .text((row) => `${fmtInt(row.value)} · ${fmtPct1(row.value / summary.analysed)}`);

      svg
        .append("text")
        .attr("class", "good-label")
        .attr("x", PAD.left)
        .attr("y", 27)
        .text("AGREEMENT SUPPORTS THE RANKING");
      svg
        .append("text")
        .attr("class", "bad-label")
        .attr("x", W - PAD.right)
        .attr("y", 27)
        .attr("text-anchor", "end")
        .text("DISAGREEMENT REQUIRES MEASUREMENT");
      xTitle(svg, "canonical stops (same denominator: 1,700)", (PAD.left + W - PAD.right) / 2, H - 15);
    },
  );
}

/** Only one small molecule has a matched native-stop atlas. A full bar means measured coverage;
 * the dashed rows are intentionally empty rather than zero-valued. */
export function atlasTherapyCoverage(
  therapies: { id: string; name: string; measured: boolean }[],
  analysed: number,
): string {
  const W = 1100;
  const H = 430;
  const PAD = { top: 70, bottom: 52, left: 230, right: 38 };
  const y = d3
    .scaleBand<string>()
    .domain(therapies.map((therapy) => therapy.name))
    .range([PAD.top, H - PAD.bottom])
    .padding(0.26);

  return draw(
    W,
    H,
    `Matched native-stop measurement coverage across six scored small molecules. G418 has measurements for ${fmtInt(analysed)} canonical stops in HEK293T cells; the other five molecules have no matched empirical atlas.`,
    (svg) => {
      axisY(svg, y, PAD.left);
      tip(
        svg
          .selectAll("rect.atlas-therapy")
          .data(therapies)
          .join("rect")
          .attr("class", (therapy) =>
            therapy.measured ? "atlas-therapy measured" : "atlas-therapy absent",
          )
          .attr("x", PAD.left + 18)
          .attr("y", (therapy) => y(therapy.name) ?? 0)
          .attr("width", W - PAD.right - PAD.left - 18)
          .attr("height", y.bandwidth())
          .attr("rx", 4)
          .attr("fill", (therapy) => (therapy.measured ? C.measured : "none"))
          .attr("stroke", (therapy) => (therapy.measured ? C.measured : C.absent))
          .attr("stroke-width", (therapy) => (therapy.measured ? 0 : 2))
          .attr("stroke-dasharray", (therapy) => (therapy.measured ? null : "7,5")),
        (therapy) =>
          therapy.measured
            ? `${therapy.name}: measured for ${fmtInt(analysed)} canonical stops in HEK293T cells`
            : `${therapy.name}: no matched empirical native-stop atlas`,
      );

      svg
        .selectAll("text.atlas-therapy-label")
        .data(therapies)
        .join("text")
        .attr("class", (therapy) =>
          therapy.measured ? "segment-label atlas-therapy-label" : "group-label atlas-therapy-label",
        )
        .attr("x", PAD.left + 34)
        .attr("y", (therapy) => (y(therapy.name) ?? 0) + y.bandwidth() / 2 + 4)
        .text((therapy) =>
          therapy.measured
            ? `MEASURED · HEK293T · ${fmtInt(analysed)} STOPS`
            : "NO MATCHED EMPIRICAL ATLAS",
        );

      svg
        .append("text")
        .attr("class", "good-label")
        .attr("x", PAD.left + 18)
        .attr("y", 28)
        .text("FILLED = MEASURED NATIVE-STOP BURDEN");
      svg
        .append("text")
        .attr("class", "bad-label")
        .attr("x", W - PAD.right)
        .attr("y", 28)
        .attr("text-anchor", "end")
        .text("DASHED = UNKNOWN, NOT ZERO");
    },
  );
}

export { SITE_LABEL, fmtInt };

const GAP_WORD = ["", "partly answered", "indirect evidence only", "no measurement at all"];
const FEAS_WORD = ["", "hard", "moderate"];

/** Every proposed experiment placed by how little evidence exists today against how buildable it is.
 * The frontier is the set nothing else beats on both, so it is the set worth starting with. */
export function experimentFrontier(rows: Experiment[]): string {
  const W = 1060;
  const H = 400;
  const PAD = { top: 56, bottom: 66, left: 92, right: 210 };
  const x = d3.scaleLinear().domain([0.5, 3.5]).range([PAD.left, W - PAD.right]);
  const y = d3.scaleLinear().domain([0.5, 2.5]).range([H - PAD.bottom, PAD.top]);

  return draw(
    W,
    H,
    `Five proposed experiments placed by evidence gap against feasibility. ${rows.filter((r) => r.on_frontier).length} sit on the frontier: nothing else is both more needed and easier to run.`,
    (svg) => {
      svg
        .append("rect")
        .attr("x", x(2.5))
        .attr("y", PAD.top)
        .attr("width", W - PAD.right - x(2.5))
        .attr("height", y(1.5) - PAD.top)
        .attr("fill", C.heldSoft)
        .attr("opacity", 0.55);
      svg
        .append("text")
        .attr("class", "good-label")
        .attr("x", W - PAD.right - 8)
        .attr("y", PAD.top + 16)
        .attr("text-anchor", "end")
        .text("MOST NEEDED, STILL BUILDABLE");

      [1, 2, 3].forEach((tick) =>
        svg
          .append("line")
          .attr("x1", x(tick))
          .attr("x2", x(tick))
          .attr("y1", PAD.top)
          .attr("y2", H - PAD.bottom)
          .attr("stroke", C.grid),
      );
      axisX(svg, x, H - PAD.bottom, [1, 2, 3], (v) => GAP_WORD[v] ?? String(v));
      svg
        .append("g")
        .attr("class", "axis")
        .attr("transform", `translate(${PAD.left},0)`)
        .call(d3.axisLeft(y).tickValues([1, 2]).tickFormat((v) => FEAS_WORD[v as number] ?? ""))
        .call((g) => g.select(".domain").remove());

      yTitle(svg, "how buildable it is", 18, (PAD.top + H - PAD.bottom) / 2);
      xTitle(svg, "how little evidence exists today →", (PAD.left + W - PAD.right) / 2, H - 14);

      rows.forEach((row, index) => {
        const jitter = index % 2 === 0 ? -11 : 11;
        const cx = x(row.evidence_gap);
        const cy = y(row.feasibility) + jitter;
        const marker = outcome(row.on_frontier, row.on_frontier ? C.held : C.absent);
        tip(
          svg
            .append("circle")
            .datum(row)
            .attr("cx", cx)
            .attr("cy", cy)
            .attr("r", 9)
            .attr("fill", marker.fill)
            .attr("stroke", marker.stroke)
            .attr("stroke-width", marker.strokeWidth),
          (d) =>
            `${d.experiment_id}: ${GAP_WORD[d.evidence_gap]}, ${FEAS_WORD[d.feasibility]} to run${
              d.on_frontier ? " — on the frontier" : ` — beaten by ${d.dominated_by}`
            }`,
        );
        svg
          .append("text")
          .attr("class", row.on_frontier ? "value-label" : "label")
          .attr("x", cx + 15)
          .attr("y", cy + 4)
          .attr("fill", row.on_frontier ? C.ink : C.muted)
          .text(row.experiment_id);
      });
    },
  );
}

/** What each experiment touches against who it would reach if the result generalised. The gap spans
 * four orders of magnitude, which is the argument for keeping the two numbers in separate columns. */
export function scopeGap(rows: Experiment[]): string {
  const W = 1060;
  const PAD = { top: 58, bottom: 60, left: 210, right: 132 };
  const H = rows.length * 46 + PAD.top + PAD.bottom;
  const max = d3.max(rows, (r) => r.potential_variants) ?? 1;
  const x = d3.scaleLog().domain([1, max]).range([PAD.left, W - PAD.right]);
  const y = d3
    .scaleBand<string>()
    .domain(rows.map((r) => r.experiment_id))
    .range([PAD.top, H - PAD.bottom])
    .padding(0.45);

  return draw(
    W,
    H,
    `Each experiment measures one variant in one cell line directly, while the population behind its generalisation runs from ${fmtInt(d3.min(rows, (r) => r.potential_variants) ?? 0)} to ${fmtInt(max)} variants.`,
    (svg) => {
      const ticks = [1, 10, 100, 1000, 10000, 70000];
      gridX(svg, x as unknown as d3.ScaleLinear<number, number>, PAD.top, H - PAD.bottom, ticks);
      axisX(svg, x as unknown as d3.ScaleLinear<number, number>, H - PAD.bottom, ticks, (v) => fmtInt(v));
      axisY(svg, y, PAD.left);

      svg
        .append("text")
        .attr("class", "good-label")
        .attr("x", PAD.left)
        .attr("y", 24)
        .text("MEASURED DIRECTLY");
      svg
        .append("text")
        .attr("class", "group-label")
        .attr("x", W - PAD.right)
        .attr("y", 24)
        .attr("text-anchor", "end")
        .text("REACHED ONLY IF IT GENERALISES");

      rows.forEach((row) => {
        const cy = (y(row.experiment_id) ?? 0) + y.bandwidth() / 2;
        svg
          .append("line")
          .attr("x1", x(1))
          .attr("x2", x(row.potential_variants))
          .attr("y1", cy)
          .attr("y2", cy)
          .attr("stroke", C.line)
          .attr("stroke-width", 2)
          .attr("stroke-dasharray", "4,4");
        tip(
          svg
            .append("circle")
            .datum(row)
            .attr("cx", x(1))
            .attr("cy", cy)
            .attr("r", 6)
            .attr("fill", C.measured),
          () => "Measured directly: one variant, one gene, one cell line",
        );
        tip(
          svg
            .append("circle")
            .datum(row)
            .attr("cx", x(row.potential_variants))
            .attr("cy", cy)
            .attr("r", 6)
            .attr("fill", C.panel)
            .attr("stroke", C.absent)
            .attr("stroke-width", 2),
          (d) => `${fmtInt(d.potential_variants)} variants sit behind this generalisation, none of them measured`,
        );
        svg
          .append("text")
          .attr("class", "value-label")
          .attr("x", W - PAD.right + 8)
          .attr("y", cy + 4)
          .text(fmtInt(row.potential_variants));
      });

      xTitle(svg, "variants, log scale", (PAD.left + W - PAD.right) / 2, H - 14);
    },
  );
}

/** The shape of the build itself: how many entries in each phase were things going wrong rather
 * than things being built. Written as a chart because the ratio is the point of the page. */
export function buildShape(
  rows: { phase: string; kind: "built" | "broke" | "decided" }[],
  phases: string[],
): string {
  const W = 1060;
  const H = 300;
  const PAD = { top: 56, bottom: 58, left: 168, right: 24 };
  const KINDS = ["broke", "decided", "built"] as const;
  const COLOUR: Record<string, string> = { broke: C.alarm, decided: C.predicted, built: C.held };
  const max = Math.max(...phases.map((p) => rows.filter((r) => r.phase === p).length));
  const x = d3.scaleLinear().domain([0, max]).range([PAD.left, W - PAD.right]);
  const y = d3
    .scaleBand<string>()
    .domain(phases)
    .range([PAD.top, H - PAD.bottom])
    .padding(0.36);

  return draw(
    W,
    H,
    `Eighteen entries across six phases: ${rows.filter((r) => r.kind === "broke").length} things that went wrong, ${rows.filter((r) => r.kind === "decided").length} rules fixed in advance, ${rows.filter((r) => r.kind === "built").length} pieces built.`,
    (svg) => {
      svg
        .append("text")
        .attr("class", "group-label")
        .attr("x", PAD.left)
        .attr("y", 24)
        .text("one square = one entry · orange = something went wrong");

      gridX(svg, x, PAD.top, H - PAD.bottom, d3.range(0, max + 1));
      axisX(svg, x, H - PAD.bottom, d3.range(0, max + 1), (v) => String(v));
      axisY(svg, y, PAD.left);

      phases.forEach((phase) => {
        const here = rows.filter((r) => r.phase === phase);
        let cursor = 0;
        KINDS.forEach((kind) => {
          here
            .filter((r) => r.kind === kind)
            .forEach(() => {
              tip(
                svg
                  .append("rect")
                  .datum({ phase, kind })
                  .attr("x", x(cursor) + 2)
                  .attr("y", y(phase) ?? 0)
                  .attr("width", Math.max(4, x(1) - x(0) - 4))
                  .attr("height", y.bandwidth())
                  .attr("rx", 3)
                  .attr("fill", COLOUR[kind]),
                (d) => `${d.phase} — ${d.kind === "broke" ? "went wrong" : d.kind}`,
              );
              cursor += 1;
            });
        });
      });

      xTitle(svg, "entries in this phase", (PAD.left + W - PAD.right) / 2, H - 12);
    },
  );
}

/** Where the variant-condition rows went. The placeholder slice is the one that matters: it looks
 * like an ordinary identifier and would have become the most common condition in the set. */
export function denominatorFlow(m: ResearchAggregate["mapping_completeness"]): string {
  const W = 1060;
  const H = 240;
  const PAD = { top: 74, bottom: 62, left: 24, right: 24 };
  const segments = [
    { key: "mapped", label: "USABLE", value: m.mapped, colour: C.held, note: "mapped to a real condition" },
    { key: "placeholder", label: "PLACEHOLDER", value: m.placeholder, colour: C.alarm, note: "condition not provided" },
    { key: "medgen_only", label: "PARTIAL", value: m.medgen_only, colour: C.caveat, note: "identifier but no name" },
    { key: "unmapped", label: "UNMAPPED", value: m.unmapped, colour: C.absent, note: "no cross-reference at all" },
  ];
  const total = segments.reduce((sum, s) => sum + s.value, 0);
  const x = d3.scaleLinear().domain([0, total]).range([PAD.left, W - PAD.right]);

  return draw(
    W,
    H,
    `Of ${fmtInt(total)} variant-condition rows, ${fmtInt(m.mapped)} map to a real condition and ${fmtInt(m.placeholder)} carry a placeholder meaning the condition was not provided.`,
    (svg) => {
      svg
        .append("text")
        .attr("class", "label")
        .attr("x", PAD.left)
        .attr("y", 26)
        .text(`Every variant-condition row accounted for · ${fmtInt(total)} total`);
      svg
        .append("text")
        .attr("class", "bad-label")
        .attr("x", W - PAD.right)
        .attr("y", 26)
        .attr("text-anchor", "end")
        .text("ORANGE WOULD HAVE BECOME THE COMMONEST CONDITION");

      let cursor = 0;
      segments.forEach((segment) => {
        const start = cursor;
        cursor += segment.value;
        tip(
          svg
            .append("rect")
            .datum(segment)
            .attr("x", x(start))
            .attr("y", PAD.top)
            .attr("width", Math.max(1, x(cursor) - x(start)))
            .attr("height", 62)
            .attr("fill", segment.colour),
          (d) => `${d.note}: ${fmtInt(d.value)} rows (${fmtPct1(d.value / total)})`,
        );
        const mid = (x(start) + x(cursor)) / 2;
        if (segment.value / total > 0.05) {
          svg
            .append("text")
            .attr("class", "segment-label")
            .attr("x", mid)
            .attr("y", PAD.top + 24)
            .attr("text-anchor", "middle")
            .text(segment.label);
          svg
            .append("text")
            .attr("class", "value-label")
            .attr("x", mid)
            .attr("y", PAD.top + 43)
            .attr("text-anchor", "middle")
            .text(fmtInt(segment.value));
          svg
            .append("text")
            .attr("class", "segment-pct")
            .attr("x", mid)
            .attr("y", PAD.top + 57)
            .attr("text-anchor", "middle")
            .text(fmtPct1(segment.value / total));
        }
        svg
          .append("text")
          .attr("class", "group-label")
          .attr("x", mid)
          .attr("y", PAD.top + 82)
          .attr("text-anchor", "middle")
          .text(segment.value / total > 0.05 ? segment.note : "");
      });
    },
  );
}
