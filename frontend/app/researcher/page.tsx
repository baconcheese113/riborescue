"use client";

import { useEffect, useMemo, useState } from "react";
import type { FrontierStep, ResearchAggregate } from "../types";

// Validated categorical hues (blue / aqua / orange), assigned in fixed order to the three frontiers.
const SERIES: { key: "variants" | "genes" | "conditions"; label: string; color: string }[] = [
  { key: "variants", label: "Variants", color: "#2a78d6" },
  { key: "genes", label: "Genes", color: "#1baf7a" },
  { key: "conditions", label: "Conditions (reach)", color: "#eb6834" },
];

const W = 360;
const H = 200;
const PAD = { l: 36, r: 14, t: 12, b: 28 };

function path(steps: FrontierStep[], maxRank: number): string {
  const x = (rank: number) => PAD.l + (rank / maxRank) * (W - PAD.l - PAD.r);
  const y = (frac: number) => PAD.t + (1 - frac) * (H - PAD.t - PAD.b);
  // Anchor the curve at (0 designs, 0 coverage) so the shape reads from the origin.
  const points = [[0, 0] as const, ...steps.map((s) => [s.rank, s.cumulative_fraction] as const)];
  return points.map(([r, f], i) => `${i ? "L" : "M"}${x(r).toFixed(1)},${y(f).toFixed(1)}`).join(" ");
}

// The three coverage curves on one axis: cumulative fraction reached against panel size. A design
// panel of size k reaches the fraction the curve shows at k. One axis, three series, direct-labelled.
function FrontierChart({ frontiers }: { frontiers: ResearchAggregate["frontiers"] }) {
  const maxRank = Math.max(1, ...SERIES.map((s) => frontiers[s.key].length));
  const x = (rank: number) => PAD.l + (rank / maxRank) * (W - PAD.l - PAD.r);
  const y = (frac: number) => PAD.t + (1 - frac) * (H - PAD.t - PAD.b);
  return (
    <figure className="chart-figure">
      <svg viewBox={`0 0 ${W} ${H}`} role="img" aria-label="Coverage frontiers by panel size" className="chart">
        {[0, 0.25, 0.5, 0.75, 1].map((f) => (
          <g key={f}>
            <line x1={PAD.l} x2={W - PAD.r} y1={y(f)} y2={y(f)} className="grid" />
            <text x={PAD.l - 6} y={y(f) + 3} className="axis-label" textAnchor="end">
              {f * 100}%
            </text>
          </g>
        ))}
        {[0, Math.ceil(maxRank / 2), maxRank].map((r) => (
          <text key={r} x={x(r)} y={H - PAD.b + 16} className="axis-label" textAnchor="middle">
            {r}
          </text>
        ))}
        <text x={(W) / 2} y={H - 2} className="axis-title" textAnchor="middle">
          designs in panel
        </text>
        {SERIES.map((s) => (
          <path key={s.key} d={path(frontiers[s.key], maxRank)} fill="none" stroke={s.color} strokeWidth={2} />
        ))}
        {SERIES.map((s) => {
          const last = frontiers[s.key];
          if (!last.length) return null;
          const end = last[last.length - 1];
          return <circle key={s.key} cx={x(end.rank)} cy={y(end.cumulative_fraction)} r={3} fill={s.color} />;
        })}
      </svg>
      <div className="legend">
        {SERIES.map((s) => {
          const steps = frontiers[s.key];
          const at = steps.find((step) => step.rank === 5) ?? steps[steps.length - 1];
          return (
            <span className="legend-item" key={s.key}>
              <span className="swatch" style={{ background: s.color }} />
              {s.label}
              {at ? <span className="legend-note"> · {(at.cumulative_fraction * 100).toFixed(0)}% at {at.rank}</span> : null}
            </span>
          );
        })}
      </div>
    </figure>
  );
}

function CompletenessBars({ counts }: { counts: Record<string, number> }) {
  const total = Object.values(counts).reduce((a, b) => a + b, 0) || 1;
  const order = ["mapped", "medgen_only", "placeholder", "unmapped"];
  return (
    <div className="bars">
      {order
        .filter((k) => k in counts)
        .map((k) => (
          <div className="bar-row" key={k}>
            <span className="bar-label">{k.replace("_", " ")}</span>
            <span className="bar-track">
              <span className="bar-fill" style={{ width: `${(counts[k] / total) * 100}%` }} />
            </span>
            <span className="bar-num">
              {counts[k].toLocaleString()} ({((counts[k] / total) * 100).toFixed(0)}%)
            </span>
          </div>
        ))}
    </div>
  );
}

export default function ResearcherPage() {
  const [data, setData] = useState<ResearchAggregate | null>(null);
  const [error, setError] = useState<string | null>(null);

  useEffect(() => {
    fetch("/riborescue_research.json")
      .then((r) => {
        if (!r.ok) throw new Error(`riborescue_research.json returned ${r.status}`);
        return r.json();
      })
      .then(setData)
      .catch((e) => setError(String(e?.message ?? e)));
  }, []);

  const top = useMemo(() => data?.condition_coverage_top.slice(0, 25) ?? [], [data]);

  if (error)
    return (
      <div className="wrap">
        <h1>Researcher</h1>
        <div className="banner" style={{ borderLeftColor: "var(--bad)" }}>
          <span className="tag" style={{ color: "var(--bad)" }}>Could not load data</span>
          <p>
            {error}. Generate it with{" "}
            <code>riborescue export-research … --out frontend/public/riborescue_research.json</code>,
            then reload.
          </p>
        </div>
      </div>
    );
  if (!data) return <div className="wrap">Loading…</div>;

  const p = data.provenance;
  const r = data.reach_denominator;
  return (
    <div className="wrap">
      <h1>Researcher</h1>
      <p className="lede">
        Where the suppressor-tRNA designs reach across the ClinVar nonsense-variant set, and how far
        each condition&apos;s variants are covered. Coverage is exact restoration — a design decodes
        the stop and reinserts the native residue — a prioritisation aid, <b>not a clinical claim</b>.
        Entities are ClinVar <b>conditions</b> (MedGen concepts): diseases, but also findings and
        broad labels, so they are counted as condition entities, not verified diseases.
      </p>

      <div className="prov">
        <span>ClinVar {p.clinvar_release}</span>
        <span>commit {p.commit || "—"}</span>
        <span>{p.qualifying_variants.toLocaleString()} qualifying variants</span>
        <span>{p.scoreable_variants.toLocaleString()} scoreable</span>
        <span>{p.condition_entities.toLocaleString()} condition entities</span>
      </div>

      <section className="panel">
        <h2>Coverage frontiers</h2>
        <p className="panel-note">
          The fewest designs covering the most variants, genes, and — as a <b>reach</b> frontier —
          condition entities. An entity is reached when a design restores at least one of its
          variants; that is weaker than covering all of it. That {r.reachable_entities.toLocaleString()}{" "}
          reachable entities (of {r.eligible_condition_entities.toLocaleString()} eligible, over{" "}
          {r.unique_variant_condition_pairs.toLocaleString()} unique variant–condition pairs) are all
          reached by the full panel is partly a closure property of the design universe, not a result.
        </p>
        <FrontierChart frontiers={data.frontiers} />
      </section>

      <section className="panel">
        <h2>Cross-reference completeness</h2>
        <p className="panel-note">
          How many variant-condition rows carry OMIM or Orphanet beyond MedGen — and these are
          identifiers only, adding no prevalence or treatment status. Roughly a quarter of rows are
          placeholder or unmapped, a material missingness for any condition-level reading; it is
          reported, not filtered.
        </p>
        <CompletenessBars counts={data.mapping_completeness} />
      </section>

      <section className="panel">
        <h2>Conditions by eligible variants</h2>
        <p className="panel-note">
          The largest condition entities by eligible nonsense-variant denominator. Three metrics are
          kept apart: <b>reach</b> (any variant covered), the <b>covered fraction</b> (covered over
          eligible), and <b>complete</b> (every eligible variant covered).
        </p>
        <div className="scroller">
          <table>
            <thead>
              <tr>
                <th>Condition entity</th>
                <th>Covered</th>
                <th>Fraction</th>
                <th>Complete</th>
                <th>Genes</th>
                <th>Refs</th>
              </tr>
            </thead>
            <tbody>
              {top.map((d) => (
                <tr key={d.medgen}>
                  <td>{d.condition_name.replace(/_/g, " ")}</td>
                  <td className="num">
                    {d.model_covered.toLocaleString()} / {d.eligible_variants.toLocaleString()}
                  </td>
                  <td className="num">{(d.covered_fraction * 100).toFixed(1)}%</td>
                  <td>
                    <span className={`pill ${d.complete ? "good" : "muted"}`}>
                      {d.complete ? "complete" : "partial"}
                    </span>
                  </td>
                  <td className="num">{d.genes}</td>
                  <td>
                    <span className={`pill ${d.mapping_completeness === "mapped" ? "good" : "muted"}`}>
                      {d.mapping_completeness === "mapped" ? "OMIM/Orphanet" : "MedGen only"}
                    </span>
                  </td>
                </tr>
              ))}
            </tbody>
          </table>
        </div>
      </section>

      <section className="panel">
        <h2>Downloads &amp; what is not here yet</h2>
        <p className="panel-note">
          <a href="/riborescue_research.json" download>
            Download the aggregate JSON
          </a>{" "}
          — frontiers, per-condition coverage, denominators, and provenance.
        </p>
        <div className="banner unavailable">
          <span className="tag">Planned</span>
          <p>
            A therapy × condition heatmap needs a therapy-level condition join, and an NMD/function
            model filter needs those layers built. They are deliberately absent rather than mocked.
          </p>
        </div>
      </section>
    </div>
  );
}
