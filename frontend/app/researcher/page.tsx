"use client";

import { useEffect, useMemo, useState } from "react";
import type { Experiment, FrontierStep, ResearchAggregate } from "../types";

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


// One programme. The plain question leads; the protocol, the statistics and the provenance sit
// behind a disclosure, because a researcher deciding what to run next reads the question first and
// the decision rule only once they care.
function ExperimentCard({ e }: { e: Experiment }) {
  return (
    <article className={`experiment${e.on_frontier ? " frontier" : ""}`}>
      <header>
        <h3>{e.question}</h3>
        {e.on_frontier ? (
          <span className="pill good">Nothing beats it on every measure</span>
        ) : (
          <span className="pill muted">Beaten by {e.dominated_by.replace(/; /g, ", ")}</span>
        )}
      </header>
      <p className="why">{e.why_it_matters}</p>

      <dl className="plain">
        <dt>What the lab does</dt><dd>{e.what_the_lab_does}</dd>
        <dt>Compared against</dt><dd>{e.comparison}</dd>
        <dt>What would count as a yes</dt><dd>{e.success_criterion}</dd>
        <dt>What an unfavourable result would mean</dt><dd>{e.if_it_fails}</dd>
        <dt>Why we cannot answer it now</dt><dd>{e.evidence_gap_reason}</dd>
      </dl>

      <div className="scope">
        <div>
          <span className="scope-label">Directly measures</span>
          <p>{e.direct_scope}</p>
        </div>
        <div>
          <span className="scope-label">
            Could reach {e.potential_variants.toLocaleString()} variants,{" "}
            {e.potential_conditions.toLocaleString()} conditions — <em>if it generalises</em>
          </span>
          <p>{e.generalisation_required}</p>
        </div>
      </div>

      <p className="counts">
        {e.open_questions_addressed > 0
          ? `Addresses ${e.open_questions_addressed} open question${e.open_questions_addressed === 1 ? "" : "s"} (${e.claims_named.replace(/; /g, ", ")})`
          : e.claims_named !== "none"
            ? `Strengthens ${e.claims_named.replace(/; /g, ", ")}, already recorded as supported`
            : "Opens a question the ledger does not yet carry"} · evidence
        today: {e.evidence_grade} · complexity: {e.complexity}
        {e.safety_relevant === "TRUE" ? " · bears on safety" : ""}
      </p>

      <details>
        <summary>
          Technical details<span className="sr-only"> for {e.experiment_id}</span>
        </summary>
        <dl className="technical">
          <dt>Assay</dt><dd>{e.assay}</dd>
          <dt>Model system</dt><dd>{e.model_system}</dd>
          <dt>Endpoint</dt><dd>{e.endpoint}</dd>
          <dt>Decision rule</dt><dd>{e.decision_rule}</dd>
          <dt>Replicates</dt><dd>{e.replicates}</dd>
          <dt>Power endpoint</dt><dd>{e.replicate_endpoint}</dd>
          <dt>Effect size assumed</dt><dd>{e.replicate_effect}</dd>
          <dt>Variance source</dt><dd>{e.replicate_variance_source}</dd>
          <dt>Design</dt><dd>{e.replicate_design}</dd>
          <dt>Alpha and power</dt><dd>{e.replicate_alpha_power}</dd>
          <dt>Method</dt><dd>{e.replicate_method}</dd>
          <dt>Provenance</dt><dd>{e.provenance}</dd>
        </dl>
      </details>
    </article>
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
  const n = data.nmd;
  return (
    <div className="wrap">
      <h1>Researcher</h1>
      <p className="lede">
        Where the suppressor-tRNA designs reach across the ClinVar nonsense-variant set, and how far
        each ClinVar condition&apos;s variants are covered. Coverage is exact restoration — a design
        decodes the stop and reinserts the native residue — a prioritisation aid, <b>not a clinical
        claim</b>. A <b>ClinVar condition</b> is a disease, syndrome, phenotype, susceptibility or
        broader clinical label associated with variants in ClinVar — so these are counted as
        conditions, not verified diseases.
      </p>

      <div className="prov">
        <span>ClinVar {p.clinvar_release}</span>
        <span>commit {p.commit || "—"}</span>
        <span>{p.qualifying_variants.toLocaleString()} qualifying variants</span>
        <span>{p.scoreable_variants.toLocaleString()} scoreable</span>
        <span>{p.condition_entities.toLocaleString()} ClinVar conditions</span>
      </div>

      <section className="panel">
        <h2>What should researchers test next?</h2>
        <p className="panel-note">
          This project&rsquo;s two research claims came back negative and one dataset could not answer
          the question it was collected for. What follows from that is not a better ranking of
          therapies but the measurements that would settle the open questions. Nothing here says a
          therapy will work.
        </p>
        <div className="banner">
          <span className="tag">How to read this</span>
          <p>
            The questions and protocols below are <strong>written by hand</strong> — a lab procedure
            cannot be generated from a table. What is <strong>computed</strong> is everything that
            makes them comparable: how many variants and conditions each could reach, which recorded
            questions it addresses, and the replicate counts our own measured variance implies. There
            is deliberately no single score: trading fewer variants for a closed question is your
            decision, not ours.
          </p>
        </div>
        <div className="experiments">
          {[...data.experiments]
            .sort(
              (a, b) =>
                Number(b.on_frontier) - Number(a.on_frontier) ||
                b.open_questions_addressed - a.open_questions_addressed ||
                b.potential_variants - a.potential_variants,
            )
            .map((e) => (
              <ExperimentCard key={e.experiment_id} e={e} />
            ))}
        </div>
      </section>

      <section className="panel">
        <h2>Coverage frontiers</h2>
        <p className="panel-note">
          The fewest designs covering the most variants, genes, and — as a <b>reach</b> frontier —
          ClinVar conditions. A condition is reached when a design restores at least one of its
          variants; that is weaker than covering all of it. That {r.reachable_entities.toLocaleString()}{" "}
          reachable conditions (of {r.eligible_condition_entities.toLocaleString()} eligible, over{" "}
          {r.unique_variant_condition_pairs.toLocaleString()} unique variant–condition pairs) are all
          reached by the full panel is partly a closure property of the design universe, not a result.
        </p>
        <FrontierChart frontiers={data.frontiers} />
      </section>

      {data.escape && (
        <section className="panel">
          <h2>Escape map: base-editing reachability</h2>
          <p className="panel-note">
            A second modality, as geometry only: whether a base editor ({data.escape.panel}) can be
            placed on the premature stop — a compatible PAM at the right distance with the
            stop-removing base in the editing window. <b>Not editing efficiency, off-target activity,
            delivery, tissue access, splice consequence, or eligibility.</b> Every variant is
            accounted for: {data.escape.total.toLocaleString()} total →{" "}
            {data.escape.scoreable.toLocaleString()} placed on a stop-forming codon (
            {data.escape.unscoreable.toLocaleString()} not placeable), and the bars are fractions of
            the placed set.
          </p>
          <div className="bars">
            <div className="bar-row">
              <span className="bar-label">exact restoration</span>
              <span className="bar-track">
                <span
                  className="bar-fill"
                  style={{ width: `${(data.escape.exact / data.escape.scoreable) * 100}%` }}
                />
              </span>
              <span className="bar-num">
                {data.escape.exact.toLocaleString()} (
                {((data.escape.exact / data.escape.scoreable) * 100).toFixed(1)}%)
              </span>
            </div>
            <div className="bar-row">
              <span className="bar-label">alternative sense</span>
              <span className="bar-track">
                <span
                  className="bar-fill alt"
                  style={{ width: `${(data.escape.alternative / data.escape.scoreable) * 100}%` }}
                />
              </span>
              <span className="bar-num">
                {data.escape.alternative.toLocaleString()} (
                {((data.escape.alternative / data.escape.scoreable) * 100).toFixed(1)}%)
              </span>
            </div>
            <div className="bar-row">
              <span className="bar-label">not editable under panel</span>
              <span className="bar-track">
                <span
                  className="bar-fill warn"
                  style={{ width: `${(data.escape.not_editable / data.escape.scoreable) * 100}%` }}
                />
              </span>
              <span className="bar-num">
                {data.escape.not_editable.toLocaleString()} (
                {((data.escape.not_editable / data.escape.scoreable) * 100).toFixed(1)}%)
              </span>
            </div>
          </div>
          <p className="panel-note">
            {data.escape.reachable.toLocaleString()} variants are reachable, of which{" "}
            {data.escape.reachable_bystander_free.toLocaleString()} by a guide with no coding
            bystander, and {data.escape.exact_bystander_free.toLocaleString()} of the{" "}
            {data.escape.exact.toLocaleString()} exact restorations are bystander-free. Every
            reachable candidate uses <b>ABE7.10</b>: reverting the first base of a stop codon is a
            T→C change, which requires editing the opposite-strand A — a constraint of the genetic
            code, since no stop codon contains a cytosine for a C→T editor to act on. The categorical
            negative class is <b>not base-editable under this panel</b>, never &ldquo;no route
            exists&rdquo; — prime editing is unevaluated and readthrough is a separate continuous
            axis.
          </p>
        </section>
      )}

      <section className="panel">
        <h2>NMD escape: two rules, and where they disagree</h2>
        <p className="panel-note">
          Whether a premature stop escapes nonsense-mediated decay decides how much transcript
          survives to rescue. The 50-nt last-junction rule (the one ClinGen PVS1 uses), against the
          fuller Lindeboom rule set that adds start-proximal and long-exon escape. This split is between{" "}
          <b>two rule sets</b>, and a disagreement is where the fuller rules turn the classification —
          not evidence of which rule is correct. The model tier below reads the same variants a third
          and fourth way.
        </p>
        <div className="bars">
          <div className="bar-row">
            <span className="bar-label">guideline escape</span>
            <span className="bar-track">
              <span className="bar-fill" style={{ width: `${n.guideline_fraction * 100}%` }} />
            </span>
            <span className="bar-num">
              {n.escape_guideline.toLocaleString()} ({(n.guideline_fraction * 100).toFixed(1)}%)
            </span>
          </div>
          <div className="bar-row">
            <span className="bar-label">full-rules escape</span>
            <span className="bar-track">
              <span className="bar-fill alt" style={{ width: `${n.full_rules_fraction * 100}%` }} />
            </span>
            <span className="bar-num">
              {n.escape_full_rules.toLocaleString()} ({(n.full_rules_fraction * 100).toFixed(1)}%)
            </span>
          </div>
          <div className="bar-row">
            <span className="bar-label">disagree</span>
            <span className="bar-track">
              <span className="bar-fill warn" style={{ width: `${n.disagree_fraction * 100}%` }} />
            </span>
            <span className="bar-num">
              {n.disagree.toLocaleString()} ({(n.disagree_fraction * 100).toFixed(1)}%)
            </span>
          </div>
        </div>
        <p className="panel-note">
          Of {n.scoreable.toLocaleString()} scoreable stops the two rules split on{" "}
          {n.disagree.toLocaleString()} — driven by the long-exon rule (
          {n.driven_by_long_exon.toLocaleString()}), the start-proximal rule (
          {n.driven_by_start_proximal.toLocaleString()}), or both ({n.driven_by_both.toLocaleString()}
          ). Every split is the fuller rules escaping a stop the guideline calls decay.
        </p>
        {n.aenmd && (
          <p className="panel-note" style={{ borderTop: "1px solid var(--edge)", paddingTop: "0.75rem" }}>
            <b>Checked against aenmd</b> (Klonowski et al. 2023), the published rule tool. On the{" "}
            {n.aenmd.both_available.toLocaleString()} stops aenmd scores ({(n.aenmd.aenmd_available_fraction * 100).toFixed(0)}%
            of the set — it drops variants overlapping a splice region, and its Ensembl v105 predates
            some MANE transcripts), the hand-rolled full rule set agrees with it on{" "}
            <b>{(n.aenmd.full_rules_vs_aenmd_agree_fraction * 100).toFixed(2)}%</b>. The{" "}
            {(n.aenmd.full_escape_aenmd_decay + n.aenmd.full_decay_aenmd_escape).toLocaleString()}{" "}
            disagreements ({n.aenmd.full_decay_aenmd_escape} aenmd-only escapes,{" "}
            {n.aenmd.full_escape_aenmd_decay} rule-only) are aenmd&apos;s own edge cases, not a different
            idea of NMD.
          </p>
        )}
        {n.nmdetective && n.nmdetective.full_rules?.separation != null && (
          <p className="panel-note" style={{ borderTop: "1px solid var(--edge)", paddingTop: "0.75rem" }}>
            <b>And against NMDetective-AI</b> (Veiner et al. 2026), the deep model — a fine-tuned
            Orthrus/Mamba encoder run here on a local GPU. It gives a continuous NMD-efficiency score,
            not a verdict, so nothing thresholds it. Over the{" "}
            {n.nmdetective.both_available.toLocaleString()} stops it scored (
            {(n.nmdetective.available_fraction * 100).toFixed(0)}% — its GENCODE v26 reference predates
            more MANE transcripts than aenmd&apos;s), the stops the full rule set calls escape carry a
            markedly lower mean efficiency than those it calls decay (
            {n.nmdetective.full_rules.escape_mean_efficiency?.toFixed(2)} vs{" "}
            {n.nmdetective.full_rules.decay_mean_efficiency?.toFixed(2)}, a separation of{" "}
            <b>{n.nmdetective.full_rules.separation.toFixed(2)}</b>). Lower efficiency means more
            escape, so the deep model moves in the same direction as the rules without being told them.
          </p>
        )}
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
          The largest ClinVar conditions by eligible nonsense-variant denominator. Three metrics are
          kept apart: <b>reach</b> (any variant covered), the <b>covered fraction</b> (covered over
          eligible), and <b>complete</b> (every eligible variant covered).
        </p>
        <div className="scroller">
          <table>
            <thead>
              <tr>
                <th>ClinVar condition</th>
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
            A therapy × condition heatmap needs a therapy-level condition join, and a function
            filter needs a function layer that does not exist yet (ADR-0018 proposes one). Both are
            deliberately absent rather than mocked. The NMD layer is built and is above.
          </p>
        </div>
      </section>
    </div>
  );
}
