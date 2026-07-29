import {
  Activity,
  ArrowRight,
  CircleDot,
  Dna,
  FlaskConical,
  Gauge,
  Microscope,
  Pill,
  ScanLine,
  ShieldCheck,
  TestTube2,
  Waves,
  type LucideIcon,
} from "lucide-react";
import { experimentFrontier, scopeGap } from "../../components/charts";
import { Chart, Figure, ScopeStrip, Section } from "../../components/ui";
import { RESEARCH, count } from "../../lib/data";
import type { Experiment } from "../types";

export const metadata = {
  title: "What would settle it · RiboRescue",
  description:
    "Five experiments that would resolve the questions this project could not answer, each with what it would measure directly and what it would only imply.",
};

const GRADE_WORD: Record<string, string> = {
  none: "No measurement exists",
  indirect: "Indirect evidence only",
  partial: "Partly answered already",
};

const PROGRAMME_FLOW: Record<string, [LucideIcon, string][]> = {
  "SUPTRNA-REPLICATED": [
    [Dna, "Suppressor tRNA"],
    [TestTube2, "Patient-stop reporter"],
    [ScanLine, "Ribosome profiling"],
    [Activity, "Readthrough signature"],
  ],
  "UPSTREAM-RESIDUE": [
    [CircleDot, "Codon–residue pairs"],
    [FlaskConical, "Synthetic reporters"],
    [Microscope, "Matched library"],
    [Waves, "Context effect"],
  ],
  "G418-INDEPENDENT": [
    [Pill, "G418 + stall control"],
    [TestTube2, "Independent cell line"],
    [ScanLine, "Ribosome profiling"],
    [ShieldCheck, "Signature repeats"],
  ],
  "SUPTRNA-SAFETY": [
    [Dna, "Suppressor tRNA"],
    [CircleDot, "Normal stop codons"],
    [ScanLine, "Ribosome profiling"],
    [Gauge, "Burden map"],
  ],
  "FUNCTION-RESCUE": [
    [Pill, "Dose response"],
    [TestTube2, "Patient-relevant cells"],
    [Microscope, "Protein activity"],
    [Activity, "Function restored"],
  ],
};

function Programme({ experiment, index }: { experiment: Experiment; index: number }) {
  const e = experiment;
  const flow = PROGRAMME_FLOW[e.experiment_id] ?? [
    [FlaskConical, "Intervention"],
    [TestTube2, "Model system"],
    [Microscope, "Assay"],
    [Gauge, "Decision"],
  ];
  return (
    <article className={e.on_frontier ? "programme" : "programme dominated"}>
      <header className="programme-head">
        <div>
          <p className="route-kicker">
            {String(index + 1).padStart(2, "0")} · {e.experiment_id}
          </p>
          <h3>{e.question}</h3>
        </div>
        <span className={`tag ${e.on_frontier ? "measured" : "absent"}`}>
          {e.on_frontier ? "start here" : "beaten by another"}
        </span>
      </header>

      <div className="programme-flow" aria-label={`Experimental flow for ${e.experiment_id}`}>
        {flow.map(([Icon, label], flowIndex) => (
          <div className="programme-flow-part" key={label}>
            <div className="programme-node">
              <Icon aria-hidden="true" />
              <span>{label}</span>
            </div>
            {flowIndex < flow.length - 1 && <ArrowRight className="programme-arrow" aria-hidden="true" />}
          </div>
        ))}
      </div>

      {/* The two numbers that must never be added together: one is measured, one is a hope. */}
      <div className="scope-split">
        <div>
          <p className="scope-label measured">Measured directly</p>
          <p className="scope-figure">{e.direct_scope}</p>
        </div>
        <div>
          <p className="scope-label absent">Reached only if it generalises</p>
          <p className="scope-figure">
            <b>{count(e.potential_variants)}</b> variants · <b>{count(e.potential_genes)}</b> genes ·{" "}
            <b>{count(e.potential_conditions)}</b> conditions
          </p>
        </div>
      </div>

      <details className="programme-detail">
        <summary>How the replicate count was derived, and the rest of the design</summary>
        <dl className="detail-grid">
          {[
            ["What the lab does", e.what_the_lab_does],
            ["Compared against", e.comparison],
            ["What counts as success", e.success_criterion],
            ["If it fails", e.if_it_fails],
            ["Why it matters", e.why_it_matters],
            ["What generalisation it needs", e.generalisation_required],
            ["Why the evidence is missing", e.evidence_gap_reason],
            ["Assay", e.assay],
            ["Model system", e.model_system],
            ["Endpoint", e.endpoint],
            ["Decision rule", e.decision_rule],
            ["Replicates", e.replicates],
            ["Replicate endpoint", e.replicate_endpoint],
            ["Effect size assumed", e.replicate_effect],
            ["Variance source", e.replicate_variance_source],
            ["Design", e.replicate_design],
            ["Alpha and power", e.replicate_alpha_power],
            ["Method", e.replicate_method],
            ["Records", e.provenance],
          ].map(([label, value]) => (
            <div key={label}>
              <dt>{label}</dt>
              <dd>{value}</dd>
            </div>
          ))}
        </dl>
      </details>

      <p className="programme-foot">
        {GRADE_WORD[e.evidence_grade] ?? e.evidence_grade} · {e.complexity} to run
        {e.safety_relevant ? " · bears on safety" : ""}
        {e.resolves !== "none" ? ` · would settle ${e.resolves}` : " · opens a new question"}
      </p>
    </article>
  );
}

export default function NextPage() {
  const experiments = RESEARCH.experiments;
  const frontier = experiments.filter((e) => e.on_frontier);
  const noEvidence = experiments.filter((e) => e.evidence_grade === "none").length;

  return (
    <div className="page">
      <header className="display" style={{ padding: "56px 0 12px" }}>
        <p className="kicker">What would settle it</p>
        <h1>
          Three answers we did not get, and <em>the experiments that would.</em>
        </h1>
        <p className="lede">
          Two negatives and a dataset that could not answer its own question. Each open question
          carries a designed experiment, down to the replicate count.
        </p>
        <ScopeStrip
          items={[
            { value: String(experiments.length), label: "designed follow-up experiments" },
            { value: String(frontier.length), label: "nothing else beats on both axes" },
            { value: String(noEvidence), label: "with no measurement in existence today" },
          ]}
        />
      </header>

      <Section
        number="01"
        kicker="Where to start"
        title={
          <>
            Most needed against <em>most buildable.</em>
          </>
        }
        intro="Worth starting if nothing else is both more urgently missing and easier to run. Three of five clear that bar."
      >
        <Figure
          kind="rule"
          title="Three experiments sit on the frontier"
          caption="Evidence gap against feasibility, for five designed programmes"
          legend={[
            { className: "held", label: "filled — on the frontier, start here" },
            { className: "hollow", label: "hollow — another programme beats it on both" },
          ]}
          note={
            <>
              Position is a design judgement, not a measurement. The hollow points are not bad
              experiments — their answers partly arrive for free if the filled ones run first.
            </>
          }
        >
          <Chart svg={experimentFrontier(experiments)} size="wide" />
        </Figure>
      </Section>

      <Section
        flat
        number="02"
        kicker="The programmes"
        title={
          <>
            What each one measures, and <em>what it only implies.</em>
          </>
        }
        intro="Each card keeps two numbers apart: what the experiment touches, and who it would reach only if the result carries."
      >
        <Figure
          kind="rule"
          title="Every programme measures one variant and implies tens of thousands"
          caption="Direct scope against the population behind its generalisation, log scale"
          legend={[
            { className: "measured", label: "filled — actually measured" },
            { className: "hollow", label: "hollow — reached only if it generalises" },
          ]}
          note={
            <>
              The dashed line is the inferential distance, and it runs to four orders of magnitude.
              Nothing on the right is evidence.
            </>
          }
        >
          <Chart svg={scopeGap(experiments)} size="wide" />
        </Figure>

        <div className="programmes">
          {experiments.map((experiment, index) => (
            <Programme key={experiment.experiment_id} experiment={experiment} index={index} />
          ))}
        </div>
      </Section>

      <footer className="provenance">
        <span>designs recorded in experiments/programs.tsv</span>
        <span>replicate counts derived per ADR-0019</span>
        <span>research use only · not medical advice</span>
      </footer>
    </div>
  );
}
