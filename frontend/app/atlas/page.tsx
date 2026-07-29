import {
  atlasAgreement,
  atlasScatter,
  atlasScope,
  atlasTherapyCoverage,
} from "../../components/charts";
import { Chart, EvidenceKey, Figure, Finding, ScopeStrip, Section } from "../../components/ui";
import { EVIDENCE, RESEARCH, WEB, count, percent } from "../../lib/data";

export const metadata = {
  title: "Where can readthrough become a liability? · RiboRescue",
  description:
    "The measured G418 native-stop atlas, its agreement with prediction and the experiment that closes the suppressor-tRNA safety gap.",
};

export default function AtlasPage() {
  const safety = EVIDENCE.safety;

  if (!safety) {
    return (
      <div className="page">
        <header className="display">
          <p className="kicker">The native-stop safety atlas</p>
          <h1>
            The measurement is <em>not available.</em>
          </h1>
          <p className="lede">
            RiboRescue shows this as an evidence gap rather than treating missing measurements as
            zero unintended readthrough.
          </p>
        </header>
      </div>
    );
  }

  const therapies = WEB.therapies.map((id) => ({
    id,
    name: WEB.therapy_names[id] ?? id,
    measured: id === "G418",
  }));
  const disagreement =
    (safety.quadrants["predicted only"] ?? 0) + (safety.quadrants["measured only"] ?? 0);
  const safetyExperiment = RESEARCH.experiments.find(
    (experiment) => experiment.experiment_id === "SUPTRNA-SAFETY",
  );

  return (
    <div className="page">
      <header className="display">
        <div className="display-head">
          <div>
            <p className="kicker">The native-stop safety atlas</p>
            <h1>
              A rescue that keeps going can create <em>a second problem.</em>
            </h1>
          </div>
          <p className="lede">
            Readthrough should cross a premature stop—not a healthy protein&apos;s normal stop. This
            atlas measures that unwanted continuation for G418.
          </p>
        </div>
        <ScopeStrip
          items={[
            {
              value: count(safety.canonical_stops_scored),
              label: "normal stops represented by the model",
            },
            {
              value: count(safety.analysed),
              label: "stops with prediction and measurement",
            },
            { value: "1 drug · 1 cell line", label: "empirical atlas coverage" },
          ]}
        />
        <EvidenceKey />
      </header>

      <Section
        number="01"
        kicker="Measurement coverage"
        title={
          <>
            Most stops could be predicted. <em>Few could be measured.</em>
          </>
        }
        intro="The model scores one MANE reference transcript per gene. Direct comparison also needs enough ribosome reads around its normal stop."
      >
        <Figure
          kind="measured"
          title="Only 1,700 of 18,887 scored normal stops have a matched measurement"
          caption="Prediction coverage versus prediction-and-measurement coverage, on the same denominator"
          verdict="BLUE CAN TEST THE MODEL · Grey is missing empirical evidence, not evidence of zero native-stop readthrough."
          note={
            <>
              The matched subset is <strong>{percent(safety.analysed / safety.canonical_stops_scored)}</strong>{" "}
              of the scored atlas. Stops outside it were not expressed or sequenced deeply enough
              for this comparison in HEK293T; they are not assumed to behave like the matched set.
            </>
          }
        >
          <Chart svg={atlasScope(safety)} size="wide" />
        </Figure>
        <Finding lead={percent(safety.analysed / safety.canonical_stops_scored)} rest="directly checked">
          All later results use these {count(safety.analysed)} stops, in one cell line.
        </Finding>
      </Section>

      <Section
        flat
        number="02"
        kicker="Prediction versus measurement"
        title={
          <>
            The model sees a signal. <em>It cannot replace the assay.</em>
          </>
        }
        intro="Each dot is one normal stop. Right means higher prediction; up means more measured continuation after G418."
      >
        <Figure
          kind="measured"
          title="Predicted and measured native-stop burden have weak-to-moderate rank agreement"
          caption={`${count(safety.analysed)} matched genes · G418 versus control · HEK293T`}
          legend={[
            { className: "atlas-both", label: "high in model and measurement" },
            { className: "predicted", label: "high in model only" },
            { className: "measured", label: "high in measurement only" },
            { className: "atlas-neither", label: "lower in both" },
          ]}
          verdict="UP AND RIGHT = MORE UNINTENDED CONTINUATION · Lower-left means lower on these two measures, not proven safe."
          verdictKind="fail"
          note={
            <>
              “High” is descriptive: a measured increase above 5 percentage points in the share of
              ribosomes found after the stop, or the highest quarter of model predictions. These
              cuts were chosen to organise the atlas,{" "}
              <strong>not as biological or clinical safety thresholds</strong>. Hover a dot for the
              gene and exact values.
            </>
          }
        >
          <Chart svg={atlasScatter(safety)} size="wide" />
        </Figure>
        <Finding lead={safety.rho.toFixed(2)} rest="rank correlation">
          1 means identical ranking; 0 means none. The 95% interval is {safety.low.toFixed(2)}–
          {safety.high.toFixed(2)}.
        </Finding>
      </Section>

      <Section
        number="03"
        kicker="Where the model disagrees"
        title={
          <>
            Four groups account for every stop. <em>407 cross the line.</em>
          </>
        }
        intro="Agreement means both values fall on the same side of descriptive cuts. Both-high still means greater burden."
      >
        <Figure
          kind="measured"
          title="Nearly one quarter are high on only one of the two measures"
          caption={`Complete partition of the same ${count(safety.analysed)} matched normal stops`}
          legend={[
            { className: "atlas-both", label: "high in both · agreement, higher burden" },
            { className: "predicted", label: "model high only · disagreement" },
            { className: "measured", label: "measurement high only · disagreement" },
            { className: "atlas-neither", label: "lower in both · agreement, not a safety call" },
          ]}
          verdict="AGREEMENT SUPPORTS RANKING · Disagreement shows where direct measurement changes the model-only story."
          note={
            <>
              The groups use the same post-analysis cuts as the scatter. They are a compact
              diagnostic, not a test of harm. In particular, “lower in both” does{" "}
              <strong>not</strong> demonstrate normal protein termination or lack of toxicity.
            </>
          }
        >
          <Chart svg={atlasAgreement(safety)} size="wide" />
        </Figure>
        <Finding lead={percent(disagreement / safety.analysed)} rest="model–assay disagreement">
          {count(disagreement)} stops change category when measurement replaces prediction.
        </Finding>
      </Section>

      <Section
        flat
        number="04"
        kicker="Highest-value next experiment"
        title={
          <>
            G418 has an atlas. <em>Every other route is an open question.</em>
          </>
        }
        intro="Only G418 has a measured atlas; five other drugs and every suppressor tRNA remain unknown."
      >
        <Figure
          kind="absent"
          title="Five of six scored small molecules have no matched native-stop atlas"
          caption="Empirical coverage by therapy; dashed rows mean unknown, never zero"
          verdict="THE SAFETY LAYER IS INCOMPLETE · All suppressor-tRNA designs also lack a qualifying native-stop measurement."
          verdictKind="fail"
          note={
            <>
              The existing result is downstream ribosome occupancy for G418 in HEK293T cells. It is
              not protein-extension abundance, cell toxicity or a result in disease-relevant
              tissue. Each new therapy needs its own treated-versus-control atlas.
            </>
          }
        >
          <Chart svg={atlasTherapyCoverage(therapies, safety.analysed)} size="wide" />
        </Figure>
        <Finding lead="3 / arm" rest="next measurement">
          Test one suppressor tRNA against its empty delivery construct: three biological replicates
          per arm, transcriptome-wide ribosome profiling, no invented pass/fail cutoff.
        </Finding>
        {safetyExperiment && (
          <p className="figure-sub">
            Direct scope · {safetyExperiment.direct_scope} · endpoint · {safetyExperiment.endpoint}
          </p>
        )}
      </Section>

      <footer className="provenance">
        <span>{EVIDENCE.provenance.dataset}</span>
        <span>G418 · HEK293T</span>
        <span>{count(safety.analysed)} matched canonical stops</span>
        <span>analysis {EVIDENCE.provenance.commit}</span>
      </footer>
    </div>
  );
}
