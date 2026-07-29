import {
  Activity,
  Archive,
  BadgeCheck,
  BrainCircuit,
  ChartNoAxesCombined,
  ChartSpline,
  Crosshair,
  Dna,
  FileJson,
  Files,
  LayoutDashboard,
  ListChecks,
  PackageCheck,
  PackageOpen,
  PanelsTopLeft,
  ScanSearch,
  Scissors,
  Shapes,
  ShieldCheck,
  Tally5,
  Waypoints,
  Workflow,
  type LucideIcon,
} from "lucide-react";
import { denominatorFlow } from "../../components/charts";
import { Chart, EvidenceKey, Figure, ScopeStrip, Section } from "../../components/ui";
import { EVIDENCE, RESEARCH, count, has } from "../../lib/data";
import { LANES } from "../../lib/tools";

export const metadata = {
  title: "How it works, and what it will not tell you · RiboRescue",
  description:
    "The pipeline behind RiboRescue, the limits on every claim it makes, and where each number comes from.",
};

const GLOSSARY: { term: string; plain: string }[] = [
  {
    term: "Nonsense variant",
    plain:
      "A single-letter change that turns an amino-acid instruction into a stop instruction, so the protein ends early.",
  },
  {
    term: "Readthrough",
    plain:
      "Persuading the ribosome to ignore that premature stop and keep going, without changing the DNA.",
  },
  {
    term: "Nonsense-mediated decay",
    plain:
      "A quality-control system that destroys messages carrying a premature stop. If the message is gone, there is nothing left to read through.",
  },
  {
    term: "Suppressor tRNA",
    plain:
      "An engineered adapter molecule built to recognise one stop codon and deliver one specific amino acid.",
  },
  {
    term: "Base editing",
    plain:
      "Rewriting a single DNA letter in place. Whether an editor can be positioned on a given stop is decidable from sequence.",
  },
  {
    term: "Ribosome profiling",
    plain:
      "Freezing cells, digesting the RNA that is not protected by a ribosome, and sequencing what survives — a snapshot of where ribosomes were.",
  },
  {
    term: "Interval",
    plain:
      "The range a number could plausibly take. Two numbers whose intervals overlap have not been shown to differ.",
  },
  {
    term: "Scoreable",
    plain:
      "A variant that could be placed on a reference transcript and given a prediction. Variants that could not are counted, not dropped.",
  },
];

const TOOL_ICON: Record<string, LucideIcon> = {
  "SRA toolkit": Archive,
  FastQC: ScanSearch,
  cutadapt: Scissors,
  STAR: Waypoints,
  samtools: Files,
  riboWaltz: Activity,
  featureCounts: Tally5,
  MultiQC: LayoutDashboard,
  pigz: PackageOpen,
  cyvcf2: FileJson,
  "Pydantic + Pandera": ShieldCheck,
  Biopython: Dna,
  bioframe: PanelsTopLeft,
  "statsmodels + SciPy": ChartNoAxesCombined,
  aenmd: ListChecks,
  "NMDetective-AI": BrainCircuit,
  BLOSUM62: Shapes,
  "Editor panel": Crosshair,
  Nextflow: Workflow,
  Pixi: PackageCheck,
  "R + caret": BadgeCheck,
  matplotlib: ChartSpline,
};

export default function MethodsPage() {
  const caveats = RESEARCH.caveats ?? {};
  const provenance = RESEARCH.provenance;

  return (
    <div className="page">
      <header className="display" style={{ padding: "56px 0 12px" }}>
        <p className="kicker">Methods and limits</p>
        <h1>
          How it works, and <em>what it will not tell you.</em>
        </h1>
        <p className="lede">
          Public data enter the pipeline below. Its limits live here.
        </p>
        <ScopeStrip
          items={[
            { value: count(provenance.qualifying_variants), label: "qualifying ClinVar nonsense variants" },
            { value: count(provenance.scoreable_variants), label: "placed on a transcript and scored" },
            { value: provenance.clinvar_release, label: "ClinVar release" },
          ]}
        />
      </header>

      <Section
        number="01"
        kicker="The pipeline"
        title={
          <>
            Two lanes. One touches raw reads, <em>one never does.</em>
          </>
        }
        intro="Two chains that share no inputs. Keeping them apart is what lets a figure say which kind of number it shows."
      >
        <div className="lanes">
          {LANES.map((lane) => (
            <section key={lane.key} className={`lane lane-${lane.key}`}>
              <header className="lane-head">
                <h3>{lane.title}</h3>
                <span className="lane-scale">{lane.scale}</span>
              </header>
              <p className="lane-intro">{lane.intro}</p>
              <ol className="lane-steps">
                {lane.tools.map((tool) => {
                  const Icon = TOOL_ICON[tool.name] ?? Shapes;
                  return (
                    <li key={tool.name}>
                      <span className="tool-mark" aria-hidden="true">
                        <Icon />
                      </span>
                      <span className="tool-copy">
                        <b>{tool.name}</b>
                        <span>{tool.role}</span>
                      </span>
                    </li>
                  );
                })}
              </ol>
            </section>
          ))}
        </div>
      </Section>

      <Section
        flat
        number="02"
        kicker="Reading a number"
        title={
          <>
            Four kinds of number, <em>never blended.</em>
          </>
        }
        intro="Every figure carries a marker for where its numbers came from. A measurement and an estimate are never combined."
      >
        <EvidenceKey />
        <div className="note-grid">
          <div>
            <b>Predictions are not results</b>
            <span>Model estimate from reporter data—not a patient, cell line, or tissue measurement.</span>
          </div>
          <div>
            <b>Counts are variants, not people</b>
            <span>Counts describe ClinVar variants, never people or population burden.</span>
          </div>
          <div>
            <b>Absence is not evidence</b>
            <span>Missing records and unreachable placements are unknowns, not proof of absence.</span>
          </div>
        </div>
      </Section>

      <Section
        number="03"
        kicker="Denominators"
        title={
          <>
            Every row accounted for, <em>including the useless ones.</em>
          </>
        }
        intro="A count is only meaningful against the set it came from. This is where the variant-condition rows went before any coverage figure was drawn."
      >
        <Figure
          kind="rule"
          title="A quarter of the rows name no condition at all"
          caption="Variant-condition rows by how completely they map to a real condition"
          note={
            <>
              The orange slice is a placeholder the database uses for &ldquo;condition not
              provided&rdquo;. It looks like an ordinary identifier, and counting it would have made
              a label meaning nothing into the most common disease in the set.
            </>
          }
        >
          <Chart svg={denominatorFlow(RESEARCH.mapping_completeness)} size="wide" />
        </Figure>
      </Section>

      <Section
        flat
        number="04"
        kicker="Stated limits"
        title={
          <>
            The caveats that <em>belong to the data.</em>
          </>
        }
        intro="Carried in the exported data itself, so a figure and its limit cannot drift apart."
      >
        <div className="caveats">
          {Object.entries(caveats).map(([key, text]) => (
            <details key={key}>
              <summary>{key.replace(/_/g, " ")}</summary>
              <p>{text}</p>
            </details>
          ))}
          {has(EVIDENCE.provenance) && (
            <details>
              <summary>confirmation scope</summary>
              <p>{EVIDENCE.provenance.scope}</p>
            </details>
          )}
        </div>
      </Section>

      <Section
        number="05"
        kicker="Words"
        title={
          <>
            The vocabulary, <em>in plain terms.</em>
          </>
        }
        intro="Defined once here rather than glossed differently on each page."
      >
        <div className="glossary">
          {GLOSSARY.map((entry) => (
            <details key={entry.term}>
              <summary>{entry.term}</summary>
              <p>{entry.plain}</p>
            </details>
          ))}
        </div>
      </Section>

      <footer className="provenance">
        <span>ClinVar {provenance.clinvar_release}</span>
        <span>contracts {EVIDENCE.contracts_version}</span>
        <span>analysis {provenance.commit}</span>
        <span>research use only · not medical advice</span>
      </footer>
    </div>
  );
}
