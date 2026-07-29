import { Bug, CheckCircle2, Scale, type LucideIcon } from "lucide-react";
import { buildShape } from "../../components/charts";
import { Chart, Figure, ScopeStrip } from "../../components/ui";
import { PHASES, TIMELINE, type Milestone } from "../../lib/timeline";

export const metadata = {
  title: "How this was built · RiboRescue",
  description:
    "The decisions and the failures behind RiboRescue, in order: what was built, what broke, and what each one cost.",
};

const KIND_WORD: Record<Milestone["kind"], string> = {
  built: "Built",
  broke: "Went wrong",
  decided: "Decided",
};

const KIND_ICON: Record<Milestone["kind"], LucideIcon> = {
  built: CheckCircle2,
  broke: Bug,
  decided: Scale,
};

function Entry({ milestone }: { milestone: Milestone }) {
  const Icon = KIND_ICON[milestone.kind];
  return (
    <article className={`tl-entry tl-${milestone.kind}`}>
      <div className="tl-rail" aria-hidden="true">
        <span className="tl-dot">
          <Icon />
        </span>
      </div>
      <div className="tl-body">
        <p className="tl-kind">
          {KIND_WORD[milestone.kind]}
          {milestone.record ? <span className="tl-record">{milestone.record}</span> : null}
        </p>
        <div className={milestone.figure ? "tl-lead has-figure" : "tl-lead"}>
          {milestone.figure && (
            <p className="tl-figure">
              <b>{milestone.figure}</b>
              <span>{milestone.figureLabel}</span>
            </p>
          )}
          <div>
            <h3>{milestone.title}</h3>
            <p className="tl-text">{milestone.body}</p>
          </div>
        </div>
      </div>
    </article>
  );
}

export default function TimelinePage() {
  const counts = {
    broke: TIMELINE.filter((m) => m.kind === "broke").length,
    decided: TIMELINE.filter((m) => m.kind === "decided").length,
    built: TIMELINE.filter((m) => m.kind === "built").length,
  };

  return (
    <div className="page">
      <header className="display" style={{ padding: "56px 0 12px" }}>
        <p className="kicker">How this was built</p>
        <h1>
          Most of the work was <em>finding out what was wrong.</em>
        </h1>
        <p className="lede">
          Every project reports what it found. This one kept what it got wrong too — several of
          those errors would have left the numbers looking perfectly reasonable.
        </p>
        <ScopeStrip
          items={[
            { value: String(counts.broke), label: "problems found and traced to a cause" },
            { value: String(counts.decided), label: "rules fixed before the answer was known" },
            { value: String(counts.built), label: "pieces that survived their own checks" },
          ]}
        />
      </header>

      <section className="section">
        <Figure
          kind="rule"
          title="Half of the entries are things that went wrong"
          caption="Eighteen recorded entries across six phases of the build"
          legend={[
            { className: "alarm", label: "went wrong" },
            { className: "predicted", label: "decided in advance" },
            { className: "held", label: "built and survived its checks" },
          ]}
          note={
            <>
              Every orange square is an error that ran without complaint until someone went looking.
              Three of them would have left the published numbers looking entirely reasonable.
            </>
          }
        >
          <Chart svg={buildShape(TIMELINE, PHASES)} size="wide" />
        </Figure>
      </section>

      {PHASES.map((phase, index) => (
        <section
          key={phase}
          className={index % 2 === 0 ? "section flat" : "section"}
          aria-labelledby={`phase-${index}`}
        >
          <div className="section-head">
            <div>
              <span className="section-number">{String(index + 1).padStart(2, "0")}</span>
              <p className="kicker">{phase}</p>
              <h2 id={`phase-${index}`}>{PHASE_TITLE[phase] ?? phase}</h2>
            </div>
            <p className="section-intro">{PHASE_INTRO[phase]}</p>
          </div>
          <div className="tl-list">
            {TIMELINE.filter((m) => m.phase === phase).map((m) => (
              <Entry key={m.title} milestone={m} />
            ))}
          </div>
        </section>
      ))}

      <footer className="provenance">
        <span>decisions recorded in docs/decisions</span>
        <span>{TIMELINE.length} entries</span>
      </footer>
    </div>
  );
}

const PHASE_TITLE: Record<string, React.ReactNode> = {
  "Setting up": (
    <>
      First, make it <em>rebuildable.</em>
    </>
  ),
  "Reading the data": (
    <>
      Public data is not <em>ready to use.</em>
    </>
  ),
  Measuring: (
    <>
      Decide the rule <em>before you look.</em>
    </>
  ),
  Predicting: (
    <>
      Reproduce someone else, <em>exactly.</em>
    </>
  ),
  "What we found": (
    <>
      Three answers, and <em>none was the one we wanted.</em>
    </>
  ),
  "Where it leads": (
    <>
      One question <em>sequence can still answer.</em>
    </>
  ),
};

const PHASE_INTRO: Record<string, string> = {
  "Setting up":
    "Before any biology: could a stranger rebuild this and get the same numbers?",
  "Reading the data":
    "Three errors in how public data described itself. Each would have run without complaint.",
  Measuring:
    "Ribosome data can be talked into showing anything, so the pass rule was frozen first.",
  Predicting:
    "Scoring seventy thousand variants is only worth doing if the scoring is checked independently.",
  "What we found":
    "Two hypotheses came back negative and a dataset could not answer its question.",
  "Where it leads":
    "Where readthrough stays uncertain, geometry can still settle a different question.",
};
