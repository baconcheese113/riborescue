import {
  Dna,
  FlaskConical,
  GitBranch,
  Map,
  Microscope,
  Orbit,
  Pill,
  Search,
  ShieldAlert,
  type LucideIcon,
} from "lucide-react";
import Link from "next/link";
import { ScopeStrip, Section } from "../components/ui";
import { EVIDENCE, RESEARCH, count } from "../lib/data";

const DOORS: {
  href: string;
  kicker: string;
  title: string;
  blurb: string;
  icon: LucideIcon;
}[] = [
  {
    href: "/evidence",
    kicker: "The instrument",
    title: "Does the measurement work?",
    blurb: "Six checks on real ribosomes—including a failed hypothesis.",
    icon: Microscope,
  },
  {
    href: "/lookup",
    kicker: "One variant",
    title: "What is known about this variant?",
    blurb: "Four evidence layers for one stop.",
    icon: Search,
  },
  {
    href: "/variants",
    kicker: "The landscape",
    title: "Which routes are open, and for how many?",
    blurb: "Three rescue routes across 70,376 variants.",
    icon: Map,
  },
  {
    href: "/atlas",
    kicker: "The cost side",
    title: "What does a readthrough drug do to normal stops?",
    blurb: "G418 measured where readthrough is unwanted.",
    icon: ShieldAlert,
  },
  {
    href: "/next",
    kicker: "The gaps",
    title: "What would settle the open questions?",
    blurb: "Five experiments for the unanswered questions.",
    icon: FlaskConical,
  },
  {
    href: "/timeline",
    kicker: "The build",
    title: "How was this put together?",
    blurb: "Decisions, failures, and what survived.",
    icon: GitBranch,
  },
];

export default function Home() {
  const escape = RESEARCH.escape;
  const provenance = RESEARCH.provenance;
  const g418 = EVIDENCE.readthrough?.g418_vs_dmso.quantities ?? [];
  const held = g418.filter((q) => q.consistent).length;

  return (
    <div className="page">
      <header className="display" style={{ padding: "60px 0 12px" }}>
        <p className="kicker">Nonsense variants and the therapies aimed at them</p>
        <h1>
          A protein that stops <em>too soon.</em>
        </h1>
        <p className="lede">
          A nonsense variant stops a protein early. RiboRescue compares three rescue routes—and
          shows where the evidence runs out.
        </p>
        <ScopeStrip
          items={[
            { value: count(provenance.qualifying_variants), label: "pathogenic nonsense variants examined" },
            { value: "3", label: "possible rescue routes compared" },
            { value: `${held} of 3`, label: "conditions the positive control met" },
          ]}
        />
      </header>

      <Section
        number="01"
        kicker="The problem"
        title={
          <>
            A gene is a sentence. This is <em>a full stop in the middle.</em>
          </>
        }
        intro="Genes are read in three-letter codons. A nonsense variant turns an amino-acid codon into an early stop."
      >
        <div className="sentence">
          <p className="sentence-strip" aria-hidden="true">
            {["AUG", "GCU", "CAG", "GGU"].map((codon) => (
              <span className="codon" key={codon}>
                {codon}
              </span>
            ))}
            <span className="codon stop">UAG</span>
            {["ACU", "GAA", "CUG", "UAA"].map((codon) => (
              <span className="codon after" key={codon}>
                {codon}
              </span>
            ))}
          </p>
          <p className="sentence-key">
            <span className="made">protein built</span>
            <strong>premature stop</strong>
            <span className="lost">never built</span>
          </p>
          <p className="sentence-note">
            The cell may also destroy the message. A therapy therefore needs both surviving RNA and
            readthrough.
          </p>
        </div>
      </Section>

      <Section
        flat
        number="02"
        kicker="The three routes"
        title={
          <>
            Three ways out, and they are <em>not equally knowable.</em>
          </>
        }
        intro="One route is uncertain, one unvalidated, and one geometrically decidable."
      >
        <div className="routes">
          <article className="route">
            <span className="route-art predicted" aria-hidden="true">
              <Pill />
            </span>
            <p className="route-kicker">Route one</p>
            <h3>Small-molecule readthrough</h3>
            <p>
              Six drugs are predicted, but their uncertainty intervals overlap almost everywhere.
              They cannot be ranked honestly.
            </p>
            <p className="route-verdict verdict-warn">Predicted, but not separable</p>
          </article>
          <article className="route">
            <span className="route-art measured" aria-hidden="true">
              <Orbit />
            </span>
            <p className="route-kicker">Route two</p>
            <h3>Suppressor tRNA</h3>
            <p>
              A custom adapter restores one amino acid at one stop. Public data do not yet validate
              it.
            </p>
            <p className="route-verdict verdict-warn">Designable, but unvalidated</p>
          </article>
          <article className="route">
            <span className="route-art rule" aria-hidden="true">
              <Dna />
            </span>
            <p className="route-kicker">Route three</p>
            <h3>Base editing</h3>
            <p>
              DNA-letter geometry is decidable from sequence: <b>{count(escape?.reachable ?? 0)}</b>{" "}
              of {count(escape?.scoreable ?? 0)} scoreable variants are reachable.
            </p>
            <p className="route-verdict verdict-good">Decidable from sequence</p>
          </article>
        </div>
      </Section>

      <Section
        number="03"
        kicker="Where to go"
        title={
          <>
            Six questions, <em>six places to look.</em>
          </>
        }
        intro="Each page answers one question, keeps its denominator, and avoids a composite treatment score."
      >
        <div className="doors">
          {DOORS.map((door) => {
            const Icon = door.icon;
            return (
              <Link key={door.href} href={door.href} className="door">
                <span className="door-art" aria-hidden="true">
                  <Icon />
                </span>
                <p className="route-kicker">{door.kicker}</p>
                <h3>{door.title}</h3>
                <p>{door.blurb}</p>
                <span className="door-go" aria-hidden="true">
                  Open →
                </span>
              </Link>
            );
          })}
        </div>
      </Section>

      <footer className="provenance">
        <span>ClinVar {provenance.clinvar_release}</span>
        <span>analysis {provenance.commit}</span>
        <span>research use only · not medical advice</span>
      </footer>
    </div>
  );
}
