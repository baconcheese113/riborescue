import {
  amenabilityThresholds,
  editingReach,
  nmdRuleExpansion,
  suppressorFrontier,
} from "../../components/charts";
import { Chart, EvidenceKey, Figure, Finding, ScopeStrip, Section } from "../../components/ui";
import { LANDSCAPE, RESEARCH, count, percent } from "../../lib/data";

export const metadata = {
  title: "Which variants can be reached? · RiboRescue",
  description:
    "The nonsense-variant landscape across predicted readthrough, base editing, transcript decay and exact-restoration suppressor tRNA designs.",
};

export default function VariantsPage() {
  const { provenance, escape, nmd, frontiers } = RESEARCH;
  const frontier = frontiers.variants;
  const firstTwo = frontier[1];
  const onePercent = LANDSCAPE.find((row) => row.readthrough_threshold === 0.01);

  return (
    <div className="page">
      <header className="display">
        <div className="display-head">
          <div>
            <p className="kicker">The variant landscape</p>
            <h1>
              One premature stop. <em>Three possible routes.</em>
            </h1>
          </div>
          <p className="lede">
            Each variant gets three separate questions: drug readthrough, DNA editing, or
            exact-restoration suppressor tRNA.
          </p>
        </div>
        <ScopeStrip
          items={[
            { value: count(provenance.scoreable_variants), label: "variants with a complete score" },
            { value: count(provenance.condition_entities), label: "mapped condition entities" },
            { value: count(frontier.length), label: "exact-restoration tRNA designs" },
          ]}
        />
        <EvidenceKey />
      </header>

      <Section
        number="01"
        kicker="Small molecules"
        title={
          <>
            A broad prediction narrows when <em>uncertainty bites.</em>
          </>
        }
        intro="Point estimates screen broadly; interval lower bounds demand stronger support. Decay and residue tolerance narrow both."
      >
        <Figure
          kind="predicted"
          title="Every stronger requirement removes candidates"
          caption={`Threshold sweep across ${count(provenance.scoreable_variants)} scoreable ClinVar nonsense variants`}
          legend={[
            { className: "predicted-soft", label: "pale — model midpoint clears" },
            { className: "predicted", label: "dark — interval lower bound also clears" },
          ]}
          verdict="DARK BARS ARE STRONGER SUPPORT · More candidates means broader predicted reach, not proven benefit."
          note={
            <>
              “All gates” requires predicted readthrough, transcript escape and a majority of tested
              amino-acid substitutions to be tolerated. It is a <strong>candidate screen</strong>,
              not evidence that a compound works in a patient.
            </>
          }
        >
          <Chart svg={amenabilityThresholds(LANDSCAPE)} size="wide" />
        </Figure>
        {onePercent && (
          <Finding lead={count(onePercent.all_conditions_lower_bound)} rest="high-confidence candidates">
            The midpoint alone leaves {count(onePercent.all_conditions)}—more reach, weaker support.
          </Finding>
        )}
      </Section>

      {escape && (
        <Section
          flat
          number="02"
          kicker="Base editing"
          title={
            <>
              Editors can reach nearly a third. <em>Reach is only geometry.</em>
            </>
          }
          intro="This is editor placement geometry—not efficiency, delivery, off-targets, or clinical eligibility."
        >
          <Figure
            kind="rule"
            title="Placement reaches 22,042 variants; removing extra nearby edits retains 18,249"
            caption={`Complete partition of ${count(escape.scoreable)} variants scoreable against ${escape.panel}`}
            verdict="GREEN MEANS A CANDIDATE PLACEMENT · Exact restoration is preferable; neither green segment demonstrates efficacy."
            note={
              <>
                Exact means the edit restores the reference amino acid. Alternative means it creates
                a different sense codon. A bystander-free placement has no additional editable base
                inside the declared activity window; it still needs experimental efficiency and
                specificity testing.
              </>
            }
          >
            <Chart svg={editingReach(escape)} size="wide" />
          </Figure>
          <Finding lead={percent(escape.reachable / escape.scoreable)} rest="geometrically reachable">
            {count(escape.exact)} restore the exact amino acid.
          </Finding>
        </Section>
      )}

      <Section
        number="03"
        kicker="Transcript survival"
        title={
          <>
            Whether the RNA survives depends on <em>which rule you ask.</em>
          </>
        }
        intro="Adding start-proximal and long-exon exceptions reveals variants whose decay prediction depends on the rule."
      >
        <Figure
          kind="rule"
          title="The full nonsense-mediated decay rule expands predicted escape from 10.1% to 30.3%"
          caption={`A complete partition of ${count(nmd.scoreable)} scoreable variants`}
          verdict="GREEN KEEPS THE READTHROUGH ROUTE OPEN · Grey predicts RNA loss; orange marks rule-sensitive uncertainty."
          note={
            <>
              The orange band is not an error bar. It is the exact set whose classification changes
              when the declared biological exceptions are added. Those{" "}
              <strong>{count(nmd.disagree)} variants</strong> are the most valuable targets for
              orthogonal transcript-stability measurements.
            </>
          }
        >
          <Chart svg={nmdRuleExpansion(nmd)} size="wide" />
        </Figure>
        <Finding lead={percent(nmd.disagree_fraction)} rest="rule-sensitive">
          One in five classifications changes. Both verdicts stay visible.
        </Finding>
      </Section>

      <Section
        flat
        number="04"
        kicker="Suppressor tRNA"
        title={
          <>
            Two designs reach a third. <em>The tail still matters.</em>
          </>
        }
        intro="Each design adds only exact-restoration variants not reached by earlier designs."
      >
        <Figure
          kind="rule"
          title="Nineteen stop-residue designs close the exact-restoration universe"
          caption="Marginal variants per design, with cumulative coverage over the same scoreable denominator"
          legend={[
            { className: "held", label: "new variants added by this design" },
            { className: "measured", label: "cumulative exact-restoration coverage" },
          ]}
          verdict="Coverage is sequence logic, not therapeutic evidence. No suppressor tRNA in this project has a qualifying ribosome-profiling measurement."
          note={
            <>
              The first design, UAG-Q, reaches {count(frontier[0]?.marginal ?? 0)} variants. UGA-R
              takes the pair to <strong>{count(firstTwo?.cumulative ?? 0)}</strong>. The final designs
              add few variants, but they are the only exact-restoration route for those stop-residue
              combinations. UAG-Q means a design for UAG stops that restores glutamine (Q).
            </>
          }
        >
          <Chart svg={suppressorFrontier(frontier)} size="wide" />
        </Figure>
        <Finding lead={percent(firstTwo?.cumulative_fraction)} rest="covered by two designs">
          Best starting portfolio, not evidence that either construct works or is safe.
        </Finding>
      </Section>

      <footer className="provenance">
        <span>ClinVar {provenance.clinvar_release}</span>
        <span>{provenance.qualifying_variants.toLocaleString("en-US")} qualifying variants</span>
        <span>analysis {provenance.commit}</span>
      </footer>
    </div>
  );
}
