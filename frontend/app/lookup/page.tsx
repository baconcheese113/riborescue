import { Lookup } from "../../components/lookup";
import { ScopeStrip } from "../../components/ui";
import { RESEARCH, count } from "../../lib/data";

export const metadata = {
  title: "Look up a variant · RiboRescue",
  description:
    "Every evidence layer for one nonsense variant: which rescue routes are open, all six compounds with their uncertainty, and whether the message survives at all.",
};

export default function LookupPage() {
  const provenance = RESEARCH.provenance;
  const escape = RESEARCH.escape;
  const nmd = RESEARCH.nmd;

  return (
    <div className="page">
      <header className="display" style={{ padding: "56px 0 12px" }}>
        <p className="kicker">One variant at a time</p>
        <h1>
          Four questions about a variant, <em>answered separately.</em>
        </h1>
        <p className="lede">
          Every layer of evidence behind one stop — including the layer that does not exist yet.
        </p>
        <ScopeStrip
          items={[
            { value: count(provenance.scoreable_variants), label: "variants, every one scored" },
            { value: count(nmd.escape_guideline), label: "whose message is predicted to survive" },
            { value: count(escape?.reachable ?? 0), label: "an editor can be placed on" },
          ]}
        />
      </header>

      <section className="section">
        <Lookup />
      </section>

      <footer className="provenance">
        <span>all {count(provenance.scoreable_variants)} scoreable variants · ClinVar {provenance.clinvar_release}</span>
        <span>research use only · not medical advice</span>
      </footer>
    </div>
  );
}
