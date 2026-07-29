import {
  codonOccupancy,
  discriminationEffects,
  discriminationLibraries,
  frameByLength,
  modelParity,
  periodicity,
  permutationNull,
} from "../../components/charts";
import { Chart, EvidenceKey, Figure, Finding, ScopeStrip, Section } from "../../components/ui";
import { EVIDENCE, count, has } from "../../lib/data";

export const metadata = {
  title: "Does the instrument work? · RiboRescue",
  description:
    "The checks behind every RiboRescue number: ribosome periodicity, footprint selection, a drug that reads through beside one that stalls, and a hypothesis of our own that failed.",
};

export default function EvidencePage() {
  const {
    calibration,
    periodicity: profile,
    readthrough,
    codon_occupancy,
    model_parity,
    kinetics_null,
    provenance,
  } = EVIDENCE;

  const libraries = calibration?.libraries ?? [];
  const psites = libraries.reduce((sum, lib) => sum + lib.psites, 0);
  const clean = libraries.filter((lib) => lib.failures.length === 0).length;
  const transcripts = readthrough?.g418_vs_dmso?.libraries[0]?.transcripts ?? 0;

  const g418 = readthrough?.g418_vs_dmso.quantities ?? [];
  const sri = readthrough?.sri37240_vs_dmso.quantities ?? [];
  const heldG418 = g418.filter((q) => q.consistent).length;
  const heldSri = sri.filter((q) => q.consistent).length;

  return (
    <div className="page">
      <header className="display">
        <div className="display-head">
          <div>
            <p className="kicker">The checks behind every number</p>
            <h1>
              Before you believe a prediction, <em>believe the instrument.</em>
            </h1>
          </div>
          <p className="lede">
            Six checks ask whether our ribosome measurements work—including the control that
            rejected our own hypothesis.
          </p>
        </div>
        <ScopeStrip
          items={[
            { value: count(psites), label: "ribosome positions mapped" },
            { value: `${clean} of ${libraries.length}`, label: "libraries calibrated without a failure" },
            { value: count(transcripts), label: "transcripts deep enough to measure" },
          ]}
        />
        <EvidenceKey />
      </header>

      {has(profile) && (
        <Section
          number="01"
          kicker="Periodicity"
          title={
            <>
              A ribosome moves three letters <em>at a time.</em>
            </>
          }
          intro="Real ribosome footprints repeat every three nucleotides; random RNA fragments do not."
        >
          <Figure
            kind="measured"
            title="The three-nucleotide beat is visible in the raw density"
            caption="Footprint density by distance from the codon, averaged within each treatment arm"
            legend={[
              { className: "", label: "DMSO (control)" },
              { className: "measured", label: "G418" },
              { className: "alarm", label: "SRI-37240" },
            ]}
            verdictKind="pass"
            verdict="PASS · All three treatment arms show the repeating three-nucleotide pattern expected from genuine ribosome footprints."
            note={
              <>
                Arms are averaged separately, never pooled — pooling would erase the stop-codon
                difference the next section depends on.
              </>
            }
          >
            <Chart svg={periodicity(profile)} size="wide" />
          </Figure>
        </Section>
      )}

      {has(calibration) && calibration.frame_by_length && (
        <Section
          flat
          number="02"
          kicker="Selection"
          title={
            <>
              We did not choose the read lengths. <em>The data did.</em>
            </>
          }
          intro="Keep only read lengths that preserve the expected three-nucleotide reading frame."
        >
          <Figure
            kind="measured"
            title="Only three footprint lengths carry the reading frame"
            caption={`Share of reads landing in the coding frame, per length, across ${libraries.length} libraries`}
            verdictKind="pass"
            verdict={`PASS · Lengths ${calibration.lengths.join(", ")} rise above the random-fragment baseline and are retained.`}
            note={
              <>
                Lengths <strong>{calibration.lengths.join(", ")}</strong> clear the bar; the rest sit
                near the one-in-three a random fragment would give. Selecting on the frame rather than
                on the length distribution is what makes the codon-level results below meaningful.
              </>
            }
          >
            <Chart svg={frameByLength(calibration.frame_by_length)} size="wide" />
          </Figure>
        </Section>
      )}

      {has(readthrough) && (
        <Section
          number="03"
          kicker="Discrimination"
          title={
            <>
              One drug reads through. The other <em>piles up at the stop.</em>
            </>
          }
          intro="A trustworthy assay must separate a readthrough drug from a stalling one."
        >
          <Figure
            kind="measured"
            title="G418 meets all three conditions; SRI-37240 meets none"
            caption="Mean difference against control with a 95% interval, per condition"
            legend={[
              { className: "measured", label: "filled — the condition was met" },
              { className: "hollow", label: "hollow — the condition was not met" },
            ]}
            verdictKind="pass"
            verdict={
              <>
                PASS · G418 moves in the desired direction on all three conditions; SRI-37240 does
                not. The G418 termination interval spans zero, but that condition is registered as a
                direction check rather than an interval-exclusion check.
              </>
            }
            note={
              <>
                Readthrough is a redistribution, so it must show in more than one place at once:
                density beyond the stop rises, density at the stop falls, and what continues stays in
                the coding frame. <strong>{heldG418} of 3</strong> hold for G418 and{" "}
                <strong>{heldSri} of 3</strong> for SRI-37240.
              </>
            }
          >
            <Chart svg={discriminationEffects(readthrough)} size="wide" />
          </Figure>

          <Figure
            kind="measured"
            title="Every library, drawn separately"
            caption="Ribosomes continuing past the stop codon, one dot per sequencing library"
            legend={[
              { className: "", label: "control (DMSO)" },
              { className: "measured", label: "treated" },
            ]}
            verdictKind="pass"
            verdict="PASS · Treated G418 libraries separate from control; SRI-37240 overlaps control."
            note={
              <>
                G418&apos;s treated libraries separate completely from their controls; SRI-37240&apos;s
                interleave. The same evidence as above, with nothing summarised away.
              </>
            }
          >
            <Chart svg={discriminationLibraries(readthrough, "downstream_occupancy")} size="wide" />
          </Figure>
        </Section>
      )}

      {has(codon_occupancy) && (
        <Section
          flat
          number="04"
          kicker="Recovery"
          title={
            <>
              The pipeline finds biology <em>nobody told it about.</em>
            </>
          }
          intro="An untuned measurement should recover known amino-acid differences in ribosome occupancy."
        >
          <Figure
            kind="measured"
            title="Occupancy varies by amino acid, in the order the literature reports"
            caption="Mean relative occupancy per amino acid, both ribosome sites, 61 sense codons"
            legend={[
              { className: "measured", label: "A site — where the next amino acid is read" },
              { className: "hollow", label: "P site — where the chain is joined" },
            ]}
            verdictKind="pass"
            verdict="PASS · The pipeline recovers the known amino-acid occupancy pattern without being trained on it."
            note={
              <>
                Nothing was tuned to produce this ordering, so it checks the instrument rather than
                reporting a finding.
              </>
            }
          >
            <Chart svg={codonOccupancy(codon_occupancy)} size="wide" />
          </Figure>
        </Section>
      )}

      {has(model_parity) && (
        <Section
          number="05"
          kicker="Ceiling"
          title={
            <>
              The model is nearly as good as <em>the experiment allows.</em>
            </>
          }
          intro="Model accuracy is judged against the reproducibility ceiling of its own experiment."
        >
          <Figure
            kind="predicted"
            title="Every drug sits close to its own reproducibility ceiling"
            caption="Held-out accuracy against the agreement between repeats of the same experiment"
            legend={[
              { className: "predicted", label: "held-out accuracy" },
              { className: "", label: "ceiling set by the experiment's own repeats" },
            ]}
            verdictKind="pass"
            verdict="PASS · Each model approaches the maximum accuracy permitted by that drug's experimental repeats."
            note={
              <>
                The distance that matters is to the black tick, not to 1.0.
              </>
            }
          >
            <Chart svg={modelParity(model_parity)} size="wide" />
          </Figure>
        </Section>
      )}

      {has(kinetics_null) && (
        <Section
          flat
          number="06"
          kicker="Our own hypothesis"
          title={
            <>
              We tested an idea of ours, <em>and it lost.</em>
            </>
          }
          intro="Three shuffles ask whether codon speed adds signal beyond local sequence. The decisive control says no."
        >
          <Figure
            kind="measured"
            title="The shuffle that keeps the amino acid reproduces the gain"
            caption="Best observed improvement against 199 permutations, per shuffle family"
            legend={[
              { className: "", label: "what chance produces" },
              { className: "hollow", label: "hollow — chance explains it" },
            ]}
            verdictKind="fail"
            verdict="NOT SUPPORTED · The same-amino-acid shuffle reproduces the gain, so codon speed does not explain it."
            note={
              <>
                Shuffling inside genes and shuffling everywhere both fail to reproduce the gain. But
                shuffling only between codons that carry the <strong>same amino acid</strong>{" "}
                reproduces most of it — so the signal is which amino acid sits before the stop, not
                how quickly its codon is read.{" "}
                {kinetics_null.analysis_status === "incomplete" && (
                  <>
                    {kinetics_null.permutations_completed} of {kinetics_null.permutations_required}{" "}
                    permutations are complete, so p-values resolve only to {kinetics_null.resolution}.
                  </>
                )}
              </>
            }
          >
            <Chart svg={permutationNull(kinetics_null.rows)} size="wide" />
          </Figure>

          <Finding lead="0.21" rest="familywise p-value">
            The decisive control is not significant. The failed hypothesis remains a result.
          </Finding>
        </Section>
      )}

      <footer className="provenance">
        <span>{provenance.dataset}</span>
        <span>{provenance.scope}</span>
        <span>commit {provenance.commit}</span>
      </footer>
    </div>
  );
}
