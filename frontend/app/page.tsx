import Link from "next/link";

// The landing page routes a visitor by what they came for, and says plainly what this is and is not
// before any number appears. Three destinations, one disclaimer, no data loaded here.
const ROLES: { href: string; title: string; blurb: string }[] = [
  {
    href: "/patient",
    title: "Patient & clinic",
    blurb:
      "Start from a gene or variant and see, for one nonsense variant, how each readthrough therapy is predicted to compare — with what is known and what is not yet modelled shown side by side.",
  },
  {
    href: "/researcher",
    title: "Researcher",
    blurb:
      "Coverage across the ClinVar nonsense-variant set: which suppressor-tRNA designs reach the most variants, genes and ClinVar conditions, and how far each condition's variants are covered — every figure with its denominator.",
  },
  {
    href: "/explorer",
    title: "Advanced data explorer",
    blurb:
      "The per-variant table: every therapy's predicted readthrough and interval, the suppressor tRNA, the NMD and residue-compatibility flags, and the native-stop safety atlas. A technical view over an example slice.",
  },
];

export default function Home() {
  return (
    <div className="wrap">
      <header className="home-hero">
        <h1>RiboRescue</h1>
        <p className="home-lede">
          Matching a patient&apos;s nonsense variant to candidate readthrough therapies. Every number
          here is a model prediction or a measured quantity with a stated scope — <b>not medical
          advice, not a claim any therapy works</b>, and not a clinical recommendation. A research
          project, shown for inspection.
        </p>
      </header>

      <div className="role-grid">
        {ROLES.map((role) => (
          <Link key={role.href} href={role.href} className="role-card">
            <h2>{role.title}</h2>
            <p>{role.blurb}</p>
            <span className="role-go">Open →</span>
          </Link>
        ))}
      </div>

      <p className="home-foot">
        Predictions rest on ClinVar, MANE and a reproduced baseline efficacy model; the safety atlas
        rests on public Ribo-seq measured for one drug in one cell line. What is built and what is
        deliberately not is stated on each page.
      </p>
    </div>
  );
}
