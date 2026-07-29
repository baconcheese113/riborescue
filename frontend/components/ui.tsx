import type { ReactNode } from "react";

export type Provenance = "measured" | "rule" | "predicted" | "absent";

const PROVENANCE_WORD: Record<Provenance, string> = {
  measured: "measured",
  rule: "rule",
  predicted: "predicted",
  absent: "absent",
};

/** Where a number came from, as a quiet pill. The modifier sets only colour; the dot inherits it. */
export function Tag({ kind, children }: { kind: Provenance; children?: ReactNode }) {
  return <span className={`tag ${kind}`}>{children ?? PROVENANCE_WORD[kind]}</span>;
}

/** Numbered, headed by the question it answers, and alternately raised or flat down the page. */
export function Section({
  number,
  kicker,
  title,
  intro,
  flat,
  children,
}: {
  number: string;
  kicker: string;
  title: ReactNode;
  intro: ReactNode;
  flat?: boolean;
  children: ReactNode;
}) {
  const id = `s${number}`;
  return (
    <section className={flat ? "section flat" : "section"} aria-labelledby={id}>
      <div className="section-head">
        <div>
          <span className="section-number">{number}</span>
          <p className="kicker">{kicker}</p>
          <h2 id={id}>{title}</h2>
        </div>
        <p className="section-intro">{intro}</p>
      </div>
      {children}
    </section>
  );
}

/** Title says what was found, sub-caption says on what denominator, note carries the limit. */
export function Figure({
  title,
  caption,
  kind,
  note,
  legend,
  verdict,
  verdictKind = "caveat",
  children,
}: {
  title: string;
  caption: string;
  kind: Provenance;
  note: ReactNode;
  legend?: { className?: string; label: string }[];
  verdict?: ReactNode;
  verdictKind?: "pass" | "caveat" | "fail";
  children: ReactNode;
}) {
  return (
    <figure className="figure">
      <div className="figure-head">
        <div>
          <h3>{title}</h3>
          <p>{caption}</p>
        </div>
        <Tag kind={kind} />
      </div>
      {children}
      {legend && (
        <ul className="chart-legend">
          {legend.map((item) => (
            <li key={item.label}>
              <i className={item.className ? `swatch ${item.className}` : "swatch"} />
              {item.label}
            </li>
          ))}
        </ul>
      )}
      {verdict && <p className={`verdict ${verdictKind}`}>{verdict}</p>}
      <figcaption className="figure-note">
        <details>
          <summary>Details &amp; limits</summary>
          <div>{note}</div>
        </details>
      </figcaption>
    </figure>
  );
}

/** Build-time SVG. The colours inside are CSS variables, so this markup serves both themes. */
export function Chart({ svg, size }: { svg: string; size?: "compact" | "wide" }) {
  return (
    <div className="chart-frame">
      <span className="chart-pan-hint" aria-hidden="true">
        Swipe chart to see the full axis →
      </span>
      <div className="chart-scroll">
        <div
          className={size ? `chart-host ${size}` : "chart-host"}
          dangerouslySetInnerHTML={{ __html: svg }}
        />
      </div>
    </div>
  );
}

/** The inverted card: one number, one clause of it in the alarm colour. */
export function Finding({ lead, rest, children }: { lead: string; rest: string; children: ReactNode }) {
  return (
    <div className="finding">
      <b>
        <span>{lead}</span>
        <em>{rest}</em>
      </b>
      <p>{children}</p>
    </div>
  );
}

export function ScopeStrip({ items }: { items: { value: string; label: string }[] }) {
  return (
    <div className="scope-strip">
      {items.map((item) => (
        <div key={item.label}>
          <b>{item.value}</b>
          <span>{item.label}</span>
        </div>
      ))}
    </div>
  );
}

/** Four shapes, not four colours: filled means we have it, hollow means we inferred it, dashed
 * means it is not there. The distinction survives greyscale and colour-vision deficiency. */
export function EvidenceKey() {
  const kinds: { kind: Provenance; title: string; body: string }[] = [
    { kind: "measured", title: "Measured", body: "Observed here" },
    { kind: "rule", title: "Rule", body: "Derived from sequence" },
    { kind: "predicted", title: "Predicted", body: "Model estimate" },
    { kind: "absent", title: "Absent", body: "Unknown, not zero" },
  ];
  return (
    <ul className="evidence-key">
      {kinds.map((k) => (
        <li key={k.kind} className={k.kind}>
          <b>{k.title}</b>
          <span>{k.body}</span>
        </li>
      ))}
    </ul>
  );
}
