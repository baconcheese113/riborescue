"use client";

import { useEffect, useMemo, useState } from "react";
import type { Variant, WebTable } from "../types";

const percent = (value: number | null) =>
  value === null ? "—" : `${(value * 100).toFixed(1)}%`;

// Every therapy for the selected variant as a card: its predicted readthrough and interval, or an
// explicit no-prediction state with the reason. The scale is shared so the bars are comparable.
function TherapyCards({ variant }: { variant: Variant }) {
  const scale = Math.max(0.001, ...variant.therapies.map((t) => t.high ?? 0));
  return (
    <div className="cards">
      {variant.therapies.map((t) => (
        <div className={t.available ? "card" : "card muted"} key={t.id}>
          <div className="card-head">
            <span className="card-name">{t.name ?? t.id}</span>
            {t.available ? (
              <span className="card-num">{percent(t.readthrough)}</span>
            ) : (
              <span className="card-na">no prediction</span>
            )}
          </div>
          {t.available ? (
            <>
              <span className="track">
                <span
                  className="ci"
                  style={{
                    left: `${((t.low ?? 0) / scale) * 100}%`,
                    width: `${(((t.high ?? 0) - (t.low ?? 0)) / scale) * 100}%`,
                  }}
                />
                <span className="point" style={{ left: `${((t.readthrough ?? 0) / scale) * 100}%` }} />
              </span>
              <span className="card-foot">lower bound {percent(t.low)}</span>
            </>
          ) : (
            <span className="card-foot">{t.reason ?? "unavailable"}</span>
          )}
        </div>
      ))}
    </div>
  );
}

// The modalities that could target this stop, each a geometric candidate rather than a therapy: is
// there a predicted readthrough, and can a base editor be placed on the stop. Neither is eligibility.
function RescueRoutes({ variant }: { variant: Variant }) {
  const e = variant.editing;
  const exact = e?.reach_class === "base_editable_exact";
  const alternative = e?.reach_class === "base_editable_alternative";
  return (
    <div className="routes">
      <div className="route">
        <div className="route-head">
          <span className="route-name">Readthrough</span>
          {variant.best ? (
            <span className="pill good">candidate</span>
          ) : (
            <span className="pill muted">no predicted readthrough</span>
          )}
        </div>
        <div className="route-body">
          {variant.best ? (
            <span className="slot-note">
              best predicted {percent(variant.best.readthrough)} with {variant.best.therapy}, lower
              bound {percent(variant.best.low)} — the therapies compared below
            </span>
          ) : (
            <span className="slot-note">no therapy has a prediction for this variant</span>
          )}
        </div>
      </div>

      <div className="route">
        <div className="route-head">
          <span className="route-name">Base editing</span>
          {e && (exact || alternative) ? (
            <span className="pill good">base-editable</span>
          ) : e ? (
            <span className="pill bad">not base-editable under this panel</span>
          ) : (
            <span className="pill muted">not evaluated</span>
          )}
        </div>
        <div className="route-body">
          {e && (exact || alternative) ? (
            <>
              <span className="slot-note">
                {exact ? (
                  <>
                    restores the <b>original {variant.residue}</b> exactly
                  </>
                ) : (
                  <>
                    removes the stop by <b>alternative sense-codon conversion</b> — not the original{" "}
                    {variant.residue}, and protein function is not established
                  </>
                )}{" "}
                · {e.editor} on the {e.strand} strand
                {e.bystander_free ? (
                  <>
                    {" "}
                    · <span className="pill good">no coding bystanders</span>
                  </>
                ) : (
                  <>
                    {" "}
                    · <span className="pill bad">has coding bystanders</span>
                  </>
                )}
              </span>
              <details className="route-more">
                <summary>Guide detail</summary>
                <dl className="route-detail">
                  <div>
                    <dt>Editor</dt>
                    <dd>{e.editor}</dd>
                  </div>
                  <div>
                    <dt>Strand</dt>
                    <dd>{e.strand}</dd>
                  </div>
                  <div>
                    <dt>Editing-window position</dt>
                    <dd>{e.window_position ?? "—"}</dd>
                  </div>
                  <div>
                    <dt>Candidate guides</dt>
                    <dd>{e.candidate_guides}</dd>
                  </div>
                </dl>
                <p className="route-caveat">
                  Geometric candidate only — whether an editor can be placed on the stop. Not editing
                  efficiency, off-target activity, product purity, delivery, tissue access, splice
                  consequence, or clinical eligibility.
                </p>
              </details>
            </>
          ) : e ? (
            <span className="slot-note">
              no BE4max or ABE7.10 guide with an NGG PAM places the stop in its editing window
            </span>
          ) : (
            <span className="slot-note">base-editing reachability was not exported for this variant</span>
          )}
        </div>
      </div>
    </div>
  );
}

// The evidence beside the therapies, each slot saying what is known and — where a layer is not built
// — that it is not, rather than leaving a blank the reader fills in optimistically.
function Evidence({ variant }: { variant: Variant }) {
  return (
    <dl className="evidence">
      <div className="slot">
        <dt>Escapes NMD</dt>
        <dd>
          {variant.nmd ? (
            <>
              <span className={`pill ${variant.nmd.escape_guideline ? "good" : "bad"}`}>
                guideline: {variant.nmd.escape_guideline ? "escape" : "decay"}
              </span>{" "}
              <span className={`pill ${variant.nmd.escape_full_rules ? "good" : "bad"}`}>
                full rules: {variant.nmd.escape_full_rules ? "escape" : "decay"}
              </span>{" "}
              {variant.nmd.aenmd?.available && (
                <>
                  <span className={`pill ${variant.nmd.aenmd.escape ? "good" : "bad"}`}>
                    aenmd: {variant.nmd.aenmd.escape ? "escape" : "decay"}
                  </span>{" "}
                </>
              )}
              {variant.nmd.nmdetective?.available && variant.nmd.nmdetective.efficiency !== null && (
                <>
                  <span className={`pill ${variant.nmd.nmdetective.efficiency <= 0 ? "good" : "bad"}`}>
                    NMDetective-AI: {variant.nmd.nmdetective.efficiency.toFixed(2)}
                  </span>{" "}
                </>
              )}
              {variant.nmd.disagree ? (
                <span className="slot-note">
                  the two rule sets <b>disagree</b> here — the fuller rules escape a stop the 50-nt
                  guideline calls decay
                  {variant.nmd.rules.long_exon ? ", via the long-exon rule" : ""}
                  {variant.nmd.rules.start_proximal ? ", via the start-proximal rule" : ""}
                </span>
              ) : (
                <span className="slot-note">the two rule sets agree here</span>
              )}
              {variant.nmd.aenmd &&
                (variant.nmd.aenmd.available ? (
                  <span className="slot-note">
                    aenmd, the published rule tool,{" "}
                    {variant.nmd.aenmd.escape === variant.nmd.escape_full_rules
                      ? "agrees with the full rule set"
                      : "differs from the full rule set"}
                  </span>
                ) : (
                  <span className="slot-note">
                    aenmd did not score this variant ({variant.nmd.aenmd.reason.replace(/_/g, " ")})
                  </span>
                ))}
              {variant.nmd.nmdetective?.available && variant.nmd.nmdetective.efficiency !== null && (
                <span className="slot-note">
                  NMDetective-AI, the deep model, scores NMD efficiency{" "}
                  {variant.nmd.nmdetective.efficiency.toFixed(2)} — higher means more decay, so this
                  variant leans {variant.nmd.nmdetective.efficiency <= 0 ? "escape" : "decay"} (a
                  continuous score, not a yes/no call)
                </span>
              )}
            </>
          ) : (
            <>
              <span className={`pill ${variant.escapes_decay ? "good" : "bad"}`}>
                {variant.escapes_decay ? "predicted escape" : "predicted decay"}
              </span>{" "}
              <span className="slot-note">rule-based (50-nt last-junction rule)</span>
            </>
          )}
        </dd>
      </div>
      <div className="slot">
        <dt>Residue compatibility</dt>
        <dd>
          {percent(variant.tolerable_insertion_share)}{" "}
          <span className="slot-note">share of insertions BLOSUM-compatible with the native residue</span>
        </dd>
      </div>
      <div className="slot">
        <dt>Suppressor tRNA</dt>
        <dd>
          {variant.suppressor ? (
            <>
              <code>{variant.suppressor.design}</code>{" "}
              <span className="slot-note">reinserts the native {variant.residue} exactly</span>
            </>
          ) : (
            <span className="slot-note">no suppressor design maps to this stop</span>
          )}
        </dd>
      </div>
      <div className="slot">
        <dt>Native-stop safety</dt>
        <dd>
          <span className="slot-note">
            measured for G418 in HEK293T only; the other therapies have no matched empirical atlas
          </span>
        </dd>
      </div>
      <div className="slot unbuilt">
        <dt>Protein function</dt>
        <dd>
          <span className="pill muted">not yet modelled</span>{" "}
          <span className="slot-note">conservation and domain context are a planned layer</span>
        </dd>
      </div>
    </dl>
  );
}

export default function PatientPage() {
  const [data, setData] = useState<WebTable | null>(null);
  const [error, setError] = useState<string | null>(null);
  const [query, setQuery] = useState("");
  const [selected, setSelected] = useState<string | null>(null);

  useEffect(() => {
    fetch("/riborescue.json")
      .then((r) => {
        if (!r.ok) throw new Error(`riborescue.json returned ${r.status}`);
        return r.json();
      })
      .then(setData)
      .catch((e) => setError(String(e?.message ?? e)));
  }, []);

  const matches = useMemo(() => {
    if (!data) return [];
    const q = query.trim().toLowerCase();
    if (!q) return data.variants.slice(0, 12);
    return data.variants
      .filter((v) =>
        [v.gene, v.stop, v.residue, v.id, String(v.protein_position ?? "")]
          .join(" ")
          .toLowerCase()
          .includes(q),
      )
      .slice(0, 24);
  }, [data, query]);

  const variant = useMemo(
    () => data?.variants.find((v) => v.id === selected) ?? null,
    [data, selected],
  );

  if (error)
    return (
      <div className="wrap">
        <h1>Patient &amp; clinic</h1>
        <div className="banner" style={{ borderLeftColor: "var(--bad)" }}>
          <span className="tag" style={{ color: "var(--bad)" }}>Could not load data</span>
          <p>{error}. Generate it with <code>pixi run export-web-example</code>, then reload.</p>
        </div>
      </div>
    );

  return (
    <div className="wrap">
      <h1>Patient &amp; clinic</h1>
      <p className="lede">
        Find a variant, then read how each readthrough therapy is predicted to compare, with the
        evidence beside it. <b>This is a research tool, not medical advice</b>, and it searches an
        example subset of variants, not the full ClinVar set.
      </p>

      <div className="controls">
        <input
          type="search"
          placeholder="Search by gene, variant, stop or residue…"
          value={query}
          onChange={(e) => setQuery(e.target.value)}
          aria-label="Search variants"
        />
      </div>

      {!data ? (
        <p>Loading…</p>
      ) : (
        <div className="patient-body">
          <ul className="hits" aria-label="Matching variants">
            {matches.map((v) => (
              <li key={v.id}>
                <button
                  className={v.id === selected ? "hit active" : "hit"}
                  onClick={() => setSelected(v.id)}
                  aria-current={v.id === selected ? "true" : undefined}
                >
                  <span className="hit-gene">{v.gene}</span>
                  <span className="hit-meta">
                    {v.stop} · {v.residue}
                    {v.protein_position ? ` · codon ${v.protein_position}` : ""}
                  </span>
                </button>
              </li>
            ))}
            {matches.length === 0 && <li className="hit-none">No variant matches “{query}”.</li>}
          </ul>

          <section className="detail-pane" aria-live="polite">
            {variant ? (
              <>
                <h2>
                  {variant.gene} <span className="detail-sub">{variant.stop} at {variant.residue}
                  {variant.protein_position ? `, codon ${variant.protein_position}` : ""}</span>
                </h2>
                <h3 className="cards-title">Possible rescue routes</h3>
                <RescueRoutes variant={variant} />
                <h3 className="cards-title">Every therapy, predicted</h3>
                <TherapyCards variant={variant} />
                <h3 className="cards-title">Evidence</h3>
                <Evidence variant={variant} />
              </>
            ) : (
              <p className="detail-empty">Select a variant to see its therapies and evidence.</p>
            )}
          </section>
        </div>
      )}
    </div>
  );
}
