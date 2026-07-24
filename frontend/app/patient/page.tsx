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
            <span className="card-name">{t.id}</span>
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

// The evidence beside the therapies, each slot saying what is known and — where a layer is not built
// — that it is not, rather than leaving a blank the reader fills in optimistically.
function Evidence({ variant }: { variant: Variant }) {
  return (
    <dl className="evidence">
      <div className="slot">
        <dt>Escapes NMD</dt>
        <dd>
          <span className={`pill ${variant.escapes_decay ? "good" : "bad"}`}>
            {variant.escapes_decay ? "predicted yes" : "predicted no"}
          </span>{" "}
          <span className="slot-note">rule-based (50-nt last-junction rule), not an ensemble yet</span>
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
