"use client";

import { useEffect, useMemo, useState } from "react";

/** The whole scored set, columnar. Fetched rather than bundled: 70,376 variants is a payload this
 * page asks for when a reader opens it, not something every other page carries. */
type Index = {
  scales: { readthrough: number; tolerance: number };
  therapies: string[];
  genes: string[];
  conditions: string[];
  residues: string[];
  reach_classes: string[];
  ids: string[];
  gene: number[];
  position: number[];
  stop: number[];
  residue: number[];
  escapes: number[];
  escapes_full: number[];
  stars: number[];
  best: number[];
  best_low: number[];
  tolerance: number[];
  extra_conditions: number[];
  condition: number[];
  reach: number[];
  arm_mid: number[];
  arm_low: number[];
  arm_high: number[];
};

const STOPS = ["UAG", "UGA", "UAA"];
const ABSENT = -1;
const ROWS_SHOWN = 60;

const THERAPY_NAME: Record<string, string> = {
  CC90009: "CC-90009",
  Clitocine: "Clitocine",
  DAP: "2,6-diaminopurine",
  G418: "G418 (geneticin)",
  SJ6986: "SJ6986",
  SRI: "SRI-41315",
};

const REACH_WORD: Record<string, string> = {
  base_editable_exact: "Reachable, restores the original amino acid",
  base_editable_alternative: "Reachable, but inserts a different amino acid",
  not_base_editable_under_panel: "No editor fits under the declared panel",
};

const pct = (v: number | null, digits = 2) => (v === null ? "—" : `${(v * 100).toFixed(digits)}%`);

function Ring({ value, tone, absent }: { value: number | null; tone: string; absent?: boolean }) {
  const r = 26;
  const c = 2 * Math.PI * r;
  const filled = absent || value === null ? 0 : Math.max(0, Math.min(1, value));
  return (
    <svg className={`ring ring-${tone}`} viewBox="0 0 64 64" aria-hidden="true">
      <circle className="ring-track" cx="32" cy="32" r={r} strokeDasharray={absent ? "3 4" : undefined} />
      {!absent && (
        <circle
          className="ring-fill"
          cx="32"
          cy="32"
          r={r}
          strokeDasharray={`${filled * c} ${c}`}
          transform="rotate(-90 32 32)"
        />
      )}
    </svg>
  );
}

function Layer({
  label,
  verdict,
  detail,
  value,
  tone,
  absent,
}: {
  label: string;
  verdict: string;
  detail: string;
  value: number | null;
  tone: string;
  absent?: boolean;
}) {
  return (
    <article className={`layer${absent ? " layer-absent" : ""}`}>
      <div className="layer-head">
        <p className="route-kicker">{label}</p>
        <Ring value={value} tone={tone} absent={absent} />
      </div>
      <p className="layer-verdict">{verdict}</p>
      <p className="layer-detail">{detail}</p>
    </article>
  );
}

type Arm = { name: string; mid: number; low: number; high: number };

function Intervals({ arms }: { arms: Arm[] }) {
  const usable = arms.filter((a) => a.mid >= 0);
  if (!usable.length) return <p className="intervals-note">No therapy was scored for this variant.</p>;
  const max = Math.max(0.00001, ...usable.map((a) => a.high));
  const scale = (v: number) => `${(v / max) * 100}%`;
  const ranked = [...usable].sort((a, b) => b.mid - a.mid);
  const [first, second] = ranked;
  const lo = Math.max(first?.low ?? 0, second?.low ?? 0);
  const hi = Math.min(first?.high ?? 0, second?.high ?? 0);
  const overlaps = Boolean(second) && hi > lo;
  return (
    <div className="intervals">
      <div className="interval-rows">
        {overlaps && (
          <div className="interval-overlap" style={{ left: scale(lo), width: scale(hi - lo) }} aria-hidden="true" />
        )}
        {ranked.map((a) => (
          <div className="interval-row" key={a.name}>
            <span className="interval-name">{a.name}</span>
            <span className="interval-track">
              <span
                className="interval-bar"
                style={{ left: scale(a.low), width: scale(a.high - a.low) }}
                title={`${a.name}: ${pct(a.mid)} (${pct(a.low)} to ${pct(a.high)})`}
              />
              <span className="interval-point" style={{ left: scale(a.mid) }} />
            </span>
            <span className="interval-value">{pct(a.mid)}</span>
          </div>
        ))}
      </div>
      <p className="intervals-note">
        {overlaps ? (
          <>
            <b>
              {first.name} and {second.name} overlap between {pct(lo)} and {pct(hi)}.
            </b>{" "}
            True for all but five of the 70,376 scored variants.
          </>
        ) : (
          <>
            <b>The leading interval clears the runner-up here.</b> True for about five variants in
            seventy thousand.
          </>
        )}
      </p>
    </div>
  );
}

const haystack = (data: Index, i: number) =>
  `${data.genes[data.gene[i]]} ${STOPS[data.stop[i]]} ${
    data.condition[i] === ABSENT ? "" : data.conditions[data.condition[i]]
  }`.toLowerCase();

export function Lookup() {
  const [data, setData] = useState<Index | null>(null);
  const [error, setError] = useState<string | null>(null);
  const [query, setQuery] = useState("");
  const [picked, setPicked] = useState(0);

  useEffect(() => {
    fetch("/riborescue_index.json")
      .then((r) => (r.ok ? r.json() : Promise.reject(new Error(`index returned ${r.status}`))))
      .then(setData)
      .catch((e) => setError(String(e?.message ?? e)));
  }, []);

  const found = useMemo(() => {
    if (!data) return { rows: [] as number[], total: 0 };
    const q = query.trim().toLowerCase();
    if (!q) return { rows: data.ids.slice(0, ROWS_SHOWN).map((_, i) => i), total: data.ids.length };
    const rows: number[] = [];
    let total = 0;
    for (let i = 0; i < data.ids.length; i += 1) {
      if (!haystack(data, i).includes(q)) continue;
      total += 1;
      if (rows.length < ROWS_SHOWN) rows.push(i);
    }
    return { rows, total };
  }, [data, query]);

  if (error) return <p className="lookup-state">Could not load the variant index: {error}</p>;
  if (!data) return <p className="lookup-state">Loading all 70,376 variants…</p>;

  const i = found.rows.includes(picked) ? picked : (found.rows[0] ?? 0);
  const S = data.scales;
  const gene = data.genes[data.gene[i]];
  const stop = STOPS[data.stop[i]];
  const residue = data.residues[data.residue[i]];
  const conditionName = data.condition[i] === ABSENT ? null : data.conditions[data.condition[i]];
  const alsoIn = data.extra_conditions[i];
  const reach = data.reach[i] === ABSENT ? null : data.reach_classes[data.reach[i]];
  const escapes = data.escapes[i] === 1;
  const escapesFull = data.escapes_full[i];
  const best = data.best[i] / S.readthrough;
  const bestLow = data.best_low[i] / S.readthrough;
  const tolerance = data.tolerance[i] / S.tolerance;
  const width = data.therapies.length;
  const arms: Arm[] = data.therapies.map((id, slot) => ({
    name: THERAPY_NAME[id] ?? id,
    mid: data.arm_mid[i * width + slot] / S.readthrough,
    low: data.arm_low[i * width + slot] / S.readthrough,
    high: data.arm_high[i * width + slot] / S.readthrough,
  }));
  const editable = Boolean(reach) && reach !== "not_base_editable_under_panel";

  return (
    <div className="lookup">
      <aside className="lookup-list">
        <label className="lookup-search">
          <span className="sr-only">Search variants</span>
          <input
            type="search"
            value={query}
            onChange={(e) => setQuery(e.target.value)}
            placeholder="Gene, disease or stop codon…"
          />
        </label>
        <p className="lookup-count">
          {found.total.toLocaleString("en-US")} match{found.total === 1 ? "" : "es"}
          {found.total > found.rows.length ? ` · first ${found.rows.length}` : ""}
        </p>
        <ul>
          {found.rows.map((row) => (
            <li key={data.ids[row]}>
              <button
                type="button"
                className={row === i ? "lookup-hit active" : "lookup-hit"}
                onClick={() => setPicked(row)}
                aria-current={row === i ? "true" : undefined}
              >
                <b>{data.genes[data.gene[row]]}</b>
                <span>
                  {STOPS[data.stop[row]]} · {data.residues[data.residue[row]]} · codon {data.position[row]}
                </span>
                <span className="lookup-cond">
                  {data.condition[row] === ABSENT
                    ? "condition not recorded"
                    : data.conditions[data.condition[row]]}
                </span>
              </button>
            </li>
          ))}
          {found.rows.length === 0 && <li className="lookup-none">Nothing matches “{query}”.</li>}
        </ul>
      </aside>

      <div className="lookup-main">
        <header className="lookup-head">
          <div>
            <p className="route-kicker">Selected variant</p>
            <h2>
              {gene} <em>{stop}</em>
            </h2>
            <p className="lookup-condition">
              {conditionName ?? "Condition not recorded in the archive"}
              {alsoIn > 0 ? ` · and ${alsoIn} other condition${alsoIn === 1 ? "" : "s"}` : ""}
            </p>
            <p className="lookup-sub">
              {residue} at codon {data.position[i]} · {data.ids[i]}
            </p>
          </div>
          <span className="tag absent">no composite score</span>
        </header>

        <section className="lookup-block">
          <h3 className="block-title">Four conditions, kept apart</h3>
          <div className="layers">
            <Layer
              label="1 · Message survives"
              verdict={escapes ? "Survives" : "Destroyed"}
              detail={escapes ? "mRNA survives to be read through." : "Nothing left for a drug to act on."}
              value={escapes ? 1 : 0}
              tone="held"
            />
            <Layer
              label="2 · Ribosome reads through"
              verdict={pct(best)}
              detail={`Lower bound ${pct(bestLow)}.`}
              value={Math.min(1, best / 0.05)}
              tone="predicted"
            />
            <Layer
              label="3 · Inserted residue tolerated"
              verdict={pct(tolerance, 0)}
              detail="Share of insertable amino acids compatible with the original."
              value={tolerance}
              tone="held"
            />
            <Layer
              label="4 · Protein still works"
              verdict="Not modelled"
              detail="A planned layer, shown as a gap."
              value={null}
              tone="absent"
              absent
            />
          </div>
        </section>

        <section className="lookup-block">
          <h3 className="block-title">Which routes are open</h3>
          <div className="routes">
            <article className="route">
              <p className="route-kicker">Route one · small molecule</p>
              <h3>{pct(best)} at best</h3>
              <p>Lower bound {pct(bestLow)}. Intervals overlap, so this is a magnitude, not a choice.</p>
              <p className="route-verdict verdict-warn">Predicted, not separable</p>
            </article>
            <article className="route">
              <p className="route-kicker">Route two · suppressor tRNA</p>
              <h3>
                {stop}-{residue}
              </h3>
              <p>
                A design reading {stop} and inserting {residue} restores the original residue exactly.
              </p>
              <p className="route-verdict verdict-warn">Designable, unvalidated</p>
            </article>
            <article className="route">
              <p className="route-kicker">Route three · base editing</p>
              <h3>{editable ? "Reachable" : "Not reachable"}</h3>
              <p>{reach ? (REACH_WORD[reach] ?? reach) : "Not scored for this variant."}</p>
              <p className={`route-verdict ${editable ? "verdict-good" : "verdict-absent"}`}>
                {editable ? "Decidable — geometry fits" : "Decidable — geometry does not fit"}
              </p>
            </article>
          </div>
        </section>

        <section className="lookup-block">
          <h3 className="block-title">All six compounds, with their uncertainty</h3>
          <Intervals arms={arms} />
        </section>

        <section className="lookup-block">
          <h3 className="block-title">Does the message survive? Asked two ways</h3>
          <ul className="nmd-rows">
            <li className={escapes ? "ok" : "no"}>
              <span className="nmd-mark" aria-hidden="true" />
              <span className="nmd-label">50-nucleotide guideline</span>
              <span className="nmd-verdict">{escapes ? "survives" : "destroyed"}</span>
              <span className="nmd-note">Last exon, or near the last junction.</span>
            </li>
            {escapesFull !== ABSENT && (
              <li className={escapesFull === 1 ? "ok" : "no"}>
                <span className="nmd-mark" aria-hidden="true" />
                <span className="nmd-label">Fuller rule set</span>
                <span className="nmd-verdict">{escapesFull === 1 ? "survives" : "destroyed"}</span>
                <span className="nmd-note">Adds start-proximal and long-exon exceptions.</span>
              </li>
            )}
          </ul>
          {escapesFull !== ABSENT && (
            <p className={`verdict ${escapesFull === (escapes ? 1 : 0) ? "pass" : "caveat"}`}>
              {escapesFull === (escapes ? 1 : 0)
                ? "Both rule sets agree here."
                : "The rule sets disagree. Both verdicts are kept rather than one being picked."}
            </p>
          )}
        </section>
      </div>
    </div>
  );
}
