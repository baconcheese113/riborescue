"use client";

import {
  createColumnHelper,
  flexRender,
  getCoreRowModel,
  getFilteredRowModel,
  getSortedRowModel,
  useReactTable,
  type SortingState,
} from "@tanstack/react-table";
import { Fragment, useEffect, useMemo, useState } from "react";
import type { SafetyAtlas, Variant, WebTable } from "./types";

// The safety atlas is the measured layer beside the predicted one, kept visibly apart: what was
// measured, for which therapy and cell line, over how many stops, and how weakly the model tracks it.
function SafetyPanel({ safety }: { safety: SafetyAtlas }) {
  const { rho, low, high } = safety.concordance;
  const groups: [string, string][] = [
    ["both", "model and measurement agree it opens"],
    ["predicted only", "model flags it, measurement does not confirm"],
    ["measured only", "opens measurably, model missed it"],
    ["neither", "neither calls it open"],
  ];
  return (
    <details className="safety">
      <summary>
        <span className="tag measured">Native-stop atlas · {safety.measured_therapy} in{" "}
          {safety.cell_line}</span>
        <span className="summary-note">
          measured over {safety.analysed.toLocaleString()} canonical stops of{" "}
          {safety.canonical_stops_scored.toLocaleString()} scored
        </span>
      </summary>
      <div className="safety-body">
        <p className="unavailable-note">
          Only {safety.measured_therapy} has a matched empirical atlas; {safety.analysed.toLocaleString()}{" "}
          of {safety.canonical_stops_scored.toLocaleString()} canonical stops were analysed, because{" "}
          {safety.unavailable_reason}.
        </p>
        <p>
          The model&apos;s prediction at canonical stops is a separate, out-of-distribution layer. Over
          the measured genes it shows {safety.concordance_label} — Spearman {rho} [{low}, {high}]. The
          two are never merged into one score.
        </p>
        <div className="quadrants">
          {groups.map(([key, label]) => (
            <div className="quad" key={key}>
              <div className="quad-n">{(safety.quadrants[key] ?? 0).toLocaleString()}</div>
              <div className="quad-key">{key}</div>
              <div className="quad-label">{label}</div>
            </div>
          ))}
        </div>
        <p className="caveat">{safety.caveat}</p>
      </div>
    </details>
  );
}

const percent = (value: number | null) =>
  value === null ? "—" : `${(value * 100).toFixed(2)}%`;

const column = createColumnHelper<Variant>();

// The largest readthrough interval upper bound decides the shared axis every therapy bar is drawn
// on, so two therapies for one variant are comparable at a glance.
function TherapyBars({ variant }: { variant: Variant }) {
  const scale = Math.max(
    0.001,
    ...variant.therapies.map((t) => t.high ?? 0),
  );
  return (
    <div className="detail">
      <h3>Every therapy, with its interval</h3>
      {variant.therapies.map((t) =>
        t.available ? (
          <div className="therapy" key={t.id}>
            <span className="name">{t.id}</span>
            <span className="num">{percent(t.readthrough)}</span>
            <span className="track">
              <span
                className="ci"
                style={{
                  left: `${((t.low ?? 0) / scale) * 100}%`,
                  width: `${(((t.high ?? 0) - (t.low ?? 0)) / scale) * 100}%`,
                }}
              />
              <span
                className="point"
                style={{ left: `${((t.readthrough ?? 0) / scale) * 100}%` }}
              />
            </span>
          </div>
        ) : (
          <div className="therapy" key={t.id}>
            <span className="name">{t.id}</span>
            <span className="unavailable-therapy" style={{ gridColumn: "2 / 4" }}>
              no prediction — {t.reason ?? "unavailable"}
            </span>
          </div>
        ),
      )}
      <div className="suppressor">
        {variant.suppressor ? (
          <>
            Suppressor tRNA <code>{variant.suppressor.design}</code> reinserts the original{" "}
            {variant.residue} — restoring the protein exactly.
          </>
        ) : (
          <>No suppressor design maps to this stop.</>
        )}
      </div>
    </div>
  );
}

export default function Page() {
  const [data, setData] = useState<WebTable | null>(null);
  const [error, setError] = useState<string | null>(null);
  const [sorting, setSorting] = useState<SortingState>([
    { id: "best", desc: true },
  ]);
  const [query, setQuery] = useState("");
  const [decayOnly, setDecayOnly] = useState(false);
  const [expanded, setExpanded] = useState<string | null>(null);

  useEffect(() => {
    // Same-origin file written by `riborescue export-web`. A leading slash so it resolves from the
    // site root whatever route the viewer is on.
    fetch("/riborescue.json")
      .then((r) => {
        if (!r.ok) throw new Error(`riborescue.json returned ${r.status}`);
        return r.json();
      })
      .then(setData)
      .catch((e) => setError(String(e?.message ?? e)));
  }, []);

  const columns = useMemo(
    () => [
      column.accessor("gene", { header: "Gene" }),
      column.accessor("protein_position", {
        header: "Codon",
        cell: (c) => c.getValue() ?? "—",
      }),
      column.accessor("stop", { header: "Stop" }),
      column.accessor("residue", { header: "Residue" }),
      column.accessor((v) => v.best.therapy, { id: "therapy", header: "Best therapy" }),
      column.accessor((v) => v.best.readthrough, {
        id: "best",
        header: "Readthrough",
        cell: (c) => (
          <span>
            {percent(c.getValue())}
            <span style={{ color: "var(--faint)" }}>
              {" "}
              ≥ {percent(c.row.original.best.low)}
            </span>
          </span>
        ),
      }),
      column.accessor("escapes_decay", {
        header: "Escapes decay",
        cell: (c) => (
          <span className={`pill ${c.getValue() ? "good" : "bad"}`}>
            {c.getValue() ? "yes" : "no"}
          </span>
        ),
      }),
      column.accessor("tolerable_insertion_share", {
        header: "Tolerable",
        cell: (c) => percent(c.getValue()),
      }),
    ],
    [],
  );

  const rows = useMemo(() => {
    if (!data) return [];
    return decayOnly ? data.variants.filter((v) => v.escapes_decay) : data.variants;
  }, [data, decayOnly]);

  const table = useReactTable({
    data: rows,
    columns,
    state: { sorting, globalFilter: query },
    onSortingChange: setSorting,
    onGlobalFilterChange: setQuery,
    getCoreRowModel: getCoreRowModel(),
    getSortedRowModel: getSortedRowModel(),
    getFilteredRowModel: getFilteredRowModel(),
  });

  if (error)
    return (
      <div className="wrap">
        <h1>Variant × therapy</h1>
        <div className="banner" style={{ borderLeftColor: "var(--bad)" }}>
          <span className="tag" style={{ color: "var(--bad)" }}>Could not load data</span>
          <p>
            {error}. Generate it with <code>pixi run export-web-example</code>, which writes{" "}
            <code>frontend/public/riborescue.json</code>, then reload.
          </p>
        </div>
      </div>
    );
  if (!data) return <div className="wrap">Loading…</div>;

  return (
    <div className="wrap">
      <h1>Variant × therapy</h1>
      <p className="lede">
        Each pathogenic nonsense variant and the readthrough therapy it looks most amenable to.
        Predictions are the amenability model&apos;s; they are not a claim any therapy works.
        This is a <b>{data.variants.length}-variant example</b> — type a gene to filter, click a row
        to see every therapy&apos;s interval and the suppressor tRNA.
      </p>

      <div className="banner">
        <span className="tag">Validation · {data.status.readthrough_control}</span>
        <p>{data.status.readthrough_detail}</p>
      </div>
      <div className="banner unavailable">
        <span className="tag">Safety · {data.status.safety_atlas}</span>
        <p>{data.status.safety_detail}</p>
      </div>

      {data.safety && <SafetyPanel safety={data.safety} />}

      <div className="controls">
        <input
          type="search"
          placeholder="Filter by gene, stop, residue…"
          value={query}
          onChange={(e) => setQuery(e.target.value)}
        />
        <label className="check">
          <input
            type="checkbox"
            checked={decayOnly}
            onChange={(e) => setDecayOnly(e.target.checked)}
          />
          escapes decay only
        </label>
      </div>

      <div className="scroller">
        <table>
          <thead>
            {table.getHeaderGroups().map((group) => (
              <tr key={group.id}>
                {group.headers.map((header) => (
                  <th key={header.id} onClick={header.column.getToggleSortingHandler()}>
                    {flexRender(header.column.columnDef.header, header.getContext())}
                    <span className="arrow">
                      {{ asc: " ↑", desc: " ↓" }[header.column.getIsSorted() as string] ?? ""}
                    </span>
                  </th>
                ))}
              </tr>
            ))}
          </thead>
          <tbody>
            {table.getRowModel().rows.map((row) => (
              <Fragment key={row.original.id}>
                <tr
                  className="parent"
                  onClick={() =>
                    setExpanded(expanded === row.original.id ? null : row.original.id)
                  }
                >
                  {row.getVisibleCells().map((cell) => (
                    <td
                      key={cell.id}
                      className={
                        typeof cell.getValue() === "number" ? "num" : undefined
                      }
                    >
                      {flexRender(cell.column.columnDef.cell, cell.getContext())}
                    </td>
                  ))}
                </tr>
                {expanded === row.original.id && (
                  <tr>
                    <td colSpan={columns.length} style={{ padding: 0 }}>
                      <TherapyBars variant={row.original} />
                    </td>
                  </tr>
                )}
              </Fragment>
            ))}
          </tbody>
        </table>
      </div>
      <p className="count">
        {table.getRowModel().rows.length} of {data.variants.length} variants shown
      </p>
    </div>
  );
}
