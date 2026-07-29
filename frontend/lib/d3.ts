import { JSDOM } from "jsdom";
import * as d3 from "d3";

/** Colours are written into the SVG as CSS custom properties rather than as hex, so one rendered
 * figure serves both themes and switching theme repaints without redrawing anything. */
export const C = {
  ink: "var(--ink)",
  muted: "var(--muted)",
  faint: "var(--faint)",
  line: "var(--line)",
  grid: "var(--grid)",
  panel: "var(--panel-2)",
  measured: "var(--blue)",
  measuredSoft: "var(--blue-soft)",
  held: "var(--teal)",
  heldSoft: "var(--teal-soft)",
  alarm: "var(--orange)",
  alarmSoft: "var(--orange-soft)",
  caveat: "var(--amber)",
  caveatSoft: "var(--amber-soft)",
  predicted: "var(--purple)",
  predictedSoft: "var(--purple-soft)",
  absent: "var(--gray)",
} as const;

const { document } = new JSDOM("").window;

export type Svg = d3.Selection<SVGSVGElement, unknown, null, undefined>;

/** A figure is drawn once at build time and shipped as markup. `label` is the whole finding in a
 * sentence with its numbers, because that is what a screen reader gets instead of the picture. */
export function draw(
  width: number,
  height: number,
  label: string,
  body: (svg: Svg) => void,
): string {
  const host = d3.select(document.body).append("div");
  const svg = host
    .append("svg")
    .attr("viewBox", `0 0 ${width} ${height}`)
    .attr("role", "img")
    .attr("aria-label", label) as unknown as Svg;
  body(svg);
  const markup = host.html();
  host.remove();
  return markup;
}

/** Native tooltips: a `title` child costs nothing, needs no script, and survives a static export. */
export function tip<T extends d3.BaseType, D>(
  selection: d3.Selection<T, D, d3.BaseType, unknown>,
  text: (d: D) => string,
) {
  selection.append("title").text(text as (d: D, i: number) => string);
  return selection;
}

/** Gridlines are an axis drawn inward with its labels and baseline removed. */
export function gridX(svg: Svg, scale: d3.ScaleLinear<number, number>, top: number, bottom: number, ticks: number[]) {
  svg
    .append("g")
    .attr("class", "grid")
    .attr("transform", `translate(0,${bottom})`)
    .call(d3.axisBottom(scale).tickValues(ticks).tickSize(-(bottom - top)).tickFormat(() => ""))
    .call((g) => g.select(".domain").remove());
}

export function axisX(
  svg: Svg,
  scale: d3.ScaleLinear<number, number>,
  y: number,
  ticks: number[],
  format: (value: number) => string,
) {
  svg
    .append("g")
    .attr("class", "axis")
    .attr("transform", `translate(0,${y})`)
    .call(
      d3
        .axisBottom(scale)
        .tickValues(ticks)
        .tickFormat((value) => format(value as number))
        .tickSizeOuter(0),
    );
}

export function axisY(svg: Svg, scale: d3.ScaleBand<string>, x: number) {
  svg
    .append("g")
    .attr("class", "axis")
    .attr("transform", `translate(${x},0)`)
    .call(d3.axisLeft(scale).tickSizeOuter(0))
    .call((g) => g.select(".domain").remove());
}

/** A mark whose condition held is filled; one whose condition did not is hollow on the panel
 * colour with its stroke kept, so the eye reads outcome before it reads hue. */
export function outcome(held: boolean, colour: string) {
  return { fill: held ? colour : C.panel, stroke: colour, strokeWidth: 2 };
}

export const fmtInt = d3.format(",");
export const fmtPct1 = d3.format(".1%");
export const fmtPct2 = d3.format(".2%");
export const fmtFixed = (digits: number) => d3.format(`.${digits}f`);
export const fmtSigned = (digits: number) => d3.format(`+.${digits}f`);
