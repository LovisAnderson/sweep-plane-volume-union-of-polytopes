import {
  rational,
  number,
  exact,
  compare,
  prepare,
  evaluate,
  interpolate,
} from "./math.js";
const $ = (id) => document.getElementById(id);
const NS = "http://www.w3.org/2000/svg";
function element(parent, tag, attrs = {}, text) {
  const e = document.createElementNS(NS, tag);
  for (const [k, v] of Object.entries(attrs)) e.setAttribute(k, v);
  if (text !== undefined) e.textContent = text;
  parent.append(e);
  return e;
}
const line = (s, x1, y1, x2, y2, attrs = {}) =>
  element(s, "line", { x1, y1, x2, y2, stroke: "#c9d5dd", ...attrs });
const fmt = (n) =>
  Number(n).toLocaleString(undefined, { maximumSignificantDigits: 6 });
function clip(points, a, lambda) {
  const out = [];
  for (let i = 0; i < points.length; i++) {
    const p = points[i],
      q = points[(i + 1) % points.length],
      u = a[0] * p[0] + a[1] * p[1] - lambda,
      v = a[0] * q[0] + a[1] * q[1] - lambda;
    if (u <= 0) out.push(p);
    if ((u < 0 && v > 0) || (u > 0 && v < 0)) {
      const t = u / (u - v);
      out.push([p[0] + t * (q[0] - p[0]), p[1] + t * (q[1] - p[1])]);
    }
  }
  return out;
}
// Renderer contract: construct from scene data, then update(lambda).
// A future 3D renderer can implement the same contract without changing controls or solving.
class Renderer2D {
  constructor(data) {
    this.svg = $("geometry");
    this.svg.replaceChildren();
    this.a = data.direction.map((x) => number(rational(x)));
    this.polys = data.polytopes.map((p) =>
      p.vertices.map((v) => v.map((x) => number(rational(x)))),
    );
    const points = this.polys.flat(),
      xs = points.map((p) => p[0]),
      ys = points.map((p) => p[1]);
    const xmin = Math.min(...xs),
      xmax = Math.max(...xs),
      ymin = Math.min(...ys),
      ymax = Math.max(...ys);
    const scale = Math.min(450 / (xmax - xmin), 290 / (ymax - ymin));
    if (!Number.isFinite(scale) || scale <= 0 || !this.a.every(Number.isFinite))
      throw Error("Coordinates exceed the display range");
    this.xy = ([x, y]) => [
      280 + (x - (xmin + xmax) / 2) * scale,
      200 - (y - (ymin + ymax) / 2) * scale,
    ];
    this.box = [
      [xmin - 25 / scale, ymin - 25 / scale],
      [xmax + 25 / scale, ymin - 25 / scale],
      [xmax + 25 / scale, ymax + 25 / scale],
      [xmin - 25 / scale, ymax + 25 / scale],
    ];
    for (const p of this.polys) this.polygon(p, { fill: "#dce4e9" });
    this.fill = element(this.svg, "g");
    // Outlines are drawn faintly to reveal the input polytopes; fill has no overlap darkening.
    const outlines = element(this.svg, "g");
    for (const p of this.polys)
      this.polygon(
        p,
        { fill: "none", stroke: "#8499a6", "stroke-width": 1 },
        outlines,
      );
    this.sweep = line(this.svg, 0, 0, 0, 0, {
      stroke: "#d15527",
      "stroke-width": 2,
    });
    element(
      this.svg,
      "text",
      { x: 35, y: 385 },
      `x: ${fmt(xmin)} … ${fmt(xmax)}     y: ${fmt(ymin)} … ${fmt(ymax)}`,
    );
  }
  polygon(p, attrs, parent = this.svg) {
    return element(parent, "polygon", {
      points: p
        .map(this.xy)
        .map((v) => v.join(","))
        .join(" "),
      ...attrs,
    });
  }
  update(lambda) {
    this.fill.replaceChildren();
    for (const p of this.polys)
      this.polygon(clip(p, this.a, lambda), { fill: "#69b9ae" }, this.fill);
    const hits = [];
    for (let i = 0; i < 4; i++) {
      const p = this.box[i],
        q = this.box[(i + 1) % 4],
        u = this.a[0] * p[0] + this.a[1] * p[1] - lambda,
        v = this.a[0] * q[0] + this.a[1] * q[1] - lambda;
      if (u === 0) hits.push(p);
      if (u * v < 0) {
        const t = u / (u - v);
        hits.push([p[0] + t * (q[0] - p[0]), p[1] + t * (q[1] - p[1])]);
      }
    }
    this.sweep.setAttribute(
      "visibility",
      hits.length >= 2 ? "visible" : "hidden",
    );
    if (hits.length >= 2) {
      const [x1, y1] = this.xy(hits[0]),
        [x2, y2] = this.xy(hits[1]);
      for (const [k, v] of Object.entries({ x1, y1, x2, y2 }))
        this.sweep.setAttribute(k, v);
    }
  }
}
class VolumeGraph {
  constructor(data, pieces, lo, hi) {
    this.svg = $("graph");
    this.svg.replaceChildren();
    const min = number(lo),
      max = number(hi),
      total = number(rational(data.total));
    if (![min, max, total].every(Number.isFinite) || max <= min || total <= 0)
      throw Error("Values exceed the display range");
    this.x = (v) => 65 + (460 * (v - min)) / (max - min);
    this.y = (v) => 345 - (305 * v) / total;
    for (let i = 0; i <= 4; i++) {
      const y = 345 - (305 * i) / 4;
      line(this.svg, 65, y, 525, y);
      element(
        this.svg,
        "text",
        { x: 57, y: y + 4, "text-anchor": "end" },
        fmt((total * i) / 4),
      );
    }
    line(this.svg, 65, 40, 65, 345);
    line(this.svg, 65, 345, 525, 345);
    element(this.svg, "text", { x: 65, y: 369 }, fmt(min));
    element(
      this.svg,
      "text",
      { x: 525, y: 369, "text-anchor": "end" },
      fmt(max),
    );
    element(this.svg, "text", { x: 295, y: 392 }, "λ");
    const samples = Array.from({ length: 501 }, (_, i) =>
      interpolate(lo, hi, i * 2),
    );
    for (const p of pieces)
      if (p.lo && compare(p.lo, lo) >= 0 && compare(p.lo, hi) <= 0)
        samples.push(p.lo);
    samples.sort(compare);
    element(this.svg, "path", {
      d: samples
        .map(
          (q, i) =>
            `${i ? "L" : "M"}${this.x(number(q))},${this.y(number(evaluate(pieces, q)))}`,
        )
        .join(" "),
      fill: "none",
      stroke: "#146c72",
      "stroke-width": 2.5,
    });
    this.mark = line(this.svg, 65, 40, 65, 345, {
      stroke: "#d15527",
      "stroke-dasharray": "5 4",
    });
    this.dot = element(this.svg, "circle", { r: 5, fill: "#d15527" });
  }
  update(lambda, area) {
    const x = this.x(lambda);
    this.mark.setAttribute("x1", x);
    this.mark.setAttribute("x2", x);
    this.dot.setAttribute("cx", x);
    this.dot.setAttribute("cy", this.y(area));
  }
}
let state;
function status(text, error = false) {
  $("status").textContent = text;
  $("status").className = error ? "error" : "";
}
function update(q) {
  if (!state) return;
  q =
    compare(q, state.lo) < 0
      ? state.lo
      : compare(q, state.hi) > 0
        ? state.hi
        : q;
  const value = evaluate(state.pieces, q),
    l = number(q);
  state.renderer.update(l);
  state.graph.update(l, number(value));
  $("lambda").value = exact(q);
  $("slider").value = Math.round(
    (1000 * (l - number(state.lo))) / (number(state.hi) - number(state.lo)),
  );
  $("area").textContent =
    `${fmt(number(value))} / ${fmt(number(rational(state.data.total)))}`;
  $("area").title = `${exact(value)} / ${state.data.total}`;
}
$("slider").addEventListener("input", () =>
  update(interpolate(state.lo, state.hi, Number($("slider").value))),
);
$("lambda").addEventListener("change", () => {
  try {
    update(rational($("lambda").value));
    status("Updated λ from the cached sweep function.");
  } catch (e) {
    status(e.message, true);
  }
});
function invalidate() {
  state = null;
  $("slider").disabled = $("lambda").disabled = true;
  status("Input changed. Compute to update the views.");
  $("geometry").replaceChildren();
  $("graph").replaceChildren();
  $("area").textContent = "—";
  $("area").title = "";
  $("minimum").textContent = $("maximum").textContent = "—";
  $("lambda").value = "";
}
$("files").addEventListener("change", invalidate);
$("direction").addEventListener("input", invalidate);
$("compute-form").addEventListener("submit", async (event) => {
  event.preventDefault();
  const controls = ["compute", "files", "direction", "initial"];
  try {
    const initial = $("initial").value.trim()
      ? rational($("initial").value)
      : null;
    const files = Array.from($("files").files);
    if (!files.length) throw Error("Choose at least one file");
    if (files.reduce((n, f) => n + f.size, 0) > 7 * 1024 * 1024)
      throw Error("Choose files totaling less than 7 MiB");
    state = null;
    $("slider").disabled = $("lambda").disabled = true;
    for (const id of controls) $(id).disabled = true;
    status("Computing the complete sweep function…");
    const response = await fetch("/api/compute", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({
        files: await Promise.all(files.map((f) => f.text())),
        direction: $("direction").value,
      }),
    });
    const data = await response.json();
    if (!response.ok) throw Error(data.error);
    if (data.dimension !== 2) throw Error("No renderer for this dimension");
    const [lo, hi] = data.lambdaRange.map(rational),
      pieces = prepare(data);
    const renderer = new Renderer2D(data),
      graph = new VolumeGraph(data, pieces, lo, hi);
    state = { data, lo, hi, pieces, renderer, graph };
    $("minimum").textContent = exact(lo);
    $("maximum").textContent = exact(hi);
    $("slider").disabled = $("lambda").disabled = false;
    update(initial || lo);
    status("Computed. Move λ freely; no further algorithm calls are needed.");
  } catch (e) {
    status(e.message, true);
  } finally {
    for (const id of controls) $(id).disabled = false;
  }
});
