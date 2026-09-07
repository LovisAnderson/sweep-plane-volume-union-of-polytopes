#!/usr/bin/env python3
"""Rebuild docs/tank.svg from the example inputs and nefvol's exact sweeps.

Requires matplotlib (only for this optional figure, not the solver/viewer).
Run from any directory after cargo build --release.
"""
import argparse
from fractions import Fraction as Q
import json
from pathlib import Path
import subprocess

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon

ROOT = Path(__file__).resolve().parent.parent


def query(binary, command, filename, direction):
    args = [str(binary), command, "--input", str(ROOT / "inputs" / filename),
            f"--direction={direction}", "--threads", "2"]
    if command == "sweep":
        args.append("--json")
    return json.loads(subprocess.check_output(args, text=True))


def evaluate(data, level):
    for piece in data["pieces"]:
        if ((piece["lo"] is None or Q(piece["lo"]) <= level)
                and (piece["hi"] is None or level < Q(piece["hi"]))):
            value = Q(0)
            for coefficient in reversed(piece["coefficients"]):
                value = value * level + Q(coefficient)
            return value
    raise ValueError("Level outside sweep pieces")


def clip_below(points, level):
    clipped = []
    for p, r in zip(points, points[1:] + points[:1]):
        if p[1] <= level:
            clipped.append(p)
        if (p[1] < level < r[1]) or (r[1] < level < p[1]):
            t = (level - p[1]) / (r[1] - p[1])
            clipped.append([p[0] + t * (r[0] - p[0]), level])
    return clipped


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, default=ROOT / "target/release/nefvol")
    parser.add_argument("--output", type=Path, default=ROOT / "docs/tank.svg", help="Output SVG path")
    args = parser.parse_args()
    section = query(args.binary.resolve(), "view-data", "fao_tank2d.ine", "0,1")
    upright = query(args.binary.resolve(), "sweep", "fao_tank3d.ine", "0,0,1")
    tilted = query(args.binary.resolve(), "sweep", "fao_tank3d.ine", "1/4,0,1")
    half = Q(211, 50)
    # Gauge at x=1.5 dm, y=1.75 dm. At roll slope 1/4, lambda = h + x/4.
    offset = Q(3, 8)
    v0, v1 = evaluate(upright, half), evaluate(tilted, half + offset)
    total = Q(upright["total"])
    assert Q(tilted["total"]) == total

    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
                         "svg.fonttype": "none", "svg.hashsalt": "nefvol-tank"})
    fig, (shape, graph) = plt.subplots(1, 2, figsize=(12, 6),
                                      gridspec_kw={"width_ratios": [1, 1.45]})
    fig.patch.set_facecolor("#f5f7fa")
    fig.suptitle("Halfway up the tank. How much fuel is left?", fontsize=19, x=.06, ha="left")
    fig.text(.06, .89, "FAO trawler aft tank · simplified interior · capacity %.3f L" % total,
             color="#536273")
    for poly in section["polytopes"]:
        points = [[float(Q(c)) * 10 for c in p] for p in poly["vertices"]]
        shape.add_patch(Polygon(points, facecolor="#e0e7ef", edgecolor="#718096"))
        wet = clip_below(points, float(half) * 10)
        if len(wet) >= 3:
            shape.add_patch(Polygon(wet, facecolor="#63b9c4", edgecolor="none"))
    shape.axhline(float(half) * 10, color="#b7522c", lw=1.5, ls="--")
    shape.plot([15], [float(half) * 10], "o", color="#b7522c")
    shape.annotate("Gauge: 42.2 cm", (15, 42.2), (25, 32),
                   arrowprops={"arrowstyle": "-", "color": "#b7522c"})
    shape.annotate("Sump makes the union nonconvex", (15, -4), (0, -24),
                   fontsize=9, arrowprops={"arrowstyle": "->", "color": "#536273"})
    shape.set(xlim=(-5, 90), ylim=(-28, 92), xlabel="Distance from inboard wall (cm)",
              ylabel="Height above main floor (cm)", title="Section through the sump · upright")
    shape.set_aspect("equal")

    heights = [Q(-3, 4) + (Q(211, 25) + Q(3, 4)) * Q(i, 400) for i in range(401)]
    for data, shift, color, label in [(upright, Q(0), "#167787", "Upright"),
                                      (tilted, offset, "#b7522c", "14.0° roll (slope 1/4)")]:
        graph.plot([float(h) * 10 for h in heights],
                   [float(evaluate(data, h + shift)) for h in heights],
                   color=color, lw=2.5, label=label)
    graph.axvline(42.2, color="#8493a5", ls="--", lw=1)
    for value, color in [(v0, "#167787"), (v1, "#b7522c")]:
        graph.plot([42.2], [float(value)], "o", color=color)
    graph.text(.05, .95, "Same gauge height: 42.2 cm\n"
               f"Upright: {float(v0):.2f} L ({float(v0 / total):.1%})\n"
               f"Rolled:  {float(v1):.2f} L ({float(v1 / total):.1%})",
               transform=graph.transAxes, va="top", linespacing=1.7,
               bbox={"facecolor": "white", "edgecolor": "#dce3ea", "pad": 10})
    graph.set(xlabel="Gauge height above main floor (cm)", ylabel="Fuel in the 3D tank (litres)",
              title="Calibration from the complete 3D sweep", ylim=(0, 420))
    graph.legend(loc="lower right", frameon=False)
    graph.grid(alpha=.2)
    for ax in (shape, graph):
        ax.spines[["top", "right"]].set_visible(False)
        ax.set_facecolor("#f5f7fa")
    fig.text(.06, .035, "Section shading shows area; the calibration uses the full 3D tank. "
             "Dimensions idealized; static liquid, no slosh or overflow model.", fontsize=9, color="#536273")
    fig.subplots_adjust(top=.79, bottom=.18, left=.07, right=.97, wspace=.3)
    fig.savefig(args.output, format="svg", metadata={"Date": None})
    # Matplotlib emits trailing spaces in path coordinates; keep the artifact diff-clean.
    args.output.write_text("\n".join(line.rstrip() for line in args.output.read_text().splitlines()) + "\n")
    print(f"Capacity: {total} L = {float(total):.3f} L")
    print(f"At h={half} dm: upright {v0} L; rolled {v1} L")
    print(f"Figure: {args.output}")


if __name__ == "__main__":
    main()
