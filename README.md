# nefvol — exact sweep-plane volume of a union of polytopes

`nefvol` computes, for a union `U = P_1 ∪ … ∪ P_k` of rational polytopes
given by inequalities and a sweep direction `a`,

    s(a, λ) = vol( U ∩ { x : ⟨a, x⟩ ≤ λ } )

as an **exact** piecewise polynomial in `λ`, plus point evaluations and the
total volume `s(a, +∞)`.  Polytopes may overlap, nest, touch or be disjoint.

It is one algorithm made of two halves that share a single per-vertex object:

* **Reverse search** (Avis–Fukuda) walks the vertices of the hyperplane
  arrangement on `∂U`, one facet at a time, each vertex exactly once, with no
  visited set — memory is the current root path only.
* **Sweep plane** (Bieri–Nef) turns the tangent cone at each such vertex into
  signed simplicial terms `α (λ − ⟨a,v⟩)_+^d` that accumulate into a knot map.

Both are driven by the central arrangement of the hyperplanes through the
vertex: its lines are the search's edges, its bases are the sweep's cones.
See [`docs/DESIGN.md`](docs/DESIGN.md) for the mathematics (in particular the
symbolic-perturbation form of the tangent-cone decomposition that makes the
kernel a single factorisation per basis) and
[`docs/REPORT.md`](docs/REPORT.md) for measurements.  The original brief is
[`SPEC.md`](SPEC.md).

## Build

Rust 1.75+ (stable).  No GMP: arithmetic is `i128` with automatic promotion
to `num-bigint`.

    cargo build --release
    cargo test --release          # unit, golden, differential, Monte Carlo, recorded lrs values

## Usage

    nefvol volume --input <file> [--input <file>…] [--direction a1,a2,…]
    nefvol sweep  --input <file> --direction a1,a2,… [--at λ]… [--json]
    nefvol plot   --input <file> --direction a1,a2,… --samples N

Common options: `--threads N`, `--budget N` (nodes per subtree task, 0 = no
splitting), `--assume-generic` (fail instead of perturbing when `a` is
orthogonal to an edge), `--oracle` (brute-force enumerator, for checking),
`--seed S`, `-v`, `--stats-json` (run statistics on stderr).

Example:

    $ nefvol sweep --input inputs/square.ine --direction 1,2 --at 1
    direction a = (1, 2)
    breakpoints (4):
      0
      1
      2
      3
    pieces (coefficients of λ^0..λ^2):
      [-inf, 0): 0
      [0, 1): 1/4·λ^2
      [1, 2): -1/4 + 1/2·λ
      [2, 3): -5/4 + 3/2·λ + -1/4·λ^2
      [3, +inf): 1
    total volume: 1
    s(1) = 1/4

### Input formats

* **lrs/cdd `.ine`** H-representation.  Each row `b a_1 … a_d` means
  `b + Σ a_i x_i ≥ 0`.  Several `begin … end` blocks in one file, or several
  `--input` files, form a union.  `linearity` rows are rejected (polytopes
  must be full-dimensional).
* **JSON**: `{"polytopes": [{"A": [[…]], "b": […]}, …]}` meaning `A x ≤ b`;
  entries may be numbers or strings such as `"7/2"`.

Each polytope must be bounded and full-dimensional; the preprocessing checks
this with exact LPs and exits with a message otherwise.  Redundant
constraints are removed automatically.

## Layout

| file | role |
|---|---|
| `src/num.rs` | hybrid integer `Z` (i128 → BigInt), rationals `Q`, points |
| `src/linalg.rs` | fraction-free `(det, adj)`, incremental echelon form |
| `src/geom.rs` | canonical hyperplanes, polytopes, shared arrangement, fast `i64` mirror |
| `src/local.rs` | **the fused kernel**: lines and signed cones at a vertex, Laurent terms |
| `src/walk.rs` | tangent-cone filter, parent rule, ratio test, exact simplex walk |
| `src/reverse.rs` | per-facet reverse search with budgeted subtree tasks |
| `src/vertex.rs` | `∂U` test, ownership, per-vertex accumulation |
| `src/sweep.rs` | knot map, pieces, evaluation, pole-cancellation check |
| `src/preprocess.rs` | vertex finding, boundedness, dimension, redundancy |
| `src/oracle.rs` | brute-force `d`-subset enumerator (test oracle) |
| `src/driver.rs` | pipeline and parallel execution |
| `src/io.rs`, `src/main.rs` | formats and CLI |
| `scripts/gen.py` | benchmark instance generator |
| `scripts/lrs_check.py` | optional: regenerates `tests/lrs_golden.txt` from an `lrs` binary |
| `scripts/bench.sh` | reproduces the report tables |

## Non-goals

Full Nef generality (complements, open faces), unbounded polyhedra,
floating-point fast paths.  The tangent-cone representation (signed cells per
basis) is local to `local.rs`, so Nef generality can be added there without
touching the accumulator.

## Interactive browser viewer (2D)

    cargo build --release
    python3 scripts/viewer.py

Open **http://127.0.0.1:8765**, choose one or more `.ine` or JSON files, enter
an unnormalized sweep direction (for example `1,2`) and optionally an initial
λ, then click **Compute**. The union and swept portion appear beside the
cumulative area graph. Drag the slider or enter a number/fraction in the λ
field to move both markers. λ is clamped to the minimum and maximum vertex
projections. Hover over the area readout for its exact rational value.

Each Compute invokes the existing algorithm once and caches the complete
piecewise polynomial in the browser. Slider movement makes no server requests;
polynomial evaluation uses exact rational arithmetic, with floating-point
conversion only for display. Changing files or direction requires Compute
again. Invalid, unbounded, lower-dimensional, and non-2D inputs show errors.

For a more complex example, load [`inputs/overlap2d_50.ine`](inputs/overlap2d_50.ine):
50 overlapping convex polygons formed by clipping random boxes with four
additional halfspaces. They form one connected union, with 289 pairs sharing
positive area. Try direction `1,2` and move λ through the full range.
Regenerate the same input with:

    python3 scripts/gen.py rand 2 4 50 42 > inputs/overlap2d_50.ine

Python 3 and a modern browser are required; no npm install or external web
assets are needed. The server binds only to localhost. Use `--port` to select
another port or `--binary` to select a nefvol executable. Stop with Ctrl+C.
Requests are limited to 8 MiB and computations to 120 seconds.

The `nefvol view-data --input … --direction 1,2` command returns presentation
JSON: dimension, coordinate-array vertices per polytope, direction, λ range,
exact polynomial pieces, and total volume. Geometry extraction and the
`Renderer2D` interface are separate from solving, graphing, and controls. A 3D
extension needs polyhedron vertices/faces and a renderer (such as Three.js)
implementing `update(lambda)`; this version does not yet render 3D.

Viewer checks: `cargo test --release --test viewer` and
`node viewer/math.test.mjs`.
For the optional browser regression check, start the viewer, make Playwright
available to Node, and run `node viewer/browser.test.cjs`. `VIEWER_URL` overrides
the server URL; `PLAYWRIGHT_MODULE` can point to an external Playwright install.
The check exercises file loading, slider synchronization, cached requests,
fraction input, direction changes, and unsupported dimensions.
