# Fuel-gauge calibration for a trawler tank

**Question: the gauge reads halfway up the main tank body. How many litres
remain, and how does the answer change when the boat heels?**

![Tank section and computed upright/rolled fuel calibration](tank.svg)

This example reconstructs a simplified aft tank from the FAO drawing
[12.8 m Ferrocement Trawler — Fuel and Water Tanks, IND 101A / 008](https://www.fao.org/fishery/services/storage/fs/fishery/images/vesseldesign/IND-101/IND-101_45.pdf)
(title block: Norwich, October 1984; revision dated 4 September 1985).
The drawing is the source of the dimensions, not a ready-made CAD dataset.
The original drawing is linked rather than bundled; the inputs and figure
are our reconstruction. No endorsement by FAO is implied.

## What is modeled

The tank body is convex. Its small sump projects below part of the floor,
making the combined interior **nonconvex**. The input is the union of two
convex polytopes sharing the sump opening. This is a small application
example, not a stress test for heavily overlapping unions.

Coordinates are **decimetres**, so the 3D volume is directly **litres**.
`x` measures distance from the inboard wall, `y` runs along the tank, and
`z=0` is the main floor. The sump extends to `z=-0.75`.

| Quantity read from the aft-tank views | Dimension used |
|---|---:|
| Length | 720 mm |
| Height above main floor | 500 + 344 = 844 mm |
| Bottom width at the two ends | 550 / 450 mm |
| Top width at the two ends | 850 / 750 mm |
| Sump size | 100 mm across × 150 mm along × 75 mm deep |
| Sump centre from inboard wall / first end | 150 / 175 mm |

The widths interpolate linearly along the length and height. In decimetres,
the body is

```text
0 ≤ y ≤ 36/5
0 ≤ z ≤ 211/25
0 ≤ x ≤ 11/2 − 5y/36 + 75z/211
```

and the sump is

```text
1 ≤ x ≤ 2,   1 ≤ y ≤ 5/2,   −3/4 ≤ z ≤ 0.
```

For an explicit nonconvexity witness, `(5, 1.75, 0.1)` and
`(1.5, 1.75, -0.7)` are inside the union, while their midpoint
`(3.25, 1.75, -0.3)` is outside it.

**Simplifications:** the listed dimensions are treated as interior dimensions;
the drawing's 3 mm plate thickness, 2 mm baffle, piping, sockets and fittings
are omitted. The sump is an open rectangular box. Its orientation and centre
positions are interpreted from the orthographic views. These choices define
the reproducible model; its computed capacity is not a certified tank capacity
and is not fitted to the drawing's nominal capacity notes.

Liquid is assumed settled under gravity in communicating, vented space.
The model computes geometric volume below a plane; it does not simulate
slosh, trapped air, or spilling through vents. Gauge placement at
`(x,y)=(1.5,1.75)` is an example assumption, not a sender shown in the drawing.

## Compute the fuel quantity

```bash
cargo build --release
./target/release/nefvol sweep --input inputs/fao_tank3d.ine \
  --direction 0,0,1 --at -0.75 --at 0 --at 4.22 --at 8.44
```

The upright sweep has the following independent geometric check, where
`h` is height above the main floor in decimetres:

```text
V(h) = 0                              h ≤ −3/4
       (3/2)(h + 3/4)                −3/4 ≤ h ≤ 0
       9/8 + 36h + (270/211)h²        0 ≤ h ≤ 211/25
       396117/1000                    h ≥ 211/25
```

The body holds 394.992 L and the sump 1.125 L, totaling **396.117 L**.
At `h=4.22 dm` (42.2 cm, halfway up the main body), it contains
**175.833 L**, or **44.39%** of capacity.

For a roll with slope `1/4` (about 14.04°), use direction `(1/4,0,1)`.
The plane is `x/4 + z = λ`. A gauge reading `h` at `x=1.5` means
**`λ = h + 0.375`**, not `λ=h`:

```bash
./target/release/nefvol sweep --input inputs/fao_tank3d.ine \
  --direction 1/4,0,1 --at 4.595 --json
```

At the same gauge height, this gives `72199101/459500` L ≈ **157.125 L**,
or **39.67%**. Ignoring this roll would overestimate the remaining fuel by
about **18.71 L**. This compares two states with the same gauge reading;
tilting a sealed quantity of fuel changes its level, not its volume.

Each direction produces the entire calibration function in one sweep.
Subsequent level readings only require polynomial evaluation. A new
orientation requires another sweep. Upright integration is simple here;
the tilted curves demonstrate why arbitrary plane directions are useful.

## Explore the shape in the existing viewer

```bash
python3 scripts/viewer.py
```

Open http://127.0.0.1:8765 and load
[`inputs/fao_tank2d.ine`](../inputs/fao_tank2d.ine). Use direction `0,1`
and initial λ `4.22`, then Compute. For the rolled section use `1/4,1`
and initial λ `4.595`. Dragging λ moves the liquid boundary.

This input is the exact `y=1.75 dm` section through the sump centre of
[`inputs/fao_tank3d.ine`](../inputs/fao_tank3d.ine). Its viewer coordinates
are `(x,z)`; the generic viewer labels them `(x,y)`. **The viewer reports
section area in dm², not fuel in litres.** The body tapers along its length
and the sump is shorter than the tank, so multiplying this section's area
by tank length does not give the volume. Use the 3D CLI for litres; the
existing browser renderer supports only 2D.

## Reproduce the figure and checks

The figure uses section vertices from `view-data` and volume curves from
two 3D `sweep --json` calls. Polynomial evaluation uses Python fractions;
only drawing coordinates are converted to floats.

```bash
python3 -m venv /tmp/nefvol-tank-plot
/tmp/nefvol-tank-plot/bin/pip install matplotlib
/tmp/nefvol-tank-plot/bin/python scripts/tank_figure.py
cargo test --release --test tank
```

Matplotlib is optional and only needed to regenerate `docs/tank.svg`.
The tests check the upright curve against integrated horizontal areas,
prove nonconvexity, check the 2D section against the 3D halfspaces, and
compare three tilted sweeps against the brute-force oracle (including
total-volume invariance).
