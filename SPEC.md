# Reverse search + sweep plane: exact volume of a union of polytopes, in Rust

## 0. The thesis — read this before anything else

This project is **one algorithm made of two halves that fuse**, not a volume
algorithm with an enumeration subroutine bolted on.

- **Reverse search** (Avis–Fukuda) enumerates the vertices of the hyperplane
  arrangement induced by the input, visiting each exactly once, with no visited
  set and `O(depth)` memory. To find a node's children it must compute the
  **local arrangement at that vertex**.
- **The sweep plane** (Bieri–Nef) computes `vol(U ∩ {⟨a,x⟩ ≤ λ})` from purely
  local data at each arrangement vertex: the **tangent cone**, decomposed into
  simplicial cones. It needs nothing else — no cells, no global vertex list.

They fuse because **both need exactly the same per-vertex object**: the central
arrangement of the hyperplanes through `v`. Its extreme rays are the reverse
search's edges; its cells are the sweep plane's tangent cone pieces. Compute it
once per vertex and feed both.

The result is a single streaming traversal:

    reverse search tree node  v
        │
        ├─ local arrangement of H_v ──┬──► extreme rays ──► children (recurse)
        │                             │
        │                             └──► cells ∩ U ──► simplicial cones
        │                                                    │
        └────────────────────────────────────────────────────┴──► knot map

Nothing is ever stored except the knot map `{ ⟨a,v⟩ : coefficient }`, which is
tiny. No vertex list, no cell list, no incidence structure, no dedup table.
Memory must stay flat as the vertex count grows. An implementation that stores
the vertices or the cells has missed the point and must be reworked.

Two independent correctness pillars the design must preserve: reverse search
gives *each vertex exactly once*, and Bieri–Nef gives *the exact contribution of
each vertex*. Test them separately (§7) before testing the composition.

## 1. Problem statement

Input: `U = P_1 ∪ … ∪ P_k`, each `P_i = { x ∈ Q^d : A_i x ≤ b_i }` a bounded,
full-dimensional rational polytope. They may overlap, nest, or be disjoint; `U`
is in general non-convex and disconnected. Let `H` be all facet hyperplanes of
all `P_i`, canonicalised to integer-primitive form and deduplicated.

Output: for a sweep direction `a ∈ Q^d`,

    s(a, λ) = vol( U ∩ { x : ⟨a, x⟩ ≤ λ } )

as an **exact** piecewise polynomial of degree `d` in `λ`, plus exact point
evaluation and the total volume `s(a, +∞)`.

References: Bieri & Nef, *A sweep-plane algorithm for computing the volume of
polyhedra represented in boolean form*, Linear Algebra Appl. 52:69–97 (1983);
Anderson & Hiller, ZIB-Report 18-37 / doi:10.1007/978-3-030-18500-8_12;
Avis & Fukuda, *A pivoting algorithm for convex hulls and vertex enumeration of
arrangements and polyhedra*, DCG 8:295–313 (1992).

## 2. The reverse search half

### 2.1 What is enumerated

The arrangement vertices lying on `∂U`. Only these contribute: if `v ∈ int(U)`
then `T_v(U) = R^d`, and contributions over any fan covering `R^d` sum to zero
because `R^d` contains a line and is therefore zero in the valuation (§3.1).
This restriction is what makes the enumeration output-sensitive, and therefore
what makes reverse search worth using instead of enumerating all `C(|H|, d)`
bases.

### 2.2 Reduce to `d−1` dimensions

`∂U ⊆ ⋃_i ∂P_i`, so run one search per facet. Fix a facet `F` of `P_i` lying in
hyperplane `h_0`. Inside `aff(h_0) ≅ R^{d−1}` take the induced hyperplanes
`G = { h ∩ h_0 : h ∈ H \ {h_0} }`, discarding non-hyperplane intersections.
Every hyperplane bounding `F` is a facet hyperplane of `P_i` and hence already
in `G`, so `F` is a union of cells of `A(G)` and the subproblem is
self-contained:

> enumerate the vertices of an arrangement of `|G|` hyperplanes in `R^{d−1}`,
> restricted to the bounded convex region `F`.

### 2.3 The local search function

Fix a generic linear objective `w` on `aff(h_0)`, symbolically perturbed so no
two vertices tie and no edge is `w`-flat. Define `parent(v)` as the neighbour
along the canonically-first `w`-improving edge at `v` (Bland-style tie-breaking
on a fixed ordering of edge directions, so it is deterministic and computable
from `v` alone, with no memory).

**Unique sink, with proof.** Let `v` be a vertex of the induced complex that is
not `w`-maximal in `F`. `F` is convex and bounded, so a direction of increase
exists at `v`; it lies in some cell `C` of the complex, and `v` is a vertex of
`C`. `C` is convex, so it has a `w`-improving edge at `v`, and edges of cells
are edges of the complex. Hence an improving edge always exists unless `v` is
the `w`-maximum; `w` strictly increases along `parent`; the search terminates
with the `w`-maximal vertex as unique sink. Root there, reached by ascent from
any starting vertex.

Working inside the convex region `F` is exactly what makes this valid. The
naive alternative — improving pivots on the vertex graph of `U` itself — fails,
because reflex vertices of a non-convex union create extra local maxima and the
local-search function becomes a forest rather than a tree. No pivot rule
repairs that; the fix is to restrict to a convex region, which is what §2.2
does.

### 2.4 Adjacency

The edges of the complex at `v` are the extreme rays of the cells of the central
arrangement of `H_v`. **This is the same object §3 needs.** Compute it once
(§4). Given an edge direction `r`, the neighbour comes from a ratio test,
`t* = min { t > 0 : v + t·r ∈ h }` over `h ∈ G` not through `v`; because `F` is
bounded with its bounding hyperplanes in `G`, `t*` always exists, so every edge
ends at another vertex. Share this code with the pivoting layer.

### 2.5 Deduplication without storage

A `∂U` vertex generally lies on facets of several polytopes and is reached by
several facet searches. Emit `v` only from the lex-smallest
`(polytope index, facet index)` among the facets containing it — a local test
from `H_v`, no global store. A vertex skipped in `F` is still found in its owner
facet `F'`, since its `d` independent tight normals in `R^d` restrict to `d−1`
independent ones in `aff(F')`.

### 2.6 Degeneracy

Many concurrent hyperplanes is the normal case here, not an exception. Break
`w`-ties by symbolic perturbation, never by random jitter. Test with inputs
where dozens of hyperplanes meet at a point.

## 3. The sweep-plane half

### 3.1 Local contribution of a simplicial cone

For `K = v + cone(u_1, …, u_d)` with independent generators and
`⟨a, u_j⟩ ≠ 0` for all `j`:

    contrib(K) = |det[u_1 … u_d]| · (λ − ⟨a,v⟩)_+^d
                 ────────────────────────────────────
                      d! · ∏_j ⟨a, u_j⟩

The product is **signed**; the term's sign is `(−1)^{#{ j : ⟨a,u_j⟩ < 0 }}`.
This is the "modulo cones containing a line" valuation identity: flipping
`u_j → −u_j` negates the class. Cones containing a line contribute zero, which
is what kills interior vertices in §2.1.

### 3.2 Global identity

    s(a, λ) = Σ_v Σ_{K ∈ Δ(T_v(U))} contrib(K)

over arrangement vertices `v`, where `T_v(U)` is the tangent cone (in general
non-convex) and `Δ` is any decomposition into simplicial cones with disjoint
interiors. Every term attaches to the single knot `⟨a,v⟩`, so the whole answer
accumulates into `{ knot : coefficient }` with
`s(a,λ) = Σ_κ α_κ (λ − κ)_+^d`. That map is Bieri–Nef's "local information
collected at every vertex", and it is the only thing the traversal stores.

### 3.3 Self-check to implement first

Unit square `[0,1]^2`, `a = (1,2)`. Vertex contributions: `λ²/4`,
`−(λ−1)_+²/4`, `−(λ−2)_+²/4`, `(λ−3)_+²/4`. At `λ=5`: `(25−16−9+4)/4 = 1`.
At `λ=1`: `1/4`. Get this passing before anything else — it pins down the sign
and orientation conventions.

### 3.4 Genericity of `a`

The condition is `⟨a, u_j⟩ ≠ 0` for every generator of every simplicial cone,
i.e. `a` not orthogonal to any edge direction. It is **not** required that
vertices have distinct `⟨a,v⟩`; coinciding knots simply merge in the
accumulator.

When `a` is non-generic the individual terms have poles that cancel in the sum.
Handle this properly rather than by jittering `a`, which would move the
breakpoints and answer a different question: use `a_ε = a + ε r` for random `r`,
carry each contribution as a truncated Laurent series in `ε` (orders `−d … 0`),
sum, take the constant term. Provide `--assume-generic` to skip it.

## 4. The fused per-vertex kernel

This is the heart of the implementation. Write it as one function.

    fn visit(v: Vertex) -> (impl Iterator<Item = Vertex>, KnotContributions)

Given `v`, compute `H_v` and the central arrangement of `H_v` **once**, then:

- **for reverse search**: its extreme rays are the edge directions; ratio-test
  each for the neighbour; keep those `u` with `parent(u) == v` as children.
- **for the sweep plane**: its full-dimensional cells, tested for membership in
  `U` by evaluating at `v + ε w` for `w` interior to the cell, triangulated into
  simplicial cones, each contributing per §3.1.

### 4.1 One factorisation yields everything

For a `d`-subset `B ⊆ H_v` with normal matrix `N` (rows `n_j`, constraints
`n_j·x ≤ β_j`) and a sign pattern:

- `v = N^{-1} β_B`
- generators `u_j = −N^{-1} e_j`, so `|det[u_1 … u_d]| = 1/|det N|`
- with `y` solving `N^T y = a`, `⟨a, u_j⟩ = −y_j`

giving `contrib = (λ − ⟨a,v⟩)_+^d / ( d! · |det N| · ∏_j (−y_j) )`.

**A single LU of `N` yields the vertex, the determinant, the edge directions and
the dual vector `y`** — every quantity both halves need. Never recompute cones
in a second pass. This is precisely the data an LP dictionary already carries,
and precisely what `lrs --printcobasis` emits as the cobasis plus `det=`; if you
ever swap in lrs or cddlib for the local work, that is the interface.

### 4.2 Caching

The same combinatorial type of `H_v` recurs constantly in structured inputs.
Cache the local cell/ray structure keyed by that type. Reuse scratch buffers for
`N`, the LU and the sign-vector enumeration so the kernel allocates nothing per
vertex.

## 5. Rust implementation notes

- **Arithmetic.** Hybrid escalation like lrs: `i64`, fall back to `i128`, then
  arbitrary precision (`rug` if GMP is acceptable, else `malachite` or
  `num-bigint`). Store coordinates as integer numerators over a common
  denominator, not per-entry rationals; normalise rarely. Knot coefficients
  cancel heavily — never floats there.
- **Parallelism** (`rayon`), two independent levels, both from §2:
  *coarse* — facet subproblems are fully independent; parallelise over
  `(polytope, facet)` first. *fine* — within one search the subtrees are
  independent; split mplrs-style: descend to an initial depth, push subtree
  roots onto a work queue, workers pull, each with a node budget after which
  unfinished subtrees return to the queue. Reverse search has no shared state,
  so the queue is the only synchronisation. Per-worker knot `HashMap`, merged at
  the end; never a shared concurrent map in the kernel.
- **Determinism.** Results must not depend on thread count. Exact arithmetic
  throughout makes summation order irrelevant — prefer that over imposing a
  canonical ordering at merge time.

## 6. CLI

    nefvol volume  --input <file> [--direction a1,a2,...]
    nefvol sweep   --input <file> --direction a1,a2,... [--at λ]... [--json]
    nefvol plot    --input <file> --direction a1,a2,... --samples N

- Input: lrs/cdd `.ine` H-representation; multiple `begin…end` blocks in one
  file mean a union, or pass several `--input` paths. Also a simple JSON schema.
- `--direction` defaults to a generic pseudo-random rational direction.
- `sweep` prints breakpoints and per-interval polynomial coefficients as exact
  rationals; `--at λ` for point evaluation.
- `--threads N`, `--assume-generic`, `--verbose`, `--stats-json`.
- Non-zero exit with a clear message on unbounded or lower-dimensional input.

## 7. Testing — the two halves separately, then together

**Reverse search alone** (no volume):

- Differential test against a brute-force `d`-subset enumerator on small
  instances: identical `∂U` vertex sets, each emitted exactly once.
- Debug assertions: `parent` strictly increases `w`; `parent(child) == v` for
  every child enumerated at `v`; the root is the unique `w`-maximum.
- Memory profile: peak RSS flat as vertex count grows by orders of magnitude.

**Sweep plane alone** (fed a known vertex/cone list):

- §3.3 self-check, then golden exact values for boxes and simplices, `d = 2..5`.
- **Invariant:** `s(a,λ) = 0` below the smallest knot and constant above the
  largest — equivalently the degree `1..d` coefficients cancel in both outer
  intervals. Strong check on the signs, nearly free.
- **Direction independence:** total volume identical for many random `a`.

**Composed:**

- Golden values on L-shapes, nested, overlapping and disjoint unions.
- Monte Carlo cross-checks on rotated, non-axis-aligned unions.
- `vol(P ∪ Q) + vol(P ∩ Q) = vol(P) + vol(Q)` on random pairs.
- Total volume against `lrs --volume` per convex piece, and against `polymake`.
- Deliberately degenerate inputs: concurrent hyperplanes, coincident facets,
  `a` orthogonal to edges, disconnected unions.

## 8. Milestones

1. Exact arithmetic layer; §3.3 self-check on a single simple polytope.
2. Sweep-plane half complete for one convex polytope, arbitrary degeneracy,
   driven by a brute-force `d`-subset vertex enumerator.
3. Extend to unions: `∂U` restriction, local cell membership, cone
   triangulation. Still brute-force enumeration. **This is the oracle.**
4. **Reverse search half**, replacing the enumerator, with the fused kernel of
   §4. Differential-test against milestone 3 before proceeding. Verify memory
   is flat.
5. Non-generic `a` via Laurent series in `ε`.
6. Parallelisation at both levels; hybrid arithmetic escalation; benchmarks.
7. CLI and I/O.

Report at each milestone: vertices visited, reverse-search tree depth and node
count, subtree size histogram, simplicial cones generated, peak memory, wall
time, speedup vs. thread count.

## 9. Non-goals

Full Nef-polyhedron generality (complements, open/closed faces), unbounded
polyhedra, floating-point fast paths. Structure the tangent-cone representation
so Nef generality can be added later without touching the accumulator.

A cell-decomposition fallback (reverse search over the *cells* of `A(H)` inside
each convex `P_i`, then triangulating cells) is easier to get right but does
strictly more work, stores more, and abandons the vertex-local character that
makes this design worth building. It is a last resort, not a starting point.
