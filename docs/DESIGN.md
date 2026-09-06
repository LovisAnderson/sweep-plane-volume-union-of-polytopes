# Design notes: how the two halves fuse

This document records the mathematics the implementation relies on, in the
form actually used by the code.  SPEC.md is the brief; this is what was built
and why it is correct.

## 1. Objects

* `H`: the deduplicated, primitive-integer facet hyperplanes `n_h·x = b_h`
  of all input polytopes, with the first non-zero entry of `n_h` positive.
* Polytope `P_i`: a list of oriented references `(h, o_ih)` meaning
  `o_ih·(n_h·x − b_h) ≤ 0`, `o ∈ {±1}`.
* `H_v ⊆ H`: hyperplanes through an arrangement vertex `v`.
* A *basis* `B ⊆ H_v`: `d` hyperplanes with independent normals, matrix
  `N_B` (rows `n_h`), `det = det N_B`, `adj = det·N_B^{-1}` (integer).

## 2. The sweep-plane half: one signed term per basis

### 2.1 The Lawrence valuation

For a simplicial cone `K = v + cone(u_1,…,u_d)` and a direction `a` with
`⟨a,u_j⟩ ≠ 0`:

    Ψ_λ(K) = |det[u_1…u_d]| · (λ − ⟨a,v⟩)_+^d / (d! · Π_j ⟨a,u_j⟩)

`Ψ_λ` extends uniquely to a valuation on polyhedra modulo those containing a
line (it is the degree-`d` part of Brion's exponential valuation).  Two facts
drive everything:

* **Brion for unions.** `[U] ≡ Σ_v [T_v(U)]` mod lines, `v` over the
  vertices of the arrangement, `T_v(U)` the tangent cone.  For a union this
  follows from inclusion–exclusion because `T_v` commutes with `∩` and `∪`.
* **Sweep = truncation.** The vertices that `{⟨a,x⟩ ≤ λ}` cuts into `U`
  have tangent cones containing the direction of the cut edge, and their
  Lawrence term carries the factor `(⟨a,v'⟩ − λ)^d = 0`.  Hence
  `vol(U ∩ {⟨a,x⟩ ≤ λ}) = Σ_{v: ⟨a,v⟩<λ} Ψ_λ(T_v(U))`, i.e. the `(·)_+`.

### 2.2 Decomposing `T_v(U)` without cells: symbolic perturbation

Instead of enumerating the cells of the central arrangement of `H_v` and
triangulating them, shift every hyperplane `h` to `n_h·x = b_h + ε^{h+1}`
(one shift per *hyperplane*, shared by every polytope using it).  For small
`ε > 0`:

* the perturbed arrangement is simple: every vertex `v_B` is the unique
  point on exactly the `d` hyperplanes of one basis `B`;
* `vol(U_ε ∩ {⟨a,x⟩ ≤ λ})` is a polynomial in `ε` that converges to the
  unperturbed value;
* at a simple vertex `v_B`, `T_{v_B}(U_ε)` is a union of some of the `2^d`
  simplicial cells `C_σ = {w : σ_k n_k·w ≤ 0, k ∈ B}`, `σ ∈ {±1}^d`, and
  `Ψ_λ(v_B + C_σ) = sgn(σ) · f(B) · (λ − ⟨a,v_B⟩)_+^d` with

      f(B) = det^d / (d! · |det| · Π_k (−Y_k)),   Y = adjᵀ a.

Cell `C_σ` lies in `U_ε` iff there is a polytope `P_i ∋ v` such that
(a) `v_B ∈ P_{i,ε}` and (b) `σ_k = o_ik` for every `k ∈ B` used by `P_i`.
Condition (a) is a **lexicographic sign test**: for each tight constraint `h`
of `P_i` not in `B`,

    n_h·v_B − b_h − ε_h = Σ_{k∈B} (R_hk / det) ε_k − ε_h,   R_hk = n_h · adj[:,k],

and the sign is that of the coefficient with the smallest hyperplane index
(`−1` if `h` itself comes first); feasible iff `o_ih · sign < 0`.  Everything
is integer arithmetic on `adj`.

Summing the sign-weighted cells gives the basis contribution
`χ_B · f(B) · (λ − ⟨a,v⟩)_+^d` with `χ_B = Σ_{σ passing} Π_k σ_k` — computed
by inclusion–exclusion over the sub-cubes of `{±1}^d` defined by the feasible
polytopes.  A sub-cube with a free coordinate has `χ = 0` (a cone containing
a line), so most bases cost only the cell test.

Taking `ε → 0` term by term (every term is a polynomial in `ε` with a fixed
cone) yields, for each unperturbed vertex `v`,

    contribution(v) = Σ_{B ⊆ H_v basis} χ_B · f(B) · (λ − ⟨a,v⟩)_+^d.

### 2.3 Why only `∂U` vertices

* `v ∉ U`: no `P_{i,ε}` contains any `v_B`, all `χ_B = 0`.
* `v ∈ int U`: the local perturbed union `W_ε` differs from `R^d` only inside
  a ball of radius `O(ε)` around `v`, so `Σ_B χ_B f(B) (λ−⟨a,v⟩)^d
  = Ψ_λ(W_ε) = −vol(hole) → 0`; being a polynomial in `ε` it is identically
  zero at `ε = 0`.

So the reverse search only needs vertices with `v ∈ U` and `v ∉ int P_i` for
all `i` — a local test.  (Note that interior vertices of the union may have
non-zero *individual* basis terms — e.g. three cones meeting at a point whose
perturbation opens a triangular hole — but their sum vanishes, so skipping the
vertex is exact.)

### 2.4 Non-generic sweep directions

If `Y_k = 0` for a basis with `χ_B ≠ 0`, `a` is orthogonal to an arrangement
edge.  With `a_ε = a + ε r` (fixed pseudo-random integer `r`) each such term
is a Laurent series in `ε` of order `≥ −p` (`p` = number of zeros), and
`(λ − ⟨a_ε,v⟩)^d` expands in powers of `ε` too.  The constant term of the
product is accumulated per knot into coefficients of `(λ−κ)^j`, `j ≤ d`
(lower degrees genuinely appear: the square swept by `(1,0)` gives `s = λ`).
Negative orders are accumulated as well and must cancel per knot; the
accumulator verifies this and reports an internal error otherwise.

## 3. The reverse-search half

Per facet `F = h_j ∩ P_i`, in the full `R^d` coordinates (no
`(d−1)`-parametrisation; the equality `h_j` is a region constraint):

* **Vertices**: points of `F` on `d` independent hyperplanes of `H`.
* **Edges** at `v`: the lines of the central arrangement of `H_v` that lie in
  `h_j` (rank-`(d−1)` subsets containing `h_j`), both directions, filtered by
  the tangent cone of `F` (`o·n_h·r ≤ 0` for the tight constraints of `P_i`).
* **Objective**: pure lexicographic order on the edge direction
  (`w = 0` symbolically perturbed by `e_1, …, e_d`): no edge is flat, ties are
  impossible, and the lex-maximal vertex of the convex bounded `F` is the
  unique sink.
* **parent(v)**: the neighbour along the lexicographically smallest
  lex-positive feasible direction.  It is a function of `v` alone, so the
  search is memoryless.
* **Children**: for each lex-*negative* feasible direction `r` at `v`, take
  the ratio-test neighbour `u` and check `parent(u)` returns along `−r`.
* **Ownership**: a `∂U` vertex is emitted only by the lex-smallest
  `(polytope, constraint position)` among the tight facets containing it.

* **Restricted view**: the search for the facets of `P_i` runs in the
  arrangement of `H_i ⊆ H`, the hyperplanes of polytopes whose bounding box
  meets `bbox(P_i)`, minus constraints whose hyperplane misses `bbox(P_i)`
  (constant sign there: violated ⇒ the polytope is dropped, satisfied ⇒ the
  constraint can never be tight).  This is exact: a vertex contributes only
  if it is a vertex of the arrangement of the hyperplanes tight for the
  polytopes containing it, and all of those lie in `H_i`; the ownership rule
  sees the same containing set with the global `(polytope, position)`.  It
  makes the per-vertex cost independent of far-away polytopes and removes
  the arrangement vertices that far hyperplanes would cut into `F`.

The per-vertex arithmetic is one scan of all hyperplanes of the view (`n_h·v − b_h` for
every `h`, in a flat `i64/i128` fast path), which feeds `H_v`, the tangent
cone, the polytope classification and every ratio test from `v`.

## 4. The fused kernel

For a visited vertex the same `H_v` drives both halves:

* `lines_at(H_v, forced = h_j)` — the rank-`(d−1)` subsets, via an
  incremental fraction-free echelon form with dependent prefixes pruned;
* `cones_at(H_v)` — all bases via the same pruned DFS; per basis one
  fraction-free Gauss–Jordan pass over `[N_B | I]` returns `(det, adj)` from
  which the edge directions (columns of `adj`), the dual vector `Y` and the
  lexicographic rows `R_h` are read.

Nothing but the knot map `{⟨a,v⟩ ↦ coefficients}` survives a vertex.

## 5. Preprocessing

Each polytope is checked and cleaned with exact simplex walks (the same
`walk` used for the root ascent, on the polytope's own hyperplanes):

1. a first vertex by pruned subset enumeration;
2. boundedness: maximise `±x_k` for every `k`;
3. full-dimensionality: minimise every constraint functional; a minimum equal
   to `b_h` is an implicit equality;
4. redundancy: maximise each constraint functional over the polytope without
   it; if the maximum does not exceed `b_h` the constraint is dropped;
5. a start vertex on every remaining facet by maximising its functional.

## 6. Parallelism and determinism

Coarse level: facet tasks are independent (`rayon` parallel iterator).  Fine
level: each task explores with a node budget; when it is exhausted every stack
frame `(vertex, next edge index)` becomes a new task returned to the pool.
Restarting from `(vertex, index)` is exact because the edge list of a vertex
is a deterministic function of the vertex.  Each task owns its accumulator;
maps are merged on return.  All arithmetic is exact, so results are
independent of the thread count and the budget.
