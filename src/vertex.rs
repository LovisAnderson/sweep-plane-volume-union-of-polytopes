//! Processing of one arrangement vertex: `∂U` test, ownership (SPEC §2.5)
//! and the sweep-plane contribution.

use crate::geom::{Problem, Side, Tight};
use crate::local::{cones_at, Containing, KernelError, KernelStats, SweepParams};
use crate::num::{Point, Z};
use crate::sweep::Accumulator;

/// Outcome of looking at a vertex.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Visit {
    /// `v ∈ int U` or `v ∉ U`: no contribution.
    NotBoundary,
    /// on `∂U` but another facet owns it.
    OwnedElsewhere,
    /// contribution accumulated.
    Emitted,
}

/// Owner facet of a boundary vertex: lex-smallest `(polytope, constraint
/// position)` over all polytopes containing `v` on their boundary.
pub fn owner(containing: &[(usize, Vec<Tight>)]) -> (usize, usize) {
    containing
        .iter()
        .map(|(p, t)| (*p, t.iter().map(|x| x.pos).min().unwrap()))
        .min()
        .unwrap()
}

/// Classify `v` against all polytopes.  `None` if `v` is interior to some
/// polytope or outside all of them; otherwise the containing list.
pub fn containing_polytopes(prob: &Problem, v: &Point) -> Option<Vec<(usize, Vec<Tight>)>> {
    containing_from_evals(prob, &prob.evals(v))
}

/// Same from precomputed hyperplane evaluations.
pub fn containing_from_evals(prob: &Problem, evals: &[Z]) -> Option<Vec<(usize, Vec<Tight>)>> {
    let mut out = Vec::new();
    for i in 0..prob.polys.len() {
        let (side, tight) = prob.classify_evals(i, evals);
        match side {
            Side::Interior => return None,
            Side::Boundary => out.push((prob.polys[i].id, tight)),
            Side::Outside => {}
        }
    }
    if out.is_empty() {
        None
    } else {
        Some(out)
    }
}

/// Full per-vertex step.  If `expect_owner` is given, the vertex is only
/// emitted when that facet owns it.
pub fn process_vertex(
    prob: &Problem,
    v: &Point,
    evals: &[Z],
    expect_owner: Option<(usize, usize)>,
    params: &SweepParams,
    acc: &mut Accumulator,
    stats: &mut KernelStats,
) -> Result<Visit, KernelError> {
    let containing = match containing_from_evals(prob, evals) {
        Some(c) => c,
        None => return Ok(Visit::NotBoundary),
    };
    if let Some(o) = expect_owner {
        if owner(&containing) != o {
            return Ok(Visit::OwnedElsewhere);
        }
    }
    let cont: Vec<Containing> = containing.iter().map(|(p, t)| Containing { poly: *p, tight: t.as_slice() }).collect();
    // A basis containing a hyperplane that is tight for no containing polytope
    // has chi = 0 (no cube fixes that coordinate), so only the hyperplanes
    // tight for some containing polytope are enumerated (ascending, as
    // required by the lexicographic cell test).
    let mut tv: Vec<usize> = containing.iter().flat_map(|(_, t)| t.iter().map(|x| x.c.h)).collect();
    tv.sort_unstable();
    tv.dedup();
    let contrib = cones_at(prob, v, &tv, &cont, params, stats)?;
    let knot = v.dot_q(params.a);
    acc.add(knot, &contrib);
    Ok(Visit::Emitted)
}
