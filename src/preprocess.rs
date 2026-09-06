//! Per-polytope preprocessing with exact simplex walks: find a vertex, check
//! boundedness and full-dimensionality, strip redundant constraints and find
//! one vertex on every facet (the start of that facet's reverse search).

use crate::geom::{Constraint, Problem, RawConstraint};
use crate::linalg::{det_adj, Echelon};
use crate::num::{Point, Z};
use crate::walk::{walk, Region, WalkError};

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum PrepError {
    NoVertex,
    Unbounded,
    LowerDimensional,
}

impl std::fmt::Display for PrepError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            PrepError::NoVertex => write!(f, "polytope has no vertex: it is empty, unbounded or not full-dimensional"),
            PrepError::Unbounded => write!(f, "polytope is unbounded"),
            PrepError::LowerDimensional => write!(f, "polytope is not full-dimensional"),
        }
    }
}

#[derive(Debug, Clone)]
pub struct Prepared {
    /// irredundant facet constraints, in input order
    pub cons: Vec<RawConstraint>,
    /// one vertex on each facet, parallel to `cons`
    pub facet_starts: Vec<Point>,
    pub removed_redundant: usize,
    /// bounding box (lo, hi) from the boundedness walks
    pub bbox: (Vec<crate::num::Q>, Vec<crate::num::Q>),
}

/// First feasible basis of `region` in the local arrangement, by pruned
/// subset enumeration.
pub fn find_vertex(prob: &Problem, region: &Region) -> Option<Point> {
    let d = prob.d;
    let m = prob.hyps.len();
    let mut ech = Echelon::new(d);
    let mut choice: Vec<usize> = Vec::new();
    let mut next = 0usize;
    loop {
        if choice.len() == d {
            let rows: Vec<&[Z]> = choice.iter().map(|&h| prob.hyps[h].n.as_slice()).collect();
            let (det, adj) = det_adj(&rows);
            if !det.is_zero() {
                let num: Vec<Z> = (0..d)
                    .map(|i| {
                        let mut s = Z::ZERO;
                        for (k, &h) in choice.iter().enumerate() {
                            s += &(&adj[i * d + k] * &prob.hyps[h].b);
                        }
                        s
                    })
                    .collect();
                let p = Point::new(num, det);
                if region.contains(prob, &p) {
                    return Some(p);
                }
            }
            let c = choice.pop().unwrap();
            ech.pop();
            next = c + 1;
            continue;
        }
        let slots = d - choice.len();
        if next + slots > m {
            let c = choice.pop()?;
            ech.pop();
            next = c + 1;
            continue;
        }
        if ech.try_push(&prob.hyps[next].n) {
            choice.push(next);
        }
        next += 1;
    }
}

fn objective(prob: &Problem, c: Constraint, sign: i8) -> Vec<Z> {
    let s = Z::from((c.o * sign) as i64);
    prob.hyps[c.h].n.iter().map(|x| x * &s).collect()
}

/// Full preprocessing of one polytope given as raw constraints.
pub fn prepare(d: usize, raw: &[RawConstraint]) -> Result<Prepared, PrepError> {
    let local = Problem::from_raw(d, &[raw.to_vec()]).map_err(|_| PrepError::NoVertex)?;
    let mut cons: Vec<Constraint> = local.polys[0].cons.clone();
    if cons.len() < d + 1 {
        // fewer than d+1 constraints cannot bound a full-dimensional polytope
        if cons.is_empty() || d == 0 {
            return Err(PrepError::Unbounded);
        }
    }
    let region = Region { cons: cons.clone(), eq: None };
    let v0 = find_vertex(&local, &region).ok_or(PrepError::NoVertex)?;

    // boundedness: max ±x_k for all k
    let mut lo = vec![crate::num::Q::zero(); d];
    let mut hi = vec![crate::num::Q::zero(); d];
    for k in 0..d {
        for s in [1i64, -1] {
            let mut w = vec![Z::ZERO; d];
            w[k] = Z::from(s);
            let u = walk(&local, &region, &v0, Some(&w)).map_err(|WalkError::Unbounded| PrepError::Unbounded)?;
            if s > 0 {
                hi[k] = u.coord(k);
            } else {
                lo[k] = u.coord(k);
            }
        }
    }

    // implicit equalities (full-dimensionality) + witnesses off each hyperplane
    let mut witness: Vec<Point> = Vec::with_capacity(cons.len());
    for &c in &cons {
        let w = objective(&local, c, -1);
        let u = walk(&local, &region, &v0, Some(&w)).map_err(|_| PrepError::Unbounded)?;
        if local.hyps[c.h].contains(&u) {
            return Err(PrepError::LowerDimensional);
        }
        witness.push(u);
    }

    // redundancy: constraint j is redundant iff max of its functional over
    // P without j does not exceed b_j
    let mut keep = vec![true; cons.len()];
    let mut removed = 0usize;
    for j in 0..cons.len() {
        let sub: Vec<Constraint> = cons.iter().enumerate().filter(|(i, _)| *i != j && keep[*i]).map(|(_, c)| *c).collect();
        let r = Region { cons: sub, eq: None };
        let w = objective(&local, cons[j], 1);
        match walk(&local, &r, &witness[j], Some(&w)) {
            Err(WalkError::Unbounded) => {}
            Ok(u) => {
                let e = local.hyps[cons[j].h].eval(&u).signum() * cons[j].o;
                if e <= 0 {
                    keep[j] = false;
                    removed += 1;
                }
            }
        }
    }
    let kept: Vec<Constraint> = cons.iter().zip(&keep).filter(|(_, k)| **k).map(|(c, _)| *c).collect();
    cons = kept;
    let region = Region { cons: cons.clone(), eq: None };

    // one vertex on each facet
    let mut facet_starts = Vec::with_capacity(cons.len());
    for &c in &cons {
        let w = objective(&local, c, 1);
        let u = walk(&local, &region, &v0, Some(&w)).map_err(|_| PrepError::Unbounded)?;
        debug_assert!(local.hyps[c.h].contains(&u), "facet maximum not tight");
        facet_starts.push(u);
    }

    let out_cons = cons
        .iter()
        .map(|c| {
            let h = &local.hyps[c.h];
            let o = Z::from(c.o as i64);
            RawConstraint {
                a: h.n.iter().map(|x| crate::num::Q::int(x * &o)).collect(),
                b: crate::num::Q::int(&h.b * &o),
            }
        })
        .collect();
    Ok(Prepared { cons: out_cons, facet_starts, removed_redundant: removed, bbox: (lo, hi) })
}
