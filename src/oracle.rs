//! Brute-force `d`-subset vertex enumerator (SPEC §8, milestone 3).  This is
//! the oracle the reverse search is differential-tested against; it stores
//! every vertex, which the real traversal never does.

use crate::geom::Problem;
use crate::linalg::{det_adj, Echelon};
use crate::local::{KernelError, KernelStats, SweepParams};
use crate::num::{Point, Z};
use crate::sweep::Accumulator;
use crate::vertex::{containing_polytopes, process_vertex};
use std::collections::HashSet;

/// All distinct vertices of the arrangement `A(H)`.
pub fn all_vertices(prob: &Problem) -> Vec<Point> {
    let d = prob.d;
    let m = prob.hyps.len();
    let mut set: HashSet<Point> = HashSet::new();
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
                set.insert(Point::new(num, det));
            }
            let c = choice.pop().unwrap();
            ech.pop();
            next = c + 1;
            continue;
        }
        let slots = d - choice.len();
        if next + slots > m {
            match choice.pop() {
                Some(c) => {
                    ech.pop();
                    next = c + 1;
                    continue;
                }
                None => break,
            }
        }
        if ech.try_push(&prob.hyps[next].n) {
            choice.push(next);
        }
        next += 1;
    }
    let mut v: Vec<Point> = set.into_iter().collect();
    v.sort_by_key(|a| a.coords());
    v
}

/// Vertices of `A(H)` lying on `∂U`.
pub fn boundary_vertices(prob: &Problem) -> Vec<Point> {
    all_vertices(prob).into_iter().filter(|v| containing_polytopes(prob, v).is_some()).collect()
}

/// Sweep function by brute force.
pub fn oracle_sweep(prob: &Problem, params: &SweepParams) -> Result<(Accumulator, KernelStats, usize), KernelError> {
    let mut acc = Accumulator::new(prob.d);
    let mut stats = KernelStats::default();
    let mut emitted = 0usize;
    for v in all_vertices(prob) {
        let evals = prob.evals(&v);
        if process_vertex(prob, &v, &evals, None, params, &mut acc, &mut stats)? == crate::vertex::Visit::Emitted {
            emitted += 1;
        }
    }
    Ok((acc, stats, emitted))
}
