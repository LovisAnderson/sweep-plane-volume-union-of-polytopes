//! Vertex-to-vertex moves shared by the reverse search, the root ascent and
//! the preprocessing LPs (SPEC §2.3–2.4).
//!
//! The objective is a linear functional `w` symbolically perturbed to
//! `(w·r, r_1, …, r_d)` compared lexicographically, so no edge is ever flat
//! and ties never occur; with `w = 0` this is the pure lexicographic order.

use crate::geom::{Constraint, Problem};
use crate::local::lines_at;
use crate::num::{dot, negate, lex_sign, Point, Z};
use std::cmp::Ordering;

/// A convex region: intersection of oriented constraints, optionally within
/// one hyperplane (`eq`), all referring to `Problem::hyps`.
#[derive(Clone, Debug)]
pub struct Region {
    pub cons: Vec<Constraint>,
    pub eq: Option<usize>,
}

impl Region {
    /// Tight constraints at `p` (which must lie in the region).
    pub fn tight(&self, prob: &Problem, p: &Point) -> Vec<Constraint> {
        self.cons.iter().copied().filter(|c| prob.hyps[c.h].contains(p)).collect()
    }
    /// Tight constraints from precomputed hyperplane evaluations.
    pub fn tight_evals(&self, evals: &[Z]) -> Vec<Constraint> {
        self.cons.iter().copied().filter(|c| evals[c.h].is_zero()).collect()
    }
    pub fn contains(&self, prob: &Problem, p: &Point) -> bool {
        if let Some(e) = self.eq {
            if !prob.hyps[e].contains(p) {
                return false;
            }
        }
        self.cons.iter().all(|c| prob.hyps[c.h].eval(p).signum() * c.o <= 0)
    }
}

/// Lexicographic sign of the perturbed objective along `r`.
#[inline]
pub fn objective_sign(w: Option<&[Z]>, r: &[Z]) -> i8 {
    if let Some(w) = w {
        let s = dot(w, r).signum();
        if s != 0 {
            return s;
        }
    }
    lex_sign(r)
}

/// Edge directions at `p` that stay inside the region: both signs of every
/// line of `H_p` (inside `eq` if present), filtered by the tangent cone.
pub fn feasible_dirs(prob: &Problem, region: &Region, hp: &[usize], tight: &[Constraint]) -> Vec<Vec<Z>> {
    let lines = lines_at(prob, hp, region.eq);
    let mut out = Vec::with_capacity(lines.len());
    for l in lines {
        let neg = negate(&l);
        if tight.iter().all(|&c| prob.dir_side(c, &l) <= 0) {
            out.push(l);
        }
        if tight.iter().all(|&c| prob.dir_side(c, &neg) <= 0) {
            out.push(neg);
        }
    }
    out
}

/// The canonical improving direction: lexicographically smallest among the
/// feasible directions with positive perturbed objective.  This *is* the
/// reverse-search `parent` rule (SPEC §2.3).
pub fn best_improving<'a>(dirs: &'a [Vec<Z>], w: Option<&[Z]>) -> Option<&'a Vec<Z>> {
    let mut best: Option<&Vec<Z>> = None;
    for r in dirs {
        if objective_sign(w, r) > 0 {
            match best {
                None => best = Some(r),
                Some(b) => {
                    if r.as_slice().cmp(b.as_slice()) == Ordering::Less {
                        best = Some(r);
                    }
                }
            }
        }
    }
    best
}

/// Ratio test: first hyperplane hit from `p` along `r` (over all hyperplanes
/// of the arrangement not through `p`).  `evals` are the hyperplane
/// evaluations at `p`.  `None` if unbounded.
pub fn ratio_step(prob: &Problem, evals: &[Z], p: &Point, r: &[Z]) -> Option<Point> {
    // t = A / (den·S) with A = b·den − n·num = −eval, S = n·r ; need t > 0
    let best: Option<(Z, Z)>; // (A, S) both made positive
    let fast = match (&prob.fast, crate::geom::small_vec(r)) {
        (Some(f), Some(r64)) => Some((f, r64)),
        _ => None,
    };
    if let Some((f, r64)) = fast {
        // all-i128 scan; the candidate comparison falls back to Z on overflow
        let d = prob.d;
        let mut bi: Option<(i128, i128)> = None;
        for h in 0..prob.hyps.len() {
            let e = match &evals[h] {
                Z::S(0) => continue,
                Z::S(e) => *e,
                Z::B(_) => return ratio_step_slow(prob, evals, p, r),
            };
            let row = &f.n[h * d..(h + 1) * d];
            let mut s: i128 = 0;
            for (a, x) in row.iter().zip(&r64) {
                s += *a as i128 * *x as i128;
            }
            if s == 0 || (e > 0) == (s > 0) {
                continue;
            }
            let (a, s) = if s < 0 { (e, -s) } else { (-e, s) };
            match bi {
                None => bi = Some((a, s)),
                Some((ba, bs)) => match (a.checked_mul(bs), ba.checked_mul(s)) {
                    (Some(l), Some(rr)) => {
                        if l < rr {
                            bi = Some((a, s));
                        }
                    }
                    _ => {
                        if &Z::S(a) * &Z::S(bs) < &Z::S(ba) * &Z::S(s) {
                            bi = Some((a, s));
                        }
                    }
                },
            }
        }
        best = bi.map(|(a, s)| (Z::S(a), Z::S(s)));
    } else {
        return ratio_step_slow(prob, evals, p, r);
    }
    let (a, s) = best?;
    // u = (num·S + A·r) / (den·S)
    let num: Vec<Z> = p.num.iter().zip(r).map(|(x, ri)| &(x * &s) + &(&a * ri)).collect();
    Some(Point::new(num, &p.den * &s))
}

fn ratio_step_slow(prob: &Problem, evals: &[Z], p: &Point, r: &[Z]) -> Option<Point> {
    let mut best: Option<(Z, Z)> = None;
    for (h, hyp) in prob.hyps.iter().enumerate() {
        let e = &evals[h];
        if e.is_zero() {
            continue;
        }
        let s = dot(&hyp.n, r);
        if s.is_zero() {
            continue;
        }
        if e.signum() == s.signum() {
            continue; // t ≤ 0
        }
        let a = -e;
        let (a, s) = if s.is_negative() { (-a, -s) } else { (a, s) };
        match &best {
            None => best = Some((a, s)),
            Some((ba, bs)) => {
                if &a * bs < ba * &s {
                    best = Some((a, s));
                }
            }
        }
    }
    let (a, s) = best?;
    // u = (num·S + A·r) / (den·S)
    let num: Vec<Z> = p.num.iter().zip(r).map(|(x, ri)| &(x * &s) + &(&a * ri)).collect();
    Some(Point::new(num, &p.den * &s))
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum WalkError {
    Unbounded,
}

/// Simplex-style ascent over arrangement vertices inside `region`, maximising
/// the perturbed objective.  Terminates at the unique lex-maximal vertex.
pub fn walk(prob: &Problem, region: &Region, start: &Point, w: Option<&[Z]>) -> Result<Point, WalkError> {
    let mut v = start.clone();
    debug_assert!(region.contains(prob, &v), "walk start outside region");
    let mut steps = 0usize;
    loop {
        let evals = prob.evals(&v);
        let hv = Problem::hv_from_evals(&evals);
        let tight = region.tight_evals(&evals);
        let dirs = feasible_dirs(prob, region, &hv, &tight);
        match best_improving(&dirs, w) {
            None => return Ok(v),
            Some(r) => match ratio_step(prob, &evals, &v, r) {
                None => return Err(WalkError::Unbounded),
                Some(u) => {
                    debug_assert!(region.contains(prob, &u));
                    v = u;
                }
            },
        }
        steps += 1;
        debug_assert!(steps < 10_000_000, "walk did not terminate");
    }
}
