//! Hyperplanes, polytopes and the shared arrangement (SPEC §1).
//!
//! Every facet hyperplane of every polytope is canonicalised to primitive
//! integer form `n·x = b` with the first non-zero entry of `n` positive and
//! deduplicated into one list `H`.  A polytope is a list of oriented references
//! `(h, o)` meaning `o·(n_h·x − b_h) ≤ 0`.

use crate::num::{dot, make_primitive, lex_sign, Point, Q, Z};
use std::collections::HashMap;

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct Hyperplane {
    pub n: Vec<Z>,
    pub b: Z,
}

impl Hyperplane {
    /// `n·x − b` scaled by `den(p)`: the sign is the side of `p`.
    #[inline]
    pub fn eval(&self, p: &Point) -> Z {
        &dot(&self.n, &p.num) - &(&self.b * &p.den)
    }
    #[inline]
    pub fn contains(&self, p: &Point) -> bool {
        self.eval(p).is_zero()
    }
    pub fn dim(&self) -> usize {
        self.n.len()
    }
}

/// `o · (n_h·x − b_h) ≤ 0`
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct Constraint {
    pub h: usize,
    pub o: i8,
}

#[derive(Clone, Debug, Default)]
pub struct Polytope {
    pub cons: Vec<Constraint>,
    /// global polytope index (differs from the position only in a local view)
    pub id: usize,
    /// global constraint position of every entry of `cons`
    pub pos: Vec<usize>,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Side {
    Outside,
    Boundary,
    Interior,
}

/// One tight constraint of a polytope at a point.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Tight {
    /// position of the constraint within `Polytope::cons`
    pub pos: usize,
    pub c: Constraint,
}

/// Flat `i64` mirror of the hyperplanes for the hot scans (SPEC §5, hybrid
/// arithmetic).  Only built when every entry is below `2^31` in magnitude,
/// so that with point coordinates below `2^62` every intermediate fits an
/// `i128` without checks.
#[derive(Clone, Debug)]
pub struct FastHyps {
    pub n: Vec<i64>,
    pub b: Vec<i64>,
}

pub const FAST_COORD_BOUND: i128 = 1 << 62;
const FAST_ENTRY_BOUND: i128 = 1 << 31;

#[derive(Clone, Debug)]
pub struct Problem {
    pub d: usize,
    pub hyps: Vec<Hyperplane>,
    pub polys: Vec<Polytope>,
    pub fast: Option<FastHyps>,
}

/// Point as small integers if every numerator and the denominator is below
/// [`FAST_COORD_BOUND`] in magnitude.
pub fn small_point(p: &Point) -> Option<(Vec<i64>, i64)> {
    let den = match &p.den {
        Z::S(x) if x.abs() < FAST_COORD_BOUND => *x as i64,
        _ => return None,
    };
    let mut num = Vec::with_capacity(p.num.len());
    for x in &p.num {
        match x {
            Z::S(v) if v.abs() < FAST_COORD_BOUND => num.push(*v as i64),
            _ => return None,
        }
    }
    Some((num, den))
}

pub fn small_vec(v: &[Z]) -> Option<Vec<i64>> {
    let mut out = Vec::with_capacity(v.len());
    for x in v {
        match x {
            Z::S(t) if t.abs() < FAST_COORD_BOUND => out.push(*t as i64),
            _ => return None,
        }
    }
    Some(out)
}

/// Raw input constraint `a·x ≤ b` in rationals.
#[derive(Clone, Debug)]
pub struct RawConstraint {
    pub a: Vec<Q>,
    pub b: Q,
}

/// Canonicalise `a·x ≤ b` into `(n, b, o)` with `n` primitive, lex-positive.
/// Returns `None` if `a = 0`.
pub fn canonicalise(a: &[Q], b: &Q) -> Option<(Hyperplane, i8)> {
    // clear denominators
    let mut den = b.den.clone();
    for c in a {
        den = &den * &c.den.div_exact(&den.gcd(&c.den));
    }
    let mut n: Vec<Z> = a.iter().map(|c| &c.num * &den.div_exact(&c.den)).collect();
    let mut bb = &b.num * &den.div_exact(&b.den);
    let s = lex_sign(&n);
    if s == 0 {
        return None;
    }
    // gcd over n and b
    let mut all = n.clone();
    all.push(bb.clone());
    make_primitive(&mut all);
    bb = all.pop().unwrap();
    n = all;
    if s < 0 {
        for x in n.iter_mut() {
            *x = -&*x;
        }
        bb = -bb;
    }
    Some((Hyperplane { n, b: bb }, s))
}

impl Problem {
    /// Build the shared arrangement from raw polytopes.  Constraints with a
    /// zero normal are rejected (`b ≥ 0` would be trivially true, `b < 0`
    /// infeasible) with an error.
    pub fn from_raw(d: usize, raw: &[Vec<RawConstraint>]) -> Result<Problem, String> {
        let mut hyps: Vec<Hyperplane> = Vec::new();
        let mut index: HashMap<Hyperplane, usize> = HashMap::new();
        let mut polys = Vec::new();
        for (pi, rc) in raw.iter().enumerate() {
            let mut cons = Vec::new();
            for (ci, c) in rc.iter().enumerate() {
                if c.a.len() != d {
                    return Err(format!("polytope {pi} constraint {ci}: expected {d} coefficients, got {}", c.a.len()));
                }
                match canonicalise(&c.a, &c.b) {
                    Some((h, o)) => {
                        let idx = match index.get(&h) {
                            Some(&i) => i,
                            None => {
                                hyps.push(h.clone());
                                index.insert(h, hyps.len() - 1);
                                hyps.len() - 1
                            }
                        };
                        let c = Constraint { h: idx, o };
                        if !cons.contains(&c) {
                            cons.push(c);
                        }
                    }
                    None => {
                        if c.b.signum() < 0 {
                            return Err(format!("polytope {pi} constraint {ci}: 0 ≤ {} is infeasible", c.b));
                        }
                        // trivially satisfied: drop
                    }
                }
            }
            let pos = (0..cons.len()).collect();
            polys.push(Polytope { cons, id: pi, pos });
        }
        let mut p = Problem { d, hyps, polys, fast: None };
        p.build_fast();
        Ok(p)
    }

    pub fn build_fast(&mut self) {
        let mut n = Vec::with_capacity(self.hyps.len() * self.d);
        let mut b = Vec::with_capacity(self.hyps.len());
        for h in &self.hyps {
            for x in &h.n {
                match x {
                    Z::S(v) if v.abs() < FAST_ENTRY_BOUND => n.push(*v as i64),
                    _ => return,
                }
            }
            match &h.b {
                Z::S(v) if v.abs() < FAST_ENTRY_BOUND => b.push(*v as i64),
                _ => return,
            }
        }
        self.fast = Some(FastHyps { n, b });
    }

    /// Raw constraints of polytope `i` (for building sub-problems).
    pub fn raw_of(&self, i: usize) -> Vec<RawConstraint> {
        self.polys[i]
            .cons
            .iter()
            .map(|c| {
                let h = &self.hyps[c.h];
                let o = Z::from(c.o as i64);
                RawConstraint {
                    a: h.n.iter().map(|x| Q::int(x * &o)).collect(),
                    b: Q::int(&h.b * &o),
                }
            })
            .collect()
    }

    /// Indices of hyperplanes through `p`, ascending.
    pub fn hyperplanes_through(&self, p: &Point) -> Vec<usize> {
        self.hyps.iter().enumerate().filter(|(_, h)| h.contains(p)).map(|(i, _)| i).collect()
    }

    /// `eval` of every hyperplane at `p`: the one scan per vertex that feeds
    /// `H_v`, the polytope classification and the ratio tests.
    pub fn evals(&self, p: &Point) -> Vec<Z> {
        if let Some(f) = &self.fast {
            if let Some((num, den)) = small_point(p) {
                let d = self.d;
                let den = den as i128;
                let mut out = Vec::with_capacity(self.hyps.len());
                for h in 0..self.hyps.len() {
                    let row = &f.n[h * d..(h + 1) * d];
                    let mut acc: i128 = 0;
                    for (a, x) in row.iter().zip(&num) {
                        acc += *a as i128 * *x as i128;
                    }
                    out.push(Z::S(acc - f.b[h] as i128 * den));
                }
                return out;
            }
        }
        self.hyps.iter().map(|h| h.eval(p)).collect()
    }

    pub fn hv_from_evals(evals: &[Z]) -> Vec<usize> {
        evals.iter().enumerate().filter(|(_, e)| e.is_zero()).map(|(i, _)| i).collect()
    }

    /// Like [`Problem::classify`] but from precomputed evaluations.
    pub fn classify_evals(&self, i: usize, evals: &[Z]) -> (Side, Vec<Tight>) {
        let mut tight = Vec::new();
        for (pos, c) in self.polys[i].cons.iter().enumerate() {
            let s = evals[c.h].signum() * c.o;
            if s > 0 {
                return (Side::Outside, Vec::new());
            }
            if s == 0 {
                tight.push(Tight { pos: self.polys[i].pos[pos], c: *c });
            }
        }
        if tight.is_empty() {
            (Side::Interior, tight)
        } else {
            (Side::Boundary, tight)
        }
    }

    /// Restricted view for the facet searches of polytope `i` (SPEC §2.2):
    /// only polytopes whose bounding box meets `bbox_i`, and of those only the
    /// constraints whose hyperplane meets `bbox_i`.  Polytope `i` comes first
    /// (local index 0).  Hyperplane order is preserved (the lexicographic
    /// perturbation only depends on relative order).
    pub fn local_view(&self, i: usize, bboxes: &[(Vec<Q>, Vec<Q>)]) -> Problem {
        let d = self.d;
        let (lo_i, hi_i) = &bboxes[i];
        let mut used = vec![false; self.hyps.len()];
        let mut polys_raw: Vec<(usize, Vec<Constraint>, Vec<usize>)> = Vec::new();
        let order = std::iter::once(i).chain((0..self.polys.len()).filter(|&j| j != i));
        'poly: for j in order {
            let (lo_j, hi_j) = &bboxes[j];
            if (0..d).any(|k| lo_j[k] > hi_i[k] || lo_i[k] > hi_j[k]) {
                continue;
            }
            let mut cons = Vec::new();
            let mut pos = Vec::new();
            for (p, c) in self.polys[j].cons.iter().enumerate() {
                let h = &self.hyps[c.h];
                let mut mn = Q::int(-&h.b);
                let mut mx = mn.clone();
                for k in 0..d {
                    let a = Q::int(h.n[k].clone());
                    let x1 = a.mul(&lo_i[k]);
                    let x2 = a.mul(&hi_i[k]);
                    let (l, u) = if x1 <= x2 { (x1, x2) } else { (x2, x1) };
                    mn = mn.add(&l);
                    mx = mx.add(&u);
                }
                if mn.signum() <= 0 && mx.signum() >= 0 {
                    cons.push(*c);
                    pos.push(p);
                } else if mn.signum() * c.o > 0 {
                    continue 'poly; // violated on all of bbox_i: P_j misses it
                }
                // else strictly satisfied on bbox_i: never tight, never violated
            }
            for c in &cons {
                used[c.h] = true;
            }
            polys_raw.push((j, cons, pos));
        }
        let mut map = vec![usize::MAX; self.hyps.len()];
        let mut hyps = Vec::new();
        for (h, u) in used.iter().enumerate() {
            if *u {
                map[h] = hyps.len();
                hyps.push(self.hyps[h].clone());
            }
        }
        let polys = polys_raw
            .into_iter()
            .map(|(id, cons, pos)| Polytope { cons: cons.iter().map(|c| Constraint { h: map[c.h], o: c.o }).collect(), id, pos })
            .collect();
        let mut p = Problem { d, hyps, polys, fast: None };
        p.build_fast();
        p
    }

    /// Classify `p` against polytope `i`; the tight list is only meaningful
    /// when the side is `Boundary` (it is also filled for `Interior`, empty).
    pub fn classify(&self, i: usize, p: &Point) -> (Side, Vec<Tight>) {
        let mut tight = Vec::new();
        for (pos, c) in self.polys[i].cons.iter().enumerate() {
            let e = self.hyps[c.h].eval(p);
            let s = e.signum() * c.o;
            if s > 0 {
                return (Side::Outside, Vec::new());
            }
            if s == 0 {
                tight.push(Tight { pos: self.polys[i].pos[pos], c: *c });
            }
        }
        if tight.is_empty() {
            (Side::Interior, tight)
        } else {
            (Side::Boundary, tight)
        }
    }

    /// Value `o·(n·r)`: `≤ 0` means direction `r` respects constraint `c`.
    #[inline]
    pub fn dir_side(&self, c: Constraint, r: &[Z]) -> i8 {
        dot(&self.hyps[c.h].n, r).signum() * c.o
    }

    pub fn total_constraints(&self) -> usize {
        self.polys.iter().map(|p| p.cons.len()).sum()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn q(x: i64) -> Q {
        Q::int(x)
    }

    #[test]
    fn canonical_dedup() {
        // x ≤ 1 and -2x ≤ -2 → same hyperplane, opposite orientations
        let raw = vec![vec![
            RawConstraint { a: vec![q(1), q(0)], b: q(1) },
            RawConstraint { a: vec![q(-2), q(0)], b: q(-2) },
            RawConstraint { a: vec![Q::parse("1/2").unwrap(), q(0)], b: Q::parse("1/2").unwrap() },
        ]];
        let p = Problem::from_raw(2, &raw).unwrap();
        assert_eq!(p.hyps.len(), 1);
        assert_eq!(p.hyps[0].n, vec![Z::from(1), Z::from(0)]);
        assert_eq!(p.hyps[0].b, Z::from(1));
        assert_eq!(p.polys[0].cons, vec![Constraint { h: 0, o: 1 }, Constraint { h: 0, o: -1 }]);
    }
}
