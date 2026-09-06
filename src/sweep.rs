//! The knot map — the only thing the traversal stores (SPEC §3.2).
//!
//! `s(a, λ) = Σ_κ Σ_j c_{κ,j} (λ − κ)_+^j`.  For generic `a` only `j = d`
//! occurs; for non-generic `a` lower degrees appear (SPEC §3.4) and negative
//! Laurent orders are kept until [`Accumulator::finish`] verifies they cancel.

use crate::local::VertexContribution;
use crate::num::{binomial, Q};
use std::collections::HashMap;

#[derive(Clone, Debug)]
struct Entry {
    top: Q,
    laurent: Option<Vec<Q>>,
}

#[derive(Clone, Debug)]
pub struct Accumulator {
    pub d: usize,
    map: HashMap<Q, Entry>,
}

impl Accumulator {
    pub fn new(d: usize) -> Self {
        Accumulator { d, map: HashMap::new() }
    }
    pub fn len(&self) -> usize {
        self.map.len()
    }
    pub fn is_empty(&self) -> bool {
        self.map.is_empty()
    }

    pub fn add(&mut self, knot: Q, c: &VertexContribution) {
        if c.is_zero() {
            return;
        }
        let e = self.map.entry(knot).or_insert_with(|| Entry { top: Q::zero(), laurent: None });
        if !c.top.is_zero() {
            e.top = e.top.add(&c.top);
        }
        if let Some(l) = &c.laurent {
            let el = e.laurent.get_or_insert_with(|| vec![Q::zero(); l.len()]);
            for (x, y) in el.iter_mut().zip(l) {
                if !y.is_zero() {
                    *x = x.add(y);
                }
            }
        }
    }

    pub fn merge(&mut self, other: Accumulator) {
        for (k, e) in other.map {
            let c = VertexContribution { d: self.d, top: e.top, laurent: e.laurent };
            self.add(k, &c);
        }
    }

    /// Collapse into a sweep function; errors if the Laurent poles did not
    /// cancel (which would indicate a bug, not bad input).
    pub fn finish(self) -> Result<SweepFunction, String> {
        let d = self.d;
        let w = d + 1;
        let mut knots: Vec<(Q, Vec<Q>)> = Vec::with_capacity(self.map.len());
        for (k, e) in self.map {
            let mut coef = vec![Q::zero(); w];
            coef[d] = e.top;
            if let Some(l) = e.laurent {
                for m in 1..=d {
                    for j in 0..=d {
                        if !l[m * w + j].is_zero() {
                            return Err(format!(
                                "internal error: pole ε^-{m} at knot {k} did not cancel (coefficient of (λ-κ)^{j} is {})",
                                l[m * w + j]
                            ));
                        }
                    }
                }
                for j in 0..=d {
                    coef[j] = coef[j].add(&l[j]);
                }
            }
            if coef.iter().any(|c| !c.is_zero()) {
                knots.push((k, coef));
            }
        }
        knots.sort_by(|a, b| a.0.cmp(&b.0));
        Ok(SweepFunction { d, knots })
    }
}

/// `s(λ) = Σ_κ Σ_j c_{κ,j} (λ − κ)_+^j`, knots ascending.
#[derive(Clone, Debug)]
pub struct SweepFunction {
    pub d: usize,
    pub knots: Vec<(Q, Vec<Q>)>,
}

/// One polynomial piece in powers of `λ`.
#[derive(Clone, Debug)]
pub struct Piece {
    pub lo: Option<Q>,
    pub hi: Option<Q>,
    /// coefficients of λ^0 … λ^d
    pub poly: Vec<Q>,
}

impl SweepFunction {
    pub fn eval(&self, lambda: &Q) -> Q {
        let mut s = Q::zero();
        for (k, coef) in &self.knots {
            if lambda <= k {
                break;
            }
            let t = lambda.sub(k);
            let mut pw = Q::one();
            for (j, c) in coef.iter().enumerate() {
                if j > 0 {
                    pw = pw.mul(&t);
                }
                if !c.is_zero() {
                    s = s.add(&c.mul(&pw));
                }
            }
        }
        s
    }

    /// Polynomial pieces in powers of `λ` over each interval between knots.
    pub fn pieces(&self) -> Vec<Piece> {
        let d = self.d;
        let mut out = Vec::new();
        let mut poly = vec![Q::zero(); d + 1];
        out.push(Piece { lo: None, hi: self.knots.first().map(|k| k.0.clone()), poly: poly.clone() });
        for (i, (k, coef)) in self.knots.iter().enumerate() {
            // add Σ_j c_j (λ − k)^j = Σ_j c_j Σ_i C(j,i) λ^i (−k)^{j−i}
            for (j, c) in coef.iter().enumerate() {
                if c.is_zero() {
                    continue;
                }
                let negk = k.neg();
                for (l, pl) in poly.iter_mut().enumerate().take(j + 1) {
                    let term = c.mul(&Q::int(binomial(j, l))).mul(&negk.pow((j - l) as u32));
                    *pl = pl.add(&term);
                }
            }
            let hi = self.knots.get(i + 1).map(|k| k.0.clone());
            out.push(Piece { lo: Some(k.clone()), hi, poly: poly.clone() });
        }
        out
    }

    /// Total volume `s(+∞)`; also verifies the tail is constant.
    pub fn total(&self) -> Result<Q, String> {
        let pieces = self.pieces();
        let last = pieces.last().unwrap();
        for (i, c) in last.poly.iter().enumerate().skip(1) {
            if !c.is_zero() {
                return Err(format!("internal error: sweep function is not constant above the last knot (λ^{i} coefficient {c})"));
            }
        }
        Ok(last.poly[0].clone())
    }

    /// Reparametrise: if this is `s(a', λ')` with `a' = scale · a`, then
    /// `s(a, λ) = s(a', scale·λ)`: knots divide by `scale`, the coefficient of
    /// `(λ−κ)^j` multiplies by `scale^j`.
    pub fn rescale(&mut self, scale: &Q) {
        if scale == &Q::one() {
            return;
        }
        for (k, coef) in self.knots.iter_mut() {
            *k = k.div(scale);
            let mut pw = Q::one();
            for (j, c) in coef.iter_mut().enumerate() {
                if j > 0 {
                    pw = pw.mul(scale);
                }
                if !c.is_zero() {
                    *c = c.mul(&pw);
                }
            }
        }
        self.knots.sort_by(|a, b| a.0.cmp(&b.0));
    }

    pub fn min_knot(&self) -> Option<&Q> {
        self.knots.first().map(|k| &k.0)
    }
    pub fn max_knot(&self) -> Option<&Q> {
        self.knots.last().map(|k| &k.0)
    }
}

pub fn eval_poly(poly: &[Q], x: &Q) -> Q {
    let mut s = Q::zero();
    for c in poly.iter().rev() {
        s = s.mul(x).add(c);
    }
    s
}

