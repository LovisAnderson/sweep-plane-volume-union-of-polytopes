//! Small constructors for test inputs and examples.

use crate::geom::RawConstraint;
use crate::num::Q;

pub fn q(x: i64) -> Q {
    Q::int(x)
}

/// Axis-aligned box `Π [lo_i, hi_i]`.
pub fn boxp(lo: &[Q], hi: &[Q]) -> Vec<RawConstraint> {
    let d = lo.len();
    let mut c = Vec::new();
    for i in 0..d {
        let mut a = vec![Q::zero(); d];
        a[i] = Q::one();
        c.push(RawConstraint { a: a.clone(), b: hi[i].clone() });
        a[i] = Q::int(-1);
        c.push(RawConstraint { a, b: lo[i].neg() });
    }
    c
}

pub fn unit_box(d: usize) -> Vec<RawConstraint> {
    boxp(&vec![Q::zero(); d], &vec![Q::one(); d])
}

pub fn int_box(lo: &[i64], hi: &[i64]) -> Vec<RawConstraint> {
    let l: Vec<Q> = lo.iter().map(|&x| q(x)).collect();
    let h: Vec<Q> = hi.iter().map(|&x| q(x)).collect();
    boxp(&l, &h)
}

/// Standard simplex `x ≥ 0, Σ x ≤ s`.
pub fn simplex(d: usize, s: Q) -> Vec<RawConstraint> {
    let mut c = Vec::new();
    for i in 0..d {
        let mut a = vec![Q::zero(); d];
        a[i] = Q::int(-1);
        c.push(RawConstraint { a, b: Q::zero() });
    }
    c.push(RawConstraint { a: vec![Q::one(); d], b: s });
    c
}

/// Generic constraint `a·x ≤ b` from integers.
pub fn con(a: &[i64], b: i64) -> RawConstraint {
    RawConstraint { a: a.iter().map(|&x| q(x)).collect(), b: q(b) }
}

/// Apply a unimodular-ish integer linear map `x ↦ M x + t` to a polytope
/// (constraints transform as `a·M^{-1}`; here we pass the inverse directly).
/// `minv` is the inverse matrix (row-major, rational).
pub fn transform(cons: &[RawConstraint], minv: &[Vec<Q>], t: &[Q]) -> Vec<RawConstraint> {
    // y = M x + t  ⇒ x = M^{-1}(y − t);  a·x ≤ b ⇒ (a M^{-1})·y ≤ b + (a M^{-1})·t
    let d = t.len();
    cons.iter()
        .map(|c| {
            let mut a2 = vec![Q::zero(); d];
            for (i, row) in minv.iter().enumerate() {
                for (j, m) in row.iter().enumerate() {
                    a2[j] = a2[j].add(&c.a[i].mul(m));
                }
            }
            let mut b2 = c.b.clone();
            for j in 0..d {
                b2 = b2.add(&a2[j].mul(&t[j]));
            }
            RawConstraint { a: a2, b: b2 }
        })
        .collect()
}
