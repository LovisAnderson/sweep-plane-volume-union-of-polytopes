use nefvol::num::Q;
use nefvol::preprocess::{prepare, PrepError};
use nefvol::shapes::*;

#[test]
fn rejects_bad_input() {
    // unbounded: x ≥ 0, y ≥ 0
    let r = prepare(2, &[con(&[-1, 0], 0), con(&[0, -1], 0)]);
    assert_eq!(r.err(), Some(PrepError::Unbounded));
    // half-space plus box in x only
    let r = prepare(2, &[con(&[-1, 0], 0), con(&[1, 0], 1)]);
    assert_eq!(r.err(), Some(PrepError::NoVertex));
    // lower-dimensional: x ≤ 1, x ≥ 1, 0 ≤ y ≤ 1
    let mut c = unit_box(2);
    c.push(con(&[-1, 0], -1));
    assert_eq!(prepare(2, &c).err(), Some(PrepError::LowerDimensional));
    // empty
    let mut c = unit_box(2);
    c.push(con(&[1, 0], -1));
    assert_eq!(prepare(2, &c).err(), Some(PrepError::NoVertex));
}

#[test]
fn strips_redundant_constraints() {
    let mut c = unit_box(2);
    c.push(con(&[1, 1], 2)); // touches only the vertex (1,1)
    c.push(con(&[1, 0], 5)); // never tight
    c.push(con(&[1, 0], 1)); // duplicate
    let p = prepare(2, &c).unwrap();
    assert_eq!(p.cons.len(), 4);
    assert_eq!(p.removed_redundant, 2);
    for (c, s) in p.cons.iter().zip(&p.facet_starts) {
        let val: Q = c.a.iter().zip(&s.coords()).fold(Q::zero(), |acc, (a, x)| acc.add(&a.mul(x)));
        assert_eq!(val, c.b);
    }
}
