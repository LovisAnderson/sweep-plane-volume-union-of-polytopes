mod common;
use common::*;
use nefvol::geom::RawConstraint;
use nefvol::num::Q;
use nefvol::shapes::*;

/// Coefficients far beyond the i64 fast-path bounds: everything runs through
/// the generic Z/BigInt paths and must still be exact.
#[test]
fn huge_coefficients_are_exact() {
    let big = Q::parse("123456789012345678901234567890").unwrap(); // ~1.2e29 > 2^62
    let lo = vec![Q::zero(), Q::zero(), Q::zero()];
    let hi = vec![big.clone(), big.clone(), big.clone()];
    let b = boxp(&lo, &hi);
    // shear x -> x + 3y so normals are (1,-3,0)-like with the huge rhs
    let minv = vec![qv(&[1, 0, 0]), qv(&[-3, 1, 0]), qv(&[0, 0, 1])];
    let sheared: Vec<RawConstraint> = transform(&b, &minv, &qv(&[0, 0, 0]));
    let v = volume(3, std::slice::from_ref(&sheared));
    assert_eq!(v, big.pow(3));
    // union with a copy translated by half the extent along y (the shear is y -> 3x + y,
    // so y-slices of the copy overlap the original half-way)
    let shift = big.div(&Q::int(2));
    let moved: Vec<RawConstraint> = sheared.iter().map(|c| RawConstraint { a: c.a.clone(), b: c.b.add(&c.a[1].mul(&shift)) }).collect();
    let vu = volume(3, &[sheared, moved]);
    assert_eq!(vu, big.pow(3).mul(&Q::parse("3/2").unwrap()));
}

#[test]
fn rational_direction_and_tiny_polytope() {
    // sweep direction with denominators; polytope with tiny rational extent
    let eps = Q::parse("1/1000000007").unwrap();
    let cube = boxp(&vec![Q::zero(); 3], &vec![eps.clone(); 3]);
    let dir = vec![Q::parse("1/3").unwrap(), Q::parse("2/7").unwrap(), Q::parse("5/11").unwrap()];
    let (f, _) = nefvol::driver::solve(3, &[cube], &dir, &Default::default()).unwrap();
    assert_eq!(f.total().unwrap(), eps.pow(3));
    // half-way knot value equals half the volume by central symmetry of the cube
    let mid = dir.iter().fold(Q::zero(), |s, a| s.add(&a.mul(&eps))).div(&Q::int(2));
    assert_eq!(f.eval(&mid), eps.pow(3).div(&Q::int(2)));
}
