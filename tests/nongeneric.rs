mod common;
use common::*;
use nefvol::driver::{solve, Options};
use nefvol::num::Q;
use nefvol::shapes::*;

#[test]
fn axis_direction_on_box_gives_linear_pieces() {
    // a = (1,0): s(λ) = λ on [0,1] for the unit square; every edge parallel to y is orthogonal to a
    let f = sweep_dir(2, &[unit_box(2)], &[1, 0]);
    assert_eq!(f.eval(&Q::parse("1/2").unwrap()), Q::parse("1/2").unwrap());
    assert_eq!(f.total().unwrap(), Q::int(1));
    let pieces = f.pieces();
    // piece on [0,1): λ
    assert_eq!(pieces[1].poly, vec![Q::zero(), Q::one(), Q::zero()]);
    // 3d: a = (1,1,0) on the unit cube: s(λ)= λ²/2 on [0,1]
    let f = sweep_dir(3, &[unit_box(3)], &[1, 1, 0]);
    assert_eq!(f.eval(&Q::parse("1/2").unwrap()), Q::parse("1/8").unwrap());
    assert_eq!(f.eval(&Q::int(1)), Q::parse("1/2").unwrap());
    assert_eq!(f.total().unwrap(), Q::int(1));
}

#[test]
fn nongeneric_union_matches_generic_total_and_monte_carlo() {
    let polys = vec![int_box(&[0, 0, 0], &[2, 2, 1]), int_box(&[1, 1, 0], &[3, 3, 2]), simplex(3, Q::int(2))];
    let generic = sweep_dir(3, &polys, &[3, 5, 7]).total().unwrap();
    for dir in [[1i64, 0, 0], [0, 1, 0], [1, 1, 0], [1, -1, 0], [0, 0, 1], [1, 1, 1]] {
        let f = sweep_dir(3, &polys, &dir);
        assert_eq!(f.total().unwrap(), generic, "dir {dir:?}");
        let af: Vec<f64> = dir.iter().map(|&x| x as f64).collect();
        for l in [0.5f64, 1.5, 2.5] {
            let ex = f.eval(&Q::parse(&l.to_string()).unwrap()).to_f64();
            let mc = monte_carlo(3, &polys, &[0.0; 3], &[3.0; 3], Some((&af, l)), 300_000, 11);
            assert!((ex - mc).abs() < 0.06, "dir {dir:?} λ={l}: exact {ex} vs MC {mc}");
        }
    }
}

#[test]
fn assume_generic_fails_loudly() {
    let opts = Options { assume_generic: true, ..Default::default() };
    let r = solve(2, &[unit_box(2)], &qv(&[1, 0]), &opts);
    assert!(r.is_err());
    let msg = format!("{}", r.err().unwrap());
    assert!(msg.contains("orthogonal"), "{msg}");
}

#[test]
fn rational_direction_is_reported_on_the_users_scale() {
    // a = (1/2, 1/3) on the unit square vs a = (3, 2): s_{a}(λ) = s_{(3,2)}(6λ)
    let dir = vec![Q::parse("1/2").unwrap(), Q::parse("1/3").unwrap()];
    let (f, _) = solve(2, &[unit_box(2)], &dir, &Options::default()).unwrap();
    let g = sweep_dir(2, &[unit_box(2)], &[3, 2]);
    for l in ["1/12", "1/4", "1/2", "2/3", "5/6"] {
        let l = Q::parse(l).unwrap();
        assert_eq!(f.eval(&l), g.eval(&l.mul(&Q::int(6))), "λ = {l}");
    }
    assert_eq!(f.max_knot().unwrap(), &Q::parse("5/6").unwrap());
}
