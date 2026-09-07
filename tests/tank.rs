mod common;

use common::*;
use nefvol::io::{parse_ine, Input};
use nefvol::num::Q;

fn tank() -> Input {
    parse_ine(include_str!("../inputs/fao_tank3d.ine")).unwrap()
}

fn q(s: &str) -> Q {
    Q::parse(s).unwrap()
}

fn contains(input: &Input, p: &[Q]) -> bool {
    input.polytopes.iter().any(|poly| {
        poly.iter().all(|c| {
            c.a.iter().zip(p).fold(Q::zero(), |s, (a, x)| s.add(&a.mul(x))) <= c.b
        })
    })
}

#[test]
fn tank_calibration_matches_integrated_horizontal_area() {
    let input = tank();
    let f = sweep_dir(input.d, &input.polytopes, &[0, 0, 1]);
    // Independent geometry: mean floor width 5 dm, length 7.2 dm,
    // widening by 3 dm over 8.44 dm height. Sump is 1 x 1.5 x 0.75 dm.
    for h in ["0", "1/10", "1", "211/50", "8", "211/25"] {
        let h = q(h);
        let expected = q("9/8").add(&Q::int(36).mul(&h))
            .add(&q("270/211").mul(&h).mul(&h));
        assert_eq!(f.eval(&h), expected, "level {h}");
    }
    assert_eq!(f.eval(&q("-1")), Q::zero());
    assert_eq!(f.eval(&q("-3/4")), Q::zero());
    assert_eq!(f.eval(&q("-3/8")), q("9/16"));
    assert_eq!(f.eval(&Q::int(9)), q("396117/1000"));
    assert_eq!(f.total().unwrap(), q("396117/1000"));
    assert_eq!(f.eval(&q("211/50")), q("175833/1000"));
}

#[test]
fn sump_makes_tank_and_viewer_section_nonconvex() {
    let input = tank();
    let section = parse_ine(include_str!("../inputs/fao_tank2d.ine")).unwrap();
    // Two interior points have a midpoint below the floor but outside the sump.
    let p = [q("5"), q("7/4"), q("1/10")];
    let r = [q("3/2"), q("7/4"), q("-7/10")];
    let m: Vec<_> = p.iter().zip(&r).map(|(a, b)| a.add(b).div(&Q::int(2))).collect();
    for (point, inside) in [(&p[..], true), (&r[..], true), (&m[..], false)] {
        assert_eq!(contains(&input, point), inside);
        assert_eq!(contains(&section, &[point[0].clone(), point[2].clone()]), inside);
    }
    // The entire halfspace section, not only the witness points, matches 3D.
    for (poly3, poly2) in input.polytopes.iter().zip(&section.polytopes) {
        let sliced: Vec<_> = poly3.iter().filter(|c| !c.a[0].is_zero() || !c.a[2].is_zero()).collect();
        assert_eq!(sliced.len(), poly2.len());
        for (c3, c2) in sliced.iter().zip(poly2) {
            assert_eq!(c2.a, vec![c3.a[0].clone(), c3.a[2].clone()]);
            assert_eq!(c2.b, c3.b.sub(&c3.a[1].mul(&q("7/4"))));
        }
    }
    let f = sweep_dir(2, &section.polytopes, &[0, 1]);
    let trapezoid = q("757/144").add(&q("1189/144")).mul(&q("211/25")).div(&Q::int(2));
    assert_eq!(f.total().unwrap(), trapezoid.add(&q("3/4")));
}

#[test]
fn tilted_tank_matches_oracle_and_preserves_capacity() {
    let input = tank();
    for direction in [[1, 0, 4], [-1, 0, 4], [1, 1, 4]] {
        assert_eq!(differential(3, &input.polytopes, &direction), q("396117/1000"));
    }
}
