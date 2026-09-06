mod common;
use common::*;
use nefvol::num::Q;
use nefvol::shapes::*;

#[test]
fn boxes_and_simplices_oracle_and_reverse() {
    for d in 1..=5 {
        let a: Vec<i64> = (0..d).map(|i| (2 * i + 3) as i64).collect();
        let v = differential(d, &[unit_box(d)], &a);
        assert_eq!(v, Q::int(1), "unit box d={d}");
        let v = differential(d, &[simplex(d, Q::int(1))], &a);
        assert_eq!(v, Q::one().div(&Q::int(nefvol::num::factorial(d))), "simplex d={d}");
        // scaled box with rational bounds
        let lo: Vec<Q> = (0..d).map(|i| Q::new(nefvol::num::Z::from(-(i as i64) - 1), nefvol::num::Z::from(2))).collect();
        let hi: Vec<Q> = (0..d).map(|i| Q::new(nefvol::num::Z::from(3 * i as i64 + 2), nefvol::num::Z::from(3))).collect();
        let expect = lo.iter().zip(&hi).fold(Q::one(), |acc, (l, h)| acc.mul(&h.sub(l)));
        let v = differential(d, &[boxp(&lo, &hi)], &a);
        assert_eq!(v, expect, "rational box d={d}");
    }
}

#[test]
fn sweep_function_of_unit_cube_matches_closed_form() {
    // cube [0,1]^3, a = (1,1,1): s(λ) = λ³/6 on [0,1]
    let f = sweep_dir(3, &[unit_box(3)], &[1, 1, 1]);
    assert_eq!(f.eval(&Q::parse("1/2").unwrap()), Q::parse("1/48").unwrap());
    // at λ=3/2: 1/2 by symmetry
    assert_eq!(f.eval(&Q::parse("3/2").unwrap()), Q::parse("1/2").unwrap());
    assert_eq!(f.total().unwrap(), Q::int(1));
    // constant tail & zero head are enforced by total(); check zero below
    assert_eq!(f.eval(&Q::int(-1)), Q::zero());
}

#[test]
fn degenerate_vertices() {
    // square pyramid: apex (0,0,1) with 4 facets → degenerate vertex
    let pyr = vec![
        con(&[0, 0, -1], 0),
        con(&[1, 0, 1], 1),
        con(&[-1, 0, 1], 1),
        con(&[0, 1, 1], 1),
        con(&[0, -1, 1], 1),
    ];
    // base [-1,1]^2 area 4, height 1 → 4/3
    assert_eq!(differential(3, std::slice::from_ref(&pyr), &[3, 5, 7]), Q::parse("4/3").unwrap());
    // octahedron |x|+|y|+|z| ≤ 1 → vol 4/3, every vertex on 4 facets
    let mut octa = Vec::new();
    for sx in [-1i64, 1] {
        for sy in [-1i64, 1] {
            for sz in [-1i64, 1] {
                octa.push(con(&[sx, sy, sz], 1));
            }
        }
    }
    assert_eq!(differential(3, &[octa.clone()], &[3, 5, 7]), Q::parse("4/3").unwrap());
    // 4d cross-polytope: 2^4/4! = 2/3
    let mut cross4 = Vec::new();
    for m in 0..16i64 {
        let s: Vec<i64> = (0..4).map(|i| if (m >> i) & 1 == 1 { -1 } else { 1 }).collect();
        cross4.push(con(&s, 1));
    }
    assert_eq!(differential(4, &[cross4], &[3, 5, 7, 11]), Q::parse("2/3").unwrap());
    // cone over a 24-gon-like base: 24 planes through the apex (0,0,2)
    let mut cone = vec![con(&[0, 0, -1], 0)];
    let n = 24i64;
    for k in 0..n {
        // integer directions on a circle-ish: use (cos,sin) approximations via Pythagorean-ish pairs
        let ang = 2.0 * std::f64::consts::PI * k as f64 / n as f64;
        let cx = (ang.cos() * 1000.0).round() as i64;
        let cy = (ang.sin() * 1000.0).round() as i64;
        // plane through apex (0,0,2): cx x + cy y + 500 z ≤ 1000  (at z=0 radius 2)
        cone.push(con(&[cx, cy, 500], 1000));
    }
    let v = differential(3, &[cone.clone()], &[3, 5, 7]);
    // base is the 24-gon circumscribed about the unit circle: area 24·tan(π/24), height 2
    let vf = v.to_f64();
    let expect = 24.0 * (std::f64::consts::PI / 24.0).tan() * 2.0 / 3.0;
    assert!((vf - expect).abs() < 1e-3, "cone volume {vf} vs {expect}");
    // pyramid ∪ octahedron (overlapping, degenerate) vs Monte Carlo
    let v = differential(3, &[pyr, octa], &[3, 5, 7]);
    let mc = monte_carlo(3, &[
        vec![con(&[0, 0, -1], 0), con(&[1, 0, 1], 1), con(&[-1, 0, 1], 1), con(&[0, 1, 1], 1), con(&[0, -1, 1], 1)],
        {
            let mut o = Vec::new();
            for sx in [-1i64, 1] { for sy in [-1i64, 1] { for sz in [-1i64, 1] { o.push(con(&[sx, sy, sz], 1)); } } }
            o
        },
    ], &[-1.0; 3], &[1.0; 3], None, 400_000, 3);
    assert!((v.to_f64() - mc).abs() < 0.02, "exact {} vs MC {mc}", v.to_f64());
}

#[test]
fn unions_golden() {
    // L-shape: [0,2]×[0,1] ∪ [0,1]×[0,2] → 3
    let l = vec![int_box(&[0, 0], &[2, 1]), int_box(&[0, 0], &[1, 2])];
    assert_eq!(differential(2, &l, &[3, 5]), Q::int(3));
    // nested: [0,3]^2 ⊇ [1,2]^2 → 9
    let nested = vec![int_box(&[0, 0], &[3, 3]), int_box(&[1, 1], &[2, 2])];
    assert_eq!(differential(2, &nested, &[3, 5]), Q::int(9));
    // disjoint: two unit squares → 2
    let disj = vec![int_box(&[0, 0], &[1, 1]), int_box(&[5, 5], &[6, 6])];
    assert_eq!(differential(2, &disj, &[3, 5]), Q::int(2));
    // adjacent (shared facet): [0,1]² ∪ [1,2]×[0,1] → 2
    let adj = vec![int_box(&[0, 0], &[1, 1]), int_box(&[1, 0], &[2, 1])];
    assert_eq!(differential(2, &adj, &[3, 5]), Q::int(2));
    // three boxes overlapping in 3d
    let three = vec![int_box(&[0, 0, 0], &[2, 2, 2]), int_box(&[1, 1, 1], &[3, 3, 3]), int_box(&[0, 2, 0], &[1, 4, 1])];
    // 8 + 8 − 1 + 2 = 17  (third box: [0,1]×[2,4]×[0,1] vol 2, overlaps first on the face y=2 only)
    assert_eq!(differential(3, &three, &[3, 5, 7]), Q::int(17));
    // union with a shared vertex only
    let touch = vec![int_box(&[0, 0], &[1, 1]), int_box(&[1, 1], &[2, 2])];
    assert_eq!(differential(2, &touch, &[3, 5]), Q::int(2));
    // four squares around the origin (interior vertex of the union)
    let four = vec![int_box(&[0, 0], &[1, 1]), int_box(&[-1, 0], &[0, 1]), int_box(&[-1, -1], &[0, 0]), int_box(&[0, -1], &[1, 0])];
    assert_eq!(differential(2, &four, &[3, 5]), Q::int(4));
    // three triangles meeting at the origin, six sectors each pair-covered (hole under perturbation)
    let tri = vec![
        vec![con(&[-1, 0], 0), con(&[0, -1], 0), con(&[1, 1], 1)],
        vec![con(&[1, 0], 0), con(&[0, -1], 0), con(&[-1, 1], 1)],
        vec![con(&[0, 1], 0), con(&[1, -1], 1), con(&[-1, -1], 1)],
    ];
    // areas: 1/2 + 1/2 + 1 = 2 (third: y ≤ 0, x−y ≤ 1, −x−y ≤ 1 : triangle (−1,0),(1,0),(0,−1) area 1)
    assert_eq!(differential(2, &tri, &[3, 5]), Q::int(2));
}

#[test]
fn inclusion_exclusion_on_random_pairs() {
    let mut rng = Rng(99);
    let mut checked = 0;
    for _ in 0..30 {
        let d = 2 + (rng.next() % 2) as usize;
        let mk = |rng: &mut Rng| {
            // random simplex-ish polytope: box ∩ 2 random halfspaces, then check it's full-dim via prepare
            let mut c = int_box(&vec![0; d], &vec![4; d]);
            for _ in 0..2 {
                let a: Vec<i64> = (0..d).map(|_| rng.int(-3, 3)).collect();
                if a.iter().all(|&x| x == 0) {
                    continue;
                }
                c.push(con(&a, rng.int(1, 8)));
            }
            let t: Vec<i64> = (0..d).map(|_| rng.int(-2, 2)).collect();
            // translate: a·(x − t) ≤ b
            c.iter().map(|rc| {
                let shift: Q = rc.a.iter().zip(&t).fold(Q::zero(), |s, (ai, ti)| s.add(&ai.mul(&Q::int(*ti))));
                nefvol::geom::RawConstraint { a: rc.a.clone(), b: rc.b.add(&shift) }
            }).collect::<Vec<_>>()
        };
        let p = mk(&mut rng);
        let q = mk(&mut rng);
        let ok = |x: &Vec<nefvol::geom::RawConstraint>| nefvol::preprocess::prepare(d, x).is_ok();
        if !ok(&p) || !ok(&q) {
            continue;
        }
        let vp = volume(d, std::slice::from_ref(&p));
        let vq = volume(d, std::slice::from_ref(&q));
        let vu = volume(d, &[p.clone(), q.clone()]);
        let mut inter = p.clone();
        inter.extend(q.clone());
        let vi = match nefvol::preprocess::prepare(d, &inter) {
            Ok(_) => volume(d, &[inter]),
            Err(_) => Q::zero(),
        };
        assert_eq!(vu.add(&vi), vp.add(&vq), "inclusion-exclusion failed");
        checked += 1;
    }
    assert!(checked >= 10, "only {checked} pairs checked");
}

#[test]
fn direction_independence_and_random_unions_vs_monte_carlo() {
    let mut rng = Rng(2024);
    for trial in 0..6 {
        let d = 2 + trial % 2;
        // rotated boxes: unimodular shear-like maps keep things integral
        let mut polys = Vec::new();
        for _ in 0..3 {
            let lo: Vec<Q> = (0..d).map(|_| Q::int(rng.int(-3, 1))).collect();
            let hi: Vec<Q> = lo.iter().map(|l| l.add(&Q::int(rng.int(1, 4)))).collect();
            let b = boxp(&lo, &hi);
            // inverse matrix of a random unit-upper-triangular-ish map: M^{-1} = I - s E_ij
            let mut minv: Vec<Vec<Q>> = (0..d).map(|i| (0..d).map(|j| if i == j { Q::one() } else { Q::zero() }).collect()).collect();
            let i = (rng.next() % d as u64) as usize;
            let j = (rng.next() % d as u64) as usize;
            if i != j {
                minv[i][j] = Q::int(rng.int(-2, 2));
            }
            let t: Vec<Q> = (0..d).map(|_| Q::int(rng.int(-1, 1))).collect();
            polys.push(transform(&b, &minv, &t));
        }
        let mut vols = Vec::new();
        for k in 0..4 {
            let dir: Vec<i64> = (0..d).map(|i| rng.int(1, 50) * (if (i + k) % 2 == 0 { 1 } else { -1 })).collect();
            let f = sweep_dir(d, &polys, &dir);
            vols.push(f.total().unwrap());
        }
        for v in &vols[1..] {
            assert_eq!(v, &vols[0], "direction dependence");
        }
        // also the differential test against the oracle
        let v = differential(d, &polys, &[7, 11, 13][..d]);
        assert_eq!(v, vols[0]);
        let mc = monte_carlo(d, &polys, &vec![-12.0; d], &vec![12.0; d], None, 300_000, 17 + trial as u64);
        let ex = v.to_f64();
        let tol = 0.05 * ex.max(1.0) + (24f64.powi(d as i32) / 300_000f64).sqrt() * 2.0;
        assert!((ex - mc).abs() < tol, "trial {trial}: exact {ex} vs MC {mc}");
    }
}

#[test]
fn sweep_partial_volumes_vs_monte_carlo() {
    // L-shape in 3d, generic direction; check s(λ) at several λ
    let polys = vec![int_box(&[0, 0, 0], &[3, 1, 1]), int_box(&[0, 0, 0], &[1, 3, 1]), int_box(&[0, 0, 0], &[1, 1, 3])];
    let dir = [2i64, 3, 5];
    let f = sweep_dir(3, &polys, &dir);
    assert_eq!(f.total().unwrap(), Q::int(7));
    let af: Vec<f64> = dir.iter().map(|&x| x as f64).collect();
    for l in [2.0f64, 4.5, 7.0, 9.5] {
        let ex = f.eval(&Q::parse(&l.to_string()).unwrap()).to_f64();
        let mc = monte_carlo(3, &polys, &[0.0; 3], &[3.0; 3], Some((&af, l)), 400_000, 5);
        assert!((ex - mc).abs() < 0.05, "λ={l}: exact {ex} vs MC {mc}");
    }
}

#[test]
fn budget_splitting_is_exact_and_deterministic() {
    use nefvol::driver::{solve, Options};
    let polys = vec![int_box(&[0, 0, 0], &[3, 1, 1]), int_box(&[0, 0, 0], &[1, 3, 1]), int_box(&[0, 0, 0], &[1, 1, 3]), simplex(3, Q::int(2))];
    let dir = qv(&[3, 5, 7]);
    let (f0, r0) = solve(3, &polys, &dir, &Options { budget: None, ..Default::default() }).unwrap();
    for b in [1u64, 2, 3, 7, 50] {
        let (f, r) = solve(3, &polys, &dir, &Options { budget: Some(b), ..Default::default() }).unwrap();
        assert_eq!(f.knots, f0.knots, "budget {b}");
        assert_eq!(r.search.nodes, r0.search.nodes, "budget {b}: node count");
        assert_eq!(r.search.emitted, r0.search.emitted, "budget {b}: emitted");
        if b <= 7 {
            assert!(r.search.tasks > r0.search.tasks, "budget {b} produced no split");
        }
    }
}
