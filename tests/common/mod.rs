#![allow(dead_code)]
use nefvol::driver::{prepare_all, Options};
use nefvol::geom::{Problem, RawConstraint};
use nefvol::local::SweepParams;
use nefvol::num::{integerize, Point, Q, Z};
use nefvol::oracle::{boundary_vertices, oracle_sweep};
use nefvol::reverse::run_facet;
use nefvol::sweep::{Accumulator, SweepFunction};
use nefvol::vertex::Visit;
use std::collections::HashMap;

pub fn zv(v: &[i64]) -> Vec<Z> {
    v.iter().map(|&x| Z::from(x)).collect()
}
pub fn qv(v: &[i64]) -> Vec<Q> {
    v.iter().map(|&x| Q::int(x)).collect()
}

/// Oracle sweep function on raw polytopes (after preprocessing).
pub fn oracle_fn(d: usize, raw: &[Vec<RawConstraint>], a: &[Z]) -> (SweepFunction, usize) {
    let prep = prepare_all(d, raw).unwrap();
    let r = nefvol::driver::perturbation(d, 7);
    let (acc, _, emitted) = oracle_sweep(&prep.prob, &SweepParams { a, r: Some(&r) }).unwrap();
    (acc.finish().unwrap(), emitted)
}

/// Reverse search sweep function plus the multiset of emitted vertices.
pub fn reverse_fn(d: usize, raw: &[Vec<RawConstraint>], a: &[Z]) -> (SweepFunction, Vec<Point>, Problem) {
    let prep = prepare_all(d, raw).unwrap();
    let r = nefvol::driver::perturbation(d, 7);
    let params = SweepParams { a, r: Some(&r) };
    let mut acc = Accumulator::new(d);
    let mut emitted = Vec::new();
    for t in &prep.tasks {
        let mut sink: Vec<(Point, Visit)> = Vec::new();
        let (a2, _) = run_facet(&prep.prob, t, &params, &mut sink).unwrap();
        acc.merge(a2);
        emitted.extend(sink.into_iter().filter(|(_, v)| *v == Visit::Emitted).map(|(p, _)| p));
    }
    (acc.finish().unwrap(), emitted, prep.prob)
}

/// Assert the reverse search emits exactly the oracle's ∂U vertices, each once,
/// and the sweep functions agree.  Returns the total volume.
pub fn differential(d: usize, raw: &[Vec<RawConstraint>], a: &[i64]) -> Q {
    let a = zv(a);
    let (fo, n_o) = oracle_fn(d, raw, &a);
    let (fr, emitted, prob) = reverse_fn(d, raw, &a);
    let expected = boundary_vertices(&prob);
    let mut counts: HashMap<Point, usize> = HashMap::new();
    for p in &emitted {
        *counts.entry(p.clone()).or_default() += 1;
    }
    for p in &expected {
        assert_eq!(counts.get(p).copied().unwrap_or(0), 1, "vertex {p} emitted {} times", counts.get(p).copied().unwrap_or(0));
    }
    assert_eq!(emitted.len(), expected.len(), "emitted {} vertices, oracle has {}", emitted.len(), expected.len());
    assert_eq!(n_o, expected.len());
    assert_eq!(fr.knots.len(), fo.knots.len(), "knot counts differ");
    for ((k1, c1), (k2, c2)) in fr.knots.iter().zip(&fo.knots) {
        assert_eq!(k1, k2);
        assert_eq!(c1, c2, "coefficients differ at knot {k1}");
    }
    let t = fr.total().unwrap();
    assert_eq!(t, fo.total().unwrap());
    t
}

pub fn volume(d: usize, raw: &[Vec<RawConstraint>]) -> Q {
    nefvol::driver::volume(d, raw, &Options::default()).unwrap()
}

pub fn sweep_dir(d: usize, raw: &[Vec<RawConstraint>], dir: &[i64]) -> SweepFunction {
    let dq = qv(dir);
    let (f, _) = nefvol::driver::solve(d, raw, &dq, &Options::default()).unwrap();
    f
}

pub fn integer_dir(dir: &[Q]) -> Vec<Z> {
    integerize(dir)
}

/// Deterministic xorshift.
pub struct Rng(pub u64);
impl Rng {
    pub fn next(&mut self) -> u64 {
        let mut s = self.0;
        s ^= s << 13;
        s ^= s >> 7;
        s ^= s << 17;
        self.0 = s;
        s
    }
    pub fn int(&mut self, lo: i64, hi: i64) -> i64 {
        lo + (self.next() % ((hi - lo + 1) as u64)) as i64
    }
    pub fn unit(&mut self) -> f64 {
        (self.next() >> 11) as f64 / (1u64 << 53) as f64
    }
}

/// Monte Carlo estimate of vol(U ∩ {a·x ≤ λ}) inside a bounding box.
pub fn monte_carlo(d: usize, raw: &[Vec<RawConstraint>], lo: &[f64], hi: &[f64], a: Option<(&[f64], f64)>, n: usize, seed: u64) -> f64 {
    let mut rng = Rng(seed | 1);
    let mut hits = 0usize;
    let polys: Vec<Vec<(Vec<f64>, f64)>> = raw
        .iter()
        .map(|p| p.iter().map(|c| (c.a.iter().map(|x| x.to_f64()).collect(), c.b.to_f64())).collect())
        .collect();
    let mut x = vec![0.0; d];
    for _ in 0..n {
        for i in 0..d {
            x[i] = lo[i] + (hi[i] - lo[i]) * rng.unit();
        }
        if let Some((av, l)) = a {
            let s: f64 = av.iter().zip(&x).map(|(p, q)| p * q).sum();
            if s > l {
                continue;
            }
        }
        let inside = polys.iter().any(|p| p.iter().all(|(av, b)| av.iter().zip(&x).map(|(p, q)| p * q).sum::<f64>() <= *b));
        if inside {
            hits += 1;
        }
    }
    let vol_box: f64 = (0..d).map(|i| hi[i] - lo[i]).product();
    vol_box * hits as f64 / n as f64
}
