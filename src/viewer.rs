//! Presentation data is separate from solving and from dimension-specific rendering.
use crate::{geom::RawConstraint, io::Input, num::Q, sweep::SweepFunction};
use serde_json::{json, Value};

/// Intersect pairs of supporting lines, retaining exactly feasible vertices.
/// This is geometry for display; it does not run the volume algorithm again.
fn polygon(cons: &[RawConstraint]) -> Vec<Vec<Q>> {
    let mut points = Vec::new();
    for (i, c) in cons.iter().enumerate() {
        for e in &cons[..i] {
            let det = c.a[0].mul(&e.a[1]).sub(&c.a[1].mul(&e.a[0]));
            if det.is_zero() { continue; }
            let p = vec![c.b.mul(&e.a[1]).sub(&c.a[1].mul(&e.b)).div(&det),
                         c.a[0].mul(&e.b).sub(&c.b.mul(&e.a[0])).div(&det)];
            if cons.iter().all(|r| r.a[0].mul(&p[0]).add(&r.a[1].mul(&p[1])) <= r.b) && !points.contains(&p) {
                points.push(p);
            }
        }
    }
    // Exact monotone-chain hull gives cyclic order without floating-point angles.
    points.sort();
    fn cross(a: &[Q], b: &[Q], c: &[Q]) -> Q {
        b[0].sub(&a[0]).mul(&c[1].sub(&a[1])).sub(&b[1].sub(&a[1]).mul(&c[0].sub(&a[0])))
    }
    let mut hull: Vec<Vec<Q>> = Vec::new();
    for p in &points {
        while hull.len() >= 2 && cross(&hull[hull.len()-2], &hull[hull.len()-1], p) <= Q::zero() { hull.pop(); }
        hull.push(p.clone());
    }
    let lower = hull.len();
    for p in points.iter().rev().skip(1) {
        while hull.len() > lower && cross(&hull[hull.len()-2], &hull[hull.len()-1], p) <= Q::zero() { hull.pop(); }
        hull.push(p.clone());
    }
    hull.pop();
    hull
}

pub fn scene(input: &Input, direction: &[Q], f: &SweepFunction) -> Result<Value, String> {
    if input.d != 2 { return Err("the viewer currently supports only 2D inputs".into()); }
    if direction.len() != 2 { return Err("a 2D sweep direction needs two coordinates".into()); }
    let polygons: Vec<_> = input.polytopes.iter().map(|c| polygon(c)).collect();
    let mut projections: Vec<Q> = polygons.iter().flatten().map(|p| direction[0].mul(&p[0]).add(&direction[1].mul(&p[1]))).collect();
    projections.sort();
    let lo = projections.first().ok_or("no vertices to display")?;
    let hi = projections.last().unwrap();
    let strings = |v: &[Q]| v.iter().map(ToString::to_string).collect::<Vec<_>>();
    let geometry: Vec<_> = polygons.iter().map(|p| json!({"vertices": p.iter().map(|v| strings(v)).collect::<Vec<_>>()})).collect();
    let pieces: Vec<_> = f.pieces().iter().map(|p| json!({"lo": p.lo.as_ref().map(ToString::to_string), "hi": p.hi.as_ref().map(ToString::to_string), "coefficients": strings(&p.poly)})).collect();
    Ok(json!({"dimension": input.d, "direction": strings(direction), "polytopes": geometry,
        "lambdaRange": [lo.to_string(), hi.to_string()], "pieces": pieces, "total": f.total()?.to_string()}))
}
