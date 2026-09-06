use nefvol::geom::Problem;
use nefvol::local::SweepParams;
use nefvol::num::{Q, Z};
use nefvol::oracle::oracle_sweep;
use nefvol::shapes::*;

fn zv(v: &[i64]) -> Vec<Z> {
    v.iter().map(|&x| Z::from(x)).collect()
}

#[test]
fn spec_3_3_unit_square() {
    let prob = Problem::from_raw(2, &[unit_box(2)]).unwrap();
    let a = zv(&[1, 2]);
    let (acc, stats, emitted) = oracle_sweep(&prob, &SweepParams { a: &a, r: None }).unwrap();
    assert_eq!(emitted, 4);
    assert_eq!(stats.bases, 4);
    let f = acc.finish().unwrap();
    let knots: Vec<(String, String)> = f.knots.iter().map(|(k, c)| (k.to_string(), c[2].to_string())).collect();
    assert_eq!(
        knots,
        vec![
            ("0".to_string(), "1/4".to_string()),
            ("1".to_string(), "-1/4".to_string()),
            ("2".to_string(), "-1/4".to_string()),
            ("3".to_string(), "1/4".to_string())
        ]
    );
    assert_eq!(f.eval(&Q::int(5)), Q::int(1));
    assert_eq!(f.eval(&Q::int(1)), Q::parse("1/4").unwrap());
    assert_eq!(f.total().unwrap(), Q::int(1));
}
