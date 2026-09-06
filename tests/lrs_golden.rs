//! Per-piece volumes recorded from `lrs` (lrslib-073, `scripts/lrs_check.py`,
//! H → V enumeration followed by `volume` on the V-representation).  Every
//! convex piece of every shipped instance must match exactly, so the suite
//! needs no lrs installation.  Regenerate `tests/lrs_golden.txt` with
//! `scripts/lrs_check.py --golden` if instances change.
use nefvol::driver::{solve, Options};
use nefvol::io::parse_auto;
use nefvol::num::Q;
use std::collections::BTreeMap;

#[test]
fn every_convex_piece_matches_lrs() {
    let root = env!("CARGO_MANIFEST_DIR");
    let golden = std::fs::read_to_string(format!("{root}/tests/lrs_golden.txt")).unwrap();
    let mut by_instance: BTreeMap<String, Vec<(usize, Q)>> = BTreeMap::new();
    for line in golden.lines().filter(|l| !l.trim().is_empty()) {
        let t: Vec<&str> = line.split_whitespace().collect();
        by_instance.entry(t[0].to_string()).or_default().push((t[1].parse().unwrap(), Q::parse(t[2]).unwrap()));
    }
    assert!(by_instance.len() >= 10);
    let mut checked = 0;
    for (name, pieces) in &by_instance {
        let text = std::fs::read_to_string(format!("{root}/inputs/{name}.ine")).unwrap();
        let inp = parse_auto(&text).unwrap();
        assert_eq!(inp.polytopes.len(), pieces.len(), "{name}: piece count");
        for (k, expected) in pieces {
            let piece = std::slice::from_ref(&inp.polytopes[*k]);
            let dir = nefvol::driver::default_direction(inp.d, 1 + *k as u64);
            let (f, _) = solve(inp.d, piece, &dir, &Options::default()).unwrap();
            assert_eq!(&f.total().unwrap(), expected, "{name} piece {k}");
            checked += 1;
        }
    }
    assert_eq!(checked, 276);
}
