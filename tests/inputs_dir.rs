//! Every small shipped input: reverse search == oracle, and the .ine/JSON
//! parsers round-trip through `to_ine`.
use nefvol::driver::{solve, Options};
use nefvol::io::{parse_auto, to_ine};

fn check(path: &str) {
    let text = std::fs::read_to_string(path).unwrap();
    let inp = parse_auto(&text).unwrap();
    let dir = nefvol::driver::default_direction(inp.d, 5);
    let (f, r) = solve(inp.d, &inp.polytopes, &dir, &Options::default()).unwrap();
    let (fo, _) = solve(inp.d, &inp.polytopes, &dir, &Options { oracle: true, ..Default::default() }).unwrap();
    assert_eq!(f.knots, fo.knots, "{path}: reverse search differs from oracle");
    assert!(r.search.emitted > 0);
    // round trip through the .ine writer
    let again = parse_auto(&to_ine(inp.d, &inp.polytopes)).unwrap();
    let (f2, _) = solve(again.d, &again.polytopes, &dir, &Options::default()).unwrap();
    assert_eq!(f2.total().unwrap(), f.total().unwrap(), "{path}: round trip");
}

#[test]
fn shipped_inputs_agree_with_oracle() {
    for f in ["square.ine", "lshape3d.ine", "union.json", "cones3d_30x3.ine", "grid3d_4.ine", "boxes3d_20.ine", "rand4d_8.ine"] {
        check(&format!("{}/inputs/{f}", env!("CARGO_MANIFEST_DIR")));
    }
}

#[test]
fn ine_errors_are_reported() {
    assert!(parse_auto("begin\n2 3 integer\n1 0\nend\n").is_err());
    assert!(parse_auto("linearity 1 1\nbegin\n1 3 integer\n1 0 0\nend\n").is_err());
    assert!(parse_auto("{\"polytopes\": [{\"A\": [[1,0]], \"b\": [1, 2]}]}").is_err());
}
