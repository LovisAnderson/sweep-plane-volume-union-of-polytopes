use nefvol::{driver::{solve, Options}, io::parse_auto, num::Q, viewer::scene};

#[test]
fn geometry_and_cached_function_agree() {
    let input = parse_auto(include_str!("../inputs/square.ine")).unwrap();
    for direction in [vec![Q::int(1), Q::int(2)], vec![Q::int(-1), Q::zero()], vec![Q::parse("1/2").unwrap(), Q::int(1)]] {
        let (f, _) = solve(input.d, &input.polytopes, &direction, &Options::default()).unwrap();
        let data = scene(&input, &direction, &f).unwrap();
        assert_eq!(data["polytopes"][0]["vertices"].as_array().unwrap().len(), 4);
        let lo = Q::parse(data["lambdaRange"][0].as_str().unwrap()).unwrap();
        let hi = Q::parse(data["lambdaRange"][1].as_str().unwrap()).unwrap();
        assert_eq!(f.eval(&lo), Q::zero());
        assert_eq!(f.eval(&hi), Q::one());
        assert_eq!(data["total"], "1");
    }
}

#[test]
fn overlapping_squares_are_a_union() {
    let input = parse_auto(r#"{"polytopes":[{"A":[[-1,0],[1,0],[0,-1],[0,1]],"b":[0,2,0,2]},{"A":[[-1,0],[1,0],[0,-1],[0,1]],"b":[-1,3,-1,3]}]}"#).unwrap();
    let direction = vec![Q::int(1), Q::int(1)];
    let (f, _) = solve(2, &input.polytopes, &direction, &Options::default()).unwrap();
    let data = scene(&input, &direction, &f).unwrap();
    assert_eq!(data["total"], "7");
    assert_eq!(data["lambdaRange"], serde_json::json!(["0", "6"]));
    assert_eq!(data["polytopes"].as_array().unwrap().len(), 2);
}
