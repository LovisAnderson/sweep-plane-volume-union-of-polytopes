//! Input formats: lrs/cdd `.ine` H-representations (several `begin…end`
//! blocks = union) and a small JSON schema.

use crate::geom::RawConstraint;
use crate::num::Q;

pub struct Input {
    pub d: usize,
    pub polytopes: Vec<Vec<RawConstraint>>,
}

pub fn parse_ine(text: &str) -> Result<Input, String> {
    let mut lines = text.lines().map(|l| l.trim()).filter(|l| !l.is_empty() && !l.starts_with('*')).peekable();
    let mut polytopes = Vec::new();
    let mut d: Option<usize> = None;
    while let Some(line) = lines.next() {
        let lower = line.to_ascii_lowercase();
        if lower.starts_with("linearity") {
            return Err("linearity rows are not supported (polytopes must be full-dimensional)".into());
        }
        if lower.starts_with("v-representation") {
            return Err("V-representations are not supported; give an H-representation".into());
        }
        if lower != "begin" {
            continue;
        }
        let header = lines.next().ok_or("unexpected end after 'begin'")?;
        let toks: Vec<&str> = header.split_whitespace().collect();
        if toks.len() < 2 {
            return Err(format!("bad header line '{header}'"));
        }
        let m: usize = toks[0].parse().map_err(|_| format!("bad row count '{}'", toks[0]))?;
        let n: usize = toks[1].parse().map_err(|_| format!("bad column count '{}'", toks[1]))?;
        if n < 2 {
            return Err("need at least 2 columns (b and one coordinate)".into());
        }
        let dd = n - 1;
        if let Some(prev) = d {
            if prev != dd {
                return Err(format!("dimension mismatch: {prev} vs {dd}"));
            }
        }
        d = Some(dd);
        let mut nums: Vec<Q> = Vec::with_capacity(m * n);
        while nums.len() < m * n {
            let l = lines.next().ok_or("unexpected end inside block")?;
            if l.eq_ignore_ascii_case("end") {
                return Err(format!("block ended after {} of {} numbers", nums.len(), m * n));
            }
            for t in l.split_whitespace() {
                nums.push(Q::parse(t).ok_or_else(|| format!("bad number '{t}'"))?);
            }
        }
        if nums.len() != m * n {
            return Err("row length mismatch".into());
        }
        match lines.next() {
            Some(l) if l.eq_ignore_ascii_case("end") => {}
            _ => return Err("expected 'end'".into()),
        }
        let mut cons = Vec::with_capacity(m);
        for i in 0..m {
            let row = &nums[i * n..(i + 1) * n];
            // b + a·x ≥ 0  ⇔  (−a)·x ≤ b
            cons.push(RawConstraint { a: row[1..].iter().map(|x| x.neg()).collect(), b: row[0].clone() });
        }
        polytopes.push(cons);
    }
    let d = d.ok_or("no 'begin' block found")?;
    Ok(Input { d, polytopes })
}

fn json_q(v: &serde_json::Value) -> Result<Q, String> {
    match v {
        serde_json::Value::Number(n) => Q::parse(&n.to_string()).ok_or_else(|| format!("bad number {n}")),
        serde_json::Value::String(s) => Q::parse(s).ok_or_else(|| format!("bad number '{s}'")),
        _ => Err(format!("expected a number, got {v}")),
    }
}

/// `{"polytopes": [{"A": [[...]], "b": [...]}, ...]}` meaning `A x ≤ b`.
pub fn parse_json(text: &str) -> Result<Input, String> {
    let v: serde_json::Value = serde_json::from_str(text).map_err(|e| e.to_string())?;
    let polys = v.get("polytopes").and_then(|p| p.as_array()).ok_or("expected top-level 'polytopes' array")?;
    let mut d: Option<usize> = None;
    let mut polytopes = Vec::new();
    for (pi, p) in polys.iter().enumerate() {
        let a = p.get("A").and_then(|x| x.as_array()).ok_or(format!("polytope {pi}: missing 'A'"))?;
        let b = p.get("b").and_then(|x| x.as_array()).ok_or(format!("polytope {pi}: missing 'b'"))?;
        if a.len() != b.len() {
            return Err(format!("polytope {pi}: A has {} rows but b has {}", a.len(), b.len()));
        }
        let mut cons = Vec::new();
        for (ri, (row, bb)) in a.iter().zip(b).enumerate() {
            let row = row.as_array().ok_or(format!("polytope {pi} row {ri}: expected array"))?;
            let av: Vec<Q> = row.iter().map(json_q).collect::<Result<_, _>>()?;
            match d {
                None => d = Some(av.len()),
                Some(dd) if dd != av.len() => return Err(format!("polytope {pi} row {ri}: dimension mismatch")),
                _ => {}
            }
            cons.push(RawConstraint { a: av, b: json_q(bb)? });
        }
        polytopes.push(cons);
    }
    Ok(Input { d: d.ok_or("no constraints")?, polytopes })
}

pub fn parse_auto(text: &str) -> Result<Input, String> {
    if text.trim_start().starts_with('{') {
        parse_json(text)
    } else {
        parse_ine(text)
    }
}

pub fn read_inputs(paths: &[String]) -> Result<Input, String> {
    let mut d: Option<usize> = None;
    let mut polytopes = Vec::new();
    for p in paths {
        let text = std::fs::read_to_string(p).map_err(|e| format!("{p}: {e}"))?;
        let inp = parse_auto(&text).map_err(|e| format!("{p}: {e}"))?;
        if let Some(dd) = d {
            if dd != inp.d {
                return Err(format!("{p}: dimension {} does not match {dd}", inp.d));
            }
        }
        d = Some(inp.d);
        polytopes.extend(inp.polytopes);
    }
    Ok(Input { d: d.ok_or("no input")?, polytopes })
}

/// Write an `.ine` file with one block per polytope.
pub fn to_ine(d: usize, polys: &[Vec<RawConstraint>]) -> String {
    let mut s = String::new();
    for (i, p) in polys.iter().enumerate() {
        s.push_str(&format!("polytope_{i}\nH-representation\nbegin\n{} {} rational\n", p.len(), d + 1));
        for c in p {
            s.push_str(&format!("{}", c.b));
            for a in &c.a {
                s.push_str(&format!(" {}", a.neg()));
            }
            s.push('\n');
        }
        s.push_str("end\n");
    }
    s
}
