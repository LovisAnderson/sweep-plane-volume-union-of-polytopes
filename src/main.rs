use clap::{Args, Parser, Subcommand};
use nefvol::driver::{self, Options};
use nefvol::io::read_inputs;
use nefvol::num::Q;
use nefvol::sweep::SweepFunction;
use std::process::ExitCode;

#[derive(Parser)]
#[command(name = "nefvol", about = "Exact sweep-plane volume of a union of polytopes (reverse search + Bieri–Nef)")]
struct Cli {
    #[command(subcommand)]
    cmd: Cmd,
}

#[derive(Args, Clone)]
struct Common {
    /// Input file(s): lrs/cdd .ine (several begin…end blocks = union) or JSON
    #[arg(long, required = true, num_args = 1..)]
    input: Vec<String>,
    /// Sweep direction a1,a2,... (rationals); default: pseudo-random generic
    #[arg(long, value_delimiter = ',', allow_hyphen_values = true)]
    direction: Option<Vec<String>>,
    /// Worker threads (default: all cores)
    #[arg(long)]
    threads: Option<usize>,
    /// Fail instead of perturbing when the direction is orthogonal to an edge
    #[arg(long)]
    assume_generic: bool,
    /// Use the brute-force d-subset enumerator instead of reverse search
    #[arg(long)]
    oracle: bool,
    /// Seed for the default direction and the symbolic perturbation
    #[arg(long, default_value_t = 0x9e3779b97f4a7c15u64)]
    seed: u64,
    /// Node budget per subtree task (0 = one task per facet, no splitting)
    #[arg(long, default_value_t = 2000)]
    budget: u64,
    #[arg(long, short)]
    verbose: bool,
    /// Print run statistics as JSON to stderr
    #[arg(long)]
    stats_json: bool,
}

#[derive(Subcommand)]
enum Cmd {
    /// Geometry and cached sweep function for the browser viewer (JSON)
    ViewData {
        #[command(flatten)]
        common: Common,
    },
    /// Total volume of the union
    Volume {
        #[command(flatten)]
        common: Common,
    },
    /// Sweep function s(a, λ) as exact piecewise polynomial
    Sweep {
        #[command(flatten)]
        common: Common,
        /// Evaluate at λ (repeatable)
        #[arg(long, allow_hyphen_values = true)]
        at: Vec<String>,
        #[arg(long)]
        json: bool,
    },
    /// Sample s(a, λ) at N points between the extreme knots (CSV)
    Plot {
        #[command(flatten)]
        common: Common,
        #[arg(long, default_value_t = 100)]
        samples: usize,
    },
}

fn peak_rss_kb() -> Option<u64> {
    let s = std::fs::read_to_string("/proc/self/status").ok()?;
    for l in s.lines() {
        if let Some(rest) = l.strip_prefix("VmHWM:") {
            return rest.trim().trim_end_matches("kB").trim().parse().ok();
        }
    }
    None
}

fn run(common: &Common) -> Result<(SweepFunction, driver::Report, Vec<Q>), Box<dyn std::error::Error>> {
    let inp = read_inputs(&common.input)?;
    run_input(common, &inp)
}

fn run_input(common: &Common, inp: &nefvol::io::Input) -> Result<(SweepFunction, driver::Report, Vec<Q>), Box<dyn std::error::Error>> {
    if let Some(t) = common.threads {
        rayon::ThreadPoolBuilder::new().num_threads(t).build_global()?;
    }
    let direction: Vec<Q> = match &common.direction {
        Some(v) => {
            let dir: Vec<Q> = v.iter().map(|s| Q::parse(s).ok_or_else(|| format!("bad direction entry '{s}'"))).collect::<Result<_, _>>()?;
            if dir.len() != inp.d {
                return Err(format!("direction has {} entries, input has dimension {}", dir.len(), inp.d).into());
            }
            dir
        }
        None => driver::default_direction(inp.d, common.seed),
    };
    let opts = Options {
        assume_generic: common.assume_generic,
        seed: common.seed,
        oracle: common.oracle,
        verbose: common.verbose,
        budget: if common.budget == 0 { None } else { Some(common.budget) },
    };
    if common.verbose {
        eprintln!("d = {}, polytopes = {}, direction = {}", inp.d, inp.polytopes.len(), fmt_vec(&direction));
    }
    let t0 = std::time::Instant::now();
    let (f, rep) = driver::solve(inp.d, &inp.polytopes, &direction, &opts)?;
    let wall = t0.elapsed().as_millis();
    if common.verbose {
        eprintln!(
            "hyperplanes = {}, facet tasks = {}, redundant removed = {}, prep {} ms, search {} ms",
            rep.hyperplanes, rep.facets, rep.redundant_removed, rep.prep_ms, rep.search_ms
        );
        let s = &rep.search;
        eprintln!(
            "vertices visited = {} (emitted {}, owned elsewhere {}, not boundary {}), bases = {}, terms = {} (+{} Laurent), max depth = {}, knots = {}",
            s.nodes, s.emitted, s.owned_elsewhere, s.not_boundary, s.kernel.bases, s.kernel.terms, s.kernel.laurent_terms, s.max_depth, f.knots.len()
        );
    }
    if common.stats_json {
        let s = &rep.search;
        let hist: Vec<String> = s.subtree_hist.iter().enumerate().filter(|(_, c)| **c > 0).map(|(b, c)| format!("\"{}\":{}", 1u64 << b, c)).collect();
        eprintln!(
            "{{\"d\":{},\"polytopes\":{},\"hyperplanes\":{},\"facet_tasks\":{},\"redundant_removed\":{},\"vertices_visited\":{},\"vertices_emitted\":{},\"owned_elsewhere\":{},\"not_boundary\":{},\"edges_tested\":{},\"bases\":{},\"terms\":{},\"laurent_terms\":{},\"max_depth\":{},\"max_stack_dirs\":{},\"knots\":{},\"subtree_histogram\":{{{}}},\"threads\":{},\"prep_ms\":{},\"search_ms\":{},\"wall_ms\":{},\"peak_rss_kb\":{}}}",
            rep.d, rep.polytopes, rep.hyperplanes, rep.facets, rep.redundant_removed, s.nodes, s.emitted, s.owned_elsewhere, s.not_boundary, s.edges_tested,
            s.kernel.bases, s.kernel.terms, s.kernel.laurent_terms, s.max_depth, s.max_stack_dirs, f.knots.len(), hist.join(","),
            rayon::current_num_threads(), rep.prep_ms, rep.search_ms, wall, peak_rss_kb().unwrap_or(0)
        );
    }
    Ok((f, rep, direction))
}

fn fmt_vec(v: &[Q]) -> String {
    let s: Vec<String> = v.iter().map(|x| x.to_string()).collect();
    format!("({})", s.join(", "))
}

fn poly_str(p: &[Q]) -> String {
    let mut parts = Vec::new();
    for (i, c) in p.iter().enumerate() {
        if c.is_zero() {
            continue;
        }
        parts.push(match i {
            0 => format!("{c}"),
            1 => format!("{c}·λ"),
            _ => format!("{c}·λ^{i}"),
        });
    }
    if parts.is_empty() {
        "0".into()
    } else {
        parts.join(" + ")
    }
}

fn main() -> ExitCode {
    let cli = Cli::parse();
    let res = match &cli.cmd {
        Cmd::ViewData { common } => (|| {
            let input = read_inputs(&common.input)?;
            if input.d != 2 { return Err("the viewer currently supports only 2D inputs".into()); }
            let (f, _, direction) = run_input(common, &input)?;
            println!("{}", nefvol::viewer::scene(&input, &direction, &f)?);
            Ok(())
        })(),
        Cmd::Volume { common } => run(common).and_then(|(f, _, _)| {
            let v = f.total()?;
            println!("{v}");
            if f.knots.is_empty() {
                Err("no boundary vertices found".into())
            } else {
                Ok(())
            }
        }),
        Cmd::Sweep { common, at, json } => run(common).and_then(|(f, rep, dir)| {
            let ats: Vec<Q> = at.iter().map(|s| Q::parse(s).ok_or_else(|| format!("bad λ '{s}'"))).collect::<Result<_, _>>()?;
            let total = f.total()?;
            let pieces = f.pieces();
            if *json {
                let knots: Vec<String> = f.knots.iter().map(|k| format!("\"{}\"", k.0)).collect();
                let ps: Vec<String> = pieces
                    .iter()
                    .map(|p| {
                        let lo = p.lo.as_ref().map(|x| format!("\"{x}\"")).unwrap_or("null".into());
                        let hi = p.hi.as_ref().map(|x| format!("\"{x}\"")).unwrap_or("null".into());
                        let c: Vec<String> = p.poly.iter().map(|x| format!("\"{x}\"")).collect();
                        format!("{{\"lo\":{lo},\"hi\":{hi},\"coefficients\":[{}]}}", c.join(","))
                    })
                    .collect();
                let evals: Vec<String> = ats.iter().map(|l| format!("{{\"lambda\":\"{l}\",\"value\":\"{}\"}}", f.eval(l))).collect();
                let dv: Vec<String> = dir.iter().map(|x| format!("\"{x}\"")).collect();
                println!(
                    "{{\"d\":{},\"direction\":[{}],\"knots\":[{}],\"pieces\":[{}],\"total\":\"{}\",\"at\":[{}],\"vertices_emitted\":{}}}",
                    rep.d,
                    dv.join(","),
                    knots.join(","),
                    ps.join(","),
                    total,
                    evals.join(","),
                    rep.search.emitted
                );
            } else {
                println!("direction a = {}", fmt_vec(&dir));
                println!("breakpoints ({}):", f.knots.len());
                for (k, _) in &f.knots {
                    println!("  {k}");
                }
                println!("pieces (coefficients of λ^0..λ^{}):", rep.d);
                for p in &pieces {
                    let lo = p.lo.as_ref().map(|x| x.to_string()).unwrap_or("-inf".into());
                    let hi = p.hi.as_ref().map(|x| x.to_string()).unwrap_or("+inf".into());
                    println!("  [{lo}, {hi}): {}", poly_str(&p.poly));
                }
                println!("total volume: {total}");
                for l in &ats {
                    println!("s({l}) = {}", f.eval(l));
                }
            }
            Ok(())
        }),
        Cmd::Plot { common, samples } => run(common).and_then(|(f, _, _)| {
            let (lo, hi) = match (f.min_knot(), f.max_knot()) {
                (Some(a), Some(b)) => (a.clone(), b.clone()),
                _ => return Err("no knots".into()),
            };
            let n = (*samples).max(2);
            println!("lambda,volume");
            for i in 0..n {
                let t = Q::new(nefvol::num::Z::from(i as i64), nefvol::num::Z::from((n - 1) as i64));
                let l = lo.add(&hi.sub(&lo).mul(&t));
                println!("{},{}", l.to_f64(), f.eval(&l).to_f64());
            }
            Ok(())
        }),
    };
    match res {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("error: {e}");
            ExitCode::from(2)
        }
    }
}
