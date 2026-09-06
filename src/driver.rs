//! End-to-end pipeline: preprocessing → shared arrangement → facet tasks →
//! parallel reverse search → sweep function.

use crate::geom::{Problem, RawConstraint};
use crate::local::SweepParams;
use crate::num::{integerize, Q, Z};
use crate::oracle::oracle_sweep;
use crate::preprocess::{prepare, PrepError};
use crate::reverse::{facet_root, run_facet, run_subtree, FacetTask, SearchStats, SubtreeTask};
use crate::sweep::{Accumulator, SweepFunction};
use rayon::prelude::*;

#[derive(Clone, Debug)]
pub struct Options {
    pub assume_generic: bool,
    pub seed: u64,
    pub oracle: bool,
    pub verbose: bool,
    /// node budget per subtree task before the remainder goes back to the
    /// work queue (`None` = one task per facet)
    pub budget: Option<u64>,
}

impl Default for Options {
    fn default() -> Self {
        Options { assume_generic: false, seed: 0x9e3779b97f4a7c15, oracle: false, verbose: false, budget: Some(2000) }
    }
}

#[derive(Debug)]
pub enum Error {
    Input(String),
    Polytope { index: usize, err: PrepError },
    Kernel(crate::local::KernelError),
    Internal(String),
}

impl std::fmt::Display for Error {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Error::Input(s) => write!(f, "input error: {s}"),
            Error::Polytope { index, err } => write!(f, "polytope {index}: {err}"),
            Error::Kernel(e) => write!(f, "{e}"),
            Error::Internal(s) => write!(f, "{s}"),
        }
    }
}
impl std::error::Error for Error {}

#[derive(Clone, Debug, Default)]
pub struct Report {
    pub d: usize,
    pub polytopes: usize,
    pub hyperplanes: usize,
    pub facets: usize,
    pub redundant_removed: usize,
    pub search: SearchStats,
    pub prep_ms: u128,
    pub search_ms: u128,
}

pub struct Prepared {
    pub prob: Problem,
    /// per-polytope restricted arrangement (SPEC §2.2), indexed like `polys`
    pub views: Vec<Problem>,
    pub tasks: Vec<FacetTask>,
    pub redundant_removed: usize,
}

/// Preprocess all polytopes and build the shared arrangement.
pub fn prepare_all(d: usize, raw: &[Vec<RawConstraint>]) -> Result<Prepared, Error> {
    if raw.is_empty() {
        return Err(Error::Input("no polytopes".into()));
    }
    let preps: Vec<_> = raw
        .par_iter()
        .enumerate()
        .map(|(i, r)| prepare(d, r).map_err(|err| Error::Polytope { index: i, err }))
        .collect::<Result<Vec<_>, _>>()?;
    let cleaned: Vec<Vec<RawConstraint>> = preps.iter().map(|p| p.cons.clone()).collect();
    let prob = Problem::from_raw(d, &cleaned).map_err(Error::Input)?;
    let mut tasks = Vec::new();
    for (i, p) in preps.iter().enumerate() {
        debug_assert_eq!(prob.polys[i].cons.len(), p.facet_starts.len());
        for (j, s) in p.facet_starts.iter().enumerate() {
            tasks.push(FacetTask { poly: i, facet: j, start: s.clone() });
        }
    }
    let bboxes: Vec<_> = preps.iter().map(|p| p.bbox.clone()).collect();
    let views: Vec<Problem> = (0..preps.len()).into_par_iter().map(|i| prob.local_view(i, &bboxes)).collect();
    Ok(Prepared { prob, views, tasks, redundant_removed: preps.iter().map(|p| p.removed_redundant).sum() })
}

/// Pseudo-random integer perturbation direction from a seed.
pub fn perturbation(d: usize, seed: u64) -> Vec<Z> {
    let mut s = seed | 1;
    (0..d)
        .map(|_| {
            s ^= s << 13;
            s ^= s >> 7;
            s ^= s << 17;
            let x = (s % (1u64 << 40)) as i64 - (1i64 << 39);
            Z::from(if x == 0 { 1 } else { x })
        })
        .collect()
}

/// Default sweep direction: generic-looking pseudo-random small integers.
pub fn default_direction(d: usize, seed: u64) -> Vec<Q> {
    let mut s = seed.wrapping_mul(0x2545F4914F6CDD1D) | 1;
    (0..d)
        .map(|_| {
            s ^= s << 13;
            s ^= s >> 7;
            s ^= s << 17;
            Q::int((s % 9973) as i64 + 1)
        })
        .collect()
}

/// Run the whole computation.
pub fn solve(d: usize, raw: &[Vec<RawConstraint>], direction: &[Q], opts: &Options) -> Result<(SweepFunction, Report), Error> {
    let t0 = std::time::Instant::now();
    let prep = prepare_all(d, raw)?;
    let prep_ms = t0.elapsed().as_millis();
    let a = integerize(direction);
    if a.iter().all(|x| x.is_zero()) {
        return Err(Error::Input("sweep direction must be non-zero".into()));
    }
    let r = perturbation(d, opts.seed);
    let params = SweepParams { a: &a, r: if opts.assume_generic { None } else { Some(&r) } };
    let t1 = std::time::Instant::now();
    let (acc, search) = if opts.oracle {
        let (acc, k, emitted) = oracle_sweep(&prep.prob, &params).map_err(Error::Kernel)?;
        let mut s = SearchStats::new();
        s.kernel = k;
        s.emitted = emitted as u64;
        (acc, s)
    } else {
        let results: Vec<(Accumulator, SearchStats)> = prep
            .tasks
            .par_iter()
            .map(|t| {
                // the facet search runs in the restricted view of its polytope,
                // where that polytope has local index 0
                let pr = &prep.views[t.poly];
                let lt = FacetTask { poly: 0, facet: t.facet, start: t.start.clone() };
                match opts.budget {
                    None => run_facet(pr, &lt, &params, &mut ()),
                    Some(b) => run_task_par(pr, facet_root(pr, &lt), &params, b),
                }
            })
            .collect::<Result<Vec<_>, _>>()
            .map_err(Error::Kernel)?;
        let mut acc = Accumulator::new(d);
        let mut stats = SearchStats::new();
        for (a, s) in results {
            acc.merge(a);
            stats.merge(&s);
        }
        (acc, stats)
    };
    let search_ms = t1.elapsed().as_millis();
    let mut f = acc.finish().map_err(Error::Internal)?;
    // the kernel used the integerised direction a_int = scale · a
    let i = a.iter().position(|x| !x.is_zero()).unwrap();
    let scale = Q::int(a[i].clone()).div(&direction[i]);
    f.rescale(&scale);
    let report = Report {
        d,
        polytopes: prep.prob.polys.len(),
        hyperplanes: prep.prob.hyps.len(),
        facets: prep.tasks.len(),
        redundant_removed: prep.redundant_removed,
        search,
        prep_ms,
        search_ms,
    };
    Ok((f, report))
}

/// Budgeted subtree exploration; leftovers are pushed to rayon's queue and
/// stolen by idle workers.  Every worker owns its accumulator; maps are only
/// merged on return (SPEC §5).
fn run_task_par(prob: &Problem, task: SubtreeTask, params: &SweepParams, budget: u64) -> Result<(Accumulator, SearchStats), crate::local::KernelError> {
    let (mut acc, mut stats, leftovers) = run_subtree(prob, &task, params, Some(budget), &mut ())?;
    if !leftovers.is_empty() {
        let subs: Vec<(Accumulator, SearchStats)> = leftovers
            .into_par_iter()
            .map(|t| run_task_par(prob, t, params, budget))
            .collect::<Result<Vec<_>, _>>()?;
        for (a, s) in subs {
            acc.merge(a);
            stats.merge(&s);
        }
    }
    Ok((acc, stats))
}

/// Total volume with a default direction.
pub fn volume(d: usize, raw: &[Vec<RawConstraint>], opts: &Options) -> Result<Q, Error> {
    let dir = default_direction(d, opts.seed);
    let (f, _) = solve(d, raw, &dir, opts)?;
    f.total().map_err(Error::Internal)
}
