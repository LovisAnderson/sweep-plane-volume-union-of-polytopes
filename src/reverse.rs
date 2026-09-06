//! The reverse-search half (SPEC §2): one depth-first traversal per facet
//! `F = h_j ∩ P_i`, enumerating every arrangement vertex in `F` exactly once
//! with no visited set.  Memory is the explicit stack of the current root
//! path only.

use crate::geom::Problem;
use crate::local::{KernelError, KernelStats, SweepParams};
use crate::num::{lex_sign, negate, Point, Z};
use crate::sweep::Accumulator;
use crate::vertex::{process_vertex, Visit};
use crate::walk::{best_improving, feasible_dirs, ratio_step, walk, Region};

#[derive(Clone, Debug)]
pub struct FacetTask {
    pub poly: usize,
    pub facet: usize,
    pub start: Point,
}

pub const HIST_BUCKETS: usize = 40;

#[derive(Clone, Debug, Default)]
pub struct SearchStats {
    pub tasks: u64,
    pub nodes: u64,
    pub emitted: u64,
    pub owned_elsewhere: u64,
    pub not_boundary: u64,
    pub edges_tested: u64,
    pub max_depth: usize,
    pub max_stack_dirs: usize,
    pub root_ascent_steps: u64,
    pub kernel: KernelStats,
    /// bucket `b` counts subtrees with `2^b ≤ size < 2^{b+1}`
    pub subtree_hist: Vec<u64>,
}

impl SearchStats {
    pub fn new() -> Self {
        SearchStats { subtree_hist: vec![0; HIST_BUCKETS], ..Default::default() }
    }
    pub fn merge(&mut self, o: &SearchStats) {
        self.tasks += o.tasks;
        self.nodes += o.nodes;
        self.emitted += o.emitted;
        self.owned_elsewhere += o.owned_elsewhere;
        self.not_boundary += o.not_boundary;
        self.edges_tested += o.edges_tested;
        self.max_depth = self.max_depth.max(o.max_depth);
        self.max_stack_dirs = self.max_stack_dirs.max(o.max_stack_dirs);
        self.root_ascent_steps += o.root_ascent_steps;
        self.kernel.bases += o.kernel.bases;
        self.kernel.terms += o.kernel.terms;
        self.kernel.laurent_terms += o.kernel.laurent_terms;
        if self.subtree_hist.len() < o.subtree_hist.len() {
            self.subtree_hist.resize(o.subtree_hist.len(), 0);
        }
        for (a, b) in self.subtree_hist.iter_mut().zip(&o.subtree_hist) {
            *a += b;
        }
    }
    fn record_subtree(&mut self, size: u64) {
        let b = (64 - size.leading_zeros()).saturating_sub(1) as usize;
        let b = b.min(HIST_BUCKETS - 1);
        self.subtree_hist[b] += 1;
    }
}

struct Frame {
    v: Point,
    evals: Vec<Z>,
    dirs: Vec<Vec<Z>>,
    idx: usize,
    subtree: u64,
}

/// Optional sink for every vertex visited (differential tests).
pub trait VertexSink {
    fn visited(&mut self, v: &Point, visit: Visit);
}
impl VertexSink for () {
    fn visited(&mut self, _v: &Point, _visit: Visit) {}
}
impl VertexSink for Vec<(Point, Visit)> {
    fn visited(&mut self, v: &Point, visit: Visit) {
        self.push((v.clone(), visit));
    }
}

/// A (sub)tree of one facet search: explore the children of `root` reached
/// through edges `start_edge..`, and their subtrees.  `visit_root` is set
/// only for the tree root itself.  Restarting from `(root, start_edge)` is
/// exact because the edge list at a vertex is a deterministic function of
/// the vertex alone (SPEC §2.3).
#[derive(Clone, Debug)]
pub struct SubtreeTask {
    pub poly: usize,
    pub facet: usize,
    pub root: Point,
    pub start_edge: usize,
    pub visit_root: bool,
}

/// Run the reverse search of one facet to completion.
pub fn run_facet<S: VertexSink>(
    prob: &Problem,
    task: &FacetTask,
    params: &SweepParams,
    sink: &mut S,
) -> Result<(Accumulator, SearchStats), KernelError> {
    let root_task = facet_root(prob, task);
    let (acc, stats, left) = run_subtree(prob, &root_task, params, None, sink)?;
    debug_assert!(left.is_empty());
    Ok((acc, stats))
}

/// Ascend from the facet's start vertex to the tree root (SPEC §2.3).
pub fn facet_root(prob: &Problem, task: &FacetTask) -> SubtreeTask {
    let facet_h = prob.polys[task.poly].cons[task.facet].h;
    let region = Region { cons: prob.polys[task.poly].cons.clone(), eq: Some(facet_h) };
    debug_assert!(region.contains(prob, &task.start), "facet start not on facet");
    let root = walk(prob, &region, &task.start, None).expect("facet region is bounded");
    SubtreeTask { poly: task.poly, facet: task.facet, root, start_edge: 0, visit_root: true }
}

/// Explore a subtree with an optional node budget (mplrs-style, SPEC §5).
/// When the budget is exhausted the unexplored remainder is returned as new
/// tasks: one per stack frame, `(vertex, next edge index)`.
pub fn run_subtree<S: VertexSink>(
    prob: &Problem,
    task: &SubtreeTask,
    params: &SweepParams,
    budget: Option<u64>,
    sink: &mut S,
) -> Result<(Accumulator, SearchStats, Vec<SubtreeTask>), KernelError> {
    let mut acc = Accumulator::new(prob.d);
    let mut stats = SearchStats::new();
    stats.tasks = 1;
    let facet_h = prob.polys[task.poly].cons[task.facet].h;
    let region = Region { cons: prob.polys[task.poly].cons.clone(), eq: Some(facet_h) };
    let owner = Some((prob.polys[task.poly].id, prob.polys[task.poly].pos[task.facet]));

    let visit = |v: &Point, evals: &[Z], acc: &mut Accumulator, stats: &mut SearchStats, sink: &mut S| -> Result<(), KernelError> {
        stats.nodes += 1;
        let r = process_vertex(prob, v, evals, owner, params, acc, &mut stats.kernel)?;
        match r {
            Visit::Emitted => stats.emitted += 1,
            Visit::OwnedElsewhere => stats.owned_elsewhere += 1,
            Visit::NotBoundary => stats.not_boundary += 1,
        }
        sink.visited(v, r);
        Ok(())
    };

    let mut stack: Vec<Frame> = Vec::new();
    let evals = prob.evals(&task.root);
    let hv = Problem::hv_from_evals(&evals);
    let tight = region.tight_evals(&evals);
    let dirs = feasible_dirs(prob, &region, &hv, &tight);
    if task.visit_root {
        debug_assert!(best_improving(&dirs, None).is_none(), "root has an improving edge");
        visit(&task.root, &evals, &mut acc, &mut stats, sink)?;
    }
    stack.push(Frame { v: task.root.clone(), evals, dirs, idx: task.start_edge, subtree: 1 });

    let mut leftovers = Vec::new();
    while let Some(fr) = stack.last_mut() {
        if fr.idx >= fr.dirs.len() {
            let done = stack.pop().unwrap();
            stats.record_subtree(done.subtree);
            if let Some(p) = stack.last_mut() {
                p.subtree += done.subtree;
            }
            continue;
        }
        if let Some(b) = budget {
            if stats.nodes >= b {
                // hand the rest of every frame back to the queue
                for f in stack.drain(..) {
                    if f.idx < f.dirs.len() {
                        leftovers.push(SubtreeTask { poly: task.poly, facet: task.facet, root: f.v, start_edge: f.idx, visit_root: false });
                    }
                }
                break;
            }
        }
        let r = fr.dirs[fr.idx].clone();
        fr.idx += 1;
        // the parent edge of a child points lex-upward, so from v it is lex-negative
        if lex_sign(&r) > 0 {
            continue;
        }
        stats.edges_tested += 1;
        let u = ratio_step(prob, &fr.evals, &fr.v, &r).expect("facet region is bounded");
        let evals_u = prob.evals(&u);
        let hu = Problem::hv_from_evals(&evals_u);
        let tight_u = region.tight_evals(&evals_u);
        let dirs_u = feasible_dirs(prob, &region, &hu, &tight_u);
        let back = negate(&r);
        let is_child = matches!(best_improving(&dirs_u, None), Some(b) if *b == back);
        if is_child {
            visit(&u, &evals_u, &mut acc, &mut stats, sink)?;
            stack.push(Frame { v: u, evals: evals_u, dirs: dirs_u, idx: 0, subtree: 1 });
            stats.max_depth = stats.max_depth.max(stack.len());
            let total_dirs: usize = stack.iter().map(|f| f.dirs.len()).sum();
            stats.max_stack_dirs = stats.max_stack_dirs.max(total_dirs);
        }
    }
    Ok((acc, stats, leftovers))
}
