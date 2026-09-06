//! The fused per-vertex kernel (SPEC §4).
//!
//! At an arrangement vertex `v` the only object ever computed is the central
//! arrangement of `H_v` (the hyperplanes through `v`), enumerated as its
//! bases `B ⊆ H_v` by a pruned subset DFS.  One fraction-free factorisation of
//! `N_B` gives:
//!
//! * the columns of `adj N_B` — the directions of the arrangement lines through
//!   `v` (the reverse search's edges, [`lines_at`]);
//! * the dual vector `Y = adjᵀ a`, so `⟨a, u_k⟩ = −σ_k Y_k / det`;
//! * the rows `R_j = n_j · adj` needed for the lexicographic cell test.
//!
//! **Cell test.**  Let every hyperplane `h` be shifted to `n_h·x = b_h + ε^{h+1}`
//! (ε → 0⁺).  The perturbed arrangement is simple, the perturbed union `U_ε`
//! has volume converging to `vol U`, and its sweep function is the sum over its
//! (simple) vertices `v_B` of Lawrence-type terms for the cells of the local
//! arrangement of `B` that lie in `U_ε`.  Cell `σ ∈ {±1}^d` of `B` (the cone
//! `{σ_k n_k·w ≤ 0}`) lies in `U_ε` iff some polytope `P_i ∋ v` has
//! `v_B ∈ P_{i,ε}` (a lexicographic sign test on `R_j` for every tight `j ∉ B`)
//! and `σ_k = o_{ik}` for all `k ∈ B ∩ J_i`.  Since the term for cell `σ`
//! is `sgn(σ)·f(B)`, the basis contributes `χ_B · f(B)` where
//! `χ_B = Σ_{σ passing} Π σ_k`, and `f(B) = det^d / (d!·|det|·Π_k(−Y_k))`
//! multiplying `(λ − ⟨a,v⟩)_+^d`.  Vertices interior to `U` sum to zero, so
//! only `∂U` vertices are ever visited.
//!
//! **Non-generic `a`.**  When some `Y_k = 0` with `χ_B ≠ 0`, the term is
//! carried as a truncated Laurent series in `ε` for `a_ε = a + ε r` (SPEC §3.4)
//! together with the expansion of `(λ − ⟨a_ε, v⟩)^d`; the accumulator keeps
//! the negative orders so their cancellation can be verified.

use crate::geom::{Problem, Tight};
use crate::linalg::{det_adj, Echelon};
use crate::num::{binomial, factorial, Point, Q, Z};

/// Statistics of one kernel invocation.
#[derive(Clone, Debug, Default)]
pub struct KernelStats {
    pub bases: u64,
    pub terms: u64,
    pub laurent_terms: u64,
}

/// Contribution of one vertex: coefficients `c[m][j]` of `ε^{-m} (λ − κ)_+^j`.
/// For generic `a` only `c[0][d]` is non-zero and `laurent` is `None`.
#[derive(Clone, Debug)]
pub struct VertexContribution {
    pub d: usize,
    pub top: Q,
    /// `(d+1)×(d+1)` row-major `[m][j]`, present only if a pole was met.
    pub laurent: Option<Vec<Q>>,
}

impl VertexContribution {
    pub fn zero(d: usize) -> Self {
        VertexContribution { d, top: Q::zero(), laurent: None }
    }
    pub fn is_zero(&self) -> bool {
        self.top.is_zero() && self.laurent.as_ref().is_none_or(|l| l.iter().all(|q| q.is_zero()))
    }
    fn laurent_mut(&mut self) -> &mut Vec<Q> {
        let n = (self.d + 1) * (self.d + 1);
        self.laurent.get_or_insert_with(|| vec![Q::zero(); n])
    }
}

/// Enumerate the lines of the central arrangement of `H_v` (given as sorted
/// hyperplane indices `hv`).  If `forced` is given, only lines lying inside
/// that hyperplane are returned (the edges of the facet subproblem, SPEC §2.2).
/// Directions are primitive, lex-positive and deduplicated.
pub fn lines_at(prob: &Problem, hv: &[usize], forced: Option<usize>) -> Vec<Vec<Z>> {
    let d = prob.d;
    let mut ech = Echelon::new(d);
    let mut out: Vec<Vec<Z>> = Vec::new();
    if d == 1 {
        // the only line is R^1 itself, unless we are confined to a point
        if forced.is_none() {
            out.push(vec![Z::ONE]);
        }
        return out;
    }
    let mut rest: Vec<usize> = Vec::with_capacity(hv.len());
    if let Some(f) = forced {
        debug_assert!(hv.contains(&f));
        let ok = ech.try_push(&prob.hyps[f].n);
        debug_assert!(ok);
        rest.extend(hv.iter().copied().filter(|&h| h != f));
    } else {
        rest.extend_from_slice(hv);
    }
    let need = d - 1;
    // iterative DFS over `rest` choosing indices in increasing order
    let mut choice: Vec<usize> = Vec::with_capacity(need);
    let mut next = 0usize;
    loop {
        if ech.len() == need {
            out.push(ech.nullspace_direction());
            // backtrack
            match choice.pop() {
                Some(c) => {
                    ech.pop();
                    next = c + 1;
                }
                None => break,
            }
            continue;
        }
        // remaining slots must be fillable
        let slots = need - ech.len();
        if next + slots > rest.len() {
            match choice.pop() {
                Some(c) => {
                    ech.pop();
                    next = c + 1;
                }
                None => break,
            }
            continue;
        }
        if ech.try_push(&prob.hyps[rest[next]].n) {
            choice.push(next);
            next += 1;
        } else {
            next += 1;
        }
    }
    out.sort();
    out.dedup();
    out
}

/// A polytope containing `v` on its boundary, with its tight constraints.
pub struct Containing<'a> {
    pub poly: usize,
    pub tight: &'a [Tight],
}

/// `χ` of a union of sub-cubes of `{±1}^d`.  A cube is `(mask, vals)`: the
/// coordinates in `mask` are fixed to `−1` where the bit of `vals` is set,
/// `+1` otherwise.  `χ(S) = Σ_{σ∈S} Π_k σ_k`; cubes with a free coordinate
/// contribute 0, so inclusion–exclusion only needs subsets whose masks cover
/// everything.
fn chi_union(cubes: &[(u64, u64)], d: usize) -> i64 {
    let full: u64 = if d >= 64 { u64::MAX } else { (1u64 << d) - 1 };
    // fast path: single cube
    if cubes.len() == 1 {
        let (m, v) = cubes[0];
        return if m == full { if v.count_ones().is_multiple_of(2) { 1 } else { -1 } } else { 0 };
    }
    fn rec(cubes: &[(u64, u64)], i: usize, mask: u64, vals: u64, size: usize, full: u64, acc: &mut i64) {
        if i == cubes.len() {
            if size > 0 && mask == full {
                let point = if vals.count_ones().is_multiple_of(2) { 1 } else { -1 };
                *acc += if size % 2 == 1 { point } else { -point };
            }
            return;
        }
        // skip cube i
        rec(cubes, i + 1, mask, vals, size, full, acc);
        // take cube i if consistent
        let (m, v) = cubes[i];
        let common = mask & m;
        if (vals & common) == (v & common) {
            rec(cubes, i + 1, mask | m, vals | (v & m), size + 1, full, acc);
        }
    }
    let mut acc = 0;
    rec(cubes, 0, 0, 0, 0, full, &mut acc);
    acc
}

/// Options for the sweep-plane half of the kernel.
pub struct SweepParams<'a> {
    /// integer sweep direction `a`
    pub a: &'a [Z],
    /// perturbation direction `r` for non-generic `a` (`None` = assume generic)
    pub r: Option<&'a [Z]>,
}

/// Errors the kernel can raise.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum KernelError {
    /// `a` is orthogonal to an edge direction and no perturbation was allowed.
    NonGeneric { vertex: String },
    /// Both `a` and `r` orthogonal to the same edge: choose another seed.
    BadPerturbation { vertex: String },
}

impl std::fmt::Display for KernelError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            KernelError::NonGeneric { vertex } => write!(f, "sweep direction is orthogonal to an arrangement edge at vertex {vertex}; drop --assume-generic"),
            KernelError::BadPerturbation { vertex } => write!(f, "perturbation direction is also orthogonal to an edge at vertex {vertex}; use another --seed"),
        }
    }
}

/// The sweep-plane half at `v`: sum over all bases `B ⊆ H_v` of the signed
/// Lawrence term (see module docs).  `containing` lists the polytopes with
/// `v` on their boundary; the caller guarantees `v ∉ int P_i` for all `i`.
pub fn cones_at(
    prob: &Problem,
    v: &Point,
    hv: &[usize],
    containing: &[Containing],
    params: &SweepParams,
    stats: &mut KernelStats,
) -> Result<VertexContribution, KernelError> {
    let d = prob.d;
    let mut contrib = VertexContribution::zero(d);
    if hv.len() < d {
        return Ok(contrib);
    }
    let dfact = factorial(d);
    let mut ech = Echelon::new(d);
    let mut choice: Vec<usize> = Vec::with_capacity(d);
    let mut next = 0usize;
    let mut cubes: Vec<(u64, u64)> = Vec::new();
    let mut y: Vec<Z> = vec![Z::ZERO; d];
    let mut yr: Vec<Z> = vec![Z::ZERO; d];
    let rho = params.r.map(|r| v.dot_q(r));

    loop {
        if choice.len() == d {
            stats.bases += 1;
            let basis: Vec<usize> = choice.iter().map(|&c| hv[c]).collect();
            let rows: Vec<&[Z]> = basis.iter().map(|&h| prob.hyps[h].n.as_slice()).collect();
            let (det, adj) = det_adj(&rows);
            debug_assert!(!det.is_zero());
            let dsign = det.signum();
            // --- cell test per containing polytope
            cubes.clear();
            'poly: for c in containing {
                let mut mask = 0u64;
                let mut vals = 0u64;
                for t in c.tight {
                    let h = t.c.h;
                    if let Some(k) = basis.iter().position(|&b| b == h) {
                        mask |= 1 << k;
                        if t.c.o < 0 {
                            vals |= 1 << k;
                        }
                    } else {
                        // lex sign of  Σ_{k∈B} (R_hk/det) ε_k − ε_h  with ε_k = ε^{k+1}
                        let mut sign: i8 = -1;
                        for (k, &bk) in basis.iter().enumerate() {
                            if bk > h {
                                break;
                            }
                            let nh = &prob.hyps[h].n;
                            let mut r = Z::ZERO;
                            for j in 0..d {
                                if !nh[j].is_zero() && !adj[j * d + k].is_zero() {
                                    r += &(&nh[j] * &adj[j * d + k]);
                                }
                            }
                            if !r.is_zero() {
                                sign = r.signum() * dsign;
                                break;
                            }
                        }
                        if t.c.o * sign >= 0 {
                            continue 'poly;
                        }
                    }
                }
                cubes.push((mask, vals));
            }
            if cubes.is_empty() {
                backtrack(&mut choice, &mut ech, &mut next);
                if choice.is_empty() && next >= hv.len() {
                    break;
                }
                continue;
            }
            let chi = chi_union(&cubes, d);
            if chi != 0 {
                // Y_k = a · adj[:,k]
                let mut zeros = 0usize;
                for k in 0..d {
                    let mut s = Z::ZERO;
                    for j in 0..d {
                        if !params.a[j].is_zero() && !adj[j * d + k].is_zero() {
                            s += &(&params.a[j] * &adj[j * d + k]);
                        }
                    }
                    if s.is_zero() {
                        zeros += 1;
                    }
                    y[k] = s;
                }
                // f = chi · det^d / (d! · |det| · Π(−Y_k))
                if zeros == 0 {
                    stats.terms += 1;
                    let mut den = &dfact * &det.abs();
                    for yk in y.iter().take(d) {
                        den = &den * &(-yk);
                    }
                    let num = &det.pow(d as u32) * &Z::from(chi);
                    contrib.top = contrib.top.add(&Q::new(num, den));
                } else {
                    let r = match params.r {
                        Some(r) => r,
                        None => return Err(KernelError::NonGeneric { vertex: v.to_string() }),
                    };
                    stats.laurent_terms += 1;
                    for k in 0..d {
                        let mut s = Z::ZERO;
                        for j in 0..d {
                            if !r[j].is_zero() && !adj[j * d + k].is_zero() {
                                s += &(&r[j] * &adj[j * d + k]);
                            }
                        }
                        yr[k] = s;
                    }
                    laurent_term(&mut contrib, d, chi, &det, &dfact, &y, &yr, zeros, rho.as_ref().unwrap(), v)?;
                }
            }
            backtrack(&mut choice, &mut ech, &mut next);
            if choice.is_empty() && next >= hv.len() {
                break;
            }
            continue;
        }
        let slots = d - choice.len();
        if next + slots > hv.len() {
            if choice.is_empty() {
                break;
            }
            backtrack(&mut choice, &mut ech, &mut next);
            continue;
        }
        if ech.try_push(&prob.hyps[hv[next]].n) {
            choice.push(next);
        }
        next += 1;
    }
    Ok(contrib)
}

#[inline]
fn backtrack(choice: &mut Vec<usize>, ech: &mut Echelon, next: &mut usize) {
    if let Some(c) = choice.pop() {
        ech.pop();
        *next = c + 1;
    }
}

/// Laurent expansion of one basis term with `p = zeros` vanishing `Y_k`.
///
/// `T(ε) = C · ε^{-p} · Π_{Y_k≠0} 1/(1 + ε Y'_k/Y_k)` with
/// `C = χ det^d / (d! |det| Π_{Y_k=0}(−Y'_k) Π_{Y_k≠0}(−Y_k))`, multiplied by
/// `(λ − κ − ε ρ)^d = Σ_j C(d,j) (λ−κ)^{d−j} (−ρ)^j ε^j`.
#[allow(clippy::too_many_arguments)]
fn laurent_term(
    contrib: &mut VertexContribution,
    d: usize,
    chi: i64,
    det: &Z,
    dfact: &Z,
    y: &[Z],
    yr: &[Z],
    p: usize,
    rho: &Q,
    v: &Point,
) -> Result<(), KernelError> {
    let mut den = dfact * &det.abs();
    for (yk, yrk) in y.iter().zip(yr) {
        if yk.is_zero() {
            if yrk.is_zero() {
                return Err(KernelError::BadPerturbation { vertex: v.to_string() });
            }
            den = &den * &(-yrk);
        } else {
            den = &den * &(-yk);
        }
    }
    let c = Q::new(&det.pow(d as u32) * &Z::from(chi), den);
    // power series  Π_{Y_k≠0} 1/(1 + ε t_k),  t_k = Y'_k/Y_k, truncated at order p
    let mut series: Vec<Q> = vec![Q::zero(); p + 1];
    series[0] = Q::one();
    for k in 0..d {
        if y[k].is_zero() {
            continue;
        }
        let t = Q::new(yr[k].clone(), y[k].clone());
        if t.is_zero() {
            continue;
        }
        // multiply by Σ_n (−t)^n ε^n
        let mut ns = vec![Q::zero(); p + 1];
        let mut pw = Q::one();
        for n in 0..=p {
            if n > 0 {
                pw = pw.mul(&t.neg());
            }
            for m in 0..=(p - n) {
                if !series[m].is_zero() {
                    ns[m + n] = ns[m + n].add(&series[m].mul(&pw));
                }
            }
        }
        series = ns;
    }
    // [ε^{-m}] T = C · series[p − m],  m = 0..p
    // coefficient of ε^{-m} (λ−κ)^{d−j}:  Σ_j [ε^{-m-j}]T · C(d,j) (−ρ)^j
    let l = contrib.laurent_mut();
    let w = d + 1;
    let mut rpow = Q::one();
    for j in 0..=d {
        if j > 0 {
            rpow = rpow.mul(&rho.neg());
        }
        let bin = Q::int(binomial(d, j));
        for m in 0..=d {
            let order = m + j; // need [ε^{-(m+j)}] T
            if order > p {
                continue;
            }
            let coef = c.mul(&series[p - order]).mul(&bin).mul(&rpow);
            if !coef.is_zero() {
                l[m * w + (d - j)] = l[m * w + (d - j)].add(&coef);
            }
        }
    }
    Ok(())
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn chi_basic() {
        // single full cube point (+,+): +1 ; (+,-): -1
        assert_eq!(chi_union(&[(0b11, 0b00)], 2), 1);
        assert_eq!(chi_union(&[(0b11, 0b10)], 2), -1);
        // free coordinate → 0
        assert_eq!(chi_union(&[(0b01, 0b00)], 2), 0);
        // two half-cubes covering the whole cube → 0 (R^2)
        assert_eq!(chi_union(&[(0b01, 0b00), (0b01, 0b01)], 2), 0);
        // L-shape: 3 of 4 quadrants: (+,+),(−,+),(−,−): χ = 1 − 1 + 1 = 1
        assert_eq!(chi_union(&[(0b11, 0b00), (0b11, 0b01), (0b11, 0b11)], 2), 1);
        // two overlapping cubes: {σ_0=+} ∪ {σ_1=+} = 3 points: (+,+),(+,-),(-,+): 1 −1 −1 = −1
        assert_eq!(chi_union(&[(0b01, 0b00), (0b10, 0b00)], 2), -1);
    }
}
