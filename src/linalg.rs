//! Exact integer linear algebra for the local kernel.
//!
//! * [`det_adj`]: one fraction-free Gauss–Jordan (Bareiss) pass over `[N | I]`
//!   returning `det N` and `adj N = det N · N^{-1}` — the single factorisation
//!   of SPEC §4.1 from which vertex, determinant, edge directions and the dual
//!   vector are all read off.
//! * [`Echelon`]: incremental fraction-free row echelon form used to prune the
//!   subset enumeration at a vertex (dependent prefixes are never extended) and
//!   to compute the direction of a line given `d−1` independent hyperplanes.

use crate::num::{canonical_direction, make_primitive, Z};

/// Returns `(det, adj)` with `adj` row-major `d×d` and `N · adj = det · I`.
/// `det` is the true determinant of `rows`; if singular returns `(0, [])`.
pub fn det_adj(rows: &[&[Z]]) -> (Z, Vec<Z>) {
    let d = rows.len();
    let w = 2 * d;
    let mut m: Vec<Z> = Vec::with_capacity(d * w);
    for (i, r) in rows.iter().enumerate() {
        debug_assert_eq!(r.len(), d);
        m.extend_from_slice(r);
        for j in 0..d {
            m.push(if i == j { Z::ONE } else { Z::ZERO });
        }
    }
    let mut sign = 1i8;
    let mut prev = Z::ONE;
    for k in 0..d {
        let p = match (k..d).find(|&i| !m[i * w + k].is_zero()) {
            Some(p) => p,
            None => return (Z::ZERO, Vec::new()),
        };
        if p != k {
            for j in 0..w {
                m.swap(p * w + j, k * w + j);
            }
            sign = -sign;
        }
        let pivot = m[k * w + k].clone();
        for i in 0..d {
            if i == k {
                continue;
            }
            let f = m[i * w + k].clone();
            if f.is_zero() {
                // row unchanged up to the common rescale by pivot/prev
                if !prev.is_one() || !pivot.is_one() {
                    for j in 0..w {
                        let x = &m[i * w + j];
                        if !x.is_zero() {
                            m[i * w + j] = (x * &pivot).div_exact(&prev);
                        }
                    }
                }
            } else {
                for j in 0..w {
                    let a = &m[i * w + j];
                    let b = &m[k * w + j];
                    let val = if b.is_zero() {
                        if a.is_zero() {
                            continue;
                        }
                        a * &pivot
                    } else if a.is_zero() {
                        -(&f * b)
                    } else {
                        &(a * &pivot) - &(&f * b)
                    };
                    m[i * w + j] = val.div_exact(&prev);
                }
            }
        }
        prev = pivot;
    }
    // Left block is now prev·I with prev = det(P N) = sign·det N; right block = prev · N^{-1}.
    let det = if sign > 0 { prev } else { -&prev };
    let mut adj = Vec::with_capacity(d * d);
    for i in 0..d {
        for j in 0..d {
            let x = &m[i * w + d + j];
            adj.push(if sign > 0 { x.clone() } else { -x });
        }
    }
    (det, adj)
}

/// Incremental fraction-free echelon form.  Rows are reduced only against
/// rows inserted *earlier*, so `pop` is O(1).
#[derive(Clone, Debug)]
pub struct Echelon {
    d: usize,
    rows: Vec<Vec<Z>>,
    pivots: Vec<usize>,
}

impl Echelon {
    pub fn new(d: usize) -> Echelon {
        Echelon { d, rows: Vec::with_capacity(d), pivots: Vec::with_capacity(d) }
    }
    pub fn len(&self) -> usize {
        self.rows.len()
    }
    pub fn is_empty(&self) -> bool {
        self.rows.is_empty()
    }
    pub fn clear(&mut self) {
        self.rows.clear();
        self.pivots.clear();
    }

    /// Try to add `row`; returns `false` (and leaves the form unchanged) if it
    /// is linearly dependent on the current rows.
    pub fn try_push(&mut self, row: &[Z]) -> bool {
        debug_assert_eq!(row.len(), self.d);
        let mut r: Vec<Z> = row.to_vec();
        for (er, &p) in self.rows.iter().zip(&self.pivots) {
            if r[p].is_zero() {
                continue;
            }
            let c = r[p].clone();
            let pv = er[p].clone();
            for j in 0..self.d {
                let a = &r[j];
                let b = &er[j];
                r[j] = if b.is_zero() {
                    if a.is_zero() {
                        continue;
                    }
                    a * &pv
                } else if a.is_zero() {
                    -(&c * b)
                } else {
                    &(a * &pv) - &(&c * b)
                };
            }
            debug_assert!(r[p].is_zero());
        }
        let p = match r.iter().position(|x| !x.is_zero()) {
            Some(p) => p,
            None => return false,
        };
        make_primitive(&mut r);
        self.rows.push(r);
        self.pivots.push(p);
        true
    }

    pub fn pop(&mut self) {
        self.rows.pop();
        self.pivots.pop();
    }

    /// Direction of the line orthogonal to all rows.  Requires `len() == d−1`.
    /// Result is primitive and lexicographically positive.
    pub fn nullspace_direction(&self) -> Vec<Z> {
        let d = self.d;
        assert_eq!(self.rows.len(), d - 1, "nullspace_direction needs rank d-1");
        let mut is_pivot = vec![false; d];
        for &p in &self.pivots {
            is_pivot[p] = true;
        }
        let free = (0..d).find(|&j| !is_pivot[j]).expect("one free column");
        let mut x: Vec<Z> = vec![Z::ZERO; d];
        let mut solved = vec![false; d];
        x[free] = Z::ONE;
        solved[free] = true;
        // rows in reverse insertion order: row i has zeros at pivots of rows < i,
        // so its unknowns are its own pivot plus already-solved columns.
        for i in (0..d - 1).rev() {
            let row = &self.rows[i];
            let p = self.pivots[i];
            let mut s = Z::ZERO;
            for j in 0..d {
                if j != p && solved[j] && !row[j].is_zero() && !x[j].is_zero() {
                    s += &(&row[j] * &x[j]);
                }
            }
            // x_p · row[p] = -s  ⇒ scale all solved by row[p], set x_p = -s
            let pv = &row[p];
            if !pv.is_one() {
                for j in 0..d {
                    if solved[j] && !x[j].is_zero() {
                        x[j] = &x[j] * pv;
                    }
                }
            }
            x[p] = -s;
            solved[p] = true;
        }
        canonical_direction(x)
    }
}

/// Determinant only (Bareiss).
pub fn det(rows: &[&[Z]]) -> Z {
    det_adj(rows).0
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::num::dot;

    fn z(v: &[i64]) -> Vec<Z> {
        v.iter().map(|&x| Z::from(x)).collect()
    }

    #[test]
    fn det_adj_2x2() {
        let r0 = z(&[2, 1]);
        let r1 = z(&[1, 3]);
        let (d, adj) = det_adj(&[&r0, &r1]);
        assert_eq!(d, Z::from(5));
        assert_eq!(adj, z(&[3, -1, -1, 2]));
    }

    #[test]
    fn det_adj_3x3_with_swap() {
        let rows = [z(&[0, 2, 1]), z(&[1, 0, 3]), z(&[4, 1, 0])];
        let rr: Vec<&[Z]> = rows.iter().map(|r| r.as_slice()).collect();
        let (d, adj) = det_adj(&rr);
        // det = 0*(0-3) - 2*(0-12) + 1*(1-0) = 24 + 1 = 25
        assert_eq!(d, Z::from(25));
        // N * adj = det I
        for (i, row) in rows.iter().enumerate() {
            for j in 0..3 {
                let col: Vec<Z> = (0..3).map(|k| adj[k * 3 + j].clone()).collect();
                let v = dot(row, &col);
                assert_eq!(v, if i == j { d.clone() } else { Z::ZERO }, "i={i} j={j}");
            }
        }
    }

    #[test]
    fn det_adj_random_identity() {
        // pseudo-random 5x5 with entries in [-9,9]
        let mut s: u64 = 12345;
        for _ in 0..20 {
            let mut rows = vec![];
            for _ in 0..5 {
                let mut r = vec![];
                for _ in 0..5 {
                    s ^= s << 13;
                    s ^= s >> 7;
                    s ^= s << 17;
                    r.push(Z::from((s % 19) as i64 - 9));
                }
                rows.push(r);
            }
            let rr: Vec<&[Z]> = rows.iter().map(|r| r.as_slice()).collect();
            let (d, adj) = det_adj(&rr);
            if d.is_zero() {
                continue;
            }
            for (i, row) in rows.iter().enumerate() {
                for j in 0..5 {
                    let col: Vec<Z> = (0..5).map(|k| adj[k * 5 + j].clone()).collect();
                    assert_eq!(dot(row, &col), if i == j { d.clone() } else { Z::ZERO });
                }
            }
        }
    }

    #[test]
    fn echelon_and_nullspace() {
        let mut e = Echelon::new(3);
        assert!(e.try_push(&z(&[1, 1, 0])));
        assert!(!e.try_push(&z(&[2, 2, 0])));
        assert!(e.try_push(&z(&[0, 1, 1])));
        let n = e.nullspace_direction();
        // orthogonal to (1,1,0) and (0,1,1): (1,-1,1)
        assert_eq!(n, z(&[1, -1, 1]));
        e.pop();
        assert!(e.try_push(&z(&[0, 0, 1])));
        assert_eq!(e.nullspace_direction(), z(&[1, -1, 0]));
    }
}
