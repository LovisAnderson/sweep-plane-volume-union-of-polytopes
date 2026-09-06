//! Exact arithmetic layer.
//!
//! `Z` is a hybrid integer: an `i128` fast path that promotes to `BigInt` on
//! overflow (the lrs-style escalation of SPEC §5).  The invariant is that the
//! `B` variant is only ever used for values that do *not* fit an `i128`, so
//! derived `Eq`/`Hash` are canonical.
//!
//! `Q` is a normalised rational over `Z`.

use num_bigint::BigInt;
use num_integer::Integer;
use num_traits::{Signed, ToPrimitive, Zero};
use std::cmp::Ordering;
use std::fmt;
use std::ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign};

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub enum Z {
    S(i128),
    B(BigInt),
}

impl Default for Z {
    fn default() -> Self {
        Z::S(0)
    }
}

impl Z {
    pub const ZERO: Z = Z::S(0);
    pub const ONE: Z = Z::S(1);

    #[inline]
    pub fn from_big(b: BigInt) -> Z {
        match b.to_i128() {
            Some(x) => Z::S(x),
            None => Z::B(b),
        }
    }

    #[inline]
    pub fn to_big(&self) -> BigInt {
        match self {
            Z::S(x) => BigInt::from(*x),
            Z::B(b) => b.clone(),
        }
    }

    #[inline]
    pub fn is_zero(&self) -> bool {
        matches!(self, Z::S(0))
    }

    #[inline]
    pub fn is_one(&self) -> bool {
        matches!(self, Z::S(1))
    }

    #[inline]
    pub fn signum(&self) -> i8 {
        match self {
            Z::S(x) => x.signum() as i8,
            Z::B(b) => {
                if b.is_positive() {
                    1
                } else {
                    -1
                }
            }
        }
    }

    #[inline]
    pub fn is_positive(&self) -> bool {
        self.signum() > 0
    }

    #[inline]
    pub fn is_negative(&self) -> bool {
        self.signum() < 0
    }

    pub fn abs(&self) -> Z {
        match self {
            Z::S(x) => match x.checked_abs() {
                Some(a) => Z::S(a),
                None => Z::from_big(BigInt::from(*x).abs()),
            },
            Z::B(b) => Z::from_big(b.abs()),
        }
    }

    pub fn gcd(&self, o: &Z) -> Z {
        match (self, o) {
            (Z::S(a), Z::S(b)) => {
                if *a == i128::MIN || *b == i128::MIN {
                    Z::from_big(BigInt::from(*a).gcd(&BigInt::from(*b)))
                } else {
                    Z::S(a.abs().gcd(&b.abs()))
                }
            }
            _ => Z::from_big(self.to_big().gcd(&o.to_big())),
        }
    }

    /// Exact division; panics in debug builds if not exact.
    pub fn div_exact(&self, o: &Z) -> Z {
        debug_assert!(!o.is_zero(), "division by zero");
        match (self, o) {
            (Z::S(a), Z::S(b)) => {
                debug_assert!(a % b == 0, "inexact division {a} / {b}");
                match a.checked_div(*b) {
                    Some(q) => Z::S(q),
                    None => Z::from_big(BigInt::from(*a) / BigInt::from(*b)),
                }
            }
            _ => {
                let (q, r) = self.to_big().div_rem(&o.to_big());
                debug_assert!(r.is_zero(), "inexact division");
                Z::from_big(q)
            }
        }
    }

    pub fn divides(&self, o: &Z) -> bool {
        match (self, o) {
            (Z::S(a), Z::S(b)) => *a != 0 && b % a == 0,
            _ => !self.is_zero() && (o.to_big() % self.to_big()).is_zero(),
        }
    }

    pub fn pow(&self, e: u32) -> Z {
        let mut r = Z::ONE;
        for _ in 0..e {
            r = &r * self;
        }
        r
    }

    pub fn to_i64(&self) -> Option<i64> {
        match self {
            Z::S(x) => i64::try_from(*x).ok(),
            Z::B(_) => None,
        }
    }

    pub fn to_f64(&self) -> f64 {
        match self {
            Z::S(x) => *x as f64,
            Z::B(b) => b.to_f64().unwrap_or(f64::NAN),
        }
    }

    pub fn parse(s: &str) -> Option<Z> {
        let s = s.trim();
        if let Ok(x) = s.parse::<i128>() {
            return Some(Z::S(x));
        }
        s.parse::<BigInt>().ok().map(Z::from_big)
    }
}

impl From<i64> for Z {
    fn from(x: i64) -> Z {
        Z::S(x as i128)
    }
}
impl From<i32> for Z {
    fn from(x: i32) -> Z {
        Z::S(x as i128)
    }
}
impl From<usize> for Z {
    fn from(x: usize) -> Z {
        Z::S(x as i128)
    }
}
impl From<i128> for Z {
    fn from(x: i128) -> Z {
        Z::S(x)
    }
}

impl PartialOrd for Z {
    fn partial_cmp(&self, o: &Z) -> Option<Ordering> {
        Some(self.cmp(o))
    }
}
impl Ord for Z {
    fn cmp(&self, o: &Z) -> Ordering {
        match (self, o) {
            (Z::S(a), Z::S(b)) => a.cmp(b),
            _ => self.to_big().cmp(&o.to_big()),
        }
    }
}

impl fmt::Display for Z {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Z::S(x) => write!(f, "{x}"),
            Z::B(b) => write!(f, "{b}"),
        }
    }
}

#[inline]
fn add_zz(a: &Z, b: &Z) -> Z {
    match (a, b) {
        (Z::S(x), Z::S(y)) => match x.checked_add(*y) {
            Some(s) => Z::S(s),
            None => Z::from_big(BigInt::from(*x) + BigInt::from(*y)),
        },
        _ => Z::from_big(a.to_big() + b.to_big()),
    }
}
#[inline]
fn sub_zz(a: &Z, b: &Z) -> Z {
    match (a, b) {
        (Z::S(x), Z::S(y)) => match x.checked_sub(*y) {
            Some(s) => Z::S(s),
            None => Z::from_big(BigInt::from(*x) - BigInt::from(*y)),
        },
        _ => Z::from_big(a.to_big() - b.to_big()),
    }
}
#[inline]
fn mul_zz(a: &Z, b: &Z) -> Z {
    match (a, b) {
        (Z::S(x), Z::S(y)) => match x.checked_mul(*y) {
            Some(s) => Z::S(s),
            None => Z::from_big(BigInt::from(*x) * BigInt::from(*y)),
        },
        _ => Z::from_big(a.to_big() * b.to_big()),
    }
}
#[inline]
fn neg_z(a: &Z) -> Z {
    match a {
        Z::S(x) => match x.checked_neg() {
            Some(s) => Z::S(s),
            None => Z::from_big(-BigInt::from(*x)),
        },
        Z::B(b) => Z::from_big(-b.clone()),
    }
}

macro_rules! impl_binop {
    ($tr:ident, $m:ident, $f:ident) => {
        impl $tr<&Z> for &Z {
            type Output = Z;
            #[inline]
            fn $m(self, o: &Z) -> Z {
                $f(self, o)
            }
        }
        impl $tr<Z> for &Z {
            type Output = Z;
            #[inline]
            fn $m(self, o: Z) -> Z {
                $f(self, &o)
            }
        }
        impl $tr<&Z> for Z {
            type Output = Z;
            #[inline]
            fn $m(self, o: &Z) -> Z {
                $f(&self, o)
            }
        }
        impl $tr<Z> for Z {
            type Output = Z;
            #[inline]
            fn $m(self, o: Z) -> Z {
                $f(&self, &o)
            }
        }
    };
}
impl_binop!(Add, add, add_zz);
impl_binop!(Sub, sub, sub_zz);
impl_binop!(Mul, mul, mul_zz);

impl Neg for &Z {
    type Output = Z;
    fn neg(self) -> Z {
        neg_z(self)
    }
}
impl Neg for Z {
    type Output = Z;
    fn neg(self) -> Z {
        neg_z(&self)
    }
}
impl AddAssign<&Z> for Z {
    fn add_assign(&mut self, o: &Z) {
        *self = add_zz(self, o);
    }
}
impl SubAssign<&Z> for Z {
    fn sub_assign(&mut self, o: &Z) {
        *self = sub_zz(self, o);
    }
}
impl MulAssign<&Z> for Z {
    fn mul_assign(&mut self, o: &Z) {
        *self = mul_zz(self, o);
    }
}

/// Dot product of two integer vectors (i128 fast path, exact fallback).
#[inline]
pub fn dot(a: &[Z], b: &[Z]) -> Z {
    debug_assert_eq!(a.len(), b.len());
    let mut acc: i128 = 0;
    for (x, y) in a.iter().zip(b) {
        match (x, y) {
            (Z::S(p), Z::S(q)) => match p.checked_mul(*q).and_then(|m| acc.checked_add(m)) {
                Some(v) => acc = v,
                None => return dot_slow(a, b),
            },
            _ => return dot_slow(a, b),
        }
    }
    Z::S(acc)
}

fn dot_slow(a: &[Z], b: &[Z]) -> Z {
    let mut s = Z::ZERO;
    for (x, y) in a.iter().zip(b) {
        if !x.is_zero() && !y.is_zero() {
            s += &(x * y);
        }
    }
    s
}

/// Gcd of all entries (non-negative; 0 if all zero).
pub fn gcd_all(v: &[Z]) -> Z {
    let mut g = Z::ZERO;
    for x in v {
        if !x.is_zero() {
            g = if g.is_zero() { x.abs() } else { g.gcd(x) };
            if g.is_one() {
                break;
            }
        }
    }
    g
}

/// Divide the vector by the gcd of its entries.  Does not touch the sign.
pub fn make_primitive(v: &mut [Z]) {
    let g = gcd_all(v);
    if !g.is_zero() && !g.is_one() {
        for x in v.iter_mut() {
            *x = x.div_exact(&g);
        }
    }
}

/// Sign of the first non-zero entry (0 for the zero vector).
pub fn lex_sign(v: &[Z]) -> i8 {
    for x in v {
        let s = x.signum();
        if s != 0 {
            return s;
        }
    }
    0
}

/// Primitive and lexicographically positive (first non-zero entry > 0).
pub fn canonical_direction(mut v: Vec<Z>) -> Vec<Z> {
    make_primitive(&mut v);
    if lex_sign(&v) < 0 {
        for x in v.iter_mut() {
            *x = -&*x;
        }
    }
    v
}

pub fn negate(v: &[Z]) -> Vec<Z> {
    v.iter().map(|x| -x).collect()
}

/// Normalised rational: `den > 0`, `gcd(num, den) = 1`.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct Q {
    pub num: Z,
    pub den: Z,
}

impl Q {
    pub fn new(num: Z, den: Z) -> Q {
        assert!(!den.is_zero(), "rational with zero denominator");
        let (mut num, mut den) = if den.is_negative() { (-num, -den) } else { (num, den) };
        let g = num.gcd(&den);
        if !g.is_one() && !g.is_zero() {
            num = num.div_exact(&g);
            den = den.div_exact(&g);
        }
        if num.is_zero() {
            den = Z::ONE;
        }
        Q { num, den }
    }
    pub fn int(x: impl Into<Z>) -> Q {
        Q { num: x.into(), den: Z::ONE }
    }
    pub fn zero() -> Q {
        Q::int(0)
    }
    pub fn one() -> Q {
        Q::int(1)
    }
    pub fn is_zero(&self) -> bool {
        self.num.is_zero()
    }
    pub fn is_integer(&self) -> bool {
        self.den.is_one()
    }
    pub fn signum(&self) -> i8 {
        self.num.signum()
    }
    pub fn neg(&self) -> Q {
        Q { num: -&self.num, den: self.den.clone() }
    }
    pub fn add(&self, o: &Q) -> Q {
        if self.den == o.den {
            return Q::new(&self.num + &o.num, self.den.clone());
        }
        Q::new(&(&self.num * &o.den) + &(&o.num * &self.den), &self.den * &o.den)
    }
    pub fn sub(&self, o: &Q) -> Q {
        self.add(&o.neg())
    }
    pub fn mul(&self, o: &Q) -> Q {
        Q::new(&self.num * &o.num, &self.den * &o.den)
    }
    pub fn div(&self, o: &Q) -> Q {
        assert!(!o.is_zero(), "division by zero rational");
        Q::new(&self.num * &o.den, &self.den * &o.num)
    }
    pub fn recip(&self) -> Q {
        Q::one().div(self)
    }
    pub fn mul_int(&self, z: &Z) -> Q {
        Q::new(&self.num * z, self.den.clone())
    }
    pub fn div_int(&self, z: &Z) -> Q {
        Q::new(self.num.clone(), &self.den * z)
    }
    pub fn pow(&self, e: u32) -> Q {
        Q { num: self.num.pow(e), den: self.den.pow(e) }
    }
    pub fn to_f64(&self) -> f64 {
        match (&self.num, &self.den) {
            (Z::S(n), Z::S(d)) => *n as f64 / *d as f64,
            _ => {
                let n = self.num.to_big();
                let d = self.den.to_big();
                // scale to avoid overflow of to_f64 on huge values
                let shift = d.bits().saturating_sub(60);
                let d2 = &d >> shift;
                let n2 = &n >> shift;
                n2.to_f64().unwrap_or(f64::NAN) / d2.to_f64().unwrap_or(f64::NAN)
            }
        }
    }
    /// Parse "a", "a/b", or a decimal like "-1.25".
    pub fn parse(s: &str) -> Option<Q> {
        let s = s.trim();
        if let Some((a, b)) = s.split_once('/') {
            let n = Z::parse(a)?;
            let d = Z::parse(b)?;
            if d.is_zero() {
                return None;
            }
            return Some(Q::new(n, d));
        }
        if let Some((ip, fp)) = s.split_once('.') {
            let neg = ip.starts_with('-');
            let ip = ip.trim_start_matches(['-', '+']);
            let ip = if ip.is_empty() { "0" } else { ip };
            let digits = format!("{ip}{fp}");
            let n = Z::parse(&digits)?;
            let d = Z::from(10i64).pow(fp.len() as u32);
            let q = Q::new(n, d);
            return Some(if neg { q.neg() } else { q });
        }
        if let Some((m, e)) = s.split_once(['e', 'E']) {
            let base = Q::parse(m)?;
            let e: i32 = e.parse().ok()?;
            let p = Q::int(10).pow(e.unsigned_abs());
            return Some(if e >= 0 { base.mul(&p) } else { base.div(&p) });
        }
        Z::parse(s).map(Q::int)
    }
}

impl PartialOrd for Q {
    fn partial_cmp(&self, o: &Q) -> Option<Ordering> {
        Some(self.cmp(o))
    }
}
impl Ord for Q {
    fn cmp(&self, o: &Q) -> Ordering {
        // both dens positive
        (&self.num * &o.den).cmp(&(&o.num * &self.den))
    }
}

impl fmt::Display for Q {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.den.is_one() {
            write!(f, "{}", self.num)
        } else {
            write!(f, "{}/{}", self.num, self.den)
        }
    }
}

/// A rational point stored as integer numerators over one positive denominator,
/// normalised so that `gcd(num_1, …, num_d, den) = 1`.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct Point {
    pub num: Vec<Z>,
    pub den: Z,
}

impl Point {
    pub fn new(mut num: Vec<Z>, den: Z) -> Point {
        assert!(!den.is_zero());
        let mut den = den;
        if den.is_negative() {
            den = -den;
            for x in num.iter_mut() {
                *x = -&*x;
            }
        }
        let mut g = den.clone();
        for x in &num {
            if g.is_one() {
                break;
            }
            if !x.is_zero() {
                g = g.gcd(x);
            }
        }
        if !g.is_one() {
            for x in num.iter_mut() {
                *x = x.div_exact(&g);
            }
            den = den.div_exact(&g);
        }
        Point { num, den }
    }
    pub fn dim(&self) -> usize {
        self.num.len()
    }
    pub fn coord(&self, i: usize) -> Q {
        Q::new(self.num[i].clone(), self.den.clone())
    }
    pub fn coords(&self) -> Vec<Q> {
        (0..self.dim()).map(|i| self.coord(i)).collect()
    }
    pub fn from_q(cs: &[Q]) -> Point {
        let mut den = Z::ONE;
        for c in cs {
            den = &den * &c.den.div_exact(&den.gcd(&c.den));
        }
        let num = cs.iter().map(|c| &c.num * &den.div_exact(&c.den)).collect();
        Point::new(num, den)
    }
    /// `<a, p>` as a rational.
    pub fn dot_q(&self, a: &[Z]) -> Q {
        Q::new(dot(a, &self.num), self.den.clone())
    }
    pub fn to_f64(&self) -> Vec<f64> {
        self.coords().iter().map(|c| c.to_f64()).collect()
    }
}

impl fmt::Display for Point {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "(")?;
        for (i, c) in self.coords().iter().enumerate() {
            if i > 0 {
                write!(f, ", ")?;
            }
            write!(f, "{c}")?;
        }
        write!(f, ")")
    }
}

/// Multiplies `a` by the LCM-free denominator scale to make a rational vector integer & primitive.
pub fn integerize(v: &[Q]) -> Vec<Z> {
    let mut den = Z::ONE;
    for c in v {
        den = &den * &c.den.div_exact(&den.gcd(&c.den));
    }
    let mut out: Vec<Z> = v.iter().map(|c| &c.num * &den.div_exact(&c.den)).collect();
    make_primitive(&mut out);
    out
}

pub fn factorial(n: usize) -> Z {
    let mut r = Z::ONE;
    for i in 2..=n {
        r = &r * &Z::from(i);
    }
    r
}

pub fn binomial(n: usize, k: usize) -> Z {
    if k > n {
        return Z::ZERO;
    }
    let mut r = Z::ONE;
    for i in 0..k {
        r = (&r * &Z::from(n - i)).div_exact(&Z::from(i + 1));
    }
    r
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn overflow_promotes_and_demotes() {
        let big = Z::S(i128::MAX);
        let s = &big + &Z::ONE;
        assert!(matches!(s, Z::B(_)));
        let back = &s - &Z::ONE;
        assert_eq!(back, big);
        assert_eq!((&s - &s), Z::ZERO);
        assert!(s > big);
    }

    #[test]
    fn rationals() {
        let a = Q::new(Z::from(6), Z::from(-4));
        assert_eq!(a, Q::new(Z::from(-3), Z::from(2)));
        assert_eq!(a.to_string(), "-3/2");
        assert_eq!(Q::parse("0.75").unwrap(), Q::new(Z::from(3), Z::from(4)));
        assert_eq!(Q::parse("-1/3").unwrap().mul(&Q::int(3)), Q::int(-1));
        assert_eq!(Q::parse("2e-2").unwrap(), Q::new(Z::from(1), Z::from(50)));
    }

    #[test]
    fn points() {
        let p = Point::new(vec![Z::from(2), Z::from(4)], Z::from(-6));
        assert_eq!(p.den, Z::from(3));
        assert_eq!(p.num, vec![Z::from(-1), Z::from(-2)]);
        assert_eq!(p.coord(1), Q::new(Z::from(-2), Z::from(3)));
    }
}
