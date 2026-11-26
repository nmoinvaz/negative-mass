//! Truncated power series in t with at most KMAX coefficients, over any lane of `Num`.
//! `n` counts the leading coefficients that may be nonzero (0: the zero series), so a constant or a linear
//! function costs little in a product.
use crate::num::Num;

pub const KMAX: usize = 32;

#[derive(Clone, Copy, Debug)]
pub struct Ser<T: Num> {
    pub c: [T; KMAX],
    pub n: usize,
}

impl<T: Num> Ser<T> {
    pub fn zero() -> Self {
        Ser { c: [T::zero(); KMAX], n: 0 }
    }
    pub fn cst(x: T) -> Self {
        let mut s = Self::zero();
        s.c[0] = x;
        s.n = 1;
        s
    }
    /// x0 + t
    pub fn var(x0: T) -> Self {
        let mut s = Self::cst(x0);
        s.c[1] = T::one();
        s.n = 2;
        s
    }
    pub fn from_slice(v: &[T]) -> Self {
        let mut s = Self::zero();
        let n = v.len().min(KMAX);
        s.c[..n].copy_from_slice(&v[..n]);
        s.n = n;
        s
    }
    pub fn get(&self, k: usize) -> T {
        if k < self.n {
            self.c[k]
        } else {
            T::zero()
        }
    }
    /// value at t
    pub fn eval(&self, t: T) -> T {
        let mut s = T::zero();
        for k in (0..self.n).rev() {
            s = s * t + self.c[k];
        }
        s
    }
}

pub fn add<T: Num>(a: &Ser<T>, b: &Ser<T>) -> Ser<T> {
    let n = a.n.max(b.n);
    let mut o = Ser::zero();
    for k in 0..n {
        o.c[k] = a.get(k) + b.get(k);
    }
    o.n = n;
    o
}

pub fn sub<T: Num>(a: &Ser<T>, b: &Ser<T>) -> Ser<T> {
    let n = a.n.max(b.n);
    let mut o = Ser::zero();
    for k in 0..n {
        o.c[k] = a.get(k) - b.get(k);
    }
    o.n = n;
    o
}

pub fn neg<T: Num>(a: &Ser<T>) -> Ser<T> {
    let mut o = *a;
    for k in 0..a.n {
        o.c[k] = -a.c[k];
    }
    o
}

/// a * x for a scalar x
pub fn scl<T: Num>(a: &Ser<T>, x: T) -> Ser<T> {
    let mut o = Ser::zero();
    for k in 0..a.n {
        o.c[k] = a.c[k] * x;
    }
    o.n = a.n;
    o
}

/// a + x for a scalar x
pub fn adds<T: Num>(a: &Ser<T>, x: T) -> Ser<T> {
    let mut o = *a;
    if o.n == 0 {
        o.n = 1;
        o.c[0] = x;
    } else {
        o.c[0] = o.c[0] + x;
    }
    o
}

pub fn mul<T: Num>(a: &Ser<T>, b: &Ser<T>) -> Ser<T> {
    let mut o = Ser::zero();
    if a.n == 0 || b.n == 0 {
        return o;
    }
    let n = (a.n + b.n - 1).min(KMAX);
    for k in 0..n {
        let lo = k.saturating_sub(b.n - 1);
        let hi = k.min(a.n - 1);
        let mut s = T::zero();
        for i in lo..=hi {
            s = s + a.c[i] * b.c[k - i];
        }
        o.c[k] = s;
    }
    o.n = n;
    o
}

/// 1/a; a.c[0] must be nonzero. The reciprocal of the leading term is formed in the lane's precision.
pub fn recip<T: Num>(a: &Ser<T>) -> Ser<T> {
    assert!(a.n >= 1, "recip of the zero series");
    let mut o = Ser::zero();
    let b0 = a.c[0].recip();
    o.c[0] = b0;
    if a.n == 1 {
        o.n = 1;
        return o;
    }
    for k in 1..KMAX {
        let mut s = T::zero();
        for i in 1..=k.min(a.n - 1) {
            s = s + a.c[i] * o.c[k - i];
        }
        o.c[k] = -(b0 * s);
    }
    o.n = KMAX;
    o
}

pub fn div<T: Num>(a: &Ser<T>, b: &Ser<T>) -> Ser<T> {
    mul(a, &recip(b))
}

pub fn powi<T: Num>(a: &Ser<T>, e: i32) -> Ser<T> {
    if e < 0 {
        return recip(&powi(a, -e));
    }
    let mut r = Ser::cst(T::one());
    let mut b = *a;
    let mut e = e as u32;
    let mut first = true;
    while e > 0 {
        if e & 1 == 1 {
            r = if first { b } else { mul(&r, &b) };
            first = false;
        }
        e >>= 1;
        if e > 0 {
            b = mul(&b, &b);
        }
    }
    if first {
        Ser::cst(T::one())
    } else {
        r
    }
}

pub fn sqrt<T: Num>(a: &Ser<T>) -> Ser<T> {
    assert!(a.n >= 1, "sqrt of the zero series");
    let mut o = Ser::zero();
    let b0 = a.c[0].sqrt();
    o.c[0] = b0;
    if a.n == 1 {
        o.n = 1;
        return o;
    }
    let inv2b0 = (b0 + b0).recip();
    for k in 1..KMAX {
        let mut s = T::zero();
        for i in 1..k {
            s = s + o.c[i] * o.c[k - i];
        }
        o.c[k] = (a.get(k) - s) * inv2b0;
    }
    o.n = KMAX;
    o
}

/// The coefficients of f(x0 + delta + t) from those of f(x0 + t), k = 0..kk, by repeated synthetic division
/// (Horner's scheme); c holds nb coefficients and is overwritten.
pub fn taylor_shift<T: Num>(c: &mut [T], delta: T, kk: usize) {
    let nb = c.len();
    for i in 0..kk.min(nb) {
        for m in (i..nb - 1).rev() {
            c[m] = c[m] + delta * c[m + 1];
        }
    }
}
