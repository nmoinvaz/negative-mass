//! Complex numbers over the crate's lanes: a pair (re, im) of lane values that implements `Num` (negative-mass #188).
//!
//! Every generic piece of the crate (series, planes, determinants, the projector, the Riccati integration) runs on `Cx<T>`
//! unchanged. A number with im = 0 behaves bit for bit as its real part: the products are (ac - bd, ad + bc) with exact
//! zeros, and division, reciprocal, square root, logarithm and exponential take the real-axis branch of the lane when im = 0,
//! so the complex lane at a real s reproduces the real lane's bits (check R2).
//!
//! `ln`, `exp`, `to_f64` and `to_dd` of a number with im != 0 panic: the R diagonal of the plane stepper is real and positive,
//! and a silent drop of an imaginary part is the failure the panics catch. Magnitudes and real parts are `mag` and `re_f64`.
use crate::num::Num;
use multidouble::dd::DD;
use std::ops::{Add, Div, Mul, Neg, Sub};
use std::sync::atomic::{AtomicBool, Ordering::Relaxed};

// the deliberate errors of the controls (`mutate=NAME`), read once by `init`
static BILINEAR: AtomicBool = AtomicBool::new(false);
static REDOT: AtomicBool = AtomicBool::new(false);
static RECIP_GENERAL: AtomicBool = AtomicBool::new(false);
static MUL_SIGN: AtomicBool = AtomicBool::new(false);

/// Reads the mutation names relevant to the complex type; call once after `mutate::set`.
pub fn init() {
    BILINEAR.store(crate::mutate::on("bilinear"), Relaxed);
    REDOT.store(crate::mutate::on("redot"), Relaxed);
    RECIP_GENERAL.store(crate::mutate::on("cx_recip_general"), Relaxed);
    MUL_SIGN.store(crate::mutate::on("cx_mul_sign"), Relaxed);
}

#[derive(Clone, Copy, Debug)]
pub struct Cx<T: Num> {
    pub re: T,
    pub im: T,
}

impl<T: Num> Cx<T> {
    pub fn new(re: T, im: T) -> Self {
        Cx { re, im }
    }
    pub fn is_real(self) -> bool {
        self.im.to_f64() == 0.0
    }
}

impl<T: Num> Add for Cx<T> {
    type Output = Cx<T>;
    fn add(self, o: Cx<T>) -> Cx<T> {
        Cx { re: self.re + o.re, im: self.im + o.im }
    }
}

impl<T: Num> Sub for Cx<T> {
    type Output = Cx<T>;
    fn sub(self, o: Cx<T>) -> Cx<T> {
        Cx { re: self.re - o.re, im: self.im - o.im }
    }
}

impl<T: Num> Neg for Cx<T> {
    type Output = Cx<T>;
    fn neg(self) -> Cx<T> {
        Cx { re: -self.re, im: -self.im }
    }
}

impl<T: Num> Mul for Cx<T> {
    type Output = Cx<T>;
    #[inline]
    fn mul(self, o: Cx<T>) -> Cx<T> {
        let (a, b, c, d) = (self.re, self.im, o.re, o.im);
        // mutation cx_mul_sign: the sign of the imaginary product in the real part
        let re = if MUL_SIGN.load(Relaxed) { a * c + b * d } else { a * c - b * d };
        Cx { re, im: a * d + b * c }
    }
}

impl<T: Num> Div for Cx<T> {
    type Output = Cx<T>;
    fn div(self, o: Cx<T>) -> Cx<T> {
        if o.is_real() && !RECIP_GENERAL.load(Relaxed) {
            return Cx { re: self.re / o.re, im: self.im / o.re };
        }
        let n = o.re * o.re + o.im * o.im;
        Cx { re: (self.re * o.re + self.im * o.im) / n, im: (self.im * o.re - self.re * o.im) / n }
    }
}

impl<T: Num> Num for Cx<T> {
    const NAME: &'static str = "cx";
    const EPS: f64 = T::EPS;
    fn zero() -> Self {
        Cx { re: T::zero(), im: T::zero() }
    }
    fn one() -> Self {
        Cx { re: T::one(), im: T::zero() }
    }
    fn from_f64(x: f64) -> Self {
        Cx { re: T::from_f64(x), im: T::zero() }
    }
    fn from_dd(x: DD) -> Self {
        Cx { re: T::from_dd(x), im: T::zero() }
    }
    fn to_dd(self) -> DD {
        assert!(self.is_real(), "Cx::to_dd of a number with a nonzero imaginary part");
        self.re.to_dd()
    }
    fn to_f64(self) -> f64 {
        assert!(self.is_real(), "Cx::to_f64 of a number with a nonzero imaginary part: use mag() or re_f64()");
        self.re.to_f64()
    }
    fn recip(self) -> Self {
        if self.is_real() && !RECIP_GENERAL.load(Relaxed) {
            return Cx { re: self.re.recip(), im: T::zero() };
        }
        let n = self.re * self.re + self.im * self.im;
        let inv = n.recip();
        Cx { re: self.re * inv, im: -(self.im * inv) }
    }
    /// principal square root (Re >= 0)
    fn sqrt(self) -> Self {
        let (a, b) = (self.re, self.im);
        if self.is_real() {
            return if a.to_f64() >= 0.0 { Cx { re: a.sqrt(), im: T::zero() } } else { Cx { re: T::zero(), im: (-a).sqrt() } };
        }
        let m = (a * a + b * b).sqrt();
        let t = ((m + a.abs()) * T::from_f64(0.5)).sqrt();
        let two_t = t + t;
        if a.to_f64() >= 0.0 {
            Cx { re: t, im: b / two_t }
        } else {
            Cx { re: b.abs() / two_t, im: if b.to_f64() < 0.0 { -t } else { t } }
        }
    }
    fn exp(self) -> Self {
        assert!(self.is_real(), "Cx::exp of a number with a nonzero imaginary part");
        Cx { re: self.re.exp(), im: T::zero() }
    }
    fn ln(self) -> Self {
        // mutation bilinear: the "norm" sqrt(sum y y) is complex; the control takes ln |z| so that it runs and its target reads a number
        if BILINEAR.load(Relaxed) {
            return Cx { re: (self.re * self.re + self.im * self.im).ln() * T::from_f64(0.5), im: T::zero() };
        }
        assert!(self.is_real() && self.re.to_f64() > 0.0, "Cx::ln only of a positive real number (the R diagonal)");
        Cx { re: self.re.ln(), im: T::zero() }
    }
    fn parse(s: &str) -> Self {
        Cx { re: T::parse(s), im: T::zero() }
    }
    fn sci(self, digits: usize) -> String {
        format!("({}, {})", self.re.sci(digits), self.im.sci(digits))
    }
    fn conj(self) -> Self {
        // mutation redot: the imaginary part of the projector dropped; bilinear: no conjugate at all
        if REDOT.load(Relaxed) {
            return Cx { re: self.re, im: T::zero() };
        }
        if BILINEAR.load(Relaxed) {
            return self;
        }
        Cx { re: self.re, im: -self.im }
    }
    fn abs2(self) -> Self {
        if BILINEAR.load(Relaxed) || self.is_real() {
            return self * self;
        }
        Cx { re: self.re * self.re + self.im * self.im, im: T::zero() }
    }
    fn mag(self) -> f64 {
        if self.is_real() {
            self.re.to_f64().abs()
        } else {
            self.re.to_f64().hypot(self.im.to_f64())
        }
    }
    fn re_f64(self) -> f64 {
        self.re.to_f64()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn bits(x: f64) -> u64 {
        x.to_bits()
    }

    #[test]
    fn real_embedding_is_exact_f64() {
        let xs = [0.3_f64, -2.75, 1e-7, 123456.789, 0.1];
        for &a in &xs {
            for &c in &xs {
                let (x, y) = (Cx::<f64>::from_f64(a), Cx::<f64>::from_f64(c));
                assert_eq!(bits((x * y).re), bits(a * c));
                assert_eq!(bits((x / y).re), bits(a / c));
                assert_eq!(bits((x + y).re), bits(a + c));
                assert_eq!(bits(x.recip().re), bits(1.0 / a));
                assert!((x * y).im == 0.0 && (x / y).im == 0.0);
            }
            if a > 0.0 {
                assert_eq!(bits(Cx::<f64>::from_f64(a).sqrt().re), bits(a.sqrt()));
                assert_eq!(bits(Cx::<f64>::from_f64(a).ln().re), bits(a.ln()));
                assert_eq!(bits(Cx::<f64>::from_f64(a).abs2().re), bits(a * a));
            }
        }
    }

    #[test]
    fn real_embedding_is_exact_dd() {
        let xs = [DD::from_f64(0.3) / DD::from_f64(7.0), DD::from_f64(-2.75), DD::from_f64(1e-7) / DD::from_f64(3.0)];
        for &a in &xs {
            for &c in &xs {
                let (x, y) = (Cx::<DD>::from_dd(a), Cx::<DD>::from_dd(c));
                let (p, q) = ((x * y).re, a * c);
                assert!(p.hi.to_bits() == q.hi.to_bits() && p.lo.to_bits() == q.lo.to_bits());
                let (p, q) = ((x / y).re, a / c);
                assert!(p.hi.to_bits() == q.hi.to_bits() && p.lo.to_bits() == q.lo.to_bits());
                let (p, q) = (x.recip().re, a.recip());
                assert!(p.hi.to_bits() == q.hi.to_bits() && p.lo.to_bits() == q.lo.to_bits());
            }
        }
    }

    #[test]
    fn sqrt_squares_back_and_divides() {
        for &(a, b) in &[(3.0_f64, 4.0), (-3.0, 4.0), (-3.0, -4.0), (3.0, -4.0), (0.0, 2.0), (-5.0, 0.0), (1e-3, 7.0)] {
            let z = Cx::<f64>::new(a, b);
            let r = z.sqrt();
            let back = r * r;
            assert!((back.re - a).abs() < 1e-14 * (a.abs() + b.abs() + 1.0) && (back.im - b).abs() < 1e-14 * (a.abs() + b.abs() + 1.0), "{a} {b}");
            assert!(r.re >= 0.0);
            let one = z / z;
            assert!((one.re - 1.0).abs() < 1e-15 && one.im.abs() < 1e-15);
            let w = z * z.recip();
            assert!((w.re - 1.0).abs() < 1e-15 && w.im.abs() < 1e-15);
        }
    }

    #[test]
    fn conj_abs2_mag() {
        let z = Cx::<f64>::new(3.0, -4.0);
        assert_eq!(z.conj().im, 4.0);
        assert_eq!(z.abs2().re, 25.0);
        assert_eq!(z.abs2().im, 0.0);
        assert_eq!(z.mag(), 5.0);
    }
}
