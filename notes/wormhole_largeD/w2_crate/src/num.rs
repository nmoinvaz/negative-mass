//! The scalar lanes: f64, double-double, quad-double and octo-double behind one trait.
use multidouble::dd::DD;
use multidouble::od::OD;
use multidouble::qd::QD;
use std::fmt::Debug;
use std::ops::{Add, Div, Mul, Neg, Sub};

/// 140 decimal digits of pi
pub const PI_STR: &str = "3.14159265358979323846264338327950288419716939937510582097494459230781640628620899862803482534211706798214808651328230664709384460955058223172535940812848111745028410270193852110555964462294895493038196";

pub trait Num:
    Copy + Add<Output = Self> + Sub<Output = Self> + Mul<Output = Self> + Div<Output = Self> + Neg<Output = Self> + Debug + Send + Sync + 'static
{
    const NAME: &'static str;
    /// unit roundoff of the lane
    const EPS: f64;
    fn zero() -> Self;
    fn one() -> Self;
    fn from_f64(x: f64) -> Self;
    fn from_dd(x: DD) -> Self;
    fn to_dd(self) -> DD;
    fn to_f64(self) -> f64;
    /// 1/x at the lane's precision
    fn recip(self) -> Self;
    fn sqrt(self) -> Self;
    fn exp(self) -> Self;
    fn ln(self) -> Self;
    /// decimal string at the lane's precision
    fn parse(s: &str) -> Self;
    fn pi() -> Self {
        Self::parse(PI_STR)
    }
    fn abs(self) -> Self {
        if self.to_f64() < 0.0 {
            -self
        } else {
            self
        }
    }
    /// a/b for small integers, at the lane's precision
    fn ratio(a: i64, b: i64) -> Self {
        Self::from_f64(a as f64) / Self::from_f64(b as f64)
    }
    /// decimal text with about `digits` significant digits
    fn sci(self, digits: usize) -> String;
}

impl Num for f64 {
    const NAME: &'static str = "f64";
    const EPS: f64 = 1.1102230246251565e-16;
    fn zero() -> f64 {
        0.0
    }
    fn one() -> f64 {
        1.0
    }
    fn from_f64(x: f64) -> f64 {
        x
    }
    fn from_dd(x: DD) -> f64 {
        x.hi + x.lo
    }
    fn to_dd(self) -> DD {
        DD::from_f64(self)
    }
    fn to_f64(self) -> f64 {
        self
    }
    fn recip(self) -> f64 {
        1.0 / self
    }
    fn sqrt(self) -> f64 {
        f64::sqrt(self)
    }
    fn exp(self) -> f64 {
        f64::exp(self)
    }
    fn ln(self) -> f64 {
        f64::ln(self)
    }
    fn parse(s: &str) -> f64 {
        s.parse().unwrap()
    }
    fn sci(self, digits: usize) -> String {
        format!("{:.*e}", digits.min(17).saturating_sub(1), self)
    }
}

impl Num for DD {
    const NAME: &'static str = "dd";
    const EPS: f64 = 1.232595164407831e-32;
    fn zero() -> DD {
        DD::ZERO
    }
    fn one() -> DD {
        DD::ONE
    }
    fn from_f64(x: f64) -> DD {
        DD::from_f64(x)
    }
    fn from_dd(x: DD) -> DD {
        x
    }
    fn to_dd(self) -> DD {
        self
    }
    fn to_f64(self) -> f64 {
        DD::to_f64(self)
    }
    fn recip(self) -> DD {
        DD::recip(self)
    }
    fn sqrt(self) -> DD {
        DD::sqrt(self)
    }
    fn exp(self) -> DD {
        DD::exp(self)
    }
    fn ln(self) -> DD {
        DD::ln(self)
    }
    fn parse(s: &str) -> DD {
        QD::parse(s).unwrap().to_dd()
    }
    fn sci(self, digits: usize) -> String {
        QD::from_dd(self).to_sci(digits.min(33))
    }
}

impl Num for QD {
    const NAME: &'static str = "qd";
    const EPS: f64 = 1.504632769052528e-64;
    fn zero() -> QD {
        QD::ZERO
    }
    fn one() -> QD {
        QD::ONE
    }
    fn from_f64(x: f64) -> QD {
        QD::from_f64(x)
    }
    fn from_dd(x: DD) -> QD {
        QD::from_dd(x)
    }
    fn to_dd(self) -> DD {
        QD::to_dd(self)
    }
    fn to_f64(self) -> f64 {
        QD::to_f64(self)
    }
    fn recip(self) -> QD {
        QD::recip(self)
    }
    fn sqrt(self) -> QD {
        QD::sqrt(self)
    }
    fn exp(self) -> QD {
        QD::exp(self)
    }
    fn ln(self) -> QD {
        QD::ln(self)
    }
    fn parse(s: &str) -> QD {
        QD::parse(s).unwrap()
    }
    fn sci(self, digits: usize) -> String {
        self.to_sci(digits.min(64))
    }
}

impl Num for OD {
    const NAME: &'static str = "od";
    const EPS: f64 = 2.3e-128;
    fn zero() -> OD {
        OD::from_f64(0.0)
    }
    fn one() -> OD {
        OD::from_f64(1.0)
    }
    fn from_f64(x: f64) -> OD {
        OD::from_f64(x)
    }
    fn from_dd(x: DD) -> OD {
        OD::from_dd(x)
    }
    fn to_dd(self) -> DD {
        OD::to_dd(self)
    }
    fn to_f64(self) -> f64 {
        OD::to_f64(self)
    }
    fn recip(self) -> OD {
        OD::recip(self)
    }
    fn sqrt(self) -> OD {
        OD::sqrt(self)
    }
    fn exp(self) -> OD {
        OD::exp(self)
    }
    fn ln(self) -> OD {
        OD::ln(self)
    }
    fn parse(s: &str) -> OD {
        OD::parse(s).unwrap()
    }
    fn sci(self, digits: usize) -> String {
        self.to_sci(digits.min(128))
    }
}
