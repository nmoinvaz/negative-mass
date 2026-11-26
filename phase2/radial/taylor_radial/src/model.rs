//! The model constants, as double-double numbers every lane embeds (so that the star in double-double, the J of the
//! perturbation in any lane and the quad-double checks solve the same discrete problem): lambda = 10240,
//! mu = lambda^2 / 4.8 and pi, each the nearest double-double.
use crate::num::{Num, PI_STR};
use multidouble::dd::DD;
use multidouble::qd::QD;

pub const LAM: f64 = 10240.0;

pub fn pi_dd() -> DD {
    QD::parse(PI_STR).unwrap().to_dd()
}

/// lambda^2 / 4.8 = lambda^2 5 / 24, in double-double
pub fn mu_dd() -> DD {
    DD::from_f64(LAM * LAM) * DD::from_f64(5.0) / DD::from_f64(24.0)
}

pub fn lam<T: Num>() -> T {
    T::from_f64(LAM)
}
pub fn mu<T: Num>() -> T {
    T::from_dd(mu_dd())
}
pub fn pi<T: Num>() -> T {
    T::from_dd(pi_dd())
}
