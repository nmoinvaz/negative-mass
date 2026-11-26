//! The system of E2 in l = s ln r (see the raw README, "The equation"):
//!   Y_l = eps w,   F w_l = A Y + B w,   F = 1 - E,  E = e^{-l},  P = e^{2 eps l},
//!   A = (1 + eps Om) P - (1/2 + eps) + eps E,   B = (1/2 + 2 eps) - 2 eps E.
//! The planted system replaces A by eps (Om - Om*) + s kappa u [F (1 + 2 kappa u) + B], u = E / (1 + kappa F), kappa =
//! min(eps, 1/1000) (radius of convergence of Y* at the throat at least ln 1001 = 6.9 > 2 pi), whose solution at Om = Om* is Y* = 1 / (1 + kappa F) (analytic at the throat, tends to a constant, so the other solution
//! e^{(1/2 + 2 eps) l} swamps it outward). Series are in tau with l = l0 + dh tau, so that coefficients stay near
//! unit size however long the step.
use crate::mutate;
use crate::num::Num;

#[derive(Clone)]
pub struct Sys<T: Num> {
    pub eps: T,
    pub s: T,
    pub planted: bool,
    pub kappa: T,
    pub om_star: T,
    /// 1/k in the lane, k = 0..=kmax (entry 0 unused)
    pub rk: Vec<T>,
}

/// Taylor coefficients (in tau) of F, A, dA/dOm and B
pub struct Coef<T> {
    pub f: Vec<T>,
    pub a: Vec<T>,
    pub ad: Vec<T>,
    pub b: Vec<T>,
}

impl<T: Num> Sys<T> {
    pub fn new(d: u64, planted: bool, kmax: usize) -> Self {
        let s = T::from_f64(2.0 * (d as f64 - 3.0));
        let rk = (0..=kmax + 2)
            .map(|k| {
                if k == 0 {
                    T::zero()
                } else if mutate::on("recip_f64") {
                    T::from_f64(1.0 / k as f64)
                } else {
                    T::one() / T::from_f64(k as f64)
                }
            })
            .collect();
        Sys { eps: T::one() / s, s, planted, kappa: if d < 503 { T::ratio(1, 1000) } else { T::one() / s }, om_star: T::ratio(5, 2), rk }
    }

    fn half() -> T {
        T::ratio(1, 2)
    }

    pub fn coef(&self, l0: T, dh: T, om: T, kk: usize) -> Coef<T> {
        let eps = self.eps;
        let two_eps = eps + eps;
        let half = Self::half();
        let e0 = (-l0).exp();
        let p0 = (two_eps * l0).exp();
        let ser = |x0: T, rate: T| -> Vec<T> {
            let mut v = vec![T::zero(); kk + 1];
            v[0] = x0;
            for k in 1..=kk {
                v[k] = v[k - 1] * rate * self.rk[k];
            }
            v
        };
        let e = ser(e0, -dh);
        let p = ser(p0, two_eps * dh);
        let mut f = vec![T::zero(); kk + 1];
        let mut a = vec![T::zero(); kk + 1];
        let mut ad = vec![T::zero(); kk + 1];
        let mut b = vec![T::zero(); kk + 1];
        f[0] = T::one() - e[0];
        for k in 1..=kk {
            f[k] = -e[k];
        }
        // eps in A, for the control wrong_eps
        let eps_a = if mutate::on("wrong_eps") { eps * T::from_f64(1.0 + 1e-6) } else { eps };
        if !self.planted {
            let c1 = T::one() + eps_a * om;
            for k in 0..=kk {
                a[k] = c1 * p[k] + eps_a * e[k];
                if mutate::on("wrong_A_sign") {
                    // the E/2 term of A with the wrong sign: -(1/2+eps)(1-E) + E/2 instead of ... - E/2
                    a[k] = c1 * p[k] + (T::one() + eps_a) * e[k];
                }
            }
            a[0] = a[0] - (half + eps_a);
        } else {
            let ka = if mutate::on("wrong_planted") { self.kappa * T::from_f64(1.0 + 1e-3) } else { self.kappa };
            // u = E / G, G = 1 + kappa F, A* = s kappa u [F (1 + 2 kappa u) + B]
            let g1 = half + two_eps;
            let mut bs = vec![T::zero(); kk + 1];
            for k in 0..=kk {
                bs[k] = -two_eps * e[k];
            }
            bs[0] = bs[0] + g1;
            let mut gs = vec![T::zero(); kk + 1];
            for k in 0..=kk {
                gs[k] = ka * f[k];
            }
            gs[0] = gs[0] + T::one();
            let r = cinv(&gs, kk);
            let u = cmul(&e, &r, kk);
            let mut t = vec![T::zero(); kk + 1];
            for k in 0..=kk {
                t[k] = T::from_f64(2.0) * ka * u[k];
            }
            t[0] = t[0] + T::one();
            let ft = cmul(&f, &t, kk);
            let mut br = vec![T::zero(); kk + 1];
            for k in 0..=kk {
                br[k] = ft[k] + bs[k];
            }
            let a_star = cmul(&u, &br, kk);
            let sk = self.s * ka;
            for k in 0..=kk {
                a[k] = sk * a_star[k];
            }
            a[0] = a[0] + eps_a * (om - self.om_star);
        }
        for k in 0..=kk {
            ad[k] = if self.planted { T::zero() } else { eps * p[k] };
        }
        if self.planted {
            ad[0] = eps;
        }
        // B = (1/2 + 2 eps) - 2 eps E; the control wrong_B drops E/2 from (1/2 + 2 eps)(1 - E) + E/2
        for k in 0..=kk {
            b[k] = if mutate::on("wrong_B") { -(half + two_eps) * e[k] } else { -two_eps * e[k] };
        }
        b[0] = b[0] + half + two_eps;
        Coef { f, a, ad, b }
    }

    /// (Y*, w*) at l for the planted system; w = s Y_l
    pub fn exact(&self, l: T) -> (T, T) {
        let e = (-l).exp();
        let g = T::one() + self.kappa * (T::one() - e);
        let y = g.recip();
        (y, -(self.s * self.kappa * e * y * y))
    }
}

/// Cauchy product of two series, orders 0..=kk
pub fn cmul<T: Num>(a: &[T], b: &[T], kk: usize) -> Vec<T> {
    (0..=kk)
        .map(|k| (0..=k).fold(T::zero(), |s, j| s + a[j] * b[k - j]))
        .collect()
}

/// Reciprocal of a series with g[0] != 0
pub fn cinv<T: Num>(g: &[T], kk: usize) -> Vec<T> {
    let i0 = g[0].recip();
    let mut r = vec![T::zero(); kk + 1];
    r[0] = i0;
    for k in 1..=kk {
        let mut s = T::zero();
        for j in 1..=k {
            s = s + g[j] * r[k - j];
        }
        r[k] = -(s * i0);
    }
    r
}
