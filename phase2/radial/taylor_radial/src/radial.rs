//! D(s2) of the radial perturbation system on a star: the regular 2-plane from r_s outward to r_m, the decaying
//! 2-plane from r_far inward to r_m from an arbitrary plane (the inward flow contracts every plane onto the pair that
//! grows inward, by e^{-(kappa_a + kappa_b)(r_far - r)}), both orthonormalized after every Taylor step; D = det[Y Z].
use crate::linalg::det;
use crate::num::Num;
use crate::plane::*;
use crate::ser::*;
use crate::table::*;
use multidouble::dd::DD;

/// exterior steps beyond r_cut: P = 0, m = M, Z = 1 / (1 - 2 M / r) and u from the exterior Riccati equation
///     u' = -u^2 - (2/r + f'/f) u + (1 - w^2/f)/f,   f = 1 - 2M/r,
/// integrated inward from r_far with u = -kappa - 1/r there, kappa = sqrt((1 - w^2/f)/f). The decaying solution is the
/// attractor of the inward flow, so the start's error is damped by e^{-2 kappa (r_far - r)}.
pub fn ext_steps<T: Num>(mass: T, omega: T, r_cut: f64, r_far: f64, tol: f64) -> Vec<StarStep<T>> {
    let one = T::one();
    let two = one + one;
    let inv = crate::plane::recips::<T>();
    let mut out = vec![];
    let rf = T::from_f64(r_far);
    let f_far = one - two * mass / rf;
    let kap = ((one - omega * omega / f_far) / f_far).sqrt();
    let mut u0 = -kap - one / rf;
    let mut x0 = r_far;
    while x0 > r_cut {
        let r = Ser::var(T::from_f64(x0));
        let rinv = recip(&r);
        let f = adds(&scl(&rinv, -two * mass), one);
        let finv = recip(&f);
        let c = add(&scl(&rinv, two), &mul(&scl(&powi(&rinv, 2), two * mass), &finv));
        let d = sub(&finv, &scl(&powi(&finv, 2), omega * omega));
        let mut uc = [T::zero(); KMAX + 1];
        uc[0] = u0;
        let mut u2 = [T::zero(); KMAX + 1];
        for k in 0..NS {
            let mut s2 = T::zero();
            for i in 0..=k {
                s2 = s2 + uc[i] * uc[k - i];
            }
            u2[k] = s2;
            let mut cu = T::zero();
            for i in 0..=k {
                cu = cu + c.get(i) * uc[k - i];
            }
            uc[k + 1] = (-(u2[k]) - cu + d.get(k)) * inv[k + 1];
        }
        let mut h = 4.0f64;
        for k in [NS - 1, NS] {
            let e = uc[k].to_f64().abs() / uc[0].to_f64().abs().max(1e-6);
            if e > 0.0 {
                h = h.min((tol / e).powf(1.0 / k as f64));
            }
        }
        h = h.min(0.4 * (x0 - 2.0 * mass.to_f64()));
        let x1 = if x0 - h <= r_cut { r_cut } else { x0 - h };
        let hh = T::from_f64(x1) - T::from_f64(x0);
        let mut s = uc[NS];
        for k in (0..NS).rev() {
            s = s * hh + uc[k];
        }
        let mut coef = [[T::zero(); KMAX + 1]; 4];
        coef[1] = uc;
        coef[2][0] = mass;
        for k in 0..=NS {
            coef[3][k] = finv.get(k);
        }
        out.push(StarStep { lo: x1, hi: x0, x0, coef });
        u0 = s;
        x0 = x1;
    }
    out.sort_by(|a, b| a.lo.partial_cmp(&b.lo).unwrap());
    out
}

pub struct Radial<T: Num> {
    pub table: JTable<T>,
    pub rs: f64,
    pub r_m: f64,
    pub r_max: f64,
    pub omega: T,
    /// a / alpha at r_s
    pub a_over_alpha: T,
    pub ctl: Ctl,
}

impl<T: Num> Radial<T> {
    pub fn new(table: JTable<T>, first: &StarStep<DD>, rs: f64, omega: T, pc: f64, ctl: Ctl) -> Self {
        let r_max = *table.x1.last().unwrap();
        // r_m: the first table step boundary where P < 1e-6 P_c
        let r_m = table.x0.iter().zip(&table.p_start).find(|(_, &p)| p < 1e-6 * pc).map(|(&x, _)| x).unwrap_or(r_max);
        let z0 = T::from_dd(first.coef[3][0]);
        let m0 = T::from_dd(first.coef[2][0]);
        let a2 = (T::one() - (T::one() + T::one()) * m0 / T::from_f64(rs)).recip();
        let a_over_alpha = (a2 * z0).sqrt();
        Radial { table, rs, r_m, r_max, omega, a_over_alpha, ctl }
    }

    fn j_matrix(&self, x0: f64, s2: T) -> [[T; 4]; 4] {
        let mut jets = Jets::<T>::new(&self.table.pattern());
        self.table.jets(x0, s2, 1, &mut jets);
        let mut m = [[T::zero(); 4]; 4];
        for en in &jets.ent {
            if en.n > 0 {
                m[en.i][en.j] = en.c[0];
            }
        }
        m
    }

    /// the regular plane's start at r_s
    pub fn regular_start(&self, s2: T) -> [[T; 4]; 2] {
        let rs = T::from_f64(self.rs);
        let j = self.j_matrix(self.rs, s2);
        let three = T::from_f64(3.0);
        let mut y = [[T::zero(); 4]; 2];
        for (c, (a0, g0)) in [(T::one(), T::zero()), (T::zero(), T::one())].into_iter().enumerate() {
            let chi_factor = if crate::mutate::on("start") { T::from_f64(2.0) } else { T::one() };
            let chi0 = rs * rs * rs * (g0 - (T::one() + T::one()) * self.omega * a0) * self.a_over_alpha / three * chi_factor;
            let f0 = j[1][0] * a0 + j[1][2] * g0 + j[1][3] * chi0;
            y[c] = [a0 + f0 * rs * rs / T::from_f64(6.0), f0 * rs / three, g0, chi0];
        }
        y
    }

    pub fn d(&self, s2: T, r_far: f64) -> T {
        let r_far = r_far.min(self.r_max);
        let mut y = self.regular_start(s2);
        let mut z = [[T::zero(); 4]; 2];
        let z0 = [[1.0, 0.3, 0.0, 0.2], [0.0, 0.5, 1.0, 0.1]];
        for c in 0..2 {
            for i in 0..4 {
                z[c][i] = T::from_f64(z0[c][i]);
            }
        }
        integrate::<T, 4, 2>(&self.table, &mut y, self.rs, self.r_m, s2, &self.ctl);
        integrate::<T, 4, 2>(&self.table, &mut z, r_far, self.r_m, s2, &self.ctl);
        let mut m = [[T::zero(); 4]; 4];
        for c in 0..2 {
            for i in 0..4 {
                m[i][c] = y[c][i];
                m[i][2 + c] = z[c][i];
            }
        }
        det(&m)
    }
}
