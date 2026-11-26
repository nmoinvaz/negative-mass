//! A star, its J table out to the far radius of a scan, its mass and Noether number, and the roots of D on a grid.
use crate::bvp::{StarCfg, StarInfo};
use crate::num::Num;
use crate::plane::Ctl;
use crate::radial::{ext_steps, Radial};
use crate::roots;
use crate::ser::*;
use crate::table::*;
use multidouble::dd::DD;

pub struct Prepared {
    pub eps: f64,
    pub cfg: StarCfg,
    pub info: StarInfo,
    pub tab: JTable<DD>,
    pub first: StarStep<DD>,
    pub r_cut: f64,
    pub r_far: f64,
    pub mass: DD,
    pub noether: DD,
    pub omega: DD,
    pub t_star: f64,
    pub t_table: f64,
}

/// kappa_min at s2 for frequency om: the smaller of sqrt(1 - (om +- sigma)^2), and at most k
pub fn kappa_min(om: f64, s2: f64) -> f64 {
    let sg = s2.max(0.0).sqrt();
    let kp = (1.0 - (om + sg).powi(2)).max(1e-12).sqrt();
    let km = (1.0 - (om - sg).powi(2)).max(1e-12).sqrt();
    kp.min(km).min((1.0 - om * om).sqrt())
}

/// 4 pi w r^2 phi0^2 a / alpha integrated analytically over the steps (series about each center)
pub fn noether(steps: &[StarStep<DD>], omega: DD, rs: f64, pc: DD, zc: DD) -> DD {
    let four_pi_w = crate::model::pi_dd().mul_f(4.0) * omega;
    let mut total = DD::ZERO;
    for st in steps {
        let (r, alpha, a, _phi0, _u) = bg_series(&st.coef, st.x0);
        let p = Ser::from_slice(&st.coef[0][..=NS]);
        let integrand = scl(&mul(&mul(&powi(&r, 2), &p), &mul(&a, &recip(&alpha))), four_pi_w);
        let t1 = DD::from_f64(st.hi) - DD::from_f64(st.x0);
        let t0 = DD::from_f64(st.lo) - DD::from_f64(st.x0);
        let (mut p1, mut p0) = (t1, t0);
        for k in 0..NS {
            total = total + integrand.get(k) * (p1 - p0) / DD::from_f64((k + 1) as f64);
            p1 = p1 * t1;
            p0 = p0 * t0;
        }
    }
    // the core 0 .. rs: (4 pi / 3) w P_c sqrt(Z_c) rs^3
    let rs3 = DD::from_f64(rs) * DD::from_f64(rs) * DD::from_f64(rs);
    total + crate::model::pi_dd().mul_f(4.0) / DD::from_f64(3.0) * omega * pc * zc.sqrt() * rs3
}

pub fn prepare(cfg: &StarCfg, info: StarInfo, t_star: f64, s2_hi: f64) -> Result<Prepared, String> {
    let t0 = std::time::Instant::now();
    if !info.converged {
        return Err(format!("star at omega {} not converged: {:?}", info.omega.to_f64(), info.fnorm));
    }
    let steps = info.steps.clone();
    let om = info.omega;
    let kk = info.k.to_f64();
    let r_t = info.r_t;
    let cut_target = r_t - 40.0 / kk;
    let r_cut = steps.iter().filter(|s| s.hi <= cut_target + 1e-9).map(|s| s.hi).fold(0.0, f64::max);
    let kmin = kappa_min(om.to_f64(), s2_hi);
    let r_far = r_t.max(r_cut + 40.0 / kmin).max(r_cut + 40.0 / kk);
    let mass = info.mass;
    let ext = ext_steps::<DD>(mass, om, r_cut, r_far, 1e-32);
    let noeth = noether(&steps, om, cfg.rs, info.pc, info.zc);
    let first = steps[0].clone();
    let mut all: Vec<StarStep<DD>> = steps.into_iter().filter(|s| s.hi <= r_cut + 1e-9).collect();
    all.extend(ext);
    let model = Model { omega: om, epsilon: cfg.eps };
    let tab = build_table(&all, &model, cfg.rs, r_far);
    let t_table = t0.elapsed().as_secs_f64();
    Ok(Prepared { eps: cfg.eps, cfg: cfg.clone(), info, tab, first, r_cut, r_far, mass, noether: noeth, omega: om, t_star, t_table })
}

pub fn prepare_k(cfg: &StarCfg, omega: DD, s2_hi: f64) -> Result<Prepared, String> {
    let t0 = std::time::Instant::now();
    let sol = crate::bvp::solve(cfg, omega, false)?;
    let t_star = t0.elapsed().as_secs_f64();
    prepare(cfg, sol.info(), t_star, s2_hi)
}

impl Prepared {
    pub fn radial<T: Num>(&self, order: usize, tol: f64) -> Radial<T> {
        let ctl = Ctl { order, tol, hmax: 4.0, kj: order };
        Radial::<T>::new(self.tab.to_lane::<T>(), &self.first, self.cfg.rs, T::from_dd(self.omega), self.info.pc.to_f64(), ctl)
    }
}

impl<T: Num> Radial<T> {
    /// the far radius for s2: r_m + factor / kappa_min, within the table
    pub fn rfar(&self, s2: f64, factor: f64) -> f64 {
        let om = self.omega.to_f64();
        (self.r_m + factor / kappa_min(om, s2)).min(self.r_max).max(self.r_m + 5.0)
    }
}

/// the declared root grid: 60 points log-spaced on [-neg_lo, -1e-7], 0, 30 log-spaced on [1e-7, 0.2 thr] and 30 linear on
/// [0.25 thr, 0.995 thr]
pub fn grid(neg_lo: f64, thr: f64, nneg: usize, npos: usize) -> Vec<f64> {
    let mut g = vec![];
    for i in 0..nneg {
        let t = i as f64 / (nneg - 1) as f64;
        g.push(-(10f64.powf(neg_lo.log10() + t * (-7.0 - neg_lo.log10()))));
    }
    g.push(0.0);
    let nlog = npos / 2;
    for i in 0..nlog {
        let t = i as f64 / (nlog - 1) as f64;
        g.push(10f64.powf(-7.0 + t * ((0.2 * thr).log10() + 7.0)));
    }
    let nlin = npos - nlog;
    for i in 0..nlin {
        let t = i as f64 / (nlin - 1) as f64;
        g.push(0.25 * thr + t * (0.995 - 0.25) * thr);
    }
    g
}

pub struct Found<T: Num> {
    pub roots: Vec<(T, usize, f64)>,
    pub vals: Vec<(f64, f64)>,
    pub evals: usize,
}

pub fn scan<T: Num>(rp: &Radial<T>, grid: &[f64], factor: f64) -> Found<T> {
    let mut f = |s: T| rp.d(s, rp.rfar(s.to_f64(), factor));
    let gt: Vec<T> = grid.iter().map(|&s| T::from_f64(s)).collect();
    let vals: Vec<T> = gt.iter().map(|&s| f(s)).collect();
    let mut found = Found { roots: vec![], vals: grid.iter().cloned().zip(vals.iter().map(|v| v.to_f64())).collect(), evals: grid.len() };
    for i in 0..grid.len() - 1 {
        if (vals[i].to_f64() < 0.0) != (vals[i + 1].to_f64() < 0.0) {
            let r = roots::find::<T>(&mut f, gt[i], gt[i + 1], vals[i], vals[i + 1]);
            found.evals += r.evals;
            found.roots.push((r.x, r.evals, r.fx.to_f64()));
        }
    }
    found
}
