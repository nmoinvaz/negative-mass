//! The phantom boson star of background.py as its own Taylor series, found by Newton shooting in the lane.
//!
//! State (P, u, m, Z) = (phi0^2, phi0'/phi0, m(r), alpha^-2): the right-hand sides of background.py's `odes` are
//! rational in these, so the Taylor recursion needs only Cauchy sums and one reciprocal. The solution is regular
//! at r = 0 and decays at infinity; the central data (P_c, Z_c) are the two unknowns:
//!   * P_c by bisection on the shot's class (u turns up: undershoot, u -> -inf: overshoot), to the last bit;
//!   * Z_c by the secant method on Z (1 - 2 m / r) - 1 at r = R_norm, where alpha is 1 at infinity.
//! Nothing from scipy's mesh enters.
use crate::num::Num;
use crate::ser::KMAX;

pub const NS: usize = 30;

#[derive(Clone, Copy, Debug)]
pub struct StarPars<T: Num> {
    pub w: T,
    pub eps: T,
    pub lam: T,
    pub mu: T,
    pub pi: T,
}

pub type Coef<T> = [[T; KMAX + 1]; 4];

#[derive(Clone)]
pub struct StarStep<T: Num> {
    pub x0: f64,
    pub x1: f64,
    pub coef: Coef<T>,
}

pub struct Run<T: Num> {
    pub steps: Vec<StarStep<T>>,
    /// +1: u turned up (P_c too small), -1: u fell to -inf (P_c too large), 0: reached r_end
    pub class: i32,
    pub x_end: f64,
    pub state_end: [T; 4],
}

fn conv<T: Num>(a: &[T], b: &[T], k: usize) -> T {
    let mut s = T::zero();
    for i in 0..=k {
        s = s + a[i] * b[k - i];
    }
    s
}

/// Taylor coefficients 0..=n of (P, u, m, Z) about x0 from the state there.
pub fn star_coeffs<T: Num>(x0: T, st: [T; 4], n: usize, p: &StarPars<T>, inv: &[T; KMAX + 1]) -> Coef<T> {
    let z = T::zero();
    let one = T::one();
    let two = one + one;
    let half = one / two;
    let nn = n + 2;
    let mut cp = vec![z; nn];
    let mut cu = vec![z; nn];
    let mut cm = vec![z; nn];
    let mut cz = vec![z; nn];
    cp[0] = st[0];
    cu[0] = st[1];
    cm[0] = st[2];
    cz[0] = st[3];
    // 1/r about x0, and the coefficients of r^2
    let r0inv = x0.recip();
    let mut rinv = vec![z; nn];
    let mut pw = r0inv;
    for (k, x) in rinv.iter_mut().enumerate() {
        *x = if k % 2 == 0 { pw } else { -pw };
        pw = pw * r0inv;
    }
    let x0sq = x0 * x0;
    let (mut bb, mut aa, mut pz, mut uu, mut uup, mut uupb) = (vec![z; nn], vec![z; nn], vec![z; nn], vec![z; nn], vec![z; nn], vec![z; nn]);
    let (mut p2, mut p3, mut rho, mut pr, mut apr, mut arho) = (vec![z; nn], vec![z; nn], vec![z; nn], vec![z; nn], vec![z; nn], vec![z; nn]);
    let (mut am1, mut nup, mut lamp, mut g, mut diff) = (vec![z; nn], vec![z; nn], vec![z; nn], vec![z; nn], vec![z; nn]);
    let w2 = p.w * p.w;
    let four = two + two;
    let fpie = four * p.pi * p.eps; // 4 pi eps
    let epie = fpie + fpie; // 8 pi eps
    let quarter_lam = p.lam / four;
    let mu6 = p.mu / T::from_f64(6.0);
    for k in 0..n {
        let d0 = if k == 0 { one } else { z };
        let mr = conv(&cm, &rinv, k);
        bb[k] = d0 - two * mr;
        if k == 0 {
            aa[0] = bb[0].recip();
        } else {
            let mut s = z;
            for i in 1..=k {
                s = s + bb[i] * aa[k - i];
            }
            aa[k] = -(s * aa[0]);
        }
        pz[k] = conv(&cp, &cz, k);
        uu[k] = conv(&cu, &cu, k);
        uup[k] = conv(&uu, &cp, k);
        uupb[k] = conv(&uup, &bb, k);
        p2[k] = conv(&cp, &cp, k);
        p3[k] = conv(&p2, &cp, k);
        let v = half * cp[k] - quarter_lam * p2[k] + mu6 * p3[k];
        let kin = half * (w2 * pz[k] + uupb[k]);
        rho[k] = kin + v;
        pr[k] = kin - v;
        let mut r2rho = x0sq * rho[k];
        if k >= 1 {
            r2rho = r2rho + two * x0 * rho[k - 1];
        }
        if k >= 2 {
            r2rho = r2rho + rho[k - 2];
        }
        let mp = fpie * r2rho;
        apr[k] = conv(&aa, &pr, k);
        arho[k] = conv(&aa, &rho, k);
        let mut rapr = x0 * apr[k];
        let mut rarho = x0 * arho[k];
        if k >= 1 {
            rapr = rapr + apr[k - 1];
            rarho = rarho + arho[k - 1];
        }
        am1[k] = aa[k] - d0;
        let am1r = conv(&am1, &rinv, k);
        nup[k] = am1r + epie * rapr;
        lamp[k] = -am1r + epie * rarho;
        g[k] = two * rinv[k] + half * nup[k] - half * lamp[k];
        let gu = conv(&g, &cu, k);
        diff[k] = (d0 - p.lam * cp[k] + p.mu * p2[k]) - w2 * cz[k];
        let adz = conv(&aa, &diff, k);
        let up = -uu[k] - gu + adz;
        let pp = two * conv(&cu, &cp, k);
        let zp = -conv(&nup, &cz, k);
        let ik = inv[k + 1];
        cp[k + 1] = pp * ik;
        cu[k + 1] = up * ik;
        cm[k + 1] = mp * ik;
        cz[k + 1] = zp * ik;
    }
    let mut out = [[z; KMAX + 1]; 4];
    for k in 0..=n {
        out[0][k] = cp[k];
        out[1][k] = cu[k];
        out[2][k] = cm[k];
        out[3][k] = cz[k];
    }
    out
}

/// The regular state at rs from the central data: leading order in r, whose error feeds the singular solution
/// only at O(rs^5).
pub fn central_state<T: Num>(pc: T, zc: T, rs: f64, p: &StarPars<T>) -> [T; 4] {
    let one = T::one();
    let two = one + one;
    let four = two + two;
    let half = one / two;
    let pc2 = pc * pc;
    let c = (one - p.lam * pc + p.mu * pc2) - p.w * p.w * zc;
    let v = half * pc - p.lam / four * pc2 + p.mu / T::from_f64(6.0) * pc2 * pc;
    let rho_c = half * p.w * p.w * pc * zc + v;
    let r = T::from_f64(rs);
    let three = two + one;
    // mutation start_term: the wrong leading coefficient of u, 1/2 instead of 1/3
    let ucoef = if crate::mutate::on("start_term") { r * c / two } else { r * c / three };
    [pc, ucoef, four * p.pi * p.eps * rho_c * r * r * r / three, zc]
}

/// Integrates from rs to r_end by Taylor steps of order NS. `store` keeps every step's coefficients.
pub fn shoot<T: Num>(pc: T, zc: T, p: &StarPars<T>, rs: f64, r_end: f64, store: bool, tol: f64, stop_on_class: bool) -> Run<T> {
    let inv = crate::plane::recips::<T>();
    let mut st = central_state(pc, zc, rs, p);
    let mut x0 = rs;
    let mut steps = Vec::new();
    let mut class = 0;
    let floor = [0.0, 1e-6, 0.0, 0.0];
    let mut taken = 0usize;
    while x0 < r_end {
        taken += 1;
        assert!(taken <= 200_000, "star: more than 200000 Taylor steps by r = {x0:e} (the shot has left the decaying solution)");
        let co = star_coeffs(T::from_f64(x0), st, NS, p, &inv);
        let mut h = 1.0f64;
        for k in [NS - 1, NS] {
            let mut e = 0.0f64;
            for f in 0..4 {
                let sc = co[f][0].to_f64().abs().max(floor[f]).max(1e-300);
                e = e.max(co[f][k].to_f64().abs() / sc);
            }
            if e > 0.0 {
                h = h.min((tol / e).powf(1.0 / k as f64));
            }
        }
        let x1 = if x0 + h >= r_end { r_end } else { x0 + h };
        let hh = T::from_f64(x1) - T::from_f64(x0);
        let mut nst = st;
        for f in 0..4 {
            let mut s = co[f][NS];
            for k in (0..NS).rev() {
                s = s * hh + co[f][k];
            }
            nst[f] = s;
        }
        if store {
            steps.push(StarStep { x0, x1, coef: co });
        }
        st = nst;
        x0 = x1;
        let u = st[1].to_f64();
        if stop_on_class {
            if x0 > 0.05 && u > 0.0 {
                class = 1;
                break;
            }
            if u < -50.0 {
                class = -1;
                break;
            }
        }
    }
    Run { steps, class, x_end: x0, state_end: st }
}

pub struct Star<T: Num> {
    pub pc: T,
    pub zc: T,
    pub steps: Vec<StarStep<T>>,
    pub r_end: f64,
}

pub struct SolveCfg {
    pub rs: f64,
    pub tol: f64,
    pub r_norm: f64,
    pub r_class: f64,
    pub verbose: bool,
}

fn bisect_pc<T: Num>(zc: T, lo0: T, hi0: T, p: &StarPars<T>, cfg: &SolveCfg) -> (T, T) {
    let (mut lo, mut hi) = (lo0, hi0);
    let half = T::from_f64(0.5);
    for _ in 0..400 {
        let mid = (lo + hi) * half;
        if ((hi - lo) / hi).to_f64().abs() < 4.0 * T::EPS {
            break;
        }
        let run = shoot(mid, zc, p, cfg.rs, cfg.r_class, false, cfg.tol, true);
        if run.class > 0 {
            lo = mid;
        } else if run.class < 0 {
            hi = mid;
        } else {
            // reached r_class without an event: the root is within the lane's resolution
            return (lo, hi);
        }
    }
    (lo, hi)
}

/// A bracket (lo: u turns up, hi: u falls) of the central amplitude at the central lapse zc.
fn bracket_pc<T: Num>(zc: T, guess: Option<T>, p: &StarPars<T>, cfg: &SolveCfg) -> (T, T) {
    let cls = |pc: T| shoot(pc, zc, p, cfg.rs, cfg.r_class, false, cfg.tol, true).class;
    if let Some(g) = guess {
        let mut w = 1e-4;
        for _ in 0..40 {
            let (lo, hi) = (g * T::from_f64(1.0 - w), g * T::from_f64(1.0 + w));
            if cls(lo) > 0 && cls(hi) < 0 {
                return (lo, hi);
            }
            w *= 2.0;
            if w > 0.5 {
                break;
            }
        }
    }
    // scan phi_c = 0.008 .. 0.03 as boson_star.find_star does
    let n = 220;
    let mut prev: Option<(T, i32)> = None;
    for i in 0..=n {
        let phic = 0.008 + (0.03 - 0.008) * i as f64 / n as f64;
        let pc = T::from_f64(phic) * T::from_f64(phic);
        let c = cls(pc);
        if let Some((pp, pcl)) = prev {
            if pcl > 0 && c < 0 {
                return (pp, pc);
            }
        }
        prev = Some((pc, c));
    }
    panic!("no bracket of the central amplitude for zc = {:?}", zc);
}

/// F2(zc): Z (1 - 2 m / r) - 1 at r_norm, with P_c the bisected root for that zc.
fn f2<T: Num>(zc: T, guess: &mut Option<T>, p: &StarPars<T>, cfg: &SolveCfg) -> (T, T) {
    let (lo, hi) = bracket_pc(zc, *guess, p, cfg);
    let (lo, hi) = bisect_pc(zc, lo, hi, p, cfg);
    let pc = (lo + hi) * T::from_f64(0.5);
    *guess = Some(pc);
    let run = shoot(pc, zc, p, cfg.rs, cfg.r_norm, false, cfg.tol, false);
    let r = T::from_f64(run.x_end);
    let val = run.state_end[3] * (T::one() - (T::one() + T::one()) * run.state_end[2] / r) - T::one();
    (pc, val)
}

/// Solves for (P_c, Z_c); returns them. `z_guess` starts the secant method.
pub fn solve_central<T: Num>(p: &StarPars<T>, z_guess: f64, cfg: &SolveCfg) -> (T, T, usize) {
    let mut guess = None;
    let z0 = T::from_f64(z_guess);
    let z1 = T::from_f64(z_guess + 0.01);
    let (_, f0) = f2(z0, &mut guess, p, cfg);
    let (mut pc1, mut f1) = f2(z1, &mut guess, p, cfg);
    let (mut za, mut zb, mut fa, mut fb) = (z0, z1, f0, f1);
    let mut iters = 2;
    for _ in 0..30 {
        if cfg.verbose {
            eprintln!("star: z {:?} f2 {:e}", zb.to_f64(), fb.to_f64());
        }
        if fb.to_f64().abs() < 50.0 * T::EPS {
            break;
        }
        let zn = zb - fb * (zb - za) / (fb - fa);
        za = zb;
        fa = fb;
        zb = zn;
        let (pc, f) = f2(zb, &mut guess, p, cfg);
        pc1 = pc;
        f1 = f;
        fb = f1;
        iters += 1;
    }
    (pc1, zb, iters)
}

/// The star's steps from rs to r_end for the central data.
pub fn build_star<T: Num>(pc: T, zc: T, p: &StarPars<T>, cfg: &SolveCfg, r_end: f64) -> Star<T> {
    let run = shoot(pc, zc, p, cfg.rs, r_end, true, cfg.tol, false);
    Star { pc, zc, steps: run.steps, r_end: run.x_end }
}

impl<T: Num> Star<T> {
    /// (P, u, m, Z) at r by Horner on the step holding r
    pub fn state(&self, r: f64) -> [T; 4] {
        let idx = match self.steps.binary_search_by(|s| s.x0.partial_cmp(&r).unwrap()) {
            Ok(i) => i,
            Err(i) => i.saturating_sub(1),
        };
        let s = &self.steps[idx];
        let h = T::from_f64(r) - T::from_f64(s.x0);
        let mut out = [T::zero(); 4];
        for f in 0..4 {
            let mut v = s.coef[f][NS];
            for k in (0..NS).rev() {
                v = v * h + s.coef[f][k];
            }
            out[f] = v;
        }
        out
    }
    /// (alpha, a, phi0, u) at r, the arguments of J_and_h
    pub fn bg(&self, r: f64) -> [T; 4] {
        let s = self.state(r);
        let rr = T::from_f64(r);
        let alpha = (s[3]).recip().sqrt();
        let a = (T::one() - (T::one() + T::one()) * s[2] / rr).recip().sqrt();
        [alpha, a, s[0].sqrt(), s[1]]
    }
}
