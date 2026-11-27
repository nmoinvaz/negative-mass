//! J(r; s2) = J0(r) + s2 J1(r) as Taylor tables along the star, and the star's steps they are built from.
//!
//! The star is the multishoot layout's solution (`StarStep` per Taylor step of every segment, series of (P, u, m, Z)
//! about the step's source point), followed by the vacuum exterior beyond r_cut (P = 0, m = M, Z = 1/(1 - 2M/r)) with u
//! from the exterior Riccati equation integrated inward. J is linear in s2 (checked at the table build at s2 = 1.5), so
//! one table serves every s2. A plane step takes the jets of J's entries about its start by shifting the table series,
//! and a table step is at most WIDTH_FRACTION of the width at which J's own series has converged to TAIL (#186's defect:
//! a shift past J's radius of convergence).
use crate::num::Num;
use crate::plane::{JetSource, Jets};
use crate::rgen::{j_series, Pars};
use crate::ser::*;
use multidouble::dd::DD;

/// coefficients kept per entry
pub const NS: usize = 30;
pub const NB: usize = NS;

#[derive(Clone)]
pub struct StarStep<T: Num> {
    pub lo: f64,
    pub hi: f64,
    /// the series' center
    pub x0: f64,
    /// coefficients 0..=NS of (P, u, m, Z)
    pub coef: [[T; KMAX + 1]; 4],
}

#[derive(Clone, Copy, Debug)]
pub struct Model {
    pub omega: DD,
    pub epsilon: f64,
}

pub fn pars<T: Num>(m: &Model, s2: T) -> Pars<T> {
    let omega = if crate::mutate::on("jomega") { T::from_dd(m.omega) * T::from_f64(1.0 + 1e-6) } else { T::from_dd(m.omega) };
    let mu = if crate::mutate::on("jmu") { crate::model::mu::<T>() * T::from_f64(1.0 + 1e-6) } else { crate::model::mu::<T>() };
    Pars { omega, epsilon: T::from_f64(m.epsilon), s2, lambda: crate::model::lam::<T>(), mu, pi: crate::model::pi::<T>() }
}

/// series of (r, alpha, a, phi0, u) about b from the star's series about its center
pub fn bg_series<T: Num>(co: &[[T; KMAX + 1]; 4], b: f64) -> (Ser<T>, Ser<T>, Ser<T>, Ser<T>, Ser<T>) {
    let r = Ser::var(T::from_f64(b));
    let p = Ser::from_slice(&co[0][..=NS]);
    let u = Ser::from_slice(&co[1][..=NS]);
    let m = Ser::from_slice(&co[2][..=NS]);
    let z = Ser::from_slice(&co[3][..=NS]);
    let rinv = recip(&r);
    let alpha = sqrt(&recip(&z));
    let bb = adds(&scl(&mul(&m, &rinv), -(T::one() + T::one())), T::one());
    let a = sqrt(&recip(&bb));
    let phi0 = if p.c[0].to_f64() == 0.0 { Ser::zero() } else { sqrt(&p) };
    (r, alpha, a, phi0, u)
}

fn width_of(sers: &[[DD; NB]], hmax: f64, tail: f64) -> f64 {
    let ok = |h: f64| -> bool {
        sers.iter().all(|c| {
            let mut sum = 0.0;
            let mut hk = 1.0;
            let mut last = 0.0f64;
            for (k, x) in c.iter().enumerate() {
                let t = x.to_f64().abs() * hk;
                sum += t;
                if k + 2 >= NB {
                    last = last.max(t);
                }
                hk *= h;
            }
            sum == 0.0 || last <= tail * sum
        })
    };
    if ok(hmax) {
        return hmax;
    }
    let (mut lo, mut hi) = (0.0, hmax);
    for _ in 0..60 {
        let mid = 0.5 * (lo + hi);
        if ok(mid) {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    lo
}

pub const TAIL: f64 = 1e-25;
pub const WIDTH_FRACTION: f64 = 0.3;

pub struct JTable<T: Num> {
    pub pattern: Vec<(usize, usize)>,
    pub x0: Vec<f64>,
    pub x1: Vec<f64>,
    pub c0: Vec<Vec<[T; NB]>>,
    pub c1: Vec<Vec<[T; NB]>>,
    /// largest relative |J(s2 = 1.5) - (J0 + 1.5 J1)| at order 0
    pub linear_defect: f64,
    /// P / P_c at each table step's start, for r_m
    pub p_start: Vec<f64>,
}

type JM = [[Ser<DD>; 4]; 4];

/// The table over [rmin, rmax] from steps sorted by position (series about `x0`, valid on [lo, hi]).
pub fn build_table(steps: &[StarStep<DD>], model: &Model, rmin: f64, rmax: f64) -> JTable<DD> {
    let mut t = JTable { pattern: vec![], x0: vec![], x1: vec![], c0: vec![], c1: vec![], linear_defect: 0.0, p_start: vec![] };
    let mut pat: Vec<(usize, usize)> = vec![];
    let mut tmp: Vec<(JM, JM, JM)> = vec![];
    for st in steps {
        if st.hi < rmin || st.lo > rmax {
            continue;
        }
        let mut b = st.lo.max(rmin);
        let top = st.hi.min(rmax);
        while b < top {
            let mut co = st.coef;
            let delta = DD::from_f64(b) - DD::from_f64(st.x0);
            for f in 0..4 {
                taylor_shift(&mut co[f][..=NS], delta, NS + 1);
            }
            let (r, alpha, a, phi0, u) = bg_series(&co, b);
            let ja = j_series(&r, &alpha, &a, &phi0, &u, &pars(model, DD::ZERO));
            let jb = j_series(&r, &alpha, &a, &phi0, &u, &pars(model, DD::from_f64(1.0)));
            let jc = j_series(&r, &alpha, &a, &phi0, &u, &pars(model, DD::from_f64(1.5)));
            let mut sers: Vec<[DD; NB]> = vec![];
            for i in 0..4 {
                for j in 0..4 {
                    if ja[i][j].n > 0 || jb[i][j].n > 0 {
                        if !pat.contains(&(i, j)) {
                            pat.push((i, j));
                        }
                        let mut a0 = [DD::ZERO; NB];
                        let mut a1 = [DD::ZERO; NB];
                        for k in 0..NB {
                            a0[k] = ja[i][j].get(k);
                            a1[k] = jb[i][j].get(k) - ja[i][j].get(k);
                        }
                        sers.push(a0);
                        sers.push(a1);
                    }
                }
            }
            let hmax = top - b;
            let h = if crate::mutate::on("wide_table") { hmax } else { (WIDTH_FRACTION * width_of(&sers, hmax / WIDTH_FRACTION, TAIL)).min(hmax) };
            let x1 = if b + h >= top || h < 1e-14 * b { top } else { b + h };
            t.x0.push(b);
            t.x1.push(x1);
            t.p_start.push(co[0][0].to_f64());
            tmp.push((ja, jb, jc));
            b = x1;
        }
    }
    pat.sort();
    for (ja, jb, jc) in &tmp {
        let mut e0 = vec![];
        let mut e1 = vec![];
        for &(i, j) in &pat {
            let mut a0 = [DD::ZERO; NB];
            let mut a1 = [DD::ZERO; NB];
            for k in 0..NB {
                a0[k] = ja[i][j].get(k);
                a1[k] = jb[i][j].get(k) - ja[i][j].get(k);
            }
            let pred = a0[0] + DD::from_f64(2.25) * a1[0];
            let _ = pred;
            // J(s2 = 1.5) = J0 + 1.5 J1
            let pred = a0[0] + DD::from_f64(1.5) * a1[0];
            let got = jc[i][j].get(0);
            let rel = ((pred - got).to_f64() / got.to_f64().abs().max(1e-300)).abs();
            if got.to_f64() != 0.0 {
                t.linear_defect = t.linear_defect.max(rel);
            }
            if crate::mutate::on("jsign") && (i, j) == (2, 3) {
                for k in 0..NB {
                    a0[k] = -a0[k];
                    a1[k] = -a1[k];
                }
            }
            e0.push(a0);
            e1.push(a1);
        }
        t.c0.push(e0);
        t.c1.push(e1);
    }
    t.pattern = pat;
    t
}

impl JTable<DD> {
    pub fn to_lane<T: Num>(&self) -> JTable<T> {
        let cv = |v: &Vec<Vec<[DD; NB]>>| -> Vec<Vec<[T; NB]>> {
            v.iter()
                .map(|e| {
                    e.iter()
                        .map(|a| {
                            let mut o = [T::zero(); NB];
                            for k in 0..NB {
                                o[k] = T::from_dd(a[k]);
                            }
                            o
                        })
                        .collect()
                })
                .collect()
        };
        JTable { pattern: self.pattern.clone(), x0: self.x0.clone(), x1: self.x1.clone(), c0: cv(&self.c0), c1: cv(&self.c1), linear_defect: self.linear_defect, p_start: self.p_start.clone() }
    }
}

impl<T: Num> JetSource<T> for JTable<T> {
    fn pattern(&self) -> Vec<(usize, usize)> {
        self.pattern.clone()
    }
    fn jets(&self, x0: f64, s2: T, kk: usize, out: &mut Jets<T>) {
        let idx = match self.x0.binary_search_by(|v| v.partial_cmp(&x0).unwrap()) {
            Ok(i) => i,
            Err(i) => i.saturating_sub(1),
        };
        assert!(x0 >= self.x0[0] * (1.0 - 1e-15) && x0 <= *self.x1.last().unwrap() * (1.0 + 1e-15), "x0 = {x0} outside the table");
        let delta = if crate::mutate::on("shift_off") { T::zero() } else { T::from_f64(x0) - T::from_f64(self.x0[idx]) };
        for (e, en) in out.ent.iter_mut().enumerate() {
            let a0 = &self.c0[idx][e];
            let a1 = &self.c1[idx][e];
            let mut tmp = [T::zero(); NB];
            for k in 0..NB {
                tmp[k] = a0[k] + s2 * a1[k];
            }
            taylor_shift(&mut tmp, delta, kk);
            en.c[..kk].copy_from_slice(&tmp[..kk]);
            en.n = kk;
        }
    }
}

// ------------------------------------------------------------------------------------------------------------
// the steps of the layout's solution
// ------------------------------------------------------------------------------------------------------------
use multishoot::lane::DV;
use multishoot::{Fl, Work};

/// The Taylor steps of every segment of the layout at z, as series of (P, u, m, Z) in double-double, ordered by position.
pub fn steps_layout<E: multishoot::Equations<1>>(sh: &multishoot::Shooting<E, 1>, z: &[Vec<f64>], g: &[DV]) -> Vec<StarStep<DD>> {
    let kt = sh.kt.min(NS);
    assert_eq!(sh.kt, NS, "the table takes the layout's order {NS}");
    let mut out: Vec<StarStep<DD>> = vec![];
    for b in 1..sh.k() {
        let pts = &sh.segs[b - 1];
        let src = sh.seg_source(b);
        let mut st: Vec<DD> = (0..E::NF).map(|f| crate::bvp::zget(z, sh.off_node(src) + f * 8)).collect();
        for w2 in pts.windows(2) {
            let (x0, x1) = (w2[0], w2[1]);
            let mut w = Work::<DV, 1>::new(E::NF, E::NAUX, kt);
            for f in 0..E::NF {
                w.c[f][0] = Fl::splat(st[f]);
            }
            sh.eq.taylor_coeffs::<DV>(x0, g, &mut w);
            let mut coef = [[DD::ZERO; KMAX + 1]; 4];
            let mut full = vec![[DD::ZERO; KMAX + 1]; E::NF];
            for f in 0..E::NF {
                for k in 0..=kt {
                    full[f][k] = w.c[f][k].get(0);
                }
            }
            for f in 0..4 {
                coef[f] = full[f];
            }
            let hh = DD::from_f64(x1) - DD::from_f64(x0);
            for f in 0..E::NF {
                let mut s = full[f][kt];
                for k in (0..kt).rev() {
                    s = s * hh + full[f][k];
                }
                st[f] = s;
            }
            out.push(StarStep { lo: x0.min(x1), hi: x0.max(x1), x0, coef });
        }
    }
    out.sort_by(|a, b| a.lo.partial_cmp(&b.lo).unwrap());
    out
}
