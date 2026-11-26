//! J(r; s^2) = J0(r) + s^2 J1(r) as Taylor tables along the star, and the vacuum jets of the exterior.
//!
//! The matrix of polar_matrix_generated.py is linear in s^2 (checked at the table's build against the matrix at
//! s^2 = 2.25), so one table per l, built once from the star's steps, serves every gamma: a step of the plane
//! integrator takes the jets of each entry about its start by shifting the table's series (J0 + s^2 J1) from the
//! start of the table step that holds it.
use crate::jgen::{j_series, Pars};
use crate::num::Num;
use crate::plane::{Jets, JetSource};
use crate::ser::*;
use crate::star::{Star, NS};
use multidouble::dd::DD;

/// coefficients kept per entry
pub const NB: usize = NS; // 0..NB-1

#[derive(Clone, Copy, Debug)]
pub struct Model {
    pub omega: f64,
    pub epsilon: f64,
    pub lambda: f64,
    pub mu: f64,
}

pub fn pars<T: Num>(m: &Model, l: usize, s: T) -> Pars<T> {
    Pars {
        omega: T::from_f64(m.omega),
        epsilon: T::from_f64(m.epsilon),
        s,
        ll: T::from_f64((l * (l + 1)) as f64),
        lambda: T::from_f64(m.lambda),
        mu: T::from_f64(m.mu),
        pi: T::pi(),
    }
}

pub struct JTable<T: Num> {
    pub l: usize,
    pub pattern: Vec<(usize, usize)>,
    /// step starts and ends
    pub x0: Vec<f64>,
    pub x1: Vec<f64>,
    /// [step][entry][k] for J0 and J1
    pub c0: Vec<Vec<[T; NB]>>,
    pub c1: Vec<Vec<[T; NB]>>,
    /// largest |J(s^2 = 2.25) - (J0 + 2.25 J1)| over entries and steps, relative to |J|, at k = 0
    pub linear_defect: f64,
}

/// The series of (alpha, a, phi0, u) about the start of a star step.
pub fn bg_series<T: Num>(co: &[[T; KMAX + 1]; 4], x0: f64) -> (Ser<T>, Ser<T>, Ser<T>, Ser<T>, Ser<T>) {
    let r = Ser::var(T::from_f64(x0));
    let p = Ser::from_slice(&co[0][..=NS]);
    let u = Ser::from_slice(&co[1][..=NS]);
    let m = Ser::from_slice(&co[2][..=NS]);
    let z = Ser::from_slice(&co[3][..=NS]);
    let rinv = recip(&r);
    let alpha = sqrt(&recip(&z));
    let b = adds(&scl(&mul(&m, &rinv), -(T::one() + T::one())), T::one());
    let a = sqrt(&recip(&b));
    let phi0 = sqrt(&p);
    (r, alpha, a, phi0, u)
}

/// Largest width h (up to hmax) at which the series of every entry still converges to `tail` relative to the sum of
/// its terms: |c_k| h^k / sum_j |c_j| h^j <= tail for the last two coefficients kept, by bisection on h.
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

/// Table steps of J along the star. A star step is cut where J's own series stops converging: J carries 1/r and 1/r^2
/// terms the star's (regular) series lack, so a table step is at most `WIDTH_FRACTION` of the width at which J's series
/// has converged to `TAIL`; the plane integrator's shift plus its own step then stay inside the radius.
pub const TAIL: f64 = 1e-25;
pub const WIDTH_FRACTION: f64 = 0.3;

pub fn build_table(star: &Star<DD>, l: usize, model: &Model, rmin: f64, rmax: f64) -> JTable<DD> {
    let mut t = JTable { l, pattern: vec![], x0: vec![], x1: vec![], c0: vec![], c1: vec![], linear_defect: 0.0 };
    let mut pat: Vec<(usize, usize)> = vec![];
    let mut tmp: Vec<([[Ser<DD>; 6]; 6], [[Ser<DD>; 6]; 6], [[Ser<DD>; 6]; 6])> = vec![];
    for st in &star.steps {
        if st.x1 < rmin || st.x0 > rmax {
            continue;
        }
        let mut b = st.x0.max(rmin * 0.999);
        while b < st.x1 {
            // the star's series re-expanded about b
            let mut co = st.coef;
            let delta = DD::from_f64(b) - DD::from_f64(st.x0);
            for f in 0..4 {
                taylor_shift(&mut co[f][..=NS], delta, NS + 1);
            }
            let (r, alpha, a, phi0, u) = bg_series(&co, b);
            let ja = j_series(&r, &alpha, &a, &phi0, &u, &pars(model, l, DD::ZERO));
            let jb = j_series(&r, &alpha, &a, &phi0, &u, &pars(model, l, DD::from_f64(1.0)));
            let jc = j_series(&r, &alpha, &a, &phi0, &u, &pars(model, l, DD::from_f64(1.5)));
            let mut sers: Vec<[DD; NB]> = vec![];
            for i in 0..6 {
                for j in 0..6 {
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
            let hmax = st.x1 - b;
            let h = if crate::mutate::on("wide_table") { hmax } else { (WIDTH_FRACTION * width_of(&sers, hmax / WIDTH_FRACTION, TAIL)).min(hmax) };
            let x1 = if b + h >= st.x1 || h < 1e-14 * b { st.x1 } else { b + h };
            t.x0.push(b);
            t.x1.push(x1);
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
            // linearity in s^2: J(s = 1.5) = J0 + 2.25 J1
            let pred = a0[0] + DD::from_f64(2.25) * a1[0];
            let got = jc[i][j].get(0);
            let rel = ((pred - got).to_f64() / got.to_f64().abs().max(1e-300)).abs();
            if got.to_f64() != 0.0 {
                t.linear_defect = t.linear_defect.max(rel);
            }
            if crate::mutate::on("jsign") && (i, j) == (3, 2) {
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
        JTable {
            l: self.l,
            pattern: self.pattern.clone(),
            x0: self.x0.clone(),
            x1: self.x1.clone(),
            c0: cv(&self.c0),
            c1: cv(&self.c1),
            linear_defect: self.linear_defect,
        }
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
        assert!(x0 >= self.x0[0] && x0 <= *self.x1.last().unwrap() + 1e-12, "x0 = {x0} outside the table");
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

/// Jets of J in the vacuum of mass M (alpha = sqrt f, a = 1/sqrt f, phi0 = u = 0), restricted to the rows and columns
/// `lo..6` (the scalar block when lo = 2), at the fixed s.
pub struct VacSource<T: Num> {
    pub mass: T,
    pub model: Model,
    pub l: usize,
    pub s: T,
    pub lo: usize,
}

impl<T: Num> VacSource<T> {
    pub fn full(&self, x0: f64) -> [[Ser<T>; 6]; 6] {
        let r = Ser::var(T::from_f64(x0));
        let f = adds(&scl(&recip(&r), -(self.mass + self.mass)), T::one());
        let alpha = sqrt(&f);
        let a = recip(&alpha);
        let z = Ser::zero();
        j_series(&r, &alpha, &a, &z, &z, &pars(&self.model, self.l, self.s))
    }
}

impl<T: Num> JetSource<T> for VacSource<T> {
    fn pattern(&self) -> Vec<(usize, usize)> {
        // the scalar block's rows and columns; the metric rows must not feed on the scalar columns (checked in jets)
        let mut v = vec![];
        for i in self.lo..6 {
            for j in self.lo..6 {
                v.push((i - self.lo, j - self.lo));
            }
        }
        v
    }
    fn jets(&self, x0: f64, _s2: T, kk: usize, out: &mut Jets<T>) {
        let j = self.full(x0);
        for en in out.ent.iter_mut() {
            let s = &j[en.i + self.lo][en.j + self.lo];
            for k in 0..kk {
                en.c[k] = s.get(k);
            }
            en.n = if s.n == 0 { 0 } else { kk.min(s.n) };
        }
    }
}
