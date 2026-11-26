//! A walk by Taylor steps along x for the implicit cubic, with the singularity locator.
//!
//! Each step takes the Taylor coefficients about the current point (the recursion of `poly`), the step from multishoot's
//! planner rule (`taylor::step_size`: the last two terms below `tol`), and the state by Horner's rule. The locator reads
//! the coefficients of each field: the ratios r_k = |c_(k-1)/c_k| of a series with its nearest singularity at distance R
//! behave as R (1 + gamma/k + delta/k^2 + ...), so a least-squares fit in 1/k over the last ten orders gives R, with gamma
//! the exponent of the singular part (F ~ (x - x*)^alpha has gamma = alpha for the field F... see `Est::alpha`). The sign
//! of c_k/c_(k-1) gives the side. The singular point is x* = x_n + dir R_n; the walk stops when the estimate is real, ahead
//! and R <= `rstop` over five consecutive steps, and the spread of those five estimates is the error bar.
use crate::poly::{abs_dd, PolyN};
use multishoot::lane::DV;
use multishoot::taylor::{all_finite, step_size};
use multishoot::{horner, Equations, Fl, Plan, Work, DD};

pub type S3 = [DD; 3];

#[derive(Clone, Copy, Debug)]
pub struct Est {
    pub r: f64,
    pub gamma: f64,
    pub resid: f64,
    pub field: usize,
    /// +1 if the singularity is at x + R, -1 at x - R
    pub side: i32,
}

impl Est {
    /// the exponent alpha of F ~ (x - x*)^alpha implied by the field's gamma: field f behaves as (x - x*)^(alpha - f)
    pub fn alpha(&self) -> f64 {
        self.gamma + self.field as f64 - 1.0
    }
}

#[derive(Clone, Debug)]
pub enum Outcome {
    /// reached the end of the interval
    End,
    /// a real singularity ahead: x* and the spread of the last five estimates
    Singular { xstar: f64, bar: f64, est: Est, x_last: f64 },
    StepLimit,
    NonFinite { x: f64 },
}

#[derive(Clone, Debug)]
pub struct Step {
    pub x: f64,
    pub s: S3,
    /// |N| / sum |terms of N| at the state
    pub n_rel: f64,
    pub est: Option<Est>,
}

#[derive(Clone, Debug)]
pub struct Walk {
    pub steps: Vec<Step>,
    pub outcome: Outcome,
}

impl Walk {
    pub fn last(&self) -> &Step {
        self.steps.last().unwrap()
    }
}

pub struct WalkOpts {
    pub plan: Plan,
    pub rstop: f64,
    pub max_steps: usize,
    /// control `plant-locator`: the side test inverted
    pub flip_side: bool,
}

impl WalkOpts {
    pub fn standard() -> Self {
        WalkOpts { plan: Plan { kt: 30, tol: 1e-32, hmax: 10.0 }, rstop: 1e-6, max_steps: 4000, flip_side: false }
    }
}

fn solve3(mut a: [[f64; 3]; 3], mut b: [f64; 3]) -> Option<[f64; 3]> {
    for c in 0..3 {
        let mut piv = c;
        for r in c + 1..3 {
            if a[r][c].abs() > a[piv][c].abs() {
                piv = r;
            }
        }
        if a[piv][c] == 0.0 {
            return None;
        }
        a.swap(c, piv);
        b.swap(c, piv);
        for r in c + 1..3 {
            let f = a[r][c] / a[c][c];
            for k in c..3 {
                a[r][k] -= f * a[c][k];
            }
            b[r] -= f * b[c];
        }
    }
    let mut x = [0.0; 3];
    for c in (0..3).rev() {
        let mut s = b[c];
        for k in c + 1..3 {
            s -= a[c][k] * x[k];
        }
        x[c] = s / a[c][c];
    }
    Some(x)
}

/// the radius estimate of one field's coefficients c[0..=kt], from the orders kt-9..=kt
pub fn fit_radius(c: &[f64], kt: usize, field: usize) -> Option<Est> {
    let k0 = kt.saturating_sub(9).max(2);
    let mut ks = Vec::new();
    let mut rs = Vec::new();
    let mut sign = 0i32;
    for k in k0..=kt {
        let q = c[k] / c[k - 1];
        if !q.is_finite() || q == 0.0 {
            return None;
        }
        let s = if q > 0.0 { 1 } else { -1 };
        if sign == 0 {
            sign = s;
        } else if s != sign {
            return None;
        }
        ks.push(k as f64);
        rs.push(1.0 / q.abs());
    }
    // least squares r = u0 + u1/k + u2/k^2
    let mut a = [[0.0f64; 3]; 3];
    let mut b = [0.0f64; 3];
    for (&k, &r) in ks.iter().zip(&rs) {
        let v = [1.0, 1.0 / k, 1.0 / (k * k)];
        for i in 0..3 {
            for j in 0..3 {
                a[i][j] += v[i] * v[j];
            }
            b[i] += v[i] * r;
        }
    }
    let u = solve3(a, b)?;
    if !(u[0].is_finite() && u[0] > 0.0) {
        return None;
    }
    let resid = ks
        .iter()
        .zip(&rs)
        .map(|(&k, &r)| (r - (u[0] + u[1] / k + u[2] / (k * k))).abs() / r)
        .fold(0.0, f64::max);
    Some(Est { r: u[0], gamma: u[1] / u[0], resid, field, side: sign })
}

/// the best radius estimate over the three fields (smallest fit residual, below 1e-3)
pub fn estimate(w: &Work<DV, 1>) -> Option<Est> {
    let kt = w.kt;
    let mut best: Option<Est> = None;
    for f in 0..3 {
        let c: Vec<f64> = (0..=kt).map(|k| w.c[f][k].get(0).hi).collect();
        if let Some(e) = fit_radius(&c, kt, f) {
            if e.resid < 1e-3 && best.is_none_or(|b| e.resid < b.resid) {
                best = Some(e);
            }
        }
    }
    best
}

fn n_rel(p: &PolyN, s: &S3, x: f64) -> f64 {
    let (v, sc) = p.eval(s[0], s[1], s[2], DD::from_f64(x));
    if sc.hi == 0.0 {
        0.0
    } else {
        (abs_dd(v) / sc).hi
    }
}

/// Walk from x0 with the state s0 to x1. `g` is the model's global. The first element of `steps` is the start.
pub fn walk<E: Equations<1>>(eq: &E, p: &PolyN, s0: S3, x0: f64, x1: f64, g: DD, opts: &WalkOpts) -> Walk {
    let kt = opts.plan.kt;
    let mut w = Work::<DV, 1>::new(E::NF, E::NAUX, kt);
    let gl = [DV::splat(g)];
    let dir: f64 = if x1 > x0 { 1.0 } else { -1.0 };
    let (mut x, mut s) = (x0, s0);
    let mut steps = vec![Step { x, s, n_rel: n_rel(p, &s, x), est: None }];
    let mut hist: Vec<(f64, f64)> = Vec::new(); // (x*, R)
    while (x1 - x) * dir > 0.0 {
        if steps.len() > opts.max_steps {
            return Walk { steps, outcome: Outcome::StepLimit };
        }
        for (f, sf) in s.iter().enumerate() {
            w.c[f][0] = Fl::splat(*sf);
        }
        eq.taylor_coeffs(x, &gl, &mut w);
        if w.c.iter().any(|cf| cf[..=kt].iter().any(|c| !all_finite(c))) {
            return Walk { steps, outcome: Outcome::NonFinite { x } };
        }
        let est = estimate(&w);
        let ahead = |e: &Est| (e.side as f64 == dir) != opts.flip_side;
        match est {
            Some(e) if ahead(&e) => hist.push((x + dir * e.r, e.r)),
            _ => hist.clear(),
        }
        steps.last_mut().unwrap().est = est;
        if hist.len() >= 5 {
            let last5 = &hist[hist.len() - 5..];
            let (xs, r) = *last5.last().unwrap();
            if r <= opts.rstop {
                let lo = last5.iter().map(|h| h.0).fold(f64::INFINITY, f64::min);
                let hi = last5.iter().map(|h| h.0).fold(f64::NEG_INFINITY, f64::max);
                return Walk {
                    steps,
                    outcome: Outcome::Singular { xstar: xs, bar: hi - lo, est: est.unwrap(), x_last: x },
                };
            }
        }
        let mut h = step_size(&w, &opts.plan);
        if !(h.is_finite() && h > 0.0) {
            return Walk { steps, outcome: Outcome::NonFinite { x } };
        }
        let rest = (x1 - x) * dir;
        if h >= rest {
            h = rest;
        } else if h > 0.6 * rest {
            h = 0.5 * rest;
        }
        let xn = if h == rest { x1 } else { x + dir * h };
        if xn == x {
            return Walk { steps, outcome: Outcome::NonFinite { x } };
        }
        let hd = DD::from_f64(xn) - DD::from_f64(x);
        let ns: Vec<Fl<DV, 1>> = (0..3).map(|f| horner(&w.c[f][..=kt], hd)).collect();
        if ns.iter().any(|f| !all_finite(f)) {
            return Walk { steps, outcome: Outcome::NonFinite { x } };
        }
        s = [ns[0].get(0), ns[1].get(0), ns[2].get(0)];
        x = xn;
        steps.push(Step { x, s, n_rel: n_rel(p, &s, x), est: None });
    }
    Walk { steps, outcome: Outcome::End }
}
