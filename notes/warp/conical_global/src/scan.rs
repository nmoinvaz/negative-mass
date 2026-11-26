//! The verdict of one (F_inf, C): the inward walk from eta = +200, and the bisection in C between verdict classes.
use crate::model::{conical, end_state, Mutation};
use crate::poly::{abs_dd, PolyN};
use crate::walk::{walk, Est, Outcome, Walk, WalkOpts, S3};
use multishoot::lane::DV;
use multishoot::{Fl, DD};

pub const ETA_MIN: f64 = 1e-6;
pub const LEND: f64 = crate::tables::LEND;

#[derive(Clone, Debug)]
pub struct Through {
    /// F_inf^-, C^- of the 2-d fit of the series at eta = -200 to the integrated state
    pub fi_minus: f64,
    pub c_minus: f64,
    /// the integrated state against the series (F_inf, C^-) at -200, F_inf kept: max(|dF|, |d(eta F')|)
    pub mismatch: f64,
    pub state_end: S3,
}

#[derive(Clone, Debug)]
pub struct Verdict {
    pub fi: f64,
    pub c: DD,
    /// "SINGULAR", "FLOOR", "THROUGH", "STEPLIMIT", "NONFINITE"
    pub class: String,
    /// kind of a singular verdict; for FLOOR "log" or "bounded"
    pub kind: String,
    pub xstar: Option<f64>,
    pub bar: Option<f64>,
    pub est: Option<Est>,
    pub x_last: f64,
    pub s_last: S3,
    /// F, F', F'' at eta = 1e-2, 1e-4, 1e-6 (FLOOR only)
    pub floor: Vec<S3>,
    pub through: Option<Through>,
    pub global: bool,
    pub steps: usize,
    /// max |N|/sum|terms| over the steps until the locator's R falls below 1e-3
    pub n_rel_max: f64,
    /// theoretical exponent 2/(x*^2 + 4) for kind F->0
    pub alpha_pred: Option<f64>,
}

pub fn lane0(f: &Fl<DV, 1>) -> DD {
    f.get(0)
}

pub fn start_state(fi: f64, c: DD, m: &Mutation) -> S3 {
    let s = end_state::<DV>(Fl::splat(DD::from_f64(fi)), Fl::splat(c), 1.0, m);
    [lane0(&s[0]), lane0(&s[1]), lane0(&s[2])]
}

/// The kind label of a singular verdict, from the last state (x_last) and the locator's estimate (labels, not criteria):
/// `F->0` if F vanishes at x* like Delta^alpha (alpha from the fit, 0 < alpha, and R F'/F = alpha to 10%); `F->inf` if
/// |F| >= 1e3 max(1, |F_inf|); `a3->0` if the F''^3 coefficient's factor F^2 - x^2 - 1 is below 1e-3 of F^2 + x^2 + 1;
/// `fold` if dN/dF'' is below 1e-6 of the sum of its terms; otherwise `other`.
pub fn conical_kind(p: &PolyN, s: &S3, x: f64, fi: f64, est: Option<&Est>) -> String {
    let (f, f1, f2) = (s[0].hi, s[1].hi, s[2].hi);
    if let Some(e) = est {
        let a = e.alpha();
        if a > 0.0 && f != 0.0 && (e.r * f1 / f - a).abs() <= 0.1 * a {
            return "F->0".into();
        }
    }
    if f.abs() >= 1e3 * fi.abs().max(1.0) {
        return "F->inf".into();
    }
    let c = p.coef_f2(s[0], s[1], DD::from_f64(x));
    let a = p.d_f2(s[0], s[1], s[2], DD::from_f64(x));
    let den = (c[0] * s[2] * s[2]).mul_f(3.0).hi.abs() + (c[1] * s[2]).mul_f(2.0).hi.abs() + c[2].hi.abs();
    let arel = if den > 0.0 { a.hi.abs() / den } else { 1.0 };
    let r3 = (f * f - x * x - 1.0).abs() / (f * f + x * x + 1.0);
    if r3 <= 1e-3 {
        return "a3->0".into();
    }
    if arel <= 1e-6 {
        return "fold".into();
    }
    let _ = f2;
    "other".into()
}

fn n_rel_until(w: &Walk) -> f64 {
    let mut m = 0.0f64;
    for st in &w.steps {
        if let Some(e) = &st.est {
            if e.r < 1e-3 {
                break;
            }
        }
        m = m.max(st.n_rel);
    }
    m
}

/// the series fit at eta = -200: (F_inf, C) from (F, F'), Newton with a difference Jacobian in double-double
pub fn fit_minus(s: &S3, fi0: f64, c0: DD, m: &Mutation) -> Option<(DD, DD)> {
    fit_end(s, fi0, c0, -1.0, m)
}

/// (F_inf, C) of the series at eta = sign 200 from (F, F')
pub fn fit_end(s: &S3, fi0: f64, c0: DD, sign: f64, m: &Mutation) -> Option<(DD, DD)> {
    let eval = |fi: DD, c: DD| {
        let st = end_state::<DV>(Fl::splat(fi), Fl::splat(c), sign, m);
        [lane0(&st[0]), lane0(&st[1])]
    };
    let (mut fi, mut c) = (DD::from_f64(fi0), c0);
    for _ in 0..30 {
        let v = eval(fi, c);
        let r = [v[0] - s[0], v[1] - s[1]];
        if abs_dd(r[0]).hi.max(abs_dd(r[1]).hi) < 1e-28 {
            return Some((fi, c));
        }
        let hf = DD::from_f64(1e-9 * (1.0 + fi.hi.abs()));
        let hc = DD::from_f64(1e-9 * (1.0 + c.hi.abs()));
        let vf = eval(fi + hf, c);
        let vc = eval(fi, c + hc);
        let j = [[(vf[0] - v[0]) / hf, (vc[0] - v[0]) / hc], [(vf[1] - v[1]) / hf, (vc[1] - v[1]) / hc]];
        let det = j[0][0] * j[1][1] - j[0][1] * j[1][0];
        if det.hi == 0.0 || !det.hi.is_finite() {
            return None;
        }
        let dfi = (j[1][1] * r[0] - j[0][1] * r[1]) / det;
        let dc = (j[0][0] * r[1] - j[1][0] * r[0]) / det;
        fi = fi - dfi;
        c = c - dc;
    }
    None
}

/// C at eta = -200 with F_inf fixed, from F alone; the mismatch max(|dF|, |d(eta F')|) after
pub fn fit_c_fixed(s: &S3, fi: f64, c0: DD, m: &Mutation) -> (DD, f64) {
    let eval = |c: DD| {
        let st = end_state::<DV>(Fl::splat(DD::from_f64(fi)), Fl::splat(c), -1.0, m);
        [lane0(&st[0]), lane0(&st[1])]
    };
    let mut c = c0;
    for _ in 0..30 {
        let v = eval(c);
        let r = v[0] - s[0];
        if abs_dd(r).hi < 1e-30 {
            break;
        }
        let h = DD::from_f64(1e-6 * (1.0 + c.hi.abs()));
        let d = (eval(c + h)[0] - v[0]) / h;
        if d.hi == 0.0 {
            break;
        }
        c = c - r / d;
    }
    let v = eval(c);
    let dp = (v[1] - s[1]).hi.abs() * LEND;
    ((c), (v[0] - s[0]).hi.abs().max(dp))
}

/// the verdict of the inward walk for (F_inf, C)
pub fn verdict(fi: f64, c: DD, m: &Mutation, opts: &WalkOpts) -> Verdict {
    let eq = conical(*m);
    let p = eq.n.clone();
    let g = DD::from_f64(fi);
    let s0 = start_state(fi, c, m);
    let mut v = Verdict {
        fi,
        c,
        class: String::new(),
        kind: String::new(),
        xstar: None,
        bar: None,
        est: None,
        x_last: LEND,
        s_last: s0,
        floor: Vec::new(),
        through: None,
        global: false,
        steps: 0,
        n_rel_max: 0.0,
        alpha_pred: None,
    };
    let targets = [1e-2, 1e-4, ETA_MIN];
    let (mut x, mut s) = (LEND, s0);
    let mut nrel = 0.0f64;
    for &t in &targets {
        let w = walk(&eq, &p, s, x, t, g, opts);
        v.steps += w.steps.len();
        nrel = nrel.max(n_rel_until(&w));
        match &w.outcome {
            Outcome::End => {
                x = t;
                s = w.last().s;
                v.floor.push(s);
            }
            Outcome::Singular { xstar, bar, est, x_last } => {
                v.class = "SINGULAR".into();
                v.xstar = Some(*xstar);
                v.bar = Some(*bar);
                v.est = Some(*est);
                v.x_last = *x_last;
                v.s_last = w.last().s;
                v.kind = conical_kind(&p, &v.s_last, *x_last, fi, Some(est));
                if v.kind == "F->0" {
                    v.alpha_pred = Some(2.0 / (xstar * xstar + 4.0));
                }
                v.n_rel_max = nrel;
                return v;
            }
            Outcome::StepLimit => {
                v.class = "STEPLIMIT".into();
                v.x_last = w.last().x;
                v.s_last = w.last().s;
                v.n_rel_max = nrel;
                return v;
            }
            Outcome::NonFinite { x: xx } => {
                v.class = "NONFINITE".into();
                v.x_last = *xx;
                v.s_last = w.last().s;
                v.n_rel_max = nrel;
                return v;
            }
        }
    }
    v.n_rel_max = nrel;
    v.x_last = ETA_MIN;
    v.s_last = s;
    let bounded = v.floor[2][1].hi.abs() <= 2.0 * v.floor[0][1].hi.abs();
    v.class = "FLOOR".into();
    v.kind = if bounded { "bounded".into() } else { "log".into() };
    if bounded {
        let w = walk(&eq, &p, s, ETA_MIN, -LEND, g, opts);
        v.steps += w.steps.len();
        match &w.outcome {
            Outcome::End => {
                let se = w.last().s;
                if let Some((fim, cm)) = fit_minus(&se, fi, c, m) {
                    let (_, mism) = fit_c_fixed(&se, fi, cm, m);
                    v.global = (fim.hi - fi).abs() <= 1e-10 && mism <= 1e-10;
                    v.through = Some(Through { fi_minus: fim.hi, c_minus: cm.hi, mismatch: mism, state_end: se });
                } else {
                    v.through = Some(Through { fi_minus: f64::NAN, c_minus: f64::NAN, mismatch: f64::NAN, state_end: se });
                }
                v.class = "THROUGH".into();
            }
            Outcome::Singular { xstar, bar, est, x_last } => {
                v.class = "SINGULAR".into();
                v.kind = conical_kind(&p, &w.last().s, *x_last, fi, Some(est));
                v.xstar = Some(*xstar);
                v.bar = Some(*bar);
                v.est = Some(*est);
                v.x_last = *x_last;
                v.s_last = w.last().s;
            }
            Outcome::StepLimit => v.class = "STEPLIMIT".into(),
            Outcome::NonFinite { x: xx } => {
                v.class = "NONFINITE".into();
                v.x_last = *xx;
            }
        }
    }
    v
}

/// the class label compared by the bisection
pub fn class_label(v: &Verdict) -> String {
    match v.class.as_str() {
        "SINGULAR" => format!("S:{}", v.kind),
        "THROUGH" => "THROUGH".into(),
        "FLOOR" => format!("FLOOR:{}", v.kind),
        c => c.into(),
    }
}

/// What one side of the forward walk from eta = 0 ended in
#[derive(Clone, Debug)]
pub struct Side {
    pub class: String,
    pub kind: String,
    pub xstar: Option<f64>,
    pub bar: Option<f64>,
    /// F_inf, C of the series fitted at +-200 when the walk reached it
    pub fit: Option<(f64, f64)>,
    pub mismatch: f64,
    pub steps: usize,
}

/// The solution through eta = 0 with F(0) = f, F'(0) = f1 and F''(0) = -f1^2/f (the root of N(0) = 0 for f^2 != 1),
/// walked to +200 and to -200.
pub fn forward(f: f64, f1: f64, m: &Mutation, opts: &WalkOpts) -> (Side, Side, f64) {
    let eq = conical(*m);
    let p = eq.n.clone();
    let (fd, f1d) = (DD::from_f64(f), DD::from_f64(f1));
    let s0: S3 = [fd, f1d, -(f1d * f1d) / fd];
    let (_, _) = (0, 0);
    let n0 = crate::poly::abs_dd(p.eval(s0[0], s0[1], s0[2], DD::ZERO).0).hi;
    let side = |to: f64| -> Side {
        let sign = if to > 0.0 { 1.0 } else { -1.0 };
        let w = walk(&eq, &p, s0, 0.0, to, fd, opts);
        let mut sd = Side { class: String::new(), kind: String::new(), xstar: None, bar: None, fit: None, mismatch: f64::NAN, steps: w.steps.len() };
        match &w.outcome {
            Outcome::End => {
                let se = w.last().s;
                sd.class = "LIMIT".into();
                if let Some((fi, c)) = fit_end(&se, f, DD::ZERO, sign, m) {
                    sd.fit = Some((fi.hi, c.hi));
                    let (_, mism) = fit_c_fixed_sign(&se, fi.hi, c, sign, m);
                    sd.mismatch = mism;
                }
            }
            Outcome::Singular { xstar, bar, est, x_last } => {
                sd.class = "SINGULAR".into();
                sd.xstar = Some(*xstar);
                sd.bar = Some(*bar);
                sd.kind = conical_kind(&p, &w.last().s, *x_last, f, Some(est));
            }
            Outcome::StepLimit => sd.class = "STEPLIMIT".into(),
            Outcome::NonFinite { .. } => sd.class = "NONFINITE".into(),
        }
        sd
    };
    (side(LEND), side(-LEND), n0)
}

fn fit_c_fixed_sign(s: &S3, fi: f64, c0: DD, sign: f64, m: &Mutation) -> (DD, f64) {
    let eval = |c: DD| {
        let st = end_state::<DV>(Fl::splat(DD::from_f64(fi)), Fl::splat(c), sign, m);
        [lane0(&st[0]), lane0(&st[1])]
    };
    let v = eval(c0);
    let dp = (v[1] - s[1]).hi.abs() * LEND;
    (c0, (v[0] - s[0]).hi.abs().max(dp))
}
