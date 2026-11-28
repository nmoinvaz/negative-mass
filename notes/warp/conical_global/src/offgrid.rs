//! offgrid: the verdict of one (F_inf, C) of notes/warp/conical_global.md section 3 on the precision-generic walker of
//! `hp.rs` (double-, quad- and octo-double), for start radii up to 819200 (issue #183, notes/warp/conical_offgrid.md).
//!
//! Differences from `scan.rs` (#167):
//! - the walker is `hp_walk` with the step variable s = (x - x0)/l (l a power of two), so the coefficients do not
//!   underflow at eta = 1e6, and with the singularity locator of `walk.rs` (`HpCfg::rstop`);
//! - the start state is the series of `tables_hp.rs` (coefficients as exact sums of 8 f64 parts, `SERIES_HP`) evaluated in
//!   the lane's precision at the first radius L = 200 * 4^k, k = 0..6, whose residual |N| / sum|terms| is at most 1e-20;
//! - the kind of a singularity, the floor test and the matching at eta = -200 are `scan.rs`'s, applied to the lane's
//!   state rounded to double-double.
use crate::hp::{hi, hp_walk, sc_of, EvalS, HpCfg, HpOut, PolyQ, Stop};
use crate::model::{conical_poly, Mutation};
use crate::scan::{conical_kind, fit_c_fixed, fit_minus};
use crate::tables_hp::{LN_L, SERIES_HP};
use crate::walk::Est;
use multishoot::multidouble::precision::Precision;
use multishoot::multidouble::simd::Native;
use multishoot::{Fl, DD};

type Sc<P> = <P as Precision<Native>>::Scalar;

/// the start radii: the smallest with residual <= 1e-20 is the start
pub const LS: [f64; 7] = [200.0, 800.0, 3200.0, 12800.0, 51200.0, 204800.0, 819200.0];
pub const ETA_MIN: f64 = 1e-6;
pub const LEND: f64 = 200.0;

pub fn div_s<P: Precision<Native>>(a: Sc<P>, b: Sc<P>) -> Sc<P> {
    let r = Fl([P::splat(a)]) / Fl([P::splat(b)]);
    P::field_get(&r, 0)
}

fn powers<P: Precision<Native>>(x: Sc<P>, n: usize) -> Vec<Sc<P>> {
    let mut v = vec![P::ONE];
    for i in 0..n {
        let t = v[i] * x;
        v.push(t);
    }
    v
}

/// (F, F', F'') of the series with logarithms at +L = LS[lk], from (F_inf, C), in the lane's precision.
/// F = F_inf + sum a_(n,j) nu^n T^j, T = ln L, nu = 1/L; G = F, p = theta G, q = theta^2 G, F' = p nu, F'' = (q - p) nu^2.
pub fn series_state<P: Precision<Native>>(fi: f64, c: f64, lk: usize, m: &Mutation) -> [Sc<P>; 3] {
    let fi_s = P::new(&[fi]);
    let c_s = P::new(&[c]);
    let mx = SERIES_HP.iter().flat_map(|t| t.mons.iter().map(|x| x.0 as usize)).max().unwrap();
    let ml = SERIES_HP.iter().flat_map(|t| t.mons.iter().map(|x| x.1 as usize)).max().unwrap();
    let ms = SERIES_HP.iter().map(|t| t.s as usize).max().unwrap();
    let nm = SERIES_HP.iter().map(|t| t.n as usize).max().unwrap();
    let jm = SERIES_HP.iter().map(|t| t.j as usize).max().unwrap();
    let (fp, cp) = (powers::<P>(fi_s, mx), powers::<P>(c_s, ml));
    let ip = powers::<P>(div_s::<P>(P::ONE, fi_s), ms);
    let nu = div_s::<P>(P::ONE, P::new(&[LS[lk]]));
    let tl = P::new(&LN_L[lk][..P::PARTS]);
    let (nup, tp) = (powers::<P>(nu, nm), powers::<P>(tl, jm));
    let tj = |e: i32| if e < 0 { P::ZERO } else { tp[e as usize] };
    let k = |x: f64| P::new(&[x]);
    let (mut g, mut p, mut q) = (fi_s, P::ZERO, P::ZERO);
    for t in SERIES_HP.iter() {
        let mut a = P::ZERO;
        for &(mm, l, co) in t.mons {
            a = a + P::new(&co[..P::PARTS]) * fp[mm as usize] * cp[l as usize];
        }
        a = a * ip[t.s as usize];
        if m.series_coeff && t.n == 4 && t.j == 1 {
            a = a * k(1.01);
        }
        let (n, j) = (t.n as f64, t.j as f64);
        let jj = t.j as i32;
        let nun = nup[t.n as usize];
        let kg = nun * tj(jj);
        let kp = nun * (k(-n) * tj(jj) + k(j) * tj(jj - 1));
        let kq = nun * (k(n * n) * tj(jj) - k(2.0 * n * j) * tj(jj - 1) + k(j * (j - 1.0)) * tj(jj - 2));
        g = g + a * kg;
        p = p + a * kp;
        q = q + a * kq;
    }
    [g, p * nu, (q - p) * nu * nu]
}

/// |N| / sum |terms| of a state at x
pub fn resid_of<P: Precision<Native>>(ev: &EvalS<P>, s: &[Sc<P>; 3], x: f64) -> f64 {
    let (v, sm, _) = ev.eval(s, P::new(&[x]));
    if hi::<P>(sm) == 0.0 {
        0.0
    } else {
        let a = if hi::<P>(v) < 0.0 { -v } else { v };
        hi::<P>(a) / hi::<P>(sm)
    }
}

#[derive(Clone, Debug)]
pub struct OgOpts {
    pub cfg: HpCfg,
    /// extra start radii above the first one that passes (the second start of the start-spread check)
    pub lk_extra: usize,
    /// the residual bar of the start rule (control `start-bar`: 1e-8)
    pub resid_bar: f64,
}

impl OgOpts {
    pub fn new(parts: usize) -> Self {
        let (kt, tol) = match parts {
            2 => (30usize, 1e-32f64),
            4 => (60, 1e-64),
            _ => (120, 1e-124),
        };
        let mut cfg = HpCfg::new(kt, tol);
        cfg.rstop = Some(1e-12);
        cfg.tiny = 1e-15;
        // the 1/k fit of the locator: the last 10 of 30 orders (double-double), the last 33 of 60 (quad-double and above)
        cfg.fit_window = if parts == 2 { 10 } else { 33 };
        cfg.project = false;
        OgOpts { cfg, lk_extra: 0, resid_bar: 1e-20 }
    }
}

#[derive(Clone, Debug)]
pub struct VerdictP {
    pub fi: f64,
    pub c: f64,
    /// SINGULAR, FLOOR, THROUGH, NOSTART, STEPLIMIT, NONFINITE, TINYSTEP
    pub class: String,
    pub kind: String,
    pub xstar: Option<f64>,
    pub bar: Option<f64>,
    pub est: Option<Est>,
    pub l: f64,
    pub resid: f64,
    /// the residuals at LS[0..=chosen index]
    pub resids: Vec<f64>,
    pub x_last: f64,
    pub s_last: [DD; 3],
    /// the state at eta = 200 (when the start is above it)
    pub s200: Option<[DD; 3]>,
    /// F, F', F'' at eta = 1e-2, 1e-4, 1e-6 (FLOOR and THROUGH)
    pub floor: Vec<[DD; 3]>,
    pub through: Option<(f64, f64, f64)>,
    pub global: bool,
    pub steps: usize,
    pub n_rel_max: f64,
    pub alpha_pred: Option<f64>,
    /// TINYSTEP only: the last steps (x, h, median radius of F'', sign changes of the last ten ratios) and the reading
    /// "shrinking" (h and the radius fall by a factor >= 1.05 at each of the last 5 steps: the walk closes on a singularity)
    /// or "stalled" (they do not)
    pub tiny: Vec<(f64, f64, f64, usize)>,
    pub tiny_kind: String,
}

fn dd3<P: Precision<Native>>(s: &[Sc<P>; 3]) -> [DD; 3] {
    [P::to_dd(s[0]), P::to_dd(s[1]), P::to_dd(s[2])]
}

fn stop_name(s: &Stop) -> &'static str {
    match s {
        Stop::End => "END",
        Stop::NonFinite => "NONFINITE",
        Stop::TinyStep => "TINYSTEP",
        Stop::StepLimit => "STEPLIMIT",
        Stop::Singular => "SINGULAR",
    }
}

/// the first index k with residual <= the bar at L = LS[k] (and the residuals up to it)
pub fn pick_start<P: Precision<Native>>(fi: f64, c: f64, m: &Mutation, o: &OgOpts, ev: &EvalS<P>) -> (Option<usize>, Vec<f64>) {
    let mut resids = Vec::new();
    for (lk, &l) in LS.iter().enumerate() {
        let s = series_state::<P>(fi, c, lk, m);
        let r = resid_of::<P>(ev, &s, l);
        resids.push(r);
        if r <= o.resid_bar && r.is_finite() {
            return (Some(lk), resids);
        }
    }
    (None, resids)
}

/// the label compared by a bisection
pub fn label(v: &VerdictP) -> String {
    match v.class.as_str() {
        "SINGULAR" => format!("S:{}", v.kind),
        "FLOOR" => format!("FLOOR:{}", v.kind),
        c => c.to_string(),
    }
}

/// the walk from the chosen start to x1 (no checkpoints): the state in parts, the steps and how it stopped
pub fn walk_to<P: Precision<Native>>(fi: f64, c: f64, x1: f64, m: &Mutation, o: &OgOpts) -> Option<(f64, HpOut<P>)> {
    let poly = PolyQ::conical(m);
    let ev = EvalS::<P>::new(&poly);
    let (pick, _) = pick_start::<P>(fi, c, m, o, &ev);
    let lk = pick? + o.lk_extra;
    if lk >= LS.len() {
        return None;
    }
    let s0 = series_state::<P>(fi, c, lk, m);
    Some((LS[lk], hp_walk::<P>(&poly, s0, LS[lk], x1, &[], &o.cfg, false)))
}

/// "shrinking" if h and the radius fall by a factor >= 1.05 at each of the last five steps, else "stalled"
pub fn tiny_reading(r: &[(f64, f64, f64, usize)]) -> String {
    if r.len() < 6 {
        return "stalled".into();
    }
    let t = &r[r.len() - 6..];
    let shrink = t.windows(2).all(|w| w[0].1 >= 1.05 * w[1].1 && w[0].2 >= 1.05 * w[1].2);
    if shrink { "shrinking".into() } else { "stalled".into() }
}

/// the verdict of the inward walk for (F_inf, C), in precision P
pub fn verdict_p<P: Precision<Native>>(fi: f64, c: f64, m: &Mutation, o: &OgOpts) -> VerdictP {
    let poly = PolyQ::conical(m);
    let ev = EvalS::<P>::new(&poly);
    let pn = conical_poly(m);
    let mut v = VerdictP {
        fi,
        c,
        class: String::new(),
        kind: String::new(),
        xstar: None,
        bar: None,
        est: None,
        l: 0.0,
        resid: f64::NAN,
        resids: Vec::new(),
        x_last: 0.0,
        s_last: [DD::ZERO; 3],
        s200: None,
        floor: Vec::new(),
        through: None,
        global: false,
        steps: 0,
        n_rel_max: 0.0,
        alpha_pred: None,
        tiny: Vec::new(),
        tiny_kind: String::new(),
    };
    // the start rule: the first radius whose residual is within the bar
    let (pick, resids) = pick_start::<P>(fi, c, m, o, &ev);
    v.resids = resids;
    let lk = match pick {
        Some(k) if k + o.lk_extra < LS.len() => k + o.lk_extra,
        _ => {
            v.class = "NOSTART".into();
            v.l = *LS.last().unwrap();
            return v;
        }
    };
    let s0 = series_state::<P>(fi, c, lk, m);
    let l0 = LS[lk];
    v.l = l0;
    v.resid = resid_of::<P>(&ev, &s0, l0);
    v.x_last = l0;
    v.s_last = dd3::<P>(&s0);
    let mut cps: Vec<f64> = Vec::new();
    if l0 > LEND {
        cps.push(LEND);
    }
    cps.extend([1e-2, 1e-4]);
    let out: HpOut<P> = hp_walk::<P>(&poly, s0, l0, ETA_MIN, &cps, &o.cfg, false);
    v.steps = out.steps;
    v.n_rel_max = out.n_rel_max;
    for cp in &out.cps {
        if cp.x == LEND {
            v.s200 = Some(dd3::<P>(&cp.s));
        } else {
            v.floor.push(dd3::<P>(&cp.s));
        }
    }
    match out.stop {
        Stop::Singular => {
            let sg = out.sing.as_ref().unwrap();
            v.class = "SINGULAR".into();
            v.xstar = Some(sg.xstar);
            v.bar = Some(sg.bar);
            v.est = Some(sg.est);
            v.x_last = sg.x_last;
            v.s_last = sg.s_last;
            v.kind = conical_kind(&pn, &sg.s_last, sg.x_last, fi, Some(&sg.est));
            if v.kind == "F->0" {
                v.alpha_pred = Some(2.0 / (sg.xstar * sg.xstar + 4.0));
            }
        }
        Stop::End => {
            let se = dd3::<P>(&out.s_end);
            v.floor.push(se);
            v.x_last = ETA_MIN;
            v.s_last = se;
            let bounded = v.floor[2][1].hi.abs() <= 2.0 * v.floor[0][1].hi.abs();
            v.class = "FLOOR".into();
            v.kind = if bounded { "bounded".into() } else { "log".into() };
            if bounded {
                let o2: HpOut<P> = hp_walk::<P>(&poly, out.s_end, ETA_MIN, -LEND, &[], &o.cfg, false);
                v.steps += o2.steps;
                match o2.stop {
                    Stop::End => {
                        let sm = dd3::<P>(&o2.s_end);
                        let cd = DD::from_f64(c);
                        match fit_minus(&sm, fi, cd, m) {
                            Some((fim, cm)) => {
                                let (_, mism) = fit_c_fixed(&sm, fi, cm, m);
                                v.global = (fim.hi - fi).abs() <= 1e-10 && mism <= 1e-10;
                                v.through = Some((fim.hi, cm.hi, mism));
                            }
                            None => v.through = Some((f64::NAN, f64::NAN, f64::NAN)),
                        }
                        v.class = "THROUGH".into();
                    }
                    Stop::Singular => {
                        let sg = o2.sing.as_ref().unwrap();
                        v.class = "SINGULAR".into();
                        v.xstar = Some(sg.xstar);
                        v.bar = Some(sg.bar);
                        v.est = Some(sg.est);
                        v.x_last = sg.x_last;
                        v.s_last = sg.s_last;
                        v.kind = conical_kind(&pn, &sg.s_last, sg.x_last, fi, Some(&sg.est));
                    }
                    other => v.class = stop_name(&other).into(),
                }
            }
        }
        other => {
            v.class = stop_name(&other).into();
            v.x_last = out.x_end;
            v.s_last = dd3::<P>(&out.s_end);
            if other == Stop::TinyStep {
                v.tiny = out.recent.clone();
                v.tiny_kind = tiny_reading(&v.tiny);
            }
        }
    }
    let _ = sc_of::<P>;
    v
}

/// the deliberate errors of the controls (Rust side); false if the name is not one of them
pub fn apply_mutation(name: &str, m: &mut Mutation, o: &mut OgOpts) -> bool {
    match name {
        "eq-coeff" => m.eq_coeff = true,
        "series-coeff" => m.series_coeff = true,
        "locator-side" => o.cfg.flip_side = true,
        "recip-dd" => o.cfg.recip_dd = true,
        "start-bar" => o.resid_bar = 1e-8,
        "n-drift" => o.cfg.skip = Some(3),
        "halving" => {
            o.cfg.kt = 14;
            o.cfg.tol = 1e-12;
        }
        "halving-coarse" => {
            o.cfg.kt = 8;
            o.cfg.tol = 1e-6;
        }
        _ => return false,
    }
    true
}

impl VerdictP {
    pub fn to_json(&self) -> String {
        use crate::json::{dd, num, opt, s3};
        let through = match &self.through {
            None => "null".to_string(),
            Some((fi, c, mism)) => format!("{{\"fi_minus\": {}, \"c_minus\": {}, \"mismatch\": {}}}", num(*fi), num(*c), num(*mism)),
        };
        let floor: Vec<String> = self.floor.iter().map(s3).collect();
        let s200 = match &self.s200 {
            None => "null".to_string(),
            Some(s) => s3(s),
        };
        let resids: Vec<String> = self.resids.iter().map(|r| num(*r)).collect();
        format!(
            "{{\"fi\": {}, \"c\": {}, \"class\": \"{}\", \"kind\": \"{}\", \"xstar\": {}, \"bar\": {}, \"est\": {}, \"l\": {}, \
             \"resid\": {}, \"resids\": [{}], \"x_last\": {}, \"s_last\": {}, \"s200\": {}, \"floor\": [{}], \"through\": {}, \
             \"global\": {}, \"steps\": {}, \"n_rel_max\": {}, \"alpha_pred\": {}, \"tiny_kind\": \"{}\", \"tiny\": [{}]}}",
            num(self.fi),
            dd(DD::from_f64(self.c)),
            self.class,
            self.kind,
            opt(self.xstar),
            opt(self.bar),
            crate::json::est(&self.est),
            num(self.l),
            num(self.resid),
            resids.join(", "),
            num(self.x_last),
            s3(&self.s_last),
            s200,
            floor.join(", "),
            through,
            self.global,
            self.steps,
            num(self.n_rel_max),
            opt(self.alpha_pred),
            self.tiny_kind,
            self.tiny.iter().map(|t| format!("[{}, {}, {}, {}]", num(t.0), num(t.1), num(t.2), t.3)).collect::<Vec<_>>().join(", ")
        )
    }
}
