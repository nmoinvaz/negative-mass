//! hp: the implicit Taylor walk of notes/warp/conical_loglog.md (#184) over multidouble's `Precision` (double-, quad-
//! and octo-double), on `multishoot::implicit` (issue #190).
//!
//! The recursion, the step (planner rule `taylor::step_from_norms`), the landing rule, the step variable s = (x - x0)/l,
//! the Newton correction (`project`) and the checks of the state come from `multishoot::implicit::Walker`. This file
//! keeps what belongs to the conical problem: N's monomials as jets (`PolyQSrc`), the singularity locator, the
//! checkpoints' states, the trace and the statistics of N. The walker's scalar is `W<P>`, the lane precision's scalar
//! with its division through the lanes. Further:
//! - every constant is a lane-precision constant: the shifted coefficients of N about the step's x0 are formed in the
//!   lane from the exact integers and the f64 x0 (`HpCfg::recip_dd` rounds them and the reciprocals to double-double,
//!   the control of the reciprocal check);
//! - the coefficients of N are exact rationals num/den, formed in the lane;
//! - `project`: one Newton correction of F'' at each step's start (N = 0 to the rounding of the current scale), against
//!   `conserve`, in which order 0 of N is not imposed and N stays constant along the flow, with the rounding of each
//!   step in it.
use crate::model::Mutation;
use crate::tables::N_TERMS;
use crate::walk::{fit_radius_w, Est, LocStats};
#[cfg(feature = "locate_lib")]
use crate::walk::Hist;
use multishoot::implicit::{
    taylor_shift, Coefficients, ImplicitControl, Jets, Mono, Real, Scalar, Stop as WStop, WalkOpts, Walker,
};
use multishoot::lane::BaseLane;
use multishoot::multidouble::precision::Precision;
use multishoot::multidouble::simd::Native;
use multishoot::{Fl, DD};
use std::collections::BTreeMap;
use std::fmt::Debug;
use std::ops::{Add, Div, Mul, Neg, Sub};

type Ln<P> = <P as Precision<Native>>::V;
type Sc<P> = <P as Precision<Native>>::Scalar;
/// (m, num, den) of the terms num/den x^m of a group
type Xp = Vec<(usize, i64, i64)>;
/// (a, b, c, terms as (m, coefficient)) of a group in the scalar
type GroupS<P> = (usize, usize, usize, Vec<(usize, Sc<P>)>);

/// N = sum (num/den) F0^a F1^b F2^c x^m with c <= 3, grouped by (a, b, c)
#[derive(Clone, Debug)]
pub struct PolyQ {
    pub groups: Vec<GroupQ>,
    pub amax: usize,
    pub bmax: usize,
    pub cmax: usize,
}

#[derive(Clone, Debug)]
pub struct GroupQ {
    pub a: usize,
    pub b: usize,
    pub c: usize,
    /// (m, num, den)
    pub xp: Vec<(usize, i64, i64)>,
}

impl PolyQ {
    pub fn new(terms: &[(u8, u8, u8, u8, i64, i64)]) -> Self {
        let mut map: BTreeMap<(usize, usize, usize), Xp> = BTreeMap::new();
        for &(a, b, c, m, num, den) in terms {
            assert!(den != 0);
            map.entry((a as usize, b as usize, c as usize)).or_default().push((m as usize, num, den));
        }
        let groups: Vec<GroupQ> = map.into_iter().map(|((a, b, c), xp)| GroupQ { a, b, c, xp }).collect();
        let amax = groups.iter().map(|g| g.a).max().unwrap_or(0);
        let bmax = groups.iter().map(|g| g.b).max().unwrap_or(0);
        let cmax = groups.iter().map(|g| g.c).max().unwrap_or(0);
        assert!(cmax <= 3);
        PolyQ { groups, amax, bmax, cmax }
    }

    /// N of `det_classify.conical_N()` (the crate's table), with the control `eq-coeff` of #167
    pub fn conical(m: &Mutation) -> Self {
        let t: Vec<_> = N_TERMS
            .iter()
            .map(|&(a, b, c, mm, co)| {
                let co = if m.eq_coeff && (a, b, c, mm) == (1, 2, 1, 2) { co + 1 } else { co };
                (a, b, c, mm, co, 1)
            })
            .collect();
        PolyQ::new(&t)
    }
}

pub fn lane_i<P: Precision<Native>>(v: i64) -> Ln<P> {
    P::splat(P::new(&[v as f64]))
}

pub fn lane_f<P: Precision<Native>>(v: f64) -> Ln<P> {
    P::splat(P::new(&[v]))
}

pub fn lane_rat<P: Precision<Native>>(num: i64, den: i64) -> Ln<P> {
    if den == 1 {
        lane_i::<P>(num)
    } else {
        lane_i::<P>(num) / lane_i::<P>(den)
    }
}

pub fn sc_of<P: Precision<Native>>(l: Ln<P>) -> Sc<P> {
    P::field_get(&Fl([l]), 0)
}

pub fn hi<P: Precision<Native>>(x: Sc<P>) -> f64 {
    P::part(x, 0)
}

pub fn abs_s<P: Precision<Native>>(x: Sc<P>) -> Sc<P> {
    if P::part(x, 0) < 0.0 || (P::part(x, 0) == 0.0 && P::part(x, 1) < 0.0) {
        -x
    } else {
        x
    }
}

fn pow_s<P: Precision<Native>>(x: Sc<P>, n: usize) -> Sc<P> {
    let mut r = P::ONE;
    for _ in 0..n {
        r = r * x;
    }
    r
}

/// N, the sum of the absolute values of its terms, and dN/dF'' at a state, in the scalar
pub struct EvalS<P: Precision<Native>> {
    groups: Vec<GroupS<P>>,
}

impl<P: Precision<Native>> EvalS<P> {
    pub fn new(p: &PolyQ) -> Self {
        let groups = p
            .groups
            .iter()
            .map(|g| {
                (g.a, g.b, g.c, g.xp.iter().map(|&(m, n, d)| (m, sc_of::<P>(lane_rat::<P>(n, d)))).collect())
            })
            .collect();
        EvalS { groups }
    }

    pub fn eval(&self, s: &[Sc<P>; 3], x: Sc<P>) -> (Sc<P>, Sc<P>, Sc<P>) {
        let (mut v, mut sm, mut d) = (P::ZERO, P::ZERO, P::ZERO);
        for (a, b, c, xp) in &self.groups {
            let base = pow_s::<P>(s[0], *a) * pow_s::<P>(s[1], *b);
            let bc = base * pow_s::<P>(s[2], *c);
            let bd = if *c >= 1 { base * pow_s::<P>(s[2], *c - 1) } else { P::ZERO };
            for (m, co) in xp {
                let xm = pow_s::<P>(x, *m);
                let t = *co * bc * xm;
                v = v + t;
                sm = sm + abs_s::<P>(t);
                if *c >= 1 {
                    d = d + *co * bd * xm * P::new(&[*c as f64]);
                }
            }
        }
        (v, sm, d)
    }
}

#[derive(Clone, Debug)]
pub struct HpCfg {
    pub kt: usize,
    pub tol: f64,
    pub project: bool,
    /// constants rounded to double-double (control of the reciprocal check)
    pub recip_dd: bool,
    /// every planned step is scaled by this (0.5 halves every step)
    pub hscale: f64,
    pub hmax: f64,
    pub max_steps: usize,
    /// control `n-drift`: F''_k left zero at this order
    pub skip: Option<usize>,
    /// the singularity locator of `walk.rs` (#183): stop when the radius estimate is real, ahead and below
    /// rstop * max(1, |x|) over five consecutive steps (relative beyond x = 1); None leaves the walk as in #184
    pub rstop: Option<f64>,
    /// false: the locator's bound is the absolute rstop (as `walk.rs`); true: rstop * max(1, |x|)
    pub rstop_rel: bool,
    /// the number of last Taylor orders the locator's 1/k fit uses (10 as `walk.rs`)
    pub fit_window: usize,
    /// the walk stops (`Stop::TinyStep`) when a step is below tiny * max(|x|, 1e-6); #184's value 1e-13, #183 uses 1e-15
    pub tiny: f64,
    /// control `plant-locator`: the locator's side test inverted
    pub flip_side: bool,
}

impl HpCfg {
    pub fn new(kt: usize, tol: f64) -> Self {
        HpCfg { kt, tol, project: true, recip_dd: false, hscale: 1.0, hmax: 1e300, max_steps: 200000, skip: None, rstop: None, rstop_rel: true, fit_window: 10, tiny: 1e-13, flip_side: false }
    }
}

#[derive(Clone, Debug, PartialEq)]
pub enum Stop {
    End,
    NonFinite,
    TinyStep,
    StepLimit,
    /// the locator found a real singularity ahead (`HpOut::sing`; only with `HpCfg::rstop`)
    Singular,
}

/// a singularity found by the locator of an `hp_walk`: x*, the spread of the last five estimates, the estimate (r in x
/// units), the last step's x and its state in double-double
#[derive(Clone, Debug)]
pub struct HpSing {
    pub xstar: f64,
    pub bar: f64,
    pub est: Est,
    pub x_last: f64,
    pub s_last: [DD; 3],
}

#[derive(Clone, Debug)]
pub struct Cp<P: Precision<Native>> {
    pub x: f64,
    pub s: [Sc<P>; 3],
    pub n_rel: f64,
    pub step: usize,
}

#[derive(Clone, Debug)]
pub struct TraceRow {
    pub x: f64,
    pub s_hi: [f64; 3],
    /// (q + p^2) / p^2 at the state, in the lane
    pub w_rel: f64,
    pub n_rel: f64,
    pub h: f64,
}

pub struct HpOut<P: Precision<Native>> {
    pub cps: Vec<Cp<P>>,
    pub trace: Vec<TraceRow>,
    pub stop: Stop,
    pub x_end: f64,
    pub s_end: [Sc<P>; 3],
    pub steps: usize,
    /// max over the step ends after the first step of n_rel after the projection (project) or at the end (conserve)
    pub n_rel_max: f64,
    /// project mode: max over the steps of n_rel of the Horner-advanced state before the projection
    pub n_rel_pre_max: f64,
    /// the singularity the locator found (only with `HpCfg::rstop`)
    pub sing: Option<HpSing>,
    /// with `HpCfg::rstop`: the last steps' (x, h, median radius of F'' from the last five Taylor ratios in x units, sign
    /// changes of the last ten ratios), oldest first (the diagnostic of a `Stop::TinyStep`)
    pub recent: Vec<(f64, f64, f64, usize)>,
    /// what the locator did (the calls of `multisolver::locate` and the f64 fit beside them; zero without `locate_lib`)
    pub loc: LocStats,
}

pub(crate) fn n_rel_of<P: Precision<Native>>(ev: &EvalS<P>, s: &[Sc<P>; 3], x: f64) -> (f64, Sc<P>, Sc<P>) {
    let (v, sm, d) = ev.eval(s, P::new(&[x]));
    let r = if hi::<P>(sm) == 0.0 { 0.0 } else { hi::<P>(abs_s::<P>(v)) / hi::<P>(sm) };
    (r, v, d)
}

/// (q + p^2)/p^2 with p = x F', q = x^2 F'' + x F', in the scalar
pub(crate) fn w_rel_of<P: Precision<Native>>(s: &[Sc<P>; 3], x: f64) -> f64 {
    let xs = P::new(&[x]);
    let p = xs * s[1];
    let q = xs * xs * s[2] + p;
    let num = q + p * p;
    let den = p * p;
    if hi::<P>(den) == 0.0 {
        return f64::NAN;
    }
    let r = Fl([P::splat(num)]) / Fl([P::splat(den)]);
    hi::<P>(P::field_get(&r, 0))
}

/// DD constant helper for the CLI (start states given as parts)
pub fn scalar_of_parts<P: Precision<Native>>(parts: &[f64]) -> Sc<P> {
    P::new(parts)
}

pub fn dd_to_sc<P: Precision<Native>>(d: DD) -> Sc<P> {
    P::from_dd(d)
}

/// the scalar of lane precision P as the library's number type
#[derive(Clone, Copy)]
pub struct W<P: Precision<Native>>(pub Sc<P>);

impl<P: Precision<Native>> Debug for W<P> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.0.fmt(f)
    }
}

impl<P: Precision<Native>> Add for W<P> {
    type Output = W<P>;
    fn add(self, o: W<P>) -> W<P> {
        W(self.0 + o.0)
    }
}

impl<P: Precision<Native>> Sub for W<P> {
    type Output = W<P>;
    fn sub(self, o: W<P>) -> W<P> {
        W(self.0 - o.0)
    }
}

impl<P: Precision<Native>> Mul for W<P> {
    type Output = W<P>;
    fn mul(self, o: W<P>) -> W<P> {
        W(self.0 * o.0)
    }
}

impl<P: Precision<Native>> Neg for W<P> {
    type Output = W<P>;
    fn neg(self) -> W<P> {
        W(-self.0)
    }
}

impl<P: Precision<Native>> Div for W<P> {
    type Output = W<P>;
    fn div(self, o: W<P>) -> W<P> {
        let r = Fl([P::splat(self.0)]) / Fl([P::splat(o.0)]);
        W(P::field_get(&r, 0))
    }
}

impl<P: Precision<Native>> Scalar for W<P> {
    fn zero() -> Self {
        W(P::ZERO)
    }
    fn one() -> Self {
        W(P::ONE)
    }
    fn from_f64(x: f64) -> Self {
        W(P::new(&[x]))
    }
    fn from_dd(x: DD) -> Self {
        W(P::from_dd(x))
    }
    fn mulf(self, s: f64) -> Self {
        W(self.0 * P::new(&[s]))
    }
    fn round_dd(self) -> Self {
        W(P::from_dd(P::to_dd(self.0)))
    }
    fn is_finite(self) -> bool {
        P::part(self.0, 0).is_finite()
    }
    /// the Cauchy sum in the lane's dot accumulator, as `field::cs` (which `hp.rs` used) forms it
    fn cauchy(a: &[Self], b: &[Self], k: usize) -> Self {
        let mut acc = <P::V as BaseLane>::acc0();
        for j in 0..=k {
            acc = <P::V as BaseLane>::acc(acc, P::splat(a[j].0), P::splat(b[k - j].0));
        }
        W(P::field_get(&Fl([<P::V as BaseLane>::fin(acc)]), 0))
    }
}

impl<P: Precision<Native>> Real for W<P> {
    const EPS: f64 = match P::PARTS {
        2 => 1.232595164407831e-32,
        4 => 1.504632769052528e-64,
        _ => 2.3e-128,
    };
    fn to_f64(self) -> f64 {
        P::part(self.0, 0)
    }
}

/// the lane's scalar as the arithmetic of `multisolver::locate`, so that the locator runs in the precision of the series
#[cfg(feature = "locate_lib")]
impl<P: Precision<Native>> multisolver::locate::Real for W<P> {
    const UNIT: f64 = match P::PARTS {
        2 => 4.930380657631324e-32,
        4 => 2.4308653429145085e-63,
        _ => 5.909106315382871e-126,
    };
    fn from_f64(x: f64) -> Self {
        W(P::new(&[x]))
    }
    fn to_f64(self) -> f64 {
        P::part(self.0, 0)
    }
    /// Newton from the f64 root: each step doubles the digits (2, 3 and 4 steps for double-, quad- and octo-double)
    fn sqrt(self) -> Self {
        let a = P::part(self.0, 0);
        if !(a > 0.0 && a.is_finite()) {
            return W(P::ZERO);
        }
        let half = W(P::new(&[0.5]));
        let mut y = W(P::new(&[a.sqrt()]));
        for _ in 0..P::PARTS / 2 + 1 {
            y = half * (y + self / y);
        }
        y
    }
}

/// the three fields' coefficients in the lane's precision (what `locate` reads); empty without `locate_lib`
pub struct Full<P: Precision<Native>> {
    #[cfg(feature = "locate_lib")]
    pub c: [Vec<W<P>>; 3],
    #[cfg(not(feature = "locate_lib"))]
    _p: std::marker::PhantomData<P>,
}

fn full_of<P: Precision<Native>>(c: &[Vec<W<P>>; 3], kt: usize) -> Full<P> {
    #[cfg(feature = "locate_lib")]
    {
        Full { c: [0, 1, 2].map(|f| c[f][..=kt].to_vec()) }
    }
    #[cfg(not(feature = "locate_lib"))]
    {
        let _ = (c, kt);
        Full { _p: std::marker::PhantomData }
    }
}

/// N's monomials as the library's source: the coefficient of each (a, b, c) is the polynomial in x of `PolyQ`, its jet the
/// Taylor shift about x0 in the lane's precision (rounded to double-double for the control `recip-dd`)
pub struct PolyQSrc<'a, P: Precision<Native>> {
    poly: &'a PolyQ,
    coefs: Vec<Vec<W<P>>>,
    recip_dd: bool,
}

impl<'a, P: Precision<Native>> PolyQSrc<'a, P> {
    pub fn new(poly: &'a PolyQ, recip_dd: bool) -> Self {
        let coefs = poly
            .groups
            .iter()
            .map(|g| {
                let deg = g.xp.iter().map(|&(m, _, _)| m).max().unwrap_or(0);
                let mut c = vec![W::<P>(P::ZERO); deg + 1];
                for &(m, num, den) in &g.xp {
                    c[m] = c[m] + W(sc_of::<P>(lane_rat::<P>(num, den)));
                }
                c
            })
            .collect();
        PolyQSrc { poly, coefs, recip_dd }
    }
}

impl<P: Precision<Native>> Coefficients<W<P>> for PolyQSrc<'_, P> {
    fn pattern(&self) -> Vec<Mono> {
        self.poly.groups.iter().map(|g| (g.a, g.b, g.c)).collect()
    }
    fn jets(&self, x0: f64, _kt: usize, out: &mut Jets<W<P>>) {
        for (g, c) in self.coefs.iter().enumerate() {
            taylor_shift(c, x0, &mut out.ent[g].c);
            out.ent[g].n = c.len();
            if self.recip_dd {
                for v in out.ent[g].c.iter_mut().take(c.len()) {
                    *v = v.round_dd();
                }
            }
        }
    }
}

fn hi_w<P: Precision<Native>>(x: W<P>) -> f64 {
    P::part(x.0, 0)
}

/// the best radius estimate of `walk::estimate` over the three fields, from the coefficients in the step variable
/// s = (x - x0)/l; r is returned in x units (times l)
fn estimate_c(c: &[Vec<f64>; 3], kt: usize, l: f64, window: usize) -> Option<Est> {
    let mut best: Option<Est> = None;
    for (f, cf) in c.iter().enumerate() {
        if let Some(mut e) = fit_radius_w(cf, kt, f, window) {
            if e.resid < 1e-3 && best.is_none_or(|b| e.resid < b.resid) {
                e.r *= l;
                best = Some(e);
            }
        }
    }
    best
}

/// Walk from x0 with the state s0 to x1, recording the state at each checkpoint (steps are cut to land on them): `hp_walk`
/// of `hp.rs` on `multishoot::implicit::Walker`.
pub fn hp_walk<P: Precision<Native>>(
    poly: &PolyQ,
    s0: [Sc<P>; 3],
    x0: f64,
    x1: f64,
    cps: &[f64],
    cfg: &HpCfg,
    want_trace: bool,
) -> HpOut<P> {
    let kt = cfg.kt;
    let ev = EvalS::<P>::new(poly);
    let src = PolyQSrc::<P>::new(poly, cfg.recip_dd);
    let mut o = WalkOpts::new(kt, cfg.tol);
    o.plan.hmax = cfg.hmax;
    o.branch = if cfg.project { multishoot::implicit::Branch::Project } else { multishoot::implicit::Branch::Conserve };
    o.norm = multishoot::implicit::StepNorm::Rel;
    o.h_scale = cfg.hscale;
    o.tiny = cfg.tiny;
    let mut w = Walker::<W<P>, _>::new(&src, x0, [W(s0[0]), W(s0[1]), W(s0[2])], o);
    if cfg.recip_dd {
        w.set_control(ImplicitControl::RecipDD);
    }
    if let Some(k) = cfg.skip {
        w.set_control(ImplicitControl::SkipOrder(k));
    }
    let dir: f64 = if x1 > x0 { 1.0 } else { -1.0 };
    let mut targets: Vec<f64> = cps.iter().copied().filter(|&c| (c - x0) * dir > 0.0 && (x1 - c) * dir > 0.0).collect();
    targets.sort_by(|a, b| (a * dir).partial_cmp(&(b * dir)).unwrap());
    targets.push(x1);
    let mut ti = 0usize;
    let mut out = HpOut::<P> {
        cps: Vec::new(),
        trace: Vec::new(),
        stop: Stop::End,
        x_end: x0,
        s_end: s0,
        steps: 0,
        n_rel_max: 0.0,
        n_rel_pre_max: 0.0,
        sing: None,
        recent: Vec::new(),
        loc: LocStats::default(),
    };
    let mut hist: Vec<(f64, f64)> = Vec::new(); // (x*, R) of the locator
    #[cfg(feature = "locate_lib")]
    let mut shadow = Hist::default();
    let mut first = true;
    loop {
        let x = w.x();
        let rest_total = (x1 - x) * dir;
        if rest_total <= 0.0 {
            break;
        }
        if out.steps >= cfg.max_steps {
            out.stop = Stop::StepLimit;
            break;
        }
        let s_pre = w.y().map(|v| v.0);
        let (n_pre, _, _) = n_rel_of::<P>(&ev, &s_pre, x);
        if !first {
            out.n_rel_pre_max = out.n_rel_pre_max.max(n_pre);
        }
        let formed = w
            .form()
            .map(|s| (s.l, [0, 1, 2].map(|f| (0..=kt).map(|k| hi_w::<P>(s.c[f][k])).collect::<Vec<f64>>()), full_of::<P>(&s.c, kt)));
        let (l, cf, full) = match formed {
            Ok(v) => v,
            Err(WStop::ZeroPivot | WStop::NonFinite) => {
                out.stop = Stop::NonFinite;
                break;
            }
            Err(e) => unreachable!("walker stop {e:?} with no pivot or drift guard"),
        };
        let s = w.y().map(|v| v.0);
        let (n_post, _, _) = n_rel_of::<P>(&ev, &s, x);
        if !first {
            out.n_rel_max = out.n_rel_max.max(if cfg.project { n_post } else { n_pre });
        }
        let _ = &full;
        if let Some(rstop) = cfg.rstop {
            #[cfg(not(feature = "locate_lib"))]
            let est = estimate_c(&cf, kt, l, cfg.fit_window);
            // `locate` in the lane's precision at the steps inside the gate, the f64 fit of the window beside it
            #[cfg(feature = "locate_lib")]
            let est = {
                let g = crate::loc::Gate {
                    h: w.planned(),
                    bound: if cfg.rstop_rel { rstop * x.abs().max(1.0) } else { rstop },
                    dir,
                    flip: cfg.flip_side,
                    x,
                    steps: out.steps,
                };
                crate::loc::step([&full.c[0], &full.c[1], &full.c[2]], l, estimate_c(&cf, kt, l, cfg.fit_window), &g, &mut out.loc, &mut shadow)
            };
            let ahead = |e: &Est| (e.side as f64 == dir) != cfg.flip_side;
            match est {
                Some(e) if ahead(&e) => hist.push((x + dir * e.r, e.r)),
                _ => hist.clear(),
            }
            if hist.len() >= 5 {
                let last5 = &hist[hist.len() - 5..];
                let (xs, r) = *last5.last().unwrap();
                if r <= if cfg.rstop_rel { rstop * x.abs().max(1.0) } else { rstop } {
                    let lo = last5.iter().map(|h| h.0).fold(f64::INFINITY, f64::min);
                    let hi_ = last5.iter().map(|h| h.0).fold(f64::NEG_INFINITY, f64::max);
                    out.stop = Stop::Singular;
                    out.sing = Some(HpSing {
                        xstar: xs,
                        bar: hi_ - lo,
                        est: est.unwrap(),
                        x_last: x,
                        s_last: [P::to_dd(s[0]), P::to_dd(s[1]), P::to_dd(s[2])],
                    });
                    break;
                }
            }
        }
        let h = w.planned();
        if !(h.is_finite() && h > 0.0) {
            out.stop = Stop::NonFinite;
            break;
        }
        if cfg.rstop.is_some() {
            let c2 = &cf[2];
            let mut rs: Vec<f64> = (kt - 4..=kt).map(|k| l / (c2[k] / c2[k - 1]).abs()).collect();
            rs.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
            let sc = (kt - 9..=kt).filter(|&k| k >= 2 && (c2[k] / c2[k - 1]) < 0.0).count();
            out.recent.push((x, h, rs[2], sc));
            if out.recent.len() > 8 {
                out.recent.remove(0);
            }
        }
        let adv = match w.advance(h, targets[ti]) {
            Ok(a) => a,
            Err(WStop::NotFormed) => unreachable!("advance follows the form of the same state"),
            Err(WStop::Tiny) => {
                out.stop = Stop::TinyStep;
                break;
            }
            Err(_) => {
                out.stop = Stop::NonFinite;
                break;
            }
        };
        if want_trace {
            out.trace.push(TraceRow {
                x,
                s_hi: [P::part(s[0], 0), P::part(s[1], 0), P::part(s[2], 0)],
                w_rel: w_rel_of::<P>(&s, x),
                n_rel: n_post,
                h: adv.h,
            });
        }
        let sn = w.y().map(|v| v.0);
        let xn = w.x();
        out.steps += 1;
        first = false;
        out.x_end = xn;
        out.s_end = sn;
        if adv.landed {
            let (nr, _, _) = n_rel_of::<P>(&ev, &sn, xn);
            if ti + 1 < targets.len() {
                out.cps.push(Cp { x: xn, s: sn, n_rel: nr, step: out.steps });
                ti += 1;
            }
        }
    }
    out.loc.step = out.steps;
    out
}

#[cfg(all(test, feature = "locate_lib"))]
mod tests {
    use super::*;
    use multidouble::precision::{Double, Octo, Quad};
    use multisolver::locate::Real as LReal;

    /// `sqrt` squares back to the argument within `units` of the lane's rounding; zero and a negative argument give zero
    fn check<P: Precision<Native>>(sqrt: impl Fn(W<P>) -> W<P>, units: f64) {
        for a in [2.0, 3.0, 0.3, 1e10, 1e-10] {
            let x = W::<P>(P::new(&[a]));
            let y = sqrt(x);
            let err = ((y * y - x) / x).0;
            assert!(P::part(err, 0).abs() <= units * <W<P> as LReal>::UNIT, "a = {a}: {:e}", P::part(err, 0));
        }
        assert_eq!(P::part(sqrt(W::<P>(P::ZERO)).0, 0), 0.0);
        assert_eq!(P::part(sqrt(W::<P>(P::new(&[-1.0]))).0, 0), 0.0);
    }

    /// the f64 root with no Newton step: the planted error the check must see
    fn f64_root<P: Precision<Native>>(x: W<P>) -> W<P> {
        let a = P::part(x.0, 0);
        W(P::new(&[if a > 0.0 { a.sqrt() } else { 0.0 }]))
    }

    #[test]
    fn newton_sqrt() {
        check::<Double>(LReal::sqrt, 4.0);
        check::<Quad>(LReal::sqrt, 4.0);
        check::<Octo>(LReal::sqrt, 4.0);
        assert!(std::panic::catch_unwind(|| check::<Double>(f64_root, 4.0)).is_err());
        assert!(std::panic::catch_unwind(|| check::<Quad>(f64_root, 4.0)).is_err());
        assert!(std::panic::catch_unwind(|| check::<Octo>(f64_root, 4.0)).is_err());
    }
}
