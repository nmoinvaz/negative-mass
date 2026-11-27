//! hp: the implicit Taylor walk of `poly.rs` and `walk.rs` written once over multidouble's `Precision`
//! (double-, quad- and octo-double), for the walks of notes/warp/conical_loglog.md (#184).
//!
//! What differs from `poly.rs`:
//! - every constant is a lane-precision constant: 1/k, the shifted coefficients of N about the step's x0 (formed in the
//!   lane from the exact integers and the f64 x0), and the division by dN/dF'' (`HpCfg::recip_dd` rounds them to
//!   double-double, the control of the reciprocal check);
//! - the coefficients of N are exact rationals num/den, formed in the lane;
//! - the step planner is the last-two-terms rule of multishoot relative to each field's own magnitude;
//! - `project`: one Newton correction of F'' at each step's start (N = 0 to the rounding of the current scale), against
//!   `conserve`, in which order 0 of N is not imposed and N stays constant along the flow, with the rounding of each
//!   step in it.
use crate::model::Mutation;
use crate::tables::N_TERMS;
use crate::walk::{fit_radius_w, Est};
use multishoot::lane::BaseLane;
use multishoot::multidouble::precision::Precision;
use multishoot::multidouble::simd::Native;
use multishoot::{cs, Fl, DD};
use std::collections::BTreeMap;

type Ln<P> = <P as Precision<Native>>::V;
type Sc<P> = <P as Precision<Native>>::Scalar;
type F1<P> = Fl<Ln<P>, 1>;

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
        let mut map: BTreeMap<(usize, usize, usize), Vec<(usize, i64, i64)>> = BTreeMap::new();
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

/// the scalar rounded to double-double and back (the control of the reciprocal check)
fn rd_l<P: Precision<Native>>(l: Ln<P>) -> Ln<P> {
    P::splat(P::from_dd(P::to_dd(sc_of::<P>(l))))
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

/// the coefficients of the polynomial sum coef x^m in h = x - x0: e_j = sum_m coef C(m, j) x0^(m - j), in the lane
fn shifted<P: Precision<Native>>(xp: &[(usize, i64, i64)], x0: Ln<P>, l: f64) -> Vec<Ln<P>> {
    let deg = xp.iter().map(|&(m, _, _)| m).max().unwrap_or(0);
    let mut pw = vec![lane_i::<P>(1)];
    for i in 1..=deg {
        pw.push(pw[i - 1] * x0);
    }
    // the step variable is s = h / l with l a power of two: e_j l^j is exact up to the lane's rounding of e_j
    let mut lp = vec![lane_i::<P>(1)];
    for i in 1..=deg {
        lp.push(lp[i - 1].mulf(l));
    }
    let mut e = vec![lane_i::<P>(0); deg + 1];
    for &(m, num, den) in xp {
        let co = lane_rat::<P>(num, den);
        let mut binom = 1.0f64;
        for (j, ej) in e.iter_mut().enumerate().take(m + 1) {
            *ej = *ej + co * lane_f::<P>(binom) * pw[m - j];
            binom = binom * (m - j) as f64 / (j + 1) as f64;
        }
    }
    for (j, ej) in e.iter_mut().enumerate() {
        *ej = *ej * lp[j];
    }
    e
}

/// N, the sum of the absolute values of its terms, and dN/dF'' at a state, in the scalar
pub struct EvalS<P: Precision<Native>> {
    groups: Vec<(usize, usize, usize, Vec<(usize, Sc<P>)>)>,
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

/// the series of one step: the fields, their powers and products, and the pieces of N, each to order k as it is formed
struct Cx<'a, P: Precision<Native>> {
    p: &'a PolyQ,
    f: [Vec<F1<P>>; 3],
    pw0: Vec<Vec<F1<P>>>,
    pw1: Vec<Vec<F1<P>>>,
    pw2: Vec<Vec<F1<P>>>,
    unit: Vec<F1<P>>,
    pairs: BTreeMap<(usize, usize), Vec<F1<P>>>,
    q: Vec<Vec<F1<P>>>,
    ex: Vec<Vec<Ln<P>>>,
}

impl<'a, P: Precision<Native>> Cx<'a, P> {
    fn new(p: &'a PolyQ, x0: f64, l: f64, kt: usize, s0: [F1<P>; 3], recip_dd: bool) -> Self {
        let len = kt + 2;
        let z = || vec![F1::<P>::zero(); len];
        let mut f = [z(), z(), z()];
        for i in 0..3 {
            f[i][0] = s0[i];
        }
        let mut unit = z();
        unit[0] = Fl([lane_i::<P>(1)]);
        let mut pairs: BTreeMap<(usize, usize), Vec<F1<P>>> = BTreeMap::new();
        for g in &p.groups {
            if g.a >= 1 && g.b >= 1 {
                pairs.entry((g.a, g.b)).or_insert_with(z);
            }
        }
        let x0l = lane_f::<P>(x0);
        let ex = p
            .groups
            .iter()
            .map(|g| {
                let e = shifted::<P>(&g.xp, x0l, l);
                if recip_dd {
                    e.into_iter().map(rd_l::<P>).collect()
                } else {
                    e
                }
            })
            .collect();
        Cx {
            p,
            f,
            pw0: (0..=p.amax).map(|_| z()).collect(),
            pw1: (0..=p.bmax).map(|_| z()).collect(),
            pw2: (0..=p.cmax.max(1)).map(|_| z()).collect(),
            unit,
            pairs,
            q: p.groups.iter().map(|g| if g.c > 0 { z() } else { Vec::new() }).collect(),
            ex,
        }
    }

    fn base(&self, g: &GroupQ) -> &[F1<P>] {
        match (g.a, g.b) {
            (0, 0) => &self.unit,
            (a, 0) => &self.pw0[a],
            (0, b) => &self.pw1[b],
            (a, b) => &self.pairs[&(a, b)],
        }
    }

    fn partial(&mut self, k: usize) {
        let p = self.p;
        for a in 1..=p.amax {
            self.pw0[a][k] = if a == 1 { self.f[0][k] } else { cs(&self.pw0[a - 1], &self.f[0], k) };
        }
        for b in 1..=p.bmax {
            self.pw1[b][k] = if b == 1 { self.f[1][k] } else { cs(&self.pw1[b - 1], &self.f[1], k) };
        }
        self.pw2[1][k] = self.f[2][k];
        for c in 2..=p.cmax {
            let v = cs(&self.pw2[c - 1], &self.f[2], k);
            self.pw2[c][k] = v;
        }
        let keys: Vec<(usize, usize)> = self.pairs.keys().copied().collect();
        for (a, b) in keys {
            let v = cs(&self.pw0[a], &self.pw1[b], k);
            self.pairs.get_mut(&(a, b)).unwrap()[k] = v;
        }
        for (gi, g) in p.groups.iter().enumerate() {
            if g.c > 0 {
                let v = cs(self.base(g), &self.pw2[g.c], k);
                self.q[gi][k] = v;
            }
        }
    }

    fn n_order(&self, k: usize) -> F1<P> {
        let mut r = F1::<P>::zero();
        for (gi, g) in self.p.groups.iter().enumerate() {
            let src: &[F1<P>] = if g.c > 0 { &self.q[gi] } else { self.base(g) };
            for (j, &e) in self.ex[gi].iter().enumerate() {
                if j <= k {
                    r = r + src[k - j].mul_t(e);
                }
            }
        }
        r
    }

    fn a0(&self) -> F1<P> {
        let mut a0 = F1::<P>::zero();
        for (gi, g) in self.p.groups.iter().enumerate() {
            if g.c == 0 {
                continue;
            }
            let pc = if g.c == 1 { Fl([lane_i::<P>(1)]) } else { self.pw2[g.c - 1][0] };
            a0 = a0 + (self.base(g)[0] * pc).mul_t(self.ex[gi][0]).mulf(g.c as f64);
        }
        a0
    }

    fn fix(&mut self, k: usize, dk: F1<P>) {
        self.f[2][k] = dk;
        let t2 = self.f[2][0] * dk;
        let corr = [F1::<P>::zero(), dk, t2.mulf(2.0), (self.f[2][0] * t2).mulf(3.0)];
        self.pw2[1][k] = dk;
        for c in 2..=self.p.cmax {
            self.pw2[c][k] = self.pw2[c][k] + corr[c];
        }
        for gi in 0..self.p.groups.len() {
            let g = &self.p.groups[gi];
            if g.c > 0 {
                let b0 = self.base(g)[0];
                self.q[gi][k] = self.q[gi][k] + b0 * corr[g.c];
            }
        }
    }
}

/// 1/k for k = 0..=kt in the lane (index 0 unused), rounded to double-double for the control
pub fn recips<P: Precision<Native>>(kt: usize, dd: bool) -> Vec<Ln<P>> {
    (0..=kt + 1)
        .map(|k| {
            if k == 0 {
                lane_i::<P>(0)
            } else {
                let r = lane_i::<P>(1) / lane_i::<P>(k as i64);
                if dd {
                    rd_l::<P>(r)
                } else {
                    r
                }
            }
        })
        .collect()
}

/// Taylor coefficients 1..=kt of (F, F', F'') about x0 in the variable s = (x - x0)/l, l a power of two, from the
/// order-0 values c[.][0]: order k of N = 0 solved for F''_k. `ik[k]` is l/k.
pub fn taylor_p<P: Precision<Native>>(p: &PolyQ, x0: f64, l: f64, c: &mut [Vec<F1<P>>; 3], ik: &[Ln<P>], cfg: &HpCfg) {
    let kt = cfg.kt;
    let mut cx = Cx::<P>::new(p, x0, l, kt, [c[0][0], c[1][0], c[2][0]], cfg.recip_dd);
    cx.partial(0);
    let mut inv = Fl([lane_i::<P>(1)]) / cx.a0();
    if cfg.recip_dd {
        inv = Fl([rd_l::<P>(inv.0[0])]);
    }
    for k in 1..=kt {
        cx.f[0][k] = cx.f[1][k - 1].mul_t(ik[k]);
        cx.f[1][k] = cx.f[2][k - 1].mul_t(ik[k]);
        cx.f[2][k] = F1::<P>::zero();
        cx.partial(k);
        let r = cx.n_order(k);
        let dk = if cfg.skip == Some(k) { F1::<P>::zero() } else { -(r * inv) };
        cx.fix(k, dk);
    }
    for (i, fi) in cx.f.iter().enumerate() {
        for k in 0..=kt {
            c[i][k] = fi[k];
        }
    }
}

fn horner_p<P: Precision<Native>>(s: &[F1<P>], h: Ln<P>) -> F1<P> {
    let mut out = F1::<P>::zero();
    for k in (0..s.len()).rev() {
        out = out.mul_t(h) + s[k];
    }
    out
}

fn fhi<P: Precision<Native>>(f: &F1<P>) -> f64 {
    hi::<P>(P::field_get(f, 0))
}

fn finite<P: Precision<Native>>(f: &F1<P>) -> bool {
    let s = P::field_get(f, 0);
    (0..P::PARTS).all(|q| P::part(s, q).is_finite())
}

/// the best radius estimate of `walk::estimate` over the three fields, from the coefficients in the step variable
/// s = (x - x0)/l; r is returned in x units (times l)
fn estimate_hp<P: Precision<Native>>(c: &[Vec<F1<P>>; 3], kt: usize, l: f64, window: usize) -> Option<Est> {
    let mut best: Option<Est> = None;
    for f in 0..3 {
        let cf: Vec<f64> = (0..=kt).map(|k| fhi::<P>(&c[f][k])).collect();
        if let Some(mut e) = fit_radius_w(&cf, kt, f, window) {
            if e.resid < 1e-3 && best.is_none_or(|b| e.resid < b.resid) {
                e.r *= l;
                best = Some(e);
            }
        }
    }
    best
}

/// the step: the k-th roots of tol |c0| / |c_k| for the last two orders and each field (absolute for a zero c0)
fn plan_h<P: Precision<Native>>(c: &[Vec<F1<P>>; 3], cfg: &HpCfg) -> f64 {
    let kt = cfg.kt;
    let mut h = cfg.hmax;
    for f in 0..3 {
        let s0 = fhi::<P>(&c[f][0]).abs();
        let scale = if s0 > 0.0 { s0 } else { 1.0 };
        for k in [kt - 1, kt] {
            let e = fhi::<P>(&c[f][k]).abs();
            if e.is_nan() {
                return f64::NAN;
            }
            if e > 0.0 {
                let hk = (cfg.tol * scale / e).powf(1.0 / k as f64);
                if hk.is_nan() {
                    return f64::NAN;
                }
                h = h.min(hk);
            }
        }
    }
    h
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

/// Walk from x0 with the state s0 to x1, recording the state at each checkpoint (steps are cut to land on them).
#[cfg(not(feature = "implicit_lib"))]
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
    let ik = recips::<P>(kt, cfg.recip_dd);
    let dir: f64 = if x1 > x0 { 1.0 } else { -1.0 };
    let mut targets: Vec<f64> = cps.iter().copied().filter(|&c| (c - x0) * dir > 0.0 && (x1 - c) * dir > 0.0).collect();
    targets.sort_by(|a, b| (a * dir).partial_cmp(&(b * dir)).unwrap());
    targets.push(x1);
    let mut ti = 0usize;
    let mut c: [Vec<F1<P>>; 3] = std::array::from_fn(|_| vec![F1::<P>::zero(); kt + 2]);
    let (mut x, mut s) = (x0, s0);
    let mut l_prev: f64 = if x0.abs() >= 1.0 { 2f64.powi(x0.abs().log2().floor() as i32) } else { 1.0 };
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
    };
    let mut hist: Vec<(f64, f64)> = Vec::new(); // (x*, R) of the locator
    let mut first = true;
    loop {
        let rest_total = (x1 - x) * dir;
        if rest_total <= 0.0 {
            break;
        }
        if out.steps >= cfg.max_steps {
            out.stop = Stop::StepLimit;
            break;
        }
        let (n_pre, nv, a0) = n_rel_of::<P>(&ev, &s, x);
        if !first {
            out.n_rel_pre_max = out.n_rel_pre_max.max(n_pre);
        }
        if cfg.project && hi::<P>(a0) != 0.0 {
            let corr = Fl([P::splat(nv)]) / Fl([P::splat(a0)]);
            s[2] = s[2] - P::field_get(&corr, 0);
        }
        let (n_post, _, _) = n_rel_of::<P>(&ev, &s, x);
        if !first {
            out.n_rel_max = out.n_rel_max.max(if cfg.project { n_post } else { n_pre });
        }
        for f in 0..3 {
            c[f][0] = Fl([P::splat(s[f])]);
        }
        let l = l_prev;
        let ik_l: Vec<Ln<P>> = ik.iter().map(|v| v.mulf(l)).collect();
        taylor_p::<P>(poly, x, l, &mut c, &ik_l, cfg);
        if c.iter().any(|cf| cf[..=kt].iter().any(|v| !finite::<P>(v))) {
            out.stop = Stop::NonFinite;
            break;
        }
        if let Some(rstop) = cfg.rstop {
            let est = estimate_hp::<P>(&c, kt, l, cfg.fit_window);
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
        let mut h = plan_h::<P>(&c, cfg) * l * cfg.hscale;
        if !(h.is_finite() && h > 0.0) {
            out.stop = Stop::NonFinite;
            break;
        }
        if cfg.rstop.is_some() {
            let cf: Vec<f64> = (0..=kt).map(|k| fhi::<P>(&c[2][k])).collect();
            let mut rs: Vec<f64> = (kt - 4..=kt).map(|k| l / (cf[k] / cf[k - 1]).abs()).collect();
            rs.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
            let sc = (kt - 9..=kt).filter(|&k| k >= 2 && (cf[k] / cf[k - 1]) < 0.0).count();
            out.recent.push((x, h, rs[2], sc));
            if out.recent.len() > 8 {
                out.recent.remove(0);
            }
        }
        let target = targets[ti];
        let rest = (target - x) * dir;
        if h >= rest {
            h = rest;
        } else if h > 0.6 * rest {
            h = 0.5 * rest;
        }
        let xn = if h == rest { target } else { x + dir * h };
        if xn == x || h < cfg.tiny * x.abs().max(1e-6) {
            out.stop = Stop::TinyStep;
            break;
        }
        let hd = (P::new(&[xn]) - P::new(&[x])) * P::new(&[1.0 / l]);
        let hl = P::splat(hd);
        let ns: Vec<F1<P>> = (0..3).map(|f| horner_p::<P>(&c[f][..=kt], hl)).collect();
        if ns.iter().any(|f| !finite::<P>(f)) {
            out.stop = Stop::NonFinite;
            break;
        }
        if want_trace {
            out.trace.push(TraceRow {
                x,
                s_hi: [hi::<P>(s[0]), hi::<P>(s[1]), hi::<P>(s[2])],
                w_rel: w_rel_of::<P>(&s, x),
                n_rel: n_post,
                h,
            });
        }
        s = [P::field_get(&ns[0], 0), P::field_get(&ns[1], 0), P::field_get(&ns[2], 0)];
        x = xn;
        // the next step's scale: the power of two nearest this step's length (the series' radius is a few of them)
        l_prev = 2f64.powi((h / cfg.hscale).log2().round() as i32).max(1e-300);
        out.steps += 1;
        first = false;
        out.x_end = x;
        out.s_end = s;
        if xn == target {
            let (nr, _, _) = n_rel_of::<P>(&ev, &s, x);
            if ti + 1 < targets.len() {
                out.cps.push(Cp { x, s, n_rel: nr, step: out.steps });
                ti += 1;
            }
        }
    }
    out
}

/// DD constant helper for the CLI (start states given as parts)
pub fn scalar_of_parts<P: Precision<Native>>(parts: &[f64]) -> Sc<P> {
    P::new(parts)
}

pub fn dd_to_sc<P: Precision<Native>>(d: DD) -> Sc<P> {
    P::from_dd(d)
}

#[cfg(feature = "implicit_lib")]
pub use crate::hp_lib::hp_walk;
