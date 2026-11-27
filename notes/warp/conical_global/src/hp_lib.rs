//! hp_lib (feature `implicit_lib`): `hp_walk` of `hp.rs` on multishoot's `implicit` facility (PR of issue #190).
//!
//! The recursion, the step (planner rule `taylor::step_from_norms`), the landing rule, the step variable s = (x - x0)/l,
//! the Newton correction (`project`) and the checks of the state come from `multishoot::implicit::Walker`; this file keeps
//! what belongs to the conical problem: N's monomials as jets (`PolyQSrc`), the singularity locator, the checkpoints'
//! states, the trace and the statistics of N, all as in `hp.rs`. The scalar of the walker is `W<P>`, the lane
//! precision's scalar with its division through the lanes, as `hp.rs` divides.
use crate::hp::{lane_rat, n_rel_of, sc_of, w_rel_of, Cp, EvalS, HpCfg, HpOut, HpSing, PolyQ, Stop, TraceRow};
use crate::walk::{fit_radius_w, Est};
use multishoot::implicit::{taylor_shift, Coefficients, ImplicitControl, Jets, Mono, Real, Scalar, Stop as WStop, WalkOpts, Walker};
use multishoot::multidouble::precision::Precision;
use multishoot::multidouble::simd::Native;
use multishoot::{Fl, DD};
use std::fmt::Debug;
use std::ops::{Add, Div, Mul, Neg, Sub};

type Sc<P> = <P as Precision<Native>>::Scalar;

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
    };
    let mut hist: Vec<(f64, f64)> = Vec::new(); // (x*, R) of the locator
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
        let formed = w.form().map(|s| (s.l, [0, 1, 2].map(|f| (0..=kt).map(|k| hi_w::<P>(s.c[f][k])).collect::<Vec<f64>>())));
        let (l, cf) = match formed {
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
        if let Some(rstop) = cfg.rstop {
            let est = estimate_c(&cf, kt, l, cfg.fit_window);
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
    out
}
