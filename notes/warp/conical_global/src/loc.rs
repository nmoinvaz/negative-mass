//! The singularity locator of multisolver 0.3.0 (`multisolver::locate`, #191 and #206) in the two walkers, in place of the f64
//! Domb-Sykes fit of `walk.rs` and `hp.rs` (issue #207). Built with the feature `locate_lib` only.
//!
//! At a step the walker hands over the three fields' Taylor coefficients in the precision of its lane (double-double in
//! `walk.rs`, the lane's `W<P>` in `hp.rs`). Each series is scaled by r^k, r a power of two near its radius (the library expects
//! a radius of order one; the scaling is exact), and goes to `locate` with its default options: the window and the degree are
//! the library's. A field that reads `Kind::Real` gives an estimate: the radius |s| r l in x units, the side from the sign of s,
//! gamma = 1 + alpha of the field (the ratios c_(k-1)/c_k = |s| (1 + (1 + alpha)/k + ...)), and the library's bar of s over |s| in
//! place of the fit residual. Of the fields the one with the smallest bar is the step's estimate. A field that reads a pair,
//! `Unresolved` or nothing gives none, as a field whose ratios change sign did for the f64 fit. The stop rule of the walkers (five
//! consecutive real estimates ahead, the last within the stop radius) is unchanged.
//!
//! **The gate.** `locate` fits 119 windows and degrees per field: about 1.2 ms in double-double (K = 30), 21 ms in quad-double
//! (K = 60) and 430 ms in octo-double (K = 120) on the Mac. It runs only at the steps where the planner's step h is at most
//! `GATE` times the stop radius. h is the radius times about 0.09 (the planner's tolerance to the power 1/K), so the gate is
//! open for the last 60 steps or so before the stop radius and for the five that the rule needs; where it is closed the step has
//! no estimate, as a step whose fit failed. The rule is checked against `locate` at every step on the `ctl` grid (raw I207-mac).
//!
//! **The shadow.** On the same steps the f64 fit of the walker runs beside `locate`, with its own five-step rule, and
//! `LocStats` keeps what it said: the step at which it would have stopped and its eta*. The walk's steps do not depend on the
//! locator, so the shadow's stop is the one the f64 walker would have made.
//!
//! **The controls** (cargo features): `locate_mut_early` replaces `locate` by `locate_at` on the ratio model with the first
//! order 2 and degree 2 (a window too early: the farther singularities act), `locate_mut_late` by the last 10 orders, degree 2
//! (a window too short: the fit of the f64 walker, through the library). `locate_ungated` is not a control: it runs `locate` at
//! every step, to check that the gate changes nothing.
use crate::walk::{Est, Hist, LocStats};
use multisolver::locate::{Kind, Location, Opts, Real};

/// the locator runs at a step when the planned step is at most this many stop radii (every step with `locate_ungated`)
#[cfg(not(feature = "locate_ungated"))]
pub const GATE: f64 = 16.0;
#[cfg(feature = "locate_ungated")]
pub const GATE: f64 = f64::INFINITY;

/// where the walker is at the step and what the gate needs
pub struct Gate {
    /// the planner's step in x units
    pub h: f64,
    /// the stop radius in x units
    pub bound: f64,
    pub dir: f64,
    /// the control `locator-side`: the side test inverted
    pub flip: bool,
    pub x: f64,
    /// the walker's steps so far
    pub steps: usize,
}

/// the power of two nearest x on the log scale, None unless x is a positive normal number whose neighbor is one too
fn pow2_near(x: f64) -> Option<f64> {
    if !(x.is_finite() && x >= f64::MIN_POSITIVE) {
        return None;
    }
    let bits = x.to_bits();
    let e = ((bits >> 52) & 0x7ff) as i32 - 1023;
    let m = f64::from_bits((bits & 0x000f_ffff_ffff_ffff) | 0x3ff0_0000_0000_0000);
    let e = if m >= std::f64::consts::SQRT_2 { e + 1 } else { e };
    if !(-1022..=1023).contains(&e) {
        return None;
    }
    Some(f64::from_bits(((e + 1023) as u64) << 52))
}

/// c_k rho^k with rho a power of two near the radius |c_(K-1)/c_K|, and rho; exact
fn scaled<T: Real>(c: &[T]) -> Option<(Vec<T>, f64)> {
    let k = c.len() - 1;
    let rho = pow2_near((c[k - 1].to_f64() / c[k].to_f64()).abs())?;
    let mut out = Vec::with_capacity(c.len());
    let mut p = 1.0f64;
    for ck in c {
        if !(p.is_finite() && (f64::MIN_POSITIVE..=1e300).contains(&p)) {
            return None;
        }
        out.push(*ck * T::from_f64(p));
        p *= rho;
    }
    Some((out, rho))
}

/// the library's reading of one series
fn call<T: Real>(c: &[T]) -> Option<Location<T>> {
    let o = Opts::default();
    #[cfg(feature = "locate_mut_early")]
    {
        multisolver::locate::locate_at(c, multisolver::locate::Model::Ratio, 2, 2, &o)
    }
    #[cfg(feature = "locate_mut_late")]
    {
        multisolver::locate::locate_at(c, multisolver::locate::Model::Ratio, c.len() - 10, 2, &o)
    }
    #[cfg(not(any(feature = "locate_mut_early", feature = "locate_mut_late")))]
    {
        multisolver::locate::locate(c, &o)
    }
}

enum Field {
    /// the estimate, the library's bar of s over |s| and the bar of alpha
    Real(Est, f64, f64),
    Pair,
    Unresolved,
    Nothing,
}

fn read_field<T: Real>(c: &[T], f: usize, l: f64) -> Field {
    let Some((sc, rho)) = scaled(c) else {
        return Field::Nothing;
    };
    let Some(loc) = call(&sc) else {
        return Field::Nothing;
    };
    match loc.kind {
        Kind::Unresolved => Field::Unresolved,
        Kind::Pair => Field::Pair,
        Kind::Real => {
            let s = loc.nearest.re.to_f64();
            let r = s.abs() * rho * l;
            let gamma = 1.0 + loc.nearest.alpha.to_f64();
            if !(r.is_finite() && r > 0.0 && gamma.is_finite() && loc.rel_bar.is_finite()) {
                return Field::Nothing;
            }
            let side = if s > 0.0 { 1 } else { -1 };
            Field::Real(Est { r, gamma, resid: loc.rel_bar, field: f, side }, loc.rel_bar, loc.nearest.alpha_bar)
        }
    }
}

/// The estimate of one step from the library, None outside the gate or with no real field. `l` is the x length of the series'
/// variable (1 in `walk.rs`), `old` the f64 fit's estimate on the same step (the shadow), `st` the walk's counts.
pub fn step<T: Real>(c: [&[T]; 3], l: f64, old: Option<Est>, g: &Gate, st: &mut LocStats, shadow: &mut Hist) -> Option<Est> {
    let ahead = |e: &Est| (e.side as f64 == g.dir) != g.flip;
    let old_ahead = old.filter(ahead);
    shadow.feed(old_ahead.map(|e| (g.x + g.dir * e.r, e.r)));
    if st.old_xstar.is_none() {
        if let Some((xs, bar)) = shadow.fired(g.bound) {
            st.old_xstar = Some(xs);
            st.old_bar = bar;
            st.old_step = g.steps;
        }
    }
    if !(g.h.is_finite() && g.h <= GATE * g.bound) {
        return None;
    }
    st.calls += 1;
    let mut best: Option<(Est, f64, f64)> = None;
    let (mut pairs, mut unres) = (0, 0);
    for (f, cf) in c.iter().enumerate() {
        match read_field(cf, f, l) {
            Field::Real(e, rb, ab) => {
                if best.is_none_or(|b| rb < b.1) {
                    best = Some((e, rb, ab));
                }
            }
            Field::Pair => pairs += 1,
            Field::Unresolved => unres += 1,
            Field::Nothing => {}
        }
    }
    let real_ahead = match &best {
        Some((e, rb, ab)) => {
            if pairs > 0 {
                st.pair_vs_real += 1;
            }
            st.rel_bar = *rb;
            st.alpha_bar = *ab;
            if ahead(e) {
                st.real += 1;
                true
            } else {
                st.behind += 1;
                false
            }
        }
        None => {
            if pairs > 0 {
                st.pair += 1;
            } else if unres > 0 {
                st.unresolved += 1;
            } else {
                st.nothing += 1;
            }
            false
        }
    };
    if old_ahead.is_some() {
        st.old_real += 1;
        if !real_ahead {
            st.unres_old_real += 1;
        }
    }
    best.map(|b| b.0)
}

#[cfg(test)]
mod tests {
    use super::*;
    use multidouble::dd::DD;
    use multidouble::qd::QD;

    #[test]
    fn pow2_near_rounds_on_the_log_scale() {
        assert_eq!(pow2_near(1.0), Some(1.0));
        assert_eq!(pow2_near(1.4), Some(1.0));
        assert_eq!(pow2_near(1.5), Some(2.0));
        assert_eq!(pow2_near(0.3), Some(0.25));
        assert_eq!(pow2_near(0.0), None);
        assert_eq!(pow2_near(f64::NAN), None);
        assert_eq!(pow2_near(1e-310), None);
    }

    /// the scaling is exact: scaled by rho^k and back, a double-double and a quad-double series come back bit for bit
    #[test]
    fn scaling_is_exact() {
        fn check<T: Real + PartialEq + std::fmt::Debug>() {
            let c: Vec<T> = (0..=30).map(|k| T::from_f64(1.0) / T::from_f64(3.0 + k as f64) * T::from_f64(0.7f64.powi(k))).collect();
            let (s, rho) = scaled(&c).unwrap();
            let mut p = 1.0;
            for (a, b) in c.iter().zip(&s) {
                assert_eq!(*b * T::from_f64(1.0 / p), *a, "scale {p}");
                p *= rho;
            }
        }
        check::<DD>();
        check::<QD>();
    }
}
