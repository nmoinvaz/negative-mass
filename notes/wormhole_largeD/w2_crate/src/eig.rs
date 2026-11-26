//! The mismatch determinant, Newton on the eigenvalue, and the contraction length.
use crate::flow::{ell_inf, flow_inward, throat, throat_series, FlowStats, LaneCfg};
use crate::num::Num;
use crate::sys::Sys;

pub struct Mismatch<T> {
    pub f: T,
    pub fp: T,
    pub tail: f64,
    pub stats: FlowStats,
    pub vin: [T; 4],
    pub vser: [T; 4],
}

/// M(Om) = Y_in w_ser - w_in Y_ser at l0 = 1 and its derivative in Om, the inward flow from (Y, w) = (1, 1) at l_inf
pub fn mismatch<T: Num>(sys: &Sys<T>, om: T, l_inf: f64, cfg: &LaneCfg) -> Mismatch<T> {
    let (vs, tail) = throat(sys, om, cfg.k0);
    let (vi, stats) = flow_inward(sys, om, T::from_f64(l_inf), T::one(), [T::one(), T::one(), T::zero(), T::zero()], cfg);
    let f = vi[0] * vs[1] - vi[1] * vs[0];
    let fp = vi[2] * vs[1] + vi[0] * vs[3] - vi[3] * vs[0] - vi[1] * vs[2];
    Mismatch { f, fp, tail, stats, vin: vi, vser: vs }
}

pub struct Solved<T> {
    pub om: T,
    pub iters: usize,
    /// log10 of |step| / max(|Om|, 1) at each iteration
    pub trace: Vec<f64>,
    pub converged: bool,
    pub l_inf: f64,
    pub tail: f64,
    pub stats: FlowStats,
    /// Om at each iteration, as f64
    pub omtrace: Vec<f64>,
    /// log10 of |M| / (max|v_in| max|v_ser|) at the last iterate
    pub m_norm: f64,
}

/// Newton from om0 until two consecutive steps are below 10^-d relative
pub fn newton<T: Num>(sys: &Sys<T>, om0: T, t: f64, cfg: &LaneCfg, d: f64, max_iter: usize) -> Solved<T> {
    let eps = sys.eps.to_f64();
    let l_inf = ell_inf(eps, om0.to_f64(), 1.0, t, sys.planted);
    let mut om = om0;
    let mut trace = vec![];
    let mut omtrace = vec![];
    let mut small = 0;
    let (mut tail, mut stats) = (0.0, FlowStats::default());
    let mut converged = false;
    for _ in 0..max_iter {
        let m = mismatch(sys, om, l_inf, cfg);
        tail = m.tail;
        stats = m.stats;
        let dom = -(m.f / m.fp);
        om = om + dom;
        omtrace.push(om.to_f64());
        let rel = dom.abs().to_f64() / om.abs().to_f64().max(1.0);
        trace.push(rel.max(1e-300).log10());
        if rel < 10f64.powf(-d) {
            small += 1;
            if small >= 2 {
                converged = true;
                break;
            }
        } else {
            small = 0;
        }
    }
    // |M| at the last iterate against the size of the two vectors
    let m = mismatch(sys, om, l_inf, cfg);
    let nv = |v: &[T; 4]| v[0].abs().to_f64().max(v[1].abs().to_f64());
    let m_norm = (m.f.abs().to_f64() / (nv(&m.vin) * nv(&m.vser))).max(1e-300).log10();
    Solved { om, iters: trace.len(), trace, converged, l_inf, tail, stats, omtrace, m_norm }
}

/// Sign changes of Y along the mode at Om: the throat series on 0 <= l <= 1 (21 points), then the inward flow's step
/// ends from l = 1 outward, the flow's overall sign fixed by continuity at l = 1
pub fn nodes<T: Num>(sys: &Sys<T>, om: T, l_inf: f64, cfg: &LaneCfg) -> usize {
    let ser = throat_series(sys, om, cfg.k0);
    let y: Vec<f64> = ser[0].iter().map(|c| c.to_f64()).collect();
    let at = |l: f64| y.iter().rev().fold(0.0, |a, c| a * l + c);
    let mut seq: Vec<f64> = (0..=20).map(|j| at(0.05 * j as f64)).collect();
    let (vi, st) = flow_inward(sys, om, T::from_f64(l_inf), T::one(), [T::one(), T::one(), T::zero(), T::zero()], cfg);
    let sgn = if vi[0].to_f64() * seq[20] >= 0.0 { 1.0 } else { -1.0 };
    for &(_, yv) in st.trace.iter().rev().skip(1) {
        seq.push(sgn * yv);
    }
    seq.windows(2).filter(|w| w[0] * w[1] < 0.0).count()
}
