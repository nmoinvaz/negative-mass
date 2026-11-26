//! w2_crate: W^2 a^2 of the Ellis wormhole's growing mode (E2 of master.py) at a given D, in f64, double-, quad- and
//! octo-double (#195).
//!   w2_crate w2 --lane=qd --D=50,100 [--W2=<guess>] [--T=<contraction>] [--mutate=NAME]
//!   w2_crate planted --lane=qd --D=5,50,1000,10000 [--mutate=NAME]
//!   w2_crate series --lane=qd --D=10000       (the throat series against the closed forms Y0, Y1, Y2 of README.md)
mod eig;
mod flow;
mod mutate;
mod num;
mod sys;

use eig::{newton, nodes};
use flow::{ell_inf, flow_inward, lane_cfg, throat, LaneCfg};
use multidouble::dd::DD;
use multidouble::od::OD;
use multidouble::qd::QD;
use num::Num;
use std::time::Instant;
use sys::Sys;

struct Args {
    lane: String,
    ds: Vec<u64>,
    w2: Vec<f64>,
    t: Option<f64>,
    off: (f64, f64),
    shift: f64,
    scan: (f64, f64, usize),
    real: bool,
}

fn lg(x: f64) -> f64 {
    if x <= 0.0 {
        -300.0
    } else {
        x.log10()
    }
}

fn rel<T: Num>(a: T, b: T) -> f64 {
    lg(((a - b) / b).abs().to_f64())
}

fn arr(v: &[f64]) -> String {
    format!("[{}]", v.iter().map(|x| format!("{x:.2}")).collect::<Vec<_>>().join(","))
}

fn arr12(v: &[f64]) -> String {
    format!("[{}]", v.iter().take(12).map(|x| format!("{x:.6}")).collect::<Vec<_>>().join(","))
}

fn run_w2<T: Num>(a: &Args) {
    let cfg = lane_cfg(&a.lane);
    for (i, &d) in a.ds.iter().enumerate() {
        let t0 = Instant::now();
        let sys = Sys::<T>::new(d, false, cfg.kt);
        let s = 2.0 * (d as f64 - 3.0);
        let eps = 1.0 / s;
        let om0 = match a.w2.get(i) {
            Some(&w) if w > 0.0 => w - s,
            _ => -4.0 * 2f64.ln() + 10.4233603787385 * eps + 73.0 * eps * eps,
        };
        let t = a.t.unwrap_or(cfg.t);
        let r = newton(&sys, T::from_f64(om0 + a.shift), t, &cfg, cfg.d_lane, 40);
        let w2 = sys.s + r.om;
        println!(
            "RESULT {{\"kind\":\"w2\",\"lane\":\"{}\",\"kt_ms\":{},\"D\":{},\"W2\":\"{}\",\"omega\":\"{}\",\"converged\":{},\"iters\":{},\"trace\":{},\"l_inf\":{},\"steps\":{},\"retries\":{},\"hmin\":{:.3},\"hmax\":{:.3},\"tail\":{:.2},\"T\":{},\"secs\":{:.3}}}",
            cfg.name,
            multishoot::Plan::optimal_kt(cfg.tol),
            d,
            w2.sci(cfg.digits as usize + 2),
            r.om.sci(cfg.digits as usize + 2),
            r.converged,
            r.iters,
            arr(&r.trace),
            r.l_inf,
            r.stats.steps,
            r.stats.retries,
            r.stats.hmin,
            r.stats.hmax,
            lg(r.tail),
            t,
            t0.elapsed().as_secs_f64()
        );
    }
}

fn run_planted<T: Num>(a: &Args) {
    let cfg = lane_cfg(&a.lane);
    for &d in &a.ds {
        let t0 = Instant::now();
        let sys = Sys::<T>::new(d, true, cfg.kt);
        let om_s = sys.om_star;
        let eps = sys.eps.to_f64();
        let t = a.t.unwrap_or(cfg.t);
        // S1: the throat series at Om*, summed at l0 = 1, against the closed form
        let (vs, tail) = throat(&sys, om_s, cfg.k0);
        let (ys, ws) = sys.exact(T::one());
        let s1y = rel(vs[0], ys);
        let s1w = rel(vs[1], ws);
        // S2: the inward flow from the exact state at l_inf
        let l_inf = ell_inf(eps, om_s.to_f64(), 1.0, t, true);
        let (yi, wi) = sys.exact(T::from_f64(l_inf));
        let (v2, st) = flow_inward(&sys, om_s, T::from_f64(l_inf), T::one(), [yi, wi, T::zero(), T::zero()], &cfg);
        let s2y = rel(v2[0], ys);
        let s2w = rel(v2[1], ws);
        // S3: the contraction from (1, 1): the direction w/Y against the closed form
        let (v3, _) = flow_inward(&sys, om_s, T::from_f64(l_inf), T::one(), [T::one(), T::one(), T::zero(), T::zero()], &cfg);
        let s3 = rel(v3[1] / v3[0], ws / ys);
        // S4: Newton from two starts; where each ends (Om, |M| against the vectors, nodes of Y)
        let s4 = |off: f64| {
            let r = newton(&sys, om_s + T::from_f64(off), t, &cfg, cfg.d_lane, 60);
            let nd = nodes(&sys, r.om, r.l_inf, &cfg);
            (rel(r.om, om_s), r.iters, r.converged, r.om.to_f64(), r.m_norm, nd, r.omtrace)
        };
        let (s4a, ia, ca, oa, ma, na, _) = s4(a.off.0);
        let (s4b, ib, cb, ob, mb, nb, tb) = s4(a.off.1);
        // S5: the shorter contraction, and its bound 100 e^{-T'} + 10^{-d} from the declared T
        let t5 = t / 2.0 + 20.0;
        let t5_decl = cfg.t_decl / 2.0 + 20.0;
        let r5 = newton(&sys, om_s + T::from_f64(a.off.0), t5, &cfg, cfg.d_lane, 60);
        let s5 = rel(r5.om, om_s);
        let bound5 = lg(100.0 * (-t5_decl).exp() + 10f64.powf(-cfg.d_lane));
        println!(
            "RESULT {{\"kind\":\"planted\",\"lane\":\"{}\",\"D\":{},\"S1y\":{:.2},\"S1w\":{:.2},\"S1tail\":{:.2},\"S2y\":{:.2},\"S2w\":{:.2},\"S3\":{:.2},\"S4a\":{:.2},\"S4b\":{:.2},\"S4iters\":[{},{}],\"S4conv\":[{},{}],\"S4om\":[{:.12},{:.12}],\"S4M\":[{:.1},{:.1}],\"S4nodes\":[{},{}],\"S4btrace\":{},\"S5\":{:.2},\"S5bound\":{:.2},\"T5\":{:.2},\"l_inf\":{},\"steps\":{},\"retries\":{},\"secs\":{:.3}}}",
            cfg.name, d, s1y, s1w, lg(tail), s2y, s2w, s3, s4a, s4b, ia, ib, ca, cb, oa, ob, ma, mb, na, nb, arr12(&tb), s5, bound5, t5_decl, l_inf, st.steps, st.retries,
            t0.elapsed().as_secs_f64()
        );
    }
}

/// M(Om) on a grid: Om, the normalized M (signed) and its derivative's sign, and the nodes of Y at that Om
fn run_scan<T: Num>(a: &Args) {
    let cfg = lane_cfg(&a.lane);
    for &d in &a.ds {
        let sys = Sys::<T>::new(d, !a.real, cfg.kt);
        let eps = sys.eps.to_f64();
        let t = a.t.unwrap_or(cfg.t);
        for i in 0..a.scan.2 {
            let om = a.scan.0 + (a.scan.1 - a.scan.0) * i as f64 / (a.scan.2 - 1).max(1) as f64;
            let l_inf = ell_inf(eps, om, 1.0, t, !a.real);
            let m = eig::mismatch(&sys, T::from_f64(om), l_inf, &cfg);
            let nv = |v: &[T; 4]| v[0].abs().to_f64().max(v[1].abs().to_f64());
            let mn = m.f.to_f64() / (nv(&m.vin) * nv(&m.vser));
            let nd = nodes(&sys, T::from_f64(om), l_inf, &cfg);
            println!(
                "RESULT {{\"kind\":\"scan\",\"lane\":\"{}\",\"D\":{},\"om\":{:.6},\"M\":{:e},\"dM\":{:e},\"nodes\":{}}}",
                cfg.name, d, om, mn, m.fp.to_f64() / (nv(&m.vin) * nv(&m.vser)), nd
            );
        }
    }
}

fn run_series<T: Num>(a: &Args) {
    let cfg = lane_cfg(&a.lane);
    for &d in &a.ds {
        let sys = Sys::<T>::new(d, false, cfg.kt);
        let eps = 1.0 / (2.0 * (d as f64 - 3.0));
        let o0 = -4.0 * 2f64.ln();
        let om = T::from_f64(o0 + 10.4233603787385 * eps + 73.0 * eps * eps);
        let (v, tail) = throat(&sys, om, cfg.k0);
        // Y = 1 + eps Y1 + eps^2 Y2, Y1 = -l, Y2 = -2 Om0 l - 3/2 l^2 at l = 1
        let ycf = 1.0 - eps + eps * eps * (-2.0 * o0 - 1.5);
        let diff = v[0].to_f64() - ycf;
        println!(
            "RESULT {{\"kind\":\"series\",\"lane\":\"{}\",\"D\":{},\"Y\":{:.15},\"Ycf\":{:.15},\"diff\":{:e},\"diff_over_eps3\":{:.4},\"tail\":{:.2}}}",
            cfg.name, d, v[0].to_f64(), ycf, diff, diff / (eps * eps * eps), lg(tail)
        );
    }
}

fn main() {
    let argv: Vec<String> = std::env::args().collect();
    let mode = argv.get(1).cloned().unwrap_or_default();
    let mut a = Args { lane: "qd".into(), ds: vec![], w2: vec![], t: None, off: (0.3, -0.4), shift: 0.0, scan: (0.0, 1.0, 11), real: false };
    for x in &argv[2..] {
        if let Some(v) = x.strip_prefix("--lane=") {
            a.lane = v.into();
        } else if let Some(v) = x.strip_prefix("--D=") {
            a.ds = v.split(',').map(|t| t.parse().unwrap()).collect();
        } else if let Some(v) = x.strip_prefix("--W2=") {
            a.w2 = v.split(',').map(|t| t.parse().unwrap()).collect();
        } else if let Some(v) = x.strip_prefix("--T=") {
            a.t = Some(v.parse().unwrap());
        } else if let Some(v) = x.strip_prefix("--scan=") {
            let p: Vec<&str> = v.split(',').collect();
            a.scan = (p[0].parse().unwrap(), p[1].parse().unwrap(), p[2].parse().unwrap());
        } else if x == "--real" {
            a.real = true;
        } else if let Some(v) = x.strip_prefix("--shift=") {
            a.shift = v.parse().unwrap();
        } else if let Some(v) = x.strip_prefix("--off=") {
            let p: Vec<f64> = v.split(',').map(|t| t.parse().unwrap()).collect();
            a.off = (p[0], p[1]);
        } else if let Some(v) = x.strip_prefix("--mutate=") {
            mutate::set(v);
        } else {
            panic!("argument {x}");
        }
    }
    macro_rules! go {
        ($f:ident) => {
            match a.lane.as_str() {
                "f64" => $f::<f64>(&a),
                "dd" => $f::<DD>(&a),
                "qd" => $f::<QD>(&a),
                "od" => $f::<OD>(&a),
                l => panic!("lane {l}"),
            }
        };
    }
    let _: Option<&LaneCfg> = None;
    match mode.as_str() {
        "w2" => go!(run_w2),
        "planted" => go!(run_planted),
        "series" => go!(run_series),
        "scan" => go!(run_scan),
        m => panic!("mode {m}"),
    }
}
