//! D at complex s (negative-mass #188): the commands `dc` (a list of s) and `dserve` (s read from stdin, one JSON line out per
//! line in), for the argument-principle counts of winding2.py. One star and one J table per process; each evaluation is
//! `polar::d_eval_s` on the complex lane `Cx<T>`.
//!
//! Options: lane (f64, dd, qd), l, pc, zc (the star), rm, r0, order, tol, hmax, kj as for `d`; `rsfix` the fixed start radius of the
//! scalar plane (winding2's R_s; default r_m + 40 / min Re k over the two k± at the s asked); `zfac` as winding2's --zfac.
//! The output D is not stripped of winding2's phase factor exp(-i Im[...]): the counter applies it in double precision.
use crate::complex::Cx;
use crate::num::Num;
use crate::plane::Ctl;
use crate::polar::{self, Cfg};
use crate::table;
use crate::{star, Args};
use multidouble::dd::DD;
use std::io::{BufRead, Write};
use std::time::Instant;

pub struct Setup<T: Num> {
    pub cfg: Cfg,
    pub tab: table::JTable<Cx<T>>,
    pub mext: Cx<T>,
    pub t_star: f64,
    pub t_table: f64,
    pub steps_star: usize,
    pub steps_table: usize,
}

/// the star and the J table of one l, as `run_d` builds them, on the complex lane
pub fn setup<T: Num>(a: &Args) -> Setup<T> {
    let m = crate::model(a);
    let l = a.u("l", 2);
    let r0 = a.f("r0", 1e-3);
    let rm = a.f("rm", 24.10980364018201);
    let order = a.u("order", 30);
    let ctl = Ctl { order, tol: a.f("tol", T::EPS * 100.0), hmax: a.f("hmax", 1.0), kj: a.u("kj", order) };
    let cfg = Cfg { l, r0, rm, ctl, model: m };
    let sp = crate::star_pars::<DD>(&m);
    let scfg = crate::solve_cfg(a);
    let (pc, zc) = crate::parse_star(a);
    let t = Instant::now();
    let st = star::build_star(pc, zc, &sp, &scfg, rm + 1.0);
    let t_star = t.elapsed().as_secs_f64();
    let t = Instant::now();
    let tab = table::build_table(&st, l, &m, r0 * 0.5, rm + 0.5);
    let t_table = t.elapsed().as_secs_f64();
    let tabl = tab.to_lane::<Cx<T>>();
    let mext = Cx::<T>::from_dd(st.state(rm)[2]);
    Setup { cfg, tab: tabl, mext, t_star, t_table, steps_star: st.steps.len(), steps_table: tab.x0.len() }
}

/// the start radius of the scalar plane at s = sre + i sim: r_m + 40 / min Re k±, k± = sqrt(1 + (s -+ i w)^2).
/// At sim = 0 the operations are those of `polar::scalar_start`, so the radius is the same bits.
pub fn scalar_start_c(rm: f64, sre: f64, sim: f64, omega: f64) -> f64 {
    let mut kmin = f64::INFINITY;
    for sign in [1.0, -1.0] {
        let d = sim - sign * omega;
        let re = 1.0 + sre * sre - d * d;
        let im = 2.0 * sre * d;
        let kr = (((re * re + im * im).sqrt() + re) / 2.0).sqrt();
        kmin = kmin.min(kr);
    }
    rm + 40.0 / kmin
}

pub fn eval_line<T: Num>(a: &Args, su: &Setup<T>, sre: f64, sim: f64) -> String {
    let cfg = &su.cfg;
    let zfac = a.f("zfac", 1.0);
    let s = Cx::new(T::from_f64(sre), T::from_f64(sim));
    let sabs = sre.hypot(sim);
    let r_s = match a.f("rsfix", f64::NAN) {
        x if x.is_nan() => scalar_start_c(cfg.rm, sre, sim, cfg.model.omega),
        x => x,
    };
    let r_far = if zfac == 1.0 { polar::far_start(cfg.rm, sabs, cfg.l) } else { cfg.rm + zfac * (40.0 + 20.0 * cfg.l as f64) / sabs + 50.0 };
    let t = Instant::now();
    let o = polar::d_eval_s(cfg, &su.tab, su.mext, s, r_s, r_far);
    let wall = t.elapsed().as_secs_f64();
    format!(
        "{{\"cmd\":\"dc\",\"lane\":\"{}\",\"l\":{},\"re\":{:?},\"im\":{:?},\"d_re\":\"{}\",\"d_im\":\"{}\",\"d_f64\":[{:e},{:e}],\"q_f64\":[{:e},{:e}],\"steps_y\":{},\"steps_q\":{},\"steps_s\":{},\"t_y\":{:.4},\"t_q\":{:.4},\"t_s\":{:.4},\"wall\":{:.4},\"r_s\":{:?},\"r_far\":{:?}}}",
        T::NAME,
        cfg.l,
        sre,
        sim,
        o.d.re.sci(130),
        o.d.im.sci(130),
        o.d.re.to_f64(),
        o.d.im.to_f64(),
        o.q.re.to_f64(),
        o.q.im.to_f64(),
        o.steps_y.steps,
        o.steps_q,
        o.steps_s.steps,
        o.t_y,
        o.t_q,
        o.t_s,
        wall,
        o.r_s,
        o.r_far
    )
}

fn setup_line<T: Num>(a: &Args, su: &Setup<T>) -> String {
    format!(
        "{{\"cmd\":\"dc_setup\",\"lane\":\"{}\",\"l\":{},\"rm\":{:?},\"r0\":{:e},\"order\":{},\"tol\":{:e},\"mext\":\"{}\",\"star_steps\":{},\"table_steps\":{},\"t_star\":{:.4},\"t_table\":{:.4},\"zfac\":{},\"rs\":\"{}\"}}",
        T::NAME,
        su.cfg.l,
        su.cfg.rm,
        su.cfg.r0,
        su.cfg.ctl.order,
        su.cfg.ctl.tol,
        su.mext.re.sci(130),
        su.steps_star,
        su.steps_table,
        su.t_star,
        su.t_table,
        a.f("zfac", 1.0),
        a.s("rsfix", "per s")
    )
}

fn dc_go<T: Num>(a: &Args) {
    let su = setup::<T>(a);
    println!("{}", setup_line(a, &su));
    let (re, im) = (a.list("re"), a.list("im"));
    assert_eq!(re.len(), im.len(), "re and im lists differ in length");
    for (x, y) in re.iter().zip(im.iter()) {
        println!("{}", eval_line(a, &su, *x, *y));
    }
}

fn serve_go<T: Num>(a: &Args) {
    let su = setup::<T>(a);
    println!("{}", setup_line(a, &su));
    std::io::stdout().flush().unwrap();
    for line in std::io::stdin().lock().lines() {
        let line = line.unwrap();
        let mut it = line.split_whitespace();
        let (Some(x), Some(y)) = (it.next(), it.next()) else { continue };
        let (x, y): (f64, f64) = (x.parse().unwrap(), y.parse().unwrap());
        println!("{}", eval_line(a, &su, x, y));
        std::io::stdout().flush().unwrap();
    }
}

pub fn dc(a: &Args) {
    match a.s("lane", "dd") {
        "f64" => dc_go::<f64>(a),
        "dd" => dc_go::<DD>(a),
        "qd" => dc_go::<multidouble::qd::QD>(a),
        x => panic!("lane {x}"),
    }
}

pub fn dserve(a: &Args) {
    match a.s("lane", "dd") {
        "f64" => serve_go::<f64>(a),
        "dd" => serve_go::<DD>(a),
        "qd" => serve_go::<multidouble::qd::QD>(a),
        x => panic!("lane {x}"),
    }
}
