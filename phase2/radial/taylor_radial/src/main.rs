//! taylor_radial: the radial boson-star modes and the phantom family by Taylor integration (negative-mass #194).
//! Usage: taylor_radial <cmd> key=value ...   (see ../../taylor_radial.py for the calls it makes)
mod bvp;
mod linalg;
mod model;
mod mutate;
mod num;
mod plane;
mod planted;
mod radial;
mod rgen;
mod roots;
mod ser;
mod star_ref;
mod table;

use multidouble::dd::DD;
use num::Num;
use std::collections::HashMap;

pub struct Args(HashMap<String, String>);

impl Args {
    fn get(&self, k: &str) -> Option<&str> {
        self.0.get(k).map(|s| s.as_str())
    }
    pub fn f(&self, k: &str, d: f64) -> f64 {
        self.get(k).map(|s| s.parse().unwrap_or_else(|_| panic!("bad {k}={s}"))).unwrap_or(d)
    }
    pub fn u(&self, k: &str, d: usize) -> usize {
        self.get(k).map(|s| s.parse().unwrap_or_else(|_| panic!("bad {k}={s}"))).unwrap_or(d)
    }
    pub fn s<'a>(&'a self, k: &str, d: &'a str) -> &'a str {
        self.get(k).unwrap_or(d)
    }
    pub fn list(&self, k: &str) -> Vec<f64> {
        self.get(k).map(|s| s.split(',').map(|x| x.parse().unwrap()).collect()).unwrap_or_default()
    }
}

/// Far-field distance in units of 1/kappa_min: the contamination of the decaying plane is e^{-2 kappa (R - r_m)},
/// below the lane's roundoff.
pub fn far_factor<T: Num>() -> f64 {
    match T::NAME {
        "f64" => 19.0,
        "dd" => 40.0,
        "qd" => 74.0,
        _ => 150.0,
    }
}

/// Taylor order of the lane: the order whose step is long at the lane's tolerance (multishoot's cost model prefers
/// about 22 at 1e-14 and 42 at 1e-30; the table's series carry 30 in f64 and dd).
pub fn default_order<T: Num>() -> usize {
    match T::NAME {
        "qd" => 56,
        _ => 30,
    }
}

fn lane_dispatch(a: &Args, cmd: &str) {
    match a.s("lane", "dd") {
        "f64" => run::<f64>(a, cmd),
        "dd" => run::<DD>(a, cmd),
        "qd" => run::<multidouble::qd::QD>(a, cmd),
        x => panic!("lane {x}"),
    }
}

fn run<T: Num>(a: &Args, cmd: &str) {
    match cmd {
        "planted" => cmd_planted::<T>(a),
        _ => panic!("unknown command {cmd}"),
    }
}

fn cmd_planted<T: Num>(a: &Args) {
    use plane::Ctl;
    let order = a.u("order", default_order::<T>());
    let ctl = Ctl { order, tol: a.f("tol", T::EPS * 100.0), hmax: a.f("hmax", 2.0), kj: a.u("kj", order) };
    let sign = if mutate::on("plant_sign") { -T::one() } else { T::one() };
    let (s_lo, s_hi, n) = (a.f("slo", -4.0), a.f("shi", 1.9), a.u("n", 60));
    let r_m = a.f("rm", 3.0);
    let kmin = (4.0 - 2.0 * s_hi).sqrt().min((3.0 - s_hi).sqrt());
    let r_far = a.f("rfar", r_m + far_factor::<T>() / kmin);
    let t0 = std::time::Instant::now();
    let mut f = |s: T| planted::d_planted::<T>(sign, s, r_m, r_far, &ctl);
    let grid: Vec<T> = (0..n).map(|i| T::from_f64(s_lo + (s_hi - s_lo) * i as f64 / (n - 1) as f64)).collect();
    let vals: Vec<T> = grid.iter().map(|&s| f(s)).collect();
    let mut found = vec![];
    for i in 0..n - 1 {
        if (vals[i].to_f64() < 0.0) != (vals[i + 1].to_f64() < 0.0) {
            let r = roots::find::<T>(&mut f, grid[i], grid[i + 1], vals[i], vals[i + 1]);
            found.push(r);
        }
    }
    let ex: Vec<(T, usize, usize)> = planted::exact::<T>();
    println!(
        "{{\"cmd\":\"planted\",\"lane\":\"{}\",\"rm\":{},\"rfar\":{},\"order\":{},\"tol\":{:e},\"grid\":{},\"nroots\":{},\"nexact\":{},\"wall\":{:.3}}}",
        T::NAME,
        r_m,
        r_far,
        order,
        ctl.tol,
        n,
        found.len(),
        ex.len(),
        t0.elapsed().as_secs_f64()
    );
    for (i, r) in found.iter().enumerate() {
        let (e, blk, nn) = ex.get(i).copied().unwrap_or((T::zero(), 9, 9));
        let rel = ((r.x - e) / e).abs();
        println!(
            "{{\"cmd\":\"planted_root\",\"i\":{},\"s\":\"{}\",\"exact\":\"{}\",\"rel\":{:e},\"block\":{},\"n\":{},\"evals\":{},\"d_at_root\":{:e}}}",
            i,
            r.x.sci(40),
            e.sci(40),
            rel.to_f64(),
            blk,
            nn,
            r.evals,
            r.fx.to_f64()
        );
    }
}

fn cmd_star(a: &Args) {
    let eps = a.f("eps", 1.0);
    let omega = DD::from_f64(a.f("omega", 0.75));
    let mut cfg = bvp::StarCfg::new(eps);
    cfg.dmax = a.f("dmax", 8.0);
    let t0 = std::time::Instant::now();
    let sol = bvp::solve(&cfg, omega, true).unwrap();
    let rep = sol.report.as_ref().unwrap();
    if std::env::var("TR_VERBOSE").is_ok() {
        eprintln!("worst rows (block, field, lane, |F| scaled by nothing): {:?}", bvp::worst_rows(&sol, 8));
    }
    println!(
        "{{\"cmd\":\"star\",\"eps\":{},\"omega\":\"{}\",\"nodes\":{},\"pc\":\"{}\",\"zc\":\"{}\",\"mass\":\"{}\",\"converged\":{},\"fnorm\":{:?},\"seconds\":{:.2}}}",
        eps,
        sol.omega().sci(32),
        sol.sh.k(),
        sol.pc().sci(32),
        sol.zc().sci(32),
        sol.mass().sci(32),
        rep.converged,
        rep.fnorm,
        t0.elapsed().as_secs_f64()
    );
}

/// A star at omega, its table, and D(s2) on a grid of s2; roots of D by sign changes
fn cmd_modes(a: &Args) {
    use multidouble::dd::DD as D;
    let eps = a.f("eps", 1.0);
    let omega = D::from_f64(a.f("omega", 0.94));
    let mut cfg = bvp::StarCfg::new(eps);
    cfg.dmax = a.f("dmax", 8.0);
    let t0 = std::time::Instant::now();
    let sol = bvp::solve(&cfg, omega, false).unwrap();
    let rep = sol.report.as_ref().unwrap();
    assert!(rep.converged, "star not converged: {:?}", rep.fnorm);
    let t_star = t0.elapsed().as_secs_f64();
    let steps = table::steps_dd(&sol);
    let kk = sol.k().to_f64();
    let r_t = sol.sh.nodes[sol.sh.k() - 1];
    // the exterior chain starts where P falls below the tail value: the last BVP step whose P(start) is below p_tail
    let r_cut = r_t - 40.0 / kk;
    let s2_hi = a.f("shi", 0.99 * (1.0 - omega.to_f64()).powi(2));
    let om = omega.to_f64();
    let kmin = (1.0 - (om + s2_hi.max(0.0).sqrt()).powi(2)).max(1e-6).sqrt().min(kk);
    let r_far = a.f("rfar", r_t.max(r_cut + 40.0 / kmin).max(r_cut + 40.0 / kk));
    let mass = sol.mass();
    let ext = radial::ext_steps::<D>(mass, sol.omega(), r_cut, r_far, 1e-32);
    let mut all: Vec<table::StarStep<D>> = steps.into_iter().filter(|s| s.hi <= r_cut + 1e-9).collect();
    // the BVP step containing r_cut is cut: keep steps wholly below, add one clipped
    let last_hi = all.last().map(|s| s.hi).unwrap_or(0.0);
    all.extend(ext.into_iter());
    let _ = last_hi;
    let model = table::Model { omega: sol.omega(), epsilon: eps };
    let tab = table::build_table(&all, &model, cfg.rs, r_far);
    let t_tab = t0.elapsed().as_secs_f64() - t_star;
    println!(
        "{{\"cmd\":\"modes_setup\",\"omega\":{},\"eps\":{},\"r_cut\":{:.3},\"r_far\":{:.3},\"table_steps\":{},\"linear_defect\":{:e},\"t_star\":{:.2},\"t_table\":{:.2}}}",
        om,
        eps,
        r_cut,
        r_far,
        tab.x0.len(),
        tab.linear_defect,
        t_star,
        t_tab
    );
    let lane = a.s("lane", "dd");
    let order = a.u("order", 30);
    let (lo, hi, n) = (a.f("slo", -0.02), a.f("shi", s2_hi), a.u("n", 20));
    match lane {
        "dd" => scan::<D>(a, &tab, &all[0], &sol, order, lo, hi, n, r_far),
        "f64" => scan::<f64>(a, &tab, &all[0], &sol, order, lo, hi, n, r_far),
        _ => panic!("lane"),
    }
}

fn scan<T: Num>(a: &Args, tab: &table::JTable<DD>, first: &table::StarStep<DD>, sol: &bvp::StarSol, order: usize, lo: f64, hi: f64, n: usize, r_far: f64) {
    let ctl = plane::Ctl { order, tol: a.f("tol", T::EPS * 100.0), hmax: a.f("hmax", 4.0), kj: order };
    let rp = radial::Radial::<T>::new(tab.to_lane::<T>(), first, 1e-8, T::from_dd(sol.omega()), sol.pc().to_f64(), ctl);
    let far_f = far_factor::<T>();
    let t0 = std::time::Instant::now();
    let grid: Vec<T> = (0..n).map(|i| T::from_f64(lo + (hi - lo) * i as f64 / (n - 1) as f64)).collect();
    let om = rp.omega.to_f64();
    let rfar_of = |s2: f64| {
        let kap = (1.0 - (om + s2.max(0.0).sqrt()).powi(2)).max(1e-6).sqrt().min((1.0 - om * om).sqrt());
        let kap2 = (1.0 - (om - s2.max(0.0).sqrt()).powi(2)).max(1e-6).sqrt();
        let km = kap.min(kap2);
        let _ = far_f;
        r_far.min(rp.r_m + far_f / km).max(rp.r_m + 5.0)
    };
    println!("{{\"cmd\":\"modes_rm\",\"r_m\":{},\"r_max\":{}}}", rp.r_m, rp.r_max);
    let mut f = |s: T| rp.d(s, rfar_of(s.to_f64()));
    let vals: Vec<T> = grid.iter().map(|&s| f(s)).collect();
    for (g, v) in grid.iter().zip(&vals) {
        println!("{{\"cmd\":\"D\",\"s2\":{:e},\"D\":{:e},\"rfar\":{:.1}}}", g.to_f64(), v.to_f64(), rfar_of(g.to_f64()));
    }
    for i in 0..n - 1 {
        if (vals[i].to_f64() < 0.0) != (vals[i + 1].to_f64() < 0.0) {
            let r = roots::find::<T>(&mut f, grid[i], grid[i + 1], vals[i], vals[i + 1]);
            println!("{{\"cmd\":\"root\",\"s2\":\"{}\",\"evals\":{},\"d\":{:e}}}", r.x.sci(30), r.evals, r.fx.to_f64());
        }
    }
    println!("{{\"cmd\":\"scan_done\",\"lane\":\"{}\",\"wall\":{:.2}}}", T::NAME, t0.elapsed().as_secs_f64());
}

fn main() {
    let argv: Vec<String> = std::env::args().collect();
    let cmd = argv.get(1).cloned().unwrap_or_default();
    let mut m = HashMap::new();
    for x in &argv[2..] {
        if let Some((k, v)) = x.split_once('=') {
            if k == "mutate" {
                mutate::set(v);
            }
            m.insert(k.to_string(), v.to_string());
        }
    }
    let a = Args(m);
    if cmd == "star" {
        cmd_star(&a);
        return;
    }
    if cmd == "modes" {
        cmd_modes(&a);
        return;
    }
    lane_dispatch(&a, &cmd);
}
