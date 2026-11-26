//! taylor_radial: the radial boson-star modes and the phantom family by Taylor integration (negative-mass #194).
//! Usage: taylor_radial <cmd> key=value ...   (see ../../taylor_radial.py for the calls it makes)
mod linalg;
mod mutate;
mod num;
mod plane;
mod planted;
mod rgen;
mod roots;
mod ser;

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
    lane_dispatch(&a, &cmd);
}
