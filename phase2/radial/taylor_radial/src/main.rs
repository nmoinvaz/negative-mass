//! taylor_radial: the radial boson-star modes and the phantom family by Taylor integration (negative-mass #194).
//! Usage: taylor_radial <cmd> key=value ...   (see ../../taylor_radial.py for the calls it makes)
mod bvp;
mod cont;
mod family;
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

/// One row: the star at omega, M, N, and the roots of D on the declared grid
fn cmd_row(a: &Args) {
    use multidouble::dd::DD as D;
    let eps = a.f("eps", 1.0);
    let omega = D::from_f64(a.f("omega", 0.94));
    let mut cfg = bvp::StarCfg::new(eps);
    cfg.dmax = a.f("dmax", 8.0);
    let thr = (1.0 - omega.to_f64()).powi(2);
    let s2_hi = a.f("shi", 0.995 * thr);
    let prep = family::prepare_k(&cfg, omega, s2_hi).unwrap();
    println!(
        "{{\"cmd\":\"star\",\"eps\":{},\"omega\":\"{}\",\"k\":\"{}\",\"pc\":\"{}\",\"zc\":\"{}\",\"alpha_c\":\"{}\",\"M\":\"{}\",\"N\":\"{}\",\"r_cut\":{:.3},\"r_far\":{:.1},\"table_steps\":{},\"linear_defect\":{:e},\"nodes\":{},\"fnorm\":{:e},\"t_star\":{:.2},\"t_table\":{:.2}}}",
        eps,
        prep.omega.sci(32),
        prep.info.k.sci(32),
        prep.info.pc.sci(32),
        prep.info.zc.sci(32),
        prep.info.alpha_c.sci(32),
        prep.mass.sci(32),
        prep.noether.sci(32),
        prep.r_cut,
        prep.r_far,
        prep.tab.x0.len(),
        prep.tab.linear_defect,
        prep.info.nodes,
        prep.info.fnorm.last().unwrap(),
        prep.t_star,
        prep.t_table
    );
    let lane = a.s("lane", "dd");
    let order = a.u("order", default_order::<f64>());
    let neg_lo = a.f("neg_lo", 4.0);
    let g = if a.get("s2list").is_some() { a.list("s2list") } else { family::grid(neg_lo, thr, a.u("nneg", 60), a.u("npos", 60)) };
    match lane {
        "dd" => report_scan::<D>(a, &prep, &g, order),
        "f64" => report_scan::<f64>(a, &prep, &g, order),
        _ => panic!("lane"),
    }
}

fn report_scan<T: Num>(a: &Args, prep: &family::Prepared, g: &[f64], order: usize) {
    let t0 = std::time::Instant::now();
    let rp = prep.radial::<T>(order, a.f("tol", T::EPS * 100.0));
    let found = family::scan(&rp, g, far_factor::<T>() * a.f("rfar_scale", 1.0));
    if a.u("dump", 0) > 0 {
        for (s, v) in &found.vals {
            println!("{{\"cmd\":\"D\",\"s2\":{:e},\"D\":{:e}}}", s, v);
        }
    }
    let signs: Vec<i32> = found.vals.iter().map(|v| if v.1 < 0.0 { -1 } else { 1 }).collect();
    let changes = signs.windows(2).filter(|w| w[0] != w[1]).count();
    for (x, ev, d) in &found.roots {
        println!("{{\"cmd\":\"root\",\"lane\":\"{}\",\"s2\":\"{}\",\"evals\":{},\"d\":{:e}}}", T::NAME, x.sci(30), ev, d);
    }
    println!(
        "{{\"cmd\":\"scan_done\",\"lane\":\"{}\",\"r_m\":{},\"order\":{},\"grid\":{},\"sign_changes\":{},\"roots\":{},\"d_evals\":{},\"wall\":{:.2}}}",
        T::NAME,
        rp.r_m,
        order,
        g.len(),
        changes,
        found.roots.len(),
        found.evals,
        t0.elapsed().as_secs_f64()
    );
}

/// Plain continuation of the star in alpha_c from a start at omega: Newton on the alpha_c-pinned layout at each step
fn cmd_follow(a: &Args) {
    use multidouble::dd::DD as D;
    let eps = a.f("eps", -1.0);
    let omega0 = D::from_f64(a.f("omega", 0.7603));
    let mut cfg = bvp::StarCfg::new(eps);
    cfg.dmax = a.f("dmax", 8.0);
    let targets = a.list("alphas");
    let h = a.f("h", 0.05);
    let t0 = std::time::Instant::now();
    let sol = bvp::solve(&cfg, omega0, false).unwrap();
    let mut info = sol.info();
    // continuation in omega (steps dw) to omega_to, then in alpha_c
    let omega_to = a.f("omega_to", omega0.to_f64());
    let dw = a.f("dw", 0.01);
    let mut wcur = omega0.to_f64();
    while (wcur - omega_to).abs() > 1e-12 {
        let wn = if (wcur - omega_to).abs() <= dw { omega_to } else { wcur - dw * (wcur - omega_to).signum() };
        let kt = (D::ONE - D::from_f64(wn) * D::from_f64(wn)).sqrt();
        let s = bvp::solve_k_from_profile(&cfg, &info.profile(), kt, a.u("verbose", 0) > 0).unwrap();
        info = s.info();
        assert!(info.converged, "omega continuation failed at {wn}: {:?}", info.fnorm);
        wcur = wn;
    }
    let mut ac = info.alpha_c;
    println!("{{\"cmd\":\"start\",\"alpha_c\":\"{}\",\"omega\":\"{}\",\"t\":{:.1}}}", ac.sci(20), info.omega.sci(20), t0.elapsed().as_secs_f64());
    let mut hcur = h;
    for &tg in &targets {
        while (tg - ac.to_f64()).abs() > 1e-13 {
            let step = (tg - ac.to_f64()).clamp(-hcur, hcur);
            let next = if (tg - ac.to_f64()).abs() <= hcur { D::from_f64(tg) } else { ac + D::from_f64(step) };
            let t1 = std::time::Instant::now();
            let prof = info.profile();
            let ni = match bvp::solve_z(&cfg, &prof, next, None, a.u("verbose", 0) > 0) {
                Ok(s) => {
                    let ni = s.info();
                    if ni.converged { Some(ni) } else { None }
                }
                Err(_) => None,
            };
            let ni = match ni {
                Some(n) => n,
                None => {
                    hcur /= 2.0;
                    if hcur < 1e-7 {
                        println!("{{\"cmd\":\"fail\",\"alpha_c\":{}}}", next.to_f64());
                        return;
                    }
                    continue;
                }
            };
            info = ni;
            ac = next;
            hcur = (hcur * 1.3).min(h);
            let n = family::noether(&info.steps, info.omega, cfg.rs, info.pc, info.zc);
            println!(
                "{{\"cmd\":\"row\",\"alpha_c\":{},\"omega\":\"{}\",\"N\":\"{}\",\"M\":\"{}\",\"k\":{:e},\"nodes\":{},\"newton\":{},\"fnorm\":{:e},\"h\":{:e},\"t\":{:.1}}}",
                ac.to_f64(),
                info.omega.sci(30),
                n.sci(30),
                info.mass.sci(30),
                info.k.to_f64(),
                info.nodes,
                info.fnorm.len(),
                info.fnorm.last().unwrap(),
                hcur,
                t1.elapsed().as_secs_f64()
            );
        }
    }
}

/// Arclength continuation in alpha_c from the star at omega through the omega fold, in chunks, rows at the targets
fn cmd_arc(a: &Args) {
    use multidouble::dd::DD as D;
    let eps = a.f("eps", -1.0);
    let omega0 = D::from_f64(a.f("omega", 0.7603));
    let mut cfg = bvp::StarCfg::new(eps);
    cfg.dmax = a.f("dmax", 8.0);
    let t0 = std::time::Instant::now();
    let sol = bvp::solve(&cfg, omega0, false).unwrap();
    let mut info = sol.info();
    let omega_to = a.f("omega_to", 0.45);
    let dw = a.f("dw", 0.003);
    let mut wcur = omega0.to_f64();
    while (wcur - omega_to).abs() > 1e-12 {
        let wn = if (wcur - omega_to).abs() <= dw { omega_to } else { wcur - dw * (wcur - omega_to).signum() };
        let kt = (D::ONE - D::from_f64(wn) * D::from_f64(wn)).sqrt();
        let s = bvp::solve_k_from_profile(&cfg, &info.profile(), kt, false).unwrap();
        info = s.info();
        wcur = wn;
    }
    println!("{{\"cmd\":\"arc_start\",\"alpha_c\":{},\"omega\":{},\"t\":{:.1}}}", info.alpha_c.to_f64(), info.omega.to_f64(), t0.elapsed().as_secs_f64());
    let targets = a.list("alphas");
    let alpha_stop = a.f("alpha_stop", *targets.last().unwrap_or(&1.3));
    let mut track = vec![];
    let mut logs = vec![];
    let res = cont::follow_family(&cfg, &info, &targets, alpha_stop, a.u("chunk", 6), a.f("ds", 0.05), a.f("ds_max", 2.0), a.u("verbose", 0) > 0, &mut track, &mut logs);
    for l in &logs {
        println!("{{\"cmd\":\"arc_chunk\",\"chunk\":{},\"points\":{},\"refused\":{},\"stalled\":{},\"alpha_end\":{},\"omega_end\":{:.10},\"folds\":{:?},\"t\":{:.1}}}", l.chunk, l.points, l.refused, l.stalled, l.alpha_end, l.omega_end, l.folds, l.t);
    }
    // the omega minimum along the path
    let wmin = track.iter().cloned().fold((0.0, 9.0), |m, p| if p.1 < m.1 { p } else { m });
    println!("{{\"cmd\":\"arc_track\",\"points\":{},\"omega_min_point\":[{},{}]}}", track.len(), wmin.0, wmin.1);
    match res {
        Ok(rows) => {
            for r in rows {
                println!(
                    "{{\"cmd\":\"row\",\"alpha_c\":{},\"omega\":\"{}\",\"N\":\"{}\",\"M\":\"{}\",\"k\":{:e},\"nodes\":{},\"fnorm\":{:e}}}",
                    r.alpha_c,
                    r.info.omega.sci(30),
                    r.n.sci(30),
                    r.info.mass.sci(30),
                    r.info.k.to_f64(),
                    r.info.nodes,
                    r.info.fnorm.last().unwrap()
                );
            }
        }
        Err(e) => println!("{{\"cmd\":\"arc_error\",\"err\":\"{}\"}}", e),
    }
    println!("{{\"cmd\":\"arc_wall\",\"t\":{:.1}}}", t0.elapsed().as_secs_f64());
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
    if cmd == "arc" {
        cmd_arc(&a);
        return;
    }
    if cmd == "follow" {
        cmd_follow(&a);
        return;
    }
    if cmd == "row" {
        cmd_row(&a);
        return;
    }
    lane_dispatch(&a, &cmd);
}
