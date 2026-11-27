//! taylor_radial: the radial boson-star modes and the phantom family by Taylor integration (negative-mass #194).
//! Usage: taylor_radial <cmd> key=value ...   (see ../../taylor_radial.py for the calls it makes)
mod bvp;
mod checks;
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
    if mutate::on("shortR") {
        return 3.0;
    }
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
    let omega = <D as Num>::parse(a.s("omega", "0.94"));
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
        for ((s, _), v) in found.vals.iter().zip(&found.val_str) {
            println!("{{\"cmd\":\"D\",\"s2\":{:e},\"D\":\"{}\"}}", s, v);
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

/// the radius where P falls to 1e-6 P_c (radial_modes' r_m rule on an exact star) and max |2m/r|
fn star_size(info: &bvp::StarInfo) -> (f64, f64) {
    let prof = info.profile();
    let pc = prof.state(info.cfg.rs)[0].to_f64();
    let mut rs = f64::NAN;
    let mut mx = 0.0f64;
    let mut r = 0.25;
    while r < prof.r_t {
        let s = prof.state(r);
        mx = mx.max((2.0 * s[2].to_f64() / r).abs());
        if rs.is_nan() && s[0].to_f64() < 1e-6 * pc {
            rs = r;
        }
        r += 0.25;
    }
    (mx, rs)
}

/// Arclength continuation in alpha_c from the star at omega through the omega fold, in chunks; at each stated alpha_c a row,
/// and (alphas in `modes`) the roots of D on the declared grid and a 4 x denser one, and (alphas in `zvariants`) the lowest
/// root under the variants of the method
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
    let modes = a.list("modes");
    let zvar = a.list("zvariants");
    let alpha_stop = a.f("alpha_stop", *targets.last().unwrap_or(&1.3));
    let mut track = vec![];
    let mut logs = vec![];
    let dir_mut = mutate::on("arcfold");
    let cfg2 = cfg.clone();
    let mut on_row = |r: &cont::FamRow| {
        let (mx, rs) = star_size(&r.info);
        println!(
            "{{\"cmd\":\"row\",\"alpha_c\":{},\"omega\":\"{}\",\"N\":\"{}\",\"M\":\"{}\",\"k\":{:e},\"max2mr\":{:.6},\"rs\":{:.4},\"nodes\":{},\"fnorm\":{:e}}}",
            r.alpha_c,
            r.info.omega.sci(30),
            r.n.sci(30),
            r.info.mass.sci(30),
            r.info.k.to_f64(),
            mx,
            rs,
            r.info.nodes,
            r.info.fnorm.last().unwrap()
        );
        let is_mode = modes.iter().any(|&x| (x - r.alpha_c).abs() < 1e-9);
        let is_z = zvar.iter().any(|&x| (x - r.alpha_c).abs() < 1e-9);
        if !(is_mode || is_z) {
            return;
        }
        let thr = (1.0 - r.info.omega.to_f64()).powi(2);
        let t1 = std::time::Instant::now();
        let prep = match family::prepare(&cfg2, clone_info(&r.info), 0.0, 0.995 * thr) {
            Ok(p) => p,
            Err(e) => {
                println!("{{\"cmd\":\"mroots_error\",\"alpha_c\":{},\"err\":\"{}\"}}", r.alpha_c, e);
                return;
            }
        };
        let t_prep = t1.elapsed().as_secs_f64();
        let (roots, ch, ng, wall) = family::scan_dd(&prep, 30, 1.232595164407831e-30, 1.0, 24, 32, 0.2);
        println!("{{\"cmd\":\"mroots\",\"alpha_c\":{},\"grid\":\"base\",\"ngrid\":{},\"sign_changes\":{},\"roots\":{:?},\"thr\":{:e},\"r_far\":{:.1},\"t_prep\":{:.1},\"wall\":{:.1}}}", r.alpha_c, ng, ch, roots, thr, prep.r_far, t_prep, wall);
        if is_mode {
            let (roots, ch, ng, wall) = family::scan_dd(&prep, 30, 1.232595164407831e-30, 1.0, 96, 128, 0.2);
            println!("{{\"cmd\":\"mroots\",\"alpha_c\":{},\"grid\":\"dense\",\"ngrid\":{},\"sign_changes\":{},\"roots\":{:?},\"thr\":{:e},\"wall\":{:.1}}}", r.alpha_c, ng, ch, roots, thr, wall);
        }
        if is_z {
            let mut variant = |name: &str, order: usize, tol: f64, scale: f64, p: &family::Prepared| {
                let (roots, ch, _ng, wall) = family::scan_dd(p, order, tol, scale, 24, 32, 0.2);
                println!("{{\"cmd\":\"zvariant\",\"alpha_c\":{},\"variant\":\"{}\",\"sign_changes\":{},\"roots\":{:?},\"wall\":{:.1}}}", r.alpha_c, name, ch, roots, wall);
            };
            variant("R x1.5", 30, 1.232595164407831e-30, 1.5, &prep);
            variant("order 24, tol 1e-22", 24, 1e-22, 1.0, &prep);
            let mut cfg3 = cfg2.clone();
            cfg3.dmax = 6.0;
            if let Ok(s3) = bvp::solve_z(&cfg3, &r.info.profile(), r.info.alpha_c, None, false) {
                let i3 = s3.info();
                if let Ok(p3) = family::prepare(&cfg3, i3, 0.0, 0.995 * thr) {
                    variant("layout dmax 6", 30, 1.232595164407831e-30, 1.0, &p3);
                }
            }
            // f64 lane
            let rp = prep.radial::<f64>(30, f64::EPSILON * 50.0);
            let g = family::grid(0.2, thr, 24, 32);
            let f = family::scan(&rp, &g, 19.0);
            println!("{{\"cmd\":\"zvariant\",\"alpha_c\":{},\"variant\":\"f64 lane\",\"roots\":{:?}}}", r.alpha_c, f.roots.iter().map(|x| x.0.sci(17)).collect::<Vec<_>>());
        }
    };
    let res = cont::follow_family(&cfg, &info, &targets, alpha_stop, a.u("chunk", 10), a.f("ds", 0.05), a.f("ds_max", 5.0), a.u("verbose", 0) > 0, &mut track, &mut logs, dir_mut, &mut on_row);
    for l in &logs {
        println!("{{\"cmd\":\"arc_chunk\",\"chunk\":{},\"points\":{},\"refused\":{},\"stalled\":{},\"alpha_end\":{},\"omega_end\":{:.10},\"folds\":{:?},\"t\":{:.1}}}", l.chunk, l.points, l.refused, l.stalled, l.alpha_end, l.omega_end, l.folds, l.t);
    }
    let wmin = track.iter().cloned().fold((0.0, 9.0), |m, p| if p.1 < m.1 { p } else { m });
    println!("{{\"cmd\":\"arc_track\",\"points\":{},\"omega_min_point\":[{},{}]}}", track.len(), wmin.0, wmin.1);
    if let Err(e) = res {
        println!("{{\"cmd\":\"arc_error\",\"err\":\"{}\"}}", e);
    }
    println!("{{\"cmd\":\"arc_wall\",\"t\":{:.1}}}", t0.elapsed().as_secs_f64());
}

/// the fold of omega: arclength in k from the star at omega_to through the minimum of omega, located by multisolver
fn cmd_fold(a: &Args) {
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
    // plain steps in alpha_c to just before the fold (pinning k is ill-posed beside its maximum)
    let a_to = a.f("alpha_to", 1.15);
    let mut ac = info.alpha_c;
    let mut hcur = 0.01f64;
    while ac.to_f64() < a_to - 1e-12 {
        let next = if a_to - ac.to_f64() <= hcur { D::from_f64(a_to) } else { ac + D::from_f64(hcur) };
        match bvp::solve_z(&cfg, &info.profile(), next, None, false) {
            Ok(sz) if sz.info().converged => {
                info = sz.info();
                ac = next;
                hcur = (hcur * 1.3).min(0.01);
            }
            _ => {
                hcur /= 2.0;
                assert!(hcur > 1e-7, "alpha_c steps failed near {}", ac.to_f64());
            }
        }
    }
    match cont::fold_search(&cfg, &info, a.u("chunk", 40), a.f("ds", 0.01), a.f("ds_max", 0.04), a.u("max_chunks", 8), a.u("verbose", 0) > 0) {
        Ok(f) => {
            let det = f.detected.map(|(x, y)| format!("[{},{}]", x.sci(20), y.sci(20))).unwrap_or("null".into());
            let (om, k) = match f.k_fold {
                Some(fk) => (format!("\"{}\"", (D::ONE - fk * fk).sqrt().sci(25)), format!("\"{}\"", fk.sci(25))),
                None => ("null".to_string(), "null".to_string()),
            };
            println!(
                "{{\"cmd\":\"fold\",\"detected_between_k\":{},\"k\":{},\"omega_min\":{},\"alpha_c\":{},\"points\":{},\"locate_corrections\":{},\"locate_error\":{},\"t\":{:.1}}}",
                det,
                k,
                om,
                f.alpha_c.map(|x| x.to_string()).unwrap_or("null".into()),
                f.points,
                f.iters,
                f.locate_err.map(|e| format!("\"{}\"", e)).unwrap_or("null".into()),
                t0.elapsed().as_secs_f64()
            );
        }
        Err(e) => println!("{{\"cmd\":\"fold_error\",\"err\":\"{}\"}}", e),
    }
}

fn clone_info(i: &bvp::StarInfo) -> bvp::StarInfo {
    bvp::StarInfo { eps: i.eps, cfg: i.cfg.clone(), steps: i.steps.clone(), k: i.k, omega: i.omega, pc: i.pc, zc: i.zc, alpha_c: i.alpha_c, mass: i.mass, r_t: i.r_t, nodes: i.nodes, fnorm: i.fnorm.clone(), converged: i.converged }
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
    match cmd.as_str() {
        "jeval" => return checks::cmd_jeval(&a),
        "gcheck" => return checks::cmd_gcheck(&a),
        "starcheck" => return checks::cmd_starcheck(&a),
        "s2" => return checks::cmd_s2(&a),
        "pd" => return checks::cmd_pd(&a),
        _ => {}
    }
    if cmd == "fold" {
        cmd_fold(&a);
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
