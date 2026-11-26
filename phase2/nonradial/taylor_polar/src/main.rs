//! taylor_polar: the polar matching determinant D(gamma) by a Taylor integrator (negative-mass #186).
//! Usage: taylor_polar <cmd> key=value ...   (see taylor_polar.py for the calls it makes)
mod jgen;
mod mutate;
mod num;
mod planted;
mod plane;
mod polar;
mod rows;
mod ser;
mod star;
mod table;

use multidouble::dd::DD;
use num::Num;
use std::collections::HashMap;
use std::time::Instant;

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

pub fn model(a: &Args) -> table::Model {
    table::Model { omega: a.f("omega", 0.7603), epsilon: a.f("eps", -1.0), lambda: a.f("lam", 10240.0), mu: a.f("mu", 10240.0f64 * 10240.0 / 4.8) }
}

pub fn star_pars<T: Num>(m: &table::Model) -> star::StarPars<T> {
    let eps = if mutate::on("star_eps") { -m.epsilon } else { m.epsilon };
    star::StarPars { w: T::from_f64(m.omega), eps: T::from_f64(eps), lam: T::from_f64(m.lambda), mu: T::from_f64(m.mu), pi: T::pi() }
}

pub fn solve_cfg(a: &Args) -> star::SolveCfg {
    star::SolveCfg { rs: a.f("rs", 1e-8), tol: a.f("startol", 1e-30), r_norm: a.f("rnorm", 32.0), r_class: a.f("rclass", 60.0), verbose: a.u("verbose", 0) > 0 }
}

fn cmd_star(a: &Args) {
    let m = model(a);
    let p = star_pars::<DD>(&m);
    let cfg = solve_cfg(a);
    let t = Instant::now();
    let (pc, mut zc, iters) = star::solve_central(&p, a.f("zguess", 0.98), &cfg);
    if mutate::on("zc1") {
        zc = DD::from_f64(1.0);
    }
    let el = t.elapsed().as_secs_f64();
    let alpha_c = zc.recip().sqrt();
    println!(
        "{{\"cmd\":\"star\",\"pc\":\"{}\",\"zc\":\"{}\",\"phic\":{:e},\"alpha_c\":{:e},\"secant_iters\":{},\"seconds\":{:.3},\"rs\":{:e},\"startol\":{:e},\"rnorm\":{}}}",
        pc.sci(34),
        zc.sci(34),
        pc.to_f64().sqrt(),
        alpha_c.to_f64(),
        iters,
        el,
        cfg.rs,
        cfg.tol,
        cfg.r_norm
    );
}

fn parse_star(a: &Args) -> (DD, DD) {
    (DD::parse(a.s("pc", "")), DD::parse(a.s("zc", "")))
}

/// state of the star at the given radii: (alpha, a, phi0, u, m) with full digits
fn cmd_starcheck(a: &Args) {
    let m = model(a);
    let p = star_pars::<DD>(&m);
    let cfg = solve_cfg(a);
    let (pc, zc) = parse_star(a);
    let rs = a.list("r");
    let rmax = rs.iter().cloned().fold(0.0, f64::max) + 0.5;
    let st = star::build_star(pc, zc, &p, &cfg, rmax);
    for r in rs {
        let s = st.state(r);
        let bg = st.bg(r);
        println!(
            "{{\"cmd\":\"starcheck\",\"r\":{:e},\"alpha\":\"{}\",\"a\":\"{}\",\"phi0\":\"{}\",\"u\":\"{}\",\"m\":\"{}\"}}",
            r,
            bg[0].sci(32),
            bg[1].sci(32),
            bg[2].sci(32),
            bg[3].sci(32),
            s[2].sci(32)
        );
    }
}

fn lane_dispatch(a: &Args) {
    match a.s("lane", "dd") {
        "f64" => run_d::<f64>(a),
        "dd" => run_d::<DD>(a),
        "qd" => run_d::<multidouble::qd::QD>(a),
        "od" => run_d::<multidouble::od::OD>(a),
        x => panic!("lane {x}"),
    }
}

fn lane_rows(a: &Args) {
    match a.s("lane", "dd") {
        "f64" => rows::run::<f64>(a),
        "dd" => rows::run::<DD>(a),
        "qd" => rows::run::<multidouble::qd::QD>(a),
        x => panic!("lane {x}"),
    }
}

fn lane_qplant(a: &Args) {
    match a.s("lane", "dd") {
        "f64" => rows::qplant::<f64>(a),
        "dd" => rows::qplant::<DD>(a),
        "qd" => rows::qplant::<multidouble::qd::QD>(a),
        x => panic!("lane {x}"),
    }
}

fn run_d<T: Num>(a: &Args) {
    let m = model(a);
    let l = a.u("l", 2);
    let r0 = a.f("r0", 1e-3);
    let rm = a.f("rm", 24.10980364018201);
    let order = a.u("order", 30);
    let ctl = plane::Ctl { order, tol: a.f("tol", T::EPS * 100.0), hmax: a.f("hmax", 1.0), kj: a.u("kj", order) };
    let cfg = polar::Cfg { l, r0, rm, ctl, model: m };
    let sp = star_pars::<DD>(&m);
    let scfg = solve_cfg(a);
    let (pc, zc) = parse_star(a);
    let t = Instant::now();
    let st = star::build_star(pc, zc, &sp, &scfg, rm + 1.0);
    let t_star = t.elapsed().as_secs_f64();
    let t = Instant::now();
    let tab = table::build_table(&st, l, &m, r0 * 0.5, rm + 0.5);
    let t_table = t.elapsed().as_secs_f64();
    let tabl = tab.to_lane::<T>();
    let mext_dd = st.state(rm)[2];
    let mext = T::from_dd(mext_dd);
    println!(
        "{{\"cmd\":\"setup\",\"kt_multishoot\":{},\"lane\":\"{}\",\"l\":{},\"r0\":{:e},\"rm\":{},\"mext\":\"{}\",\"star_steps\":{},\"table_steps\":{},\"t_star\":{:.4},\"t_table\":{:.4},\"linear_defect\":{:e}}}",
        multishoot::Plan::optimal_kt(ctl.tol),
        T::NAME,
        l,
        r0,
        rm,
        mext.sci(130),
        st.steps.len(),
        tab.x0.len(),
        t_star,
        t_table,
        tab.linear_defect
    );
    for g in a.list("gammas") {
        let t = Instant::now();
        let o = polar::d_eval(&cfg, &tabl, mext, g);
        let wall = t.elapsed().as_secs_f64();
        println!(
            "{{\"cmd\":\"d\",\"lane\":\"{}\",\"l\":{},\"gamma\":{},\"d\":\"{}\",\"d_f64\":{:e},\"q\":{:e},\"q_int\":{:e},\"qint_res\":{:e},\"steps_y\":{},\"steps_q\":{},\"steps_s\":{},\"hmin_y\":{:e},\"hmax_y\":{:e},\"t_y\":{:.4},\"t_q\":{:.4},\"t_s\":{:.4},\"wall\":{:.4},\"order\":{},\"tol\":{:e},\"kj\":{},\"r_far\":{:e},\"r_s\":{:e},\"lnvol_y\":{:e}}}",
            T::NAME,
            l,
            g,
            o.d.sci(130),
            o.d.to_f64(),
            o.q.to_f64(),
            o.q_int.to_f64(),
            o.qint_res,
            o.steps_y.steps,
            o.steps_q,
            o.steps_s.steps,
            o.steps_y.hmin,
            o.steps_y.hmax,
            o.t_y,
            o.t_q,
            o.t_s,
            wall,
            order,
            ctl.tol,
            ctl.kj,
            o.r_far,
            o.r_s,
            o.lnvol_y.to_f64()
        );
    }
}

/// J's entries at order 0 for given background values, as the translation check against J_and_h
fn cmd_jcheck(a: &Args) {
    let m = model(a);
    let l = a.u("l", 2);
    let s = a.f("s", 1.0);
    let (r, al, aa, ph, u) = (a.f("r", 2.0), a.f("alpha", 1.0), a.f("a", 1.0), a.f("phi0", 0.01), a.f("u", -0.1));
    let ser = |x: f64| ser::Ser::<DD>::cst(DD::from_f64(x));
    let j = jgen::j_series(&ser(r), &ser(al), &ser(aa), &ser(ph), &ser(u), &table::pars(&m, l, DD::from_f64(s)));
    for i in 0..6 {
        for jj in 0..6 {
            let v = if mutate::on("jsign") && (i, jj) == (3, 2) { -j[i][jj].get(0) } else { j[i][jj].get(0) };
            println!("{{\"cmd\":\"jcheck\",\"i\":{},\"j\":{},\"v\":\"{}\"}}", i, jj, v.sci(32));
        }
    }
}

/// the jets of (r, alpha, a, phi0, u) and of J about the start of the star's step holding r, as decimal strings
fn cmd_jjets(a: &Args) {
    let m = model(a);
    let l = a.u("l", 2);
    let s = a.f("s", 1.0);
    let sp = star_pars::<DD>(&m);
    let cfg = solve_cfg(a);
    let (pc, zc) = parse_star(a);
    let r = a.f("r", 2.0);
    let st = star::build_star(pc, zc, &sp, &cfg, r + 1.0);
    let step = st.steps.iter().find(|x| x.x0 <= r && r <= x.x1).unwrap();
    let (rs, al, aa, ph, u) = table::bg_series(&step.coef, step.x0);
    let mut j = jgen::j_series(&rs, &al, &aa, &ph, &u, &table::pars(&m, l, DD::from_f64(s)));
    if mutate::on("jsign") {
        for k in 0..ser::KMAX {
            j[3][2].c[k] = -j[3][2].c[k];
        }
    }
    let nb = table::NB;
    let arr = |x: &ser::Ser<DD>| -> String { (0..nb).map(|k| format!("\"{}\"", x.get(k).sci(32))).collect::<Vec<_>>().join(",") };
    println!("{{\"cmd\":\"jjets\",\"center\":{:e},\"x1\":{:e},\"nb\":{},\"l\":{},\"s\":{},", step.x0, step.x1, nb, l, s);
    println!("\"r\":[{}],\"alpha\":[{}],\"a\":[{}],\"phi0\":[{}],\"u\":[{}],", arr(&rs), arr(&al), arr(&aa), arr(&ph), arr(&u));
    let mut ents = vec![];
    for i in 0..6 {
        for jj in 0..6 {
            if j[i][jj].n > 0 {
                ents.push(format!("\"{},{}\":[{}]", i, jj, arr(&j[i][jj])));
            }
        }
    }
    println!("\"J\":{{{}}}}}", ents.join(","));
}

/// 1/k at the lane's precision, for the recip control
fn cmd_recipcheck(a: &Args) {
    fn go<T: Num>() {
        let inv = plane::recips::<T>();
        let v: Vec<String> = (1..=ser::KMAX).map(|k| format!("\"{}\"", inv[k].sci(130))).collect();
        println!("{{\"cmd\":\"recipcheck\",\"lane\":\"{}\",\"inv\":[{}]}}", T::NAME, v.join(","));
    }
    match a.s("lane", "dd") {
        "f64" => go::<f64>(),
        "dd" => go::<DD>(),
        "qd" => go::<multidouble::qd::QD>(),
        _ => go::<multidouble::od::OD>(),
    }
}

/// debug: the table's order-0 jets against J evaluated directly from the star at the same r
fn cmd_tabcheck(a: &Args) {
    let m = model(a);
    let l = a.u("l", 2);
    let sp = star_pars::<DD>(&m);
    let cfg = solve_cfg(a);
    let (pc, zc) = parse_star(a);
    let st = star::build_star(pc, zc, &sp, &cfg, 6.0);
    let tab = table::build_table(&st, l, &m, 5e-4, 5.5);
    let s = a.f("s", 1.0);
    for r in a.list("r") {
        let mut jets = plane::Jets::<DD>::new(&tab.pattern);
        plane::JetSource::jets(&tab, r, DD::from_f64(s) * DD::from_f64(s), 30, &mut jets);
        let bg = st.bg(r);
        let c = |x: DD| ser::Ser::<DD>::cst(x);
        let d = jgen::j_series(&c(DD::from_f64(r)), &c(bg[0]), &c(bg[1]), &c(bg[2]), &c(bg[3]), &table::pars(&m, l, DD::from_f64(s)));
        let mut worst = 0.0f64;
        let mut at = (0, 0);
        let ymax = jets.ent.iter().map(|e| d[e.i][e.j].get(0).to_f64().abs()).fold(0.0f64, f64::max);
        let mut norm = 0.0f64;
        for e in &jets.ent {
            norm = norm.max((e.c[0] - d[e.i][e.j].get(0)).to_f64().abs() / ymax);
            let y = d[e.i][e.j].get(0).to_f64();
            let rel = ((e.c[0] - d[e.i][e.j].get(0)).to_f64() / y.abs().max(1e-300)).abs();
            if y != 0.0 && rel > worst {
                worst = rel;
                at = (e.i, e.j);
            }
        }
        println!("{{\"cmd\":\"tabcheck\",\"r\":{:e},\"worst\":{:e},\"norm\":{:e},\"at\":\"{},{}\"}}", r, worst, norm, at.0, at.1);
        if a.u("dump", 0) > 0 {
            for e in &jets.ent {
                println!("{{\"cmd\":\"tabdump\",\"r\":{:e},\"i\":{},\"j\":{},\"table\":\"{}\",\"direct\":\"{}\",\"alpha\":\"{}\",\"a\":\"{}\",\"phi0\":\"{}\",\"u\":\"{}\"}}", r, e.i, e.j, e.c[0].sci(33), d[e.i][e.j].get(0).sci(33), bg[0].sci(33), bg[1].sci(33), bg[2].sci(33), bg[3].sci(33));
            }
        }
    }
}

fn main() {
    let argv: Vec<String> = std::env::args().collect();
    let cmd = argv.get(1).cloned().unwrap_or_default();
    let mut map = HashMap::new();
    for kv in &argv[2.min(argv.len())..] {
        if let Some((k, v)) = kv.split_once('=') {
            map.insert(k.to_string(), v.to_string());
        }
    }
    if let Some(m) = map.get("mutate") {
        mutate::set(m);
    }
    let args = Args(map);
    // the generated J keeps hundreds of series on the stack
    let h = std::thread::Builder::new()
        .stack_size(1 << 30)
        .spawn(move || match cmd.as_str() {
            "star" => cmd_star(&args),
            "starcheck" => cmd_starcheck(&args),
            "d" => lane_dispatch(&args),
            "jcheck" => cmd_jcheck(&args),
            "jjets" => cmd_jjets(&args),
            "tabcheck" => cmd_tabcheck(&args),
            "recipcheck" => cmd_recipcheck(&args),
            "planted" => planted::run(&args),
            "rows" => lane_rows(&args),
            "rmplant" => rows::rmplant(&args),
            "qplant" => lane_qplant(&args),
            _ => eprintln!("commands: star starcheck d jcheck planted rows rmplant qplant"),
        })
        .unwrap();
    h.join().unwrap();
}
