//! Commands that serve the checks of taylor_radial.py: J at a point (G1), the table's shifted jets against a direct series
//! (G2), the star at radii (S1), the two Taylor recursions of the star at one state (S2), and the planted D at given s (P2).
use crate::family;
use crate::num::Num;
use crate::ser::*;
use crate::table::*;
use crate::Args;
use multidouble::dd::DD;

fn lane_run(a: &Args, f32: impl Fn(&Args), fdd: impl Fn(&Args), fqd: impl Fn(&Args)) {
    match a.s("lane", "dd") {
        "f64" => f32(a),
        "dd" => fdd(a),
        "qd" => fqd(a),
        x => panic!("lane {x}"),
    }
}

/// J(r; s2) at a point: the entries as strings in the lane. Arguments are decimal strings parsed in the lane.
fn jeval<T: Num>(a: &Args) {
    let g = |k: &str| T::parse(a.s(k, "0"));
    let r = Ser::var(g("r"));
    let (alpha, am, phi0, u) = (Ser::cst(g("alpha")), Ser::cst(g("a")), Ser::cst(g("phi0")), Ser::cst(g("u")));
    let p = crate::rgen::Pars { omega: g("omega"), epsilon: g("eps"), s2: g("s2"), lambda: crate::model::lam::<T>(), mu: crate::model::mu::<T>(), pi: crate::model::pi::<T>() };
    let j = crate::rgen::j_series(&r, &alpha, &am, &phi0, &u, &p);
    for i in 0..4 {
        for k in 0..4 {
            println!("{{\"cmd\":\"J\",\"i\":{},\"j\":{},\"v\":\"{}\"}}", i, k, j[i][k].get(0).sci(60));
        }
    }
    println!("{{\"cmd\":\"J_const\",\"mu\":\"{}\",\"pi\":\"{}\"}}", crate::model::mu::<T>().sci(60), crate::model::pi::<T>().sci(60));
}

pub fn cmd_jeval(a: &Args) {
    lane_run(a, |a| jeval::<f64>(a), |a| jeval::<DD>(a), |a| jeval::<multidouble::qd::QD>(a));
}

/// the table's shifted jets at several x0 against the direct series of J from the star's own series about x0
pub fn cmd_gcheck(a: &Args) {
    let eps = a.f("eps", 1.0);
    let omega = DD::from_f64(a.f("omega", 0.94));
    let mut cfg = crate::bvp::StarCfg::new(eps);
    cfg.dmax = a.f("dmax", 8.0);
    let prep = family::prepare_k(&cfg, omega, 0.0).unwrap();
    let tab = &prep.tab;
    let model = Model { omega: prep.omega, epsilon: eps };
    let s2 = DD::from_f64(a.f("s2", 0.003));
    let xs = a.list("x0");
    let mut worst = 0.0f64;
    for &x0 in &xs {
        // the star step holding x0 and its series about x0
        let st = prep.info.steps.iter().find(|s| s.lo <= x0 && x0 <= s.hi).expect("x0 outside the star's steps");
        let mut co = st.coef;
        let delta = DD::from_f64(x0) - DD::from_f64(st.x0);
        for f in 0..4 {
            taylor_shift(&mut co[f][..=NS], delta, NS + 1);
        }
        let (r, alpha, am, phi0, u) = bg_series(&co, x0);
        let direct = crate::rgen::j_series(&r, &alpha, &am, &phi0, &u, &pars(&model, s2));
        let mut jets = crate::plane::Jets::<DD>::new(&crate::plane::JetSource::<DD>::pattern(tab));
        crate::plane::JetSource::<DD>::jets(tab, x0, s2, NB, &mut jets);
        let mut wx = 0.0f64;
        for en in &jets.ent {
            for k in 0..NB {
                let d = direct[en.i][en.j].get(k);
                let t = en.c[k];
                let scale = d.to_f64().abs().max(1e-300);
                if scale > 0.0 {
                    wx = wx.max(((t - d).to_f64() / scale).abs().min(1.0));
                }
            }
        }
        // the first 12 coefficients only: the table's high orders carry the shift's amplification of series noise
        let mut w12 = 0.0f64;
        for en in &jets.ent {
            for k in 0..12 {
                let d = direct[en.i][en.j].get(k);
                let t = en.c[k];
                let scale = d.to_f64().abs().max(1e-300);
                w12 = w12.max(((t - d).to_f64() / scale).abs().min(1.0));
            }
        }
        worst = worst.max(w12);
        println!("{{\"cmd\":\"G2\",\"x0\":{},\"worst_all\":{:e},\"worst_first12\":{:e}}}", x0, wx, w12);
    }
    println!("{{\"cmd\":\"G2_done\",\"worst_first12\":{:e}}}", worst);
}

/// (alpha, a, phi0, u, m) of the star at radii, in the lane's digits
pub fn cmd_starcheck(a: &Args) {
    let eps = a.f("eps", 1.0);
    let omega = DD::parse(a.s("omega", "0.75"));
    let mut cfg = crate::bvp::StarCfg::new(eps);
    cfg.dmax = a.f("dmax", 8.0);
    let sol = crate::bvp::solve(&cfg, omega, false).unwrap();
    let info = sol.info();
    let prof = info.profile();
    println!("{{\"cmd\":\"star_central\",\"pc\":\"{}\",\"zc\":\"{}\",\"k\":\"{}\",\"omega\":\"{}\",\"mass\":\"{}\"}}", info.pc.sci(32), info.zc.sci(32), info.k.sci(32), info.omega.sci(32), info.mass.sci(32));
    for r in a.list("r") {
        let s = prof.state(r);
        let rr = DD::from_f64(r);
        let alpha = s[3].recip().sqrt();
        let am = (DD::ONE - DD::from_f64(2.0) * s[2] / rr).recip().sqrt();
        println!("{{\"cmd\":\"starcheck\",\"r\":{},\"alpha\":\"{}\",\"a\":\"{}\",\"phi0\":\"{}\",\"u\":\"{}\",\"m\":\"{}\"}}", r, alpha.sci(32), am.sci(32), s[0].sqrt().sci(32), s[1].sci(32), s[2].sci(32));
    }
}

/// S2: the Taylor coefficients of (P, u, m, Z) about a node state by multishoot's recursion (Equations::taylor_coeffs on the
/// double-double lane) and by the crate's own Num-lane recursion of star_ref.rs
pub fn cmd_s2(a: &Args) {
    use multishoot::lane::DV;
    use multishoot::{Equations, Fl, Work};
    let eps = a.f("eps", 1.0);
    let omega = DD::parse(a.s("omega", "0.75"));
    let mut cfg = crate::bvp::StarCfg::new(eps);
    cfg.dmax = a.f("dmax", 8.0);
    let sol = crate::bvp::solve(&cfg, omega, false).unwrap();
    let k = sol.k();
    let x0 = a.f("x0", 4.5);
    let node = sol.sh.nodes.iter().position(|&x| (x - x0).abs() < 1e-9).expect("x0 is not a node");
    let st: [DD; 4] = [0, 1, 2, 3].map(|f| crate::bvp::zget(&sol.z, crate::bvp::zi(&sol.sh, node, f, 0)));
    let kt = 30;
    let mut w = Work::<DV, 1>::new(4, 0, kt);
    for f in 0..4 {
        w.c[f][0] = Fl::splat(st[f]);
    }
    let g = [DV::splat(k)];
    sol.sh.eq.taylor_coeffs::<DV>(x0, &g, &mut w);
    let p = crate::star_ref::StarPars::<DD> { w: sol.omega(), eps: DD::from_f64(eps), lam: DD::from_f64(crate::model::LAM), mu: crate::model::mu_dd(), pi: crate::model::pi_dd() };
    let inv = crate::plane::recips::<DD>();
    let ref_c = crate::star_ref::star_coeffs::<DD>(DD::from_f64(x0), st, kt, &p, &inv);
    let mut worst = 0.0f64;
    for f in 0..4 {
        for kk in 0..=kt {
            let a1 = w.c[f][kk].get(0);
            let a2 = ref_c[f][kk];
            let sc = a2.to_f64().abs().max(1e-300);
            let rel = ((a1 - a2).to_f64() / sc).abs();
            if kk <= 28 && a2.to_f64().abs() > 1e-250 {
                worst = worst.max(rel);
            }
        }
    }
    println!("{{\"cmd\":\"S2\",\"x0\":{},\"worst_rel\":{:e}}}", x0, worst);
}

/// D of the planted problem at the given s in the lane
fn pd<T: Num>(a: &Args) {
    use crate::plane::Ctl;
    let order = a.u("order", crate::default_order::<T>());
    let ctl = Ctl { order, tol: a.f("tol", T::EPS * 100.0), hmax: a.f("hmax", 2.0), kj: order };
    let sign = if crate::mutate::on("plant_sign") { -T::one() } else { T::one() };
    let r_m = a.f("rm", 3.0);
    let r_far = a.f("rfar", 3.0 + crate::far_factor::<T>() / (4.0f64 - 2.0 * 1.9).sqrt().min(1.0));
    for s in a.s("slist", "0").split(',') {
        let st = T::parse(s);
        let d = crate::planted::d_planted::<T>(sign, st, r_m, r_far, &ctl);
        println!("{{\"cmd\":\"pd\",\"s\":\"{}\",\"d\":\"{}\"}}", s, d.sci(60));
    }
}

pub fn cmd_pd(a: &Args) {
    lane_run(a, |a| pd::<f64>(a), |a| pd::<DD>(a), |a| pd::<multidouble::qd::QD>(a));
}

/// S3: the star's central data (P_c, Z_c) and its mass M in the quad-double lane (the crate's own Taylor shooting, `star_ref.rs`,
/// normalizing alpha = 1 at `rnorm`) against the double-double boundary-value solution (multishoot)
pub fn cmd_s3(a: &Args) {
    use crate::star_ref::{shoot, solve_central, SolveCfg, StarPars};
    use multidouble::qd::QD;
    let eps = a.f("eps", 1.0);
    let omega = DD::parse(a.s("omega", "0.75"));
    let mut cfg = crate::bvp::StarCfg::new(eps);
    cfg.dmax = a.f("dmax", 8.0);
    let t0 = std::time::Instant::now();
    let info = crate::bvp::solve(&cfg, omega, false).unwrap().info();
    let t_dd = t0.elapsed().as_secs_f64();
    let p = StarPars::<QD> { w: QD::from_dd(omega), eps: QD::from_f64(eps), lam: QD::from_f64(crate::model::LAM), mu: crate::model::mu::<QD>(), pi: crate::model::pi::<QD>() };
    let r_norm = a.f("rnorm", 60.0);
    let scfg = SolveCfg { rs: cfg.rs, tol: a.f("tol", 1e-62), r_norm, r_class: a.f("rclass", r_norm + 30.0), verbose: a.u("verbose", 0) > 0 };
    let t1 = std::time::Instant::now();
    let (pc, zc, it) = solve_central::<QD>(&p, 1.0 - 0.02 * eps, &scfg);
    let run = shoot::<QD>(pc, zc, &p, cfg.rs, r_norm, false, scfg.tol, false);
    let m_qd = run.state_end[2];
    let t_qd = t1.elapsed().as_secs_f64();
    let rel = |q: QD, d: DD| ((q - QD::from_dd(d)) / q).abs().to_f64();
    let (e_pc, e_zc, e_m) = (rel(pc, info.pc), rel(zc, info.zc), rel(m_qd, info.mass));
    println!(
        "{{\"cmd\":\"S3\",\"eps\":{},\"omega\":\"{}\",\"rnorm\":{},\"tol\":{:e},\"secant_iters\":{},\"pc_qd\":\"{}\",\"pc_dd\":\"{}\",\"zc_qd\":\"{}\",\"zc_dd\":\"{}\",\"M_qd\":\"{}\",\"M_dd\":\"{}\",\"rel_pc\":{:e},\"rel_zc\":{:e},\"rel_M\":{:e},\"worst\":{:e},\"t_dd\":{:.1},\"t_qd\":{:.1}}}",
        eps,
        omega.sci(20),
        r_norm,
        scfg.tol,
        it,
        pc.sci(40),
        info.pc.sci(40),
        zc.sci(40),
        info.zc.sci(40),
        m_qd.sci(40),
        info.mass.sci(40),
        e_pc,
        e_zc,
        e_m,
        e_pc.max(e_zc).max(e_m),
        t_dd,
        t_qd
    );
}
