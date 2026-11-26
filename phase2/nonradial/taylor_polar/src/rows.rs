//! The rows of convergence.py from the crate (negative-mass #187): the five variants of one (l, gamma), the r_m rule of
//! `PolarProblem`, and the planted cases for the two pieces of new code (the r_m grid rule, `interior_q`).
//!
//! One star (and one for the R = 90 variant) and one J table per l serve every gamma; each (gamma, variant) is one `d_eval`
//! and prints one JSON line, flushed, so a stopped run resumes at the next gamma.
use crate::num::Num;
use crate::plane::Ctl;
use crate::polar::{self, Cfg};
use crate::star::{self, Star};
use crate::table::{self, JTable, Model};
use crate::Args;
use multidouble::dd::DD;
use std::io::Write;
use std::time::Instant;

/// points of `np.linspace(r0, R, 20000)` in PolarProblem's r_m rule
pub const NGRID: usize = 20000;

/// numpy's `linspace(a, b, n)[i]`: i * step + a with step = (b - a) / (n - 1), and exactly b at the last point
pub fn linspace_at(a: f64, b: f64, n: usize, i: usize) -> f64 {
    if i == n - 1 {
        b
    } else {
        (i as f64) * ((b - a) / ((n - 1) as f64)) + a
    }
}

/// `rr[np.argmax(phi(rr) < cut)]` for rr = linspace(g0, rbig, NGRID): the first grid point with phi < cut, and its index
pub fn grid_rm(phi: &dyn Fn(f64) -> f64, g0: f64, rbig: f64, cut: f64) -> Option<(usize, f64)> {
    for i in 0..NGRID {
        let r = linspace_at(g0, rbig, NGRID, i);
        if phi(r) < cut {
            // mutation rm_off1: the index one too low
            let i = if crate::mutate::on("rm_off1") { i - 1 } else { i };
            return Some((i, linspace_at(g0, rbig, NGRID, i)));
        }
    }
    None
}

/// One variant of convergence.py's list, as the crate runs it.
#[derive(Clone, Copy, Debug)]
pub struct Var {
    pub name: &'static str,
    /// the interior start radius
    pub r0: f64,
    /// where the r_m grid starts, and its end R
    pub g0: f64,
    pub rbig: f64,
    pub cut: f64,
    pub order: usize,
    /// None: the lane's default (100 u)
    pub tol: Option<f64>,
}

pub fn variant(name: &str) -> Var {
    let base = Var { name: "base", r0: 1e-3, g0: 1e-3, rbig: 60.0, cut: 1e-9, order: 30, tol: None };
    match name {
        "base" => base,
        // the analogue of the rtol variant: a looser step control of the same integrator (the crate has no scipy rtol)
        "tol" => Var { name: "tol", tol: Some(if crate::mutate::on("tol_loose") { 1e-1 } else { 1e-12 }), ..base },
        "phi" => Var { name: "phi", cut: if crate::mutate::on("rm_cut") { 1e-11 } else { 1e-12 }, ..base },
        "r90" => Var { name: "r90", rbig: 90.0, ..base },
        "r0" => Var { name: "r0", r0: 3e-3, g0: if crate::mutate::on("r0_grid") { 1e-3 } else { 3e-3 }, ..base },
        // order 24 at tol 1e-22, #186's V3
        "o24" => Var { name: "o24", order: 24, tol: Some(1e-22), ..base },
        x => panic!("variant {x}"),
    }
}

struct StarSet<T: Num> {
    st: Star<DD>,
    tab: JTable<T>,
    t_star: f64,
    t_table: f64,
}

pub fn run<T: Num>(a: &Args) {
    let m: Model = crate::model(a);
    let l = a.u("l", 2);
    let names: Vec<String> = a.s("variants", "base,tol,phi,r90,r0").split(',').map(|s| s.to_string()).collect();
    let vars: Vec<Var> = names.iter().map(|n| variant(n)).collect();
    let sp = crate::star_pars::<DD>(&m);
    let scfg = crate::solve_cfg(a);
    // one star serves every variant: its central data do not depend on the classification radius (60 or 90: the same
    // 33 digits), so the R = 90 variant differs from base by its r_m grid only
    let (pc, zc) = (DD::parse(a.s("pc", "")), DD::parse(a.s("zc", "")));
    let t = Instant::now();
    let st = star::build_star(pc, zc, &sp, &scfg, 36.0);
    let t_star = t.elapsed().as_secs_f64();
    let mut rms: Vec<(usize, f64)> = vec![];
    let mut setups = vec![];
    for v in &vars {
        let phi = |r: f64| st.bg(r)[2].to_f64();
        let (i, rm) = grid_rm(&phi, v.g0, v.rbig, v.cut).expect("no r_m");
        // the margin of the grid rule: how far phi at r_m and at the point before it lie from the cut, in ln
        let prev = linspace_at(v.g0, v.rbig, NGRID, i.saturating_sub(1));
        rms.push((i, rm));
        setups.push(format!(
            "{{\"variant\":\"{}\",\"rm\":{:e},\"rm_i\":{},\"cut\":{:e},\"ln_phi_over_cut_at_rm\":{:e},\"ln_phi_over_cut_before\":{:e}}}",
            v.name,
            rm,
            i,
            v.cut,
            (phi(rm) / v.cut).ln(),
            (phi(prev) / v.cut).ln()
        ));
    }
    // J's table from the smallest start radius to the largest r_m + 0.5
    let rmax = rms.iter().map(|r| r.1).fold(0.0, f64::max) + 0.5;
    let rmin = vars.iter().map(|v| v.r0).fold(f64::INFINITY, f64::min) * 0.5;
    let t = Instant::now();
    let tab = table::build_table(&st, l, &m, rmin, rmax);
    let set = StarSet { st, tab: tab.to_lane::<T>(), t_star, t_table: t.elapsed().as_secs_f64() };
    println!(
        "{{\"cmd\":\"rowsetup\",\"lane\":\"{}\",\"l\":{},\"variants\":[{}],\"star_steps\":{},\"t_star\":{:.3},\"t_table\":{:.3}}}",
        T::NAME,
        l,
        setups.join(","),
        set.st.steps.len(),
        set.t_star,
        set.t_table
    );
    std::io::stdout().flush().unwrap();
    for g in a.list("gammas") {
        for (v, &(ri, rm)) in vars.iter().zip(&rms) {
            let ctl = Ctl { order: v.order, tol: v.tol.unwrap_or(T::EPS * 100.0), hmax: 1.0, kj: v.order };
            let cfg = Cfg { l, r0: v.r0, rm, ctl, model: m };
            let mext = T::from_dd(set.st.state(rm)[2]);
            let t = Instant::now();
            let o = polar::d_eval(&cfg, &set.tab, mext, g);
            let wall = t.elapsed().as_secs_f64();
            println!(
                "{{\"cmd\":\"row\",\"variant\":\"{}\",\"lane\":\"{}\",\"l\":{},\"gamma\":{},\"d\":\"{}\",\"d_f64\":{:e},\"q_ext\":\"{}\",\"q_ext_f64\":{:e},\"q_int\":\"{}\",\"q_int_f64\":{:e},\"qint_res\":{:e},\"rm\":{:e},\"rm_i\":{},\"mext\":\"{}\",\"steps_y\":{},\"steps_q\":{},\"steps_s\":{},\"wall\":{:.4},\"order\":{},\"tol\":{:e}}}",
                v.name,
                T::NAME,
                l,
                g,
                o.d.sci(130),
                o.d.to_f64(),
                o.q.sci(130),
                o.q.to_f64(),
                o.q_int.sci(130),
                o.q_int.to_f64(),
                o.qint_res,
                rm,
                ri,
                mext.sci(130),
                o.steps_y.steps,
                o.steps_q,
                o.steps_s.steps,
                wall,
                v.order,
                ctl.tol
            );
            std::io::stdout().flush().unwrap();
        }
    }
}

/// planted r_m: phi(r) = a exp(-kappa r) through the grid rule, for the control of the rule against numpy's linspace/argmax
pub fn rmplant(a: &Args) {
    let (amp, kap) = (a.f("a", 0.02), a.f("kappa", 0.65));
    let phi = move |r: f64| amp * (-kap * r).exp();
    let (g0, rbig, cut) = (a.f("g0", 1e-3), a.f("R", 60.0), a.f("cut", 1e-9));
    let (i, rm) = grid_rm(&phi, g0, rbig, cut).expect("no r_m");
    println!("{{\"cmd\":\"rmplant\",\"i\":{},\"rm\":{:e}}}", i, rm);
}

fn dy(i: i64) -> f64 {
    i as f64 / 8.0
}

/// planted planes for `interior_q`: dyadic Y (3 x 6), S (2 x 4) with Y[2:] c0 in span S for c0 = (1/2, -1/4, 3/8);
/// case "near": column 3 of Y is column 2 plus 2^-27 w (nearly degenerate planes). The python check builds the same
/// numbers with exact rationals.
pub fn qplant_inputs<T: Num>(near: bool) -> ([[T; 6]; 3], [[T; 4]; 2]) {
    let mut y = [[0.0f64; 6]; 3];
    for c in 0..3 {
        for i in 0..6 {
            y[c][i] = dy(((3 * i as i64 + 5 * c as i64 + 1) % 7) - 3);
        }
    }
    if near {
        for i in 0..6 {
            let w = (((i as i64 + 2) % 3) - 1) as f64;
            y[2][i] = y[1][i] + w / (1u64 << 27) as f64;
        }
    }
    let c0 = [0.5, -0.25, 0.375];
    let mut s = [[0.0f64; 4]; 2];
    for i in 0..4 {
        s[0][i] = (0..3).map(|c| y[c][2 + i] * c0[c]).sum();
        s[1][i] = (((2 * i as i64 + 1) % 5) - 2) as f64 / 4.0;
    }
    let cv = |x: f64| T::from_f64(x);
    let mut yo = [[T::zero(); 6]; 3];
    let mut so = [[T::zero(); 4]; 2];
    for c in 0..3 {
        for i in 0..6 {
            yo[c][i] = cv(y[c][i]);
        }
    }
    for k in 0..2 {
        for i in 0..4 {
            so[k][i] = cv(s[k][i]);
        }
    }
    (yo, so)
}

pub fn qplant<T: Num>(a: &Args) {
    let near = a.s("case", "well") == "near";
    let (y, s) = qplant_inputs::<T>(near);
    let (mext, rm, l) = (T::from_f64(1.0 / 64.0), T::from_f64(24.125), a.u("l", 3));
    let (q, res) = polar::interior_q(&y, &s, mext, rm, l);
    println!("{{\"cmd\":\"qplant\",\"lane\":\"{}\",\"case\":\"{}\",\"q_int\":\"{}\",\"q_int_f64\":{:e},\"res\":{:e}}}", T::NAME, if near { "near" } else { "well" }, q.sci(130), q.to_f64(), res);
}
