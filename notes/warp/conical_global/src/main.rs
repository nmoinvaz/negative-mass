//! conical_global: command line for conical_global.py. Commands print JSON lines on stdout.
//!   walk FI C [X KT TOL HMAX]   state at X by the inward walk from 200 (default: the verdict)
//!   series FI C                 end states at +-200 and N there
//!   ntable FILE                 N at points "f0 f1 f2 eta" (rationals), value and scale as [hi, lo]
//!   scan FI_LIST C_LIST         one verdict per pair
//!   bisect FI CLO CHI           bisection in C between two verdict classes
//!   plant NAME                  blowup | linear | cubic
//! Any argument --mutate=NAME turns on a deliberate error of the crate.
use conical_global::connect::{connect, ConnectIn};
use conical_global::json;
use conical_global::model::{conical, conical_poly, end_state, Implicit, Mutation, NoEnds};
use conical_global::planted;
use conical_global::poly::{abs_dd, PolyN};
use conical_global::scan::{class_label, conical_kind, lane0, start_state, verdict};
use conical_global::walk::{walk, Outcome, WalkOpts};
use multishoot::lane::DV;
use multishoot::{Fl, Plan, DD};

fn parse_dd(s: &str) -> DD {
    if let Some((a, b)) = s.split_once('/') {
        DD::from_f64(a.parse().unwrap()) / DD::from_f64(b.parse().unwrap())
    } else {
        DD::from_f64(s.parse().unwrap())
    }
}

fn main() {
    let mut args: Vec<String> = std::env::args().skip(1).collect();
    let mut m = Mutation::default();
    args.retain(|a| {
        if let Some(n) = a.strip_prefix("--mutate=") {
            m = Mutation::from_name(n).unwrap_or_else(|| panic!("unknown mutation {n}"));
            false
        } else {
            true
        }
    });
    let cmd = args.first().cloned().unwrap_or_default();
    let mut opts = WalkOpts::standard();
    opts.flip_side = m.plant_locator;
    match cmd.as_str() {
        "walk" => {
            let fi: f64 = args[1].parse().unwrap();
            let c = parse_dd(&args[2]);
            if args.len() == 3 {
                println!("{}", json::verdict(&verdict(fi, c, &m, &opts)));
                return;
            }
            let x: f64 = args[3].parse().unwrap();
            let kt: usize = args[4].parse().unwrap();
            let tol: f64 = args[5].parse().unwrap();
            let hmax: f64 = args[6].parse().unwrap();
            opts.plan = Plan { kt, tol, hmax };
            let eq = conical(m);
            let w = walk(&eq, &eq.n, start_state(fi, c, &m), 200.0, x, DD::from_f64(fi), &opts);
            let ok = matches!(w.outcome, Outcome::End);
            println!(
                "{{\"end\": {}, \"steps\": {}, \"x\": {}, \"state\": {}}}",
                ok,
                w.steps.len(),
                json::num(w.last().x),
                json::s3(&w.last().s)
            );
        }
        "series" => {
            let fi = parse_dd(&args[1]);
            let c = parse_dd(&args[2]);
            let p = conical_poly(&m);
            let mut out = Vec::new();
            for sign in [1.0, -1.0] {
                let s = end_state::<DV>(Fl::splat(fi), Fl::splat(c), sign, &m);
                let st = [lane0(&s[0]), lane0(&s[1]), lane0(&s[2])];
                let (v, sc) = p.eval(st[0], st[1], st[2], DD::from_f64(200.0 * sign));
                out.push(format!(
                    "{{\"sign\": {}, \"state\": {}, \"n\": {}, \"scale\": {}, \"n_rel\": {}}}",
                    sign,
                    json::s3(&st),
                    json::dd(v),
                    json::dd(sc),
                    json::num((abs_dd(v) / sc).hi)
                ));
            }
            println!("[{}]", out.join(", "));
        }
        "ntable" => {
            let p = conical_poly(&m);
            let txt = std::fs::read_to_string(&args[1]).unwrap();
            for line in txt.lines().filter(|l| !l.trim().is_empty()) {
                let v: Vec<DD> = line.split_whitespace().map(parse_dd).collect();
                let (val, sc) = p.eval(v[0], v[1], v[2], v[3]);
                println!("{{\"n\": {}, \"scale\": {}}}", json::dd(val), json::dd(sc));
            }
        }
        "scan" => {
            let fis: Vec<f64> = args[1].split(',').map(|s| s.parse().unwrap()).collect();
            let cs: Vec<DD> = args[2].split(',').map(parse_dd).collect();
            for &fi in &fis {
                for &c in &cs {
                    println!("{}", json::verdict(&verdict(fi, c, &m, &opts)));
                }
            }
        }
        "forward" => {
            // forward F0 F1 ...: pairs f,f1 as lists
            let fs: Vec<f64> = args[1].split(',').map(|s| s.parse().unwrap()).collect();
            let f1s: Vec<f64> = args[2].split(',').map(|s| s.parse().unwrap()).collect();
            for &f in &fs {
                for &f1 in &f1s {
                    let (a, b, n0) = conical_global::scan::forward(f, f1, &m, &opts);
                    let sd = |s: &conical_global::scan::Side| {
                        let fit = s.fit.map(|(x, c)| format!("[{}, {}]", json::num(x), json::num(c))).unwrap_or("null".into());
                        format!("{{\"class\": \"{}\", \"kind\": \"{}\", \"xstar\": {}, \"bar\": {}, \"fit\": {}, \"mismatch\": {}, \"steps\": {}}}", s.class, s.kind, json::opt(s.xstar), json::opt(s.bar), fit, json::num(s.mismatch), s.steps)
                    };
                    println!("{{\"f\": {}, \"f1\": {}, \"n0\": {}, \"plus\": {}, \"minus\": {}}}", json::num(f), json::num(f1), json::num(n0), sd(&a), sd(&b));
                }
            }
        }
        "fwdstate" => {
            // fwdstate F0 F1 X: the state at X of the walk from eta = 0 (diagnostic)
            let (f, f1, x): (f64, f64, f64) = (args[1].parse().unwrap(), args[2].parse().unwrap(), args[3].parse().unwrap());
            let eq = conical(m);
            let (fd, f1d) = (DD::from_f64(f), DD::from_f64(f1));
            let s0 = [fd, f1d, -(f1d * f1d) / fd];
            let w = walk(&eq, &eq.n, s0, 0.0, x, fd, &opts);
            println!("{{\"outcome\": \"{:?}\", \"x\": {}, \"state\": {}, \"n_rel\": {}, \"steps\": {}}}", w.outcome, json::num(w.last().x), json::s3(&w.last().s), json::num(w.last().n_rel), w.steps.len());
        }
        "bisect" => {
            let fi: f64 = args[1].parse().unwrap();
            let (mut lo, mut hi) = (parse_dd(&args[2]), parse_dd(&args[3]));
            let vlo = verdict(fi, lo, &m, &opts);
            let vhi = verdict(fi, hi, &m, &opts);
            let (cl, ch) = (class_label(&vlo), class_label(&vhi));
            let mut n = 0;
            while (hi - lo).hi.abs() > 1e-13 * (1.0 + lo.hi.abs().max(hi.hi.abs())) && n < 80 {
                let mid = (lo + hi).mul_f(0.5);
                let vm = verdict(fi, mid, &m, &opts);
                if class_label(&vm) == cl {
                    lo = mid;
                } else {
                    hi = mid;
                }
                n += 1;
            }
            let (vl, vh) = (verdict(fi, lo, &m, &opts), verdict(fi, hi, &m, &opts));
            println!(
                "{{\"fi\": {}, \"class_lo\": \"{}\", \"class_hi\": \"{}\", \"final_lo\": \"{}\", \"final_hi\": \"{}\", \"iters\": {}, \"lo\": {}, \"hi\": {}, \"v_lo\": {}, \"v_hi\": {}}}",
                json::num(fi),
                cl,
                ch,
                class_label(&vl),
                class_label(&vh),
                n,
                json::dd(lo),
                json::dd(hi),
                json::verdict(&vl),
                json::verdict(&vh)
            );
        }
        "plant" => match args[1].as_str() {
            "blowup" => plant_blowup(&m),
            "linear" => plant_connect(false, &m),
            "cubic" => plant_connect(true, &m),
            "walk" => plant_walk(),
            other => panic!("plant {other}: blowup | linear | cubic"),
        },
        _ => eprintln!("commands: walk series ntable scan bisect plant"),
    }
}

/// B3: the locator on F'' = 2 F^3 and F F'' + F'^2 = 0, inward and outward from two starts each
fn plant_blowup(m: &Mutation) {
    let mut opts = WalkOpts::standard();
    opts.flip_side = m.plant_locator;
    let xstar = planted::dd(7321.0, 10000.0);
    let cases: [(&str, PolyN, &str); 2] = [("cube", planted::blowup_cube(), "F->inf"), ("sqrt", planted::blowup_sqrt(), "F->0")];
    for (name, poly, kind) in cases {
        let eq = Implicit { n: poly.clone(), ends: NoEnds, skip_order: None };
        for dir in [-1.0f64, 1.0] {
            // dir = -1: walk toward smaller x, the singularity at x*, the solution on x > x*
            for d in [DD::from_f64(0.9), DD::from_f64(0.35)] {
                let x0d = if dir < 0.0 { xstar + d } else { xstar - d };
                let s = if name == "cube" {
                    let f = DD::ONE / d;
                    let f1 = if dir < 0.0 { -f * f } else { f * f };
                    [f, f1, (f * f * f).mul_f(2.0)]
                } else {
                    let c = planted::dd(17.0, 10.0);
                    let f = (c * d).sqrt();
                    let f1 = if dir < 0.0 { c / f.mul_f(2.0) } else { -c / f.mul_f(2.0) };
                    let f2 = -(f1 * f1) / f;
                    [f, f1, f2]
                };
                let x0 = x0d.hi;
                let x1 = if dir < 0.0 { -50.0 } else { 50.0 };
                let w = walk(&eq, &poly, s, x0, x1, DD::ZERO, &opts);
                match &w.outcome {
                    Outcome::Singular { xstar: xs, bar, est, x_last } => println!(
                        "{{\"case\": \"{}\", \"dir\": {}, \"d\": {}, \"xstar_true\": {}, \"xstar\": {}, \"err\": {}, \"bar\": {}, \"kind\": \"{}\", \"est\": {}}}",
                        name,
                        dir,
                        json::num(d.hi),
                        json::dd(xstar),
                        json::num(*xs),
                        json::num((DD::from_f64(*xs) - xstar).hi.abs()),
                        json::num(*bar),
                        conical_kind(&poly, &w.last().s, *x_last, 1.0, Some(est)).replace("other", kind),
                        json::est(&Some(*est))
                    ),
                    o => println!("{{\"case\": \"{}\", \"dir\": {}, \"d\": {}, \"outcome\": \"{:?}\"}}", name, dir, json::num(d.hi), o),
                }
            }
        }
    }
}

/// B1, B2: the matching on the planted equation with the exact solution alpha + (a + b x)/(1 + x^2), a = 13/10,
/// b = 2/5, alpha = 7/10, from a perturbed start
fn plant_connect(cubic: bool, m: &Mutation) {
    let (alpha, a, b) = (planted::dd(7.0, 10.0), planted::dd(13.0, 10.0), planted::dd(2.0, 5.0));
    let pert: f64 = std::env::var("PLANT_PERT").ok().and_then(|v| v.parse().ok()).unwrap_or(0.2);
    let nodes = vec![-40.0, -10.0, -3.0, -1.0, 1.0, 3.0, 10.0, 40.0];
    let states: Vec<_> = nodes
        .iter()
        .enumerate()
        .map(|(i, &x)| {
            let mut s = planted::exact(alpha, a, b, DD::from_f64(x));
            for (f, v) in s.iter_mut().enumerate() {
                let sg = if (i + f) % 2 == 0 { 1.0 } else { -1.0 };
                *v = *v + DD::from_f64(pert * (abs_dd(*v).hi * 0.02 * sg * (1.0 + i as f64 / 8.0) + 0.001 * sg));
            }
            s
        })
        .collect();
    let inp = ConnectIn {
        nodes: nodes.clone(),
        match_node: 3,
        plan: Plan { kt: 30, tol: 1e-32, hmax: 5.0 },
        node_states: states,
        left: [a.mul_f(1.0 + 0.1 * pert), DD::from_f64(0.05 * pert)],
        right: b.mul_f(1.0 - 0.1 * pert),
        global: alpha + DD::from_f64(0.1),
        pin: alpha,
        tol: 1e-28,
    };
    let out = if cubic {
        connect(planted::planted_cubic(alpha, a, b, *m), &inp)
    } else {
        connect(planted::planted_linear(alpha, a, b, *m), &inp)
    };
    match out {
        Err(e) => println!("{{\"case\": \"{}\", \"error\": \"{}\"}}", if cubic { "cubic" } else { "linear" }, e),
        Ok(o) => {
            let mut node_err = 0.0f64;
            for (i, &x) in nodes.iter().enumerate() {
                let ex = planted::exact(alpha, a, b, DD::from_f64(x));
                for f in 0..3 {
                    node_err = node_err.max(abs_dd(o.nodes[i][f] - ex[f]).hi);
                }
            }
            let fl: Vec<String> = o.fnorm.iter().map(|x| json::num(*x)).collect();
            println!(
                "{{\"case\": \"{}\", \"converged\": {}, \"fnorm\": [{}], \"err_cm\": {}, \"err_e\": {}, \"err_cp\": {}, \"err_global\": {}, \"err_nodes\": {}, \"steps\": {:?}}}",
                if cubic { "cubic" } else { "linear" },
                o.converged,
                fl.join(", "),
                json::num(abs_dd(o.left[0] - a).hi),
                json::num(abs_dd(o.left[1]).hi),
                json::num(abs_dd(o.right - b).hi),
                json::num(abs_dd(o.global - alpha).hi),
                json::num(node_err),
                o.steps
            );
        }
    }
}

/// diagnostic: walk along the exact planted solution (linear and cubic) and report the largest deviation
fn plant_walk() {
    let (alpha, a, b) = (planted::dd(7.0, 10.0), planted::dd(13.0, 10.0), planted::dd(2.0, 5.0));
    for cubic in [false, true] {
        let poly = if cubic { planted::cubic_poly(alpha) } else { planted::linear_poly(alpha) };
        let eq = Implicit { n: poly.clone(), ends: NoEnds, skip_order: None };
        let mut opts = WalkOpts::standard();
        opts.plan = Plan { kt: 30, tol: 1e-32, hmax: 5.0 };
        let s0 = planted::exact(alpha, a, b, DD::from_f64(-40.0));
        let w = walk(&eq, &poly, s0, -40.0, 40.0, alpha, &opts);
        let mut err = 0.0f64;
        let mut nrel = 0.0f64;
        for st in &w.steps {
            let ex = planted::exact(alpha, a, b, DD::from_f64(st.x));
            for f in 0..3 {
                err = err.max(abs_dd(st.s[f] - ex[f]).hi);
            }
            nrel = nrel.max(st.n_rel);
        }
        println!("{{\"cubic\": {}, \"steps\": {}, \"max_err\": {}, \"n_rel_max\": {}, \"outcome\": \"{:?}\"}}", cubic, w.steps.len(), json::num(err), json::num(nrel), w.outcome);
    }
}
