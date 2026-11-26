//! loglog: command line for conical_loglog.py (#184). Walks of the implicit equation N = 0 over double-, quad- and
//! octo-double (`hp.rs`), conical or planted, printing JSON lines on stdout.
//!
//!   loglog run --poly NAME --prec dd|qd|od --mode project|conserve (--fwd F F1 | --x0 X --s0 "p,..;p,..;p,..")
//!              --to X1 [--cp a,b,c] [--kt K] [--tol T] [--hscale S] [--trace] [--max-steps N] [--mutate=NAME]
//!   NAME: conical | pl2 | pl3 | pl4 | pl5 (the planted equations of conical_loglog.py, see `planted_poly`)
//!   --fwd F F1: start at eta = 0 with F(0) = F, F'(0) = F1, F''(0) = -F1^2/F (the forward map of #167)
//!   --s0: the start state as three groups of f64 parts (hi first), sums exact in the precision
//! Mutations: eq-coeff (N of the table), recip-dd (constants rounded to double-double), no-project, plant-c, n-drift.
use conical_global::hp::{
    abs_s, hi, hp_walk, lane_f, lane_i, sc_of, scalar_of_parts, HpCfg, HpOut, PolyQ, Stop,
};
use conical_global::model::Mutation;
use multishoot::multidouble::precision::{Double, Octo, Precision, Quad};
use multishoot::multidouble::simd::Native;
use multishoot::Fl;
use std::collections::HashMap;

fn gcd(a: i128, b: i128) -> i128 {
    if b == 0 {
        a.abs()
    } else {
        gcd(b, a % b)
    }
}

type Term = (u8, u8, u8, u8, i64, i64);

fn terms_of(p: &PolyQ) -> Vec<Term> {
    p.groups
        .iter()
        .flat_map(|g| g.xp.iter().map(move |&(m, n, d)| (g.a as u8, g.b as u8, g.c as u8, m as u8, n, d)))
        .collect()
}

fn mulq(a: &PolyQ, b: &PolyQ) -> PolyQ {
    let (x, y) = (terms_of(a), terms_of(b));
    let mut out: Vec<Term> = Vec::new();
    for &(a1, b1, c1, m1, n1, d1) in &x {
        for &(a2, b2, c2, m2, n2, d2) in &y {
            let (n, d) = (n1 as i128 * n2 as i128, d1 as i128 * d2 as i128);
            let g = gcd(n, d).max(1);
            out.push((a1 + a2, b1 + b2, c1 + c2, m1 + m2, (n / g) as i64, (d / g) as i64));
        }
    }
    PolyQ::new(&out)
}

fn addq(a: &PolyQ, b: &PolyQ) -> PolyQ {
    let mut t = terms_of(a);
    t.extend(terms_of(b));
    PolyQ::new(&t)
}

fn scaleq(a: &PolyQ, num: i64) -> PolyQ {
    PolyQ::new(&terms_of(a).into_iter().map(|(a1, b1, c1, m, n, d)| (a1, b1, c1, m, n * num, d)).collect::<Vec<_>>())
}

/// the planted equations; c of pl3 is 1/3 (1003/3000 under the control plant-c)
fn planted_poly(name: &str, m_plant_c: bool) -> PolyQ {
    match name {
        // q + p^2 = 0: eta^2 F'' + eta F' + eta^2 F'^2
        "pl2" => PolyQ::new(&[(0, 0, 1, 2, 1, 1), (0, 1, 0, 1, 1, 1), (0, 2, 0, 2, 1, 1)]),
        // p^2 w' (2 p^2 - w') eta^8 with w' = q + p^2 + 2 c s (1 - p) + c^2 s^2, s = eta^-2:
        // W = eta^4 w' = eta^6 F'' + eta^5 F' + eta^6 F'^2 + 2c eta^2 - 2c eta^3 F' + c^2; B = 2 eta^6 F'^2 - W
        "pl3" => {
            let (cn, cd) = if m_plant_c { (1003i64, 3000i64) } else { (1, 3) };
            let (c2n, c2d) = (cn * cn, cd * cd);
            let w = PolyQ::new(&[
                (0, 0, 1, 6, 1, 1),
                (0, 1, 0, 5, 1, 1),
                (0, 2, 0, 6, 1, 1),
                (0, 0, 0, 2, 2 * cn, cd),
                (0, 1, 0, 3, -2 * cn, cd),
                (0, 0, 0, 0, c2n, c2d),
            ]);
            let b = addq(&PolyQ::new(&[(0, 2, 0, 6, 2, 1)]), &scaleq(&w, -1));
            let p2 = PolyQ::new(&[(0, 2, 0, 2, 1, 1)]);
            mulq(&mulq(&p2, &w), &b)
        }
        // (1 + eta^2) F'' + 4 eta F' + 2 (F - alpha), alpha = 7/10
        "pl4" => PolyQ::new(&[
            (0, 0, 1, 0, 1, 1),
            (0, 0, 1, 2, 1, 1),
            (0, 1, 0, 1, 4, 1),
            (1, 0, 0, 0, 2, 1),
            (0, 0, 0, 0, -7, 5),
        ]),
        // eta^6 F'^4 F'' + 4 eta^5 F'^5 (= p^4 (q + 3p))
        "pl5" => PolyQ::new(&[(0, 4, 1, 6, 1, 1), (0, 5, 0, 5, 4, 1)]),
        other => panic!("unknown planted equation {other}"),
    }
}

fn parts_json<P: Precision<Native>>(x: P::Scalar) -> String {
    let v: Vec<String> = (0..P::PARTS).map(|q| format!("{:e}", P::part(x, q))).collect();
    format!("[{}]", v.join(", "))
}

fn s_json<P: Precision<Native>>(s: &[P::Scalar; 3]) -> String {
    format!("\"F\": {}, \"F1\": {}, \"F2\": {}", parts_json::<P>(s[0]), parts_json::<P>(s[1]), parts_json::<P>(s[2]))
}

fn parse_parts(s: &str) -> Vec<f64> {
    s.split(',').map(|t| t.trim().parse().unwrap()).collect()
}

fn go<P: Precision<Native>>(a: &HashMap<String, Vec<String>>, flag: &dyn Fn(&str) -> bool, m: &Mutation, plant_c: bool) {
    let get = |k: &str| a.get(k).map(|v| v.as_slice());
    let poly_name = get("poly").map(|v| v[0].clone()).unwrap_or_else(|| "conical".into());
    let poly = if poly_name == "conical" { PolyQ::conical(m) } else { planted_poly(&poly_name, plant_c) };
    let dd = matches!(P::PARTS, 2);
    let (kt_def, tol_def) = match P::PARTS {
        2 => (30usize, 1e-32f64),
        4 => (60, 1e-64),
        _ => (120, 1e-124),
    };
    let _ = dd;
    let kt: usize = get("kt").map(|v| v[0].parse().unwrap()).unwrap_or(kt_def);
    let tol: f64 = get("tol").map(|v| v[0].parse().unwrap()).unwrap_or(tol_def);
    let mut cfg = HpCfg::new(kt, tol);
    let mode = get("mode").map(|v| v[0].clone()).unwrap_or_else(|| "project".into());
    cfg.project = mode == "project" && !flag("no-project");
    cfg.recip_dd = flag("recip-dd");
    cfg.skip = if m.n_drift { Some(3) } else { None };
    if let Some(v) = get("hscale") {
        cfg.hscale = v[0].parse().unwrap();
    }
    if let Some(v) = get("max-steps") {
        cfg.max_steps = v[0].parse().unwrap();
    }
    let to: f64 = get("to").map(|v| v[0].parse().unwrap()).expect("--to");
    let cps: Vec<f64> = get("cp").map(|v| v[0].split(',').map(|t| t.parse().unwrap()).collect()).unwrap_or_default();
    let (x0, s0): (f64, [P::Scalar; 3]) = if let Some(v) = get("fwd") {
        let (f, f1): (f64, f64) = (v[0].parse().unwrap(), v[1].parse().unwrap());
        let f2 = lane_f::<P>(f1) * lane_f::<P>(f1) / lane_f::<P>(f);
        (0.0, [sc_of::<P>(lane_f::<P>(f)), sc_of::<P>(lane_f::<P>(f1)), -sc_of::<P>(f2)])
    } else {
        let x0: f64 = get("x0").map(|v| v[0].parse().unwrap()).expect("--x0");
        let g: Vec<Vec<f64>> = get("s0").expect("--s0")[0].split(';').map(parse_parts).collect();
        (x0, [scalar_of_parts::<P>(&g[0]), scalar_of_parts::<P>(&g[1]), scalar_of_parts::<P>(&g[2])])
    };
    let want_trace = flag("trace");
    let o: HpOut<P> = hp_walk::<P>(&poly, s0, x0, to, &cps, &cfg, want_trace);
    for c in &o.cps {
        println!(
            "{{\"cp\": {:e}, \"step\": {}, \"n_rel\": {:e}, {}}}",
            c.x,
            c.step,
            c.n_rel,
            s_json::<P>(&c.s)
        );
    }
    for t in &o.trace {
        println!(
            "{{\"tr\": [{:e}, {:e}, {:e}, {:e}, {:e}, {:e}, {:e}]}}",
            t.x, t.s_hi[0], t.s_hi[1], t.s_hi[2], t.w_rel, t.n_rel, t.h
        );
    }
    let stop = match o.stop {
        Stop::End => "End",
        Stop::NonFinite => "NonFinite",
        Stop::TinyStep => "TinyStep",
        Stop::StepLimit => "StepLimit",
    };
    println!(
        "{{\"end\": true, \"stop\": \"{}\", \"x\": {:e}, \"steps\": {}, \"n_rel_max\": {:e}, \"n_rel_pre_max\": {:e}, \"kt\": {}, \"tol\": {:e}, \"project\": {}, {}}}",
        stop,
        o.x_end,
        o.steps,
        o.n_rel_max,
        o.n_rel_pre_max,
        kt,
        tol,
        cfg.project,
        s_json::<P>(&o.s_end)
    );
    let _ = (abs_s::<P>, hi::<P>, lane_i::<P>, Fl::<<P as Precision<Native>>::V, 1>::zero);
}

fn main() {
    let raw: Vec<String> = std::env::args().skip(1).collect();
    let mut m = Mutation::default();
    let mut flags: Vec<String> = Vec::new();
    let mut a: HashMap<String, Vec<String>> = HashMap::new();
    let mut i = 0;
    let cmd = raw.first().cloned().unwrap_or_default();
    while i < raw.len() {
        let t = &raw[i];
        if let Some(n) = t.strip_prefix("--mutate=") {
            match n {
                "eq-coeff" | "n-drift" => m = Mutation::from_name(n).unwrap(),
                "recip-dd" | "no-project" | "plant-c" => flags.push(n.to_string()),
                _ => panic!("unknown mutation {n}"),
            }
            i += 1;
        } else if t == "--trace" {
            flags.push("trace".into());
            i += 1;
        } else if let Some(k) = t.strip_prefix("--") {
            let mut v = Vec::new();
            i += 1;
            while i < raw.len() && !raw[i].starts_with("--") {
                v.push(raw[i].clone());
                i += 1;
            }
            a.insert(k.to_string(), v);
        } else {
            i += 1;
        }
    }
    let flag = |n: &str| flags.iter().any(|f| f == n);
    if cmd != "run" {
        eprintln!("commands: run");
        return;
    }
    let prec = a.get("prec").map(|v| v[0].clone()).unwrap_or_else(|| "dd".into());
    let plant_c = flag("plant-c");
    match prec.as_str() {
        "dd" => go::<Double>(&a, &flag, &m, plant_c),
        "qd" => go::<Quad>(&a, &flag, &m, plant_c),
        "od" => go::<Octo>(&a, &flag, &m, plant_c),
        p => panic!("--prec {p}"),
    }
}
