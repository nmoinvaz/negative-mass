//! offgrid: command line for conical_offgrid.py (#183). Verdicts of the inward walk from a series start at radius up to
//! 819200, on the walker of `hp.rs`; JSON lines on stdout.
//!
//!   offgrid scan --prec dd|qd|od --mode conserve|project --fi A,B --c X,Y [--rstop R] [--lk N] [--kt K] [--tol T]
//!                [--hscale S] [--mutate=NAME]      one verdict per (fi, c): the start is the first of L = 200 4^k with
//!                                                  residual <= 1e-20, plus N radii (the second start of the spread check)
//!   offgrid series --prec P --fi F --c C --lk K     the series state at L = 200 4^K as sums of parts, and its residual
//!   offgrid walk --prec P --fi F --c C --to X       the state at X of the walk from the chosen start
//!   offgrid bisect --prec P --fi F --c CLO,CHI      bisection in C between two verdict labels
//!   offgrid ntable --prec P --file FILE             N and its scale at rational points
//! Mutations: eq-coeff, series-coeff, locator-side, recip-dd, start-bar, n-drift, halving.
use conical_global::hp::{hi, lane_rat, sc_of, EvalS, PolyQ};
use conical_global::model::Mutation;
use conical_global::offgrid::{apply_mutation, label, resid_of, series_state, verdict_p, walk_to, OgOpts, LS};
use multishoot::multidouble::precision::{Double, Octo, Precision, Quad};
use multishoot::multidouble::simd::Native;
use std::collections::HashMap;

fn parts_json<P: Precision<Native>>(x: P::Scalar) -> String {
    let v: Vec<String> = (0..P::PARTS).map(|q| format!("{:e}", P::part(x, q))).collect();
    format!("[{}]", v.join(", "))
}

fn go<P: Precision<Native>>(cmd: &str, a: &HashMap<String, String>, muts: &[String]) {
    let mut m = Mutation::default();
    let mut o = OgOpts::new(P::PARTS);
    for n in muts {
        assert!(apply_mutation(n, &mut m, &mut o), "unknown mutation {n}");
    }
    let get = |k: &str| a.get(k).cloned();
    if let Some(v) = get("mode") {
        o.cfg.project = v == "project";
    }
    if let Some(v) = get("rstop") {
        o.cfg.rstop = Some(v.parse().unwrap());
    }
    if get("rabs").is_some() {
        o.cfg.rstop_rel = false;
    }
    if let Some(v) = get("lk") {
        o.lk_extra = v.parse().unwrap();
    }
    if let Some(v) = get("kt") {
        o.cfg.kt = v.parse().unwrap();
    }
    if let Some(v) = get("tol") {
        o.cfg.tol = v.parse().unwrap();
    }
    if let Some(v) = get("hscale") {
        o.cfg.hscale = v.parse().unwrap();
    }
    if cmd == "ntable" {
        // N and the sum of |terms| at the rational points of a file ("a/b a/b a/b a/b" per line), in parts
        let ev = EvalS::<P>::new(&PolyQ::conical(&m));
        let txt = std::fs::read_to_string(get("file").expect("--file")).unwrap();
        let rat = |t: &str| -> P::Scalar {
            let (n, d) = t.split_once('/').unwrap_or((t, "1"));
            sc_of::<P>(lane_rat::<P>(n.parse().unwrap(), d.parse().unwrap()))
        };
        for line in txt.lines().filter(|l| !l.trim().is_empty()) {
            let v: Vec<P::Scalar> = line.split_whitespace().map(rat).collect();
            let (val, sm, _) = ev.eval(&[v[0], v[1], v[2]], v[3]);
            println!("{{\"n\": {}, \"scale\": {}}}", parts_json::<P>(val), parts_json::<P>(sm));
        }
        return;
    }
    let fis: Vec<f64> = get("fi").expect("--fi").split(',').map(|s| s.parse().unwrap()).collect();
    let cs: Vec<f64> = get("c").expect("--c").split(',').map(|s| s.parse().unwrap()).collect();
    match cmd {
        "scan" => {
            for &fi in &fis {
                for &c in &cs {
                    println!("{}", verdict_p::<P>(fi, c, &m, &o).to_json());
                }
            }
        }
        "series" => {
            let lk: usize = get("lk").expect("--lk").parse().unwrap();
            let ev = EvalS::<P>::new(&PolyQ::conical(&m));
            for &fi in &fis {
                for &c in &cs {
                    let s = series_state::<P>(fi, c, lk, &m);
                    let r = resid_of::<P>(&ev, &s, LS[lk]);
                    println!(
                        "{{\"fi\": {:e}, \"c\": {:e}, \"l\": {:e}, \"resid\": {:e}, \"F\": {}, \"F1\": {}, \"F2\": {}}}",
                        fi,
                        c,
                        LS[lk],
                        r,
                        parts_json::<P>(s[0]),
                        parts_json::<P>(s[1]),
                        parts_json::<P>(s[2])
                    );
                }
            }
        }
        "walk" => {
            // the state at --to of the walk from the chosen start (no locator unless --rstop)
            let x1: f64 = get("to").expect("--to").parse().unwrap();
            for &fi in &fis {
                for &c in &cs {
                    match walk_to::<P>(fi, c, x1, &m, &o) {
                        None => println!("{{\"fi\": {:e}, \"c\": {:e}, \"class\": \"NOSTART\"}}", fi, c),
                        Some((l, out)) => println!(
                            "{{\"fi\": {:e}, \"c\": {:e}, \"l\": {:e}, \"stop\": \"{:?}\", \"x\": {:e}, \"steps\": {}, \"F\": {}, \"F1\": {}, \"F2\": {}}}",
                            fi, c, l, out.stop, out.x_end, out.steps, parts_json::<P>(out.s_end[0]), parts_json::<P>(out.s_end[1]), parts_json::<P>(out.s_end[2])
                        ),
                    }
                }
            }
        }
        "bisect" => {
            // bisection in C between clo and chi (the two ends must carry different labels) to relative width 1e-12
            let (mut lo, mut hi_) = (cs[0], cs[1]);
            let fi = fis[0];
            let vlo = verdict_p::<P>(fi, lo, &m, &o);
            let vhi = verdict_p::<P>(fi, hi_, &m, &o);
            let (cl, ch) = (label(&vlo), label(&vhi));
            let mut n = 0;
            while (hi_ - lo).abs() > 1e-12 * (1.0 + lo.abs().max(hi_.abs())) && n < 200 && cl != ch {
                let mid = 0.5 * (lo + hi_);
                if label(&verdict_p::<P>(fi, mid, &m, &o)) == cl {
                    lo = mid;
                } else {
                    hi_ = mid;
                }
                n += 1;
            }
            let (vl, vh) = (verdict_p::<P>(fi, lo, &m, &o), verdict_p::<P>(fi, hi_, &m, &o));
            println!(
                "{{\"fi\": {:e}, \"class_lo\": \"{}\", \"class_hi\": \"{}\", \"final_lo\": \"{}\", \"final_hi\": \"{}\", \"iters\": {}, \"lo\": {:e}, \"hi\": {:e}, \"v_lo\": {}, \"v_hi\": {}}}",
                fi, cl, ch, label(&vl), label(&vh), n, lo, hi_, vl.to_json(), vh.to_json()
            );
        }
        _ => eprintln!("commands: scan series walk bisect ntable"),
    }
    let _ = (hi::<P>, sc_of::<P>);
}

fn main() {
    let raw: Vec<String> = std::env::args().skip(1).collect();
    let cmd = raw.first().cloned().unwrap_or_default();
    let mut a: HashMap<String, String> = HashMap::new();
    let mut muts: Vec<String> = Vec::new();
    let mut i = 1;
    while i < raw.len() {
        if let Some(n) = raw[i].strip_prefix("--mutate=") {
            muts.push(n.to_string());
            i += 1;
        } else if let Some(k) = raw[i].strip_prefix("--") {
            a.insert(k.to_string(), raw.get(i + 1).cloned().unwrap_or_default());
            i += 2;
        } else {
            i += 1;
        }
    }
    match a.get("prec").map(|s| s.as_str()).unwrap_or("dd") {
        "dd" => go::<Double>(&cmd, &a, &muts),
        "qd" => go::<Quad>(&cmd, &a, &muts),
        "od" => go::<Octo>(&cmd, &a, &muts),
        p => panic!("--prec {p}"),
    }
}
