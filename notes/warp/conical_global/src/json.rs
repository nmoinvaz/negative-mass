//! A minimal JSON writer: numbers, strings, arrays and objects built as text. Doubles go out as [hi, lo] pairs where the
//! low part matters.
use crate::scan::Verdict;
use crate::walk::{Est, LocStats, S3};
use multishoot::DD;

pub fn num(x: f64) -> String {
    if x.is_nan() || x.is_infinite() {
        "null".into()
    } else {
        format!("{x:e}")
    }
}

pub fn dd(x: DD) -> String {
    format!("[{}, {}]", num(x.hi), num(x.lo))
}

pub fn s3(s: &S3) -> String {
    format!("[{}, {}, {}]", dd(s[0]), dd(s[1]), dd(s[2]))
}

pub fn opt(x: Option<f64>) -> String {
    x.map(num).unwrap_or_else(|| "null".into())
}

pub fn est(e: &Option<Est>) -> String {
    match e {
        None => "null".into(),
        Some(e) => format!(
            "{{\"r\": {}, \"gamma\": {}, \"resid\": {}, \"field\": {}, \"side\": {}, \"alpha\": {}}}",
            num(e.r),
            num(e.gamma),
            num(e.resid),
            e.field,
            e.side,
            num(e.alpha())
        ),
    }
}

/// the locator's counts and the f64 fit's stop beside it; written with the feature `locate_lib` only
pub fn loc(l: &LocStats) -> String {
    format!(
        "{{\"calls\": {}, \"real\": {}, \"behind\": {}, \"pair\": {}, \"unresolved\": {}, \"nothing\": {}, \"pair_vs_real\": {}, \"old_real\": {}, \
         \"unres_old_real\": {}, \"step\": {}, \"old_step\": {}, \"old_xstar\": {}, \"old_bar\": {}, \"rel_bar\": {}, \"alpha_bar\": {}}}",
        l.calls,
        l.real,
        l.behind,
        l.pair,
        l.unresolved,
        l.nothing,
        l.pair_vs_real,
        l.old_real,
        l.unres_old_real,
        l.step,
        l.old_step,
        opt(l.old_xstar),
        num(l.old_bar),
        num(l.rel_bar),
        num(l.alpha_bar)
    )
}

/// `, "loc": {...}` with the feature, nothing without it (the output of the f64 locator is unchanged)
pub fn loc_key(l: &LocStats) -> String {
    if cfg!(feature = "locate_lib") {
        format!(", \"loc\": {}", loc(l))
    } else {
        String::new()
    }
}

pub fn verdict(v: &Verdict) -> String {
    let through = match &v.through {
        None => "null".to_string(),
        Some(t) => format!(
            "{{\"fi_minus\": {}, \"c_minus\": {}, \"mismatch\": {}, \"state_end\": {}}}",
            num(t.fi_minus),
            num(t.c_minus),
            num(t.mismatch),
            s3(&t.state_end)
        ),
    };
    let floor: Vec<String> = v.floor.iter().map(s3).collect();
    format!(
        "{{\"fi\": {}, \"c\": {}, \"class\": \"{}\", \"kind\": \"{}\", \"xstar\": {}, \"bar\": {}, \"est\": {}, \"x_last\": {}, \
         \"s_last\": {}, \"floor\": [{}], \"through\": {}, \"global\": {}, \"steps\": {}, \"n_rel_max\": {}, \"alpha_pred\": {}, \"start_l\": {}, \"n_rel_start\": {}, \"n_rel_200\": {}{}}}",
        num(v.fi),
        dd(v.c),
        v.class,
        v.kind,
        opt(v.xstar),
        opt(v.bar),
        est(&v.est),
        num(v.x_last),
        s3(&v.s_last),
        floor.join(", "),
        through,
        v.global,
        v.steps,
        num(v.n_rel_max),
        opt(v.alpha_pred),
        num(v.start_l),
        num(v.n_rel_start),
        num(v.n_rel_200),
        loc_key(&v.loc)
    )
}
