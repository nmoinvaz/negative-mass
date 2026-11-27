//! Pseudo-arclength continuation of the star family by multisolver, on multishoot's layout as a family
//! (`Shooting::in_global`): the alpha_c-pinned layout (StarZ) follows the phantom branch in alpha_c; the k-pinned layout
//! (StarK) follows it in k = sqrt(1 - w^2), whose fold (w_min) it passes and locates.
use crate::bvp::*;
use multidouble::dd::DD;
use multishoot::multisolver::arclength::{continuation, locate_fold, start, ArcOpts, FoldOpts};
use multishoot::multisolver::arclength::Bordered;
use multishoot::multisolver::newton::{exact, Exact, NewtonOpts};
use multishoot::multisolver::prec::parts;
use multishoot::{Double, Family, Plan, Shooting};

pub struct PathPoint {
    pub z: Vec<Vec<f64>>,
    pub p: DD,
}

pub struct ArcRun {
    pub points: Vec<PathPoint>,
    pub ds: Vec<f64>,
    pub folds: Vec<usize>,
    pub refused: usize,
    pub stalled: bool,
}

fn exact_f<Fa: Family<Double> + ?Sized>(b: &Bordered<'_, Double, Fa>, w: &[&[f64]]) -> Result<Exact<Double>, String> {
    exact::<Double, _>(b, w, false)
}

fn opts(ds: f64, ds_max: f64, max_points: usize) -> ArcOpts {
    ArcOpts { ds, ds_min: 1e-9, ds_max, max_points, newton: NewtonOpts { tol: 1e-26, max_iter: 25, scaling: multishoot::multisolver::newton::Scaling::Iterates { floor: 1e-250 }, ..Default::default() }, ..Default::default() }
}

/// The StarZ layout of `sol` followed in alpha_c until `stop(alpha_c)`, at most `max_points` points. `dir` = +1 / -1.
pub fn follow_z(sol: &StarSolZ, dir: f64, ds: f64, ds_max: f64, max_points: usize, stop: impl Fn(f64) -> bool) -> Result<ArcRun, String> {
    let fam = sol.sh.in_global(0, 0)?;
    let (zf, p) = fam.family_z::<Double>(&parts(&sol.z));
    let factor = exact_f;
    let st = start::<Double, _, _>(&fam, &parts(&zf), p, dir, false, factor)?;
    let path = continuation::<Double, _, _>(&fam, st, &opts(ds, ds_max, max_points), factor, |pt| stop(Double_hi(pt.p)))?;
    Ok(ArcRun {
        points: path
            .points
            .iter()
            .map(|pt| PathPoint { z: fam.layout_z::<Double>(&parts(&pt.z), pt.p).unwrap(), p: DD::new(<Double as multishoot::Precision<multishoot::Native>>::part(pt.p, 0), <Double as multishoot::Precision<multishoot::Native>>::part(pt.p, 1)) })
            .collect(),
        ds: path.ds.clone(),
        folds: path.folds.clone(),
        refused: path.refused,
        stalled: path.stalled,
    })
}

#[allow(non_snake_case)]
fn Double_hi(p: <Double as multishoot::Precision<multishoot::Native>>::Scalar) -> f64 {
    <Double as multishoot::Precision<multishoot::Native>>::to_dd(p).hi
}

/// the StarK layout followed in k from `sol` until `stop`, with the fold of k (the minimum of omega) located
pub struct KRun {
    pub points: Vec<PathPoint>,
    pub folds: Vec<usize>,
    pub fold_k: Option<DD>,
    pub fold_z: Option<Vec<Vec<f64>>>,
    pub fold_iters: usize,
    /// k at the two path points around the first sign change of t_k, and locate_fold's error when it failed
    pub fold_between: Option<(DD, DD)>,
    pub locate_err: Option<String>,
    pub last: multishoot::multisolver::arclength::Point<Double>,
    pub last_z: Vec<Vec<f64>>,
}

/// The StarK layout followed in k from `sol` (from `start_pt` when given, with its tangent: no orientation lost beside the
/// fold), with the fold of k located when the path crosses it
pub fn follow_k(sol: &StarSol, start_pt: Option<multishoot::multisolver::arclength::Point<Double>>, dir: f64, ds: f64, ds_max: f64, max_points: usize) -> Result<KRun, String> {
    let fam = sol.sh.in_global(0, 0)?;
    let factor = exact_f;
    let st = match start_pt {
        Some(p) => p,
        None => {
            let (zf, p) = fam.family_z::<Double>(&parts(&sol.z));
            start::<Double, _, _>(&fam, &parts(&zf), p, dir, false, factor)?
        }
    };
    let mut o = opts(ds, ds_max, max_points);
    o.newton.tol = 1e-23;
    let path = continuation::<Double, _, _>(&fam, st, &o, factor, |_| false)?;
    let to_pp = |pt: &multishoot::multisolver::arclength::Point<Double>| PathPoint {
        z: fam.layout_z::<Double>(&parts(&pt.z), pt.p).unwrap(),
        p: <Double as multishoot::Precision<multishoot::Native>>::to_dd(pt.p),
    };
    let lastp = path.points.last().unwrap().clone();
    let mut run = KRun { points: path.points.iter().map(to_pp).collect(), folds: path.folds.clone(), fold_k: None, fold_z: None, fold_iters: 0, fold_between: None, locate_err: None, last_z: fam.layout_z::<Double>(&parts(&lastp.z), lastp.p).unwrap(), last: lastp };
    if let Some(&i) = path.folds.first() {
        let kk = |pt: &multishoot::multisolver::arclength::Point<Double>| <Double as multishoot::Precision<multishoot::Native>>::to_dd(pt.p);
        run.fold_between = Some((kk(&path.points[i]), kk(&path.points[i + 1])));
        let fopts = FoldOpts { tol: 1e-14, max_iter: 60 };
        match locate_fold::<Double, _, _>(&fam, &path.points[i], path.ds[i], &fopts, &o.newton, factor) {
            Ok((pt, its)) => {
                run.fold_k = Some(<Double as multishoot::Precision<multishoot::Native>>::to_dd(pt.p));
                run.fold_z = Some(fam.layout_z::<Double>(&parts(&pt.z), pt.p).unwrap());
                run.fold_iters = its;
            }
            Err(e) => run.locate_err = Some(e),
        }
    }
    Ok(run)
}

/// One row of the phantom family: alpha_c, omega, N, M at a star solved at a stated alpha_c
pub struct FamRow {
    pub alpha_c: f64,
    pub info: StarInfo,
    pub n: DD,
}

pub struct ChunkLog {
    pub chunk: usize,
    pub points: usize,
    pub refused: usize,
    pub stalled: bool,
    pub alpha_end: f64,
    pub omega_end: f64,
    pub folds: Vec<usize>,
    pub t: f64,
}

/// The phantom family in alpha_c by arclength, in chunks of `chunk` points; after each chunk the star at its last point
/// is re-solved on a fresh layout (new nodes and Taylor plan), and every stated alpha_c reached in the chunk is solved from
/// the chunk's nearest path point. `rows` gets the stars at the targets, `track` every path point's (alpha_c, omega).
pub fn follow_family(
    cfg: &StarCfg,
    start: &StarInfo,
    targets: &[f64],
    alpha_stop: f64,
    chunk: usize,
    ds: f64,
    ds_max: f64,
    verbose: bool,
    track: &mut Vec<(f64, f64)>,
    logs: &mut Vec<ChunkLog>,
    reverse: bool,
    on_row: &mut dyn FnMut(&FamRow),
) -> Result<(), String> {
    let mut cur = start.alpha_c;
    let mut info = StarInfo { eps: start.eps, cfg: start.cfg.clone(), steps: start.steps.clone(), k: start.k, omega: start.omega, pc: start.pc, zc: start.zc, alpha_c: start.alpha_c, mass: start.mass, r_t: start.r_t, nodes: start.nodes, fnorm: start.fnorm.clone(), converged: true };
    let mut next_target = 0usize;
    // targets at or below the start's alpha_c are solved from the start
    let mut chunk_no = 0usize;
    let mut dsc = ds;
    while cur.to_f64() < alpha_stop {
        let t0 = std::time::Instant::now();
        chunk_no += 1;
        if reverse && chunk_no > 3 {
            return Err(format!("reversed continuation: alpha_c fell to {} after 3 chunks and never reached {}", cur.to_f64(), alpha_stop));
        }
        let sz = solve_z(cfg, &info.profile(), cur, None, false)?;
        {
            let fi = sz.info();
            let fl = *fi.fnorm.last().unwrap();
            if verbose {
                eprintln!("chunk {chunk_no} start: alpha_c {:.6} fnorm {:?}", cur.to_f64(), fi.fnorm.iter().map(|v| format!("{v:.1e}")).collect::<Vec<_>>());
            }
            if fl > 1e-18 {
                let w = worst_rows_z(&sz, 4);
                let nodes: Vec<f64> = sz.sh.nodes.clone();
                return Err(format!("chunk {chunk_no}: start at alpha_c {} not converged: {:?} worst rows {:?} x_node {} nodes {:?}", cur.to_f64(), fi.fnorm, w, sz.x_node, nodes));
            }
        }
        let run = follow_z(&sz, if reverse { -1.0 } else { 1.0 }, dsc, ds_max, chunk, |_| false)?;
        let last = run.points.last().unwrap();
        for pt in &run.points {
            let k = zget(&pt.z, 8);
            track.push((pt.p.to_f64(), (DD::ONE - k * k).sqrt().to_f64()));
        }
        // targets reached in this chunk
        while next_target < targets.len() && targets[next_target] <= last.p.to_f64() {
            let tg = targets[next_target];
            // the nearest point of this chunk (by alpha_c)
            let best = run.points.iter().min_by(|a, b| (a.p.to_f64() - tg).abs().partial_cmp(&(b.p.to_f64() - tg).abs()).unwrap()).unwrap();
            let prof = sz.info_at(&best.z).profile();
            let st = solve_z(cfg, &prof, DD::from_f64(tg), None, false)?;
            let ti = st.info();
            if *ti.fnorm.last().unwrap() > 1e-24 {
                return Err(format!("target {tg}: not converged {:?}", ti.fnorm));
            }
            let n = crate::family::noether(&ti.steps, ti.omega, cfg.rs, ti.pc, ti.zc);
            let row = FamRow { alpha_c: tg, info: ti, n };
            on_row(&row);
            next_target += 1;
        }
        let ai = sz.info_at(&last.z);
        cur = ai.alpha_c;
        info = ai;
        let dss = run.ds.iter().cloned().fold(0.0, f64::max);
        dsc = (dss * 0.5).max(ds).min(ds_max);
        logs.push(ChunkLog { chunk: chunk_no, points: run.points.len(), refused: run.refused, stalled: run.stalled, alpha_end: cur.to_f64(), omega_end: info.omega.to_f64(), folds: run.folds.clone(), t: t0.elapsed().as_secs_f64() });
        if verbose {
            eprintln!("chunk {chunk_no}: {} points, refused {}, stalled {}, alpha_c {:.6} omega {:.8} t {:.1}s", run.points.len(), run.refused, run.stalled, cur.to_f64(), info.omega.to_f64(), t0.elapsed().as_secs_f64());
        }
        if run.points.len() < 2 {
            return Err(format!("chunk {chunk_no}: no progress at alpha_c {}", cur.to_f64()));
        }
    }
    Ok(())
}

pub struct FoldResult {
    pub detected: Option<(DD, DD)>,
    pub k_fold: Option<DD>,
    pub alpha_c: Option<f64>,
    pub points: usize,
    pub iters: usize,
    pub locate_err: Option<String>,
}

/// The fold of the k-pinned family (the minimum of omega): arclength in k in chunks from the star `start`, the layout re-planned
/// between chunks from the last point (which, with its tangent, starts the next chunk: pinning k beside its maximum is
/// ill-posed). The sign change of t_k is the multisolver's detection; `locate_fold` then refines it, or reports why it could not.
pub fn fold_search(cfg: &StarCfg, start: &StarInfo, chunk: usize, ds: f64, ds_max: f64, max_chunks: usize, verbose: bool) -> Result<FoldResult, String> {
    let mut sol = solve_k_from_profile(cfg, &start.profile(), start.k, false)?;
    let mut pt = None;
    let mut points = 0usize;
    for c in 0..max_chunks {
        let run = follow_k(&sol, pt.take(), 1.0, ds, ds_max, chunk)?;
        points += run.points.len();
        if verbose {
            let lk = run.points.last().unwrap().p;
            eprintln!("fold chunk {c}: {} points, folds {:?}, k {:.10} omega {:.10}", run.points.len(), run.folds, lk.to_f64(), (DD::ONE - lk * lk).sqrt().to_f64());
        }
        if run.fold_between.is_some() {
            let alpha_c = run.fold_z.as_ref().map(|fz| zget(fz, 8).recip().sqrt().to_f64());
            return Ok(FoldResult { detected: run.fold_between, k_fold: run.fold_k, alpha_c, points, iters: run.fold_iters, locate_err: run.locate_err });
        }
        sol.z = run.last_z.clone();
        let plan = Plan { kt: cfg.kt, tol: cfg.plan_tol, hmax: 4.0 };
        sol.sh.plan_segments(&sol.z[0], &sol.z[1], &plan);
        set_row_scale(&mut sol.sh, &sol.z);
        pt = Some(run.last);
    }
    Err("no sign change of t_k within the chunks".to_string())
}
