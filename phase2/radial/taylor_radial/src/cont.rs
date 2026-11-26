//! Pseudo-arclength continuation of the star family by multisolver, on multishoot's layout as a family
//! (`Shooting::in_global`): the alpha_c-pinned layout (StarZ) follows the phantom branch in alpha_c; the k-pinned layout
//! (StarK) follows it in k = sqrt(1 - w^2), whose fold (w_min) it passes and locates.
use crate::bvp::*;
use multidouble::dd::DD;
use multishoot::multisolver::arclength::{continuation, locate_fold, start, ArcOpts, FoldOpts};
use multishoot::multisolver::arclength::Bordered;
use multishoot::multisolver::newton::{exact, Exact, NewtonOpts};
use multishoot::multisolver::prec::parts;
use multishoot::{Double, Family, Shooting};

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
    ArcOpts { ds, ds_min: 1e-9, ds_max, max_points, newton: NewtonOpts { tol: 1e-26, max_iter: 25, ..Default::default() }, ..Default::default() }
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
}

pub fn follow_k(sol: &StarSol, dir: f64, ds: f64, ds_max: f64, max_points: usize, stop: impl Fn(f64, usize) -> bool) -> Result<KRun, String> {
    let fam = sol.sh.in_global(0, 0)?;
    let (zf, p) = fam.family_z::<Double>(&parts(&sol.z));
    let factor = exact_f;
    let st = start::<Double, _, _>(&fam, &parts(&zf), p, dir, false, factor)?;
    let o = opts(ds, ds_max, max_points);
    let mut count = 0usize;
    let path = continuation::<Double, _, _>(&fam, st, &o, factor, |pt| {
        count += 1;
        stop(Double_hi(pt.p), count)
    })?;
    let to_pp = |pt: &multishoot::multisolver::arclength::Point<Double>| PathPoint {
        z: fam.layout_z::<Double>(&parts(&pt.z), pt.p).unwrap(),
        p: <Double as multishoot::Precision<multishoot::Native>>::to_dd(pt.p),
    };
    let mut run = KRun { points: path.points.iter().map(to_pp).collect(), folds: path.folds.clone(), fold_k: None, fold_z: None, fold_iters: 0 };
    if let Some(&i) = path.folds.first() {
        let fopts = FoldOpts { tol: 1e-14, max_iter: 60 };
        let (pt, its) = locate_fold::<Double, _, _>(&fam, &path.points[i], path.ds[i], &fopts, &o.newton, factor)?;
        run.fold_k = Some(<Double as multishoot::Precision<multishoot::Native>>::to_dd(pt.p));
        run.fold_z = Some(fam.layout_z::<Double>(&parts(&pt.z), pt.p).unwrap());
        run.fold_iters = its;
    }
    let _ = Shooting::<StarK, 1>::new;
    let _: Option<&dyn Family<Double>> = None;
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
) -> Result<Vec<FamRow>, String> {
    let mut rows: Vec<FamRow> = vec![];
    let mut cur = start.alpha_c;
    let mut info = StarInfo { eps: start.eps, cfg: start.cfg.clone(), steps: start.steps.clone(), k: start.k, omega: start.omega, pc: start.pc, zc: start.zc, alpha_c: start.alpha_c, mass: start.mass, r_t: start.r_t, nodes: start.nodes, fnorm: start.fnorm.clone(), converged: true };
    let mut next_target = 0usize;
    // targets at or below the start's alpha_c are solved from the start
    let mut chunk_no = 0usize;
    let mut dsc = ds;
    while cur.to_f64() < alpha_stop {
        let t0 = std::time::Instant::now();
        chunk_no += 1;
        let sz = solve_z(cfg, &info.profile(), cur, None, false)?;
        {
            let fi = sz.info();
            let fl = *fi.fnorm.last().unwrap();
            if fl > 1e-18 {
                let w = worst_rows_z(&sz, 4);
                let nodes: Vec<f64> = sz.sh.nodes.clone();
                return Err(format!("chunk {chunk_no}: start at alpha_c {} not converged: {:?} worst rows {:?} x_node {} nodes {:?}", cur.to_f64(), fi.fnorm, w, sz.x_node, nodes));
            }
        }
        let run = follow_z(&sz, 1.0, dsc, ds_max, chunk, |_| false)?;
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
            rows.push(FamRow { alpha_c: tg, info: ti, n });
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
    Ok(rows)
}
