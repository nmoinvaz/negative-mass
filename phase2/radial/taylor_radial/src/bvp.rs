//! The star of background.py as a multiple-shooting boundary-value problem on multishoot's layout.
//!
//! Fields (P, u, m, Z) = (phi0^2, phi0'/phi0, m(r), alpha^-2), as `star_ref.rs`. The regular start at the first node
//! rs is the leading-order series of the central data (P_c, Z_c) (its error feeds the singular solution at O(rs^5));
//! the far end R carries u(R) = -k - 1/R and Z(R) (1 - 2 m(R)/R) = 1, with the data (P_R, m_R). The global is
//! k = sqrt(1 - w^2), so that every equation is polynomial (w^2 = 1 - k^2) and the tail condition is linear in it.
//! One extra row pins k. The 8 lanes of the layout all hold the same star (POINTWISE).
use crate::mutate;
use multishoot::lane::Lane;
use multishoot::{cs, Equations, ExtraRow, Fl, Shooting, State, Work, DD};

pub const NF: usize = 4;

pub struct StarK {
    pub eps: f64,
}

fn c<T: Lane>(x: f64) -> Fl<T, 1> {
    Fl::splat(DD::from_f64(x))
}

/// The Taylor recursion of (P, u, m, Z) about x0: orders 1..=kt from order 0 in `w.c` (the recursion of
/// star_ref.rs's `star_coeffs` over fields). w2 = w^2 as a field.
fn recursion<T: Lane>(eps: f64, x0: f64, w2: Fl<T, 1>, w: &mut Work<T, 1>) {
    let n = w.kt;
    let nn = n + 2;
    let z = || Fl::<T, 1>::zero();
    let one = c::<T>(1.0);
    let two = c::<T>(2.0);
    let half = c::<T>(0.5);
    let lam = c::<T>(crate::model::LAM);
    let quarter_lam = c::<T>(crate::model::LAM / 4.0);
    let mu = Fl::<T, 1>::splat(crate::model::mu_dd());
    let mu6 = mu / c::<T>(6.0);
    let pi = Fl::<T, 1>::splat(crate::model::pi_dd());
    let fpie = pi.mulf(4.0 * eps);
    let epie = pi.mulf(8.0 * eps);
    let x0f = c::<T>(x0);
    let r0inv = one / x0f;
    let x0sq = x0f * x0f;
    let mut rinv = vec![z(); nn];
    let mut pw = r0inv;
    for (k, x) in rinv.iter_mut().enumerate() {
        *x = if k % 2 == 0 { pw } else { -pw };
        pw = pw * r0inv;
    }
    let (mut bb, mut aa, mut pz, mut uu, mut uup, mut uupb) = (vec![z(); nn], vec![z(); nn], vec![z(); nn], vec![z(); nn], vec![z(); nn], vec![z(); nn]);
    let (mut p2, mut p3, mut rho, mut pr, mut apr, mut arho) = (vec![z(); nn], vec![z(); nn], vec![z(); nn], vec![z(); nn], vec![z(); nn], vec![z(); nn]);
    let (mut am1, mut nup, mut lamp, mut g, mut diff) = (vec![z(); nn], vec![z(); nn], vec![z(); nn], vec![z(); nn], vec![z(); nn]);
    for k in 0..n {
        let d0 = if k == 0 { one } else { z() };
        let mr = cs(&w.c[2], &rinv, k);
        bb[k] = d0 - two * mr;
        if k == 0 {
            aa[0] = one / bb[0];
        } else {
            let mut s = z();
            for i in 1..=k {
                s = s + bb[i] * aa[k - i];
            }
            aa[k] = -(s * aa[0]);
        }
        pz[k] = cs(&w.c[0], &w.c[3], k);
        uu[k] = cs(&w.c[1], &w.c[1], k);
        uup[k] = cs(&uu, &w.c[0], k);
        uupb[k] = cs(&uup, &bb, k);
        p2[k] = cs(&w.c[0], &w.c[0], k);
        p3[k] = cs(&p2, &w.c[0], k);
        let v = half * w.c[0][k] - quarter_lam * p2[k] + mu6 * p3[k];
        let kin = half * (w2 * pz[k] + uupb[k]);
        rho[k] = kin + v;
        pr[k] = kin - v;
        let mut r2rho = x0sq * rho[k];
        if k >= 1 {
            r2rho = r2rho + two * x0f * rho[k - 1];
        }
        if k >= 2 {
            r2rho = r2rho + rho[k - 2];
        }
        let mp = fpie * r2rho;
        apr[k] = cs(&aa, &pr, k);
        arho[k] = cs(&aa, &rho, k);
        let mut rapr = x0f * apr[k];
        let mut rarho = x0f * arho[k];
        if k >= 1 {
            rapr = rapr + apr[k - 1];
            rarho = rarho + arho[k - 1];
        }
        am1[k] = aa[k] - d0;
        let am1r = cs(&am1, &rinv, k);
        nup[k] = am1r + epie * rapr;
        lamp[k] = -am1r + epie * rarho;
        g[k] = two * rinv[k] + half * nup[k] - half * lamp[k];
        let gu = cs(&g, &w.c[1], k);
        diff[k] = (d0 - lam * w.c[0][k] + mu * p2[k]) - w2 * w.c[3][k];
        let adz = cs(&aa, &diff, k);
        let up = -uu[k] - gu + adz;
        let pp = two * cs(&w.c[1], &w.c[0], k);
        let zp = -cs(&nup, &w.c[3], k);
        let ik = T::splat(DD::ONE) / T::splat(DD::from_f64((k + 1) as f64));
        w.c[0][k + 1] = pp.mul_t(ik);
        w.c[1][k + 1] = up.mul_t(ik);
        w.c[2][k + 1] = mp.mul_t(ik);
        w.c[3][k + 1] = zp.mul_t(ik);
    }
}

impl Equations<1> for StarK {
    const NF: usize = 4;
    const NL: usize = 2;
    const NR: usize = 2;
    const NG: usize = 1;
    const POINTWISE: bool = true;
    type Prep = ();

    fn taylor_coeffs<T: Lane>(&self, x0: f64, g: &[T], w: &mut Work<T, 1>) {
        let k = Fl([g[0]]);
        let w2 = c::<T>(1.0) - k * k;
        recursion(self.eps, x0, w2, w);
    }

    /// (P_c, u = r (1 - lam P_c + mu P_c^2 - w^2 Z_c) / 3, m = (4 pi / 3) eps rho_c r^3, Z_c) at r = x
    fn left_state<T: Lane>(&self, x: f64, data: &[Fl<T, 1>], g: &[T], _prep: &()) -> State<T, 1> {
        let (pc, zc) = (data[0], data[1]);
        let k = Fl([g[0]]);
        let w2 = c::<T>(1.0) - k * k;
        let one = c::<T>(1.0);
        let lam = c::<T>(crate::model::LAM);
        let mu = Fl::<T, 1>::splat(crate::model::mu_dd());
        let pi = Fl::<T, 1>::splat(crate::model::pi_dd());
        let pc2 = pc * pc;
        let cc = (one - lam * pc + mu * pc2) - w2 * zc;
        let v = c::<T>(0.5) * pc - c::<T>(crate::model::LAM / 4.0) * pc2 + (mu / c::<T>(6.0)) * pc2 * pc;
        let rho_c = c::<T>(0.5) * w2 * pc * zc + v;
        let r = c::<T>(x);
        let three = c::<T>(3.0);
        let ucoef = if mutate::on("start_term") { r * cc / c::<T>(2.0) } else { r * cc / three };
        let m = pi.mulf(4.0 * self.eps) * rho_c * r * r * r / three;
        vec![pc, ucoef, m, zc]
    }

    /// (P_R, u = -k - 1/R, m_R, Z = 1 / (1 - 2 m_R / R)) at r = x = R
    fn right_state<T: Lane>(&self, x: f64, data: &[Fl<T, 1>], g: &[T], _prep: &()) -> State<T, 1> {
        let (pr, mr) = (data[0], data[1]);
        let k = Fl([g[0]]);
        let r = c::<T>(x);
        let one = c::<T>(1.0);
        let inv_r = one / r;
        let u = if mutate::on("tailk") { -k + inv_r } else { -k - inv_r };
        let zr = if mutate::on("znorm") { one } else { one / (one - c::<T>(2.0) * mr * inv_r) };
        vec![pr, u, mr, zr]
    }
}

/// The layout with its nodes and the pin of k, no steps planned yet
pub fn layout(eps: f64, nodes: Vec<f64>, kt: usize, k_target: DD) -> Shooting<StarK, 1> {
    let eps = if mutate::on("eps") { -eps } else { eps };
    let mut sh = Shooting::new(StarK { eps }, nodes, kt);
    let og = sh.off_global();
    sh.extra.push(ExtraRow { terms: vec![(og, DD::ONE)], rhs: k_target });
    sh
}

/// Index helpers: node i, field f, lane j
pub fn zi(sh: &Shooting<StarK, 1>, i: usize, f: usize, j: usize) -> usize {
    sh.off_node(i) + f * 8 + j
}

pub fn zset(z: &mut [Vec<f64>], idx: usize, v: DD) {
    z[0][idx] = v.hi;
    z[1][idx] = v.lo;
}

pub fn zget(z: &[Vec<f64>], idx: usize) -> DD {
    DD::new(z[0][idx], z[1][idx])
}

// ------------------------------------------------------------------------------------------------------------
// the start of Newton's method from the crate's own shot (star_ref.rs), and the solve
// ------------------------------------------------------------------------------------------------------------
use crate::num::Num;
use crate::star_ref::{march, solve_central, SolveCfg, StarPars};
use multishoot::multisolver::newton::{NewtonOpts, NewtonReport, Scaling};
use multishoot::{Double, Plan};

#[derive(Clone, Debug)]
pub struct StarCfg {
    pub eps: f64,
    pub rs: f64,
    /// largest spacing of nodes inside the star (the tail's is 10 / k at most)
    pub dmax: f64,
    pub kt: usize,
    /// Taylor tolerance of the planned steps
    pub plan_tol: f64,
    /// P below which the tail is vacuum (the BVP ends 40 / k beyond)
    pub p_tail: f64,
    /// P / P_c at which the node layout turns from outward to inward segments
    pub x_ratio: f64,
}

impl StarCfg {
    pub fn new(eps: f64) -> Self {
        StarCfg { eps, rs: 1e-8, dmax: 8.0, kt: 30, plan_tol: 1e-32, p_tail: 1e-44, x_ratio: 1e-9 }
    }
}

pub struct StarSol {
    pub cfg: StarCfg,
    pub sh: Shooting<StarK, 1>,
    pub z: Vec<Vec<f64>>,
    /// index of the node where the layout turns to inward segments
    pub x_node: usize,
    pub report: Option<NewtonReport>,
}

pub fn make_nodes(rs: f64, r_x: f64, r_t: f64, k: f64, dmax: f64) -> Vec<f64> {
    let rnd = |x: f64| (x * 8.0).round() / 8.0;
    let r_x = rnd(r_x);
    let r_t = rnd(r_t);
    let mut v = vec![rs];
    for f in [1.0, 2.5, 4.5] {
        if f < r_x - 1.0 {
            v.push(f);
        }
    }
    let mut cur = *v.last().unwrap();
    let cap_in = dmax.min(2.5 / k);
    while cur + cap_in < r_x - 0.5 * cap_in {
        cur = rnd(cur + cap_in);
        v.push(cur);
    }
    if *v.last().unwrap() < r_x {
        v.push(r_x);
    }
    let mut cur = r_x;
    loop {
        let sp = dmax.max(0.5 * (cur - r_x)).min(10.0 / k);
        if cur + 1.5 * sp >= r_t {
            break;
        }
        cur = rnd(cur + sp);
        v.push(cur);
    }
    v.push(r_t);
    v
}

/// The nodes, the node states and the data from the crate's own shot at frequency omega.
pub fn guess(cfg: &StarCfg, omega: DD) -> (Vec<f64>, Vec<Vec<f64>>, usize, Shooting<StarK, 1>) {
    use multidouble::dd::DD as D;
    let k = (D::ONE - omega * omega).sqrt();
    let p = StarPars::<D> { w: omega, eps: D::from_f64(cfg.eps), lam: D::from_f64(crate::model::LAM), mu: crate::model::mu_dd(), pi: crate::model::pi_dd() };
    let scfg = SolveCfg { rs: cfg.rs, tol: 1e-30, r_norm: 32.0, r_class: 60.0, verbose: false };
    let zg = 1.0 - 0.02 * cfg.eps;
    let t_dbg = std::time::Instant::now();
    let (pc, zc, _it) = solve_central::<D>(&p, zg, &scfg);
    if std::env::var("TR_VERBOSE").is_ok() {
        eprintln!("guess: solve_central {:.1}s pc {:e} zc {:e}", t_dbg.elapsed().as_secs_f64(), pc.to_f64(), zc.to_f64());
    }
    // where the star falls to x_ratio P_c: outward states at integer radii, one radius at a time (the shot's error mode
    // takes over beyond that)
    let mut st: Vec<[D; 4]> = vec![];
    let mut pre: Vec<f64> = vec![];
    let mut cur = crate::star_ref::central_state(pc, zc, cfg.rs, &p);
    let mut x_from = cfg.rs;
    let mut r_x = 60.0;
    for i in 1..=60 {
        let r = i as f64;
        let s1 = march::<D>(cur, x_from, &[r], &p, 1e-30);
        cur = s1[0];
        x_from = r;
        st.push(cur);
        pre.push(r);
        if cur[0].to_f64() < cfg.x_ratio * pc.to_f64() {
            r_x = r;
            break;
        }
    }
    if std::env::var("TR_VERBOSE").is_ok() {
        eprintln!("guess: pre-march done {:.1}s r_x {}", t_dbg.elapsed().as_secs_f64(), r_x);
    }
    let kf = k.to_f64();
    let r_x = (r_x * 8.0).ceil() / 8.0;
    let px = st[(r_x as usize).saturating_sub(1).min(st.len() - 1)][0].to_f64();
    let r_cut = r_x + (px / cfg.p_tail).ln() / (2.0 * kf);
    let r_t = r_cut + 40.0 / kf;
    let nodes = make_nodes(cfg.rs, r_x, r_t, kf, cfg.dmax);
    let x_node = nodes.iter().position(|&x| (x - (r_x * 8.0).round() / 8.0).abs() < 1e-9).unwrap_or(nodes.len() - 1);
    // outward states at nodes <= r_x
    if std::env::var("TR_VERBOSE").is_ok() {
        eprintln!("guess: nodes {:?} x_node {} r_cut {:.2} r_t {:.2}", nodes, x_node, r_cut, r_t);
    }
    let out_t: Vec<f64> = nodes[..=x_node].to_vec();
    let out = march::<D>(crate::star_ref::central_state(pc, zc, cfg.rs, &p), cfg.rs, &out_t[1..], &p, 1e-30);
    // inward from r_t: P scaled to match at r_x
    let rr = *nodes.last().unwrap();
    let sx = out.last().copied().unwrap_or(crate::star_ref::central_state(pc, zc, cfg.rs, &p));
    let rx = nodes[x_node];
    let m_x = sx[2];
    let p_guess = sx[0] * (-(D::from_f64(2.0) * k * D::from_f64(rr - rx))).exp() * D::from_f64(rx * rx / (rr * rr));
    let zr = (D::ONE - D::from_f64(2.0) * m_x / D::from_f64(rr)).recip();
    let u_r = -k - D::ONE / D::from_f64(rr);
    let in_t: Vec<f64> = nodes[x_node..nodes.len() - 1].iter().rev().cloned().collect();
    let inn = march::<D>([p_guess, u_r, m_x, zr], rr, &in_t, &p, 1e-30);
    // inn[last] is the state at r_x: rescale P
    let scale = sx[0] / inn.last().unwrap()[0];
    let mut sh = layout(cfg.eps, nodes.clone(), cfg.kt, k);
    sh.match_node = x_node;
    let mut z = vec![vec![0.0; sh.len()]; 2];
    let big_k = nodes.len();
    for j in 0..8 {
        zset(&mut z, j, pc);
        zset(&mut z, 8 + j, zc);
        for i in 0..big_k {
            let s: [D; 4] = if i == 0 {
                crate::star_ref::central_state(pc, zc, cfg.rs, &p)
            } else if i <= x_node {
                out[i - 1]
            } else if i == big_k - 1 {
                [p_guess * scale, u_r, m_x, zr]
            } else {
                // nodes after x_node: inn is in descending order of in_t = nodes[x_node..K-1] reversed
                let idx = (big_k - 2) - i; // in_t[idx] = nodes[K-2-idx]
                let mut s = inn[idx];
                s[0] = s[0] * scale;
                s
            };
            for f in 0..4 {
                zset(&mut z, zi(&sh, i, f, j), s[f]);
            }
        }
        let o = sh.off_right();
        zset(&mut z, o + j, p_guess * scale);
        zset(&mut z, o + 8 + j, m_x);
    }
    zset(&mut z, sh.off_global(), k);
    (nodes, z, x_node, sh)
}

/// power-of-two row factors that bring each block's P row to O(1): the residual of P is tiny in the tail
pub fn set_row_scale(sh: &mut Shooting<StarK, 1>, z: &[Vec<f64>]) {
    let rows = sh.rows();
    let mut rs = vec![1.0; rows];
    for b in 0..=sh.k() {
        let target = sh.block_target(b);
        let pval = z[0][target].abs();
        if pval > 0.0 {
            let e = -(pval.log2().round());
            let f = 2f64.powi(e as i32);
            for j in 0..8 {
                rs[sh.row_block(b) + j] = f;
            }
        }
    }
    sh.row_scale = rs;
}


/// Newton, and when it stops short of the tolerance a second pass with the P rows scaled by the iterate rather than by the guess
fn newton_twice<E: Equations<1>>(sh: &mut Shooting<E, 1>, z: &mut Vec<Vec<f64>>, opts: &NewtonOpts) -> Result<NewtonReport, String> {
    let mut o2 = *opts;
    if mutate::on("nonewton") {
        o2.max_iter = 1;
    }
    let opts = &o2;
    let mut rep = sh.newton::<Double>(z, opts)?;
    if !rep.converged {
        let mut rs2 = vec![1.0; sh.rows()];
        for b in 0..=sh.k() {
            let pval = z[0][sh.block_target(b)].abs();
            if pval > 0.0 {
                let f = 2f64.powi(-(pval.log2().round() as i32));
                for j in 0..8 {
                    rs2[sh.row_block(b) + j] = f;
                }
            }
        }
        sh.row_scale = rs2;
        let r2 = sh.newton::<Double>(z, opts)?;
        let mut all = rep.fnorm.clone();
        all.extend(r2.fnorm.iter().cloned());
        rep = r2;
        rep.fnorm = all;
    }
    Ok(rep)
}

/// Solves the BVP at omega from the crate's shot; Newton in double-double.
pub fn solve(cfg: &StarCfg, omega: DD, verbose: bool) -> Result<StarSol, String> {
    let (_nodes, mut z, x_node, mut sh) = guess(cfg, omega);
    sh.kt = cfg.kt;
    let plan = Plan { kt: cfg.kt, tol: cfg.plan_tol, hmax: 4.0 };
    sh.plan_segments(&z[0], &z[1], &plan);
    set_row_scale(&mut sh, &z);
    let opts = NewtonOpts { tol: 1e-28, max_iter: 25, scaling: Scaling::Iterates { floor: 1e-250 }, ..Default::default() };
    if std::env::var("TR_VERBOSE").is_ok() {
        eprintln!("solve: K {} nodes, len {}, planned segments steps {:?}", sh.k(), sh.len(), sh.segs.iter().map(|s| s.len() - 1).collect::<Vec<_>>());
    }
    let rep = newton_twice(&mut sh, &mut z, &opts)?;
    if verbose {
        eprintln!("newton: converged {} fnorm {:?}", rep.converged, rep.fnorm.iter().map(|v| format!("{v:.1e}")).collect::<Vec<_>>());
    }
    Ok(StarSol { cfg: cfg.clone(), sh, z, x_node, report: Some(rep) })
}

impl StarSol {
    pub fn k(&self) -> DD {
        zget(&self.z, self.sh.off_global())
    }
    pub fn omega(&self) -> DD {
        let k = self.k();
        (DD::ONE - k * k).sqrt()
    }
    pub fn pc(&self) -> DD {
        zget(&self.z, 0)
    }
    pub fn zc(&self) -> DD {
        zget(&self.z, 8)
    }
    pub fn mass(&self) -> DD {
        zget(&self.z, self.sh.off_right() + 8)
    }
}

/// the rows of the largest scaled residual: (block, field, lane, value)
pub fn worst_rows(sol: &StarSol, n: usize) -> Vec<(usize, usize, usize, f64)> {
    use multishoot::multisolver::prec::parts;
    let f = sol.sh.residual::<Double>(&parts(&sol.z), false).unwrap();
    let mut v: Vec<(usize, usize, usize, f64)> = vec![];
    for b in 0..=sol.sh.k() {
        for fl in 0..NF {
            for j in 0..8 {
                let i = sol.sh.row_block(b) + fl * 8 + j;
                v.push((b, fl, j, f[i].to_f64().abs()));
            }
        }
    }
    v.sort_by(|a, b| b.3.partial_cmp(&a.3).unwrap());
    v.truncate(n);
    v
}

// ------------------------------------------------------------------------------------------------------------
// the alpha_c-pinned form: fields (P, u, m, Z, k) with k' = 0, the global alpha_c
// ------------------------------------------------------------------------------------------------------------

/// Fields (P, u, m, Z, k): k = sqrt(1 - w^2) is carried as a constant field, the central lapse alpha_c is the global
/// (Z_c = alpha_c^-2). Left data (P_c, k_c), right data (P_R, m_R, k_R); NL + NR = NF, and the one global has one extra row
/// (alpha_c = target), so the layout is square and `Shooting::in_global` takes alpha_c as the family's parameter.
pub struct StarZ {
    pub eps: f64,
}

impl Equations<1> for StarZ {
    const NF: usize = 5;
    const NL: usize = 2;
    const NR: usize = 3;
    const NG: usize = 1;
    const POINTWISE: bool = true;
    type Prep = ();

    fn taylor_coeffs<T: Lane>(&self, x0: f64, _g: &[T], w: &mut Work<T, 1>) {
        let k = w.c[4][0];
        let w2 = c::<T>(1.0) - k * k;
        recursion(self.eps, x0, w2, w);
        for j in 1..w.c[4].len() {
            w.c[4][j] = Fl::zero();
        }
    }

    fn left_state<T: Lane>(&self, x: f64, data: &[Fl<T, 1>], g: &[T], _prep: &()) -> State<T, 1> {
        let (pc, kc) = (data[0], data[1]);
        let one = c::<T>(1.0);
        let alpha = Fl([g[0]]);
        let zc = one / (alpha * alpha);
        let w2 = one - kc * kc;
        let lam = c::<T>(crate::model::LAM);
        let mu = Fl::<T, 1>::splat(crate::model::mu_dd());
        let pi = Fl::<T, 1>::splat(crate::model::pi_dd());
        let pc2 = pc * pc;
        let cc = (one - lam * pc + mu * pc2) - w2 * zc;
        let v = c::<T>(0.5) * pc - c::<T>(crate::model::LAM / 4.0) * pc2 + (mu / c::<T>(6.0)) * pc2 * pc;
        let rho_c = c::<T>(0.5) * w2 * pc * zc + v;
        let r = c::<T>(x);
        let three = c::<T>(3.0);
        let m = pi.mulf(4.0 * self.eps) * rho_c * r * r * r / three;
        vec![pc, r * cc / three, m, zc, kc]
    }

    fn right_state<T: Lane>(&self, x: f64, data: &[Fl<T, 1>], _g: &[T], _prep: &()) -> State<T, 1> {
        let (pr, mr, kr) = (data[0], data[1], data[2]);
        let one = c::<T>(1.0);
        let inv_r = one / c::<T>(x);
        let u = if mutate::on("tailk") { -kr + inv_r } else { -kr - inv_r };
        let zr = if mutate::on("znorm") { one } else { one / (one - c::<T>(2.0) * mr * inv_r) };
        vec![pr, u, mr, zr, kr]
    }
}

pub fn layout_z(eps: f64, nodes: Vec<f64>, kt: usize, alpha_target: DD) -> Shooting<StarZ, 1> {
    let eps = if mutate::on("eps") { -eps } else { eps };
    let mut sh = Shooting::new(StarZ { eps }, nodes, kt);
    let og = sh.off_global();
    sh.extra.push(ExtraRow { terms: vec![(og, DD::ONE)], rhs: alpha_target });
    sh
}

pub fn zi5(sh: &Shooting<StarZ, 1>, i: usize, f: usize, j: usize) -> usize {
    sh.off_node(i) + f * 8 + j
}

/// A star's profile on [r_s, r_t] by its Taylor steps, with a Yukawa tail beyond r_t: the guess's source for a new layout
pub struct Profile {
    pub steps: Vec<crate::table::StarStep<DD>>,
    pub rs: f64,
    pub r_t: f64,
    pub k: DD,
    pub eps: f64,
}

impl Profile {
    pub fn state(&self, r: f64) -> [DD; 4] {
        if r > self.r_t {
            let last = self.state(self.r_t);
            let p = last[0] * (-(DD::from_f64(2.0) * self.k * DD::from_f64(r - self.r_t))).exp() * DD::from_f64(self.r_t * self.r_t / (r * r));
            let zz = (DD::ONE - DD::from_f64(2.0) * last[2] / DD::from_f64(r)).recip();
            return [p, -self.k - DD::ONE / DD::from_f64(r), last[2], zz];
        }
        let idx = match self.steps.binary_search_by(|s| s.lo.partial_cmp(&r).unwrap()) {
            Ok(i) => i,
            Err(i) => i.saturating_sub(1),
        };
        let s = &self.steps[idx];
        let h = DD::from_f64(r) - DD::from_f64(s.x0);
        let mut out = [DD::ZERO; 4];
        for f in 0..4 {
            let mut v = s.coef[f][crate::table::NS];
            for k in (0..crate::table::NS).rev() {
                v = v * h + s.coef[f][k];
            }
            out[f] = v;
        }
        out
    }
}

pub struct StarSolZ {
    pub cfg: StarCfg,
    pub sh: Shooting<StarZ, 1>,
    pub z: Vec<Vec<f64>>,
    pub x_node: usize,
    pub report: Option<NewtonReport>,
}

/// The layout of a star near `prof`'s: nodes from the profile's own P (as `guess` does from the shot)
pub fn nodes_from_profile(cfg: &StarCfg, prof: &Profile) -> (Vec<f64>, usize) {
    let pc = prof.state(cfg.rs)[0].to_f64();
    let mut r_x = 1.0;
    let mut px = pc;
    let mut r = 1.0;
    while r < prof.r_t {
        let p = prof.state(r)[0].to_f64();
        if p < cfg.x_ratio * pc {
            r_x = r;
            px = p;
            break;
        }
        r += 1.0;
    }
    if px == pc {
        r_x = prof.r_t.min(60.0);
        px = prof.state(r_x)[0].to_f64();
    }
    let r_x = (r_x * 8.0).ceil() / 8.0;
    let kf = prof.k.to_f64();
    let r_cut = r_x + (px / cfg.p_tail).ln().max(0.0) / (2.0 * kf);
    let r_t = r_cut + 40.0 / kf;
    let nodes = make_nodes(cfg.rs, r_x, r_t, kf, cfg.dmax);
    let x_node = nodes.iter().position(|&x| (x - r_x).abs() < 1e-9).unwrap_or(nodes.len() - 1);
    (nodes, x_node)
}

/// Solves the alpha_c-pinned star at alpha_target from a profile (nodes from `nodes`, or from the profile's size)
pub fn solve_z(cfg: &StarCfg, prof: &Profile, alpha_target: DD, nodes_in: Option<(Vec<f64>, usize)>, verbose: bool) -> Result<StarSolZ, String> {
    let (nodes, x_node) = nodes_in.unwrap_or_else(|| nodes_from_profile(cfg, prof));
    let mut sh = layout_z(cfg.eps, nodes.clone(), cfg.kt, alpha_target);
    sh.match_node = x_node;
    let mut z = vec![vec![0.0; sh.len()]; 2];
    let big_k = nodes.len();
    let s0 = prof.state(cfg.rs);
    let sr = prof.state(*nodes.last().unwrap());
    for j in 0..8 {
        zset(&mut z, j, s0[0]);
        zset(&mut z, 8 + j, prof.k);
        for i in 0..big_k {
            let s = if i == 0 { s0 } else { prof.state(nodes[i]) };
            for f in 0..4 {
                zset(&mut z, zi5(&sh, i, f, j), s[f]);
            }
            zset(&mut z, zi5(&sh, i, 4, j), prof.k);
        }
        let o = sh.off_right();
        zset(&mut z, o + j, sr[0]);
        zset(&mut z, o + 8 + j, sr[2]);
        zset(&mut z, o + 16 + j, prof.k);
    }
    zset(&mut z, sh.off_global(), alpha_target);
    let plan = Plan { kt: cfg.kt, tol: cfg.plan_tol, hmax: 4.0 };
    sh.plan_segments(&z[0], &z[1], &plan);
    // row scale: P rows to O(1)
    let rows = sh.rows();
    let mut rs = vec![1.0; rows];
    for b in 0..=sh.k() {
        let target = sh.block_target(b);
        let pval = z[0][target].abs();
        if pval > 0.0 {
            let f = 2f64.powi(-(pval.log2().round() as i32));
            for j in 0..8 {
                rs[sh.row_block(b) + j] = f;
            }
        }
    }
    sh.row_scale = rs;
    let opts = NewtonOpts { tol: 1e-28, max_iter: 25, scaling: Scaling::Iterates { floor: 1e-250 }, ..Default::default() };
    let rep = newton_twice(&mut sh, &mut z, &opts)?;
    if verbose {
        eprintln!("newton_z: converged {} fnorm {:?}", rep.converged, rep.fnorm.iter().map(|v| format!("{v:.1e}")).collect::<Vec<_>>());
    }
    let sol = StarSolZ { cfg: cfg.clone(), sh, z, x_node, report: Some(rep) };
    if std::env::var("TR_VERBOSE").is_ok() && !sol.report.as_ref().unwrap().converged {
        eprintln!("solve_z not converged at alpha {}: worst rows {:?} nodes {:?}", alpha_target.to_f64(), worst_rows_z(&sol, 4), sol.sh.nodes);
    }
    Ok(sol)
}

impl StarSolZ {
    pub fn alpha_c(&self) -> DD {
        zget(&self.z, self.sh.off_global())
    }
    pub fn k(&self) -> DD {
        zget(&self.z, 8)
    }
    pub fn omega(&self) -> DD {
        let k = self.k();
        (DD::ONE - k * k).sqrt()
    }
    pub fn pc(&self) -> DD {
        zget(&self.z, 0)
    }
    pub fn zc(&self) -> DD {
        let a = self.alpha_c();
        (a * a).recip()
    }
    pub fn mass(&self) -> DD {
        zget(&self.z, self.sh.off_right() + 8)
    }
}

/// what the radial problem needs of a solved star
pub struct StarInfo {
    pub eps: f64,
    pub cfg: StarCfg,
    pub steps: Vec<crate::table::StarStep<DD>>,
    pub k: DD,
    pub omega: DD,
    pub pc: DD,
    pub zc: DD,
    pub alpha_c: DD,
    pub mass: DD,
    pub r_t: f64,
    pub nodes: usize,
    pub fnorm: Vec<f64>,
    pub converged: bool,
}

impl StarSol {
    pub fn info(&self) -> StarInfo {
        self.info_at(&self.z)
    }

    /// the info of another state z of the same layout (a point of an arclength path in k)
    pub fn info_at(&self, z: &[Vec<f64>]) -> StarInfo {
        let k = zget(z, self.sh.off_global());
        let g = [multishoot::lane::DV::splat(k)];
        let rep = self.report.as_ref().unwrap();
        let zc = zget(z, 8);
        StarInfo {
            eps: self.cfg.eps,
            cfg: self.cfg.clone(),
            steps: crate::table::steps_layout(&self.sh, z, &g),
            k,
            omega: (DD::ONE - k * k).sqrt(),
            pc: zget(z, 0),
            zc,
            alpha_c: zc.recip().sqrt(),
            mass: zget(z, self.sh.off_right() + 8),
            r_t: *self.sh.nodes.last().unwrap(),
            nodes: self.sh.k(),
            fnorm: rep.fnorm.clone(),
            converged: rep.converged,
        }
    }
}

impl StarSolZ {
    pub fn info(&self) -> StarInfo {
        self.info_at(&self.z)
    }

    /// the info of another state z of the same layout (a point of an arclength path); the Newton report is this solve's
    pub fn info_at(&self, z: &[Vec<f64>]) -> StarInfo {
        let alpha = zget(z, self.sh.off_global());
        let g = [multishoot::lane::DV::splat(alpha)];
        let rep = self.report.as_ref().unwrap();
        let k = zget(z, 8);
        StarInfo {
            eps: self.cfg.eps,
            cfg: self.cfg.clone(),
            steps: crate::table::steps_layout(&self.sh, z, &g),
            k,
            omega: (DD::ONE - k * k).sqrt(),
            pc: zget(z, 0),
            zc: (alpha * alpha).recip(),
            alpha_c: alpha,
            mass: zget(z, self.sh.off_right() + 8),
            r_t: *self.sh.nodes.last().unwrap(),
            nodes: self.sh.k(),
            fnorm: rep.fnorm.clone(),
            converged: rep.converged,
        }
    }
}

impl StarInfo {
    pub fn profile(&self) -> Profile {
        Profile { steps: self.steps.clone(), rs: self.cfg.rs, r_t: self.r_t, k: self.k, eps: self.eps }
    }
}

/// Solves the k-pinned star (omega target) from a profile of a neighboring star (continuation in omega)
pub fn solve_k_from_profile(cfg: &StarCfg, prof: &Profile, k_target: DD, verbose: bool) -> Result<StarSol, String> {
    let (nodes, x_node) = nodes_from_profile(cfg, prof);
    let mut sh = layout(cfg.eps, nodes.clone(), cfg.kt, k_target);
    sh.match_node = x_node;
    let mut z = vec![vec![0.0; sh.len()]; 2];
    let big_k = nodes.len();
    let s0 = prof.state(cfg.rs);
    let sr = prof.state(*nodes.last().unwrap());
    for j in 0..8 {
        zset(&mut z, j, s0[0]);
        zset(&mut z, 8 + j, s0[3]);
        for i in 0..big_k {
            let s = if i == 0 { s0 } else { prof.state(nodes[i]) };
            for f in 0..4 {
                zset(&mut z, zi(&sh, i, f, j), s[f]);
            }
        }
        let o = sh.off_right();
        zset(&mut z, o + j, sr[0]);
        zset(&mut z, o + 8 + j, sr[2]);
    }
    zset(&mut z, sh.off_global(), k_target);
    let plan = Plan { kt: cfg.kt, tol: cfg.plan_tol, hmax: 4.0 };
    sh.plan_segments(&z[0], &z[1], &plan);
    set_row_scale(&mut sh, &z);
    let opts = NewtonOpts { tol: 1e-28, max_iter: 25, scaling: Scaling::Iterates { floor: 1e-250 }, ..Default::default() };
    let rep = newton_twice(&mut sh, &mut z, &opts)?;
    if verbose {
        eprintln!("newton_k: converged {} fnorm {:?}", rep.converged, rep.fnorm.iter().map(|v| format!("{v:.1e}")).collect::<Vec<_>>());
    }
    Ok(StarSol { cfg: cfg.clone(), sh, z, x_node, report: Some(rep) })
}

/// the rows of the largest scaled residual of a StarZ solution: (block, field, lane, value)
pub fn worst_rows_z(sol: &StarSolZ, n: usize) -> Vec<(usize, usize, usize, f64)> {
    use multishoot::multisolver::prec::parts;
    let f = sol.sh.residual::<Double>(&parts(&sol.z), false).unwrap();
    let mut v: Vec<(usize, usize, usize, f64)> = vec![];
    for b in 0..=sol.sh.k() {
        for fl in 0..5 {
            for j in 0..1 {
                let i = sol.sh.row_block(b) + fl * 8 + j;
                v.push((b, fl, j, f[i].to_f64().abs()));
            }
        }
    }
    v.sort_by(|a, b| b.3.partial_cmp(&a.3).unwrap());
    v.truncate(n);
    v
}
