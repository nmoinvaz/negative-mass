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
        let zr = one / (one - c::<T>(2.0) * mr * inv_r);
        vec![pr, u, mr, zr]
    }
}

/// The layout with its nodes and the pin of k, no steps planned yet
pub fn layout(eps: f64, nodes: Vec<f64>, kt: usize, k_target: DD) -> Shooting<StarK, 1> {
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
use multishoot::multisolver::newton::{NewtonOpts, NewtonReport};
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

/// Solves the BVP at omega from the crate's shot; Newton in double-double.
pub fn solve(cfg: &StarCfg, omega: DD, verbose: bool) -> Result<StarSol, String> {
    let (_nodes, mut z, x_node, mut sh) = guess(cfg, omega);
    sh.kt = cfg.kt;
    let plan = Plan { kt: cfg.kt, tol: cfg.plan_tol, hmax: 4.0 };
    sh.plan_segments(&z[0], &z[1], &plan);
    set_row_scale(&mut sh, &z);
    let opts = NewtonOpts { tol: 1e-27, max_iter: 25, ..Default::default() };
    if std::env::var("TR_VERBOSE").is_ok() {
        eprintln!("solve: K {} nodes, len {}, planned segments steps {:?}", sh.k(), sh.len(), sh.segs.iter().map(|s| s.len() - 1).collect::<Vec<_>>());
    }
    let rep = sh.newton::<Double>(&mut z, &opts)?;
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
