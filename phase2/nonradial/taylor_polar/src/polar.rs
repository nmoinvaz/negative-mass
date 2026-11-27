//! D(s) = det[Y1 Y2 Y3 Z1 Z2 Z3] at r_m, with every 3-plane orthonormalized, as polar_solver.py builds it.
//!
//!   interior   3 regular solutions from r0 with the flat data r^l, Taylor steps with Gram-Schmidt after each;
//!   exterior   the decaying Zerilli solution through the Riccati variable q = f Z'/Z (inward from r_far) and the
//!              Moncrief map; the two decaying free-KG solutions from the contracting eigenspace of the frozen
//!              scalar block at R_s, integrated inward with the same stepper.
//! Gram-Schmidt with positive R diagonal replaces LAPACK's QR and its phase bookkeeping: the determinant of the
//! orthonormalized planes is then the real D of polar_solver.py (see plane.rs).
use crate::num::Num;
use crate::plane::*;
use crate::ser::*;
use crate::table::*;
use std::time::Instant;

#[derive(Clone, Debug)]
pub struct Cfg {
    pub l: usize,
    pub r0: f64,
    pub rm: f64,
    pub ctl: Ctl,
    pub model: Model,
}

#[derive(Clone, Debug)]
pub struct DOut<T: Num> {
    pub d: T,
    pub q: T,
    /// the interior Zerilli log-derivative f Z'/Z at r_m (`interior_q`), and the residual of its null vector
    pub q_int: T,
    pub qint_res: f64,
    pub r_far: f64,
    pub r_s: f64,
    pub steps_y: Stats,
    pub steps_q: usize,
    pub steps_s: Stats,
    pub t_y: f64,
    pub t_q: f64,
    pub t_s: f64,
    pub lnvol_y: T,
}

pub fn det_n<T: Num, const N: usize>(a: &[[T; N]; N]) -> T {
    let mut m = *a;
    let mut det = T::one();
    for c in 0..N {
        let mut p = c;
        for r in c + 1..N {
            if m[r][c].mag() > m[p][c].mag() {
                p = r;
            }
        }
        if m[p][c].mag() == 0.0 {
            return T::zero();
        }
        if p != c {
            m.swap(p, c);
            det = -det;
        }
        det = det * m[c][c];
        let inv = m[c][c].recip();
        for r in c + 1..N {
            let f = m[r][c] * inv;
            for k in c..N {
                m[r][k] = m[r][k] - f * m[c][k];
            }
        }
    }
    det
}

pub fn inv_n<T: Num, const N: usize>(a: &[[T; N]; N]) -> [[T; N]; N] {
    let mut m = *a;
    let mut e = [[T::zero(); N]; N];
    for i in 0..N {
        e[i][i] = T::one();
    }
    for c in 0..N {
        let mut p = c;
        for r in c + 1..N {
            if m[r][c].mag() > m[p][c].mag() {
                p = r;
            }
        }
        m.swap(p, c);
        e.swap(p, c);
        let inv = m[c][c].recip();
        for k in 0..N {
            m[c][k] = m[c][k] * inv;
            e[c][k] = e[c][k] * inv;
        }
        for r in 0..N {
            if r != c {
                let f = m[r][c];
                for k in 0..N {
                    m[r][k] = m[r][k] - f * m[c][k];
                    e[r][k] = e[r][k] - f * e[c][k];
                }
            }
        }
    }
    e
}

fn fro<T: Num, const N: usize>(a: &[[T; N]; N]) -> f64 {
    let mut s = 0.0;
    for r in a {
        for x in r {
            s += x.mag() * x.mag();
        }
    }
    s.sqrt()
}

/// The projector (I - sign M)/2 onto the invariant subspace of M with Re lambda < 0, by the scaled Newton
/// iteration X <- (mu X + X^-1 / mu)/2 for the matrix sign function.
pub fn contracting_projector<T: Num>(m: &[[T; 4]; 4]) -> [[T; 4]; 4] {
    let mut x = *m;
    let half = T::from_f64(0.5);
    let mut scaled = true;
    for _ in 0..200 {
        let xi = inv_n(&x);
        let mu = if scaled { (fro(&xi) / fro(&x)).sqrt() } else { 1.0 };
        let (mu_t, imu) = (T::from_f64(mu), T::from_f64(1.0 / mu));
        let mut xn = x;
        let mut diff = 0.0;
        for i in 0..4 {
            for j in 0..4 {
                xn[i][j] = half * (mu_t * x[i][j] + imu * xi[i][j]);
                diff += (xn[i][j] - x[i][j]).mag().powi(2);
            }
        }
        let rel = diff.sqrt() / fro(&xn);
        x = xn;
        if rel < 1e-2 {
            scaled = false;
        }
        if rel < 1e-3 * T::EPS.sqrt() * 1e-3 && !scaled {
            // one more unscaled step doubles the digits again
            let xi = inv_n(&x);
            for i in 0..4 {
                for j in 0..4 {
                    x[i][j] = half * (x[i][j] + xi[i][j]);
                }
            }
            break;
        }
    }
    let mut p = [[T::zero(); 4]; 4];
    for i in 0..4 {
        for j in 0..4 {
            p[i][j] = half * ((if i == j { T::one() } else { T::zero() }) - x[i][j]);
        }
    }
    p
}

fn lam_l<T: Num>(l: usize) -> T {
    T::from_f64(((l - 1) * (l + 2)) as f64) / T::from_f64(2.0)
}

/// the Zerilli potential V(r) at one r
fn zerilli_v<T: Num>(r: T, mass: T, lam: T) -> T {
    let one = T::one();
    let f = one - (mass + mass) / r;
    let three_m = mass * T::from_f64(3.0);
    let n1 = lam * lam * r * r * ((lam + one) * r + three_m) + T::from_f64(9.0) * mass * mass * (lam * r + mass);
    let d = r * r * r * (lam * r + three_m) * (lam * r + three_m);
    (one + one) * f * n1 / d
}

/// q = f Z'/Z at rm of the decaying Zerilli solution: the Riccati equation dq/dr = (s^2 + V - q^2)/f integrated
/// inward from the asymptotic start, in r (the series about x0 has radius about x0).
pub fn zerilli_q<T: Num>(gamma: T, mass: T, l: usize, r_start: f64, rm: f64, order: usize, tol: f64) -> (T, usize) {
    let one = T::one();
    let lam = lam_l::<T>(l);
    let s2 = gamma * gamma;
    let mut q = -(s2 + zerilli_v(T::from_f64(r_start), mass, lam)).sqrt();
    // the decaying branch q ~ -s (polar_solver's rule for a complex s; never true for a real gamma > 0)
    if (q / gamma).re_f64() > 0.0 {
        q = -q;
    }
    let inv = recips::<T>();
    let three_m = mass * T::from_f64(3.0);
    let nine_m2 = T::from_f64(9.0) * mass * mass;
    let two = one + one;
    let mut x0 = r_start;
    let mut steps = 0usize;
    while x0 > rm {
        let r = Ser::var(T::from_f64(x0));
        let rinv = recip(&r);
        let f = adds(&scl(&rinv, -(mass + mass)), one);
        let a = recip(&f);
        let lr = scl(&r, lam);
        let n1 = add(&scl(&mul(&mul(&r, &r), &adds(&scl(&r, lam + one), three_m)), lam * lam), &scl(&adds(&lr, mass), nine_m2));
        let d1 = mul(&powi(&r, 3), &powi(&adds(&lr, three_m), 2));
        let w = scl(&mul(&n1, &recip(&d1)), two);
        let mut qc = [T::zero(); KMAX + 1];
        let mut qq = [T::zero(); KMAX + 1];
        qc[0] = q;
        for k in 0..order {
            qq[k] = {
                let mut t = T::zero();
                for i in 0..=k {
                    t = t + qc[i] * qc[k - i];
                }
                t
            };
            // coefficient k of A (s^2 - q^2) + W
            let mut s = a.get(k) * s2 + w.get(k);
            for i in 0..=k {
                s = s - a.get(i) * qq[k - i];
            }
            qc[k + 1] = s * inv[k + 1];
        }
        let q0 = qc[0].mag().max(1e-300);
        let mut h = x0 * 0.5;
        for k in [order - 1, order] {
            let e = qc[k].mag() / q0;
            if e > 0.0 {
                h = h.min((tol / e).powf(1.0 / k as f64));
            }
        }
        let x1 = if x0 - h <= rm { rm } else { x0 - h };
        let hh = T::from_f64(x1) - T::from_f64(x0);
        let mut v = qc[order];
        for k in (0..order).rev() {
            v = v * hh + qc[k];
        }
        q = v;
        x0 = x1;
        steps += 1;
    }
    (q, steps)
}

/// (H1/s, K) = T (Z, dZ/dr) at r, Schwarzschild of mass M
pub fn moncrief<T: Num>(r: T, mass: T, l: usize) -> [[T; 2]; 2] {
    let lam = lam_l::<T>(l);
    let one = T::one();
    let f = one - (mass + mass) / r;
    let three_m = mass * T::from_f64(3.0);
    let a = (lam * r * r - three_m * lam * r - three_m * mass) / ((r - mass - mass) * (lam * r + three_m));
    let c = (lam * (lam + one) * r * r + three_m * lam * r + T::from_f64(6.0) * mass * mass) / (r * r * (lam * r + three_m));
    [[a, r], [c, if crate::mutate::on("moncrief") { -f } else { f }]]
}

/// q_int = f Z'/Z at r_m of the interior solution that also matches the decaying scalar exterior
/// (`polar_solver.interior_zerilli`). y: the interior 3-plane at r_m (rows H1/s, K, F, F', G/s, G'/s); ys: the exterior scalar
/// 2-plane (its 4 scalar rows). The coefficients c with Y[2:] c in span S are the null vector of the 4 x 5 matrix [Y[2:] | -S],
/// n_j = (-1)^j det(M without column j) (the scipy path takes an SVD null vector instead). (H1/s, K) = Y[:2] c, then
/// (Z, Z') = T^-1 (H1/s, K) with T the Moncrief map. Returns q_int and max |M n| / max |n|.
pub fn interior_q<T: Num>(y: &[[T; 6]; 3], ys: &[[T; 4]; 2], mext: T, rm: T, l: usize) -> (T, f64) {
    let mut m = [[T::zero(); 5]; 4];
    for i in 0..4 {
        for c in 0..3 {
            m[i][c] = y[c][2 + i];
        }
        for k in 0..2 {
            m[i][3 + k] = -ys[k][i];
        }
    }
    let mut n = [T::zero(); 5];
    for j in 0..5 {
        let mut sub = [[T::zero(); 4]; 4];
        for i in 0..4 {
            let mut cc = 0;
            for k in 0..5 {
                if k != j {
                    sub[i][cc] = m[i][k];
                    cc += 1;
                }
            }
        }
        let d = det_n(&sub);
        // mutation qint_cofsign: the cofactor signs are not alternated
        n[j] = if j % 2 == 0 || crate::mutate::on("qint_cofsign") { d } else { -d };
    }
    let nmax = n.iter().map(|x| x.mag()).fold(0.0, f64::max);
    let mut res = 0.0f64;
    for i in 0..4 {
        let mut s = T::zero();
        for k in 0..5 {
            s = s + m[i][k] * n[k];
        }
        res = res.max(s.mag());
    }
    let mut met = [T::zero(); 2];
    for c in 0..3 {
        met[0] = met[0] + y[c][0] * n[c];
        met[1] = met[1] + y[c][1] * n[c];
    }
    let mon = moncrief(rm, mext, l);
    let det = mon[0][0] * mon[1][1] - mon[0][1] * mon[1][0];
    let z = (mon[1][1] * met[0] - mon[0][1] * met[1]) / det;
    let zp = (mon[0][0] * met[1] - mon[1][0] * met[0]) / det;
    // mutation qint_f: the factor f of q = f Z'/Z left out
    let f = if crate::mutate::on("qint_f") { T::one() } else { T::one() - (mext + mext) / rm };
    (f * zp / z, res / nmax)
}

/// the start radius of the scalar plane, R_s = r_m + 40 / min Re sqrt(1 + (s -+ i w)^2)
pub fn scalar_start(rm: f64, gamma: f64, omega: f64) -> f64 {
    // sqrt(1 + (g - i w)^2) = sqrt(1 + g^2 - w^2 - 2 i g w): real part sqrt((|z| + Re z)/2)
    let re = 1.0 + gamma * gamma - omega * omega;
    let im = -2.0 * gamma * omega;
    let kr = (((re * re + im * im).sqrt() + re) / 2.0).sqrt();
    rm + 40.0 / kr
}

pub fn far_start(rm: f64, gamma: f64, l: usize) -> f64 {
    if crate::mutate::on("qstart") {
        return rm + 1.0;
    }
    rm + 40.0 / gamma.abs() + 20.0 * l as f64 / gamma.abs() + 50.0
}

pub fn d_eval<T: Num>(cfg: &Cfg, table: &JTable<T>, mext: T, gamma: f64) -> DOut<T> {
    let r_s = scalar_start(cfg.rm, gamma, cfg.model.omega);
    let r_far = far_start(cfg.rm, gamma, cfg.l);
    d_eval_s(cfg, table, mext, T::from_f64(gamma), r_s, r_far)
}

/// D at the growth rate s (a lane value, complex on a complex lane), the scalar plane started at r_s and the Riccati
/// integration of the Zerilli solution at r_far.
pub fn d_eval_s<T: Num>(cfg: &Cfg, table: &JTable<T>, mext: T, g: T, r_s: f64, r_far: f64) -> DOut<T> {
    let l = cfg.l;
    let s2 = g * g;
    // interior
    let t0 = Instant::now();
    let r0 = T::from_f64(cfg.r0);
    let mut rl1 = T::one(); // r0^(l-1)
    for _ in 1..l {
        rl1 = rl1 * r0;
    }
    let rl = rl1 * r0;
    let lt = T::from_f64(l as f64);
    let mut y = [[T::zero(); 6]; 3];
    y[0][0] = rl * r0 / lt;
    y[0][1] = rl;
    y[1][2] = rl;
    y[1][3] = lt * rl1;
    y[2][4] = rl;
    y[2][5] = lt * rl1;
    let (lnvol_y, steps_y) = integrate::<T, 6, 3>(table, &mut y, cfg.r0, cfg.rm, s2, &cfg.ctl);
    let t_y = t0.elapsed().as_secs_f64();
    // exterior: metric
    let t0 = Instant::now();
    let mext = if crate::mutate::on("mext0") { T::zero() } else { mext };
    let rmt = T::from_f64(cfg.rm);
    let (q, steps_q) = zerilli_q(g, mext, l, r_far, cfg.rm, cfg.ctl.order, cfg.ctl.tol);
    let f = T::one() - (mext + mext) / rmt;
    let mon = moncrief(rmt, mext, l);
    let zv = [T::one(), q / f];
    let m = [mon[0][0] * zv[0] + mon[0][1] * zv[1], mon[1][0] * zv[0] + mon[1][1] * zv[1]];
    let t_q = t0.elapsed().as_secs_f64();
    // exterior: scalar plane
    let t0 = Instant::now();
    let vac = VacSource { mass: mext, model: cfg.model, l, s: g, lo: 2 };
    // the metric rows must not feed on the scalar columns in vacuum
    {
        let jf = vac.full(r_s);
        for i in 0..2 {
            for j in 2..6 {
                for k in 0..4 {
                    assert!(jf[i][j].get(k).mag() == 0.0, "J_vac[{i}][{j}] couples the blocks");
                }
            }
        }
    }
    let mut jets = Jets::<T>::new(&vac.pattern());
    vac.jets(r_s, s2, 1, &mut jets);
    let mut js = [[T::zero(); 4]; 4];
    for en in &jets.ent {
        if en.n > 0 {
            js[en.i][en.j] = en.c[0];
        }
    }
    let p = contracting_projector(&js);
    let mut ys = [[T::zero(); 4]; 2];
    for i in 0..4 {
        let (c0, c1) = if crate::mutate::on("swapcols") { (2, 0) } else { (0, 2) };
        ys[0][i] = p[i][c0];
        ys[1][i] = p[i][c1];
    }
    let (_, steps_s) = integrate::<T, 4, 2>(&vac, &mut ys, r_s, cfg.rm, s2, &cfg.ctl);
    let t_s = t0.elapsed().as_secs_f64();
    let (q_int, qint_res) = interior_q(&y, &ys, mext, rmt, l);
    // assemble
    let mut z = [[T::zero(); 6]; 3];
    z[0][0] = m[0];
    z[0][1] = m[1];
    for i in 0..4 {
        z[1][2 + i] = ys[0][i];
        z[2][2 + i] = ys[1][i];
    }
    gram_schmidt(&mut z);
    let mut mat = [[T::zero(); 6]; 6];
    for c in 0..3 {
        for i in 0..6 {
            mat[i][c] = y[c][i];
            mat[i][3 + c] = z[c][i];
        }
    }
    let d = det_n(&mat);
    DOut { d, q, q_int, qint_res, r_far, r_s, steps_y, steps_q, steps_s, t_y, t_q, t_s, lnvol_y }
}
