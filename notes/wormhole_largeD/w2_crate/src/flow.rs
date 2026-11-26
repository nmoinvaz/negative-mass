//! Taylor series of the system: the throat's Frobenius series, regular-point series, and the adaptive inward flow.
//! A state is (Y, w, dY/dOm, dw/dOm).
use crate::mutate;
use crate::num::Num;
use crate::sys::{Coef, Sys};

pub type V4<T> = [T; 4];

pub struct LaneCfg {
    pub name: &'static str,
    pub digits: f64,
    pub kt: usize,
    pub tol: f64,
    pub k0: usize,
    pub t: f64,
    /// T as declared, before a control changes it (the bound of S5 uses this one)
    pub t_decl: f64,
    /// stated digits (planted)
    pub d_lane: f64,
}

pub fn lane_cfg(name: &str) -> LaneCfg {
    let mut c = match name {
        "f64" => LaneCfg { name: "f64", digits: 16.0, kt: 24, tol: 1e-17, k0: 20, t: 51.0, t_decl: 51.0, d_lane: 11.0 },
        "dd" => LaneCfg { name: "dd", digits: 32.0, kt: 48, tol: 1e-33, k0: 46, t: 88.0, t_decl: 88.0, d_lane: 27.0 },
        "qd" => LaneCfg { name: "qd", digits: 64.0, kt: 96, tol: 1e-65, k0: 90, t: 162.0, t_decl: 162.0, d_lane: 56.0 },
        "od" => LaneCfg { name: "od", digits: 128.0, kt: 180, tol: 1e-130, k0: 170, t: 311.0, t_decl: 311.0, d_lane: 110.0 },
        _ => panic!("lane {name}"),
    };
    if mutate::on("short_series") {
        c.k0 /= 2;
    }
    if mutate::on("short_contraction") {
        c.t /= 4.0;
    }
    if mutate::on("step_loose") {
        c.tol = 1e-10;
    }
    c
}

fn mag<T: Num>(x: T) -> f64 {
    x.abs().to_f64()
}

/// The even mode's Frobenius series about l = 0 (Y(0) = 1): coefficients 0..=kk of Y, w, dY/dOm, dw/dOm
pub fn throat_series<T: Num>(sys: &Sys<T>, om: T, kk: usize) -> [Vec<T>; 4] {
    let c = sys.coef(T::zero(), T::one(), om, kk + 1);
    let eps = sys.eps;
    let z = T::zero();
    let mut y = vec![z; kk + 1];
    let mut w = vec![z; kk + 1];
    let mut yd = vec![z; kk + 1];
    let mut wd = vec![z; kk + 1];
    y[0] = T::one();
    w[0] = -(c.a[0] / c.b[0]);
    wd[0] = -(c.ad[0] / c.b[0]);
    for m in 1..=kk {
        y[m] = eps * w[m - 1] * sys.rk[m];
        yd[m] = eps * wd[m - 1] * sys.rk[m];
        let mut s = c.a[0] * y[m];
        let mut sd = c.a[0] * yd[m] + c.ad[0] * y[m];
        for j in 1..=m {
            s = s + c.a[j] * y[m - j] + c.b[j] * w[m - j];
            sd = sd + c.a[j] * yd[m - j] + c.ad[j] * y[m - j] + c.b[j] * wd[m - j];
        }
        for j in 1..m {
            let fac = c.f[j + 1] * T::from_f64((m - j) as f64);
            s = s - fac * w[m - j];
            sd = sd - fac * wd[m - j];
        }
        let inv = (T::from_f64(m as f64) - c.b[0]).recip();
        w[m] = s * inv;
        wd[m] = sd * inv;
    }
    [y, w, yd, wd]
}

/// The series summed at l = 1: the state and the size of the last two terms relative to it
pub fn throat<T: Num>(sys: &Sys<T>, om: T, kk: usize) -> (V4<T>, f64) {
    let ser = throat_series(sys, om, kk);
    let z = T::zero();
    let sum = |v: &Vec<T>| v.iter().fold(z, |a, &x| a + x);
    let st = [sum(&ser[0]), sum(&ser[1]), sum(&ser[2]), sum(&ser[3])];
    let scale = mag(st[0]).max(mag(st[1]));
    let (y, w) = (&ser[0], &ser[1]);
    let tail = [y[kk], w[kk], y[kk - 1], w[kk - 1]].iter().map(|&x| mag(x)).fold(0.0, f64::max) / scale;
    (st, tail)
}

/// Series of the state about a regular point l0 (F(l0) != 0), in tau with l = l0 + dh tau, orders 0..=kk
pub fn regular<T: Num>(sys: &Sys<T>, c: &Coef<T>, dh: T, v0: V4<T>, kk: usize) -> [Vec<T>; 4] {
    let z = T::zero();
    let mut y = vec![z; kk + 1];
    let mut w = vec![z; kk + 1];
    let mut yd = vec![z; kk + 1];
    let mut wd = vec![z; kk + 1];
    y[0] = v0[0];
    w[0] = v0[1];
    yd[0] = v0[2];
    wd[0] = v0[3];
    let inv_f0 = c.f[0].recip();
    let de = dh * sys.eps;
    for k in 0..kk {
        y[k + 1] = de * w[k] * sys.rk[k + 1];
        yd[k + 1] = de * wd[k] * sys.rk[k + 1];
        let mut s = z;
        let mut sd = z;
        for j in 0..=k {
            s = s + c.a[j] * y[k - j] + c.b[j] * w[k - j];
            sd = sd + c.a[j] * yd[k - j] + c.ad[j] * y[k - j] + c.b[j] * wd[k - j];
        }
        s = s * dh;
        sd = sd * dh;
        for j in 1..=k {
            let fac = c.f[j] * T::from_f64((k + 1 - j) as f64);
            s = s - fac * w[k + 1 - j];
            sd = sd - fac * wd[k + 1 - j];
        }
        w[k + 1] = s * sys.rk[k + 1] * inv_f0;
        wd[k + 1] = sd * sys.rk[k + 1] * inv_f0;
    }
    [y, w, yd, wd]
}

#[derive(Default, Clone, Debug)]
pub struct FlowStats {
    pub steps: usize,
    pub retries: usize,
    pub hmin: f64,
    pub hmax: f64,
    /// (l, Y) at the end of each accepted step, from l_from downward
    pub trace: Vec<(f64, f64)>,
}

/// Inward flow from l_from down to l_to in adaptive steps (the step is the largest whose last two terms are below
/// tol times the state's size, multishoot's plan_path rule)
pub fn flow_inward<T: Num>(sys: &Sys<T>, om: T, l_from: T, l_to: T, v0: V4<T>, cfg: &LaneCfg) -> (V4<T>, FlowStats) {
    let kk = cfg.kt;
    let mut l = l_from;
    let mut v = v0;
    let mut h = 1.0f64;
    let mut st = FlowStats { hmin: f64::INFINITY, ..Default::default() };
    let hcap = 30.0;
    while mag(l - l_to) > 0.0 {
        let rem = (l - l_to).to_f64();
        let mut hs = h.min(hcap);
        loop {
            let last = hs >= 0.95 * rem;
            let hl = if last { l - l_to } else { T::from_f64(hs) };
            let hf = hl.to_f64();
            let c = sys.coef(l, -hl, om, kk);
            let ser = regular(sys, &c, -hl, v, kk);
            let scale = v.iter().map(|&x| mag(x)).fold(0.0, f64::max);
            let mut fac = f64::INFINITY;
            let mut e_last = 0.0f64;
            for m in [kk - 1, kk] {
                let e = ser.iter().map(|s| mag(s[m])).fold(0.0, f64::max);
                e_last = e_last.max(e);
                if e > 0.0 {
                    fac = fac.min((cfg.tol * scale / e).powf(1.0 / m as f64));
                }
            }
            if fac >= 1.0 {
                // accept
                let mut nv = [T::zero(); 4];
                for i in 0..4 {
                    nv[i] = ser[i].iter().fold(T::zero(), |a, &x| a + x);
                }
                v = nv;
                l = if last { l_to } else { l - hl };
                st.trace.push((l.to_f64(), v[0].to_f64()));
                st.steps += 1;
                st.hmin = st.hmin.min(hf);
                st.hmax = st.hmax.max(hf);
                h = hf * fac.min(3.0) * 0.8;
                break;
            }
            st.retries += 1;
            hs = hf * fac * 0.9;
            assert!(hs > 1e-8, "step collapsed at l = {}", l.to_f64());
        }
    }
    (v, st)
}

/// l_inf on a 0.25 grid: the first l >= l0 with the integral of sigma from l0 at least t, sigma(l) =
/// sqrt(1/4 + 4 eps (1 + eps Om) e^{2 eps l}), the separation of the two solutions' exponents
pub fn ell_inf(eps: f64, om: f64, l0: f64, t: f64, planted: bool) -> f64 {
    // the planted system has no e^{2 eps l} term at its eigenvalue: its exponents differ by 1/2 or more
    let sig = |l: f64| if planted { 0.5 } else { (0.25 + 4.0 * eps * (1.0 + eps * om) * (2.0 * eps * l).exp()).sqrt() };
    let mut l = l0;
    let mut acc = 0.0;
    while acc < t {
        acc += 0.25 * 0.5 * (sig(l) + sig(l + 0.25));
        l += 0.25;
    }
    l
}
