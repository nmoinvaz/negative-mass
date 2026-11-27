//! Planted complex systems with known solutions, run through the same `integrate` as the real problem (negative-mass #188).
//!
//! `cp1`: P1 of multishoot's plane tests with a rotating phase. J = P diag(24 + i f1, 8 + i f2, -8 + i f3, -24 + i f4) P^-1, P the
//! unimodular integer matrix of that test (cond 96), so J has Gaussian-integer entries. The plane P[:, :2] is invariant:
//! Y(x) = P[:, :2] diag(e^(lambda_1 x), e^(lambda_2 x)), and the positive-diagonal orthonormal basis of Y(x) is Q0 diag(e^(i f_1 x),
//! e^(i f_2 x)), Q0 the one of P[:, :2]. `cx=0` runs the same system on the real lane (f = 0) for the bit-for-bit comparison.
//!
//! `czserve`: a planted determinant with zeros. J = A + i B + i diag(c) on 6 x 6 (below), the regular plane Y from 0 to L, the
//! decaying plane Z inward from L + Lz, both from constant start data except Y's third column, Y0(s)[:, 2] = zt + p(s) u with
//! p(s) = (s - s0)(s - conj s0), u a real column, and zt = exp(-J (L + Lz)) Z0[:, 0] (given by the caller, who computes it in mpmath):
//! det[Y(L) | Z(L)] = c p(s) for a constant c, so D has exactly the zeros of p. The output has det[Q_Y | Q_Z] and the two log-volumes.
use crate::complex::Cx;
use crate::num::Num;
use crate::planted::{Y0, Z0};
use crate::plane::*;
use crate::polar::det_n;
use crate::ser::KMAX;
use crate::Args;
use multidouble::dd::DD;
use std::io::{BufRead, Write};

pub const P: [[f64; 4]; 4] = [[1., 1., 0., 1.], [1., 2., 1., 1.], [0., 1., 2., 1.], [1., 1., 1., 3.]];
pub const PINV: [[f64; 4]; 4] = [[7., -4., 3., -2.], [-4., 3., -2., 1.], [3., -2., 2., -1.], [-2., 1., -1., 1.]];
pub const G: [f64; 4] = [24., 8., -8., -24.];
pub const F: [f64; 4] = [7., -3., 5., -11.];

/// a lane value from (re, im); the real lanes take only im = 0
pub trait Embed: Num {
    fn from_c(re: f64, im: f64) -> Self;
    /// decimal strings of the real and imaginary parts
    fn parts(self, digits: usize) -> (String, String);
}

macro_rules! real_embed {
    ($t:ty) => {
        impl Embed for $t {
            fn from_c(re: f64, im: f64) -> Self {
                assert!(im == 0.0, "a real lane needs f = 0");
                <$t as Num>::from_f64(re)
            }
            fn parts(self, digits: usize) -> (String, String) {
                (self.sci(digits), "0".to_string())
            }
        }
    };
}
real_embed!(f64);
real_embed!(DD);
real_embed!(multidouble::qd::QD);

impl<T: Num> Embed for Cx<T> {
    fn from_c(re: f64, im: f64) -> Self {
        Cx::new(T::from_f64(re), T::from_f64(im))
    }
    fn parts(self, digits: usize) -> (String, String) {
        (self.re.sci(digits), self.im.sci(digits))
    }
}

/// constant J of P1c: J_ij = sum_k P_ik (g_k + i f_k) PINV_kj, as (re, im) in f64 (integers)
pub fn j_p1c(f: &[f64; 4]) -> [[(f64, f64); 4]; 4] {
    let mut j = [[(0.0, 0.0); 4]; 4];
    for (i, row) in j.iter_mut().enumerate() {
        for (jj, e) in row.iter_mut().enumerate() {
            let (mut re, mut im) = (0.0, 0.0);
            for k in 0..4 {
                re += P[i][k] * G[k] * PINV[k][jj];
                im += P[i][k] * f[k] * PINV[k][jj];
            }
            *e = (re, im);
        }
    }
    j
}

struct Const<T: Num, const N: usize> {
    j: [[(f64, f64); N]; N],
    mk: fn(f64, f64) -> T,
}

impl<T: Num, const N: usize> JetSource<T> for Const<T, N> {
    fn pattern(&self) -> Vec<(usize, usize)> {
        let mut v = vec![];
        for i in 0..N {
            for j in 0..N {
                v.push((i, j));
            }
        }
        v
    }
    fn jets(&self, _x0: f64, _s2: T, _kk: usize, out: &mut Jets<T>) {
        for en in out.ent.iter_mut() {
            en.c = [T::zero(); KMAX];
            en.c[0] = (self.mk)(self.j[en.i][en.j].0, self.j[en.i][en.j].1);
            en.n = 1;
        }
    }
}

fn cp1_go<T: Embed>(a: &Args) {
    let x = a.f("x", 12.0);
    let fl = a.list("f");
    let f: [f64; 4] = if fl.len() == 4 { [fl[0], fl[1], fl[2], fl[3]] } else { F };
    let order = a.u("order", 30);
    let ctl = Ctl { order, tol: a.f("tol", T::EPS * 100.0), hmax: a.f("hmax", 1.0), kj: a.u("kj", order) };
    let src = Const::<T, 4> { j: j_p1c(&f), mk: T::from_c };
    let mut y = [[T::zero(); 4]; 2];
    for (c, col) in y.iter_mut().enumerate() {
        for (i, v) in col.iter_mut().enumerate() {
            *v = T::from_f64(P[i][c]);
        }
    }
    let t = std::time::Instant::now();
    let (lnvol, st) = integrate::<T, 4, 2>(&src, &mut y, 0.0, x, T::zero(), &ctl);
    let wall = t.elapsed().as_secs_f64();
    let (lre, lim) = lnvol.parts(130);
    println!(
        "{{\"cmd\":\"cp1\",\"lane\":\"{}\",\"x\":{},\"f\":{:?},\"lnvol_re\":\"{}\",\"lnvol_im\":\"{}\",\"steps\":{},\"hmin\":{:e},\"order\":{},\"tol\":{:e},\"wall\":{:.4}}}",
        T::NAME,
        x,
        f,
        lre,
        lim,
        st.steps,
        st.hmin,
        order,
        ctl.tol,
        wall
    );
    for (c, col) in y.iter().enumerate() {
        let (re, im): (Vec<String>, Vec<String>) = col.iter().map(|v| v.parts(130)).map(|(r, i)| (format!("\"{r}\""), format!("\"{i}\""))).unzip();
        println!("{{\"cmd\":\"cp1_q\",\"col\":{},\"re\":[{}],\"im\":[{}]}}", c, re.join(","), im.join(","));
    }
}

/// `cp1 lane= cx=1|0 x= f=a,b,c,d`: the planted complex plane (cx = 0: the real lane, f must be 0)
pub fn cp1(a: &Args) {
    let cx = a.u("cx", 1) == 1;
    match (a.s("lane", "dd"), cx) {
        ("f64", true) => cp1_go::<Cx<f64>>(a),
        ("dd", true) => cp1_go::<Cx<DD>>(a),
        ("qd", true) => cp1_go::<Cx<multidouble::qd::QD>>(a),
        ("f64", false) => cp1_go::<f64>(a),
        ("dd", false) => cp1_go::<DD>(a),
        ("qd", false) => cp1_go::<multidouble::qd::QD>(a),
        (x, _) => panic!("lane {x}"),
    }
}

// ---------------------------------------------------------------------------------------------------------------
// the planted determinant
// ---------------------------------------------------------------------------------------------------------------

/// J of the planted determinant: (re, im) of entry (i, j)
pub fn j_cz() -> [[(f64, f64); 6]; 6] {
    const C: [f64; 6] = [2.0, -3.0, 1.0, 4.0, -2.0, 3.0];
    let mut j = [[(0.0, 0.0); 6]; 6];
    for (i, row) in j.iter_mut().enumerate() {
        for (jj, e) in row.iter_mut().enumerate() {
            let re = (((7 * i + 3 * jj + 1) % 11) as f64 - 5.0) / 8.0;
            let mut im = (((3 * i + 5 * jj + 2) % 9) as f64 - 4.0) / 4.0;
            if i == jj {
                im += C[i];
            }
            *e = (re, im);
        }
    }
    j
}

fn cz_go<T: Num>(a: &Args) {
    let l = a.f("L", 2.0);
    let lz = a.f("Lz", 2.0);
    let s0 = (a.f("s0re", 0.75), a.f("s0im", 0.5));
    let order = a.u("order", 30);
    let ctl = Ctl { order, tol: a.f("tol", T::EPS * 100.0), hmax: a.f("hmax", 1.0), kj: a.u("kj", order) };
    let parse_list = |k: &str| -> Vec<Cx<T>> { a.s(k, "").split(',').filter(|x| !x.is_empty()).map(Cx::<T>::parse).collect() };
    let (zt_re, zt_im) = (parse_list("zt_re"), parse_list("zt_im"));
    assert!(zt_re.len() == 6 && zt_im.len() == 6, "zt_re and zt_im need six numbers each");
    let zt: Vec<Cx<T>> = zt_re.iter().zip(zt_im.iter()).map(|(r, i)| Cx::new(r.re, i.re)).collect();
    let jm = j_cz();
    let src = Const::<Cx<T>, 6> { j: jm, mk: |re, im| Cx::new(T::from_f64(re), T::from_f64(im)) };
    let s0c = Cx::new(T::from_f64(s0.0), T::from_f64(s0.1));
    let s0b = s0c.conj();
    println!(
        "{{\"cmd\":\"cz_setup\",\"lane\":\"{}\",\"L\":{},\"Lz\":{},\"s0\":[{:?},{:?}],\"order\":{},\"tol\":{:e},\"j_re\":{:?},\"j_im\":{:?}}}",
        T::NAME,
        l,
        lz,
        s0.0,
        s0.1,
        order,
        ctl.tol,
        jm.iter().flatten().map(|e| e.0).collect::<Vec<_>>(),
        jm.iter().flatten().map(|e| e.1).collect::<Vec<_>>()
    );
    std::io::stdout().flush().unwrap();
    let conj_s = crate::mutate::on("conj_s");
    for line in std::io::stdin().lock().lines() {
        let line = line.unwrap();
        let mut it = line.split_whitespace();
        let (Some(x), Some(y)) = (it.next(), it.next()) else { continue };
        let (x, y): (f64, f64) = (x.parse().unwrap(), y.parse().unwrap());
        // mutation conj_s: D evaluated at conj(s)
        let s = Cx::new(T::from_f64(x), T::from_f64(if conj_s { -y } else { y }));
        let p = (s - s0c) * (s - s0b);
        let mut yv = [[Cx::<T>::zero(); 6]; 3];
        for c in 0..2 {
            for i in 0..6 {
                yv[c][i] = Cx::from_f64(Y0[c][i]);
            }
        }
        for i in 0..6 {
            yv[2][i] = zt[i] + p * Cx::from_f64(Y0[2][i]);
        }
        let mut zv = [[Cx::<T>::zero(); 6]; 3];
        for c in 0..3 {
            for i in 0..6 {
                zv[c][i] = Cx::from_f64(Z0[c][i]);
            }
        }
        let (lny, _) = integrate::<Cx<T>, 6, 3>(&src, &mut yv, 0.0, l, Cx::zero(), &ctl);
        let (lnz, _) = integrate::<Cx<T>, 6, 3>(&src, &mut zv, l + lz, l, Cx::zero(), &ctl);
        let mut mat = [[Cx::<T>::zero(); 6]; 6];
        for c in 0..3 {
            for i in 0..6 {
                mat[i][c] = yv[c][i];
                mat[i][3 + c] = zv[c][i];
            }
        }
        let d = det_n(&mat);
        println!(
            "{{\"cmd\":\"cz\",\"re\":{:?},\"im\":{:?},\"d_re\":\"{}\",\"d_im\":\"{}\",\"d_f64\":[{:e},{:e}],\"ln_y\":\"{}\",\"ln_z\":\"{}\"}}",
            x,
            y,
            d.re.sci(130),
            d.im.sci(130),
            d.re.to_f64(),
            d.im.to_f64(),
            lny.re.sci(130),
            lnz.re.sci(130)
        );
        std::io::stdout().flush().unwrap();
    }
}

/// `czserve lane= zt_re=... zt_im=... s0re= s0im= L= Lz=`: reads `re im` lines, answers with det[Q_Y | Q_Z] and the log-volumes
pub fn czserve(a: &Args) {
    match a.s("lane", "dd") {
        "f64" => cz_go::<f64>(a),
        "dd" => cz_go::<DD>(a),
        "qd" => cz_go::<multidouble::qd::QD>(a),
        x => panic!("lane {x}"),
    }
}
