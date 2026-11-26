//! Planted systems with known solutions, run through the same `integrate` as the real problem.
//!
//! case=const: y' = A y with A = blockdiag([[0, 1], [g_k^2, 0]]) + C, g_k = gamma + (0, 1/2, 3/4)_k and
//! C_ij = (((7 i + 3 j + 1) mod 11) - 5) / 8, which has the solution expm(A r) y0 growing like e^(gamma r); the
//! growing and decaying planes are nearly parallel (angles about 2/gamma), so D falls like gamma^-3 as the real one. The regular plane Y is integrated outward from 0 to L,
//! the decaying plane Z inward from L + Lz to L, and D = det[Y Z] of the orthonormalized planes is compared with
//! mpmath's expm in taylor_polar.py. Starts: Y0 and Z0 below.
//! case=linear: A(r) = A0 + r A1 with A0 the constant matrix and A1_ij = (((5 i + 2 j + 3) mod 7) - 3) / 16: no closed form;
//! the reference is mpmath's odefun (also Taylor, but arbitrary precision and no shared code).
use crate::num::Num;
use crate::plane::*;
use crate::ser::KMAX;
use crate::Args;
use multidouble::dd::DD;

pub fn c_entry(i: usize, j: usize) -> f64 {
    (((7 * i + 3 * j + 1) % 11) as f64 - 5.0) / 8.0
}

pub fn a1_entry(i: usize, j: usize) -> f64 {
    (((5 * i + 2 * j + 3) % 7) as f64 - 3.0) / 16.0
}

pub const DELTA: [f64; 3] = [0.0, 0.5, 0.75];
pub const Y0: [[f64; 6]; 3] = [[1.0, 0.0, 0.0, 0.5, 0.0, 0.125], [0.0, 1.0, 0.0, 0.25, 0.375, 0.0], [0.0, 0.0, 1.0, 0.0, 0.125, 0.375]];
pub const Z0: [[f64; 6]; 3] = [[0.25, 0.0, 0.125, 1.0, 0.0, 0.0], [0.125, 0.375, 0.0, 0.0, 1.0, 0.0], [0.0, 0.5, 0.25, 0.0, 0.0, 1.0]];

pub struct Planted<T: Num> {
    pub gamma: T,
    pub linear: bool,
}

impl<T: Num> JetSource<T> for Planted<T> {
    fn pattern(&self) -> Vec<(usize, usize)> {
        let mut v = vec![];
        for i in 0..6 {
            for j in 0..6 {
                v.push((i, j));
            }
        }
        v
    }
    fn jets(&self, x0: f64, _s2: T, kk: usize, out: &mut Jets<T>) {
        for en in out.ent.iter_mut() {
            let (i, j) = (en.i, en.j);
            let mut a0 = T::from_f64(c_entry(i, j));
            // three blocks [[0, 1], [g^2, 0]] with g = gamma + delta_k: growth e^(g r), decay e^(-g r)
            if i / 2 == j / 2 {
                let g = self.gamma + T::from_f64(DELTA[i / 2]);
                if j == i + 1 {
                    a0 = a0 + T::one();
                } else if i == j + 1 {
                    a0 = a0 + g * g;
                }
            }
            en.c = [T::zero(); KMAX];
            if self.linear {
                // A0 + (x0 + t) A1
                let a1 = T::from_f64(a1_entry(i, j));
                en.c[0] = a0 + T::from_f64(x0) * a1;
                en.c[1] = if crate::mutate::on("lin_drop") { T::zero() } else { a1 };
                en.n = 2.min(kk);
            } else {
                en.c[0] = a0;
                en.n = 1;
            }
        }
    }
}

fn dump<T: Num, const C: usize>(name: &str, q: &[[T; 6]; C]) {
    for (c, col) in q.iter().enumerate() {
        let v: Vec<String> = col.iter().map(|x| format!("\"{}\"", x.sci(130))).collect();
        println!("{{\"cmd\":\"planted_q\",\"name\":\"{}\",\"col\":{},\"v\":[{}]}}", name, c, v.join(","));
    }
}

pub fn run(a: &Args) {
    match a.s("lane", "dd") {
        "f64" => go::<f64>(a),
        "dd" => go::<DD>(a),
        "qd" => go::<multidouble::qd::QD>(a),
        "od" => go::<multidouble::od::OD>(a),
        x => panic!("lane {x}"),
    }
}

fn go<T: Num>(a: &Args) {
    let gamma = a.f("gamma", 10.0);
    let l_end = a.f("L", 1.0);
    let lz = a.f("Lz", 1.0);
    let order = a.u("order", 30);
    let ctl = Ctl { order, tol: a.f("tol", T::EPS * 100.0), hmax: a.f("hmax", 1.0), kj: a.u("kj", order) };
    let src = Planted::<T> { gamma: T::from_f64(gamma), linear: a.s("case", "const") == "linear" };
    let to = |m: &[[f64; 6]; 3]| {
        let mut y = [[T::zero(); 6]; 3];
        for c in 0..3 {
            for i in 0..6 {
                y[c][i] = T::from_f64(m[c][i]);
            }
        }
        y
    };
    let (mut y, mut z) = (to(&Y0), to(&Z0));
    let t = std::time::Instant::now();
    let (lny, sty) = integrate::<T, 6, 3>(&src, &mut y, 0.0, l_end, T::zero(), &ctl);
    let (lnz, stz) = integrate::<T, 6, 3>(&src, &mut z, l_end + lz, l_end, T::zero(), &ctl);
    let wall = t.elapsed().as_secs_f64();
    let mut mat = [[T::zero(); 6]; 6];
    for c in 0..3 {
        for i in 0..6 {
            mat[i][c] = y[c][i];
            mat[i][3 + c] = z[c][i];
        }
    }
    let d = crate::polar::det_n(&mat);
    println!(
        "{{\"cmd\":\"planted\",\"case\":\"{}\",\"lane\":\"{}\",\"gamma\":{},\"L\":{},\"Lz\":{},\"d\":\"{}\",\"lnvol_y\":\"{}\",\"lnvol_z\":\"{}\",\"steps_y\":{},\"steps_z\":{},\"hmin\":{:e},\"order\":{},\"tol\":{:e},\"wall\":{:.4}}}",
        a.s("case", "const"),
        T::NAME,
        gamma,
        l_end,
        lz,
        d.sci(130),
        lny.sci(130),
        lnz.sci(130),
        sty.steps,
        stz.steps,
        sty.hmin,
        order,
        ctl.tol,
        wall
    );
    dump("y", &y);
    dump("z", &z);
}
