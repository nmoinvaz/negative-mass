//! P1, P2: a planted radial eigenproblem with a known spectrum, run through the same plane integrator, matching
//! determinant and root finder as the star's.
//!
//! Two decoupled Poschl-Teller problems, linear in the eigenvalue s,
//!     u_i'' = [b_i - c_i s - lam_i (lam_i + 1) sech^2 r] u_i,   u_i(0) = 0,   i = 1, 2,
//! have the bound states s = (b_i - (lam_i - n)^2) / c_i for odd n < lam_i (the odd states of the line problem), all
//! in closed form. y = S(r) y~ mixes them, S = (I + g1 e12)(I + g2 e21), g1 = 2/5 / (1 + (r/8)^2), g2 = 3/10 + r/20, so
//! J = S J~ S^-1 + S' S^-1 has no zero entry of its own and the two planes mix. The regular plane starts exactly at
//! r = 0 (S(0) times the columns (u, u') = (0, 1) of each block); the decaying plane starts from a fixed plane at R, as
//! the star's does. Parameters: lam_1 = 18/5 (b_1 = 4, c_1 = 2), lam_2 = 27/10 (b_2 = 3, c_2 = 1).
use crate::linalg::det;
use crate::num::Num;
use crate::plane::*;
use crate::ser::*;

type M4<T> = [[Ser<T>; 4]; 4];

pub const LAM: [(i64, i64); 2] = [(18, 5), (27, 10)];
pub const B: [i64; 2] = [4, 3];
pub const C: [i64; 2] = [2, 1];

/// the exact eigenvalues: (s, block, n), s = (b - (lam - n)^2) / c, as exact rationals p/q
pub fn exact_rational() -> Vec<(i64, i64, usize, usize)> {
    let mut v = vec![];
    for blk in 0..2 {
        let (ln, ld) = LAM[blk];
        let mut n = 1;
        while n * ld < ln {
            // kappa = lam - n = (ln - n ld)/ld ; s = (b - kappa^2)/c = (b ld^2 - (ln - n ld)^2) / (c ld^2)
            let k = ln - n * ld;
            let num = B[blk] * ld * ld - k * k;
            let den = C[blk] * ld * ld;
            v.push((num, den, blk, n as usize));
            n += 2;
        }
    }
    v.sort_by(|a, b| (a.0 as f64 / a.1 as f64).partial_cmp(&(b.0 as f64 / b.1 as f64)).unwrap());
    v
}

pub fn exact<T: Num>() -> Vec<(T, usize, usize)> {
    exact_rational().into_iter().map(|(p, q, blk, n)| (T::ratio(p, q), blk, n)).collect()
}

fn zero_m<T: Num>() -> M4<T> {
    [[Ser::zero(); 4]; 4]
}

fn ident<T: Num>() -> M4<T> {
    let mut m = zero_m::<T>();
    for i in 0..4 {
        m[i][i] = Ser::cst(T::one());
    }
    m
}

fn matmul<T: Num>(a: &M4<T>, b: &M4<T>) -> M4<T> {
    let mut o = zero_m::<T>();
    for i in 0..4 {
        for j in 0..4 {
            let mut s = Ser::zero();
            for k in 0..4 {
                if a[i][k].n > 0 && b[k][j].n > 0 {
                    s = add(&s, &mul(&a[i][k], &b[k][j]));
                }
            }
            o[i][j] = s;
        }
    }
    o
}

fn deriv<T: Num>(a: &Ser<T>) -> Ser<T> {
    let mut o = Ser::zero();
    for k in 0..a.n.saturating_sub(1) {
        o.c[k] = a.c[k + 1] * T::from_f64((k + 1) as f64);
    }
    o.n = a.n.saturating_sub(1);
    o
}

/// sech^2 about x0 as a series in t = r - x0
fn sech2<T: Num>(x0: T) -> Ser<T> {
    let e = x0.exp();
    let ei = e.recip();
    let mut c = Ser::zero();
    let mut fact_inv = T::one();
    for k in 0..KMAX {
        if k > 0 {
            fact_inv = fact_inv * T::from_f64(k as f64).recip();
        }
        let sgn = if k % 2 == 0 { ei } else { -ei };
        c.c[k] = (e + sgn) * T::from_f64(0.5) * fact_inv;
    }
    c.n = KMAX;
    recip(&powi(&c, 2))
}

pub struct Planted<T: Num> {
    pub s: T,
    pub sign: T,
}

impl<T: Num> Planted<T> {
    /// the series of J(r; s) about x0
    pub fn j_at(&self, x0: f64, s: T) -> M4<T> {
        let x0t = T::from_f64(x0);
        let r = Ser::var(x0t);
        let sh = sech2(x0t);
        let mut jt = zero_m::<T>();
        for blk in 0..2 {
            let lam = T::ratio(LAM[blk].0, LAM[blk].1);
            let b = T::from_f64(B[blk] as f64);
            let c = T::from_f64(C[blk] as f64);
            jt[2 * blk][2 * blk + 1] = Ser::cst(T::one());
            // q = b - c s - sign lam (lam + 1) sech^2
            let well = scl(&sh, self.sign * lam * (lam + T::one()));
            jt[2 * blk + 1][2 * blk] = adds(&neg(&well), b - c * s);
        }
        let one = T::one();
        let g1 = scl(&recip(&adds(&scl(&powi(&r, 2), T::ratio(1, 64)), one)), T::ratio(2, 5));
        let g2 = adds(&scl(&r, T::ratio(1, 20)), T::ratio(3, 10));
        let (g1p, g2p) = (deriv(&g1), deriv(&g2));
        let unit = |i: usize, j: usize, g: &Ser<T>| -> M4<T> {
            let mut m = zero_m::<T>();
            m[i][j] = *g;
            m
        };
        let addm = |a: &M4<T>, b: &M4<T>| -> M4<T> {
            let mut o = zero_m::<T>();
            for i in 0..4 {
                for j in 0..4 {
                    o[i][j] = add(&a[i][j], &b[i][j]);
                }
            }
            o
        };
        let negm = |a: &M4<T>| -> M4<T> {
            let mut o = zero_m::<T>();
            for i in 0..4 {
                for j in 0..4 {
                    o[i][j] = neg(&a[i][j]);
                }
            }
            o
        };
        let id = ident::<T>();
        // S = A B, A = I + g1 e_{1,2}, B = I + g2 e_{2,1}; the indices are (row, column)
        let a_m = addm(&id, &unit(1, 2, &g1));
        let b_m = addm(&id, &unit(2, 1, &g2));
        let ai = addm(&id, &negm(&unit(1, 2, &g1)));
        let bi = addm(&id, &negm(&unit(2, 1, &g2)));
        let s_m = matmul(&a_m, &b_m);
        let si = matmul(&bi, &ai);
        // S' = A' B + A B'
        let sp = addm(&matmul(&unit(1, 2, &g1p), &b_m), &matmul(&a_m, &unit(2, 1, &g2p)));
        let j = addm(&matmul(&matmul(&s_m, &jt), &si), &matmul(&sp, &si));
        j
    }

    pub fn s_at0(&self) -> [[T; 4]; 4] {
        let sm = self.j_s_only(0.0);
        let mut o = [[T::zero(); 4]; 4];
        for i in 0..4 {
            for j in 0..4 {
                o[i][j] = sm[i][j].get(0);
            }
        }
        o
    }

    fn j_s_only(&self, x0: f64) -> M4<T> {
        let r = Ser::var(T::from_f64(x0));
        let one = T::one();
        let g1 = scl(&recip(&adds(&scl(&powi(&r, 2), T::ratio(1, 64)), one)), T::ratio(2, 5));
        let g2 = adds(&scl(&r, T::ratio(1, 20)), T::ratio(3, 10));
        let mut a_m = ident::<T>();
        a_m[1][2] = g1;
        let mut b_m = ident::<T>();
        b_m[2][1] = g2;
        matmul(&a_m, &b_m)
    }
}

impl<T: Num> JetSource<T> for Planted<T> {
    fn pattern(&self) -> Vec<(usize, usize)> {
        let mut v = vec![];
        for i in 0..4 {
            for j in 0..4 {
                v.push((i, j));
            }
        }
        v
    }
    fn jets(&self, x0: f64, s2: T, kk: usize, out: &mut Jets<T>) {
        let j = self.j_at(x0, s2);
        for en in out.ent.iter_mut() {
            let s = &j[en.i][en.j];
            for k in 0..kk {
                en.c[k] = s.get(k);
            }
            en.n = if s.n == 0 { 0 } else { kk.min(s.n) };
        }
    }
}

/// the regular and the decaying planes of the planted problem at s, orthonormalized, and D = det[Y Z] at r_m
pub fn d_planted<T: Num>(sign: T, s: T, r_m: f64, r_far: f64, ctl: &Ctl) -> T {
    let src = Planted { s, sign };
    // regular: S(0) (0, 1, 0, 0)^T and S(0) (0, 0, 0, 1)^T
    let s0 = src.s_at0();
    let mut y = [[T::zero(); 4]; 2];
    for i in 0..4 {
        y[0][i] = s0[i][1];
        y[1][i] = s0[i][3];
    }
    // decaying: a fixed generic plane at R
    let mut z = [[T::zero(); 4]; 2];
    let z0 = [[1.0, 0.3, 0.0, 0.2], [0.0, 0.5, 1.0, 0.1]];
    for c in 0..2 {
        for i in 0..4 {
            z[c][i] = T::from_f64(z0[c][i]);
        }
    }
    integrate::<T, 4, 2>(&src, &mut y, 0.0, r_m, s, ctl);
    integrate::<T, 4, 2>(&src, &mut z, r_far, r_m, s, ctl);
    let mut m = [[T::zero(); 4]; 4];
    for c in 0..2 {
        for i in 0..4 {
            m[i][c] = y[c][i];
            m[i][2 + c] = z[c][i];
        }
    }
    det(&m)
}
