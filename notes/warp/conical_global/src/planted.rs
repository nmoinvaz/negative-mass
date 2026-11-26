//! Planted cases: equations with known solutions, run through the same recursion, walk and matching code.
//!
//! linear: (1 + x^2) F'' + 4 x F' + 2 (F - alpha) = 0, whose solutions with N = 0 are F = alpha + (a + b x)/(1 + x^2);
//! the ends are its truncated expansions in 1/x at |x| = L. cubic: M N_p with M = 1 + F'^2 + F''^2, the same solutions,
//! cubic in F''. blow-ups: F'' = 2 F^3 (F = 1/(x* - x), or 1/(x - x*)) and F F'' + F'^2 = 0 (F^2 linear in x).
use crate::model::{Ends, Implicit, Mutation};
use crate::poly::PolyN;
use multishoot::lane::Lane;
use multishoot::{Fl, State, DD};

pub fn dd(num: f64, den: f64) -> DD {
    DD::from_f64(num) / DD::from_f64(den)
}

pub fn linear_poly(alpha: DD) -> PolyN {
    let one = DD::ONE;
    PolyN::new(&[
        (0, 0, 1, 0, one),
        (0, 0, 1, 2, one),
        (0, 1, 0, 1, DD::from_f64(4.0)),
        (1, 0, 0, 0, DD::from_f64(2.0)),
        (0, 0, 0, 0, -(alpha.mul_f(2.0))),
    ])
}

pub fn cubic_poly(alpha: DD) -> PolyN {
    let one = DD::ONE;
    let m = PolyN::new(&[(0, 0, 0, 0, one), (0, 2, 0, 0, one), (0, 0, 2, 0, one)]);
    linear_poly(alpha).mul(&m)
}

pub fn blowup_cube() -> PolyN {
    PolyN::new(&[(0, 0, 1, 0, DD::ONE), (3, 0, 0, 0, DD::from_f64(-2.0))])
}

pub fn blowup_sqrt() -> PolyN {
    PolyN::new(&[(1, 0, 1, 0, DD::ONE), (0, 2, 0, 0, DD::ONE)])
}

/// (F, F', F'') of alpha + (a + b x)/(1 + x^2)
pub fn exact(alpha: DD, a: DD, b: DD, x: DD) -> [DD; 3] {
    let v = DD::ONE + x * x;
    let u = a + b * x;
    let f = alpha + u / v;
    let f1 = b / v - (x * u).mul_f(2.0) / (v * v);
    let f2 = -(x * b).mul_f(4.0) / (v * v) - u.mul_f(2.0) / (v * v) + (x * x * u).mul_f(8.0) / (v * v * v);
    [f, f1, f2]
}

pub struct PlantedEnds {
    pub a_star: DD,
    pub b_star: DD,
    pub l: f64,
    /// terms j = 1..=n of the expansion
    pub n: usize,
    pub m: Mutation,
}

impl PlantedEnds {
    /// (F, F', F'') of alpha + sum_(j=1..n) (-1)^(j+1) (a x^(-2j) + b x^(-2j+1)), termwise
    fn series<T: Lane>(&self, x: f64, a: Fl<T, 1>, b: Fl<T, 1>, alpha: T) -> State<T, 1> {
        let xd = DD::from_f64(x);
        let ix = DD::ONE / xd;
        let (mut f, mut f1, mut f2) = (Fl([alpha]), Fl::<T, 1>::zero(), Fl::<T, 1>::zero());
        for j in 1..=self.n {
            let sg = if j % 2 == 1 { 1.0 } else { -1.0 };
            for (c, m) in [(a, 2 * j), (b, 2 * j - 1)] {
                let mf = m as f64;
                let p = ix.powi(m as i32);
                f = f + c.muls(p.mul_f(sg));
                f1 = f1 + c.muls((p * ix).mul_f(-mf * sg));
                f2 = f2 + c.muls((p * ix * ix).mul_f(mf * (mf + 1.0) * sg));
            }
        }
        vec![f, f1, f2]
    }
}

impl Ends for PlantedEnds {
    fn left<T: Lane>(&self, x: f64, data: &[Fl<T, 1>], g: &[T]) -> State<T, 1> {
        assert!((x + self.l).abs() < 1e-9);
        let mut s = self.series(x, data[0], Fl::splat(self.b_star), g[0]);
        if self.m.match_sign {
            s[1] = -s[1];
        }
        s[2] = s[2] + data[1];
        s
    }
    fn right<T: Lane>(&self, x: f64, data: &[Fl<T, 1>], g: &[T]) -> State<T, 1> {
        assert!((x - self.l).abs() < 1e-9);
        self.series(x, Fl::splat(self.a_star), data[0], g[0])
    }
}

pub fn planted_linear(alpha: DD, a: DD, b: DD, m: Mutation) -> Implicit<PlantedEnds> {
    Implicit { n: linear_poly(alpha), ends: PlantedEnds { a_star: a, b_star: b, l: 40.0, n: 16, m }, skip_order: None }
}

pub fn planted_cubic(alpha: DD, a: DD, b: DD, m: Mutation) -> Implicit<PlantedEnds> {
    Implicit { n: cubic_poly(alpha), ends: PlantedEnds { a_star: a, b_star: b, l: 40.0, n: 16, m }, skip_order: None }
}
