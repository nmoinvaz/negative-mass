//! Roots of a scalar function in the lane: the Illinois method with bisection safeguarding down to a relative
//! bracket of `STAGE1`, then secant steps to the lane's roundoff.

use crate::num::Num;

pub struct Root<T: Num> {
    pub x: T,
    pub fx: T,
    pub evals: usize,
    /// width of the bracket Illinois stopped at, and the last secant step
    pub bracket: f64,
    pub last_step: f64,
}

const STAGE1: f64 = 1e-9;

/// A root in [a, b] of f, f(a) f(b) < 0. The secant stage stops when a step is within 4 eps |x| or when a step grows
/// instead of shrinking (the noise floor of f), and returns the iterate of the smallest |f|.
pub fn find<T: Num>(f: &mut dyn FnMut(T) -> T, mut a: T, mut b: T, mut fa: T, mut fb: T) -> Root<T> {
    let mut evals = 0;
    let mut side = 0i32;
    let (mut x, mut fx) = (a, fa);
    let mut last_w = (b - a).abs().to_f64();
    let mut it = 0;
    loop {
        let w = (b - a).abs().to_f64();
        let scale = a.to_f64().abs().max(b.to_f64().abs()).max(1e-300);
        if w <= STAGE1 * scale || it > 80 {
            break;
        }
        let use_bisect = it % 3 == 2 && w > 0.5 * last_w;
        if it % 3 == 0 {
            last_w = w;
        }
        x = if use_bisect { (a + b) * T::from_f64(0.5) } else { (a * fb - b * fa) / (fb - fa) };
        let (lo, hi) = (a.to_f64().min(b.to_f64()), a.to_f64().max(b.to_f64()));
        if !(x.to_f64() > lo && x.to_f64() < hi) {
            x = (a + b) * T::from_f64(0.5);
        }
        fx = f(x);
        evals += 1;
        it += 1;
        if fx.to_f64() == 0.0 {
            return Root { x, fx, evals, bracket: w, last_step: 0.0 };
        }
        if (fx.to_f64() < 0.0) == (fb.to_f64() < 0.0) {
            b = x;
            fb = fx;
            if side == -1 {
                fa = fa * T::from_f64(0.5);
            }
            side = -1;
        } else {
            a = x;
            fa = fx;
            if side == 1 {
                fb = fb * T::from_f64(0.5);
            }
            side = 1;
        }
    }
    let bracket = (b - a).abs().to_f64();
    // stage 2: secant from the bracket's ends (their true values, not the halved ones: recompute the two ends' f)
    let (mut x0, mut f0) = (a, f(a));
    let (mut x1, mut f1) = (b, f(b));
    evals += 2;
    if f0.abs().to_f64() < f1.abs().to_f64() {
        std::mem::swap(&mut x0, &mut x1);
        std::mem::swap(&mut f0, &mut f1);
    }
    // now (x1, f1) is the better; iterate
    let (mut best_x, mut best_f) = (x1, f1);
    let mut prev_step = f64::INFINITY;
    let mut last_step = 0.0;
    for _ in 0..14 {
        if f1.to_f64() == f0.to_f64() {
            break;
        }
        let xn = x1 - f1 * (x1 - x0) / (f1 - f0);
        let step = (xn - x1).abs().to_f64();
        let scale = xn.to_f64().abs().max(1e-300);
        x0 = x1;
        f0 = f1;
        x1 = xn;
        f1 = f(x1);
        evals += 1;
        last_step = step;
        if f1.abs().to_f64() <= best_f.abs().to_f64() {
            best_x = x1;
            best_f = f1;
        }
        if step <= 4.0 * T::EPS * scale || f1.to_f64() == 0.0 || step > prev_step {
            break;
        }
        prev_step = step;
    }
    Root { x: best_x, fx: best_f, evals, bracket, last_step }
}
