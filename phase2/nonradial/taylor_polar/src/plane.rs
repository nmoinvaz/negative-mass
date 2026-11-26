//! The Taylor integrator for y' = J(r) y on a k-plane, with orthonormalization after each step.
//!
//! One step from x0: the jets of J's entries about x0 give the coefficients
//!   y_(k+1) = (1/(k+1)) sum_(m<=k) J_m y_(k-m),
//! the step h is the largest whose last two terms stay below `tol` (multishoot's planner rule,
//! `taylor::plan_path`: e_(kt-1) h^(kt-1) <= tol and e_kt h^kt <= tol), and the state is the Horner sum.
//! Columns are orthonormalized by classical Gram-Schmidt with one reorthogonalization after every step, so the
//! diagonal of R is positive and the determinant needs no phase.
use crate::num::Num;
use crate::ser::KMAX;

/// A nonzero entry's jets about the current step's start.
#[derive(Clone)]
pub struct Entry<T: Num> {
    pub i: usize,
    pub j: usize,
    pub n: usize,
    pub c: [T; KMAX],
}

/// The nonzero pattern and the jets of an R x R matrix.
pub struct Jets<T: Num> {
    pub ent: Vec<Entry<T>>,
}

impl<T: Num> Jets<T> {
    pub fn new(pattern: &[(usize, usize)]) -> Self {
        Jets {
            ent: pattern.iter().map(|&(i, j)| Entry { i, j, n: 0, c: [T::zero(); KMAX] }).collect(),
        }
    }
}

/// Where the jets come from: the star's table, the vacuum formulas, or a planted matrix.
pub trait JetSource<T: Num> {
    /// the nonzero pattern, fixed for the source
    fn pattern(&self) -> Vec<(usize, usize)>;
    /// coefficients 0..kk of every entry about x0, at s^2 = s2
    fn jets(&self, x0: f64, s2: T, kk: usize, out: &mut Jets<T>);
}

#[derive(Clone, Copy, Debug)]
pub struct Ctl {
    /// Taylor order of the step (at most KMAX - 2)
    pub order: usize,
    /// bound on the last two terms, absolute (the columns have unit norm)
    pub tol: f64,
    /// largest step
    pub hmax: f64,
    /// jets used, at most `order` (the recursion needs J_m for m < order)
    pub kj: usize,
}

#[derive(Clone, Copy, Debug, Default)]
pub struct Stats {
    pub steps: usize,
    pub hmin: f64,
    pub hmax: f64,
}

/// 1/k for k = 1..=KMAX at the lane's precision (a quotient formed in the lane, not widened from double-double)
pub fn recips<T: Num>() -> [T; KMAX + 1] {
    let mut v = [T::zero(); KMAX + 1];
    for (k, x) in v.iter_mut().enumerate().skip(1) {
        *x = if crate::mutate::on("recip_dd") { T::from_dd(multidouble::dd::DD::ONE / multidouble::dd::DD::from_f64(k as f64)) } else { T::one() / T::from_f64(k as f64) };
    }
    v
}

/// Gram-Schmidt of the columns in place; returns the positive diagonal of R.
#[cfg(not(feature = "lib_plane"))]
pub fn gram_schmidt<T: Num, const R: usize, const C: usize>(y: &mut [[T; R]; C]) -> [T; C] {
    let mut diag = [T::zero(); C];
    for c in 0..C {
        for _pass in 0..2 {
            for p in 0..c {
                let mut d = T::zero();
                for i in 0..R {
                    d = d + y[p][i] * y[c][i];
                }
                for i in 0..R {
                    y[c][i] = y[c][i] - d * y[p][i];
                }
            }
        }
        let mut s = T::zero();
        for i in 0..R {
            s = s + y[c][i] * y[c][i];
        }
        let nrm = s.sqrt();
        let inv = nrm.recip();
        for i in 0..R {
            y[c][i] = y[c][i] * inv;
        }
        diag[c] = nrm;
    }
    diag
}

/// Integrates the plane y (C columns of R rows) from ra to rb, either direction, at s^2 = s2.
/// Returns the sum of ln R_ii over all steps (the log of the plane's volume growth) and the step statistics.
#[cfg(not(feature = "lib_plane"))]
pub fn integrate<T: Num, const R: usize, const C: usize>(
    src: &dyn JetSource<T>,
    y: &mut [[T; R]; C],
    ra: f64,
    rb: f64,
    s2: T,
    ctl: &Ctl,
) -> (T, Stats) {
    let n = ctl.order;
    assert!(n + 2 <= KMAX && ctl.kj <= n);
    let inv = recips::<T>();
    let mut jets = Jets::<T>::new(&src.pattern());
    let mut by_row: Vec<Vec<usize>> = vec![vec![]; R];
    for (e, en) in jets.ent.iter().enumerate() {
        by_row[en.i].push(e);
    }
    let dir = if rb >= ra { 1.0 } else { -1.0 };
    let max_s: f64 = std::env::var("TP_MAX_SECONDS").ok().and_then(|v| v.parse().ok()).unwrap_or(f64::INFINITY);
    let t_start = std::time::Instant::now();
    let mut x0 = ra;
    let mut lnvol = T::zero();
    let mut st = Stats { steps: 0, hmin: f64::INFINITY, hmax: 0.0 };
    // coefficient arrays: co[c][i][k]
    let mut co = vec![[[T::zero(); KMAX]; R]; C];
    let d0 = gram_schmidt(y);
    for d in d0 {
        lnvol = lnvol + d.ln();
    }
    while (rb - x0) * dir > 0.0 {
        if st.steps % 1024 == 1023 && t_start.elapsed().as_secs_f64() > max_s {
            panic!("integrate: more than {max_s} s by x = {x0:e} after {} steps", st.steps);
        }
        src.jets(x0, s2, ctl.kj, &mut jets);
        for c in 0..C {
            for i in 0..R {
                co[c][i][0] = y[c][i];
            }
        }
        for k in 0..n {
            for c in 0..C {
                for i in 0..R {
                    let mut s = T::zero();
                    for &e in &by_row[i] {
                        let en = &jets.ent[e];
                        let mm = k.min(en.n.saturating_sub(1));
                        if en.n == 0 {
                            continue;
                        }
                        for m in 0..=mm {
                            s = s + en.c[m] * co[c][en.j][k - m];
                        }
                    }
                    co[c][i][k + 1] = s * inv[k + 1];
                }
            }
        }
        // largest step with the last two terms below tol
        let mut h = ctl.hmax;
        for k in [n - 1, n] {
            let mut e = 0.0f64;
            for c in 0..C {
                for i in 0..R {
                    e = e.max(co[c][i][k].to_f64().abs());
                }
            }
            if !e.is_finite() {
                panic!("non-finite Taylor coefficient at x = {x0}, order {k}");
            }
            if e > 0.0 {
                h = h.min((ctl.tol / e).powf(1.0 / k as f64));
            }
        }
        let rem = (rb - x0).abs();
        let hs = if h >= rem { rem } else { h };
        let x1 = if hs >= rem { rb } else { x0 + dir * hs };
        let hh = T::from_f64(x1) - T::from_f64(x0);
        for c in 0..C {
            for i in 0..R {
                // mutation lowsum: the Taylor sum cut at order 8 while the step assumes order n
                let top = if crate::mutate::on("lowsum") { 8 } else { n };
                let mut s = co[c][i][top];
                for k in (0..top).rev() {
                    s = s * hh + co[c][i][k];
                }
                y[c][i] = s;
            }
        }
        if !crate::mutate::on("noqr") {
            let d = gram_schmidt(y);
            for dd in d {
                lnvol = lnvol + dd.ln();
            }
        }
        st.steps += 1;
        st.hmin = st.hmin.min((x1 - x0).abs());
        st.hmax = st.hmax.max((x1 - x0).abs());
        x0 = x1;
    }
    (lnvol, st)
}

/// Gram-Schmidt of the columns in place by multishoot's `plane::orthonormalize`; returns the positive diagonal of R.
#[cfg(feature = "lib_plane")]
pub fn gram_schmidt<T: Num, const R: usize, const C: usize>(y: &mut [[T; R]; C]) -> [T; C] {
    let mut flat: Vec<T::Lib> = y.iter().flatten().map(|x| x.to_lib()).collect();
    let mut diag = vec![<T::Lib as multishoot::plane::Real>::zero(); C];
    multishoot::plane::orthonormalize(&mut flat, R, C, &mut diag);
    for (c, col) in y.iter_mut().enumerate() {
        for (i, x) in col.iter_mut().enumerate() {
            *x = T::from_lib(flat[c * R + i]);
        }
    }
    std::array::from_fn(|c| T::from_lib(diag[c]))
}

/// The crate's `JetSource` seen as multishoot's: the jets at s^2 = s2 into the crate's own `Jets`, copied across. It also
/// keeps `integrate`'s wall-time guard (`TP_MAX_SECONDS`), checked every 1024 steps.
#[cfg(feature = "lib_plane")]
struct Bridge<'a, T: Num> {
    src: &'a dyn JetSource<T>,
    s2: T,
    scratch: std::cell::RefCell<Jets<T>>,
    calls: std::cell::Cell<usize>,
    t_start: std::time::Instant,
    max_s: f64,
}

#[cfg(feature = "lib_plane")]
impl<T: Num> multishoot::plane::JetSource<T::Lib> for Bridge<'_, T> {
    fn pattern(&self) -> Vec<(usize, usize)> {
        self.src.pattern()
    }
    fn jets(&self, x0: f64, kj: usize, out: &mut multishoot::plane::Jets<T::Lib>) {
        let calls = self.calls.get();
        self.calls.set(calls + 1);
        if calls % 1024 == 1023 && self.t_start.elapsed().as_secs_f64() > self.max_s {
            panic!("integrate: more than {} s by x = {x0:e} after {calls} steps", self.max_s);
        }
        let mut sc = self.scratch.borrow_mut();
        self.src.jets(x0, self.s2, kj, &mut sc);
        for (e, en) in out.ent.iter_mut().enumerate() {
            let n = sc.ent[e].n.min(kj);
            for k in 0..n {
                en.c[k] = sc.ent[e].c[k].to_lib();
            }
            en.n = n;
        }
    }
}

/// Integrates the plane y (C columns of R rows) from ra to rb, either direction, at s^2 = s2, by multishoot's
/// `plane::Plane::advance` (the mutations `lowsum` and `noqr` of the crate's own integrator do not apply here).
/// Returns the sum of ln R_ii over all steps (the log of the plane's volume growth) and the step statistics.
#[cfg(feature = "lib_plane")]
pub fn integrate<T: Num, const R: usize, const C: usize>(
    src: &dyn JetSource<T>,
    y: &mut [[T; R]; C],
    ra: f64,
    rb: f64,
    s2: T,
    ctl: &Ctl,
) -> (T, Stats) {
    assert!(ctl.order + 2 <= KMAX && ctl.kj <= ctl.order);
    for m in ["lowsum", "noqr"] {
        if crate::mutate::on(m) {
            eprintln!("lib_plane: the mutation {m} is not applied by multishoot's stepper");
        }
    }
    let max_s: f64 = std::env::var("TP_MAX_SECONDS").ok().and_then(|v| v.parse().ok()).unwrap_or(f64::INFINITY);
    let bridge = Bridge {
        src,
        s2,
        scratch: std::cell::RefCell::new(Jets::<T>::new(&src.pattern())),
        calls: std::cell::Cell::new(0),
        t_start: std::time::Instant::now(),
        max_s,
    };
    let basis: Vec<T::Lib> = y.iter().flatten().map(|x| x.to_lib()).collect();
    let mut pl = multishoot::plane::Plane::<T::Lib>::new(R, C, &basis);
    let plan = multishoot::Plan { kt: ctl.order, tol: ctl.tol, hmax: ctl.hmax };
    let st = pl.advance(&bridge, ra, rb, &plan, ctl.kj);
    for (c, col) in y.iter_mut().enumerate() {
        for (i, x) in col.iter_mut().enumerate() {
            *x = T::from_lib(pl.column(c)[i]);
        }
    }
    (T::from_lib(pl.ln_vol()), Stats { steps: st.steps, hmin: st.hmin, hmax: st.hmax })
}
