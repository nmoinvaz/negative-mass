//! A polynomial N(F0, F1, F2; x) = sum coef F0^a F1^b F2^c x^m with c <= 3, and the Taylor recursion of the implicit
//! equation N = 0 for the state (F, F', F'').
//!
//! About a step's start x0 the fields are series in h = x - x0. The recursion takes order k >= 1 of N as zero and solves
//! it for the coefficient F''_k: F''_k enters the order k of N linearly, with the factor A0 = dN/dF'' at the order-0 values,
//! and F_k, F'_k are known from F'_(k-1), F''_(k-1). The order-0 term of N is not imposed, so N is constant along the
//! flow (N' = 0) and the data's N = 0 carries through. One division, by A0, per step; no hand-differentiated F'''.
use multishoot::lane::Lane;
use multishoot::{cs, Fl, Work, DD};
use std::collections::BTreeMap;

/// the monomials F0^a F1^b F2^c with a common polynomial in x: sum over m of coef x^m
#[derive(Clone, Debug)]
pub struct Group {
    pub a: usize,
    pub b: usize,
    pub c: usize,
    pub xp: Vec<(usize, DD)>,
}

#[derive(Clone, Debug)]
pub struct PolyN {
    pub groups: Vec<Group>,
    pub amax: usize,
    pub bmax: usize,
    pub cmax: usize,
}

pub fn abs_dd(x: DD) -> DD {
    if x.hi < 0.0 || (x.hi == 0.0 && x.lo < 0.0) {
        -x
    } else {
        x
    }
}

impl PolyN {
    pub fn new(terms: &[(u8, u8, u8, u8, DD)]) -> Self {
        let mut map: BTreeMap<(usize, usize, usize), Vec<(usize, DD)>> = BTreeMap::new();
        for &(a, b, c, m, co) in terms {
            map.entry((a as usize, b as usize, c as usize)).or_default().push((m as usize, co));
        }
        let groups: Vec<Group> = map.into_iter().map(|((a, b, c), xp)| Group { a, b, c, xp }).collect();
        let amax = groups.iter().map(|g| g.a).max().unwrap_or(0);
        let bmax = groups.iter().map(|g| g.b).max().unwrap_or(0);
        let cmax = groups.iter().map(|g| g.c).max().unwrap_or(0);
        assert!(cmax <= 3, "PolyN: degree {cmax} in F'' (at most 3)");
        PolyN { groups, amax, bmax, cmax }
    }

    pub fn from_ints(terms: &[(u8, u8, u8, u8, i64)]) -> Self {
        let t: Vec<_> = terms.iter().map(|&(a, b, c, m, co)| (a, b, c, m, DD::from_f64(co as f64))).collect();
        Self::new(&t)
    }

    /// the product of two polynomials
    pub fn mul(&self, o: &PolyN) -> PolyN {
        let tl = |p: &PolyN| -> Vec<(u8, u8, u8, u8, DD)> {
            p.groups
                .iter()
                .flat_map(|g| g.xp.iter().map(move |&(m, co)| (g.a as u8, g.b as u8, g.c as u8, m as u8, co)))
                .collect()
        };
        let (x, y) = (tl(self), tl(o));
        let mut out = Vec::new();
        for &(a, b, c, m, co) in &x {
            for &(a2, b2, c2, m2, co2) in &y {
                out.push((a + a2, b + b2, c + c2, m + m2, co * co2));
            }
        }
        PolyN::new(&out)
    }

    pub fn add(&self, o: &PolyN) -> PolyN {
        let mut t = Vec::new();
        for p in [self, o] {
            for g in &p.groups {
                for &(m, co) in &g.xp {
                    t.push((g.a as u8, g.b as u8, g.c as u8, m as u8, co));
                }
            }
        }
        PolyN::new(&t)
    }

    /// N and the sum of the absolute values of its terms at (f0, f1, f2, x), in double-double
    pub fn eval(&self, f0: DD, f1: DD, f2: DD, x: DD) -> (DD, DD) {
        let (mut v, mut s) = (DD::ZERO, DD::ZERO);
        for g in &self.groups {
            let base = f0.powi(g.a as i32) * f1.powi(g.b as i32) * f2.powi(g.c as i32);
            for &(m, co) in &g.xp {
                let t = co * base * x.powi(m as i32);
                v = v + t;
                s = s + abs_dd(t);
            }
        }
        (v, s)
    }

    /// dN/dF'' in double-double
    pub fn d_f2(&self, f0: DD, f1: DD, f2: DD, x: DD) -> DD {
        let mut v = DD::ZERO;
        for g in self.groups.iter().filter(|g| g.c >= 1) {
            let base = f0.powi(g.a as i32) * f1.powi(g.b as i32) * f2.powi(g.c as i32 - 1);
            for &(m, co) in &g.xp {
                v = v + co * base * x.powi(m as i32) * DD::from_f64(g.c as f64);
            }
        }
        v
    }

    /// the coefficient of F''^3, F''^2, F'', 1 at (f0, f1, x)
    pub fn coef_f2(&self, f0: DD, f1: DD, x: DD) -> [DD; 4] {
        let mut v = [DD::ZERO; 4];
        for g in &self.groups {
            let base = f0.powi(g.a as i32) * f1.powi(g.b as i32);
            for &(m, co) in &g.xp {
                v[g.c] = v[g.c] + co * base * x.powi(m as i32);
            }
        }
        [v[3], v[2], v[1], v[0]]
    }
}

/// the coefficients of the polynomial sum coef x^m in h = x - x0: e_j = sum_m coef C(m, j) x0^(m - j)
fn shifted(xp: &[(usize, DD)], x0: DD) -> Vec<DD> {
    let deg = xp.iter().map(|&(m, _)| m).max().unwrap_or(0);
    let mut e = vec![DD::ZERO; deg + 1];
    for &(m, co) in xp {
        let mut binom = 1.0f64;
        for (j, ej) in e.iter_mut().enumerate().take(m + 1) {
            *ej = *ej + co * DD::from_f64(binom) * x0.powi((m - j) as i32);
            binom = binom * (m - j) as f64 / (j + 1) as f64;
        }
    }
    e
}

/// the series of one step: the fields, their powers and products, and the pieces of N, each to order k as it is formed
struct Cx<'a, T: Lane> {
    p: &'a PolyN,
    f: [Vec<Fl<T, 1>>; 3],
    pw0: Vec<Vec<Fl<T, 1>>>,
    pw1: Vec<Vec<Fl<T, 1>>>,
    pw2: Vec<Vec<Fl<T, 1>>>,
    unit: Vec<Fl<T, 1>>,
    pairs: BTreeMap<(usize, usize), Vec<Fl<T, 1>>>,
    q: Vec<Vec<Fl<T, 1>>>,
    ex: Vec<Vec<DD>>,
}

impl<'a, T: Lane> Cx<'a, T> {
    fn new(p: &'a PolyN, x0: f64, kt: usize, s0: [Fl<T, 1>; 3]) -> Self {
        let len = kt + 2;
        let z = || vec![Fl::<T, 1>::zero(); len];
        let mut f = [z(), z(), z()];
        for i in 0..3 {
            f[i][0] = s0[i];
        }
        let mut unit = z();
        unit[0] = Fl::splat(DD::ONE);
        let mut pairs: BTreeMap<(usize, usize), Vec<Fl<T, 1>>> = BTreeMap::new();
        for g in &p.groups {
            if g.a >= 1 && g.b >= 1 {
                pairs.entry((g.a, g.b)).or_insert_with(z);
            }
        }
        let x0d = DD::from_f64(x0);
        Cx {
            p,
            f,
            pw0: (0..=p.amax).map(|_| z()).collect(),
            pw1: (0..=p.bmax).map(|_| z()).collect(),
            pw2: (0..=p.cmax.max(1)).map(|_| z()).collect(),
            unit,
            pairs,
            q: p.groups.iter().map(|g| if g.c > 0 { z() } else { Vec::new() }).collect(),
            ex: p.groups.iter().map(|g| shifted(&g.xp, x0d)).collect(),
        }
    }

    fn base(&self, g: &Group) -> &[Fl<T, 1>] {
        base_of(g, &self.unit, &self.pw0, &self.pw1, &self.pairs)
    }

    /// orders k of the powers, products and the F''-power products, with F''_k as it stands (zero for k >= 1)
    fn partial(&mut self, k: usize) {
        let p = self.p;
        for a in 1..=p.amax {
            self.pw0[a][k] = if a == 1 { self.f[0][k] } else { cs(&self.pw0[a - 1], &self.f[0], k) };
        }
        for b in 1..=p.bmax {
            self.pw1[b][k] = if b == 1 { self.f[1][k] } else { cs(&self.pw1[b - 1], &self.f[1], k) };
        }
        self.pw2[1][k] = self.f[2][k];
        for c in 2..=p.cmax {
            let v = cs(&self.pw2[c - 1], &self.f[2], k);
            self.pw2[c][k] = v;
        }
        let keys: Vec<(usize, usize)> = self.pairs.keys().copied().collect();
        for (a, b) in keys {
            let v = cs(&self.pw0[a], &self.pw1[b], k);
            self.pairs.get_mut(&(a, b)).unwrap()[k] = v;
        }
        for (gi, g) in p.groups.iter().enumerate() {
            if g.c > 0 {
                let v = cs(self.base(g), &self.pw2[g.c], k);
                self.q[gi][k] = v;
            }
        }
    }

    /// order k of N from the current series
    fn n_order(&self, k: usize) -> Fl<T, 1> {
        let mut r = Fl::<T, 1>::zero();
        for (gi, g) in self.p.groups.iter().enumerate() {
            let src: &[Fl<T, 1>] = if g.c > 0 { &self.q[gi] } else { self.base(g) };
            for (j, &e) in self.ex[gi].iter().enumerate() {
                if j <= k {
                    r = r + src[k - j].muls(e);
                }
            }
        }
        r
    }

    /// dN/dF'' at order 0
    fn a0(&self) -> Fl<T, 1> {
        let mut a0 = Fl::<T, 1>::zero();
        for (gi, g) in self.p.groups.iter().enumerate() {
            if g.c == 0 {
                continue;
            }
            let pc = if g.c == 1 { Fl::splat(DD::ONE) } else { self.pw2[g.c - 1][0] };
            a0 = a0 + (self.base(g)[0] * pc).muls(self.ex[gi][0]).mulf(g.c as f64);
        }
        a0
    }

    /// set F''_k = dk and add the part of the F'' powers and their products that holds it
    fn fix(&mut self, k: usize, dk: Fl<T, 1>) {
        self.f[2][k] = dk;
        let t2 = self.f[2][0] * dk;
        let corr = [Fl::zero(), dk, t2.mulf(2.0), (self.f[2][0] * t2).mulf(3.0)];
        self.pw2[1][k] = dk;
        for c in 2..=self.p.cmax {
            self.pw2[c][k] = self.pw2[c][k] + corr[c];
        }
        for gi in 0..self.p.groups.len() {
            let g = &self.p.groups[gi];
            if g.c > 0 {
                let b0 = self.base(g)[0];
                self.q[gi][k] = self.q[gi][k] + b0 * corr[g.c];
            }
        }
    }
}

/// Taylor coefficients 1..=kt of (F, F', F'') about x0 from the order-0 values in `w.c`, solving order k of N = 0 for
/// F''_k. `skip` leaves F''_k = 0 at that order (a control: N then drifts).
#[cfg(not(feature = "implicit_lib"))]
pub fn implicit_taylor<T: Lane>(p: &PolyN, x0: f64, w: &mut Work<T, 1>, skip: Option<usize>) {
    let kt = w.kt;
    let mut cx = Cx::new(p, x0, kt, [w.c[0][0], w.c[1][0], w.c[2][0]]);
    cx.partial(0);
    let inv = Fl::<T, 1>::splat(DD::ONE) / cx.a0();
    for k in 1..=kt {
        let ik = DD::ONE / DD::from_f64(k as f64);
        cx.f[0][k] = cx.f[1][k - 1].muls(ik);
        cx.f[1][k] = cx.f[2][k - 1].muls(ik);
        cx.f[2][k] = Fl::zero();
        cx.partial(k);
        let r = cx.n_order(k);
        let dk = if skip == Some(k) { Fl::zero() } else { -(r * inv) };
        cx.fix(k, dk);
    }
    for i in 0..=kt {
        for (c, fc) in cx.f.iter().enumerate() {
            w.c[c][i] = fc[i];
        }
    }
}

fn base_of<'a, T: Lane>(
    g: &Group,
    unit: &'a [Fl<T, 1>],
    pw0: &'a [Vec<Fl<T, 1>>],
    pw1: &'a [Vec<Fl<T, 1>>],
    pairs: &'a BTreeMap<(usize, usize), Vec<Fl<T, 1>>>,
) -> &'a [Fl<T, 1>] {
    match (g.a, g.b) {
        (0, 0) => unit,
        (a, 0) => &pw0[a],
        (0, b) => &pw1[b],
        (a, b) => &pairs[&(a, b)],
    }
}

#[cfg(feature = "implicit_lib")]
pub use crate::poly_lib::implicit_taylor;
