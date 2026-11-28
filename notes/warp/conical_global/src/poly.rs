//! A polynomial N(F0, F1, F2; x) = sum coef F0^a F1^b F2^c x^m with c <= 3, and the Taylor recursion of the implicit
//! equation N = 0 for the state (F, F', F'') on `multishoot::implicit`.
//!
//! About a step's start x0 the fields are series in h = x - x0. The recursion takes order k >= 1 of N as zero and solves
//! it for the coefficient F''_k: F''_k enters the order k of N linearly, with the factor A0 = dN/dF'' at the order-0 values,
//! and F_k, F'_k are known from F'_(k-1), F''_(k-1). The order-0 term of N is not imposed, so N is constant along the
//! flow (N' = 0) and the data's N = 0 carries through. One division, by A0, per step; no hand-differentiated F'''.
use multishoot::implicit::{taylor_shift, Coefficients, ImplicitControl, Jets, Mono, Recursion};
use multishoot::lane::Lane;
use multishoot::{Fl, Work, DD};
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

/// N = sum coef F0^a F1^b F2^c x^m as the library's source: the jet of each group's polynomial in x by Taylor shift
struct PolyNSrc<'a>(&'a PolyN);

impl<T: Lane> Coefficients<Fl<T, 1>> for PolyNSrc<'_> {
    fn pattern(&self) -> Vec<Mono> {
        self.0.groups.iter().map(|g| (g.a, g.b, g.c)).collect()
    }
    fn jets(&self, x0: f64, _kt: usize, out: &mut Jets<Fl<T, 1>>) {
        for (gi, g) in self.0.groups.iter().enumerate() {
            let deg = g.xp.iter().map(|&(m, _)| m).max().unwrap_or(0);
            let mut c = vec![Fl::<T, 1>::zero(); deg + 1];
            for &(m, co) in &g.xp {
                c[m] = c[m] + Fl::splat(co);
            }
            taylor_shift(&c, x0, &mut out.ent[gi].c);
            out.ent[gi].n = c.len();
        }
    }
}

/// Taylor coefficients 1..=kt of (F, F', F'') about x0 from the order-0 values in `w.c`, solving order k of N = 0 for
/// F''_k. `skip` leaves F''_k = 0 at that order (a control: N then drifts).
pub fn implicit_taylor<T: Lane>(p: &PolyN, x0: f64, w: &mut Work<T, 1>, skip: Option<usize>) {
    let kt = w.kt;
    let src = PolyNSrc(p);
    let pat = Coefficients::<Fl<T, 1>>::pattern(&src);
    let mut rec = Recursion::<Fl<T, 1>>::new(&pat, kt);
    if let Some(k) = skip {
        rec.set_control(ImplicitControl::SkipOrder(k));
    }
    let mut jets = Jets::<Fl<T, 1>>::new(pat.len(), kt);
    src.jets(x0, kt, &mut jets);
    let s = rec.form(&jets, 1.0, [w.c[0][0], w.c[1][0], w.c[2][0]]);
    for i in 0..=kt {
        for c in 0..3 {
            w.c[c][i] = s.c[c][i];
        }
    }
}
