//! poly_lib (feature `implicit_lib`): `implicit_taylor` of `poly.rs` on multishoot's `implicit::Recursion`, on the field
//! `Fl<T, 1>` of any lane (values, duals for the matching's Jacobian), with N's monomials from `PolyN` as jets.
use crate::poly::PolyN;
use multishoot::implicit::{taylor_shift, Coefficients, ImplicitControl, Jets, Mono, Recursion};
use multishoot::lane::Lane;
use multishoot::{Fl, Work};

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
