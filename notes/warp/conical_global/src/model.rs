//! The equations on multishoot's trait: the implicit cubic N = 0 with state (F, F', F''), and the end maps.
//!
//! Data: left (C-, e), right (C+), global F_inf. The left state at -200 and the right state at +200 are the series of
//! `tables::SERIES` at |eta| = 200 (the mirror at the left), and e shifts F'' at the left end only: the matching returns
//! e = 0, which absorbs the one equation that N = 0 at both ends makes redundant.
use crate::poly::{implicit_taylor, PolyN};
use crate::tables::{LEND, LN_LEND, N_TERMS, SERIES};
use multishoot::lane::Lane;
use multishoot::{Equations, Fl, State, Work, DD};

/// deliberate errors, one per control; all off by default
#[derive(Clone, Copy, Debug, Default)]
pub struct Mutation {
    pub eq_coeff: bool,
    pub series_coeff: bool,
    pub match_sign: bool,
    pub plant_locator: bool,
    pub n_drift: bool,
}

impl Mutation {
    pub fn from_name(s: &str) -> Option<Mutation> {
        let mut m = Mutation::default();
        match s {
            "eq-coeff" => m.eq_coeff = true,
            "series-coeff" => m.series_coeff = true,
            "match-sign" => m.match_sign = true,
            "plant-locator" => m.plant_locator = true,
            "n-drift" => m.n_drift = true,
            _ => return None,
        }
        Some(m)
    }
}

/// N of det_classify.conical_N() as a polynomial
pub fn conical_poly(m: &Mutation) -> PolyN {
    let mut t = N_TERMS.to_vec();
    if m.eq_coeff {
        // (1, 2, 1, 2, -6) -> -5
        for term in t.iter_mut() {
            if (term.0, term.1, term.2, term.3) == (1, 2, 1, 2) {
                term.4 += 1;
            }
        }
    }
    PolyN::from_ints(&t)
}

/// G = F, p = eta F', q = eta^2 F'' + eta F' of the series at |eta| = 200 from (F_inf, C). The series is
/// F_inf + sum a_(n,j) nu^n T^j, T = ln eta, nu = 1/eta, theta = eta d/d eta.
pub fn series_gpq<T: Lane>(fi: Fl<T, 1>, c: Fl<T, 1>, m: &Mutation) -> [Fl<T, 1>; 3] {
    series_gpq_at(fi, c, LEND, m)
}

/// `series_gpq` at |eta| = l
pub fn series_gpq_at<T: Lane>(fi: Fl<T, 1>, c: Fl<T, 1>, l: f64, m: &Mutation) -> [Fl<T, 1>; 3] {
    let ln_l = if l == LEND { LN_LEND } else { DD::from_f64(l).ln() };
    let mx = SERIES.iter().flat_map(|t| t.mons.iter().map(|x| x.0 as usize)).max().unwrap();
    let ml = SERIES.iter().flat_map(|t| t.mons.iter().map(|x| x.1 as usize)).max().unwrap();
    let ms = SERIES.iter().map(|t| t.s as usize).max().unwrap();
    let one = Fl::<T, 1>::splat(DD::ONE);
    let inv = one / fi;
    let pw = |b: Fl<T, 1>, n: usize| {
        let mut v = vec![one];
        for i in 0..n {
            v.push(v[i] * b);
        }
        v
    };
    let (fp, cp, ip) = (pw(fi, mx), pw(c, ml), pw(inv, ms));
    let nu = DD::ONE / DD::from_f64(l);
    let tj = |e: i32| if e < 0 { DD::ZERO } else { ln_l.powi(e) };
    let (mut g, mut p, mut q) = (fi, Fl::zero(), Fl::zero());
    for t in SERIES.iter() {
        let mut a = Fl::<T, 1>::zero();
        for &(mm, l, co) in t.mons {
            a = a + (fp[mm as usize] * cp[l as usize]).muls(co);
        }
        a = a * ip[t.s as usize];
        if m.series_coeff && t.n == 4 && t.j == 1 {
            a = a.muls(DD::from_f64(1.01));
        }
        let (n, j) = (t.n as f64, t.j as f64);
        let nun = nu.powi(t.n as i32);
        let jj = t.j as i32;
        let kg = nun * tj(jj);
        let kp = nun * (DD::from_f64(-n) * tj(jj) + DD::from_f64(j) * tj(jj - 1));
        let kq = nun
            * (DD::from_f64(n * n) * tj(jj) - DD::from_f64(2.0 * n * j) * tj(jj - 1)
                + DD::from_f64(j * (j - 1.0)) * tj(jj - 2));
        g = g + a.muls(kg);
        p = p + a.muls(kp);
        q = q + a.muls(kq);
    }
    [g, p, q]
}

/// the state (F, F', F'') at eta = +-200 from the series; sign = +1 at +200, -1 at -200 (F' flips)
pub fn end_state<T: Lane>(fi: Fl<T, 1>, c: Fl<T, 1>, sign: f64, m: &Mutation) -> State<T, 1> {
    end_state_at(fi, c, sign, LEND, m)
}

/// `end_state` at |eta| = l
pub fn end_state_at<T: Lane>(fi: Fl<T, 1>, c: Fl<T, 1>, sign: f64, l: f64, m: &Mutation) -> State<T, 1> {
    let [g, p, q] = series_gpq_at(fi, c, l, m);
    let nu = DD::ONE / DD::from_f64(l);
    let f1 = p.muls(nu);
    let f2 = (q - p).muls(nu * nu);
    let f1 = if sign < 0.0 && !m.match_sign { -f1 } else { f1 };
    vec![g, f1, f2]
}

/// the end maps of a model: data at each end to the state at the end node
pub trait Ends: Sync {
    fn left<T: Lane>(&self, x: f64, data: &[Fl<T, 1>], g: &[T]) -> State<T, 1>;
    fn right<T: Lane>(&self, x: f64, data: &[Fl<T, 1>], g: &[T]) -> State<T, 1>;
}

pub struct SeriesEnds {
    pub m: Mutation,
}

impl Ends for SeriesEnds {
    fn left<T: Lane>(&self, x: f64, data: &[Fl<T, 1>], g: &[T]) -> State<T, 1> {
        assert!((x + LEND).abs() < 1e-9, "left end at {x}");
        let mut s = end_state(Fl([g[0]]), data[0], -1.0, &self.m);
        s[2] = s[2] + data[1];
        s
    }
    fn right<T: Lane>(&self, x: f64, data: &[Fl<T, 1>], g: &[T]) -> State<T, 1> {
        assert!((x - LEND).abs() < 1e-9, "right end at {x}");
        end_state(Fl([g[0]]), data[0], 1.0, &self.m)
    }
}

/// no ends: for walks that never use them
pub struct NoEnds;

impl Ends for NoEnds {
    fn left<T: Lane>(&self, _x: f64, _d: &[Fl<T, 1>], _g: &[T]) -> State<T, 1> {
        unreachable!("no left end")
    }
    fn right<T: Lane>(&self, _x: f64, _d: &[Fl<T, 1>], _g: &[T]) -> State<T, 1> {
        unreachable!("no right end")
    }
}

/// N = 0 with state (F, F', F''): `n` is N, `ends` the maps from data to the end states
pub struct Implicit<E: Ends> {
    pub n: PolyN,
    pub ends: E,
    pub skip_order: Option<usize>,
}

impl<E: Ends> Equations<1> for Implicit<E> {
    const NF: usize = 3;
    const NL: usize = 2;
    const NR: usize = 1;
    const NG: usize = 1;
    const POINTWISE: bool = true;
    type Prep = ();

    fn taylor_coeffs<T: Lane>(&self, x0: f64, _g: &[T], w: &mut Work<T, 1>) {
        implicit_taylor(&self.n, x0, w, self.skip_order);
    }

    fn left_state<T: Lane>(&self, x: f64, data: &[Fl<T, 1>], g: &[T], _prep: &()) -> State<T, 1> {
        self.ends.left(x, data, g)
    }

    fn right_state<T: Lane>(&self, x: f64, data: &[Fl<T, 1>], g: &[T], _prep: &()) -> State<T, 1> {
        self.ends.right(x, data, g)
    }
}

pub fn conical(m: Mutation) -> Implicit<SeriesEnds> {
    Implicit { n: conical_poly(&m), ends: SeriesEnds { m }, skip_order: if m.n_drift { Some(3) } else { None } }
}
