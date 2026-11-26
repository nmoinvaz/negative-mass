//! The matching: a multiple-shooting layout on the nodes, Newton's method in double-double, and the pin of the global.
//!
//! Unknowns: the left data (C-, e), the states at the nodes, the right data (C+) and the global (F_inf, or alpha of a
//! planted equation), 8 identical lanes (the problem is pointwise, so the lanes are independent). One extra row pins the
//! global. Segments left of `match_node` run outward from the left end, the others inward from the right end.
use crate::walk::S3;
use multishoot::lane::DV;
use multishoot::{Double, Equations, ExtraRow, Fl, Plan, Shooting, DD};
use multisolver::NewtonOpts;

pub struct ConnectIn {
    pub nodes: Vec<f64>,
    pub match_node: usize,
    pub plan: Plan,
    pub node_states: Vec<S3>,
    pub left: [DD; 2],
    pub right: DD,
    pub global: DD,
    pub pin: DD,
    pub tol: f64,
}

pub struct ConnectOut {
    pub converged: bool,
    pub fnorm: Vec<f64>,
    pub left: [DD; 2],
    pub right: DD,
    pub global: DD,
    pub nodes: Vec<S3>,
    pub steps: Vec<usize>,
}

pub fn connect<E: Equations<1>>(eq: E, inp: &ConnectIn) -> Result<ConnectOut, String> {
    let mut sh = Shooting::new(eq, inp.nodes.clone(), inp.plan.kt);
    sh.match_node = inp.match_node;
    sh.extra.push(ExtraRow { terms: vec![(sh.off_global(), DD::ONE)], rhs: inp.pin });
    if !sh.is_square() {
        return Err(format!("not square: {} rows, {} unknowns", sh.rows(), sh.len()));
    }
    let n = sh.len();
    let mut z = vec![vec![0.0; n]; 2];
    let set = |z: &mut Vec<Vec<f64>>, o: usize, v: DD| {
        for j in 0..8 {
            z[0][o + j] = v.hi;
            z[1][o + j] = v.lo;
        }
    };
    set(&mut z, 0, inp.left[0]);
    set(&mut z, 8, inp.left[1]);
    for (i, s) in inp.node_states.iter().enumerate() {
        for f in 0..3 {
            set(&mut z, sh.off_node(i) + f * 8, s[f]);
        }
    }
    set(&mut z, sh.off_right(), inp.right);
    let og = sh.off_global();
    z[0][og] = inp.global.hi;
    z[1][og] = inp.global.lo;
    sh.plan_segments(&z[0], &z[1], &inp.plan);
    let opts = NewtonOpts { tol: inp.tol, ..Default::default() };
    let rep = sh.newton::<Double>(&mut z, &opts)?;
    let get = |z: &Vec<Vec<f64>>, o: usize| Fl::<DV, 1>::load(&z[0][o..o + 8], &z[1][o..o + 8]).get(0);
    let nodes = (0..sh.k())
        .map(|i| {
            let o = sh.off_node(i);
            [get(&z, o), get(&z, o + 8), get(&z, o + 16)]
        })
        .collect();
    Ok(ConnectOut {
        converged: rep.converged,
        fnorm: rep.fnorm,
        left: [get(&z, 0), get(&z, 8)],
        right: get(&z, sh.off_right()),
        global: DD { hi: z[0][og], lo: z[1][og] },
        nodes,
        steps: sh.segs.iter().map(|s| s.len() - 1).collect(),
    })
}
