// Throughput of a stand-in for the planned dual-null cell update (#39 design; NOT the evolution code).
// rustc -C opt-level=3 -C target-cpu=native kernel_bench.rs -o kernel_bench && ./kernel_bench [N] [reps]
// Same arithmetic as kernel_bench.py: predictor-corrector diamond update, fields r, f, ghost, normal, n = 2, row sweep.
use std::time::Instant;

const NN: f64 = 2.0;
const E1: f64 = -1.0;
const E2: f64 = 1.0;

fn sweep(n: usize, h: f64, rp: &[f64], ap: &[f64], rm: &[f64], am: &[f64]) -> f64 {
    let mut old = vec![[0.0f64; 4]; n + 1];
    let mut new = vec![[0.0f64; 4]; n + 1];
    for i in 0..=n {
        old[i][0] = rp[i];
        old[i][2] = ap[i];
    }
    let h2 = h * h;
    for j in 0..n {
        new[0] = [rm[j + 1], 0.0, am[j + 1], 0.0];
        for i in 0..n {
            let s = old[i];
            let e = old[i + 1];
            let w = new[i];
            let f = |c: [f64; 4], u: [f64; 4], v: [f64; 4]| -> [f64; 4] {
                let emf = (-c[1]).exp();
                let rpm = -(NN - 1.0) * (u[0] * v[0] + 0.5 * emf) / c[0];
                let spm = E1 * u[2] * v[2] + E2 * u[3] * v[3];
                let fpm = (NN / c[0]) * rpm + spm;
                let k = -NN / (2.0 * c[0]);
                [rpm, fpm, k * (u[0] * v[2] + v[0] * u[2]), k * (u[0] * v[3] + v[0] * u[3])]
            };
            let mut c = [0.0; 4];
            let mut u = [0.0; 4];
            let mut v = [0.0; 4];
            for q in 0..4 {
                c[q] = 0.5 * (e[q] + w[q]);
                u[q] = (e[q] - s[q]) / h;
                v[q] = (w[q] - s[q]) / h;
            }
            let g = f(c, u, v);
            let mut nh = [0.0; 4];
            for q in 0..4 {
                nh[q] = e[q] + w[q] - s[q] + h2 * g[q];
            }
            for q in 0..4 {
                c[q] = 0.25 * (s[q] + e[q] + w[q] + nh[q]);
                u[q] = (e[q] + nh[q] - s[q] - w[q]) / (2.0 * h);
                v[q] = (w[q] + nh[q] - s[q] - e[q]) / (2.0 * h);
            }
            let g = f(c, u, v);
            for q in 0..4 {
                new[i + 1][q] = e[q] + w[q] - s[q] + h2 * g[q];
            }
        }
        std::mem::swap(&mut old, &mut new);
    }
    old[n][0]
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let n: usize = args.get(1).map(|s| s.parse().unwrap()).unwrap_or(1000);
    let reps: usize = args.get(2).map(|s| s.parse().unwrap()).unwrap_or(3);
    let h = 0.001;
    let x: Vec<f64> = (0..=n).map(|i| i as f64 * h).collect();
    let rp: Vec<f64> = x.iter().map(|v| (1.0 + (v / 2f64.sqrt()).powi(2)).sqrt()).collect();
    let ap: Vec<f64> = x.iter().map(|v| 2f64.sqrt() * (v / 2f64.sqrt()).atan()).collect();
    let rm = rp.clone();
    let am: Vec<f64> = ap.iter().map(|v| -v).collect();
    let load = std::fs::read_to_string("/proc/loadavg").unwrap();
    println!("N = {}, cells = {}, h = {}, load average at start {}, one thread", n, n * n, h, load.split_whitespace().take(3).collect::<Vec<_>>().join(" "));
    let mut ts = vec![];
    let mut val = 0.0;
    for _ in 0..reps {
        let t0 = Instant::now();
        val = sweep(n, h, &rp, &ap, &rm, &am);
        ts.push(t0.elapsed().as_secs_f64());
    }
    ts.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let cells = (n * n) as f64;
    println!("rust   best {:8.3} s  median {:8.3} s  -> {:7.1} ns/cell ({:6.2} Mcell/s)  value {:.6}", ts[0], ts[ts.len() / 2], 1e9 * ts[0] / cells, cells / ts[0] / 1e6, val);
}
