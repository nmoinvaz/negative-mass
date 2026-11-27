//! Small dense linear algebra over a lane: the determinant of a 4 x 4 matrix.
use crate::num::Num;

/// det of an n x n matrix by Gaussian elimination with partial pivoting
pub fn det<T: Num, const N: usize>(m: &[[T; N]; N]) -> T {
    let mut a = *m;
    let mut d = T::one();
    for c in 0..N {
        let mut p = c;
        for r in c + 1..N {
            if a[r][c].to_f64().abs() > a[p][c].to_f64().abs() {
                p = r;
            }
        }
        if a[p][c].to_f64() == 0.0 {
            return T::zero();
        }
        if p != c {
            a.swap(p, c);
            d = -d;
        }
        d = d * a[c][c];
        let inv = a[c][c].recip();
        for r in c + 1..N {
            let f = a[r][c] * inv;
            for k in c..N {
                a[r][k] = a[r][k] - f * a[c][k];
            }
        }
    }
    d
}

/// inverse of an n x n matrix by Gauss-Jordan with partial pivoting
pub fn inv<T: Num, const N: usize>(m: &[[T; N]; N]) -> [[T; N]; N] {
    let mut a = *m;
    let mut b = [[T::zero(); N]; N];
    for i in 0..N {
        b[i][i] = T::one();
    }
    for c in 0..N {
        let mut p = c;
        for r in c + 1..N {
            if a[r][c].to_f64().abs() > a[p][c].to_f64().abs() {
                p = r;
            }
        }
        a.swap(p, c);
        b.swap(p, c);
        let iv = a[c][c].recip();
        for k in 0..N {
            a[c][k] = a[c][k] * iv;
            b[c][k] = b[c][k] * iv;
        }
        for r in 0..N {
            if r != c {
                let f = a[r][c];
                if f.to_f64() != 0.0 {
                    for k in 0..N {
                        a[r][k] = a[r][k] - f * a[c][k];
                        b[r][k] = b[r][k] - f * b[c][k];
                    }
                }
            }
        }
    }
    b
}

pub fn matmul<T: Num, const N: usize>(a: &[[T; N]; N], b: &[[T; N]; N]) -> [[T; N]; N] {
    let mut o = [[T::zero(); N]; N];
    for i in 0..N {
        for j in 0..N {
            let mut s = T::zero();
            for k in 0..N {
                s = s + a[i][k] * b[k][j];
            }
            o[i][j] = s;
        }
    }
    o
}

/// The matrix sign function of a matrix with no eigenvalue on the imaginary axis, by the scaled Newton iteration
/// S <- (g S + (g S)^-1) / 2, g = |det S|^(-1/N). Its eigenvalues are +1 on the invariant subspace of the eigenvalues with
/// positive real part and -1 on the other. Returns the sign and the number of iterations.
pub fn matsign<T: Num, const N: usize>(a: &[[T; N]; N]) -> ([[T; N]; N], usize) {
    let mut s = *a;
    let half = T::from_f64(0.5);
    for it in 0..80 {
        let d = det(&s).abs().to_f64();
        let g = if it < 8 { T::from_f64(d.powf(-1.0 / N as f64)) } else { T::one() };
        let si = inv(&s);
        let mut n = [[T::zero(); N]; N];
        let mut diff = 0.0f64;
        let mut big = 0.0f64;
        for i in 0..N {
            for j in 0..N {
                n[i][j] = (g * s[i][j] + si[i][j] * g.recip()) * half;
                diff = diff.max((n[i][j] - s[i][j]).to_f64().abs());
                big = big.max(n[i][j].to_f64().abs());
            }
        }
        s = n;
        if it >= 8 && diff <= 8.0 * T::EPS * big {
            // one more step squares the error away
            return (s, it + 1);
        }
    }
    (s, 80)
}
