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
