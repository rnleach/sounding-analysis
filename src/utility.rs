// FIXME: do not export this from the crate. Used in sonde currently, but generating profiles can
// be done in this crate.
/// Bisection algorithm for finding the root of an equation given values bracketing a root. Used
/// when drawing moist adiabats.
pub fn find_root(f: &Fn(f64) -> f64, mut low_val: f64, mut high_val: f64) -> f64 {
    use std::f64;
    const MAX_IT: usize = 50;
    const EPS: f64 = 1.0e-10;

    if low_val > high_val {
        ::std::mem::swap(&mut low_val, &mut high_val);
    }

    let mut f_low = f(low_val);
    // let mut f_high = f(high_val);

    let mut mid_val = (high_val - low_val) / 2.0 + low_val;
    let mut f_mid = f(mid_val);
    for _ in 0..MAX_IT {
        if f_mid * f_low > 0.0 {
            low_val = mid_val;
            f_low = f_mid;
        } else {
            high_val = mid_val;
            // f_high = f_mid;
        }

        if (high_val - low_val).abs() < EPS {
            break;
        }

        mid_val = (high_val - low_val) / 2.0 + low_val;
        f_mid = f(mid_val);
    }

    mid_val
}

// TODO: unit conversions.

#[cfg(test)]
pub mod test_tools {
    pub fn approx_equal(val1: f64, val2: f64, eps: f64) -> bool {
        assert!(eps > 0.0);

        (val1 - val2).abs() < eps
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use utility::test_tools::*;

    #[test]
    fn test_find_root() {
        assert!(approx_equal(
            1.0,
            find_root(&|x| x * x - 1.0, 2.0, 0.0),
            1.0e-10
        ));
        assert!(approx_equal(
            -1.0,
            find_root(&|x| x * x - 1.0, -2.0, 0.0),
            1.0e-10
        ));
    }
}
