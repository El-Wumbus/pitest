// https://en.wikipedia.org/wiki/Gauss%E2%80%93Legendre_algorithm
fn gl_approx(iters: usize) -> f64 {
    let mut a = 1.0_f64;
    let mut b = 1.0_f64 / 2.0f64.sqrt();
    let mut p = 1.0_f64;
    let mut t = 1.0_f64 / 4.0_f64;

    for _ in 0..iters {
        let new_a = (a + b) / 2.0_f64;
        t = t - (p * (new_a - a).powf(2.0_f64));
        b = (a * b).sqrt();
        a = new_a;
        p *= 2.0;
    }

    (a + b).powf(2.0) / (4.0_f64 * t)
}

/// A recursive, pure-functional implementation of the [Gauss–Legendre algorithm](https://en.wikipedia.org/wiki/Gauss%E2%80%93Legendre_algorithm)
#[inline]
fn pure_gl_approx(iter: usize) -> f64 {
    pure_gl_approx_impl(iter, (1.0, 1.0 / 2.0f64.sqrt(), 1.0, 1.0 / 4.0))
}

fn pure_gl_approx_impl(iter: usize, (a, b, p, t): (f64, f64, f64, f64)) -> f64 {
    if iter == 0 {
        return (a + b).powf(2.0) / (4.0_f64 * t);
    }

    let new_a = (a + b) / 2.0_f64;
    let t = t - (p * (new_a - a).powf(2.0_f64));
    let b = (a * b).sqrt();
    let a = new_a;
    let p = p * 2.0;

    pure_gl_approx_impl(iter - 1, (a, b, p, t))
}

/// An implementation of the [Bailey–Borwein–Plouffe formula](https://en.wikipedia.org/wiki/Bailey%E2%80%93Borwein%E2%80%93Plouffe_formula)
/// `iters` is a measure of *hexadecimal* digit precsision.
fn bbp(iters: usize) -> f64 {
    let mut s = 0_f64;
    for k in 0..iters {
        let k = k as f64;
        s += (1.0 / 16.0_f64.powf(k))
            * ((4.0 / (8.0 * k + 1.0))
                - (2.0 / (8.0 * k + 4.0))
                - (1.0 / (8.0 * k + 5.0))
                - (1.0 / (8.0 * k + 6.0)));
    }
    s
}

/// A pure-functional implementation of the [Bailey–Borwein–Plouffe formula](https://en.wikipedia.org/wiki/Bailey%E2%80%93Borwein%E2%80%93Plouffe_formula)
/// `iters` is a measure of *hexadecimal* digit precsision.
fn pure_bbp(iters: usize) -> f64 {
    pure_bbp_impl(iters, 0.0)
}

fn pure_bbp_impl(iters: usize, k: f64) -> f64 {
    if iters == 0 {
        return 0.0;
    }

    (1.0 / 16.0_f64.powf(k))
        * ((4.0 / (8.0 * k + 1.0))
            - (2.0 / (8.0 * k + 4.0))
            - (1.0 / (8.0 * k + 5.0))
            - (1.0 / (8.0 * k + 6.0)))
        + pure_bbp_impl(iters - 1, k + 1.0)
}

fn main() {
    println!("std consts PI: {}", std::f64::consts::PI);
    for i in 1..=4 {
        let x = gl_approx(i);
        let y = pure_gl_approx(i);

        println!("Gauss-Legemdre approx ([Imparative] {i} iters): {x}",);
        println!("Gauss-Legemdre approx [Pure functional] ({i} iters): {y}",);
    }

    for i in 1..=20 {
        let x = bbp(i);
        let y = pure_bbp(i);

        println!("BBP ([Imparative] {i} iters): {x}",);
        println!("BBP [Pure functional] ({i} iters): {y}",);
        if y == std::f64::consts::PI {
            break;
        }
    }
}
