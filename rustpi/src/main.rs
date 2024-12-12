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
    pure_gl_approx_impl(iter, (1.0, 1.0/2.0f64.sqrt(), 1.0, 1.0 / 4.0))
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

fn main() {
    println!("std consts PI: {}", std::f64::consts::PI);
    for i in 1..=4  {
        println!("[Imperative] Gauss-Legemdre approx ({i} iters): {}", gl_approx(i));
    }
    
    for i in 1..=4  {
        println!("[Pure functional] Gauss-Legemdre approx ({i} iters): {}", pure_gl_approx(i));
    }
}
