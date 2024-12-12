module Main where
import Text.Printf

-- An implementation of the Gauss-Legendre algorithm
-- https://en.wikipedia.org/wiki/Gauss%E2%80%93Legendre_algorithm
gl_approx :: (Floating n) => Int -> n
gl_approx iters = gl_approx_impl iters (1.0, sqrt (1.0 / 2.0), 1.0, 1.0 / 4.0)

gl_approx_impl :: (Floating n) => Int -> (n, n, n, n) -> n
gl_approx_impl 0 (a, b, p, t) = (a + b) ^^ 2 / (4.0 * t)
gl_approx_impl i (a, b, p, t) =
    gl_approx_impl (i - 1) (a', b', p', t')
    where
    a' = (a + b) / 2.0
    t' = t - (p * ((a' - a) ^^ 2))
    b' = sqrt (a * b)
    p' = 2.0 * p


-- An implementation of the [Bailey–Borwein–Plouffe formula](https://en.wikipedia.org/wiki/Bailey%E2%80%93Borwein%E2%80%93Plouffe_formula)
-- `iters` is a measure of *hexadecimal* digit precsision.
bbp :: (Floating n) => Int -> n
bbp iters = bbp_impl iters 0 

bbp_impl :: (Floating n) => Int -> n -> n
bbp_impl 0 _ = 0.0
bbp_impl iters k =
    (1 / (16 ** k)) * (a - b - c - d) + next 
    where
    a = 4.0 / (8.0 * k + 1.0)
    b = 2.0 / (8.0 * k + 4.0)
    c = 1.0 / (8.0 * k + 5.0)
    d = 1.0 / (8.0 * k + 6.0)
    next = bbp_impl (iters - 1) (k + 1.0)

main = printf "GL: %f\nBBP: %f\n" (gl_approx 4 :: Double) (bbp 11 :: Double )
