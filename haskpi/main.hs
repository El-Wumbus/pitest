module Main where

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


main = putStrLn . show $ gl_approx 4
