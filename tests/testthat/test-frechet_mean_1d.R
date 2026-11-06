test_that("the mean of a sample and its translate is the midpoint translate", {
  x <- c(0.3, -1.2, 2.5, 0.9, -0.4)
  expect_equal(frechet_mean_1d(list(x, x + 2), n_grid = length(x)),
               sort(x) + 1)
})

test_that("the mean of a single sample is that sample", {
  x <- c(5, 1, 3, 2)
  expect_equal(frechet_mean_1d(list(x), n_grid = 4), sort(x))
})

test_that("descent with unit step reaches the closed form in one step", {
  set.seed(4)
  samples <- list(rnorm(80, -2), rnorm(60, 2), rexp(100))
  init <- rnorm(50)
  fit <- frechet_mean_1d_descent(samples, init = init, step = 1)
  expect_equal(fit$mean, frechet_mean_1d(samples, n_grid = 50))
  expect_true(fit$converged)
  expect_length(fit$path, 2)
})

test_that("smaller steps decrease the gradient norm geometrically", {
  set.seed(5)
  samples <- list(rnorm(40, -3), rnorm(40, 3))
  fit <- frechet_mean_1d_descent(samples, step = 0.5, max_iter = 10, tol = 0)
  ratios <- fit$grad_norm[-1] / fit$grad_norm[-length(fit$grad_norm)]
  expect_equal(ratios, rep(0.5, length(ratios)))
})

test_that("invalid input is rejected", {
  expect_error(frechet_mean_1d(list()), "non-empty")
  expect_error(frechet_mean_1d(list(1:3), n_grid = 0), "n_grid")
  expect_error(frechet_mean_1d_descent(list(1:3), step = 2), "step")
})
