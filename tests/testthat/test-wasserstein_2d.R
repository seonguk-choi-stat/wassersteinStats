test_that("a translated cloud is matched to itself at distance |v|", {
  set.seed(6)
  x <- matrix(rnorm(40), ncol = 2)
  y <- sweep(x, 2, c(3, 4), "+")
  expect_equal(wasserstein_2d(x, x), 0)
  expect_equal(wasserstein_2d(x, y), 5)
  expect_equal(transport_map_2d(x, y), y)
})

test_that("the assignment is optimal on a small example", {
  x <- cbind(c(0, 1), c(0, 0))
  y <- cbind(c(1, 0), c(1, 1))
  expect_equal(transport_map_2d(x, y), y[c(2, 1), ])
  expect_equal(wasserstein_2d(x, y), 1)
})

test_that("the assignment agrees with brute force", {
  set.seed(7)
  x <- matrix(runif(10), ncol = 2)
  y <- matrix(runif(10), ncol = 2)
  perms <- function(v) {
    if (length(v) <= 1) return(list(v))
    do.call(c, lapply(seq_along(v), function(i)
      lapply(perms(v[-i]), function(r) c(v[i], r))))
  }
  costs <- vapply(perms(1:5), function(s) mean(rowSums((x - y[s, ])^2)),
                  numeric(1))
  expect_equal(wasserstein_2d(x, y)^2, min(costs))
})

test_that("unequal sizes are rejected", {
  expect_error(wasserstein_2d(matrix(0, 2, 2), matrix(0, 3, 2)), "same number")
  expect_error(wasserstein_2d(matrix(0, 2, 3), matrix(0, 2, 3)), "two columns")
})

test_that("the mean of symmetric translates is the base cloud", {
  set.seed(8)
  base <- matrix(rnorm(30), ncol = 2)
  samples <- list(sweep(base, 2, c(2, 0), "+"), sweep(base, 2, c(-2, 0), "+"))
  fit <- frechet_mean_2d(samples, init = base)
  expect_equal(fit$mean, base)
  expect_true(fit$converged)
})
