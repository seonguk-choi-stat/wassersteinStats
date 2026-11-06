test_that("wasserstein_1d is zero on identical samples and symmetric", {
  set.seed(1)
  x <- rnorm(30)
  y <- rexp(45)
  expect_equal(wasserstein_1d(x, x), 0)
  expect_equal(wasserstein_1d(x, y), wasserstein_1d(y, x))
})

test_that("a translation by c has distance |c| for every p", {
  x <- c(1, 4, 2, 8, 5)
  for (p in c(1, 2, 3)) {
    expect_equal(wasserstein_1d(x, x + 2, p = p), 2)
  }
})

test_that("equal sizes reduce to the sorted matching", {
  set.seed(2)
  x <- rnorm(50)
  y <- rgamma(50, 2)
  expect_equal(wasserstein_1d(x, y, p = 2),
               sqrt(mean((sort(x) - sort(y))^2)))
  expect_equal(wasserstein_1d(x, y, p = 1), mean(abs(sort(x) - sort(y))))
})

test_that("unequal sizes are computed exactly", {
  # Q_x = 0 on (0, 1/2], 1 on (1/2, 1]; Q_y = 0, 1/2, 1 on thirds.
  # |Q_x - Q_y| = 1/2 on (1/3, 2/3) and 0 elsewhere.
  x <- c(0, 1)
  y <- c(0, 0.5, 1)
  expect_equal(wasserstein_1d(x, y, p = 1), 1 / 6)
  expect_equal(wasserstein_1d(x, y, p = 2), sqrt(1 / 12))
})

test_that("normal laws match the closed form for W2", {
  u <- (1:4000 - 0.5) / 4000
  w <- wasserstein_1d(qnorm(u, 0, 1), qnorm(u, 2, 3))
  expect_equal(w, sqrt(2^2 + (3 - 1)^2), tolerance = 1e-2)
})

test_that("invalid input is rejected", {
  expect_error(wasserstein_1d(numeric(0), 1), "non-empty")
  expect_error(wasserstein_1d(c(1, NA), 1), "finite")
  expect_error(wasserstein_1d(1, 2, p = 0.5), "p")
})

test_that("transport_map_1d is the sorted matching and keeps input order", {
  expect_equal(transport_map_1d(c(3, 1, 2), c(10, 30, 20)), c(30, 10, 20))
  # ties in the source are split across distinct targets
  expect_setequal(transport_map_1d(c(1, 1), c(5, 6)), c(5, 6))
})

test_that("exp_map_1d inverts log_map_1d and traces the geodesic", {
  set.seed(3)
  gamma <- rnorm(40)
  mu <- rnorm(40, 3, 2)
  v <- log_map_1d(gamma, mu)
  expect_equal(sort(exp_map_1d(gamma, v)), sort(mu))
  expect_equal(exp_map_1d(gamma, v, t = 0), gamma)
  # constant speed: W2(gamma, geodesic(t)) = t * W2(gamma, mu)
  mid <- exp_map_1d(gamma, v, t = 0.3)
  expect_equal(wasserstein_1d(gamma, mid), 0.3 * wasserstein_1d(gamma, mu))
  expect_error(exp_map_1d(gamma, v[-1]), "same length")
})
