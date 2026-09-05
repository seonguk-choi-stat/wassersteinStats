test_that("wasserstein1d computes correct distances", {
  x <- c(1, 2, 3, 4, 5)
  
  # Distance to itself is 0
  expect_equal(wasserstein1d(x, x), 0)
  
  # Distance to a shifted sample
  y <- x + 2
  expect_equal(wasserstein1d(x, y, p = 2), 2)
})