test_that("every family is standardized to the requested mean and sd", {
  set.seed(9)
  for (family in distribution_families) {
    z <- sample_distribution(family, 20000, mean = 1, sd = 2)
    expect_equal(mean(z), 1, tolerance = 0.05)
    expect_equal(stats::sd(z), 2, tolerance = 0.05)
  }
  expect_error(sample_distribution("Cauchy", 10), "Unknown")
  expect_error(sample_distribution("Normal", 10, sd = 0), "sd")
})
