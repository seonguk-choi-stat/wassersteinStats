test_that("the point-cloud module fits a barycenter of the right shape", {
  testServer(mod_cloud_2d_server, {
    session$setInputs(k = 3, n = 20, spread = 3, step = 1, seed = 1,
                      iter = 0, matched = "Cloud 1")
    expect_length(clouds(), 3)
    expect_equal(dim(fit()$mean), c(20L, 2L))
    expect_length(distances(), 3)
    # the fixed-point iteration does not increase the Frechet functional
    first <- vapply(clouds(), wasserstein_2d, numeric(1), x = fit()$path[[1]])
    last <- vapply(clouds(), wasserstein_2d, numeric(1), x = fit()$mean)
    expect_lte(mean(last^2), mean(first^2) + 1e-10)
  })
})
