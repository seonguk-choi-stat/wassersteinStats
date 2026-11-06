test_that("the Frechet module's descent reaches the closed-form barycenter", {
  testServer(mod_frechet_1d_server, {
    session$setInputs(scenario = "Shifted normals", k = 3, n = 100,
                      step = 1, seed = 1, iter = 0)
    expect_length(samples(), 3)
    expect_true(fit()$converged)
    expect_equal(fit()$mean, barycenter())

    session$setInputs(step = 0.3, iter = 2)
    expect_equal(current(), fit()$path[[3]])
    expect_match(output$status, "Showing iteration 2")
  })
})
