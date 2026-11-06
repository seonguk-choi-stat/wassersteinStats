test_that("the transport module computes the distance and the geodesic", {
  testServer(mod_transport_1d_server, {
    session$setInputs(
      src_family = "Normal", src_mean = 0, src_sd = 1,
      tgt_family = "Normal", tgt_mean = 2, tgt_sd = 1,
      n = 200, p = "2", t = 0.25, seed = 1
    )
    expect_length(samples()$x, 200)
    expect_equal(distance(), wasserstein_1d(samples()$x, samples()$y))
    # the geodesic moves at constant speed
    expect_equal(wasserstein_1d(samples()$x, interpolated()),
                 0.25 * distance())
    expect_match(output$distance, "W2\\(source, target\\)")

    session$setInputs(p = "1")
    expect_equal(distance(), wasserstein_1d(samples()$x, samples()$y, p = 1))
  })
})
