test_that("MultivDriftTest returns the expected structure", {
  local_seed(1)
  p1 <- matrix(rnorm(40*3), 40); p2 <- matrix(rnorm(40*3, 0.5), 40)
  G <- cov(p1)
  out <- MultivDriftTest(p1, p2, G, Ne = 1000, generations = 100, iterations = 200)
  expect_named(out, c("magnitude of morphological change expected through drift", "observed magnitude"))
  expect_length(out[[1]], 2)
  expect_named(out[[1]], c("2.5%", "97.5%"))
  expect_length(out[[2]], 2)
  expect_true(all(out[[1]] >= 0) && out[[1]][1] < out[[1]][2])
  # drift of the expected size: wide G gives larger magnitudes
  big <- MultivDriftTest(p1, p2, G * 100, Ne = 1000, generations = 100, iterations = 200)
  expect_true(big[[1]][2] > out[[1]][2])
})

test_that("MultivDriftTest works with a single trait", {
  local_seed(2)
  p1 <- matrix(rnorm(30), 30); p2 <- matrix(rnorm(30, 1), 30)
  out <- MultivDriftTest(p1, p2, matrix(var(p1[, 1])), Ne = 500, generations = 50, iterations = 100)
  expect_length(out[[1]], 2)
  expect_true(all(is.finite(unlist(out))))
  # scalar G and data frames also work
  out2 <- MultivDriftTest(data.frame(x = p1[, 1]), data.frame(x = p2[, 1]), var(p1[, 1]),
                          Ne = 500, generations = 50, iterations = 100)
  expect_true(all(is.finite(unlist(out2))))
})
