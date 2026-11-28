test_that("LModularity is defunct", {
  expect_error(LModularity(diag(3)), "removed in evolqg 1.0")
  expect_error(LModularity(diag(3)), "igraph::cluster_optimal")
  expect_s3_class(tryCatch(LModularity(), error = identity), "defunctError")
})
