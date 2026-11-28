test_that("DeltaZCorr returns one correlation per skewer for two matrices",
{
  local_seed(41)
  x <- RandomMatrix(6); y <- RandomMatrix(6)
  skewers <- matrix(rnorm(6 * 4), 6, 4, dimnames = list(NULL, paste0("s", 1:4)))
  res <- DeltaZCorr(x, y, skewers)
  expect_equal(names(res), paste0("s", 1:4))
  manual <- vapply(1:4, function(k) {
    a <- x %*% skewers[, k]; b <- y %*% skewers[, k]
    sum(a * b) / sqrt(sum(a^2) * sum(b^2))
  }, 1)
  expect_equal(res, manual, ignore_attr = TRUE)
  # a plain vector is one skewer
  expect_equal(DeltaZCorr(x, y, skewers[, 1]), manual[1], ignore_attr = TRUE)
  expect_length(DeltaZCorr(x, y, skewers[, 1]), 1)
  expect_error(DeltaZCorr(x, y, skewers[1:5, ]), "one row per trait")
  expect_error(DeltaZCorr(x, RandomMatrix(5), skewers), "same dimension")
})

test_that("DeltaZCorr on a list and a matrix has one column per skewer",
{
  local_seed(42)
  mats <- RandomMatrix(6, 4); y <- RandomMatrix(6)
  skewers <- matrix(rnorm(18), 6, 3)
  res <- DeltaZCorr(mats, y, skewers)
  expect_s3_class(res, "data.frame")
  expect_equal(dim(res), c(4, 3))
  expect_false(anyNA(res))
  expect_equal(unlist(res[2, ]), DeltaZCorr(mats[[2]], y, skewers), ignore_attr = TRUE)
  named <- mats; names(named) <- letters[1:4]
  res.named <- DeltaZCorr(named, y, skewers)
  expect_equal(res.named$.id, letters[1:4])
  expect_equal(dim(res.named), c(4, 4))
  one <- DeltaZCorr(mats, y, skewers[, 1])
  expect_equal(dim(one), c(4, 1))
  expect_equal(one[[1]], res[[1]])
})

test_that("DeltaZCorr on a list returns one matrix per skewer",
{
  local_seed(43)
  mats <- RandomMatrix(6, 4)
  names(mats) <- letters[1:4]
  skewers <- matrix(rnorm(18), 6, 3, dimnames = list(NULL, c("p", "q", "r")))
  res <- DeltaZCorr(mats, skewers = skewers)
  expect_equal(dim(res), c(4, 4, 3))
  expect_equal(dimnames(res), list(letters[1:4], letters[1:4], c("p", "q", "r")))
  expect_equal(res[3, 2, ], DeltaZCorr(mats[[3]], mats[[2]], skewers))
  expect_equal(res[4, 1, "q"], DeltaZCorr(mats[[4]], mats[[1]], skewers)[["q"]])
  expect_true(all(apply(res, 3, function(m) all(m[upper.tri(m, diag = TRUE)] == 0))))
  # skewer 3 and up are no longer dropped, and skewer 2 is not labelled as probabilities
  expect_false(isTRUE(all.equal(res[, , 2], res[, , 3])))

  one <- DeltaZCorr(mats, skewers = skewers[, 1])
  expect_type(one, "list")
  expect_equal(names(one), c("correlations", "probabilities"))
  expect_equal(one$correlations, res[, , 1])
  expect_equal(dimnames(one$correlations), list(letters[1:4], letters[1:4]))
  expect_error(DeltaZCorr(mats, skewers = skewers[1:3, ]), "one row per trait")
})
