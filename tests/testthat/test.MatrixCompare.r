test_that("MatrixCompare returns one row per method",
{
  local_seed(51)
  x <- RandomMatrix(6, 1, 1, 10); y <- RandomMatrix(6, 1, 10, 20)
  res <- MatrixCompare(x, y)
  expect_s3_class(res, "data.frame")
  expect_equal(res$.id, c("RS", "Mantel", "Krz", "PCA"))
  expect_equal(dim(res), c(4, 3))
  expect_equal(res$correlation[3], KrzCor(x, y)[[1]])
  expect_equal(res$correlation[4], PCAsimilarity(x, y)[[1]])
  expect_equal(res$correlation[2], MantelCor(cov2cor(x), cov2cor(y))[[1]], tolerance = 1e-12)
  expect_true(is.na(res$probability[3]) && is.na(res$probability[4]))
  expect_equal(names(MatrixCompare(x, y, id = "method"))[3], "method")
})
