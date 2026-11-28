test_that("Matrix Reimann Distance",
{
  local_seed(10, "3.5.0")
  cor.matrix.1 <- RandomMatrix(10, LKJ = FALSE)
  cor.matrix.2 <- RandomMatrix(10, LKJ = FALSE)
  results <- RiemannDist(cor.matrix.1, cor.matrix.2)
  results_rec <- RiemannDist(cor.matrix.2, cor.matrix.1)
  results_ident <- RiemannDist(cor.matrix.1, cor.matrix.1)
  expect_type(results, "double")
  expect_true(results >=  0)
  expect_equal(results, results_rec, tolerance = 10e-4)
  expect_equal(results_ident, 0)
}
)

test_that("Matrix Distribution Overlap Distance",
{
  cor.matrix.1 <- RandomMatrix(10, LKJ = FALSE)
  cor.matrix.2 <- RandomMatrix(10, LKJ = FALSE)
  results <- OverlapDist(cor.matrix.1, cor.matrix.2)
  results_rec <- OverlapDist(cor.matrix.2, cor.matrix.1)
  results_ident <- OverlapDist(cor.matrix.1, cor.matrix.1)
  expect_type(results, "double")
  expect_true(results >=  0)
  expect_equal(results, results_rec, tolerance = 10e-4)
  expect_equal(results_ident, 0)
}
)

