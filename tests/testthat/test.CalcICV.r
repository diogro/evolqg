test_that("CalcICV throws error",
{
  expect_warning(expect_error(CalcICV(array(1:100, c(10, 10))),
                              "covariance matrix must be symmetric."), "deprecated")
  w <- capture_warnings(CalcICV(RandomMatrix(10)))
  expect_match(w, "Matrix appears to be a correlation matrix! Only covariance matrices should be used for ICV.",
               all = FALSE)
})
test_that("CalcICV returns sensible results",
{
  cov.matrix <- cov(iris[,1:4])
  eVals <- eigen(cov.matrix, only.values = TRUE)$values
  ident <- diag(rep(2, 10))
  suppressWarnings({
    expect_equal(CalcICV(ident), 0)
    expect_equal(CalcICV(cov.matrix), sd(eVals)/mean(eVals))
  })
})
test_that("CalcICV recognises correlation matrices with rounding error",
{
  cm <- cov2cor(cov(iris[,1:4])); diag(cm) <- diag(cm) + 1e-12
  expect_match(capture_warnings(CalcICV(cm)), "correlation matrix", all = FALSE)
})
