test_that("CalcEigenVar throws error",
          {
            expect_error(CalcEigenVar(array(1:100, c(10, 10))),
                         "covariance matrix must be symmetric.")
          })
test_that("CalcEigenVar returns sensible results",
          {
            cov.matrix <- cov(iris[,1:4])
            eVals <- eigen(cov.matrix, only.values = TRUE)$values
            ident <- diag(rep(2, 10))
            expect_equal(CalcEigenVar(ident), 0)
            expect_equal(CalcEigenVar(cov.matrix, sd = FALSE, rel = FALSE),
                         sum((eVals - mean(eVals))^2)/length(eVals))
            expect_equal(CalcEigenVar(cov.matrix, sd = TRUE, rel = FALSE),
                         sqrt(sum((eVals - mean(eVals))^2)/length(eVals)))
          }
)
test_that("CalcEigenVar returns NA with a warning when sd = TRUE and the variance is below its expectation",
          {
            ident <- diag(4)
            expect_warning(out <- CalcEigenVar(ident, sd = TRUE, sample = 10), "Returning NA")
            expect_true(is.na(out))
            expect_warning(out <- CalcEigenVar(ident, sd = TRUE, rel = FALSE, sample = 10), "Returning NA")
            expect_true(is.na(out))
            # no warning when the corrected variance is positive
            expect_no_warning(CalcEigenVar(cov(iris[,1:4]), sd = TRUE, sample = 150))
            expect_no_warning(CalcEigenVar(cov(iris[,1:4]), sample = 10))
          })
test_that("CalcEigenVar with keep.positive uses the number of positive eigenvalues",
          {
            m <- diag(c(3, 1, 0, 0))
            e <- c(3, 1)
            expect_equal(CalcEigenVar(m, rel = FALSE), mean((e - mean(e))^2))
            # maximum for 2 eigenvalues with trace 4: (2 - 1) * 4^2 / 2^2
            expect_equal(CalcEigenVar(m), mean((e - mean(e))^2)/4)
          })
