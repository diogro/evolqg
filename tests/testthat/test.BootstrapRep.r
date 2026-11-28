test_that("BootstrapRep returns sensible results",
          {
              local_seed(1)
              ind.data <- iris[,1:4]
              expect_lte(BootstrapRep(ind.data, RandomSkewers, iterations = 50), 1)
              expect_lte(BootstrapRep(ind.data, MatrixCor, iterations = 50, correlation = TRUE), 1)
              expect_lte(BootstrapRep(ind.data, KrzCor, iterations = 50), 1)
              expect_lte(BootstrapRep(ind.data, KrzCor, iterations = 50, correlation = TRUE), 1)
              expect_gte(BootstrapRep(ind.data, RandomSkewers, iterations = 50), 0.98)
              expect_gte(BootstrapRep(ind.data, RandomSkewers, iterations = 50, correlation = TRUE), 0.98)
              expect_gte(BootstrapRep(ind.data, MatrixCor, iterations = 50, correlation = TRUE), 0.98)
              expect_gte(BootstrapRep(ind.data, KrzCor, iterations = 50), 0.98)
              expect_gte(BootstrapRep(ind.data, KrzCor, iterations = 50, correlation = TRUE), 0.98)
              expect_error(BootstrapRep(cov(ind.data), RandomSkewers, iterations = 50), "input appears to be a matrix, use residuals.")
              expect_error(BootstrapRep(cor(ind.data), RandomSkewers, iterations = 50), "input appears to be a matrix, use residuals.")
          }
)
test_that("BootstrapRep passes ... to ComparisonFunc",
          {
              local_seed(2)
              seen <- NULL
              f <- function(x, y, extra) { seen <<- extra; 1 }
              expect_equal(BootstrapRep(iris[,1:4], f, iterations = 3, extra = "hello"), 1)
              expect_equal(seen, "hello")
          })
test_that("BootstrapR2 returns sensible results",
{
  local_seed(3)
  ind.data <- iris[,1:4]
  expect_length(BootstrapR2(ind.data), 1000)
  expect_length(BootstrapR2(ind.data, 10), 10)
  results <- BootstrapR2(ind.data, 10)
  expect_true(all(results < 1 & results > 0))
  expect_error(BootstrapR2(RandomMatrix(10), 10), "input appears to be a matrix, use residuals.")
})
test_that("BootstrapR2 does not warn about deprecated functions",
          expect_no_warning(BootstrapR2(iris[,1:4], 5)))

test_that("BootstrapStat returns a data frame with the documented columns",
          {
              local_seed(4)
              out <- BootstrapStat(iris[,1:4], 6, function(x, y) c(a = y[1, 1], b = 2), cov)
              expect_s3_class(out, "data.frame")
              expect_named(out, c("X1", "a", "b"))
              expect_s3_class(out$X1, "factor")
              expect_equal(levels(out$X1), as.character(1:6))
              out <- BootstrapStat(iris[,1:4], 3, function(x, y) y[1, 1], cov)
              expect_named(out, c("X1", "V1"))
              expect_equal(nrow(BootstrapStat(iris[,1:4], 0, function(x, y) 1, cov)), 0)
              expect_equal(nrow(BootstrapStat(iris[,1:4], 1, function(x, y) 1, cov)), 1)
          })
test_that("BootstrapStat gives the same results for data frames and matrices",
          {
              local_seed(5)
              a <- BootstrapStat(iris[,1:4], 10, function(x, y) y[1, 2], cov)
              local_seed(5)
              b <- BootstrapStat(as.matrix(iris[,1:4]), 10, function(x, y) y[1, 2], cov)
              expect_identical(a, b)
          })
test_that("BootstrapStat passes data frames to custom StatFunc",
          {
              local_seed(6)
              classes <- character()
              f <- function(d) { classes <<- c(classes, class(d)[1]); colMeans(d) }
              BootstrapStat(iris[,1:4], 3, function(x, y) y[1], f)
              expect_true(all(classes == "data.frame"))
          })
test_that("BootstrapStat fails instead of hanging on NA or errors",
          {
              setTimeLimit(elapsed = 60, transient = TRUE)
              withr::defer(setTimeLimit())
              d <- iris[,1:4]; d[3, 2] <- NA
              expect_error(BootstrapStat(d, 3, function(x, y) 1, cov), "missing values")
              expect_error(BootstrapRep(d, RandomSkewers, iterations = 3), "missing values")
              expect_error(BootstrapStat(iris[,1:4], 2, function(x, y) stop("boom"), cov), "boom")
              expect_error(BootstrapStat(iris[,1:4], 2, function(x, y) NA_real_, cov), "NA")
              expect_error(BootstrapRep(iris[,1:4], function(x, y) 1, iterations = 2, nonexistent.arg = 1), "unused argument")
          })
test_that("BootstrapStat redraws resamples that fail and warns once",
          {
              local_seed(7)
              calls <- 0
              f <- function(x, y) { calls <<- calls + 1; if (calls == 2) stop("transient") else 1 }
              expect_warning(out <- BootstrapStat(iris[,1:4], 4, f, cov), "redrawn")
              expect_equal(nrow(out), 4)
              expect_false(anyNA(out$V1))
          })
test_that("BootstrapStat in parallel equals serial",
          {
              skip_if_not_installed("doParallel")
              skip_on_os("windows")
              doParallel::registerDoParallel(cores = 2)
              withr::defer(doParallel::stopImplicitCluster())
              f <- function(x, y) PCAsimilarity(x, y)[1]
              local_seed(8)
              serial <- BootstrapStat(iris[,1:4], 20, f, cov)
              local_seed(8)
              par <- BootstrapStat(iris[,1:4], 20, f, cov, parallel = TRUE)
              expect_identical(serial, par)
          })
