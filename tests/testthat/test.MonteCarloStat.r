test_that("MonteCarloR2 returns sensible results",
          {
            local_seed(1)
            expect_length(MonteCarloR2(diag(rep(2, 10)), 10), 1000)
            expect_length(MonteCarloR2(diag(rep(2, 10)), 10, 10), 10)
            results <- MonteCarloR2(diag(rep(2, 10)), 10)
            expect_true(all(results < 1 & results > 0))
            expect_warning(MonteCarloR2(RandomMatrix(10), 10),
                           "Matrix appears to be a correlation matrix! Only covariance matrices should be used in parametric resampling.")
          }
)
test_that("MonteCarloR2 does not warn about deprecated functions",
          expect_no_warning(MonteCarloR2(diag(rep(2, 5)), 10, 5)))

test_that("MonteCarloRep returns sensible results",
          {
            local_seed(2)
            expect_warning(MonteCarloRep(RandomMatrix(10), 30, RandomSkewers, iterations = 1),
                           "Matrix appears to be a correlation matrix! Only covariance matrices should be used in parametric resampling.")
            cov.matrix <- RandomMatrix(10, 1, 1, 10)
            expect_lte(MonteCarloRep(cov.matrix, sample.size = 30, RandomSkewers, iterations = 50), 1)
            expect_lte(MonteCarloRep(cov.matrix, 30, MatrixCor, correlation = TRUE, iterations = 50), 1)
            expect_lte(MonteCarloRep(cov.matrix, 30, KrzCor, iterations = 50), 1)
            expect_lte(MonteCarloRep(cov.matrix, 30, KrzCor, correlation = TRUE, iterations = 50), 1)
          }
)
test_that("MonteCarloStat throws error",
{
  expect_error(MonteCarloStat(array(1:100, c(10, 10)), 10, 10, RandomSkewers, cov),
               "covariance matrix must be symmetric.")
})
test_that("MonteCarloStat output has the documented shape",
          {
            local_seed(3)
            cm <- RandomMatrix(5, 1, 1, 10)
            out <- MonteCarloStat(cm, 20, 7, function(x, y) c(a = y[1, 1], b = 1), cov)
            expect_s3_class(out, "data.frame")
            expect_named(out, c("X1", "a", "b"))
            expect_equal(levels(out$X1), as.character(1:7))
            out <- MonteCarloStat(cm, 20, 7, function(x, y) y[1, 1], cov)
            expect_named(out, c("X1", "V1"))
          })
test_that("MonteCarloStat reproduces the mvtnorm stream",
          {
            skip_if_not_installed("mvtnorm")
            cm <- RandomMatrix(5, 1, 1, 10)
            local_seed(4)
            out <- MonteCarloStat(cm, 20, 5, function(x, y) y, function(x) x[1, 1])
            local_seed(4)
            ref <- vapply(1:5, function(i) mvtnorm::rmvnorm(20, sigma = cm, method = "svd")[1, 1], 1)
            expect_equal(out$V1, ref)
          })
test_that("MonteCarloStat fails instead of hanging on NA or errors",
          {
            setTimeLimit(elapsed = 60, transient = TRUE)
            withr::defer(setTimeLimit())
            cm <- RandomMatrix(5, 1, 1, 10)
            bad <- cm; bad[1, 2] <- bad[2, 1] <- NA
            expect_error(MonteCarloStat(bad, 20, 3, function(x, y) 1, cov), "missing values")
            expect_error(MonteCarloStat(cm, 20, 2, function(x, y) stop("boom"), cov), "boom")
            expect_error(MonteCarloStat(cm, 20, 2, function(x, y) NA_real_, cov), "NA")
            expect_error(MonteCarloRep(cm, 20, function(x, y) 1, iterations = 2, nonexistent.arg = 1), "unused argument")
          })
test_that("MonteCarloStat redraws failed populations and warns once",
          {
            local_seed(5)
            calls <- 0
            f <- function(x, y) { calls <<- calls + 1; if (calls == 2) stop("transient") else 1 }
            expect_warning(out <- MonteCarloStat(RandomMatrix(5, 1, 1, 10), 20, 4, f, cov), "redrawn")
            expect_equal(nrow(out), 4)
          })
test_that("MonteCarloStat in parallel equals serial",
          {
            skip_if_not_installed("doParallel")
            skip_on_os("windows")
            doParallel::registerDoParallel(cores = 2)
            withr::defer(doParallel::stopImplicitCluster())
            cm <- RandomMatrix(5, 1, 1, 10)
            f <- function(x, y) PCAsimilarity(x, y)[1]
            local_seed(6)
            serial <- MonteCarloStat(cm, 20, 20, f, cov)
            local_seed(6)
            par <- MonteCarloStat(cm, 20, 20, f, cov, parallel = TRUE)
            expect_identical(serial, par)
          })
test_that("MonteCarloRep returns correct results for scalled crappy matrices",
          {
           skip_if_not_installed("mvtnorm")
           local_seed(7)
           mat <- RandomMatrix(35, 1, 1 ,10)
           crappy_mat <- var(mvtnorm::rmvnorm(10, sigma = mat))
           rep_c <- MonteCarloRep(crappy_mat, 10, PCAsimilarity, iterations = 100)
           rep_sc <- MonteCarloRep(100000*crappy_mat, 10, PCAsimilarity, iterations = 100)
           expect_lt(abs(rep_c - rep_sc), 0.05)
          }
)
