check_rarefaction <- function(results, n = 30, reps = 5) {
  expect_type(results, "list")
  expect_equal(names(results), as.character(3:n))
  expect_length(results[[1]], reps)
}
test_that("RarefactionRandomSkewers returns sensible results",
          {
            local_seed(1)
            ind.data <- matrix(rnorm(30*10), 30, 10)
            check_rarefaction(Rarefaction(ind.data, RandomSkewers, num.reps = 5))
          }
)
test_that("RarefactionMantelCor returns sensible results",
          {
            local_seed(2)
            ind.data <- matrix(rnorm(30*10), 30, 10)
            check_rarefaction(Rarefaction(ind.data, MatrixCor, correlation = TRUE, num.reps = 5))
          }
)
test_that("RarefactionKrzCor returns sensible results",
          {
            local_seed(3)
            ind.data <- matrix(rnorm(30*10), 30, 10)
            check_rarefaction(Rarefaction(ind.data, KrzCor, num.reps = 5))
            check_rarefaction(Rarefaction(ind.data, KrzCor, num.reps = 5, correlation = TRUE))
          }
)
test_that("RarefactionPCAsimilarity returns sensible results",
{
  local_seed(4)
  ind.data <- matrix(rnorm(30*10), 30, 10)
  check_rarefaction(Rarefaction(ind.data, PCAsimilarity, num.reps = 5, correlation = TRUE))
}
)
test_that("Rarefaction passes replace to RarefactionStat",
          {
            local_seed(5)
            ind.data <- matrix(rnorm(12*4), 12, 4)
            seen <- list()
            f <- function(d) { seen[[length(seen) + 1]] <<- nrow(unique(d)); cov(d) }
            # with replacement, a sample as large as the data has repeated rows
            RarefactionStat(ind.data, f, function(x, y) 1, num.reps = 20, replace = TRUE)
            n.unique <- unlist(seen)[-1]
            expect_true(any(n.unique < 12))
            seen <- list()
            RarefactionStat(ind.data, f, function(x, y) 1, num.reps = 5)
            expect_true(all(unlist(seen)[-1] == rep(3:12, each = 5)))
            expect_no_error(Rarefaction(ind.data, KrzCor, num.reps = 2, replace = TRUE))
          })
test_that("Rarefaction is deterministic given a seed and equals the plain loop", {
  ind.data <- matrix(rnorm(15*4), 15, 4)
  local_seed(6)
  a <- RarefactionStat(ind.data, cov, function(x, y) KrzCor(x, y)[1], num.reps = 3)
  local_seed(6)
  b <- RarefactionStat(ind.data, cov, function(x, y) KrzCor(x, y)[1], num.reps = 3)
  expect_identical(a, b)
})
test_that("RarefactionStat passes ... to ComparisonFunc",
          {
            ind.data <- matrix(rnorm(8*3), 8, 3)
            out <- RarefactionStat(ind.data, cov, function(x, y, k) k, k = 7, num.reps = 2)
            expect_true(all(unlist(out) == 7))
          })
test_that("RarefactionStat does not hang when StatFunc always warns",
          {
            setTimeLimit(elapsed = 60, transient = TRUE)
            withr::defer(setTimeLimit())
            ind.data <- matrix(rnorm(6*3), 6, 3)
            # warns for every sample, also on the observed data
            always <- function(x) { warning("always"); cov(x) }
            expect_warning(out <- RarefactionStat(ind.data, always, function(x, y) 1, num.reps = 2),
                           "observed data")
            expect_length(out, 4)
            # warns only for small samples: kept after the retries run out
            small <- function(x) { if (nrow(x) < 5) warning("small"); cov(x) }
            expect_warning(out <- RarefactionStat(ind.data, small, function(x, y) 1, num.reps = 2),
                           "kept")
            expect_length(out, 4)
            # constant column with correlation: the observed statistic is NA
            ind.data[, 1] <- 1
            expect_error(Rarefaction(ind.data, KrzCor, correlation = TRUE, num.reps = 2), "missing values")
          })
test_that("RarefactionStat does not warn about deprecated functions with CalcR2-free code",
          expect_no_warning(RarefactionStat(matrix(rnorm(30), 10, 3), cor,
                                            function(x, y) mean(y[lower.tri(y)]^2), num.reps = 2)))
test_that("RarefactionStat needs at least 3 individuals",
          expect_error(RarefactionStat(matrix(rnorm(4), 2, 2), cov, function(x, y) 1), "at least 3"))
test_that("RarefactionStat works on data frames",
          {
            local_seed(7)
            a <- RarefactionStat(iris[1:20, 1:4], cov, function(x, y) KrzCor(x, y)[1], num.reps = 2)
            local_seed(7)
            b <- RarefactionStat(as.matrix(iris[1:20, 1:4]), cov, function(x, y) KrzCor(x, y)[1], num.reps = 2)
            expect_identical(a, b)
          })
test_that("Rarefaction in parallel equals serial",
          {
            skip_if_not_installed("doParallel")
            skip_on_os("windows")
            doParallel::registerDoParallel(cores = 2)
            withr::defer(doParallel::stopImplicitCluster())
            ind.data <- matrix(rnorm(20*4), 20, 4)
            local_seed(8)
            serial <- Rarefaction(ind.data, KrzCor, num.reps = 3)
            local_seed(8)
            par <- Rarefaction(ind.data, KrzCor, num.reps = 3, parallel = TRUE)
            expect_identical(serial, par)
          })
test_that("PlotRarefaction returns a ggplot with the expected data",
          {
            local_seed(9)
            res <- Rarefaction(matrix(rnorm(15*4), 15, 4), KrzCor, num.reps = 3)
            p <- PlotRarefaction(res, "Krz")
            expect_s3_class(p, "ggplot")
            expect_named(p$data, c("avg.corr", "L1"))
            expect_equal(nrow(p$data), 3 * 13)
            expect_equal(sort(unique(p$data$L1)), 3:15)
            expect_equal(p$data$avg.corr, unlist(res, use.names = FALSE))
            expect_no_error(ggplot2::ggplot_build(p))
            expect_no_warning(ggplot2::ggplot_build(p))
          })
