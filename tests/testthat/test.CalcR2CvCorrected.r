test_that("CalcR2CvCorrected returns sensible results on residuals",
          {
              local_seed(1)
              results <- CalcR2CvCorrected(iris[,1:4], iterations = 100)
              expect_length(results, 3)
              expect_type(results, "list")
              expect_type(results[[1]], "double")
              expect_type(results[[2]], "list")
              expect_s3_class(results[[3]], "data.frame")
              expect_equal(names(results), c("adjusted.integration.index", "models", "dist"))
              expect_equal(colnames(results[[3]]), c("r2", "eVals_cv", "mean_cv"))
              expect_true(all(results[[3]][,1] > 0 & results[[3]][,1] < 1))
              expect_length(results[[1]], 2)
              expect_length(results[[2]], 2)
              expect_equal(dim(results[[3]]), c(100, 3))
          })

test_that("CalcR2CvCorrected returns sensible results on models",
          {
              local_seed(2)
              iris.lm <- lm(as.matrix(iris[,1:4]) ~ iris[,5])
              results <- CalcR2CvCorrected(iris.lm, iterations = 100)
              expect_length(results, 3)
              expect_type(results, "list")
              expect_type(results[[1]], "double")
              expect_type(results[[2]], "list")
              expect_true(is.matrix(results[[3]]))
              expect_equal(names(results), c("adjusted.integration.index", "models", "dist"))
              expect_equal(colnames(results[[3]]), c("r2", "eVals_cv", "mean_cv"))
              expect_true(all(results[[3]][,1] > 0 & results[[3]][,1] < 1))
              expect_equal(dim(results[[3]]), c(100, 3))
          })

test_that("CalcR2CvCorrected.lm works with more than one predictor",
          {
              local_seed(3)
              d <- iris[,1:4]
              d$sex <- factor(rep(c("a", "b"), 75))
              d$sp <- iris$Species
              m <- lm(as.matrix(d[,1:4]) ~ sp + sex, data = d)
              results <- CalcR2CvCorrected(m, iterations = 30)
              expect_equal(dim(results$dist), c(30, 3))
              expect_false(anyNA(results$dist))
              # the groups are the combinations of the factors
              m1 <- lm(as.matrix(d[,1:4]) ~ sp:sex, data = d)
              expect_equal(dim(CalcR2CvCorrected(m1, iterations = 5)$dist), c(5, 3))
          })

test_that("CalcR2CvCorrected.lm with one factor matches a direct computation",
          {
              local_seed(4)
              m <- lm(as.matrix(iris[,1:4]) ~ iris[,5])
              out <- CalcR2CvCorrected(m, iterations = 3)
              local_seed(4)
              cv <- function(x) sd(x)/mean(x)
              resp <- as.matrix(iris[,1:4]); n <- 150
              s <- sample(n, replace = TRUE)
              by.group <- lapply(split(seq_len(n), iris$Species[s]), function(i) apply(resp[s,][i,], 2, cv))
              ns <- table(iris$Species[s])
              mcv <- sum(rowMeans(do.call(rbind, by.group)) * ns / sum(ns))
              expect_equal(unname(out$dist[1, "mean_cv"]), mcv)
          })

test_that("CalcR2CvCorrected.lm does not hang when a group is too small",
          {
              setTimeLimit(elapsed = 60, transient = TRUE)
              withr::defer(setTimeLimit())
              # a constant trait makes the mean cv undefined for every sample
              d <- as.matrix(iris[,1:4]); d[, 1] <- 0
              expect_error(CalcR2CvCorrected(lm(d ~ iris$Species), iterations = 2), "at least two individuals")
          })

test_that("CalcR2CvCorrected.lm in parallel equals serial",
          {
              skip_if_not_installed("doParallel")
              skip_on_os("windows")
              doParallel::registerDoParallel(cores = 2)
              withr::defer(doParallel::stopImplicitCluster())
              m <- lm(as.matrix(iris[,1:4]) ~ iris[,5])
              local_seed(5)
              serial <- CalcR2CvCorrected(m, iterations = 20)
              local_seed(5)
              par <- CalcR2CvCorrected(m, iterations = 20, parallel = TRUE)
              expect_identical(serial$dist, par$dist)
          })
test_that("CalcR2CvCorrected does not warn about deprecated functions",
          expect_no_warning(CalcR2CvCorrected(iris[,1:4], iterations = 5)))
