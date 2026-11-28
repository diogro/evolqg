test_that("CalculateMatrix returns correct result",
          {
            data(iris)
            iris.lm = lm(as.matrix(iris[,1:4])~iris[,5])
            cov.matrix = var(iris.lm$residuals)*((dim(iris.lm$residuals)[1]-1)/iris.lm$df.residual)
            expect_equal(CalculateMatrix(iris.lm), cov.matrix, tolerance = 1e-14)
            expect_equal(CalculateMatrix(lm(as.matrix(iris[,1:4])~1)), cov(iris[,1:4]), tolerance = 1e-14)
          }
)

test_that("CalculateMatrix does not depend on the factor contrasts",
          {
            lm.default <- lm(as.matrix(iris[,1:4])~iris[,5])
            old <- options(contrasts=c("contr.sum","contr.poly"))
            on.exit(options(old))
            lm.sum <- lm(as.matrix(iris[,1:4])~iris[,5])
            expect_equal(CalculateMatrix(lm.default), CalculateMatrix(lm.sum))
          }
)

test_that("CalculateMatrix does not re-center the residuals of models without intercept",
          {
            local_seed(1)
            x <- rnorm(30)
            y <- cbind(a = 2 * x + rnorm(30), b = -x + rnorm(30) + 1)
            m <- lm(y ~ 0 + x)
            res <- residuals(m)
            expect_false(isTRUE(all.equal(colMeans(res), c(a = 0, b = 0))))
            expect_equal(CalculateMatrix(m), crossprod(res) / m$df.residual)
            expect_equal(dimnames(CalculateMatrix(m)), list(c("a", "b"), c("a", "b")))
          }
)

test_that("CalculateMatrix handles a univariate model",
          {
            m <- lm(Sepal.Length ~ Species, iris)
            res <- CalculateMatrix(m)
            expect_true(is.matrix(res))
            expect_equal(dim(res), c(1, 1))
            expect_equal(drop(res), sum(residuals(m)^2) / m$df.residual)
          }
)

test_that("BayesianCalculateMatrix returns MAP, MLE and posterior samples",
          {
            m <- lm(as.matrix(iris[,1:4]) ~ Species, iris)
            res <- BayesianCalculateMatrix(m)
            expect_named(res, c("MAP", "MLE"))
            expect_equal(res$MLE, CalculateMatrix(m))
            local_seed(5)
            res <- BayesianCalculateMatrix(m, samples = 20, nu = 6)
            expect_named(res, c("MAP", "MLE", "P", "Ps"))
            expect_s3_class(res$Ps, "mcmc_sample")
            expect_equal(dim(res$Ps), c(20, 4, 4))
            expect_equal(res$P, apply(unclass(res$Ps), 2:3, median), ignore_attr = TRUE)
            expect_named(dimnames(res$P), c("X1", "X2"))
            expect_null(names(dimnames(res$MAP)))
            expect_output(print(res$Ps), "MCMC sample of 20 matrices of dimension 4x4")
            S_x <- crossprod(residuals(m))
            S_0 <- diag(diag(S_x / m$df.residual) * 6)
            expect_equal(unname(res$MAP), unname((S_0 + S_x) / (6 + m$df.residual)))
          }
)

test_that("BayesianCalculateMatrix keeps its random stream",
          {
            ## reference values computed with evolqg 0.4-3 (MCMCpack::riwish)
            m <- lm(as.matrix(iris[,1:4]) ~ Species, iris)
            local_seed(5)
            res <- BayesianCalculateMatrix(m, samples = 20, nu = 6)
            expect_equal(unname(res$P[1, ]), c(0.27614745, 0.10078031, 0.17146138, 0.04588523),
                         tolerance = 1e-7)
          }
)

test_that("BayesianCalculateMatrix handles a univariate model",
          {
            m <- lm(Sepal.Length ~ Species, iris)
            local_seed(1)
            res <- BayesianCalculateMatrix(m, samples = 10)
            expect_equal(dim(res$MAP), c(1, 1))
            expect_equal(dim(res$Ps), c(10, 1, 1))
            expect_true(all(res$Ps > 0))
          }
)
