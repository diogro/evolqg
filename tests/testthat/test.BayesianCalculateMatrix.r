# Tests for BayesianCalculateMatrix live in test.CalculateMatrix.r
test_that("print.mcmc_sample prints the median matrix",
          {
            local_seed(1)
            m <- lm(as.matrix(iris[,1:3]) ~ Species, iris)
            res <- BayesianCalculateMatrix(m, samples = 5)
            expect_output(print(res$Ps), "Median matrix")
          }
)
