test_that("CalcR2 returns correct value for correlation matrices",
          {
            cor.matrix <- cor(iris[,1:4])
            cov.matrix <- cov(iris[,1:4])
            ident <- diag(rep(1, 10))
            expected <- mean((cor.matrix[lower.tri(cor.matrix)])^2)
            expect_warning(r.ident <- CalcR2(ident), "deprecated")
            expect_equal(r.ident, 0)
            expect_warning(r.cor <- CalcR2(cor.matrix), "deprecated")
            expect_equal(r.cor, expected)
            expect_warning(r.cov <- CalcR2(cov.matrix), "deprecated")
            expect_equal(r.cov, expected)
          }
)
