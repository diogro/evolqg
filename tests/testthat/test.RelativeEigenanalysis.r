test_that("RelativeEigenanalysis returns real eigenvalues that solve the generalised problem",
          {
            local_seed(1)
            for (i in 1:10) {
              x <- cov(matrix(rnorm(35 * 30), 35)); y <- cov(matrix(rnorm(35 * 30), 35))
              res <- RelativeEigenanalysis(x, y)
              expect_type(res$values, "double")
              expect_false(is.complex(res$vectors))
              expect_equal(x %*% solve(y) %*% res$vectors, res$vectors %*% diag(res$values),
                           tolerance = 1e-6)
              expect_equal(sqrt(colSums(res$vectors^2)), rep(1, 30))
            }
          }
)

test_that("RelativeEigenanalysis agrees with eigen of the product and the symmetric version",
          {
            x <- cov(iris[,1:4]); y <- cov(iris[iris$Species == "setosa", 1:4])
            res <- RelativeEigenanalysis(x, y)
            ref <- eigen(x %*% solve(y))
            expect_equal(res$values, ref$values)
            expect_equal(abs(res$vectors), abs(ref$vectors))
            sym <- RelativeEigenanalysis(x, y, symmetric = TRUE)
            expect_equal(sym$values, ref$values)
            expect_s3_class(sym, "eigen")
          }
)
