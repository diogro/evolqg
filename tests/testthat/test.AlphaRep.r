test_that("AlphaRep returns correct values",
          {
              cor.matrix <- cor(iris[,1:4])
              tam <- nrow(iris)
              vec <- cor.matrix[lower.tri(cor.matrix)]
              var.erro <- (1 - mean(vec)^2)/(tam - 2)
              var.vec <- var(vec)
              iris.rep <- (var.vec - var.erro)/var.vec
              expect_equal(AlphaRep(cor.matrix, tam), iris.rep)
              expect_equal(AlphaRep(cor.matrix, tam, method = "cheverud"), iris.rep)
              expect_error(AlphaRep(cov(iris[,1:4])), "Matrices do not appear to be correlation matrices.")
          }
)

test_that("AlphaRep accepts correlation matrices with rounding error on the diagonal",
          {
              cor.matrix <- cor(iris[,1:4])
              diag(cor.matrix) <- diag(cor.matrix) + 1e-12
              expect_no_error(AlphaRep(cor.matrix, 150))
          })

test_that("AlphaRep methods use the documented formulas",
          {
              cor.matrix <- cor(iris[,1:4])
              N <- 40
              r <- cor.matrix[lower.tri(cor.matrix)]
              alpha <- function(e) (var(r) - e)/var(r)
              expect_equal(AlphaRep(cor.matrix, N, "squared"), alpha((1 - mean(r^2))/(N - 2)))
              expect_equal(AlphaRep(cor.matrix, N, "exact"), alpha(mean((1 - r^2)^2)/(N - 1)))
              expect_error(AlphaRep(cor.matrix, N, "other"))
          })
