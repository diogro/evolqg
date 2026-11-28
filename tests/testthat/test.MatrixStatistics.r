test_that("Autonomy returns correct results",
          {
              local_seed(1)
              data(iris)
              cov.matrix <- cov(iris[,1:4])
              beta <- Normalize(rnorm(4))
              auto <- ((t (beta) %*% solve (cov.matrix) %*% beta)^(-1)) / (t (beta) %*% cov.matrix %*% beta)
              expect_equal(Autonomy(cov.matrix, beta), as.numeric(auto))
              local_seed(2)
              rmatrix <- RandomMatrix(40, 1, 1, 10)
              ev <- eigen(rmatrix)
              ev$values[35:40] <- 0
              singular.matrix <- ev$vectors %*% diag(ev$values) %*% t(ev$vectors)
              expect_warning(Autonomy(singular.matrix), 
                             "matrix is singular, can't compute autonomy directly. Using nearPD, results could be wrong")
          }
)
test_that("ConditionalEvolvability returns correct restults",
          {
              local_seed(1)
              data(iris)
              cov.matrix <- cov(iris[,1:4])
              beta <- Normalize(rnorm(4))
              cond.evol <- (t (beta) %*% solve (cov.matrix) %*% beta)^(-1)
              expect_equal(ConditionalEvolvability(cov.matrix, beta), as.numeric(cond.evol))
              local_seed(2)
              rmatrix <- RandomMatrix(40, 1, 1, 10)
              ev <- eigen(rmatrix)
              ev$values[35:40] <- 0
              singular.matrix <- ev$vectors %*% diag(ev$values) %*% t(ev$vectors)
              expect_warning(ConditionalEvolvability(singular.matrix), 
                             "matrix is singular, can't compute conditional evolvability directly. Using nearPD, results could be wrong")
          }
)
test_that("Constraints returns correct results",
          {
              data(iris)
              cov.matrix <- cov(iris[,1:4])
              beta = Normalize(rnorm(4))
              const = abs (t (Normalize (eigen (cov.matrix)$vectors[,1])) %*% Normalize (cov.matrix %*% beta))
              expect_equal(Constraints(cov.matrix, beta), as.numeric(const))
          }
)
test_that("Evolvability returns correct results",
          {
              local_seed(1)
              data(iris)
              cov.matrix <- cov(iris[,1:4])
              beta <- Normalize(rnorm(4))
              evol <- t (beta) %*% cov.matrix %*% beta
              expect_equal(Evolvability(cov.matrix, beta), as.numeric(evol))
          }
)
test_that("Flexibility returns correct result",
         {
           local_seed(1)
           data(iris)
           beta <- Normalize(rnorm(4))
           cov.matrix <- cov(iris[,1:4])
           flex <- t (beta) %*% cov.matrix %*% beta / Norm (cov.matrix %*% beta)
           expect_equal(Flexibility(cov.matrix, beta), as.numeric(flex))
           expect_true(Flexibility(cov.matrix, beta) <=  1)
           expect_true(Flexibility(cov.matrix, beta) >= -1)
         }
)

test_that("Pc1Percent returns correct results",
          {
            local_seed(1)
            cov.matrix = cov(matrix(rnorm(30*10), 30, 10))
            pc1 <- eigen(cov.matrix)$values[1]/sum(diag(cov.matrix))
            expect_equal(Pc1Percent(cov.matrix), pc1)
            expect_true(Pc1Percent(cov.matrix) <=  1)
            expect_true(Pc1Percent(cov.matrix) > 0)
          }
)
test_that("Respondability returns correct result",
         {
           local_seed(1)
           data(iris)
           beta <- Normalize(rnorm(4))
           cov.matrix <- cov(iris[,1:4])
           repond <- Norm(cov.matrix %*% beta)
           expect_equal(Respondability(cov.matrix, beta), repond)
         }
)

test_that("MeanMatrixStatistics excludes legacy integration metrics by default",
          {
            stats <- MeanMatrixStatistics(cov(iris[,1:4]))
            expect_false("MeanSquaredCorrelation" %in% names(stats))
            expect_false("ICV" %in% names(stats))
            expect_true("pc1.percent" %in% names(stats))
            expect_true("EigenSd" %in% names(stats))
          }
)

test_that("MeanMatrixStatistics can include R2 and ICV on demand",
          {
            cov.matrix <- cov(iris[,1:4])
            stats <- suppressWarnings(MeanMatrixStatistics(cov.matrix, include.R2 = TRUE, include.ICV = TRUE))
            expect_equal(as.numeric(stats["MeanSquaredCorrelation"]), mean(cov2cor(cov.matrix)[lower.tri(cov2cor(cov.matrix))]^2))
            expect_equal(as.numeric(stats["ICV"]), sd(eigen(cov.matrix, only.values = TRUE)$values) / mean(eigen(cov.matrix, only.values = TRUE)$values))
          }
)

singular_matrix <- function() {
  local_seed(2, .local_envir = parent.frame())
  ev <- eigen(RandomMatrix(40, 1, 1, 10))
  ev$values[35:40] <- 0
  ev$vectors %*% diag(ev$values) %*% t(ev$vectors)
}

test_that("Autonomy and ConditionalEvolvability use nearPD when the matrix is singular",
          {
            S <- singular_matrix()
            beta <- apply(matrix(rnorm(40 * 25), 40), 2, Normalize)
            nearpd <- as.matrix(Matrix::nearPD(S)$mat)
            cond.direct <- 1 / colSums(beta * solve(nearpd, beta))
            auto.direct <- cond.direct / colSums(beta * (nearpd %*% beta))
            expect_warning(cond <- ConditionalEvolvability(S, beta), "singular")
            expect_warning(auto <- Autonomy(S, beta), "singular")
            expect_true(all(cond > 0))
            expect_true(all(auto > 0))
            expect_equal(cond, cond.direct)
            expect_equal(auto, auto.direct)
          }
)

test_that("Autonomy and ConditionalEvolvability do not fail on exactly singular matrices",
          {
            local_seed(3)
            X <- matrix(rnorm(5 * 8), 8, 5)
            C <- tcrossprod(X)  # rank 5, 8 x 8
            expect_warning(cond <- ConditionalEvolvability(C, iterations = 20), "singular")
            expect_true(all(is.finite(cond)) && all(cond > 0))
          }
)

test_that("MeanMatrixStatistics works with a user supplied beta.mat",
          {
            local_seed(4)
            beta <- apply(matrix(rnorm(4 * 50), 4), 2, Normalize)
            expect_output(res <- MeanMatrixStatistics(cov(iris[,1:4]), beta.mat = beta, full.results = TRUE),
                          "critical value:  [0-9.]+")
            expect_equal(nrow(res$dist), 50)
            expect_false(anyNA(res$mean))
          }
)

test_that("MeanMatrixStatistics works for one trait and through foreach",
          {
            local_seed(5)
            one <- MeanMatrixStatistics(matrix(2, 1, 1), iterations = 20)
            expect_equal(one[["evolvability"]], 2)
            foreach::registerDoSEQ()
            local_seed(6)
            seq.res <- MeanMatrixStatistics(cov(iris[,1:4]), iterations = 50, full.results = TRUE)
            local_seed(6)
            par.res <- suppressWarnings(MeanMatrixStatistics(cov(iris[,1:4]), iterations = 50,
                                                             full.results = TRUE, parallel = TRUE))
            expect_equal(par.res, seq.res)
          }
)

test_that("MeanMatrixStatistics keeps its output structure",
          {
            local_seed(7)
            res <- MeanMatrixStatistics(cov(iris[,1:4]), iterations = 30, full.results = TRUE)
            expect_named(res, c("dist", "mean"))
            expect_equal(dim(res$dist), c(30, 7))
            expect_equal(colnames(res$dist)[7], "null.dist")
            expect_named(res$mean, c("pc1.percent", "EigenSd", "respondability", "evolvability",
                                     "conditional.evolvability", "autonomy", "flexibility", "constraints"))
          }
)
