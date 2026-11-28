test_that("RandomSkewers returns correct results on pairs of matrices",
          {
            expect = c(correlation = 0.954184140494324, probability = 0.00561, correlation_sd = 0.0553934050116506)
            data(dentus)
            cor.matrix.1 <- cov(dentus[dentus$species=="A", 1:4])
            cor.matrix.2 <- cov(dentus[dentus$species=="B", 1:4])
            results <- RandomSkewers(cor.matrix.1, cor.matrix.2, num.vectors = 100000)
            expect_equal(results, expect, tolerance = 1e-2)
            expect_type(results, "double")
            expect_equal(length(results), 3)
            expect_true(results[1] <=  1)
            expect_true(results[1] >= -1)
            expect_true(results[2] <=  1)
            expect_true(results[2] >=  0)
            expect_true(results[3] >=  0)
          }
)

test_that("RandomSkewers returns correct results on lists",
          {
            mat.list <- lapply(as.list(1:10), function(x) RandomMatrix(10))
            rep.vec <- runif(length(mat.list), 0.8, 0.9)
            local_seed(42, "3.5.0")
            results.list <- RandomSkewers(mat.list)
            expect_type(results.list, "list")
            results <- results.list[[1]]
            expect_equal(sum(is.na(results)), 0)
            probabilities <- results.list[[2]]
            expect_equal(sum(is.na(probabilities)), 0)
            expect_type(results.list, "list")
            local_seed(42, "3.5.0")
            results.list.2 <- RandomSkewers(mat.list, repeat.vector = rep.vec)
            results.2 <- results.list.2[[1]]
            expect_equal(sum(is.na(results.2)), 0)
            probabilities.2 <- results.list.2[[2]]
            expect_equal(dim(results), c(length(mat.list),length(mat.list)))
            lower  <- results[lower.tri(results)]
            lower.bool = sapply(lower, function(x) isTRUE(x > -1 & x < 1))
            expect_equal(sum(lower.bool), length(lower.bool))
            upper  <- results[upper.tri(results, diag = TRUE)]
            upper.bool = sapply(upper, function(x) isTRUE(x == 0))
            expect_equal(sum(upper.bool), length(upper.bool))
            lower.2  <- results.2[lower.tri(results.2)]
            upper.2  <- t(results.2)[lower.tri(results.2)]
            expect_equal(lower.2, lower)
            expect_equal(diag(results.2), rep.vec, ignore_attr = TRUE)
            expect_equal(sum(upper.2 < lower.2), 0)
          }
)

test_that("RandomSkewers returns correct results on lists + matrices",
          {
            mat.list <- lapply(as.list(1:10), function(x) RandomMatrix(10))
            y.matrix <- RandomMatrix(10)
            local_seed(42, "3.5.0")
            results <- RandomSkewers(mat.list, y.matrix)
            expect_s3_class(results, "data.frame")
            expect_equal(sum(is.na(results)), 0)
            expect_equal(dim(results), c(length(mat.list), 3))
            names(mat.list) <- seq_along(mat.list)
            named.results <- RandomSkewers(mat.list, y.matrix)
            expect_s3_class(named.results, "data.frame")
            expect_equal(sum(is.na(named.results)), 0)
            expect_equal(dim(named.results), c(length(mat.list), 4))
            expect_equal(named.results[,".id"], names(mat.list))
          }
)
test_that("MantelCor returns correct results on non-landmark data",
          {
            cor.matrix.1 <- RandomMatrix(15)
            cor.matrix.2 <- RandomMatrix(15)
            results <- MantelCor(cor.matrix.1, cor.matrix.2)
            expect_equal(results[1], cor(cor.matrix.1[lower.tri(diag(15))],
                                              cor.matrix.2[lower.tri(diag(15))]), ignore_attr = TRUE)

          })

test_that("MantelCor returns correct results on landmark data",
          {
            cor.matrix.1 <- RandomMatrix(12)
            cor.matrix.2 <- RandomMatrix(12)
            results <- MantelCor(cor.matrix.1, cor.matrix.2, landmark.dim = 2)
            expect_equal(results[1], cor(cor.matrix.1[evolqg:::lower.tri.land(diag(12), 2)],
                                              cor.matrix.2[evolqg:::lower.tri.land(diag(12), 2)]), ignore_attr = TRUE)
            results <- MantelCor(cor.matrix.1, cor.matrix.2, landmark.dim = 3)
            expect_equal(results[1], cor(cor.matrix.1[evolqg:::lower.tri.land(diag(12), 3)],
                                              cor.matrix.2[evolqg:::lower.tri.land(diag(12), 3)]), ignore_attr = TRUE)
            expect_error(MantelCor(cor.matrix.1, cor.matrix.2, landmark.dim = 4))
          })

test_that("MantelCor returns correct results on list",
          {
            mat.list <- RandomMatrix(10, 10)
            rep.vec <- runif(length(mat.list), 0.8, 0.9)
            results.list <- MantelCor(mat.list)
            results <- results.list[[1]]
            expect_equal(results[2,1], MantelCor(mat.list[[1]],
                                                      mat.list[[2]])[1], ignore_attr = TRUE)
            expect_equal(results[8,5], MantelCor(mat.list[[8]],
                                                      mat.list[[5]])[1], ignore_attr = TRUE)
            expect_equal(sum(is.na(results)), 0)
            probabilities <- results.list[[2]]
            expect_equal(sum(is.na(probabilities)), 0)
            expect_type(results.list, "list")
            results.list.2 <- MantelCor(mat.list, repeat.vector = rep.vec)
            results.2 <- results.list.2[[1]]
            expect_equal(sum(is.na(results.2)), 0)
            probabilities.2 <- results.list.2[[2]]
            expect_equal(dim(results), c(length(mat.list),length(mat.list)))
            lower <- results[lower.tri(results)]
            lower.bool = sapply(lower, function(x) isTRUE(x > -1 & x < 1))
            expect_equal(sum(lower.bool), length(lower.bool))
            upper  <- results[upper.tri(results, diag = TRUE)]
            upper.bool = sapply(upper, function(x) isTRUE(x == 0))
            expect_equal(sum(upper.bool), length(upper.bool))
            lower.2 <- results.2[lower.tri(results.2)]
            upper.2 <- t(results.2)[lower.tri(results.2)]
            expect_equal(lower.2, lower)
            expect_equal(diag(results.2), rep.vec, ignore_attr = TRUE)
            expect_equal(sum(abs(upper.2) < abs(lower.2)), 0)
          }
)

test_that("MantelCor returns correct results on list for landmark data",
          {
            mat.list <- RandomMatrix(12, 10)
            rep.vec <- runif(length(mat.list), 0.8, 0.9)
            results.list <- MantelCor(mat.list, landmark.dim = 2)
            results <- results.list[[1]]
            expect_equal(results[2,1], MantelCor(mat.list[[1]],
                                                      mat.list[[2]], landmark.dim = 2)[1], ignore_attr = TRUE)
            expect_equal(results[8,5], MantelCor(mat.list[[8]],
                                                      mat.list[[5]], landmark.dim = 2)[1], ignore_attr = TRUE)
            results.list <- MantelCor(mat.list, landmark.dim = 3)
            results <- results.list[[1]]
            expect_equal(results[2,1], MantelCor(mat.list[[1]],
                                                      mat.list[[2]], landmark.dim = 3)[1], ignore_attr = TRUE)
            expect_equal(results[8,5], MantelCor(mat.list[[8]],
                                                      mat.list[[5]], landmark.dim = 3)[1], ignore_attr = TRUE)
            expect_equal(sum(is.na(results)), 0)
            probabilities <- results.list[[2]]
            expect_equal(sum(is.na(probabilities)), 0)
            expect_type(results.list, "list")
            results.list.2 <- MantelCor(mat.list, repeat.vector = rep.vec, landmark.dim = 3)
            results.2 <- results.list.2[[1]]
            expect_equal(sum(is.na(results.2)), 0)
            probabilities.2 <- results.list.2[[2]]
            expect_equal(dim(results), c(length(mat.list),length(mat.list)))
            lower <- results[lower.tri(results)]
            lower.bool = sapply(lower, function(x) isTRUE(x > -1 & x < 1))
            expect_equal(sum(lower.bool), length(lower.bool))
            upper  <- results[upper.tri(results, diag = TRUE)]
            upper.bool = sapply(upper, function(x) isTRUE(x == 0))
            expect_equal(sum(upper.bool), length(upper.bool))
            lower.2 <- results.2[lower.tri(results.2)]
            upper.2 <- t(results.2)[lower.tri(results.2)]
            expect_equal(lower.2, lower)
            expect_equal(diag(results.2), rep.vec, ignore_attr = TRUE)
            expect_equal(sum(abs(upper.2) < abs(lower.2)), 0)
          }
)

test_that("MantelCor returns correct results on lists + matrices",
          {
            mat.list <- lapply(as.list(1:10), function(x) RandomMatrix(11))
            y.matrix <- RandomMatrix(11)
            local_seed(42, "3.5.0")
            results <- MantelCor(mat.list, y.matrix)
            local_seed(42, "3.5.0")
            expect_equal(unlist(results[1,]), MantelCor(y.matrix, mat.list[[1]]), ignore_attr = TRUE)
            expect_s3_class(results, "data.frame")
            expect_equal(sum(is.na(results)), 0)
            expect_equal(dim(results), c(length(mat.list), 2))
            names(mat.list) <- seq_along(mat.list)
            named.results <- MantelCor(mat.list, y.matrix)
            expect_s3_class(named.results, "data.frame")
            expect_equal(sum(is.na(named.results)), 0)
            expect_equal(dim(named.results), c(length(mat.list), 3))
            expect_equal(named.results[,".id"], names(mat.list))
          }
)

test_that("MantelCor returns correct results on lists + matrices for landmark data",
          {
            mat.list <- lapply(as.list(1:10), function(x) RandomMatrix(12))
            y.matrix <- RandomMatrix(12)
            local_seed(42, "3.5.0")
            results <- MantelCor(mat.list, y.matrix, landmark.dim = 2)
            local_seed(42, "3.5.0")
            expect_equal(unlist(results[1,]), MantelCor(y.matrix, mat.list[[1]], landmark.dim = 2), ignore_attr = TRUE)
            expect_s3_class(results, "data.frame")
            expect_equal(sum(is.na(results)), 0)
            expect_equal(dim(results), c(length(mat.list), 2))
            names(mat.list) <- seq_along(mat.list)
            named.results <- MantelCor(mat.list, y.matrix)
            expect_s3_class(named.results, "data.frame")
            expect_equal(sum(is.na(named.results)), 0)
            expect_equal(dim(named.results), c(length(mat.list), 3))
            expect_equal(named.results[,".id"], names(mat.list))
          }
)

test_that("KrzCor returns correct results",
          {
              cov.matrix.1 <- RandomMatrix(10)
              cov.matrix.2 <- RandomMatrix(10)
              ret.dim = round(dim(cov.matrix.1)[1]/2 - 1)
              EigenVectors <- function (x) return (eigen(x)$vectors[,1:ret.dim])
              A <- EigenVectors (cov.matrix.1)
              B <- EigenVectors (cov.matrix.2)
              S <- t(A) %*% B %*% t(B) %*% A
              SL <- sum (eigen(S)$values) / ret.dim
              expect_equal(KrzCor(cov.matrix.1, cov.matrix.2), SL, ignore_attr = TRUE)
              expect_equal(KrzCor(cov.matrix.1, cov.matrix.2, 10), 1, ignore_attr = TRUE)
              expect_equal(KrzCor(x  <- RandomMatrix(11), x), 1, ignore_attr = TRUE)
              expect_true(KrzCor(cov.matrix.1, cov.matrix.2) <= 1)
              expect_true(KrzCor(cov.matrix.1, cov.matrix.2) > 0)
          }
)

test_that("KrzCor returns correct results",
          {
            mat.list <- lapply(as.list(1:10), function(x) RandomMatrix(10))
            rep.vec <- runif(length(mat.list), 0.8, 0.9)
            results <- KrzCor(mat.list)
            expect_equal(results[2,1], KrzCor(mat.list[[1]], mat.list[[2]]), ignore_attr = TRUE)
            results.2 <- KrzCor(mat.list, repeat.vector = rep.vec)
            expect_equal(sum(is.na(results)), 0)
            expect_equal(sum(is.na(results.2)), 0)
            expect_equal(dim(results), c(length(mat.list),length(mat.list)))
            lower  <- results[lower.tri(results)]
            lower.bool = sapply(lower, function(x) isTRUE(x > 0 & x < 1))
            expect_equal(sum(lower.bool), length(lower.bool))
            upper  <- results[upper.tri(results, diag = TRUE)]
            upper.bool = sapply(upper, function(x) isTRUE(x == 0))
            expect_equal(sum(upper.bool), length(upper.bool))
            lower.2  <- results.2[lower.tri(results.2)]
            upper.2  <- t(results.2)[lower.tri(results.2)]
            expect_equal(lower.2, lower)
            expect_equal(diag(results.2), rep.vec, ignore_attr = TRUE)
            expect_equal(sum(upper.2 < lower.2), 0)
          }
)

test_that("KrzCor returns correct results on lists + matrices",
          {
            mat.list <- lapply(as.list(1:10), function(x) RandomMatrix(13))
            y.matrix <- RandomMatrix(13)
            results <- KrzCor(mat.list, y.matrix)
            expect_equal(results[1], KrzCor(mat.list[[1]], y.matrix), ignore_attr = TRUE)
            expect_type(results, "double")
            expect_equal(sum(is.na(results)), 0)
            expect_equal(length(results), length(mat.list))
            names(mat.list) <- seq_along(mat.list)
            named.results <- KrzCor(mat.list, y.matrix)
            expect_s3_class(named.results, "data.frame")
            expect_equal(sum(is.na(named.results)), 0)
            expect_equal(dim(named.results), c(length(mat.list), 2))
            expect_equal(named.results[,".id"], names(mat.list))
          }
)

test_that("KrzProjection returns correct results on matrices",
          {
            local_seed(17, "3.5.0")
            cov.matrix.1 <- RandomMatrix(10, LKJ = FALSE)
            cov.matrix.2 <- RandomMatrix(10, LKJ = FALSE)
            expect_type(KrzProjection(cov.matrix.1, cov.matrix.2), "list")
            expect_equal(length(KrzProjection(cov.matrix.1, cov.matrix.2)), 2)
            expect_equal(length(KrzProjection(cov.matrix.1,
                                             cov.matrix.2,
                                             ret.dim.1 = 3)[[2]]), 3)
            expect_equal(length(KrzProjection(cov.matrix.1,
                                             cov.matrix.2,
                                             ret.dim.1 = 5)[[2]]), 5)
            expect_equal(length(KrzProjection(cov.matrix.1,
                                             cov.matrix.2,
                                             ret.dim.2 = 10)[[2]]), 4)
            expect_equal(KrzProjection(cov.matrix.1,
                                      cov.matrix.2,
                                      ret.dim.1 = 10, ret.dim.2 = 10)[[1]], 1)
            per.PC <- KrzProjection(cov.matrix.1, cov.matrix.2, ret.dim.1 = 10, ret.dim.2 = 10)[[2]]
            ones = sapply(per.PC, function(x) isTRUE(all.equal(x, 1)))
            expect_equal(sum(ones), 10)
          }
)

test_that("KrzProjection returns correct results on lists",
          {
            local_seed(170, "3.5.0")
            mat.list <- RandomMatrix(10, 10, LKJ = FALSE)
            expect_equal(dim(KrzProjection(mat.list)), c(length(mat.list), length(mat.list)))
            expect_equal(KrzProjection(mat.list)[1,2], KrzProjection(mat.list[[1]], mat.list[[2]])[[1]])
            expect_equal(KrzProjection(mat.list)[2,1], KrzProjection(mat.list[[2]], mat.list[[1]])[[1]])
            expect_equal(length(KrzProjection(mat.list, full.results = TRUE)), length(mat.list))
            expect_type(KrzProjection(mat.list, full.results = TRUE), "list")
          }
)

test_that("KrzProjection returns correct results on lists and matrices",
          {
            mat.list <- RandomMatrix(10, 10)
            results <- KrzProjection(mat.list, mat.list[[1]])
            test.results <- KrzProjection(mat.list)[,1]
            expect_s3_class(results, "data.frame")
            expect_equal(results[[1]], test.results, ignore_attr = TRUE)
          }
)

test_that("PCAsimilarity returns correct results",
          {
            cov.matrix.1 <- RandomMatrix(10)
            cov.matrix.2 <- RandomMatrix(10)
            eg.x <- eigen(cov.matrix.1)
            eg.y <- eigen(cov.matrix.2)
            eg.x.values <- eg.x$values
            eg.y.values <- eg.y$values
            eg.x.vectors <- eg.x$vectors
            eg.y.vectors <- eg.y$vectors

            total_var <- eg.x.values %*% eg.y.values

            SL = sum((eg.x.values %o% eg.y.values) * ((t(eg.x.vectors) %*% (eg.y.vectors))**2))/total_var
            expect_equal(PCAsimilarity(cov.matrix.1, cov.matrix.2), SL, ignore_attr = TRUE)
            expect_equal(PCAsimilarity(x  <- RandomMatrix(11), x), 1, ignore_attr = TRUE)
            expect_true(PCAsimilarity(cov.matrix.1, cov.matrix.2) <= 1)
            expect_true(PCAsimilarity(cov.matrix.1, cov.matrix.2) > 0)
          }
)

test_that("PCAsimilarity returns correct results",
          {
            mat.list <- lapply(as.list(1:10), function(x) RandomMatrix(10))
            rep.vec <- runif(length(mat.list), 0.8, 0.9)
            results <- PCAsimilarity(mat.list)
            expect_equal(results[2,1], PCAsimilarity(mat.list[[1]], mat.list[[2]]), ignore_attr = TRUE)
            results.2 <- PCAsimilarity(mat.list, repeat.vector = rep.vec)
            expect_equal(sum(is.na(results)), 0)
            expect_equal(sum(is.na(results.2)), 0)
            expect_equal(dim(results), c(length(mat.list),length(mat.list)))
            lower  <- results[lower.tri(results)]
            lower.bool = sapply(lower, function(x) isTRUE(x > 0 & x < 1))
            expect_equal(sum(lower.bool), length(lower.bool))
            upper  <- results[upper.tri(results, diag = TRUE)]
            upper.bool = sapply(upper, function(x) isTRUE(x == 0))
            expect_equal(sum(upper.bool), length(upper.bool))
            lower.2  <- results.2[lower.tri(results.2)]
            upper.2  <- t(results.2)[lower.tri(results.2)]
            expect_equal(lower.2, lower)
            expect_equal(diag(results.2), rep.vec, ignore_attr = TRUE)
            expect_equal(sum(upper.2 < lower.2), 0)
          }
)

test_that("PCAsimilarity returns correct results on lists + matrices",
          {
            mat.list <- lapply(as.list(1:10), function(x) RandomMatrix(13))
            y.matrix <- RandomMatrix(13)
            results <- PCAsimilarity(mat.list, y.matrix)
            expect_equal(results[1], PCAsimilarity(mat.list[[1]], y.matrix), ignore_attr = TRUE)
            expect_type(results, "double")
            expect_equal(sum(is.na(results)), 0)
            expect_equal(length(results), length(mat.list))
            names(mat.list) <- seq_along(mat.list)
            named.results <- PCAsimilarity(mat.list, y.matrix)
            expect_s3_class(named.results, "data.frame")
            expect_equal(sum(is.na(named.results)), 0)
            expect_equal(dim(named.results), c(length(mat.list), 2))
            expect_equal(named.results[,".id"], names(mat.list))
          }
)

# Reference values computed with evolqg 0.4-3 on the covariance matrices of the
# five dentus species (4 traits). The faster list methods must reproduce them.
dentus.covs <- function() {
  data(dentus, envir = environment())
  lapply(split(dentus[, 1:4], dentus$species), cov)
}
lower.part <- function(m) m[lower.tri(m)]
old.ref <- function(x) matrix(x, 5, 5)

test_that("list methods reproduce the values of evolqg 0.4-3",
          {
            covs <- dentus.covs()
            cors <- lapply(covs, cov2cor)
            krz <- structure(c(0, 0.9875424144, 0.8226835597, 0.751784972, 0.5889151745,
                               0, 0, 0.8124445851, 0.74235615, 0.6169841037, 0, 0, 0, 0.9468612189,
                               0.6767582496, 0, 0, 0, 0, 0.8036694634, 0, 0, 0, 0, 0), dim = c(5L, 5L))
            pca <- structure(c(0, 0.9834567252, 0.9438200613, 0.980587512, 0.8498867416,
                               0, 0, 0.9153103824, 0.9500547756, 0.9083482774, 0, 0, 0, 0.9528737037,
                               0.8522773712, 0, 0, 0, 0, 0.847976552, 0, 0, 0, 0, 0), dim = c(5L, 5L))
            matcor <- structure(c(0, 0.8841707227, 0.5770244019, 0.3920545785, -0.209902053,
                                  0, 0, 0.2655182769, -0.0730655279, 0.1245137034, 0, 0, 0, 0.6646935662,
                                  -0.2866739463, 0, 0, 0, 0, -0.5744876086, 0, 0, 0, 0, 0), dim = c(5L, 5L))
            kproj <- structure(c(0.5942705482, 0.4450025778, 0.4118931125, 0.7248974186,
                                 0.2932141237, 0.5676125908, 0.465902149, 0.3629071037, 0.6835473828,
                                 0.3527861042, 0.5030527348, 0.3474841184, 0.4865810855, 0.6241237591,
                                 0.2419547813, 0.5875669319, 0.434370101, 0.4142118544, 0.7331678536,
                                 0.2764561122, 0.4302607448, 0.4058530515, 0.2907047066, 0.500486356,
                                 0.4049835365), dim = c(5L, 5L))
            expect_equal(unname(KrzCor(covs, ret.dim = 2)), krz, tolerance = 1e-9)
            expect_equal(unname(PCAsimilarity(covs)), pca, tolerance = 1e-9)
            expect_equal(unname(MatrixCor(cors)), matcor, tolerance = 1e-9)
            expect_equal(unname(KrzProjection(covs)), kproj, tolerance = 1e-9)
            expect_equal(KrzCor(covs, covs[[1]])$krz,
                         c(1, 0.9551417154, 0.84650457, 0.9887195886, 0.7240149224),
                         tolerance = 1e-9)
            expect_equal(dimnames(KrzCor(covs, ret.dim = 2)), list(names(covs), names(covs)))
          }
)

test_that("list methods agree with the pairwise methods",
          {
            local_seed(11)
            mats <- RandomMatrix(6, 5)
            cors <- lapply(mats, cov2cor)
            pairs <- function(f, ...) {
              out <- matrix(0, 5, 5)
              for (i in 1:4) for (j in (i + 1):5) out[j, i] <- f(mats[[j]], mats[[i]], ...)
              out
            }
            low <- lower.tri(diag(5))
            expect_equal(KrzCor(mats)[low], pairs(KrzCor)[low], ignore_attr = TRUE)
            expect_equal(PCAsimilarity(mats)[low], pairs(PCAsimilarity)[low], ignore_attr = TRUE)
            expect_equal(MatrixCor(cors)[low], pairs(MatrixCor)[low], ignore_attr = TRUE)
            expect_equal(unname(KrzProjection(mats)),
                         outer(1:5, 1:5, Vectorize(function(i, j) KrzProjection(mats[[i]], mats[[j]])[[1]])))
            expect_equal(KrzCor(mats, mats[[2]]), vapply(mats, KrzCor, 1, mats[[2]]), ignore_attr = TRUE)
            expect_equal(MatrixCor(cors, cors[[2]])[[1]], vapply(cors, MatrixCor, 1, cors[[2]]),
                         ignore_attr = TRUE)
          }
)

test_that("RandomSkewers on lists shares the random vectors and is reproducible",
          {
            mats <- RandomMatrix(6, 4)
            local_seed(5)
            a <- RandomSkewers(mats, num.vectors = 2000)
            local_seed(5)
            b <- RandomSkewers(mats, num.vectors = 2000)
            expect_equal(a, b)
            expect_equal(names(a), c("correlations", "probabilities"))
            pair <- RandomSkewers(mats[[2]], mats[[1]], num.vectors = 20000)
            expect_equal(a$correlations[2, 1], pair[["correlation"]], tolerance = 0.05)
            expect_true(all(a$probabilities[lower.tri(a$probabilities)] >= 0))
            res <- RandomSkewers(mats, mats[[1]], num.vectors = 500)
            expect_equal(dim(res), c(4, 3))
            expect_equal(res$correlation[1], 1)
            expect_error(RandomSkewers(mats, RandomMatrix(5), num.vectors = 10), "do not match")
            expect_error(RS(mats[[1]], RandomMatrix(5), 10), "same dimensions")
            expect_error(RS(mats[[1]], mats[[2]], 0), "positive")
          }
)

test_that("mcmc_sample methods return one value per sample",
          {
            local_seed(21)
            n <- 4
            make.sample <- function() {
              a <- array(0, c(n, 5, 5))
              for (i in seq_len(n)) a[i, , ] <- RandomMatrix(5)
              class(a) <- "mcmc_sample"
              a
            }
            s1 <- make.sample(); s2 <- make.sample()
            one <- function(s, i) unclass(s)[i, , ]
            rs <- RandomSkewers(s1, s2, num.vectors = 200)
            expect_length(rs, n)
            expect_true(all(abs(rs) <= 1))
            krz <- KrzCor(s1, s2)
            expect_equal(krz, vapply(seq_len(n), function(i) KrzCor(one(s1, i), one(s2, i)), 1),
                         ignore_attr = TRUE)
            pca <- PCAsimilarity(s1, s2)
            expect_equal(pca, vapply(seq_len(n), function(i) PCAsimilarity(one(s1, i), one(s2, i)), 1),
                         ignore_attr = TRUE)
            mc <- MatrixCor(s1, s2)
            expect_length(mc, n)
            expect_equal(mc, vapply(seq_len(n), function(i) MatrixCor(cov2cor(one(s1, i)), cov2cor(one(s2, i))), 1),
                         ignore_attr = TRUE)
            expect_equal(MantelCor(s1, s2), mc)
            expect_equal(MantelCor(s1, s2, parallel = FALSE), mc)
            # against a single matrix: one row per sample, and correlations
            y <- RandomMatrix(5)
            single <- MatrixCor(s1, y)
            expect_s3_class(single, "data.frame")
            expect_equal(single$correlations,
                         vapply(seq_len(n), function(i) MatrixCor(cov2cor(one(s1, i)), cov2cor(y)), 1))
            expect_equal(nrow(KrzCor(s1, y)), n)
            expect_equal(nrow(PCAsimilarity(s1, y)), n)
            expect_equal(nrow(RandomSkewers(s1, y, num.vectors = 100)), n)
            short <- s2[1:2, , ]; class(short) <- "mcmc_sample"
            expect_error(KrzCor(s1, short), "same size")
          }
)

test_that("MantelCor with withinLandmark excludes the diagonal and uses all other elements",
          {
            local_seed(3)
            x <- cov2cor(RandomMatrix(6)); y <- cov2cor(RandomMatrix(6))
            res <- MantelCor(x, y, permutations = 20, landmark.dim = 2, withinLandmark = TRUE)
            expect_equal(res[["rho"]], cor(x[lower.tri(x)], y[lower.tri(y)]))
            # a diagonal in the mask would make the correlation differ
            full <- lower.tri(x, diag = TRUE)
            expect_false(isTRUE(all.equal(res[["rho"]], cor(x[full], y[full]))))
            between <- MantelCor(x, y, permutations = 20, landmark.dim = 2)
            expect_equal(between[["rho"]], cor(x[evolqg:::lower.tri.land(x, 2)], y[evolqg:::lower.tri.land(x, 2)]))
            expect_equal(evolqg:::CreateWithinLandMat(2, 2),
                         kronecker(diag(2), matrix(1, 2, 2)) == 1)
          }
)

test_that("MantelCor landmark p-value is never 0 and matches (k + 1)/(permutations + 1)",
          {
            local_seed(4)
            x <- cov2cor(RandomMatrix(6)); 
            res <- MantelCor(x, x, permutations = 50, landmark.dim = 3)
            expect_equal(res[["rho"]], 1)
            expect_gt(res[["Probability"]], 0)
            expect_equal(res[["Probability"]] * 51, round(res[["Probability"]] * 51))
          }
)

test_that("MantelCor warns on non correlation matrices only",
          {
            local_seed(6)
            x <- cov2cor(RandomMatrix(5)); y <- RandomMatrix(5)
            expect_warning(MantelCor(x, x, permutations = 10), NA)
            expect_warning(MantelCor(x, y * 5, permutations = 10), "correlation matrices")
            # diagonal equal to 1 only up to rounding error is accepted
            x2 <- x; diag(x2) <- 1 + 1e-12
            expect_warning(MantelCor(x, x2, permutations = 10), NA)
            expect_warning(MatrixCor(x, x2), NA)
            expect_warning(MatrixCor(x, y * 5), "correlation matrices")
          }
)

test_that("KrzCor and KrzProjection work for 2 and 3 traits",
          {
            local_seed(8)
            for (n in 2:3) {
              random.cov <- function() crossprod(matrix(rnorm(n * n), n))
              a <- random.cov(); b <- random.cov()
              expect_true(is.finite(KrzCor(a, b)))
              expect_equal(KrzCor(a, b), KrzCor(a, b, ret.dim = 1))
              expect_equal(KrzCor(a, a)[[1]], 1)
              kp <- KrzProjection(a, b)
              expect_true(is.finite(kp$total.variation))
              expect_length(kp$per.PC, 1)
              expect_equal(kp, KrzProjection(a, b, ret.dim.1 = 1, ret.dim.2 = 1))
              ml <- list(a, b, random.cov())
              expect_true(all(is.finite(KrzCor(ml))))
              expect_true(all(is.finite(KrzProjection(ml))))
            }
            # the default for other sizes is unchanged
            expect_equal(evolqg:::.default_ret_dim(10), 4)
            expect_equal(evolqg:::.default_ret_dim(11), 4)
            expect_equal(evolqg:::.default_ret_dim(4), 1)
          }
)

test_that("comparison methods handle partially named, duplicated and single lists",
          {
            local_seed(9)
            mats <- RandomMatrix(5, 4)
            partial <- mats; names(partial) <- c("a", "", "c", "")
            res <- KrzCor(partial)
            expect_equal(dimnames(res), list(c("a", "2", "c", "4"), c("a", "2", "c", "4")))
            dup <- mats; names(dup) <- c("a", "a", "b", "a")
            expect_equal(rownames(KrzCor(dup)), c("a", "a.1", "b", "a.2"))
            expect_equal(rownames(KrzProjection(dup)), c("a", "a.1", "b", "a.2"))
            expect_equal(colnames(KrzProjection(dup)), c("a", "a.1", "b", "a.2"))
            expect_equal(KrzProjection(dup)[2, 1], KrzProjection(mats[[2]], mats[[1]])[[1]])
            expect_equal(names(KrzProjection(partial, full.results = TRUE)), c("a", "2", "c", "4"))
            expect_equal(rownames(KrzProjection(partial)), c("a", "2", "c", "4"))
            expect_equal(rownames(PCAsimilarity(dup)), c("a", "a.1", "b", "a.2"))
            single <- KrzCor(mats[1])
            expect_equal(dim(single), c(1, 1))
            expect_equal(dim(KrzProjection(mats[1])), c(1, 1))
            expect_equal(dim(MantelCor(lapply(mats[1], cov2cor), permutations = 10)[[1]]), c(1, 1))
            expect_equal(diag(KrzCor(mats[1], repeat.vector = 0.9)), 0.9, ignore_attr = TRUE)
            expect_error(KrzCor(list()), "at least one")
            # a data frame keeps working for a named list
            expect_s3_class(KrzCor(partial, mats[[1]]), "data.frame")
          }
)

test_that("SingleComparisonMap checks both dimensions and that y is a matrix",
          {
            mats <- list(diag(3), diag(3))
            expect_error(evolqg:::SingleComparisonMap(mats, 1:3, function(x, y) 1), "do not match")
            expect_error(evolqg:::SingleComparisonMap(mats, matrix(1, 3, 2), function(x, y) 1), "do not match")
            expect_error(evolqg:::SingleComparisonMap(mats, diag(4), function(x, y) 1), "do not match")
            expect_s3_class(evolqg:::SingleComparisonMap(mats, diag(3), function(x, y) c(a = 1)), "data.frame")
          }
)

test_that("parallel = TRUE works through a registered backend",
          {
            skip_if_not_installed("foreach")
            skip_if_not_installed("doParallel")
            local_seed(12)
            mats <- RandomMatrix(5, 4)
            doParallel::registerDoParallel(cores = 2)
            withr::defer(doParallel::stopImplicitCluster())
            expect_equal(KrzCor(mats, parallel = TRUE), KrzCor(mats))
            expect_equal(KrzProjection(mats, parallel = TRUE), KrzProjection(mats))
            expect_equal(KrzProjection(mats, mats[[1]], parallel = TRUE), KrzProjection(mats, mats[[1]]))
            expect_equal(dim(RandomSkewers(mats, num.vectors = 50, parallel = TRUE)[[1]]), c(4, 4))
          }
)
