dentus_vcv <- function() {
  data(dentus, envir = environment())
  aperm(daply(dentus, .(species), function(x) cov(x[,-5])), c(2, 3, 1))
}

test_that("EigenTensorDecomposition returns the expected structure",
          {
            d <- dentus_vcv()
            etd <- EigenTensorDecomposition(d)
            expect_named(etd, c("mean", "mean.sqrt", "values", "matrices", "projection"))
            expect_equal(dim(etd$matrices), c(4, 4, 4))
            expect_equal(dimnames(etd$matrices)[[3]], paste0("ET", 1:4))
            expect_equal(dim(etd$projection), c(5, 4))
            expect_equal(rownames(etd$projection), LETTERS[1:5])
            expect_equal(etd$mean.sqrt %*% etd$mean.sqrt, etd$mean, ignore_attr = TRUE)
            expect_null(EigenTensorDecomposition(d, FALSE)$projection)
          }
)

test_that("EigenTensorDecomposition accepts lists",
          {
            d <- dentus_vcv()
            etd.list <- EigenTensorDecomposition(lapply(1:5, function(i) d[,,i]))
            expect_equal(etd.list$values, EigenTensorDecomposition(d)$values)
          }
)

test_that("Eigentensor projections are PCA-like scores in log space",
          {
            local_seed(1)
            etd <- EigenTensorDecomposition(simplify2array(RandomMatrix(5, 40, 1, 3)))
            proj <- etd$projection
            expect_equal(unname(apply(proj, 2, var)), etd$values)
            cors <- cor(proj)
            expect_lt(max(abs(cors - diag(ncol(cors)))), 1e-8)
            expect_lt(max(abs(colMeans(proj))), 1e-8)
            et <- matrix(etd$matrices, 25)
            expect_equal(crossprod(et), diag(ncol(et)), ignore_attr = TRUE)
            expect_true(all(diff(etd$values) <= 0))
          }
)

test_that("EigenTensorDecomposition works with two matrices",
          {
            d <- dentus_vcv()
            etd <- EigenTensorDecomposition(d[,,1:2])
            expect_equal(dim(etd$matrices), c(4, 4, 1))
            expect_equal(dim(etd$projection), c(2, 1))
            expect_length(etd$values, 1)
          }
)

test_that("ProjectMatrix and RevertMatrix round-trip",
          {
            d <- dentus_vcv()
            etd <- EigenTensorDecomposition(d)
            expect_equal(ProjectMatrix(d[,,2], etd), etd$projection[2, ], tolerance = 1e-8)
            expect_equal(RevertMatrix(etd$projection[3, ], etd, scaled = FALSE), d[,,3], tolerance = 1e-8, ignore_attr = TRUE)
            scores <- c(1, -0.5, 0.2, 0.1)
            rebuilt <- RevertMatrix(scores, etd, scaled = FALSE)
            expect_true(isSymmetric(rebuilt))
            expect_equal(unname(ProjectMatrix(rebuilt, etd)), scores, tolerance = 1e-8)
            expect_equal(ProjectMatrix(RevertMatrix(scores, etd, scaled = TRUE), etd) / sqrt(etd$values),
                         scores, tolerance = 1e-8, ignore_attr = TRUE)
            expect_equal(RevertMatrix(0, etd), etd$mean, tolerance = 1e-8, ignore_attr = TRUE)
            expect_equal(dimnames(rebuilt), dimnames(etd$matrices)[1:2])
          }
)

test_that("RevertMatrix pads and truncates scores",
          {
            etd <- EigenTensorDecomposition(dentus_vcv())
            expect_equal(RevertMatrix(1, etd), RevertMatrix(c(1, 0, 0, 0), etd))
            expect_warning(RevertMatrix(rep(1, 6), etd), "More scores")
          }
)
