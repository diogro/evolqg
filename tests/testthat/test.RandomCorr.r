test_that("RandomCorr returns sensible results",
          {
            local_seed(1)
            mat.list <- lapply(as.list(1:100), function(x) RandomMatrix(10))
            mean.corr <- mean(MatrixCor(mat.list)[upper.tri(array(0, c(100, 100)))])
            expect_equal(dim(mat.list[[1]]), c(10, 10))
            expect_equal(sum(diag(mat.list[[1]])), 10)
            expect_true(abs(mean.corr) < 0.01)
          }
)

test_that("RandomMatrix returns sensible results",
          {
            mat.list <- RandomMatrix(10, 10)
            expect_equal(dim(mat.list[[1]]), c(10, 10))
            expect_equal(length(mat.list), 10)
            expect_equal(sum(diag(mat.list[[1]])), 10)
            rand.mat <- RandomMatrix(10, 1, 1, 2)
            expect_true(sum(diag(rand.mat)) > 10)
            expect_true(sum(diag(rand.mat)) < 20)
            mat.list.fixed.var <- RandomMatrix(3, 2, variance = c(1, 2, 3))
            expect_equal(diag(mat.list.fixed.var[[1]]), diag(mat.list.fixed.var[[2]]))
            expect_equal(diag(mat.list.fixed.var[[1]]), c(1, 2, 3))
          }
)

test_that("RandomMatrix works for one and two traits",
          {
            local_seed(1)
            expect_equal(RandomMatrix(1), matrix(1, 1, 1))
            expect_equal(RandomMatrix(1, 1, 2, 2), matrix(2, 1, 1))
            for (LKJ in c(FALSE, TRUE)) {
              m <- RandomMatrix(2, LKJ = LKJ)
              expect_equal(dim(m), c(2, 2))
              expect_equal(diag(m), c(1, 1))
              expect_true(isSymmetric(m))
              expect_true(all(eigen(m)$values > 0))
              expect_length(RandomMatrix(2, 3, LKJ = LKJ), 3)
            }
          }
)

test_that("RandomMatrix gives positive definite matrices for both methods",
          {
            local_seed(2)
            for (n in c(3, 6, 12)) {
              for (LKJ in c(FALSE, TRUE)) {
                m <- RandomMatrix(n, LKJ = LKJ)
                expect_equal(diag(m), rep(1, n))
                expect_gt(min(eigen(m, symmetric = TRUE)$values), -1e-10)
              }
            }
          }
)

test_that("The partial correlation method keeps its random stream",
          {
            ## reference values computed with evolqg 0.4-3
            local_seed(1)
            m <- RandomMatrix(4)
            expect_equal(m[lower.tri(m)],
                         c(-0.149986846186699, 0.817026716360302, -0.2557522, -0.589663622323468,
                           0.14570672, -0.46898268), tolerance = 1e-6)
          }
)

test_that("LKJ marginals are Beta(eta - 1 + n/2) on (-1, 1) and not capped at 0.5",
          {
            local_seed(3)
            for (pars in list(c(4, 1), c(4, 3))) {
              n <- pars[1]; eta <- pars[2]
              r <- replicate(3000, { m <- RandomMatrix(n, LKJ = TRUE, shape = eta); c(m[1, 2], m[1, 3], m[2, 3], m[1, 4]) })
              a <- eta - 1 + n / 2
              for (idx in 1:4) {
                x <- (r[idx, ] + 1) / 2
                expect_gt(suppressWarnings(ks.test(x, "pbeta", a, a))$p.value, 1e-3)
              }
            }
            r1 <- replicate(2000, RandomMatrix(3, LKJ = TRUE, shape = 1)[1, 2])
            expect_gt(max(abs(r1)), 0.9)
          }
)
