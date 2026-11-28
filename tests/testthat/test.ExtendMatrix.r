test_that("ExtendMatrix returns correct results",
          {
            local_seed(43, "3.5.0")
            cov.matrix = RandomMatrix(11)
            ret.dim = sample(2:11, 1)
            ext.matrix = ExtendMatrix(cov.matrix, ret.dim = ret.dim)
            expect_type(ext.matrix, "list")
            last.eval = eigen(cov.matrix)$values[ret.dim]
            eVals = eigen(ext.matrix[[1]])$values
            extended = sapply(eVals, function(x) isTRUE(all.equal(x, last.eval)))
            expect_equal(sum(extended), 11-ret.dim+1)
            expect_equal(sum(eVals>0), length(eVals))
            expect_warning(ExtendMatrix(cov.matrix[1:9, 1:9], ret.dim = 5), "matrix is too small")
          }
)

test_that("ExtendMatrix guards against small matrices and empty extensions",
          {
            local_seed(1)
            expect_error(ExtendMatrix(RandomMatrix(6)), "too small")
            expect_error(ExtendMatrix(RandomMatrix(3)), "too small")
            cov.matrix <- RandomMatrix(11)
            expect_warning(res <- ExtendMatrix(cov.matrix, var.cut.off = 0), "no eigenvalues were extended")
            expect_equal(res$ExtMat, cov.matrix, tolerance = 1e-10)
            expect_equal(dim(res$var.grad), c(5, 2))
          }
)
