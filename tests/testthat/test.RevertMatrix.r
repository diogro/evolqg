# RevertMatrix and ProjectMatrix are tested together with
# EigenTensorDecomposition in test.EigenTensorDecomposition.r
test_that("frobenius and maximum norm helpers work",
          {
            expect_equal(frobenius.prod(diag(2), matrix(1:4, 2)), 5)
            expect_equal(maximum.norm(c(-3, 2)), 3)
            expect_equal(frobenius.norm(matrix(1, 2, 2)), 2)
          }
)
