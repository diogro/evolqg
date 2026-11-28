test_that("Normalize returns correct results",
          {
            beta <- rep(1, 10)
            expect_equal(Normalize(beta), rep(1/sqrt(10), 10))
            expect_equal(Norm(Normalize(beta)), 1)
          }
)

test_that("Norm returns correct results",
          {
            beta <- rep(1, 10)
            expect_equal(Norm(beta), sqrt(sum((rep(1, 10))^2)))
            expect_equal(Norm(rep(1/sqrt(10), 10)), 1)
          }
)
