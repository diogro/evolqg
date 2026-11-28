# The C++ routines in src/fast_utils.cpp replace R loops. These tests keep the
# old R loops as references and check that results and the random number
# stream are unchanged.

ref_rand_corr <- function(num.traits, ke = 10^-3) { # R loop of RandCorr() before C++
  random.corr = array(0., dim = c(num.traits, num.traits))
  b = array(0., dim = c(num.traits, num.traits))
  b[lower.tri(b, diag = TRUE)] = 1.
  random.corr[2:num.traits, 1] = -1 + 2*round(runif(num.traits-1)*10^8)/10^8
  b[2:num.traits, 1] = random.corr[2:num.traits, 1]
  for (i in 2:num.traits) b[i, 2:i] = sqrt(1 - random.corr[i, 1]^2)
  for (i in seq_len(num.traits)[-(1:2)]) {
    for (j in 2:(i-1)) {
      b1 = b[j, 1:(j-1)] %*% b[i, 1:(j-1)]
      b2 = b[j, j]*b[i, j]
      z = b1 + b2; y = b1 - b2
      if (b2 < ke) random.corr[i, j] = b1
      else random.corr[i, j] = y + (z - y)*round(runif(1)*10^8)/10^8
      cosinv = c((random.corr[i, j] - b1)/b2)
      if (is.finite(cosinv)) {
        if (cosinv > 1) b[i, (j+1):num.traits] = 0
        else if (cosinv < -1) { b[i, j] = -b[i, j]; b[i, (j+1):num.traits] = 0 }
        else { b[i, j] = b[i, j]*cosinv; b[i, (j+1):num.traits] = b[i, (j+1):num.traits]*sqrt(1 - cosinv^2) }
      }
    }
  }
  random.corr + t(random.corr) + diag(rep(1, num.traits))
}

ref_rand_lkj <- function(k, shape) { # R loop of RandLKJ() before C++
  P <- matrix(0, k, k); rc <- diag(k)
  for (m in seq_len(k - 1)) {
    beta <- shape + (k - 1 - m) / 2
    for (i in (m + 1):k) {
      P[m, i] <- 2 * rbeta(1, beta, beta) - 1
      p <- P[m, i]
      if (m > 1) for (l in rev(seq_len(m - 1))) p <- p * sqrt((1 - P[l, i]^2) * (1 - P[l, m]^2)) + P[l, i] * P[l, m]
      rc[m, i] <- rc[i, m] <- p
    }
  }
  rc
}

test_that("C++ partial-correlation sampler matches the R loop and its random stream", {
  for (p in c(2, 3, 5, 12, 40)) {
    local_seed(p)
    a <- ref_rand_corr(p); after.a <- runif(1)
    local_seed(p)
    b <- .rand_corr_partial_cpp(as.integer(p), 10^-3); after.b <- runif(1)
    expect_equal(b, a, tolerance = 1e-12)
    expect_identical(after.b, after.a)
  }
})

test_that("C++ LKJ sampler matches the R loop and its random stream", {
  for (p in c(2, 3, 6, 25)) for (eta in c(1, 2.5)) {
    local_seed(p)
    a <- ref_rand_lkj(p, eta); after.a <- runif(1)
    local_seed(p)
    b <- .rand_lkj_cpp(as.integer(p), eta); after.b <- runif(1)
    expect_equal(b, a, tolerance = 1e-12)
    expect_identical(after.b, after.a)
  }
})

test_that("C++ Mantel permutation correlations match cor()", {
  local_seed(3)
  x <- cov2cor(crossprod(matrix(rnorm(200), 20)))
  y <- cov2cor(crossprod(matrix(rnorm(200), 20)))
  low <- lower.tri(x)
  perms <- t(replicate(50, sample.int(10)))
  expect_equal(.mantel_perm_cor(x, y[low], low, perms),
               apply(perms, 1, function(p) cor(x[p, p][low], y[low])), tolerance = 1e-12)
  # any mask, as in the landmark test
  mask <- low & !diag(10)[, c(2:10, 1)]
  expect_equal(.mantel_perm_cor(x, y[mask], mask, perms),
               apply(perms, 1, function(p) cor(x[p, p][mask], y[mask])), tolerance = 1e-12)
  expect_equal(.mantel_perm_cor(x, x[low], low, matrix(1:10, 1)), 1)
  expect_length(.mantel_perm_cor(x, y[low], low, matrix(integer(0), 0, 10)), 0)
})
