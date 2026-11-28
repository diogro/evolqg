# Tiny synthetic landmark set: 8 landmarks in 3D and 5 in 2D
make_shapes <- function(k, m, n, sd = 0.05) {
  ref <- matrix(rnorm(k * m), k, m)
  arr <- array(NA_real_, c(k, m, n), dimnames = list(NULL, NULL, paste0("s", seq_len(n))))
  for (i in seq_len(n)) arr[, , i] <- ref + matrix(rnorm(k * m, sd = sd), k, m)
  list(ref = ref, arr = arr)
}
tess3 <- rbind(c(1, 2, 3, 4), c(2, 3, 4, 5), c(5, 6, 7, 8), c(1, 4, 6, 8))
tess2 <- rbind(c(1, 2, 3), c(2, 3, 4), c(3, 4, 5))

# central finite-difference jacobian of the TPS function
fd_jacobian <- function(f, x, h = 1e-6) {
  vapply(seq_along(x), function(j) {
    e <- replace(numeric(length(x)), j, h)
    as.numeric(f(x + e) - f(x - e)) / (2 * h)
  }, numeric(length(f(x))))
}
centroids <- function(Q, tess) t(apply(tess, 1, function(r) colMeans(Q[r, , drop = FALSE])))

test_that("TPS interpolates the target shape", {
  local_seed(1)
  s <- make_shapes(8, 3, 2)
  tp <- TPS(s$arr[, , 1], s$ref)
  fitted <- t(apply(s$ref, 1, function(p) tp$theta.p(p)))
  expect_equal(fitted, s$arr[, , 1], tolerance = 1e-8, ignore_attr = TRUE)
  s2 <- make_shapes(5, 2, 2)
  tp2 <- TPS(s2$arr[, , 1], s2$ref)
  fitted2 <- t(apply(s2$ref, 1, function(p) tp2$theta.p(p)))
  expect_equal(fitted2, s2$arr[, , 1], tolerance = 1e-8, ignore_attr = TRUE)
})

test_that("JacobianArray matches finite differences in 3D and 2D", {
  local_seed(2)
  for (case in list(list(k = 8, m = 3, tess = tess3), list(k = 5, m = 2, tess = tess2))) {
    s <- make_shapes(case$k, case$m, 2)
    tp <- TPS(s$arr[, , 1], s$ref)
    jac <- JacobianArray(tp, case$tess)
    expect_equal(dim(jac), c(case$m, case$m, nrow(case$tess)))
    cent <- centroids(s$ref, case$tess)
    for (i in seq_len(nrow(cent))) {
      fd <- fd_jacobian(tp$theta.p, cent[i, ])
      expect_equal(jac[, , i], fd, tolerance = 1e-6, ignore_attr = TRUE)
    }
  }
})

test_that("JacobianArray warns when numerical differentiation arguments are given", {
  local_seed(3)
  s <- make_shapes(8, 3, 2)
  tp <- TPS(s$arr[, , 1], s$ref)
  expect_warning(JacobianArray(tp, tess3, method = "Richardson"), "ignored")
})

test_that("Center2MeanJacobianFast returns mean-zero log determinants", {
  local_seed(4)
  A <- array(rnorm(3 * 3 * 7, sd = 0.1), c(3, 3, 7)) + array(diag(3), c(3, 3, 7))
  out <- Center2MeanJacobianFast(A)
  ld <- apply(A, 3, function(x) log(abs(det(x))))
  expect_equal(out, ld - mean(ld))
  expect_equal(mean(out), 0)
})

test_that("Rotate2MidlineMatrix returns an orthogonal matrix", {
  local_seed(5)
  X <- matrix(rnorm(8 * 3), 8, 3)
  R <- Rotate2MidlineMatrix(X, 1:4)
  expect_equal(dim(R), c(3, 3))
  expect_equal(crossprod(R), diag(3), tolerance = 1e-10)
})

test_that("LocalShapeVariables gives mean-zero local variables", {
  local_seed(6)
  s <- make_shapes(8, 3, 6)
  cs <- apply(s$arr, 3, function(x) sqrt(sum(sweep(x, 2, colMeans(x))^2)))
  out <- LocalShapeVariables(gpa = s$arr, cs = cs, tesselation = tess3, verbose = FALSE)
  expect_equal(names(out), c("tps", "jacobians", "local", "reference", "cs", "info"))
  expect_equal(dim(out$local), c(6, nrow(tess3)))
  expect_equal(dim(out$jacobians), c(3, 3, 6, nrow(tess3)))
  expect_equal(unname(colMeans(out$local)), rep(0, nrow(tess3)), tolerance = 1e-12)
  expect_equal(out$info, paste0("s", 1:6))
  expect_equal(out$reference, rowMeans(s$arr, dims = 2), ignore_attr = TRUE)
  # local variables are the centred log determinants of the jacobians
  ld <- apply(out$jacobians[, , , 2], 3, function(x) log(abs(det(x))))
  expect_equal(unname(out$local[, 2]), unname(ld - mean(ld)))
  expect_message(LocalShapeVariables(gpa = s$arr, cs = cs, tesselation = tess3), "tps done")
})

test_that("LocalShapeVariables input checks and cs from gpa", {
  local_seed(7)
  s <- make_shapes(8, 3, 4)
  expect_error(LocalShapeVariables(cs = 1:4, tesselation = tess3), "gpa must be supplied")
  expect_error(LocalShapeVariables(tesselation = tess3), "Supply either")
  expect_warning(out <- LocalShapeVariables(gpa = s$arr, tesselation = tess3, verbose = FALSE),
                 "computed from gpa")
  expect_length(out$cs, 4)
  out2 <- LocalShapeVariables(gpa = s$arr, landmarks = s$arr, tesselation = tess3, verbose = FALSE)
  expect_length(out2$cs, 4)
})

test_that("LocalShapeVariables runs the GPA when only landmarks are given", {
  skip_if_not_installed("Morpho")
  local_seed(8)
  s <- make_shapes(8, 3, 5)
  out <- LocalShapeVariables(landmarks = s$arr, tesselation = tess3, verbose = FALSE)
  expect_equal(dim(out$local), c(5, nrow(tess3)))
  expect_equal(unname(out$cs), unname(apply(s$arr, 3, Morpho::cSize)))
  expect_equal(unname(colMeans(out$local)), rep(0, nrow(tess3)), tolerance = 1e-12)
})
