# a covariance matrix with three modules and a few random extra hypotheses
.mint_example <- function() {
  modules <- matrix(c(rep(c(1, 0, 0), each = 5),
                      rep(c(0, 1, 0), each = 5),
                      rep(c(0, 0, 1), each = 5)), 15)
  cor.hypot <- CreateHypotMatrix(modules)[[4]]
  hypot.mask <- matrix(as.logical(cor.hypot), 15, 15)
  mod.cor <- matrix(NA, 15, 15)
  mod.cor[ hypot.mask] <- runif(length(mod.cor[ hypot.mask]), 0.8, 0.9) # within-modules
  mod.cor[!hypot.mask] <- runif(length(mod.cor[!hypot.mask]), 0.3, 0.4) # between-modules
  diag(mod.cor) <- 1
  mod.cor <- (mod.cor + t(mod.cor))/2 # correlation matrices should be symmetric
  list(modules = modules, mod.cor = mod.cor)
}

test_that("MINT returns correct results", {
  local_seed(42, "3.5.0")
  ex <- .mint_example()
  modules <- ex$modules; mod.cor <- ex$mod.cor
  hypothetical.modules <- cbind(modules, matrix(sample(c(1, 0), 4*15, replace = TRUE), 15, 4))
  colnames(hypothetical.modules) <- letters[1:7]
  results <- MINT(mod.cor, hypothetical.modules)
  expect_true(all(results[[2]] == modules))

  random_var <- runif(15, 1, 10)
  root <- chol(sqrt(outer(random_var, random_var)) * mod.cor)
  out <- vapply(1:50, function(i) {
    mod.cov <- cov(.rmvnorm(50, sigma = crossprod(root)))
    hypothetical.modules <- cbind(modules, matrix(sample(c(1, 0), 4*15, replace = TRUE), 15, 4))
    MINT(mod.cov, hypothetical.modules)[[1]][1, 1]
  }, character(1))
  close_enough <- sum(grepl("1_2_3|2_3|1_2|1_3", out, perl = TRUE))/50
  expect_true(close_enough > 0.95)

  xdata <- .rmvnorm(200, sigma = crossprod(root))
  out_jack <- JackKnifeMINT(xdata, hypothetical.modules, 100)
  expect_equal(out_jack[1, 1], "a_b_c")
  expect_true(out_jack[1, 5] > 0.95)

  expect_error(JackKnifeMINT(xdata, hypothetical.modules[, 1]))
  expect_error(JackKnifeMINT(cov(xdata), hypothetical.modules))
})

test_that("JackKnifeMINT upper bounds are above the estimate", {
  local_seed(3)
  ex <- .mint_example()
  hyp <- cbind(ex$modules, matrix(sample(c(1, 0), 4*15, replace = TRUE), 15, 4))
  X <- matrix(rnorm(200*15), 200) %*% chol(ex$mod.cor)
  jack <- JackKnifeMINT(X, hyp, n = 40)
  # the upper limit is the 97.5% quantile (it used to be the 9.75% quantile)
  expect_true(all(jack$upper.raw >= jack$raw.gamma))
  expect_true(all(jack$lower.raw <= jack$raw.gamma))
  expect_true(all(jack$upper.corrected >= jack$lower.corrected))
  expect_true(mean(jack$upper.corrected >= jack$corrected.gamma) > 0.9)
})

test_that("MINT works with underscores in hypothesis names", {
  local_seed(1)
  ex <- .mint_example()
  hyp <- ex$modules
  colnames(hyp) <- c("mod_1", "mod_2", "mod_3")
  res <- MINT(ex$mod.cor, hyp)
  expect_true(all(res$Modularity_hypothesis == ex$modules[, 1:2]))
  expect_equal(colnames(res$Modularity_hypothesis), colnames(hyp)[1:2])
  expect_true("mod_1_mod_2_mod_3" %in% res$gamma_rank$.id)
  expect_equal(res$gamma_rank$.id[1], "mod_1_mod_2")
})

test_that("MINT handles the null hypothesis winning", {
  local_seed(1)
  ex <- .mint_example()
  fake <- function(cor.hypot.list, ...)
    data.frame(.id = names(cor.hypot.list), .n = seq_along(cor.hypot.list),
               rank = seq_along(cor.hypot.list),
               corrected.gamma = c(-1, rep(0, length(cor.hypot.list) - 1)),
               raw.gamma = 1)
  testthat::local_mocked_bindings(SingleMINT = fake)
  expect_warning(res <- MINT(ex$mod.cor, ex$modules), "null hypothesis")
  expect_equal(res$gamma_rank$.id[1], "null")
  expect_equal(unname(res$Modularity_hypothesis), diag(15))
})

test_that("MINT significance columns", {
  local_seed(5)
  ex <- .mint_example()
  S <- 4 * ex$mod.cor * CreateHypotMatrix(ex$modules)[[4]] # exactly modular covariance matrix
  hyp <- ex$modules
  expect_warning(plain <- MINT(S, hyp), NA)
  res <- MINT(S, hyp, significance = TRUE, sample.size = 100, iterations = 50)$gamma_rank
  expect_equal(names(res), c(".id", ".n", "rank", "corrected.gamma", "raw.gamma", "sig", "p.value"))
  expect_equal(res[, 1:5], plain$gamma_rank[, 1:5])
  expect_true(all(res$sig >= 0 & res$sig <= 1))
  expect_true(all(res$p.value >= 1/51 & res$p.value <= 1))
  # the true hypothesis is consistent with the data (high p), wrong ones are not
  expect_equal(res$p.value[res$.id == "1_2_3"], 1)
  expect_equal(res$sig[res$.id == "1_2_3"], 0)
  expect_equal(res$p.value[res$.id == "null"], 1/51)
  expect_equal(res$sig[res$.id == "null"], 1)
  expect_error(MINT(S, hyp, significance = TRUE), "sample.size")
})

test_that("MINT gamma is the squared distance", {
  A <- matrix(c(1, .5, .5, 1), 2); B <- diag(2)
  expect_equal(.mint_gamma(A, B), 0.5)
  expect_equal(MINT(A, c(1, 1)), 0)
})
