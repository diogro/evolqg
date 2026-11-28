test_that("KrzSubspace returns correct results",{
  local_seed(10)
  kr.subspace_aguirre_bench <- function(Gs, vec){
    if (dim(Gs)[[1]] != dim(Gs)[[2]]){
      stop("G array must be of order n x n x m x MCMCsamp")
    }
    if (is.na(dim(Gs)[4])) {
      stop("There are no MCMCsamples")
    }
    n <- dim(Gs)[[1]]
    m <- dim(Gs)[[3]]
    MCMCsamp <- dim(Gs)[[4]] 
    if(length(vec) != m){stop("vec must have length = m")}
    h <- function (g, v){
      AA <- array(, c(n, n, m))  
      for (k in 1:m){
        g.vec <- eigen(g[,,k])$vectors[,1:(v[k])] 
        AA[,,k] <- g.vec %*% t(g.vec)
      }
      H <- apply(AA, 1:2, sum)
      list(H = H, AA = AA)
    }
    #internal function to calculate AA and H
    MCMC.H <- array(, c(n, n, MCMCsamp))
    dimnames(MCMC.H) <- list(dimnames(Gs)[[1]], dimnames(Gs)[[1]], dimnames(Gs)[[4]])      
    MCMC.AA <- array(, c(n, n, m, MCMCsamp))
    dimnames(MCMC.AA) <- list(dimnames(Gs)[[1]], dimnames(Gs)[[1]], dimnames(Gs)[[3]], dimnames(Gs)[[4]])
    for (i in 1:MCMCsamp){
      kr <- h(Gs[,,,i], v = vec)
      MCMC.H[,,i] <- kr$H
      MCMC.AA[,,,i] <- kr$AA
    }	
    #calculate AA and H for the ith MCMC sample of the G array		
    avH <- apply(MCMC.H, 1:2, mean)
    rownames(avH) <- dimnames(Gs)[[1]]
    colnames(avH) <- dimnames(Gs)[[1]]
    #calculate the posterior mean H
    avAA <- apply(MCMC.AA, 1:3, mean)
    dimnames(avAA) <- list(dimnames(Gs)[[1]], dimnames(Gs)[[1]], dimnames(Gs)[[3]])
    #calculate the posterior mean AA
    avH.vec <- eigen(avH)$vectors
    #eigenanalysis of posterior mean H	
    proj<- function(a, b) t(b) %*% a %*% b
    #internal function to do projection
    avH.theta <- matrix(, n, m)
    for (i in 1:n){
      for (i in 1:n){
        avH.theta[i,] <- acos(sqrt(apply(avAA, 3, proj, b = avH.vec[,i]))) * (180/pi)
      }
    }
    #angles between the eigenvectors posterior mean H and the posterior mean subspaces of each population
    MCMC.H.val <- matrix(, MCMCsamp, n)
    colnames(MCMC.H.val) <- paste("h", 1:n, sep="")
    for (i in 1:n){
      MCMC.H.val[,i] <- apply(MCMC.H, 3, proj, b = avH.vec[,i])
    }
    #posterior distribution of the genetic variance for the eigenvectors of posterior mean H 
    MCMC.H.theta <- array(, c(n, m, MCMCsamp))
    rownames(MCMC.H.theta) <- paste("h", 1:n, sep="")
    colnames(MCMC.H.theta) <- dimnames(Gs)[[3]]
    for(i in 1:n){
      for(j in 1:MCMCsamp){
        MCMC.H.theta[i,,j] <- acos(sqrt(apply(MCMC.AA[,,,j], 3, proj, b = avH.vec[,i]))) * (180/pi)
      }
    }
    #posterior distribution of the angles between the eigenvectors of posterior mean H and the MCMC samples of the subspaces of each population
    list(avAA = avAA, avH = avH, MCMC.AA = MCMC.AA, MCMC.H = MCMC.H, MCMC.H.val = MCMC.H.val, MCMC.H.theta = MCMC.H.theta)
  }
  cov.matrices = aperm(aaply(1:10, 1, function(x) laply(RandomMatrix(6, 40, variance = runif(6,1, 10)), 
                                                        identity)), c(3, 4, 1, 2))
  Hs = llply(alply(cov.matrices, 4, function(x) alply(x, 3)), function(x) KrzSubspace(x, 3)$H)
  avgH = Reduce("+", Hs)/length(Hs)
  avgH.vec <- eigen(avgH)$vectors
  MCMC.H.val = laply(Hs, function(mat) diag(t(avgH.vec) %*% mat %*% avgH.vec))
  x = kr.subspace_aguirre_bench(cov.matrices, vec = rep(3, 10))
  expect_equal(.hpd_interval(x$MCMC.H.val), .hpd_interval(MCMC.H.val), ignore_attr = TRUE)
})

iris_lms <- function(cols = 1:4) {
  lapply(split(iris, iris$Species), function(data) lm(as.matrix(data[cols]) ~ 1))
}

test_that("KrzSubspaceBootstrap uses the requested number of dimensions", {
  lm_models <- iris_lms()
  local_seed(1)
  result_k2 <- KrzSubspaceBootstrap(lm_models, rep = 1, MCMCsamples = 2, k = 2)
  local_seed(1)
  result_k3 <- KrzSubspaceBootstrap(lm_models, rep = 1, MCMCsamples = 2, k = 3)

  expect_equal(dim(result_k2$observed), c(2, 4))
  expect_equal(dim(result_k2$random), c(2, 4))
  expect_false(isTRUE(all.equal(result_k2$observed, result_k3$observed)))
  expect_false(isTRUE(all.equal(result_k2$random, result_k3$random)))
})

test_that("KrzSubspaceBootstrap default k is 3, also with 4 traits", {
  lm_models <- iris_lms()
  local_seed(2)
  def <- KrzSubspaceBootstrap(lm_models, rep = 2, MCMCsamples = 5)
  local_seed(2)
  k3 <- KrzSubspaceBootstrap(lm_models, rep = 2, MCMCsamples = 5, k = 3)
  expect_equal(def, k3)
  expect_equal(dim(def$observed), c(5, 4))
  expect_equal(dim(def$random), c(10, 4))
  expect_error(KrzSubspaceBootstrap(lm_models, rep = 1, MCMCsamples = 2, k = 5), "larger than")
})

test_that("KrzSubspaceBootstrap observed values match a direct computation", {
  lm_models <- iris_lms()
  local_seed(3)
  out <- KrzSubspaceBootstrap(lm_models, rep = 1, MCMCsamples = 4, k = 2)
  local_seed(3)
  Ps <- lapply(lm_models, function(m) BayesianCalculateMatrix(m, samples = 4)$Ps)
  Hs <- lapply(1:4, function(s) KrzSubspace(lapply(Ps, function(P) P[s, , ]), 2)$H)
  avgH <- Reduce("+", Hs) / 4
  v <- eigen(avgH)$vectors
  ref <- t(sapply(Hs, function(H) diag(t(v) %*% H %*% v)))
  expect_equal(unname(out$observed), unname(ref))
})

test_that("KrzSubspaceBootstrap null is a permutation of the residuals with observed df", {
  lm_models <- lapply(split(iris, iris$Species), function(d) lm(as.matrix(d[1:4]) ~ Sepal.Length, data = transform(d, Sepal.Length = seq_len(nrow(d)))))
  seen <- list()
  testthat::local_mocked_bindings(
    BayesianCalculateMatrix = function(linear.m, samples = NULL, ...) {
      seen[[length(seen) + 1]] <<- linear.m
      p <- ncol(linear.m$residuals)
      list(Ps = array(rep(diag(p), each = samples), c(samples, p, p)) + 0)
    })
  local_seed(4)
  KrzSubspaceBootstrap(lm_models, rep = 1, MCMCsamples = 2, k = 2)
  null <- seen[-(1:3)]
  expect_length(null, 3)
  expect_equal(sapply(null, function(m) m$df.residual), sapply(lm_models, function(m) m$df.residual), ignore_attr = TRUE)
  pooled <- do.call(rbind, lapply(lm_models, residuals))
  used <- do.call(rbind, lapply(null, function(m) m$residuals))
  expect_equal(used[order(used[, 1]), ], pooled[order(pooled[, 1]), ], ignore_attr = TRUE)
})

test_that("KrzSubspaceBootstrap works in parallel with a registered backend", {
  skip_if_not_installed("doParallel")
  skip_on_cran()
  skip_on_os("windows")
  doParallel::registerDoParallel(cores = 2)
  on.exit(foreach::registerDoSEQ())
  out <- KrzSubspaceBootstrap(iris_lms(), rep = 2, MCMCsamples = 3, parallel = TRUE)
  expect_equal(dim(out$random), c(6, 4))
})

make_covs <- function(p, n = 3) {
  local_seed(20, .local_envir = parent.frame())
  lapply(seq_len(n), function(i) cov(matrix(rnorm(3 * p * p), 3 * p, p)))
}

test_that("KrzSubspace handles k = 1, small p and the default k", {
  for (p in c(2, 3, 4, 6)) {
    covs <- make_covs(p)
    out <- KrzSubspace(covs, 1)
    expect_equal(dim(out$k_eVec_H), c(p, 1))
    expect_equal(dim(out$k_angles), c(3, 2))
    def <- KrzSubspace(covs)
    expect_gte(ncol(def$k_eVec_H), 1)
    expect_equal(ncol(def$k_eVec_H), max(1, round(p / 2 - 1)))
    expect_true(all(is.finite(as.matrix(def$k_angles[-1]))))
  }
  expect_equal(names(KrzSubspace(make_covs(4), 2)$k_angles), c(".id", "e1", "e2"))
  expect_equal(KrzSubspace(make_covs(4), 2)$k_angles$.id, c("m1", "m2", "m3"))
  named <- setNames(make_covs(4), c("a", "b", "c"))
  expect_equal(KrzSubspace(named, 2)$k_angles$.id, c("a", "b", "c"))
})

test_that("KrzSubspace errors clearly when k > p", {
  expect_error(KrzSubspace(make_covs(4), 5), "cannot be larger")
  expect_error(KrzSubspace(make_covs(4), 0), "positive integer")
})

test_that("KrzSubspace gives zero angles for identical matrices", {
  m <- make_covs(5, 1)[[1]]
  out <- KrzSubspace(list(m, m, m), 3)
  expect_false(anyNA(out$k_angles))
  expect_equal(as.matrix(out$k_angles[-1]), matrix(0, 3, 3, dimnames = list(NULL, c("e1", "e2", "e3"))),
               tolerance = 1e-3, ignore_attr = TRUE)
  expect_equal(out$k_eVals_H, c(3, 3, 3))
})

test_that("KrzSubspaceDataFrame and PlotKrzSubspace work", {
  local_seed(5)
  out <- KrzSubspaceBootstrap(iris_lms(), rep = 2, MCMCsamples = 10, k = 2)
  df <- KrzSubspaceDataFrame(out)
  expect_equal(dim(df), c(8, 5))
  expect_equal(df$mean, (df$lower + df$upper) / 2)
  df1 <- KrzSubspaceDataFrame(out, n = 1)
  expect_equal(dim(df1), c(2, 5))
  expect_equal(levels(df1$PC), "PC1")
  plt <- PlotKrzSubspace(df)
  expect_s3_class(plt, "ggplot")
  expect_no_error(ggplot2::ggplot_build(plt))
})
