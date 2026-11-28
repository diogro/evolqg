# Internal helpers shared across evolqg. None of these are exported.
#
# Several of these replace functions that used to come from other packages
# (mvtnorm, MCMCpack, coda, vegan, expm). The replacements were written to
# give the same results, including the random number stream, so seeded
# analyses keep reproducing; tests/testthat/test.utils.r pins this against
# reference values computed with the original packages.

# Parallel-aware lapply --------------------------------------------------------

# Runs FUN over X. With parallel = TRUE it uses foreach %dopar% on whatever
# backend the user registered (doParallel, doMC, doFuture, ...), which is what
# plyr's .parallel = TRUE did. Without a registered backend foreach runs
# sequentially and warns, as before.
.evq_lapply <- function(X, FUN, ..., parallel = FALSE) {
  FUN <- match.fun(FUN)
  if (!isTRUE(parallel) || length(X) == 0L)
    return(lapply(X, FUN, ...))
  if (!requireNamespace("foreach", quietly = TRUE))
    stop("parallel = TRUE requires the 'foreach' package and a registered ",
         "parallel backend (e.g. doParallel::registerDoParallel()).", call. = FALSE)
  `%dopar%` <- foreach::`%dopar%`
  worker <- function(element) FUN(element, ...)
  element <- NULL # silence R CMD check note about the foreach iterator
  out <- foreach::foreach(element = X) %dopar% worker(element)
  names(out) <- names(X)
  out
}

# Matrix checks ----------------------------------------------------------------

# TRUE if all diagonal elements are 1 up to floating point error.
.is_cor_matrix <- function(x, tol = sqrt(.Machine$double.eps))
  isTRUE(all(abs(diag(as.matrix(x)) - 1) < tol))

# Mean squared correlation; the computation behind the deprecated CalcR2(),
# without its deprecation warning.
.r2 <- function(c.matrix) {
  cor.matrix <- cov2cor(c.matrix)
  mean(cor.matrix[lower.tri(cor.matrix)]^2)
}

# Functions of symmetric matrices ----------------------------------------------

# f applied to the eigenvalues of a symmetric matrix: V f(L) V'.
# Much faster than the general algorithms in expm for symmetric input.
.sym_matfun <- function(x, f) {
  x <- as.matrix(x)
  x <- (x + t(x)) / 2
  e <- eigen(x, symmetric = TRUE)
  out <- e$vectors %*% (f(e$values) * t(e$vectors))
  dimnames(out) <- dimnames(x)
  out
}
.sym_sqrtm <- function(x) .sym_matfun(x, sqrt)
.sym_inv_sqrtm <- function(x) .sym_matfun(x, function(l) 1 / sqrt(l))
.sym_logm <- function(x) .sym_matfun(x, log)
.sym_expm <- function(x) .sym_matfun(x, exp)

# Multivariate normal ----------------------------------------------------------

# Matrix root R with t(R) %*% R == sigma, computed the same way as
# mvtnorm::rmvnorm so that draws match it exactly for the same seed.
.mvn_root <- function(sigma, method = c("eigen", "svd", "chol")) {
  method <- match.arg(method)
  sigma <- as.matrix(sigma)
  if (!isSymmetric(sigma, tol = sqrt(.Machine$double.eps), check.attributes = FALSE))
    stop("sigma must be a symmetric matrix")
  psd.tol <- sqrt(.Machine$double.eps)
  not.psd <- function(values) any(values < -psd.tol * abs(values[1]))
  if (method == "eigen") {
    dec <- eigen(sigma, symmetric = TRUE)
    if (not.psd(dec$values)) warning("sigma is numerically not positive semidefinite")
    t(dec$vectors %*% (t(dec$vectors) * sqrt(pmax(dec$values, 0))))
  } else if (method == "svd") {
    dec <- svd(sigma)
    if (not.psd(dec$d)) warning("sigma is numerically not positive semidefinite")
    t(dec$v %*% (t(dec$u) * sqrt(pmax(dec$d, 0))))
  } else {
    root <- chol(sigma, pivot = TRUE)
    root[, order(attr(root, "pivot")), drop = FALSE]
  }
}

# n draws from N(mean, sigma), one per row. Pass a precomputed `root` (from
# .mvn_root) to avoid refactoring sigma on every call inside loops.
.rmvnorm <- function(n, mean = NULL, sigma, method = c("eigen", "svd", "chol"),
                     root = NULL) {
  if (is.null(root)) root <- .mvn_root(sigma, match.arg(method))
  p <- ncol(root)
  if (is.null(mean)) mean <- rep(0, p)
  if (length(mean) != p) stop("mean and sigma have non-conforming size")
  out <- matrix(rnorm(n * p), nrow = n, byrow = TRUE) %*% root
  out <- sweep(out, 2, mean, "+")
  colnames(out) <- names(mean)
  out
}

# Multivariate normal density for the rows of x.
.dmvnorm <- function(x, mean = NULL, sigma, log = FALSE) {
  if (is.vector(x)) x <- matrix(x, ncol = length(x))
  p <- ncol(x)
  if (is.null(mean)) mean <- rep(0, p)
  upper <- chol(sigma)
  z <- backsolve(upper, t(x) - mean, transpose = TRUE)
  log.dens <- -sum(log(diag(upper))) - 0.5 * p * log(2 * pi) - 0.5 * colSums(z^2)
  names(log.dens) <- rownames(x)
  if (log) log.dens else exp(log.dens)
}

# Inverse Wishart --------------------------------------------------------------

# One draw from the inverse Wishart with `df` degrees of freedom and scale S,
# through the Bartlett decomposition of the Wishart(df, S^-1). The order of
# the random draws (chi-squares on the diagonal, then normals filling the
# upper triangle column by column) matches MCMCpack::riwish, so seeded
# results are unchanged.
.riwish <- function(df, S) {
  S <- as.matrix(S)
  p <- nrow(S)
  if (df < p) stop("inverse Wishart degrees of freedom must be at least the dimension")
  scale.root <- chol(solve(S))
  bartlett <- matrix(0, p, p)
  diag(bartlett) <- sqrt(rchisq(p, df:(df - p + 1)))
  upper <- upper.tri(bartlett)
  bartlett[upper] <- rnorm(sum(upper))
  solve(crossprod(bartlett %*% scale.root))
}

# Highest posterior density intervals ------------------------------------------

# Shortest interval containing `prob` of the samples, for each column of x.
# Same result and layout as coda::HPDinterval(coda::as.mcmc(x)).
.hpd_interval <- function(x, prob = 0.95) {
  x <- as.matrix(x)
  if (is.null(colnames(x))) colnames(x) <- paste0("var", seq_len(ncol(x)))
  sorted <- apply(x, 2, sort)
  if (!is.matrix(sorted)) stop("need more than one sample to compute an HPD interval")
  n <- nrow(sorted)
  k <- ncol(sorted)
  gap <- max(1, min(n - 1, round(n * prob)))
  starts <- seq_len(n - gap)
  widths <- sorted[starts + gap, , drop = FALSE] - sorted[starts, , drop = FALSE]
  best <- apply(widths, 2, which.min)
  out <- cbind(sorted[cbind(best, seq_len(k))], sorted[cbind(best + gap, seq_len(k))])
  dimnames(out) <- list(colnames(x), c("lower", "upper"))
  attr(out, "Probability") <- gap / n
  out
}

# Mantel test ------------------------------------------------------------------

# All permutations of 1:n in lexicographic order, one per row.
.all_permutations <- function(n) {
  if (n == 1L) return(matrix(1L, 1, 1))
  smaller <- .all_permutations(n - 1L)
  do.call(rbind, lapply(seq_len(n), function(first)
    cbind(first, smaller + (smaller >= first), deparse.level = 0)))
}

# Permutation Mantel test of the lower triangles of two symmetric matrices,
# permuting rows and columns of x. Follows the choices vegan::mantel makes
# (through the permute package) so p-values match it for the same seed:
# when fewer than 5040 distinct permutations exist, or more permutations are
# requested than exist, the full set is enumerated (without the identity) and,
# if larger than requested, a random subset of rows is used; otherwise each
# permutation is drawn with sample.int().
#
# Supply `perm.matrix` (one permutation per row) to reuse the same
# permutations across several tests.
.mantel_test <- function(x, y, permutations = 999, perm.matrix = NULL) {
  x <- as.matrix(x)
  y <- as.matrix(y)
  n <- nrow(x)
  low <- lower.tri(x)
  y.low <- y[low]
  statistic <- cor(x[low], y.low)
  if (is.null(perm.matrix))
    perm.matrix <- .mantel_permutations(n, permutations)
  if (nrow(perm.matrix) == 0L)
    return(list(statistic = statistic, signif = NA_real_, perm = numeric(0),
                permutations = 0L))
  perm <- apply(perm.matrix, 1, function(p) cor(x[p, p][low], y.low))
  signif <- (sum(perm >= statistic - sqrt(.Machine$double.eps)) + 1) / (length(perm) + 1)
  list(statistic = statistic, signif = signif, perm = perm,
       permutations = length(perm))
}

.mantel_permutations <- function(n, permutations) {
  if (permutations < 1) return(matrix(integer(0), 0, n))
  n.distinct <- factorial(n) - 1 # the identity is not a permutation
  if (n.distinct < 5040.5 || permutations + 0.5 > n.distinct) {
    all.perms <- .all_permutations(n)[-1, , drop = FALSE]
    if (nrow(all.perms) > permutations)
      all.perms <- all.perms[sample.int(nrow(all.perms), permutations), , drop = FALSE]
    all.perms
  } else {
    t(vapply(seq_len(permutations), function(i) sample.int(n, n), integer(n)))
  }
}
