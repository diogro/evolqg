#' Compare matrices via Krzanowski Correlation
#'
#' Calculates covariance matrix correlation via Krzanowski Correlation
#'
#' @param cov.x Single covariance matrix or list of covariance matrices.
#' If single matrix is supplied, it is compared to cov.y.
#' If list is supplied and no cov.y is suplied, all matrices
#' are compared to each other.
#' If cov.y is supplied, all matrices in list are compared to it.
#' @param cov.y First argument is compared to cov.y.
#' Optional if cov.x is a list.
#' @param ret.dim number of retained dimensions in the comparison,
#' default for nxn matrix is n/2-1 (rounded, and at least 1)
#' @param repeat.vector Vector of repeatabilities for correlation correction.
#' @param parallel if TRUE and a list is passed, computations are done in parallel. Some foreach back-end must be registered, like doParallel or doMC.
#' @param ... additional arguments passed to other methods
#' @return If cov.x and cov.y are passed, returns Krzanowski correlation
#'
#' If cov.x is a list and cov.y is passed, same as above, but for all matrices in cov.x.
#'
#' If only a list is passed to cov.x, a matrix of Krzanowski correlation
#' values.
#' If repeat.vector is passed, comparison matrix is corrected above
#' diagonal and repeatabilities returned in diagonal.
#' @export
#' @rdname KrzCor
#' @references Krzanowski, W. J. (1979). Between-Groups Comparison of Principal
#' Components. Journal of the American Statistical Association, 74(367),
#' 703. doi:10.2307/2286995
#' @author Diogo Melo, Guilherme Garcia
#' @seealso \code{\link{RandomSkewers}},\code{\link{KrzProjection}},\code{\link{MantelCor}}
#' @examples
#' c1 <- RandomMatrix(10, 1, 1, 10)
#' c2 <- RandomMatrix(10, 1, 1, 10)
#' c3 <- RandomMatrix(10, 1, 1, 10)
#' KrzCor(c1, c2)
#'
#' KrzCor(list(c1, c2, c3))
#'\donttest{
#' reps <- unlist(lapply(list(c1, c2, c3), MonteCarloRep, 10, KrzCor, iterations = 10))
#' KrzCor(list(c1, c2, c3), repeat.vector = reps)
#'
#' c4 <- RandomMatrix(10)
#' KrzCor(list(c1, c2, c3), c4)
#' }
#'
#' \dontrun{
#' #Multiple threads can be used with some foreach backend library, like doParallel or doMC
#' library(doParallel)
#' registerDoParallel(cores = 2)
#' KrzCor(list(c1, c2, c3), parallel = TRUE)
#' }
#' @keywords matrixcomparison
#' @keywords matrixcorrelation
#' @keywords Krzanowski
KrzCor <- function (cov.x, cov.y, ...) UseMethod("KrzCor")

#' @rdname KrzCor
#' @method KrzCor default
#' @export
KrzCor.default <- function (cov.x, cov.y, ret.dim = NULL, ...) {
  if (is.null(ret.dim))
    ret.dim <- .default_ret_dim(nrow(cov.x))
  eVec.x <- eigen(cov.x, symmetric = TRUE)$vectors[, seq_len(ret.dim), drop = FALSE]
  eVec.y <- eigen(cov.y, symmetric = TRUE)$vectors[, seq_len(ret.dim), drop = FALSE]
  c("krz" = .krz_from_vectors(eVec.x, eVec.y))
}

# Default number of retained dimensions: n/2 - 1 rounded, but at least 1
# (the plain rule gives 0 for 2 and 3 traits).
.default_ret_dim <- function(n) max(1, round(n/2 - 1))

# Krzanowski correlation from matrices of retained eigenvectors (one per column).
.krz_from_vectors <- function(vec.x, vec.y)
  sum(crossprod(vec.x, vec.y)^2)/ncol(vec.x)

# First ret.dim eigenvectors of every matrix in a list, computed once.
.eigenvector_list <- function(matrices, ret.dim) {
  lapply(matrices, function(m) {
    if (is.null(ret.dim)) ret.dim <- .default_ret_dim(nrow(m))
    eigen(m, symmetric = TRUE)$vectors[, seq_len(ret.dim), drop = FALSE]
  })
}

#' @rdname KrzCor
#' @method KrzCor list
#' @export
KrzCor.list <- function (cov.x, cov.y = NULL,
                         ret.dim = NULL, repeat.vector = NULL,
                         parallel = FALSE, ...) {
  krz <- function(x, y) c("krz" = .krz_from_vectors(x, y), NA)
  if (is.null (cov.y)) {
    # eigenvectors are computed once per matrix, not once per pair
    output <- ComparisonMap(.eigenvector_list(cov.x, ret.dim), krz,
                            repeat.vector = repeat.vector,
                            parallel = parallel)
    output <- output[[1]]
  } else{
    .check_same_dim(cov.x, cov.y)
    output <- SingleComparisonMap(.eigenvector_list(cov.x, ret.dim),
                                  .eigenvector_list(list(cov.y), ret.dim)[[1]],
                                  krz, parallel = parallel, check.dim = FALSE)
    output <- output[,-length(output)]
  }
  return(output)
}

#' @rdname KrzCor
#' @method KrzCor mcmc_sample
#' @export
KrzCor.mcmc_sample <- function (cov.x, cov.y, ret.dim = NULL, parallel = FALSE, ...)
{
  if (inherits(cov.y, "mcmc_sample")) {
    n <- dim(cov.x)[1]
    if(dim(cov.y)[1] != n) stop("samples must be of same size")
    output <- .evq_lapply(seq_len(n), function(i) KrzCor(cov.x[i,,],
                                                         cov.y[i,,],
                                                         ret.dim = ret.dim)[["krz"]],
                          parallel = parallel)
    output <- as.numeric(unlist(output))
  } else{
    output <- SingleComparisonMap(.mcmc_to_list(cov.x), cov.y,
                                  function(x, y) KrzCor(x, y, ret.dim),
                                  parallel = parallel)
  }
  return(output)
}

# The samples of an mcmc_sample array as a list of matrices named "1", "2", ...
# (like plyr::alply(x, 1) did).
.mcmc_to_list <- function(x, f = identity) {
  out <- lapply(seq_len(dim(x)[1]), function(i) f(x[i,,]))
  names(out) <- as.character(seq_along(out))
  out
}
