#' Compare matrices using PCA similarity factor
#'
#' @param cov.x Single covariance matrix or list of covariance matrices.
#' If cov.x is a single matrix, it is compared to cov.y.
#' If cov.x is a list and no cov.y is supplied, all matrices
#' are compared to each other.
#' If cov.x is a list and cov.y is supplied, all matrices in cov.x are compared to cov.y.
#' @param ... additional arguments passed to other methods
#' @param cov.y First argument is compared to cov.y.
#' @param ret.dim number of retained dimensions in the comparison. Defaults to all.
#' @param repeat.vector Vector of repeatabilities for correlation correction.
#' @param parallel if TRUE computations are done in parallel. Some foreach back-end must be registered, like doParallel or doMC.
#' @return Ratio of projected variance to total variance
#' @references Singhal, A. and Seborg, D. E. (2005), Clustering multivariate time-series data. J. Chemometrics, 19: 427-438. doi: 10.1002/cem.945
#' @author Edgar Zanella Alvarenga
#' @seealso \code{\link{KrzProjection}},\code{\link{KrzCor}},\code{\link{RandomSkewers}},\code{\link{MantelCor}}
#' @rdname PCAsimilarity
#' @export
#' @examples
#' c1 <- RandomMatrix(10)
#' c2 <- RandomMatrix(10)
#' PCAsimilarity(c1, c2)
#'
#' m.list <- RandomMatrix(10, 3)
#' PCAsimilarity(m.list)
#'
#' PCAsimilarity(m.list, c1)
#' @keywords matrixcomparison
#' @keywords matrixcorrelation
#' @keywords Krzanowski
#' @keywords PCA

PCAsimilarity <- function(cov.x, cov.y, ...) UseMethod("PCAsimilarity")

#' @rdname PCAsimilarity
#' @method PCAsimilarity default
#' @export
PCAsimilarity.default <- function(cov.x, cov.y, ret.dim = NULL, ...) {
  if (is.null(ret.dim))
    ret.dim <- nrow(cov.x)
  c(PCAsimilarity = .pca_sim_from_eigen(.eigen_trimmed(cov.x, ret.dim),
                                        .eigen_trimmed(cov.y, ret.dim)))
}

# Eigenvalues and eigenvectors of the first ret.dim components.
.eigen_trimmed <- function(m, ret.dim) {
  eg <- eigen(m, symmetric = TRUE)
  keep <- seq_len(ret.dim)
  list(values = eg$values[keep], vectors = eg$vectors[, keep, drop = FALSE])
}

.pca_sim_from_eigen <- function(eg.x, eg.y)
  sum((eg.x$values %o% eg.y$values) * crossprod(eg.x$vectors, eg.y$vectors)^2) /
    sum(eg.x$values * eg.y$values)

#' @rdname PCAsimilarity
#' @method PCAsimilarity list
#' @export
PCAsimilarity.list <- function (cov.x, cov.y = NULL, ...,
                         repeat.vector = NULL, parallel = FALSE) {
  ret.dim <- list(...)[["ret.dim"]]
  # eigen decompositions are computed once per matrix, not once per pair
  eigen.list <- function(matrices)
    lapply(matrices, function(m)
      .eigen_trimmed(m, if (is.null(ret.dim)) nrow(m) else ret.dim))
  pca.sim <- function(x, y) c(PCAsimilarity = .pca_sim_from_eigen(x, y), NA)
  if (is.null (cov.y)) {
    output <- ComparisonMap(eigen.list(cov.x), pca.sim,
                         repeat.vector = repeat.vector,
                         parallel = parallel)
    output <- output[[1]]
  } else{
    .check_same_dim(cov.x, cov.y)
    output <- SingleComparisonMap(eigen.list(cov.x), eigen.list(list(cov.y))[[1]],
                                  pca.sim, parallel = parallel, check.dim = FALSE)
    output <- output[,-length(output)]
  }
  return(output)
}

#' @rdname PCAsimilarity
#' @method PCAsimilarity mcmc_sample
#' @export
PCAsimilarity.mcmc_sample <- function (cov.x, cov.y, ..., parallel = FALSE)
{
  if (inherits(cov.y, "mcmc_sample")) {
    n <- dim(cov.x)[1]
    if(dim(cov.y)[1] != n) stop("samples must be of same size")
    output <- .evq_lapply(seq_len(n), function(i) PCAsimilarity(cov.x[i,,],
                                                                cov.y[i,,], ...)[["PCAsimilarity"]],
                          parallel = parallel)
    output <- as.numeric(unlist(output))
  } else{
    output <- SingleComparisonMap(.mcmc_to_list(cov.x), cov.y,
                                  function(x, y) PCAsimilarity(x, y, ...),
                                  parallel = parallel)
  }
  return(output)
}
