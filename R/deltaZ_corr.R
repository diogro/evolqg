#' Compare matrices via the correlation between response vectors
#'
#' Compares the expected response to selection for two matrices for a specific set of
#' selection gradients (not random gradients like in the RandomSkewers method)
#'
#' @param cov.x Single covariance matrix or list of covariance matrices.
#' If single matrix is supplied, it is compared to cov.y.
#' If list is supplied and no cov.y is supplied, all matrices
#' are compared.
#' If cov.y is supplied, all matrices in list are compared to it.
#' @param cov.y First argument is compared to cov.y.
#' Optional if cov.x is a list.
#' @param skewers matrix of column vectors to be used as gradients (k columns).
#' A single vector is treated as a matrix with one column. Must have as many rows as the matrices have traits.
#' @param parallel if TRUE computations are done in parallel. Some foreach back-end must be registered, like doParallel or doMC.
#' @param ... additional arguments passed to other methods.
#' @return
#' For two matrices, a vector with the correlation between the expected responses of the two matrices
#' for each supplied vector, named with the column names of \code{skewers}.
#'
#' For a list and a single matrix \code{cov.y}, a data.frame with one row per matrix in the list
#' (and a \code{.id} column if the list is named) and one column per skewer.
#'
#' For a list only (all pairwise comparisons): with one skewer, a list with the matrices
#' \code{correlations} and \code{probabilities}, as in the other comparison methods
#' (the probabilities are not defined for this method and are \code{NA}). With more than one skewer,
#' an array of dimensions n x n x k, with dimnames (matrix names, matrix names, skewer names), where
#' the correlation between matrices i and j (j > i) for each skewer is in \code{[j, i, ]};
#' the diagonal and upper triangle are 0.
#' @export
#' @useDynLib evolqg, .registration = TRUE
#' @rdname DeltaZCorr
#' @references Cheverud, J. M., and Marroig, G. (2007). Comparing covariance matrices:
#' Random skewers method compared to the common principal components model.
#' Genetics and Molecular Biology, 30, 461-469.
#' @author Diogo Melo, Guilherme Garcia
#' @seealso \code{\link{KrzCor}},\code{\link{MantelCor}},\code{\link{RandomSkewers}}
#' @examples
#' x <- RandomMatrix(10, 1, 1, 10)
#' y <- RandomMatrix(10, 1, 1, 10)
#'
#' n_skewers = 10
#' skewers = matrix(rnorm(10*n_skewers), 10, n_skewers)
#' DeltaZCorr(x, y, skewers)
#'
#' # all pairs in a list, one 10 x 10 matrix of correlations per skewer
#' DeltaZCorr(list(x, y, RandomMatrix(10)), skewers = skewers)
#'
#' @keywords matrixcomparison
#' @keywords matrixcorrelation
#' @keywords randomskewers
DeltaZCorr <- function(cov.x, cov.y, skewers, ...) UseMethod("DeltaZCorr")

# Skewers as a matrix with one skewer per column, checked against the number of traits.
.check_skewers <- function(skewers, num.traits) {
  skewers <- as.matrix(skewers)
  if (!is.numeric(skewers) || ncol(skewers) < 1)
    stop("skewers must be a numeric matrix with at least one column")
  if (nrow(skewers) != num.traits)
    stop("skewers must have one row per trait (", num.traits, "), but has ", nrow(skewers))
  skewers
}

#' @rdname DeltaZCorr
#' @export
DeltaZCorr.default <- function (cov.x, cov.y, skewers, ...) {
  cov.x <- as.matrix(cov.x)
  cov.y <- as.matrix(cov.y)
  if (!identical(dim(cov.x), dim(cov.y)) || nrow(cov.x) != ncol(cov.x))
    stop("cov.x and cov.y must be square matrices of the same dimension")
  skewers <- .check_skewers(skewers, nrow(cov.x))
  output <- as.numeric(delta_z_corr(cov.x, cov.y, ncol(skewers), skewers))
  names(output) <- colnames(skewers)
  return(output)
}

#' @rdname DeltaZCorr
#' @method DeltaZCorr list
#' @export
DeltaZCorr.list <- function (cov.x, cov.y = NULL, skewers, parallel = FALSE, ...)
  {
    if (length(cov.x) < 1) stop("cov.x must contain at least one matrix")
    skewers <- .check_skewers(skewers, nrow(cov.x[[1]]))
    num.skewers <- ncol(skewers)
    if (is.null (cov.y)) {
      if (num.skewers == 1) {
        output <- ComparisonMap(cov.x,
                                function(x, y) c(DeltaZCorr(x, y, skewers), NA),
                                parallel = parallel)
      } else {
        n.matrix <- length(cov.x)
        matrix.names <- .fix_names(names(cov.x), n.matrix)
        names(cov.x) <- matrix.names
        pair.values <- .pairwise_values(cov.x, function(x, y) DeltaZCorr(x, y, skewers),
                                        parallel = parallel)
        skewer.names <- colnames(skewers)
        if (is.null(skewer.names)) skewer.names <- paste0("V", seq_len(num.skewers))
        output <- array(0, c(n.matrix, n.matrix, num.skewers),
                        dimnames = list(matrix.names, matrix.names, skewer.names))
        for (p in seq_along(pair.values$values))
          output[pair.values$pairs[p, 2], pair.values$pairs[p, 1], ] <- pair.values$values[[p]]
      }
    } else{
      .check_same_dim(cov.x, cov.y)
      output <- SingleComparisonMap(cov.x, cov.y,
                                    function(x, y) DeltaZCorr(x, y, skewers),
                                    parallel = parallel)
    }
    return(output)
  }
