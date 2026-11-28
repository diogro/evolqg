#' Compare matrices via RandomSkewers
#'
#' Calculates covariance matrix correlation via random skewers
#'
#' @param cov.x Single covariance matrix or list of covariance matrices.
#' If single matrix is supplied, it is compared to cov.y.
#' If list is supplied and no cov.y is supplied, all matrices
#' are compared.
#' If cov.y is supplied, all matrices in list are compared to it.
#' @param cov.y First argument is compared to cov.y.
#' Optional if cov.x is a list.
#' @param num.vectors Number of random vectors used in comparison.
#' @param repeat.vector Vector of repeatabilities for correlation correction.
#' @param parallel if TRUE computations are done in parallel. Some foreach back-end must be registered, like doParallel or doMC.
#' @param ... additional arguments passed to other methods.
#' @return
#' If cov.x and cov.y are passed, returns average value
#' of response vectors correlation ('correlation'), significance ('probability') and standard deviation
#' of response vectors correlation ('correlation_sd')
#'
#' If cov.x and cov.y are passed, same as above, but for all matrices in cov.x.
#'
#' If only a list is passed to cov.x, a matrix of RandomSkewers average
#' values and probabilities of all comparisons. All pairwise comparisons in a
#' list share the same random vectors, which are drawn once per call.
#' If repeat.vector is passed, comparison matrix is corrected above
#' diagonal and repeatabilities returned in diagonal.
#' @export
#' @useDynLib evolqg, .registration = TRUE
#' @rdname RandomSkewers
#' @references Cheverud, J. M., and Marroig, G. (2007). Comparing covariance matrices:
#' Random skewers method compared to the common principal components model.
#' Genetics and Molecular Biology, 30, 461-469.
#' @author Diogo Melo, Guilherme Garcia
#' @seealso \code{\link{KrzCor}},\code{\link{MantelCor}},\code{\link{DeltaZCorr}}
#' @examples
#' c1 <- RandomMatrix(10, 1, 1, 10)
#' c2 <- RandomMatrix(10, 1, 1, 10)
#' c3 <- RandomMatrix(10, 1, 1, 10)
#' RandomSkewers(c1, c2)
#'
#' RandomSkewers(list(c1, c2, c3))
#'\donttest{
#' reps <- unlist(lapply(list(c1, c2, c3), MonteCarloRep, sample.size = 10,
#'                                         RandomSkewers, num.vectors = 100,
#'                                         iterations = 10))
#' RandomSkewers(list(c1, c2, c3), repeat.vector = reps)
#'
#' c4 <- RandomMatrix(10)
#' RandomSkewers(list(c1, c2, c3), c4)
#' }
#' \dontrun{
#' #Multiple threads can be used with some foreach backend library, like doParallel or doMC
#' library(doParallel)
#' registerDoParallel(cores = 2)
#' RandomSkewers(list(c1, c2, c3), parallel = TRUE)
#' }
#'
#' @keywords matrixcomparison
#' @keywords matrixcorrelation
#' @keywords randomskewers
RandomSkewers <- function(cov.x, cov.y, ...) UseMethod("RandomSkewers")

#' @rdname RandomSkewers
#' @export
RandomSkewers.default <- function (cov.x, cov.y, num.vectors = 10000, ...) {
  output <- as.numeric(RS(cov.x, cov.y, num.vectors))
  names(output) <- c("correlation","probability","correlation_sd")
  return(output)
}

#' @rdname RandomSkewers
#' @method RandomSkewers list
#' @export
RandomSkewers.list <- function (cov.x, cov.y = NULL, num.vectors = 10000, repeat.vector = NULL, parallel = FALSE, ...)
{
  if (length(cov.x) < 1) stop("cov.x must contain at least one matrix")
  if (num.vectors < 1) stop("num.vectors must be positive")
  traits <- nrow(cov.x[[1]])
  if (is.null(cov.y)) {
    .check_same_dim(cov.x, cov.x[[1]])
  } else {
    .check_same_dim(cov.x, cov.y)
  }
  # The random vectors (and the null distribution of their correlations) are
  # drawn once and shared by all comparisons; the response of each matrix to
  # them is computed once, with columns normalised to unit length.
  skewers <- matrix(rnorm(traits * num.vectors), traits, num.vectors)
  base.vector <- rnorm(traits)
  null.dist <- .unit_columns(skewers)
  null.dist <- as.numeric(crossprod(base.vector/sqrt(sum(base.vector^2)), null.dist))
  responses <- function(matrices)
    lapply(matrices, function(m) .unit_columns(m %*% skewers))
  rs.pair <- function(x, y) {
    cosines <- colSums(x * y)
    avg <- mean(cosines)
    c(correlation = avg, probability = mean(null.dist > avg), correlation_sd = sd(cosines))
  }
  if (is.null (cov.y)) {
    output <- ComparisonMap(responses(cov.x), rs.pair,
                            repeat.vector = repeat.vector,
                            parallel = parallel)
  } else{
    output <- SingleComparisonMap(responses(cov.x), responses(list(cov.y))[[1]], rs.pair,
                                  parallel = parallel, check.dim = FALSE)
  }
  return(output)
}

# Scales each column to unit length.
.unit_columns <- function(x) x / rep(sqrt(colSums(x^2)), each = nrow(x))

#' @rdname RandomSkewers
#' @method RandomSkewers mcmc_sample
#' @export
RandomSkewers.mcmc_sample <- function (cov.x, cov.y, num.vectors = 10000, parallel = FALSE, ...)
{
  if (inherits(cov.y, "mcmc_sample")) {
    n <- dim(cov.x)[1]
    if(dim(cov.y)[1] != n) stop("samples must be of same size")
    output <- .evq_lapply(seq_len(n), function(i) RandomSkewers(cov.x[i,,],
                                                                cov.y[i,,],
                                                                num.vectors = num.vectors)[["correlation"]],
                          parallel = parallel)
    output <- as.numeric(unlist(output))
  } else{
    output <- SingleComparisonMap(.mcmc_to_list(cov.x), cov.y,
                                  function(x, y) RandomSkewers(x, y, num.vectors),
                                  parallel = parallel)
  }
  return(output)
}
