#' Compare matrices via Mantel Correlation
#'
#' Calculates correlation matrix correlation and significance via Mantel test.
#'
#' @param cor.x Single correlation matrix or list of correlation matrices.
#'
#' If single matrix is supplied, it is compared to cor.y.
#'
#' If list is supplied and no cor.y is supplied, all matrices
#' are compared.
#'
#' If cor.y is supplied, all matrices in list are compared to it.
#' @param cor.y First argument is compared to cor.y.
#' Optional if cor.x is a list.
#' @param permutations Number of permutations used in significance calculation.
#' @param repeat.vector Vector of repeatabilities for correlation correction.
#' @param parallel if TRUE computations are done in parallel. Some foreach back-end must be registered, like doParallel or doMC.
#' @param ... additional arguments passed to other methods
#' @param landmark.dim Used if permutations should be performed maintaining landmark structure in geometric morphometric data. Either 2 for 2d data or 3 for 3d data. Default is NULL for non geometric morphomotric data.
#' @param withinLandmark Logical. If TRUE within-landmark correlations are also used in the calculation of matrix correlation (all lower triangle elements are used). Only used if landmark.dim is passed, default is FALSE.
#' @param mod Set TRUE to use mantel in testing modularity hypothesis. Should only be used in MantelModTest.
#' @return If cor.x and cor.y are passed, returns matrix Pearson
#' correlation coefficient and significance via Mantel permutations.
#'
#' If cor.x is a list of matrices and cor.y is passed, same as above, but for all matrices in cor.x.
#'
#' If only cor.x is passed, a matrix of MantelCor average
#' values and probabilities of all comparisons.
#' If repeat.vector is passed, comparison matrix is corrected above
#' diagonal and repeatabilities returned in diagonal.
#' @details MantelCor computes the Pearson correlation between the elements
#' of the lower triangles (excluding the diagonal) of the two matrices. Its
#' significance comes from a permutation test: the rows and columns of
#' \code{cor.x} are permuted together, the correlation with \code{cor.y} is
#' recomputed for each permutation, and these correlations form the null
#' distribution. The p-value (reported as \code{Probability}) is
#' (sum(null >= observed) + 1)/(permutations + 1), so it is never 0.
#' Values equal to the observed correlation up to floating point error count as
#' equal. If there are fewer than 5040 distinct permutations (7 traits or fewer),
#' or more permutations are requested than exist, all distinct permutations
#' (without the identity) are enumerated instead of sampled; if there are more of
#' them than \code{permutations}, a random subset of that size is used, and the
#' number of permutations in the formula above is the number actually used.
#'
#' If \code{landmark.dim} is given, the traits are grouped in landmarks of 2 or 3
#' consecutive coordinates, and the permutations shuffle whole landmarks
#' (keeping the coordinates of a landmark together); these permutations are
#' always sampled at random. By default (\code{withinLandmark = FALSE}) only
#' correlations between traits of different landmarks are used, both for the
#' observed correlation and for the null distribution. With
#' \code{withinLandmark = TRUE} all elements of the lower triangle are used,
#' including those between coordinates of the same landmark (the diagonal is
#' always excluded).
#'
#' If \code{cor.x} and \code{cor.y} are \code{mcmc_sample} objects, the i-th
#' sample of each is compared, and a vector with one value per sample is returned
#' (both are converted to correlations with \code{cov2cor}).
#' @note If the significance is not needed, MatrixCor provides the
#' correlation and skips the permutations, so it is much faster.
#' @export
#' @rdname MantelCor
#' @references Mantel, N. (1967). The detection of disease clustering and a generalized regression approach. Cancer Research 27(2): 209-220.
#' @author Diogo Melo, Guilherme Garcia
#' @seealso \code{\link{KrzCor}},\code{\link{RandomSkewers}},\code{\link{RandomSkewers}},\code{\link{TestModularity}}, \code{\link{MantelModTest}}
#' @examples
#' c1 <- RandomMatrix(10, 1, 1, 10)
#' c2 <- RandomMatrix(10, 1, 1, 10)
#' c3 <- RandomMatrix(10, 1, 1, 10)
#' MantelCor(cov2cor(c1), cov2cor(c2))
#'
#' cov.list <- list(c1, c2, c3)
#' cor.list <- lapply(list(c1, c2, c3), cov2cor)
#'
#' MantelCor(cor.list)
#'
#'# For repeatabilities we can use MatrixCor, which skips the significance calculation
#' reps <- unlist(lapply(cov.list, MonteCarloRep, 10, MatrixCor, correlation = TRUE))
#' MantelCor(cor.list, repeat.vector = reps)
#'
#' c4 <- RandomMatrix(10)
#' MantelCor(cor.list, c4)
#'
#' \dontrun{
#' #Multiple threads can be used with some foreach backend library, like doParallel or doMC
#' library(doParallel)
#' registerDoParallel(cores = 2)
#' MantelCor(cor.list, parallel = TRUE)
#' }
#' @keywords matrixcomparison
#' @keywords matrixcorrelation
#' @keywords randomskewers
MantelCor <- function (cor.x, cor.y, ...) UseMethod("MantelCor")

#' @rdname MantelCor
#' @method MantelCor default
#' @export
MantelCor.default <- function (cor.x, cor.y, permutations = 1000, ...,
                               landmark.dim = NULL, withinLandmark = FALSE, mod = FALSE) {
  if(!mod && !(.is_cor_matrix(cor.x) && .is_cor_matrix(cor.y)))
    warning("Matrices do not appear to be correlation matrices. Use with caution.")
  if(!is.null(landmark.dim)){
    if(!any(landmark.dim == c(2, 3))) stop("landmark.dim should be either 2 or 3 dimensions")
    if(withinLandmark) cor.mask <- lower.tri(cor.x)
    else cor.mask <- lower.tri.land(cor.x, landmark.dim)
    correlation <- cor(cor.x [cor.mask],
                       cor.y [cor.mask])
    null_vector <- vector("numeric", permutations)
    num.traits <- dim(cor.x)[1]
    n.land <- num.traits/landmark.dim
    for(perm in seq_len(permutations)){
      current.perm <- rep(sample.int(n.land) * landmark.dim, each = landmark.dim) - (landmark.dim-1):0
      null_vector[perm] <- cor(cor.x[cor.mask],
                               cor.y[current.perm, current.perm][cor.mask])
    }
    prob <- (sum(null_vector >= correlation - sqrt(.Machine$double.eps)) + 1)/(permutations + 1)
  } else{
    mantel.output <- .mantel_test(cor.x, cor.y, permutations = permutations)
    correlation <- mantel.output$statistic
    prob <- mantel.output$signif
  }
  output <- c(correlation, prob)
  names(output) <- c("rho", "Probability")
  return (output)
}

# Logical matrix, TRUE for pairs of traits that belong to the same landmark.
CreateWithinLandMat <- function(num.land, land.dim){
  kronecker(diag(num.land), matrix(1, land.dim, land.dim)) == 1
}

lower.tri.land <- function(x, landmark.dim = NULL){
  num.land <- dim(x)[1] / landmark.dim
  !CreateWithinLandMat(num.land, landmark.dim) & lower.tri(x)
}

#' @rdname MantelCor
#' @method MantelCor list
#' @export
MantelCor.list <- function (cor.x, cor.y = NULL,
                            permutations = 1000, repeat.vector = NULL,
                            parallel = FALSE, ...)
{
  if (is.null (cor.y)) {
    output <- ComparisonMap(cor.x,
                         function(x, cor.y) MantelCor(x, cor.y, permutations, ...),
                         repeat.vector = repeat.vector,
                         parallel = parallel)
  } else{
    output <- SingleComparisonMap(cor.x, cor.y,
                               function(x, y) MantelCor(y,
                                                        x,
                                                        permutations, ...),
                               parallel = parallel)
  }
  return(output)
}

#' @rdname MantelCor
#' @method MantelCor mcmc_sample
#' @export
MantelCor.mcmc_sample <- function (cor.x, cor.y, ..., parallel = FALSE)
{
  MatrixCor(cor.x, cor.y, parallel = parallel)
}

#' @export
#' @rdname MantelCor
MatrixCor <- function (cor.x, cor.y, ...) UseMethod("MatrixCor")

#' @rdname MantelCor
#' @method MatrixCor default
#' @export
MatrixCor.default <- function (cor.x, cor.y, ...)
{
  if(!(.is_cor_matrix(cor.x) && .is_cor_matrix(cor.y)))
    warning("Matrices do not appear to be correlation matrices. Use with caution.")
  c("correlations" = cor(cor.x[lower.tri(cor.x)], cor.y[lower.tri(cor.y)], ...))
}

#' @rdname MantelCor
#' @method MatrixCor list
#' @export
MatrixCor.list <- function (cor.x, cor.y = NULL,
                            permutations = 1000, repeat.vector = NULL,
                            parallel = FALSE, ...)
{
  if (length(cor.x) < 1) stop("cor.x must contain at least one matrix")
  if (!all(vapply(cor.x, .is_cor_matrix, logical(1))) ||
      (!is.null(cor.y) && !.is_cor_matrix(cor.y)))
    warning("Matrices do not appear to be correlation matrices. Use with caution.")
  # lower triangles stacked in columns: all correlations in a single call to cor()
  lower.triangles <- function(matrices)
    vapply(matrices, function(m) m[lower.tri(m)], numeric(nrow(matrices[[1]]) * (nrow(matrices[[1]]) - 1) / 2))
  if (is.null (cor.y)) {
    .check_same_dim(cor.x, cor.x[[1]])
    n.matrix <- length(cor.x)
    matrix.names <- .fix_names(names(cor.x), n.matrix)
    all.cors <- matrix(cor(matrix(lower.triangles(cor.x), ncol = n.matrix)), n.matrix, n.matrix)
    correlations <- matrix(0, n.matrix, n.matrix)
    correlations[upper.tri(correlations)] <- t(all.cors)[upper.tri(all.cors)]
    output <- .assemble_comparison(correlations, matrix(0, n.matrix, n.matrix),
                                   matrix.names, repeat.vector)[[1]]
  } else{
    .check_same_dim(cor.x, cor.y)
    cors <- cor(matrix(lower.triangles(cor.x), ncol = length(cor.x)),
                cor.y[lower.tri(cor.y)])
    output <- .list_to_df(lapply(cors, function(r) c(correlations = r)), names(cor.x))
  }
  return(output)
}

#' @rdname MantelCor
#' @method MatrixCor mcmc_sample
#' @export
MatrixCor.mcmc_sample <- function (cor.x, cor.y, ..., parallel = FALSE)
{
  if (inherits(cor.y, "mcmc_sample")) {
    n <- dim(cor.x)[1]
    if(dim(cor.y)[1] != n) stop("samples must be of same size")
    output <- .evq_lapply(seq_len(n), function(i) MatrixCor(cov2cor(cor.x[i,,]),
                                                            cov2cor(cor.y[i,,]))[["correlations"]],
                          parallel = parallel)
    output <- as.numeric(unlist(output))
  } else{
    output <- SingleComparisonMap(.mcmc_to_list(cor.x, cov2cor), cov2cor(cor.y),
                                  function(x, y) MatrixCor(x, y),
                                  parallel = parallel)
  }
  return(output)
}
