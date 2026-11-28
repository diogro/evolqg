#' Compare matrices via Modified Krzanowski Correlation
#'
#' Calculates the modified Krzanowski correlation between matrices,
#' projecting the variance in each principal  components of the first
#' matrix in to the ret.dim.2 components of the second matrix.
#'
#' @param cov.x Single covariance matrix ou list of covariance matrices.
#' If cov.x is a single matrix is supplied, it is compared to cov.y.
#' If cov.x is a list of matrices is supplied and no cov.y is supplied, all matrices
#' are compared between each other.
#' If cov.x is a list of matrices and a single cov.y matrix is supplied, all matrices in list are compared to it.
#' @param ... additional arguments passed to other methods
#' @param cov.y First argument is compared to cov.y.
#' If cov.x is a list, every element in cov.x is projected in cov.y.
#' @param ret.dim.1 number of retained dimensions for first matrix in comparison, default for nxn matrix is n/2-1 (rounded, and at least 1)
#' @param ret.dim.2 number of retained dimensions for second matrix in comparison, default for nxn matrix is n/2-1 (rounded, and at least 1)
#' @param parallel if TRUE computations are done in parallel. Some foreach back-end must be registered, like doParallel or doMC.
#' @param full.results if FALSE returns only total variance, if TRUE also per PC variance.
#' @return Ratio of projected variance to total variance, and ratio of projected total in each PC
#' @references Krzanowski, W. J. (1979). Between-Groups Comparison of Principal
#' Components. Journal of the American Statistical Association, 74(367),
#' 703. doi:10.2307/2286995
#' @author Diogo Melo, Guilherme Garcia
#' @seealso \code{\link{RandomSkewers}},\code{\link{MantelCor}}
#' @rdname KrzProjection
#' @export
#' @examples
#' c1 <- RandomMatrix(10)
#' c2 <- RandomMatrix(10)
#' KrzProjection(c1, c2)
#'
#'\donttest{
#' m.list <- RandomMatrix(10, 3)
#' KrzProjection(m.list)
#' KrzProjection(m.list, full.results = TRUE)
#' KrzProjection(m.list, ret.dim.1 = 5, ret.dim.2 = 4)
#' KrzProjection(m.list, ret.dim.1 = 4, ret.dim.2 = 5)
#'
#' KrzProjection(m.list, c1)
#' KrzProjection(m.list, c1, full.results = TRUE)
#' }
#' \dontrun{
#' #Multiple threads can be used with some foreach backend library, like doParallel or doMC
#' library(doParallel)
#' registerDoParallel(cores = 2)
#' KrzProjection(m.list, parallel = TRUE)
#' }
#' @keywords matrixcomparison
#' @keywords matrixcorrelation
#' @keywords Krzanowski

KrzProjection <- function(cov.x, cov.y, ...) UseMethod("KrzProjection")

#' @rdname KrzProjection
#' @method KrzProjection default
#' @export
KrzProjection.default <- function (cov.x, cov.y, ret.dim.1 = NULL, ret.dim.2 = NULL, ...) {
  num.traits <- dim(cov.x)[1]
  if (is.null(ret.dim.1))
    ret.dim.1 <- .default_ret_dim(num.traits)
  if (is.null(ret.dim.2))
    ret.dim.2 <- .default_ret_dim(num.traits)
  .krz_projection(eigen(cov.x, symmetric = TRUE), eigen(cov.y, symmetric = TRUE),
                  ret.dim.1, ret.dim.2)
}

# Projection of the variance in the first ret.dim.1 PCs of the first matrix onto
# the first ret.dim.2 PCs of the second one, from their eigen decompositions.
.krz_projection <- function(eigen.1, eigen.2, ret.dim.1, ret.dim.2) {
  keep.1 <- seq_len(ret.dim.1)
  eVal.1 <- eigen.1$values
  eVar.1 <- sweep(eigen.1$vectors[, keep.1, drop = FALSE], 2, sqrt(eVal.1[keep.1]), "*")
  ProjectionNorms <- colSums(crossprod(eigen.2$vectors[, seq_len(ret.dim.2), drop = FALSE],
                                       eVar.1)^2)
  names(ProjectionNorms) <- keep.1
  list(total.variation = sum(ProjectionNorms)/sum(eVal.1),
       per.PC = ProjectionNorms/eVal.1[keep.1])
}

#' @rdname KrzProjection
#' @method KrzProjection list
#' @export
KrzProjection.list <- function(cov.x, cov.y = NULL,
                               ret.dim.1 = NULL, ret.dim.2 = NULL,
                               parallel = FALSE, full.results = FALSE, ...){
  n.matrix <- length(cov.x)
  if (n.matrix < 1) stop("cov.x must contain at least one matrix")
  dims <- function(m) {
    n <- nrow(m)
    c(if (is.null(ret.dim.1)) .default_ret_dim(n) else ret.dim.1,
      if (is.null(ret.dim.2)) .default_ret_dim(n) else ret.dim.2)
  }
  # one eigen decomposition per matrix
  eigen.x <- lapply(cov.x, eigen, symmetric = TRUE)
  if(is.null(cov.y)){
    matrix.names <- .fix_names(names(cov.x), n.matrix)
    names(eigen.x) <- matrix.names
    # only the outer loop (one reference matrix at a time) is parallel
    CompareToNProj <- function(n) {
      k <- dims(cov.x[[n]])
      res <- lapply(eigen.x, function(x) .krz_projection(x, eigen.x[[n]], k[1], k[2]))
      if (full.results) res else vapply(res, function(r) r$total.variation, numeric(1))
    }
    comparisons.proj <- .evq_lapply(seq_len(n.matrix), CompareToNProj, parallel = parallel)
    if(full.results){
      names(comparisons.proj) <- matrix.names
      return(comparisons.proj)
    }
    # rows are the projected matrices, columns the matrices projected onto
    out <- do.call(cbind, comparisons.proj)
    dimnames(out) <- list(matrix.names, matrix.names)
    return(out)
  } else{
    .check_same_dim(cov.x, cov.y)
    eigen.y <- eigen(cov.y, symmetric = TRUE)
    k <- dims(cov.y)
    comparisons.proj <- .evq_lapply(eigen.x, function(x) .krz_projection(x, eigen.y, k[1], k[2]),
                                    parallel = parallel)
    if(full.results){
      return(comparisons.proj)
    } else{
      return(data.frame(V1 = vapply(comparisons.proj, function(r) r$total.variation, numeric(1)),
                        row.names = NULL))
    }
  }
}
