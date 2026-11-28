#' Generic Single Comparison Map functions for creating parallel list methods
#' Internal functions for making efficient comparisons.
#' @param matrix.list list of matrices being compared
#' @param y.mat single matrix being compared to list
#' @param MatrixCompFunc Function used to compare pair of matrices, must output a vector: comparisons and probabilities
#' @param ... Additional arguments to MatrixCompFunc
#' @param parallel if TRUE computations are done in parallel. Some foreach back-end must be registered, like doParallel or doMC.
#' @param check.dim if TRUE (default) check that every matrix in the list has the same dimensions as y.mat.
#' @return Data frame with one row per matrix in the list, and a \code{.id} column if the list is named.
#' @author Diogo Melo
#' @seealso \code{\link{MantelCor}}, \code{\link{KrzCor}},\code{\link{RandomSkewers}}
SingleComparisonMap  <- function(matrix.list, y.mat, MatrixCompFunc, ..., parallel = FALSE,
                                 check.dim = TRUE){
  if (check.dim) .check_same_dim(matrix.list, y.mat)
  dots <- list(...)
  one.matrix <- function(x) do.call(MatrixCompFunc, c(list(x, y.mat), dots))
  res <- .evq_lapply(matrix.list, one.matrix, parallel = parallel)
  .list_to_df(res, names(matrix.list))
}

.check_same_dim <- function(matrix.list, y.mat) {
  if (!is.matrix(y.mat))
    stop("Matrices on list and single matrice dimension do not match: single matrix is not a matrix")
  same <- vapply(matrix.list, function(x) {
    d <- dim(x)
    length(d) == 2L && all(d == dim(y.mat))
  }, logical(1))
  if (!all(same))
    stop("Matrices on list and single matrice dimension do not match")
  invisible(TRUE)
}

# Stacks a list of equal-length vectors into a data frame, one row per element,
# with a character .id column when ids are given. Column names come from the
# names of the first vector, empty ones become V<position>.
.list_to_df <- function(res, ids = NULL) {
  mat <- do.call(rbind, unname(res))
  if (is.null(mat)) mat <- matrix(numeric(0), 0, 0)
  nm <- colnames(mat)
  if (is.null(nm)) nm <- rep("", ncol(mat))
  empty <- is.na(nm) | nm == ""
  nm[empty] <- paste0("V", seq_len(ncol(mat)))[empty]
  colnames(mat) <- nm
  rownames(mat) <- NULL
  out <- as.data.frame(mat, stringsAsFactors = FALSE)
  if (!is.null(ids)) out <- data.frame(.id = as.character(ids), out, stringsAsFactors = FALSE,
                                       check.names = FALSE)
  out
}
