#' Generic Comparison Map functions for creating parallel list methods
#' Internal functions for making eficient comparisons.
#' @param matrix.list list of matrices being compared
#' @param MatrixCompFunc Function used to compare pair of matrices, must output a vector: comparisons and probabilities
#' @param ... Aditional arguments to MatrixCompFunc
#' @param repeat.vector Vector of repeatabilities for correlation correction.
#' @param parallel if TRUE computations are done in parallel. Some foreach backend must be registered, like doParallel or doMC.
#' @return Matrix of comparisons, matrix of probabilities.
#' @details Unnamed lists get the names 1, 2, ...; empty names are filled
#' with their position and duplicated names are made unique with
#' \code{\link{make.unique}}. A list with a single matrix returns 1x1 matrices.
#' @author Diogo Melo
#' @seealso \code{\link{MantelCor}}, \code{\link{KrzCor}},\code{\link{RandomSkewers}}
ComparisonMap <- function (matrix.list, MatrixCompFunc, ..., repeat.vector = NULL, parallel = FALSE){
  n.matrix <- length(matrix.list)
  if (n.matrix < 1) stop("matrix.list must contain at least one matrix")
  matrix.names <- .fix_names(names(matrix.list), n.matrix)
  names(matrix.list) <- matrix.names
  pair.values <- .pairwise_values(matrix.list, MatrixCompFunc, ..., parallel = parallel)
  values <- do.call(rbind, lapply(pair.values$values, function(v) v[1:2]))
  probabilities <- array (0, c(n.matrix, n.matrix))
  correlations <- probabilities
  if (n.matrix > 1) {
    correlations[pair.values$pairs] <- values[, 1]
    probabilities[pair.values$pairs] <- values[, 2]
  }
  .assemble_comparison(correlations, probabilities, matrix.names, repeat.vector)
}

# Names for a list of matrices: NULL becomes 1, 2, ..., empty names are filled
# with their position, and duplicates are made unique.
.fix_names <- function(nms, n) {
  if (is.null(nms)) nms <- rep("", n)
  nms <- as.character(nms)
  empty <- is.na(nms) | nms == ""
  nms[empty] <- seq_len(n)[empty]
  make.unique(nms, sep = ".")
}

# Applies f(matrix.list[[j]], matrix.list[[i]]) to every pair i < j, one
# reference matrix i at a time (so random numbers are drawn in the same order
# as in earlier versions), and returns the list of
# results with the row/column indices (i, j) of each pair. Only the pairs are
# distributed to the parallel workers.
.pairwise_values <- function(matrix.list, f, ..., parallel = FALSE) {
  n.matrix <- length(matrix.list)
  pairs <- which(upper.tri(matrix(0, n.matrix, n.matrix)), arr.ind = TRUE)
  pairs <- pairs[order(pairs[, 1], pairs[, 2]), , drop = FALSE]
  dots <- list(...)
  one.pair <- function(p) {
    do.call(f, c(list(matrix.list[[pairs[p, 2]]], matrix.list[[pairs[p, 1]]]), dots))
  }
  values <- .evq_lapply(seq_len(nrow(pairs)), one.pair, parallel = parallel)
  list(pairs = pairs, values = values)
}

# Repeatability correction and final layout of ComparisonMap. `correlations`
# and `probabilities` hold the comparisons in their upper triangles.
.assemble_comparison <- function(correlations, probabilities, matrix.names,
                                 repeat.vector = NULL) {
  if (!is.null (repeat.vector)) {
    repeat.matrix <- sqrt(outer(repeat.vector, repeat.vector))
    correlations[lower.tri(correlations)] <- t(correlations/repeat.matrix)[lower.tri(correlations)]
    diag (correlations) <- repeat.vector
  }
  rownames (correlations) <- matrix.names
  colnames (correlations) <- matrix.names
  dimnames (probabilities) <- dimnames (correlations)
  list ('correlations' = t(correlations), 'probabilities' = t(probabilities))
}
