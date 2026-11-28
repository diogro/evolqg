#' Calculate mean values for various matrix statistics
#'
#' Calculates: Mean Squared Correlation, ICV, Autonomy, ConditionalEvolvability, Constraints, Evolvability, Flexibility, Pc1Percent, Respondability.
#' @aliases Autonomy ConditionalEvolvability Constraints Evolvability Flexibility Pc1Percent Respondability
#' @param cov.matrix A covariance matrix
#' @param beta.mat Matrix of direction vectors used to calculate the stochastic statistics. If NULL, a random matrix of normalized vectors is generated.
#' @param iterations Number of random vectors to be used in calculating the stochastic statistics when beta.mat is NULL.
#' @param full.results If TRUE, full distribution of statistics will be returned.
#' @param parallel if TRUE computations are done in parallel. Some foreach backend must be registered, like doParallel or doMC.
#' @param include.R2 Logical. If TRUE, include the legacy mean squared correlation metric in the returned output.
#' @param include.ICV Logical. If TRUE, include the legacy ICV metric in the returned output.
#' @details Autonomy and ConditionalEvolvability need the inverse of the covariance matrix. If the
#' matrix is not positive definite (or is numerically singular) a warning is issued and the nearest
#' positive definite matrix (\code{Matrix::nearPD}) is used instead, so the results are approximate.
#' @return dist Full distribution of stochastic statistics, only if full.resuts == TRUE
#' @return mean Mean value for all statistics
#' @export
#' @references Hansen, T. F., and Houle, D. (2008). Measuring and comparing evolvability
#' and constraint in multivariate characters. Journal of evolutionary
#' biology, 21(5), 1201-19. doi:10.1111/j.1420-9101.2008.01573.x
#' @author Diogo Melo Guilherme Garcia
#' @examples
#' cov.matrix <- cov(iris[,1:4])
#' MeanMatrixStatistics(cov.matrix)
#'
#' \dontrun{
#' #Multiple threads can be used with some foreach backend library, like doParallel
#' library(doParallel)
#' registerDoParallel(cores = 2)
#' MeanMatrixStatistics(cov.matrix, parallel = TRUE)
#' }
#' @keywords Autonomy
#' @keywords ConditionalEvolvability
#' @keywords Constraints
#' @keywords Evolvability
#' @keywords Flexibility
#' @keywords Pc1Percent
#' @keywords Respondability
MeanMatrixStatistics <- function (cov.matrix, iterations = 1000, full.results = FALSE, parallel = FALSE,
                                beta.mat = NULL, include.R2 = FALSE, include.ICV = FALSE) {
  matrix.stat.functions = list ('respondability' = Respondability,
                                'evolvability' = Evolvability,
                                'conditional.evolvability' = ConditionalEvolvability,
                                'autonomy' = Autonomy,
                                'flexibility' = Flexibility,
                                'constraints' = Constraints)
  num.traits <- dim (cov.matrix) [1]
  if (is.null(beta.mat))
    beta.mat <- .random_directions(num.traits, iterations)
  iso.vec <- Normalize (rep(1, num.traits))
  null.dist <- abs (t (iso.vec) %*% beta.mat)
  null.dist <- sort (null.dist)
  crit.value <- null.dist [round (0.95 * ncol(beta.mat))]
  if(full.results) cat ('critical value: ', crit.value, '\n')
  MatrixStatisticsMap <- function (CurrentFunc) CurrentFunc(cov.matrix = cov.matrix, beta.mat = beta.mat)
  stat.dist <- do.call(cbind, .evq_lapply(matrix.stat.functions, MatrixStatisticsMap, parallel = parallel))
  stat.dist <- cbind (stat.dist, null.dist)
  rownames (stat.dist) <- as.character(seq_len(nrow(stat.dist))) # as in earlier versions
  colnames (stat.dist) <- c('respondability',
                            'evolvability',
                            'conditional.evolvability',
                            'autonomy',
                            'flexibility',
                            'constraints',
                            'null.dist')
  stat.mean <- colMeans (stat.dist[,-7, drop = FALSE])
  integration <- c (Pc1Percent (cov.matrix),
                   CalcEigenVar(cov.matrix, sd=TRUE, rel=FALSE))
  names (integration) <- c ('pc1.percent', 'EigenSd')
  if (include.R2) {
    integration <- c (integration, "MeanSquaredCorrelation" = .r2(cov.matrix))
  }
  if (include.ICV) {
    integration <- c (integration, "ICV" = suppressWarnings(CalcICV(cov.matrix)))
  }
  stat.mean <- c (integration, stat.mean)
  if(full.results)
    return (list ('dist' = stat.dist, 'mean' = stat.mean))
  else
    return (stat.mean)
}

#' @export
#' @rdname MeanMatrixStatistics
#' @importFrom Matrix nearPD
Autonomy <- function (cov.matrix, beta.mat = NULL, iterations = 1000){
  num.traits <- dim (cov.matrix) [1]
  if(is.null(beta.mat))
    beta.mat <- .random_directions(num.traits, iterations)
  beta.mat <- as.matrix(beta.mat)
  dec <- .chol_or_nearpd(cov.matrix, "autonomy")
  inv.quad <- colSums(backsolve(dec$chol, beta.mat, transpose = TRUE)^2)
  (1 / inv.quad) / colSums(beta.mat * (dec$cov.matrix %*% beta.mat))
}
#' @export
#' @rdname MeanMatrixStatistics
#' @importFrom Matrix nearPD
ConditionalEvolvability <- function (cov.matrix, beta.mat = NULL, iterations = 1000){
  num.traits <- dim (cov.matrix) [1]
  if(is.null(beta.mat))
    beta.mat <- .random_directions(num.traits, iterations)
  beta.mat <- as.matrix(beta.mat)
  dec <- .chol_or_nearpd(cov.matrix, "conditional evolvability")
  1 / colSums(backsolve(dec$chol, beta.mat, transpose = TRUE)^2)
}
#' @export
#' @rdname MeanMatrixStatistics
Constraints <- function  (cov.matrix, beta.mat = NULL, iterations = 1000){
  num.traits <- dim (cov.matrix) [1]
  if(is.null(beta.mat))
    beta.mat <- .random_directions(num.traits, iterations)
  beta.mat <- as.matrix(beta.mat)
  PC1 <- eigen(cov.matrix, symmetric = TRUE)$vectors[,1]
  CB <- cov.matrix %*% beta.mat
  abs(as.vector(crossprod(PC1, CB)) / sqrt(colSums(CB^2)))
}
#' @export
#' @rdname MeanMatrixStatistics
Evolvability <- function (cov.matrix, beta.mat = NULL, iterations = 1000){
  num.traits <- dim (cov.matrix) [1]
  if(is.null(beta.mat))
    beta.mat <- .random_directions(num.traits, iterations)
  beta.mat <- as.matrix(beta.mat)
  colSums(beta.mat * (cov.matrix %*% beta.mat))
}
#' @export
#' @rdname MeanMatrixStatistics
Flexibility <- function (cov.matrix, beta.mat = NULL, iterations = 1000){
  num.traits <- dim (cov.matrix) [1]
  if(is.null(beta.mat))
    beta.mat <- .random_directions(num.traits, iterations)
  beta.mat <- as.matrix(beta.mat)
  CB <- cov.matrix %*% beta.mat
  colSums(beta.mat * CB) / sqrt(colSums(CB^2))
}
#' @export
#' @rdname MeanMatrixStatistics
Pc1Percent <- function (cov.matrix) {
  values <- eigen(cov.matrix, symmetric = TRUE, only.values = TRUE)$values
  values[1] / sum(values)
}
#' @export
#' @rdname MeanMatrixStatistics
Respondability <- function (cov.matrix, beta.mat = NULL, iterations = 1000) {
  num.traits <- dim (cov.matrix) [1]
  if(is.null(beta.mat))
    beta.mat <- .random_directions(num.traits, iterations)
  beta.mat <- as.matrix(beta.mat)
  sqrt(colSums((cov.matrix %*% beta.mat)^2))
}

# Random unit-length directions, one per column. Same random stream as
# normalising columns of rnorm(), and always a matrix, even for one trait.
.random_directions <- function(num.traits, iterations) {
  beta.mat <- matrix(rnorm(num.traits * iterations), num.traits, iterations)
  beta.mat / rep(sqrt(colSums(beta.mat^2)), each = num.traits)
}

# Cholesky factor of a covariance matrix, falling back to the nearest positive
# definite matrix when the matrix is singular (or numerically so).
.chol_or_nearpd <- function(cov.matrix, what) {
  cov.matrix <- as.matrix(cov.matrix)
  R <- tryCatch(chol(cov.matrix), error = function(cond) NULL)
  if (is.null(R) || min(diag(R)) <= sqrt(.Machine$double.eps) * max(diag(R))) {
    warning("matrix is singular, can't compute ", what,
            " directly. Using nearPD, results could be wrong")
    cov.matrix <- as.matrix(nearPD(cov.matrix)$mat)
    R <- chol(cov.matrix)
  }
  list(cov.matrix = cov.matrix, chol = R)
}

