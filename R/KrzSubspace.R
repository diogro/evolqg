#' Krzanowski common subspaces analysis
#' 
#' Calculates the subspace most similar across a set of covariance matrices.
#' 
#' @param cov.matrices list of covariance matrices
#' @param k number of dimensions to be retained in calculating the subspace.
#' Must be between 1 and the number of traits. Default is \code{max(1, round(p/2 - 1))}, where p is the number of traits.
#' @return H shared space matrix
#' @return k_eVals_H eigen values for shared space matrix, maximum value for each is the number of matrices, representing a fully shared direction
#' @return k_eVec_H eigen vectors of shared space matrix (a p x k matrix)
#' @return k_angles data frame with the angles (in degrees) between each population subspace and each eigen vector of shared space matrix
#' @note can be used to implement the Bayesian comparison from Aguirre et al. 2014
#' @references Aguirre, J. D., E. Hine, K. McGuigan, and M. W. Blows. "Comparing G: multivariate analysis of genetic variation in multiple populations." Heredity 112, no. 1 (2014): 21-29.
#' @export
#' @examples 
#' data(dentus)
#' dentus.matrices = dlply(dentus, .(species), function(x) cov(x[-5]))
#' KrzSubspace(dentus.matrices, k = 2)
#' 
#' \dontrun{
#' # The method in Aguirre et al. 2014 can de implemented using this function as follows:
#' 
#' #Random input data with dimensions traits x traits x populations x MCMCsamples:
#' cov.matrices = aperm(aaply(1:10, 1, function(x) laply(RandomMatrix(6, 40,                    
#'                                                       variance = runif(6,1, 10)), 
#'                            identity)), 
#'                      c(3, 4, 1, 2))
#'    
#' Hs = alply(cov.matrices, 4, function(x) alply(x, 3)) |> llply(function(x) KrzSubspace(x, 3)$H)
#' avgH = Reduce("+", Hs)/length(Hs)
#' avgH.vec <- eigen(avgH)$vectors
#' MCMC.H.val = laply(Hs, function(mat) diag(t(avgH.vec) %*% mat %*% avgH.vec))
#' 
#' # confidence intervals for variation in shared subspace directions
#' library(coda)
#' HPDinterval(as.mcmc(MCMC.H.val))
#' }
KrzSubspace <- function(cov.matrices, k = NULL){
  p <- nrow(cov.matrices[[1]])
  if (is.null(k))
    k <- max(1, round(p/2 - 1))
  .check_krz_k(k, p)
  AA_T <- lapply(cov.matrices, function(x) {
    A <- eigen(x, symmetric = TRUE)$vectors[, seq_len(k), drop = FALSE]
    tcrossprod(A)
  })
  H <- Reduce("+", AA_T)
  eigen_H <- eigen(H, symmetric = TRUE)
  b <- eigen_H$vectors[, seq_len(k), drop = FALSE]
  angles <- t(vapply(AA_T, function(AA_Ti) {
    # clamp: rounding can put the cosine slightly outside [0, 1]
    acos(sqrt(pmin(1, pmax(0, colSums(b * (AA_Ti %*% b)))))) * (180/pi)
  }, numeric(k)))
  angles <- matrix(angles, nrow = length(AA_T), ncol = k)
  ids <- names(cov.matrices)
  if (is.null(ids)) ids <- paste0("m", seq_along(cov.matrices))
  angles <- data.frame(ids, angles, stringsAsFactors = FALSE)
  colnames(angles) <- c(".id", paste0("e", seq_len(k)))
  list(H = H,
       k_eVals_H = eigen_H$values[seq_len(k)],
       k_eVec_H = b,
       k_angles = angles)
}

.check_krz_k <- function(k, p) {
  if (!is.numeric(k) || length(k) != 1 || is.na(k) || k < 1 || k != round(k))
    stop("k must be a single positive integer.", call. = FALSE)
  if (k > p)
    stop("k (", k, ") cannot be larger than the number of traits (", p, ").", call. = FALSE)
  invisible(TRUE)
}

# Shared subspace matrix H only, without the angles computed by KrzSubspace().
.krz_H <- function(cov.matrices, k) {
  Reduce("+", lapply(cov.matrices, function(x)
    tcrossprod(eigen(x, symmetric = TRUE)$vectors[, seq_len(k), drop = FALSE])))
}
