#' Mean Covariance Matrix
#'
#' Estimate geometric mean for a set of covariance matrices
#'
#' @param matrix.array k x k x m array of covariance matrices, with k traits and m matrices
#' @param tol relative tolerance for the stopping rule: iteration stops when the spectral norm of the
#' mean of the log-mapped matrices (the Riemannian gradient) falls below \code{tol} times the Frobenius norm of the current mean
#' @details The mean is the Karcher (geometric, Riemannian) mean, computed by the Richardson-like
#' iteration of Bini and Iannazzo (2013). For two matrices it is the exact geometric mean
#' \eqn{A^{1/2}(A^{-1/2} B A^{-1/2})^{1/2} A^{1/2}}. A single matrix is returned as is.
#' @return geometric mean covariance matrix
#' 
#' @rdname MeanMatrix
#' @references Bini, D. A., Iannazzo, B. 2013. Computing the Karcher Mean of Symmetric
#' Positive Definite Matrices. Linear Algebra and Its Applications, 16th ILAS Conference
#' Proceedings, Pisa 2010, 438 (4): 1700-1710. doi:10.1016/j.laa.2011.08.052. 
#' @author Guilherme Garcia, Diogo Melo
#'
#' @examples
#' matrices <- simplify2array(RandomMatrix(5, 10, 1, 3))
#' MeanMatrix(matrices)
#' 
#' @seealso \code{\link{EigenTensorDecomposition}}, \code{\link{RiemannDist}}
#'
#' @export
#' 
MeanMatrix <- function (matrix.array, tol = 1e-10)
{
  A <- matrix.array
  k <- dim(A)[1]
  m <- dim(A)[3]
  dn <- dimnames(A)[1:2]
  if(is.null(dn) || all(vapply(dn, is.null, TRUE))) dn <- NULL
  if(m == 1)
    return(matrix(A[, , 1], k, k, dimnames = dn))

  A.chol <- lapply(seq_len(m), function(i) chol(A[, , i]))
  X <- rowMeans(A, dims = 2)
  X.next <- X
  niold <- Inf
  max.iter <- 1000L

  for(iter in seq_len(max.iter))
  {
    X.chol <- chol(X)
    X.chinv <- backsolve(X.chol, diag(k))
    S <- matrix(0, k, k)
    beta <- 0
    gamma <- 0
    for(i in seq_len(m))
    {
      Z <- A.chol[[i]] %*% X.chinv
      eig <- eigen(crossprod(Z), symmetric = TRUE)
      V <- eig$values
      S <- S + eig$vectors %*% (log(V) * t(eig$vectors))
      ## dynamic choice of theta
      ch <- max(V) / min(V)
      dh <- if(ch - 1 < 1e-8) 1 else log(ch) / (ch - 1)
      beta <- beta + dh
      gamma <- gamma + ch * dh
    }
    S <- (S + t(S)) / 2
    theta <- 2 / (gamma + beta)

    eig.S <- eigen(S, symmetric = TRUE)
    Z <- sqrt(exp(eig.S$values * theta)) * t(eig.S$vectors) %*% X.chol
    X.next <- crossprod(Z)
    ni <- max(abs(eig.S$values))

    if(ni > niold)
    {
      ## the gradient grew: keep the previous iterate, which was better
      X.next <- X
      break
    }
    if(ni < frobenius.norm(X.next) * tol)
      break
    niold <- ni
    X <- X.next
  }
  if(iter == max.iter)
    warning("MeanMatrix did not converge in ", max.iter, " iterations")
  X.next <- (X.next + t(X.next)) / 2
  if(!is.null(dn)) dimnames(X.next) <- dn
  X.next
}
