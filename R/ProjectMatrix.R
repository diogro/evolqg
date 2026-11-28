#' Project Covariance Matrix
#'
#' This function projects a given covariance matrix into the basis provided by an eigentensor decomposition.
#'
#' @param matrix A symmetric covariance matrix for k traits
#' @param etd Eigentensor decomposition of m covariance matrices for k traits
#' (obtained from \code{\link{EigenTensorDecomposition}})
#' 
#' @return Vector of scores of given covariance matrix onto eigentensor basis. The scores are the
#' Frobenius products of the eigentensors with the log-centered matrix, as in the \code{projection}
#' of \code{\link{EigenTensorDecomposition}}, and \code{\link{RevertMatrix}} is its inverse.
#' 
#' @references Basser P. J., Pajevic S. 2007. Spectral decomposition of a 4th-order 
#' covariance tensor: Applications to diffusion tensor MRI. Signal Processing. 87:220-236.
#' @references Hine E., Chenoweth S. F., Rundle H. D., Blows M. W. 2009. Characterizing 
#' the evolution of genetic variance using genetic covariance tensors. Philosophical 
#' transactions of the Royal Society of London. Series B, Biological sciences. 364:1567-78.
#' 
#' @author Guilherme Garcia, Diogo Melo
#' @seealso \code{\link{EigenTensorDecomposition}}, \code{\link{RevertMatrix}}
#' 
#' @examples
#' # this function is useful for projecting posterior samples for a set of 
#' # covariance matrices onto the eigentensor decomposition done 
#' # on their estimated means
#' \donttest{
#' data(dentus)
#' 
#' dentus.models <- dlply(dentus, .(species), lm, 
#'                        formula = cbind(humerus, ulna, femur, tibia) ~ 1)
#'
#' dentus.matrices <- llply(dentus.models, BayesianCalculateMatrix, samples = 100)
#'
#' dentus.post.vcv <- laply(dentus.matrices, function (L) L $ Ps)
#' dentus.post.vcv <- aperm(dentus.post.vcv, c(3, 4, 1, 2))
#' 
#' dentus.mean.vcv <- aaply(dentus.post.vcv, 3, MeanMatrix)
#' dentus.mean.vcv <- aperm(dentus.mean.vcv, c(2, 3, 1))
#' 
#' dentus.mean.etd <- EigenTensorDecomposition(dentus.mean.vcv)
#' dentus.mean.proj <- data.frame('species' = LETTERS [1:5], dentus.mean.etd $ projection)
#' 
#' dentus.post.proj <- adply(dentus.post.vcv, c(3, 4), ProjectMatrix, etd = dentus.mean.etd)
#' colnames(dentus.post.proj) [1:2] <- c('species', 'sample')
#' levels(dentus.post.proj $ species) <- LETTERS[1:5]
#'
#' require(ggplot2)
#' ggplot() +
#'   geom_point(aes(x = ET1, y = ET2, color = species), 
#'      data = dentus.mean.proj, shape = '+', size = 8) +
#'   geom_point(aes(x = ET1, y = ET2, color = species), 
#'      data = dentus.post.proj, shape = '+', size = 3) +
#'   theme_bw()
#' }
#'
#' @export
#' 
ProjectMatrix <-
  function (matrix, etd)
  {
    mean.is <- .sym_matfun (etd $ mean.sqrt, function (l) 1 / l)
    log.center.mat <- .sym_logm (mean.is %*% matrix %*% mean.is)
    k <- nrow(log.center.mat)
    projection <- as.vector (crossprod (matrix (etd $ matrices, k * k), as.vector (log.center.mat)))
    names (projection) <- dimnames (etd $ matrices) [[3]]
    projection
  }
