#' Calculate Mahalonabis distance for many vectors
#'
#' Calculates the Mahalanobis distance between a list of species mean, using a global covariance matrix
#' @param means list or array of species means being compared. array must have means in the rows.
#' @param cov.matrix a single covariance matrix defining the scale (or metric tensor) to be used in the distance calculation.
#' @param parallel Deprecated and ignored. Distances are computed in a single vectorised step.
#' @details The Mahalanobis distance between two means x and y is
#' d(x, y) = sqrt((x - y)' S^-1 (x - y)), where S is the covariance matrix
#' given in \code{cov.matrix} and S^-1 is its inverse. It is the euclidean
#' distance after rescaling the trait space by S. The function computes this
#' distance (not squared) for all pairs of means, using the same shared
#' covariance matrix S for every pair.
#' @return returns a matrix of species-species (non-squared) Mahalanobis distances. If there is a single mean, a 1 x 1 matrix with a zero is returned.
#' @author Diogo Melo
#' @export
#' @references Mahalanobis, P. C. (1936). On the generalised distance in statistics. Proceedings of the National Institute of Sciences of India 2: 49-55.
#' @seealso \code{\link{mahalanobis}}
#' @examples
#' mean.1 <- colMeans(matrix(rnorm(30*10), 30, 10))
#' mean.2 <- colMeans(matrix(rnorm(30*10), 30, 10))
#' mean.3 <- colMeans(matrix(rnorm(30*10), 30, 10))
#' mean.list <- list(mean.1, mean.2, mean.3)
#'
#' # If cov.matrix is the identity, calculated distance is euclidian:
#' euclidian <- MultiMahalanobis(mean.list, diag(10))
#' # Using a matrix with half the variance will give twice the distance between each mean:
#' half.euclidian  <- MultiMahalanobis(mean.list, diag(10)/2)
#'
#' # Other covariance matrices will give different distances, measured in the scale of the matrix
#' non.euclidian <- MultiMahalanobis(mean.list, RandomMatrix(10))
#'
#' #Input can be an array with means in each row
#' mean.array = array(1:36, c(9, 4))
#' mat = RandomMatrix(4)
#' MultiMahalanobis(mean.array, mat)
#'
MultiMahalanobis <- function (means, cov.matrix, parallel = FALSE) {
  if(is.data.frame(means) || (!is.array(means) && !is.list(means)))
    stop("means must be in a list or an array.")
  if(!isSymmetric(cov.matrix)) stop("covariance matrix must be symmetric.")
  if(is.array(means)){
    mean.names <- rownames(means)
    mean.array <- as.matrix(means)
  } else {
    mean.names <- names(means)
    mean.array <- do.call(rbind, means)
  }
  num.mean <- nrow(mean.array)
  if(is.null(mean.names)) mean.names <- as.character(seq_len(num.mean))
  # With U'U = S^-1, (x - y)' S^-1 (x - y) = |U (x - y)|^2, so Mahalanobis
  # distances are euclidean distances between the transformed means
  transformed <- mean.array %*% t(chol(solve(cov.matrix)))
  distances <- as.matrix(dist(unname(transformed)))
  dimnames(distances) <- list(mean.names, mean.names)
  return (distances)
}
