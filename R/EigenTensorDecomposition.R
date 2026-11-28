#' Eigentensor Decomposition 
#' 
#' This function performs eigentensor decomposition on a set of covariance matrices.
#'
#' @param matrices k x k x m array of m covariance matrices with k traits;
#' @param return.projection Should we project covariance matrices into estimated eigentensors? Defaults to TRUE
#' @param ... additional arguments for methods
#' @return List with the following components:
#' @return mean mean covariance matrices used to center the sample (obtained from \code{\link{MeanMatrix}})
#' @return mean.sqrt square root of mean matrix (saved for use in other functions, 
#' such as \code{\link{ProjectMatrix}} and \code{\link{RevertMatrix}})
#' @return values vector of ordered eigenvalues associated with eigentensors;
#' @return matrices array of eigentensor in matrix form;
#' @return projection matrix of unstandardized projected covariance matrices over eigentensors.
#' 
#' @details The number of estimated eigentensors is the minimum between the number of data points (m) minus one and 
#' the number of independent variables (k(k + 1)/2), in a similar manner to the usual 
#' principal component analysis.
#'
#' All computations are done in the affine-invariant (log-Euclidean) space around the Riemannian mean
#' of the matrices: each matrix A is mapped to \eqn{\log(M^{-1/2} A M^{-1/2})}, with M the mean from
#' \code{\link{MeanMatrix}}, and the covariance tensor Sigma of these log-centered matrices
#' is decomposed. The eigentensors are therefore orthonormal (Frobenius norm) symmetric matrices in
#' this space, the eigenvalues are the variances of the projections, and the projections are
#' uncorrelated PCA-like scores (they have mean zero, since the mean is the Karcher mean).
#' \code{\link{ProjectMatrix}} and \code{\link{RevertMatrix}} move between matrices and scores in
#' this same space. (Versions before 1.0 built Sigma from the raw matrix elements, but projected
#' log-centered matrices, so values, tensors and scores did not correspond to each other.)
#' 
#' @references Basser P. J., Pajevic S. 2007. Spectral decomposition of a 4th-order 
#' covariance tensor: Applications to diffusion tensor MRI. Signal Processing. 87:220-236.
#' @references Hine E., Chenoweth S. F., Rundle H. D., Blows M. W. 2009. Characterizing 
#' the evolution of genetic variance using genetic covariance tensors. Philosophical 
#' transactions of the Royal Society of London. Series B, Biological sciences. 364:1567-78.
#' 
#' @author Guilherme Garcia, Diogo Melo
#' @seealso \code{\link{ProjectMatrix}}, \code{\link{RevertMatrix}}
#' 
#' @examples 
#' data(dentus)
#'
#' dentus.vcv <- daply (dentus, .(species), function(x) cov(x[,-5]))
#'
#' dentus.vcv <- aperm(dentus.vcv, c(2, 3, 1))
#'
#' dentus.etd <- EigenTensorDecomposition(dentus.vcv, TRUE)
#'
#' # Plot some results
#' oldpar <- par(mfrow = c(1,2))  
#' plot(dentus.etd $ values, pch = 20, type = 'b', ylab = 'Eigenvalue')
#' plot(dentus.etd $ projection [, 1:2], pch = 20, 
#'      xlab = 'Eigentensor 1', ylab = 'Eigentensor 2')
#' text(dentus.etd $ projection [, 1:2],
#'      labels = rownames (dentus.etd $ projection), pos = 2)
#' par(oldpar)  
#' 
#' \donttest{
#' # we can also deal with posterior samples of covariance matrices using plyr
#' 
#' dentus.models <- dlply(dentus, .(species), 
#'                        lm, formula = cbind(humerus, ulna, femur, tibia) ~ 1)
#'
#' dentus.matrices <- llply(dentus.models, BayesianCalculateMatrix, samples = 100)
#'
#' dentus.post.vcv <- laply(dentus.matrices, function (L) L $ Ps)
#'
#' dentus.post.vcv <- aperm(dentus.post.vcv, c(3, 4, 1, 2))
#'
#' # this will perform one eigentensor decomposition for each set of posterior samples
#' dentus.post.etd <- alply(dentus.post.vcv, 4, EigenTensorDecomposition)
#' 
#' # which would allow us to observe the posterior 
#' # distribution of associated eigenvalues, for example
#' dentus.post.eval <- laply (dentus.post.etd, function (L) L $ values)
#' 
#' boxplot(dentus.post.eval, xlab = 'Index', ylab = 'Value', 
#'         main = 'Posterior Eigenvalue Distribution')
#' }
#'
#' @rdname EigenTensorDecomposition 
#' @export
#'
EigenTensorDecomposition <- function (matrices, return.projection = TRUE, ...) UseMethod('EigenTensorDecomposition')

#' @rdname EigenTensorDecomposition
#' @method EigenTensorDecomposition list
#' @export
EigenTensorDecomposition.list <- function (matrices, return.projection = TRUE, ...)
  {
    k <- nrow(matrices[[1]])
    list2array <- array(unlist(matrices, use.names = FALSE), c(k, k, length(matrices)),
                        dimnames = list(rownames(matrices[[1]]), colnames(matrices[[1]]),
                                        names(matrices)))
    EigenTensorDecomposition(list2array, return.projection, ...)
  }

#' @rdname EigenTensorDecomposition
#' @method EigenTensorDecomposition default
#' @export
EigenTensorDecomposition.default <- function (matrices, return.projection = TRUE, ...)
  {
    n.traits <- dim(matrices) [1]
    n.matrix <- dim(matrices) [3]
    mean.matrix <- MeanMatrix(matrices, ...)
    mean.sqrt <- .sym_sqrtm(mean.matrix)
    mean.is <- .sym_inv_sqrtm(mean.matrix)

    ### log-centered matrices, k x k x m
    lce.array <- array(0, dim(matrices))
    for (i in seq_len(n.matrix))
      lce.array[, , i] <- .sym_logm(mean.is %*% matrices[, , i] %*% mean.is)

    ### Build Sigma in the log-centered space. Each matrix is a vector with the
    ### diagonal followed by the lower triangle times sqrt(2), so that the
    ### Euclidean product of vectors is the Frobenius product of matrices.
    lower <- lower.tri(diag(n.traits))
    n.vars <- n.traits + sum(lower)
    vec.lce <- vapply(seq_len(n.matrix),
                      function (i) c(diag(matrix(lce.array[, , i], n.traits)),
                                     sqrt(2) * lce.array[, , i][lower]),
                      numeric(n.vars))
    vec.lce <- t(matrix(vec.lce, n.vars, n.matrix))
    Sigma <- var (vec.lce)

    n.tensor <- min (n.matrix - 1, dim (Sigma) [1])
    eigen.dec <- eigen (Sigma, symmetric = TRUE)
    et.vectors <- eigen.dec $ vectors [, seq_len(n.tensor), drop = FALSE]
    eigen.matrices <- array(0, c(n.traits, n.traits, n.tensor))
    for (j in seq_len(n.tensor))
      {
        eigen.mat <- matrix(0, n.traits, n.traits)
        eigen.mat [lower] <- et.vectors [-seq_len(n.traits), j] / sqrt (2)
        eigen.mat <- eigen.mat + t(eigen.mat)
        diag (eigen.mat) <- et.vectors [seq_len(n.traits), j]
        eigen.matrices [, , j] <- eigen.mat
      }
    dimnames (eigen.matrices) <-
      list (rownames (matrices),
            colnames (matrices),
            paste ('ET', seq_len(n.tensor), sep = ''))
    
    out <- list ('mean' = mean.matrix,
                 'mean.sqrt' = mean.sqrt,
                 'values' = eigen.dec $ values [seq_len(n.tensor)],
                 'matrices' = eigen.matrices)
   
    if (return.projection)
      {
        projection <- crossprod (matrix (lce.array, n.traits^2), matrix (eigen.matrices, n.traits^2))
        dn <- list (dimnames (matrices) [[3]], dimnames (eigen.matrices) [[3]])
        if (!is.null (names (dimnames (matrices))))
          names (dn) <- c (names (dimnames (matrices)) [3], "")
        dimnames (projection) <- dn
        out $ projection <- projection
      }
    return (out)
  }
