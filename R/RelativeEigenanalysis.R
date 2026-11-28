#' Relative Eigenanalysis
#'
#' Computes relative eigenvalues and eigenvectors between a pair of covariance matrices.
#'
#' @param cov.x covariance matrix
#' @param cov.y covariance matrix
#' @param symmetric Logical. If TRUE, computes symmetric eigenanalysis.
#' @return an object of class \code{eigen} (a list with \code{values} and \code{vectors}).
#' The values are always real, and ordered decreasingly.
#' @details With \code{symmetric = FALSE} the eigenvalues and eigenvectors of \code{cov.x \%*\% solve(cov.y)}
#' are returned; the eigenvectors are scaled to unit length. The computation uses a Cholesky whitening of
#' \code{cov.y} and symmetric eigenanalysis, so the eigenvalues are real and numerically stable, even for
#' many traits (the direct eigenanalysis of the non-symmetric product can return complex values).
#' With \code{symmetric = TRUE} the eigenanalysis is of \code{cov.y^(-1/2) \%*\% cov.x \%*\% cov.y^(-1/2)},
#' with unit-length eigenvectors in the whitened space. The sign of each eigenvector is arbitrary.
#'
#' @export
#'
#' 
#' @rdname RelativeEigenanalysis
#'
#' @references Bookstein, F. L., and P. Mitteroecker, P. "Comparing Covariance Matrices by
#' Relative Eigenanalysis, with Applications to Organismal Biology." Evolutionary Biology 41, no. 2
#' (June 1, 2014): 336-350. doi:10.1007/s11692-013-9260-5.
#'
#' @author Guilherme Garcia, Diogo Melo
#'
#' @examples
#' data(dentus)
#' dentus.vcv <- dlply(dentus, .(species), function(df) var(df[, -5]))
#'
#' dentus.eigrel <- RelativeEigenanalysis(dentus.vcv [[1]], dentus.vcv[[5]])
#' 
RelativeEigenanalysis <- function(cov.x, cov.y, symmetric = FALSE)
{
    if(symmetric)
    {
        out.elem <- .sym_inv_sqrtm(cov.y)
        matrix.prod <- out.elem %*% cov.x %*% out.elem
        matrix.prod <- (matrix.prod + t(matrix.prod)) / 2
        out <- eigen(matrix.prod, symmetric = TRUE)
    }
    else
    {
        ## cov.x %*% solve(cov.y) = L M L^-1, with cov.y = L L' and M = L^-1 cov.x L^-T symmetric,
        ## so eigenvectors are L u, u being the eigenvectors of M
        L <- t(chol(cov.y))
        M <- forwardsolve(L, t(forwardsolve(L, cov.x)))
        M <- (M + t(M)) / 2
        out <- eigen(M, symmetric = TRUE)
        vectors <- L %*% out$vectors
        out$vectors <- unname(vectors / rep(sqrt(colSums(vectors^2)), each = nrow(vectors)))
    }
    out
}
