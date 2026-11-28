#' Calculate Covariance Matrix from a linear model fitted with lm()
#'
#' Calculates the within-group (residual) covariance matrix of a multivariate linear model: the cross
#' product of the model residuals divided by the residual degrees of freedom. This is the unbiased
#' estimator, not the maximum likelihood estimator, which would divide by the number of observations.
#' For a univariate model a 1 x 1 matrix is returned.
#'
#' @details Covariance matrices of a sample are usually wanted after removing the effect of known
#' sources of variation, such as sex, population, age or size means, so that the matrix describes the
#' variation around those fixed effects and not the differences between their levels. Fitting these
#' effects with \code{lm()} (with a matrix of traits as response) and taking the covariance
#' matrix of the residuals does exactly that.
#'
#' The estimator used is \code{crossprod(residuals(linear.m)) / linear.m$df.residual}. The residuals
#' are not re-centered (a model without intercept is handled correctly), and the denominator is the
#' residual degrees of freedom, \eqn{n - r}, with n observations and r the number of estimated
#' coefficients (rank of the design matrix) for each trait: every estimated fixed effect uses up one
#' degree of freedom, which makes the estimator unbiased. For an intercept-only model
#' (\code{lm(y ~ 1)}) this is \eqn{n - 1}, and the result is the usual \code{cov(y)}.
#'
#' Use \code{cov2cor} on the result to obtain the corresponding correlation matrix. The residuals
#' (and so the matrix) do not depend on the contrasts used for factors, as long as the model includes
#' the same fixed effects.
#'
#' For more detailed models, the \pkg{stanimal} package (\url{https://github.com/diogro/stanimal}, not on
#' CRAN; install with \code{pak::pkg_install("diogro/stanimal/package")}, it needs CmdStan through
#' \pkg{cmdstanr}) has \code{stanimal::lm_phen_cov()}, which fits a fully Bayesian phenotypic covariance
#' matrix after fixed effects in Stan, with an LKJ prior on the correlations (or a factor-structured
#' covariance), missing trait values allowed, and posterior draws of P, corrP and the fixed effects.
#' Use it when a full posterior, regularisation for many traits, or incomplete data matter.
#' \code{stanimal::lm_phen_cov_laplace()} is a fast Laplace approximation of the same model, and
#' \code{stanimal::lmm_animal()} estimates the additive genetic matrix G (and the residual E) from a
#' pedigree or genomic relationship matrix, for when G and not P is the target.
#' @param linear.m Linear model adjusted for original data, fitted with \code{lm()}.
#'
#' @return Estimated covariance matrix.
#' @author Diogo Melo, Fabio Machado
#' @export
#' @examples
#' data(iris)
#' iris.lm = lm(as.matrix(iris[,1:4])~iris[,5])
#' cov.matrix <- CalculateMatrix(iris.lm)
#' 
#' #To obtain a corrlation matrix, use:
#' cor.matrix <- cov2cor(cov.matrix)
#' @keywords covariancematrix
CalculateMatrix <- function(linear.m){
  residuals <- as.matrix(residuals(linear.m))
  cov.matrix <- crossprod(residuals) / linear.m$df.residual
  return (cov.matrix)
}

#' Calculate Covariance Matrix from a linear model fitted with lm() using different estimators
#'
#' Calculates the residual covariance matrix of a multivariate linear model with a regularized
#' (shrinkage) Wishart prior, and, if requested, samples from the resulting posterior distribution
#' and their element-wise median.
#'
#' @details Let \eqn{S_x} be the cross product of the residuals of the model (\code{crossprod(residuals)}),
#' \eqn{N} the residual degrees of freedom (as in \code{\link{CalculateMatrix}}) and p the number of traits.
#' The prior on the covariance matrix is an inverse Wishart with \eqn{\nu} degrees of freedom (\code{nu})
#' and scale matrix \eqn{S_0} (\code{S_0}). By default \eqn{\nu = p} and \eqn{S_0} is a diagonal matrix with
#' the observed variances times \eqn{\nu} on the diagonal and zero covariances, so the prior pulls the
#' covariances towards zero and keeps the variances at their observed values. A large \code{nu}
#' is a strong prior, and with the default \eqn{\nu = p} the prior carries as much weight as p extra observations,
#' which can be too strong for small samples.
#' \itemize{
#'  \item \code{MLE} is \eqn{S_x / N}, the residual cross product divided by the residual degrees of
#'  freedom. Despite the name this is not the maximum likelihood estimator (which would divide by the
#'  number of observations): it is the unbiased estimator, equal to \code{\link{CalculateMatrix}}.
#'  \item \code{MAP} is \eqn{(S_0 + S_x) / (\nu + N)}, the sample estimate shrunk towards the prior, a
#'  weighted average of the prior scale and the \code{MLE}. It is the estimate to use when the sample
#'  is small. Note that this is not the mode of the posterior distribution that is sampled
#'  below, which for k traits is \eqn{(S_0 + S_x) / (\nu + N + k + 1)}.
#'  \item \code{Ps} (only when \code{samples} is given) are \code{samples} draws from the posterior
#'  distribution, which is an inverse Wishart with \eqn{\nu + N} degrees of freedom and scale
#'  \eqn{S_0 + S_x}, returned as a samples x p x p array. The draws are independent and use the random
#'  number stream, so set a seed for reproducibility. They can be used to propagate the uncertainty of
#'  the matrix estimate to any downstream statistic (for example, one eigentensor decomposition
#'  per sample).
#'  \item \code{P} (only when \code{samples} is given) is the element-wise median of the posterior
#'  samples, a point estimate of the matrix.
#' }
#' Sampling from the posterior needs \eqn{\nu + N \geq p} (an error otherwise), and a
#' non-singular \eqn{S_0 + S_x}. With the default prior (\eqn{\nu = p}, diagonal \eqn{S_0}) the first holds
#' always and the second holds whenever every trait has some variance, so "large enough" only matters
#' for the quality of the estimates: with fewer observations than traits \eqn{S_x} is singular and the
#' posterior is dominated by the prior. With a custom, low-rank \code{S_0} or a small \code{nu}, a sample
#' at least as large as the number of traits is needed. \code{MAP} and \code{MLE} are always computed.
#'
#' For more detailed models, the \pkg{stanimal} package (\url{https://github.com/diogro/stanimal}, not on
#' CRAN; install with \code{pak::pkg_install("diogro/stanimal/package")}, it needs CmdStan through
#' \pkg{cmdstanr}) has \code{stanimal::lm_phen_cov()}, a fully Bayesian phenotypic covariance matrix
#' after fixed effects, fitted in Stan with an LKJ prior on the correlations (or a factor-structured
#' covariance), missing trait values allowed, and posterior draws of P, corrP and the fixed effects.
#' By contrast, the posterior here is the conjugate inverse Wishart for the residual covariance: quick,
#' but with a less flexible prior. Use \code{stanimal::lm_phen_cov()} when a full posterior, regularisation
#' for many traits, or incomplete data matter; \code{stanimal::lm_phen_cov_laplace()} is a fast Laplace
#' approximation of the same model, and \code{stanimal::lmm_animal()} estimates the additive genetic matrix
#' G (and the residual E) from a pedigree or genomic relationship matrix, for when G is the target.
#' @param linear.m Linear model adjusted for original data
#' @param samples number of samples to be generated from the posterior. See Details for the sample size needed
#' @param ... additional arguments, currently ignored
#' @param nu degrees of freedom in prior distribution, defaults to the number of traits (this can be a too strong prior)
#' @param S_0 cross product matrix of the prior. Default is to use the observed variances and zero covariance
#'
#' @return Estimated covariance matrices and posterior samples
#' @author Diogo Melo, Fabio Machado
#' @references Murphy, K. P. (2012). Machine learning: a probabilistic perspective. MIT press.
#' @references Schafer, J., e Strimmer, K. (2005). A shrinkage approach to large-scale covariance matrix estimation and implications for functional genomics. Statistical applications in genetics and molecular biology, 4(1).
#' @export
#' @examples
#' data(iris)
#' iris.lm = lm(as.matrix(iris[,1:4])~iris[,5])
#' matrices <- BayesianCalculateMatrix(iris.lm, nu = 0.1, samples = 100)
#'
#' @keywords covariancematrix
#' @importFrom stats median residuals
BayesianCalculateMatrix <- function(linear.m, samples = NULL, ..., nu = NULL, S_0 = NULL){
  x <- as.matrix(residuals(linear.m))
  S_x = t(x) %*% x
  N = linear.m$df.residual
  p = dim(x)[2]
  if(is.null(nu)) nu = p
  if(is.null(S_0)){
    S_0 = diag(diag(S_x/N) * nu, nrow = p)
  }
  S_N <- S_0 + S_x
  nu_N <- nu + N
  MAP <- (S_0 + S_x) / (nu + N)
  MLE <- S_x/N
  if(!is.null(samples)){
   S_sample <- array(0, c(samples, p, p),
                     dimnames = list(NULL, rownames(S_N), colnames(S_N)))
   for(i in seq_len(samples)) S_sample[i, , ] <- .riwish(nu_N, S_N)
   median.P <- .median_matrix(S_sample)
   class(S_sample) <- "mcmc_sample"
   return(list(MAP = MAP, 
               MLE = MLE, 
               P = median.P, 
               Ps = S_sample))
  }
  else return(list(MAP = MAP, MLE = MLE))
}

# Element-wise median over the first dimension of a samples x k x k array.
# Dimnames are named X1 and X2, as they always were in the output.
.median_matrix <- function(x) {
  out <- apply(x, 2:3, median)
  k <- dim(x)[2:3]
  dn <- dimnames(x)
  if (is.null(dn)) dn <- list(NULL, NULL, NULL)
  dimnames(out) <- list(X1 = if (is.null(dn[[2]])) as.character(seq_len(k[1])) else dn[[2]],
                        X2 = if (is.null(dn[[3]])) as.character(seq_len(k[2])) else dn[[3]])
  out
}

#' @export
#' @method print mcmc_sample
#' @aliases CalculateMatrix
print.mcmc_sample <- function(x, ...) {
  cat(paste0("MCMC sample of ", dim(x)[1], 
             " matrices of dimension ", dim(x)[2], 
             "x", dim(x)[3]), "\n")
  cat("Median matrix:\n")
  print(.median_matrix(unclass(x)))
}
