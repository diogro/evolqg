#' Parametric population samples with covariance or correlation matrices
#'
#' Using a multivariate normal model, random populations are generated
#' using the supplied covariance matrix. A statistic is calculated on the
#' random population and compared to the statistic calculated on the
#' original matrix. 
#'
#' @param cov.matrix Covariance matrix.
#' @param sample.size Size of the random populations
#' @param iterations Number of random populations
#' @param ComparisonFunc Comparison functions for the calculated statistic
#' @param StatFunc Function for calculating the statistic
#' @param parallel if TRUE computations are done in parallel. Some foreach back-end must be registered, like doParallel or doMC.
#' @details Since this function uses multivariate normal model to generate populations, only covariance matrices should be used.
#'
#' The random populations are simulated in the main process, in the same order whether or not
#' \code{parallel = TRUE}; only the statistics are computed in parallel.
#' Populations whose comparison fails or is NA are redrawn up to 100 times, then an error is raised;
#' a single warning reports that redraws happened.
#' @return a data frame with one row per random population: a factor column \code{X1} with the
#'   population number and the comparison values (named as the result of \code{ComparisonFunc}, or V1, V2, ...).
#' @author Diogo Melo, Guilherme Garcia
#' @seealso \code{\link{BootstrapRep}}, \code{\link{AlphaRep}}
#' @rdname MonteCarloStat
#' @export
#' @examples
#' cov.matrix <- RandomMatrix(5, 1, 1, 10)
#'
#' MonteCarloStat(cov.matrix, sample.size = 30, iterations = 50,
#'                ComparisonFunc = function(x, y) PCAsimilarity(x, y)[1],
#'                StatFunc = cov)
#'
#' #Calculating R2 confidence intervals
#' r2.dist <- MonteCarloR2(RandomMatrix(10, 1, 1, 10), 30)
#' quantile(r2.dist)
#' 
#' \dontrun{
#' #Multiple threads can be used with some foreach backend library, like doParallel
#' library(doParallel)
#' registerDoParallel(cores = 2)
#' MonteCarloStat(cov.matrix, sample.size = 30, iterations = 100,
#'                ComparisonFunc = function(x, y) KrzCor(x, y)[1],
#'                StatFunc = cov,
#'                parallel = TRUE)
#' }
#' @keywords parametricsampling
#' @keywords montecarlo
#' @keywords repeatability
MonteCarloStat <- function (cov.matrix, sample.size, iterations,
                            ComparisonFunc, StatFunc,
                            parallel = FALSE) {
  if(any(is.na(cov.matrix))) stop("covariance matrix has missing values.")
  if(!isSymmetric(cov.matrix)) stop("covariance matrix must be symmetric.")
  if(.is_cor_matrix(cov.matrix)) warning("Matrix appears to be a correlation matrix! Only covariance matrices should be used in parametric resampling.")
  # the matrix root is fixed, so it is computed once; the random draws are
  # the same as those of mvtnorm::rmvnorm(method = "svd")
  root <- .mvn_root(cov.matrix, "svd")
  # all populations are simulated here, in order, so that parallel and serial
  # runs use the same populations
  populations <- lapply(seq_len(iterations),
                        function(i) .rmvnorm(sample.size, root = root))
  worker <- function(x)
    .compute_with_retry(x,
                        function(pop) ComparisonFunc(cov.matrix, StatFunc(pop)),
                        function() .rmvnorm(sample.size, root = root))
  .collect_comparisons(.evq_lapply(populations, worker, parallel = parallel), "Monte Carlo")
}

#' Parametric repeatabilities with covariance or correlation matrices
#'
#' Using a multivariate normal model, random populations are generated
#' using the suplied covariance matrix. A statistic is calculated on the
#' random population and compared to the statistic calculated on the
#' original matrix. 
#'
#' @param cov.matrix Covariance matrix.
#' @param sample.size Size of the random populations.
#' @param ComparisonFunc comparison function.
#' @param ... Aditional arguments passed to ComparisonFunc.
#' @param iterations Number of random populations.
#' @param correlation If TRUE, correlation matrix is used, else covariance matrix. MantelCor and MatrixCor should always uses correlation matrix.
#' @param parallel If is TRUE and list is passed, computations are done in parallel. Some foreach backend must be registered, like doParallel or doMC.
#' @details Since this function uses multivariate normal model to generate populations, only covariance matrices should be used, even when computing repeatabilities for covariances matrices.
#' @return returns the mean repeatability, or mean value of comparisons from samples to original statistic.
#' @author Diogo Melo Guilherme Garcia
#' @seealso \code{\link{BootstrapRep}}, \code{\link{AlphaRep}}
#' @rdname MonteCarloRep
#' @export
#' @examples
#' cov.matrix <- RandomMatrix(5, 1, 1, 10)
#'
#' MonteCarloRep(cov.matrix, sample.size = 30, RandomSkewers, iterations = 20)
#' 
#'\donttest{
#' MonteCarloRep(cov.matrix, sample.size = 30, RandomSkewers, num.vectors = 100, 
#'               iterations = 20, correlation = TRUE)
#' MonteCarloRep(cov.matrix, sample.size = 30, MatrixCor, correlation = TRUE)
#' MonteCarloRep(cov.matrix, sample.size = 30, KrzCor, iterations = 20)
#' MonteCarloRep(cov.matrix, sample.size = 30, KrzCor, correlation = TRUE)
#'
#'
#' #Creating repeatability vector for a list of matrices
#' mat.list <- RandomMatrix(5, 3, 1, 10)
#' sapply(mat.list, MonteCarloRep, 30, KrzCor, correlation = TRUE)
#'}
#'
#' \dontrun{
#' #Multiple threads can be used with some foreach backend library, like doMC or doParallel
#' library(doParallel)
#' registerDoParallel(cores = 2)
#' MonteCarloRep(cov.matrix, 30, RandomSkewers, iterations = 100, parallel = TRUE)
#' }
#' @keywords parametricsampling
#' @keywords montecarlo
#' @keywords repeatability
MonteCarloRep <- function(cov.matrix,
                          sample.size,
                          ComparisonFunc,
                          ...,
                          iterations = 1000, 
                          correlation = FALSE, 
                          parallel = FALSE){
  if(correlation)  {StatFunc <- cov; c2v <- cov2cor
  } else {StatFunc <- cov; c2v <- function(x) x}
  repeatability <- MonteCarloStat(cov.matrix, sample.size, iterations,
                                  ComparisonFunc = function(x, y) ComparisonFunc(c2v(x), 
                                                                                 c2v(y), ...),
                                  StatFunc = StatFunc,
                                  parallel = parallel)
  return(mean(repeatability[,2]))
}

#' R2 confidence intervals by parametric sampling
#'
#' Using a multivariate normal model, random populations are generated
#' using the suplied covariance matrix. R2 is calculated on all the
#' random population, provinding a distribution based on the original matrix.
#'
#' @param cov.matrix Covariance matrix.
#' @param sample.size Size of the random populations
#' @param iterations Number of random populations
#' @param parallel if TRUE computations are done in parallel. Some foreach backend must be registered, like doParallel or doMC.
#' @details Since this function uses multivariate normal model to generate populations, only covariance matrices should be used.
#' @return returns a vector with the R2 for all populations
#' @author Diogo Melo Guilherme Garcia
#' @seealso \code{\link{BootstrapRep}}, \code{\link{AlphaRep}}
#' @export
#' @examples
#' r2.dist <- MonteCarloR2(RandomMatrix(10, 1, 1, 10), 30)
#' quantile(r2.dist)
#' @keywords parametricsampling
#' @keywords montecarlo
#' @keywords repeatability
MonteCarloR2 <- function (cov.matrix, sample.size, iterations = 1000, parallel = FALSE) {
  it.r2 <- MonteCarloStat(cov.matrix, sample.size, iterations,
                          ComparisonFunc = function(x, y) y,
                          StatFunc = function(x) .r2(cor(x)),
                          parallel = parallel)
  return (it.r2[,2])
}

