#' Bootstrap analysis via resampling
#'
#'
#'   Calculates the repeatability of the covariance matrix of the supplied data
#'   via bootstrap resampling
#'
#'   Samples with replacement are taken from the full population, a statistic calculated
#'   and compared to the full population statistic. 
#'
#' @param ind.data Matrix of residuals or individual measurements
#' @param ComparisonFunc comparison function
#' @param iterations Number of resamples to take
#' @param sample.size Size of resamples, default is the same size as ind.data 
#' @param correlation If TRUE, correlation matrix is used, else covariance matrix. 
#' @param parallel if TRUE computations are done in parallel. Some foreach backend must be registered, like doParallel or doMC.
#' @param ... Additional arguments passed to ComparisonFunc.
#' @return returns the mean repeatability, that is, the mean value of comparisons from samples to original statistic.
#' @details Resamples are drawn in the main process, in the same order whether or not
#'   \code{parallel = TRUE}, so the resampled populations (and the result, for deterministic
#'   comparison functions) do not depend on the number of workers.
#'   Resamples for which the comparison fails or returns NA are redrawn, up to 100 times per resample,
#'   after which an error is raised; a single warning reports that redraws happened.
#' @author Diogo Melo, Guilherme Garcia
#' @seealso \code{\link{MonteCarloStat}}, \code{\link{AlphaRep}}
#' @export
#' @examples
#' BootstrapRep(iris[,1:4], MantelCor, iterations = 5, correlation = TRUE)
#'              
#' BootstrapRep(iris[,1:4], RandomSkewers, iterations = 50)
#' 
#' BootstrapRep(iris[,1:4], KrzCor, iterations = 50, correlation = TRUE)
#' 
#' BootstrapRep(iris[,1:4], PCAsimilarity, iterations = 50)
#' 
#' #Multiple threads can be used with some foreach backend library, like doMC or doParallel
#' #library(doParallel)
#' ##Windows:
#' #cl <- makeCluster(2)
#' #registerDoParallel(cl)
#' ##Mac and Linux:
#' #registerDoParallel(cores = 2)
#' #BootstrapRep(iris[,1:4], PCAsimilarity,
#' #             iterations = 5,
#' #             parallel = TRUE)
#' @keywords bootstrap
#' @keywords repetabilities


BootstrapRep <- function(ind.data,
                         ComparisonFunc,
                         iterations = 1000, 
                         sample.size = dim (ind.data)[1],
                         correlation = FALSE, 
                         parallel = FALSE,
                         ...){
  if(correlation)  {StatFunc <- cor; c2v <- cov2cor
  } else {StatFunc <- cov; c2v <- function(x) x}
  repeatability <- BootstrapStat(ind.data, iterations,
                                 ComparisonFunc = function(x, y) ComparisonFunc(c2v(x), 
                                                                                c2v(y), ...),
                                 StatFunc = StatFunc,
                                 sample.size = sample.size,
                                 parallel = parallel)
  return(mean(repeatability[,2]))
}


#' Non-Parametric population samples and statistic comparison
#'
#' Random populations are generated via ressampling 
#' using the suplied population. A statistic is calculated on the
#' random population and compared to the statistic calculated on the
#' original population. 
#'
#' @param ind.data Matrix of residuals or indiviual measurments
#' @param iterations Number of resamples to take
#' @param ComparisonFunc comparison function
#' @param StatFunc Function for calculating the statistic
#' @param sample.size Size of ressamples, default is the same size as ind.data 
#' @param parallel if TRUE computations are done in parallel. Some foreach backend must be registered, like doParallel or doMC.
#' @details The row indices of all resamples are drawn in the main process, in the same
#'   order whether or not \code{parallel = TRUE}; only the statistics are computed in parallel.
#'   If \code{ind.data} is a data frame and \code{StatFunc} is \code{cov} or \code{cor} it is converted to
#'   a matrix once; with any other \code{StatFunc} the resamples are data frames.
#'   Resamples whose comparison fails or is NA are redrawn up to 100 times, then an error is raised.
#' @return a data frame with one row per resample: a factor column \code{X1} with the resample number
#'   and the comparison values.
#' @author Diogo Melo, Guilherme Garcia
#' @seealso \code{\link{BootstrapRep}}, \code{\link{AlphaRep}}
#' @rdname BootstrapStat
#' @export
#' @examples
#' cov.matrix <- RandomMatrix(5, 1, 1, 10)
#'
#' BootstrapStat(iris[,1:4], iterations = 50,
#'                ComparisonFunc = function(x, y) PCAsimilarity(x, y)[1],
#'                StatFunc = cov)
#'
#' #Calculating R2 confidence intervals
#' r2.dist <- BootstrapR2(iris[,1:4], 30)
#' quantile(r2.dist)
#' 
#' #Multiple threads can be used with some foreach backend library, like doMC or doParallel
#' #library(doParallel)
#' ##Windows:
#' #cl <- makeCluster(2)
#' #registerDoParallel(cl)
#' ##Mac and Linux:
#' #registerDoParallel(cores = 2)
#' #BootstrapStat(iris[,1:4], iterations = 100,
#' #               ComparisonFunc = function(x, y) KrzCor(x, y)[1],
#' #               StatFunc = cov,
#' #               parallel = TRUE)
#' @keywords parametricsampling
#' @keywords montecarlo
#' @keywords repeatability
BootstrapStat <- function (ind.data, 
                           iterations,
                           ComparisonFunc, 
                           StatFunc, 
                           sample.size = dim (ind.data)[1],
                           parallel = FALSE){
  if(isSymmetric(as.matrix(ind.data))) stop("input appears to be a matrix, use residuals.")
  c.matrix <- StatFunc(ind.data)
  if(any(is.na(c.matrix)))
    stop("The statistic calculated on the observed data has missing values. Check the data for NAs.")
  n_ind <- dim (ind.data)[1]
  data <- .resampling_data(ind.data, StatFunc)
  # all random draws happen here, in the main process and in iteration order,
  # so that parallel and serial runs use the same populations
  indices <- lapply(seq_len(iterations), function(i) sample.int(n_ind, sample.size, TRUE))
  worker <- function(idx)
    .compute_with_retry(idx,
                        function(i) ComparisonFunc(c.matrix, StatFunc(data[i, , drop = FALSE])),
                        function() sample.int(n_ind, sample.size, TRUE))
  .collect_comparisons(.evq_lapply(indices, worker, parallel = parallel), "bootstrap")
}

# Internal helpers shared by BootstrapStat, MonteCarloStat and RarefactionStat -------

# Convert data frames to matrices once, before resampling, when it is known
# that StatFunc (cov or cor) gives the same result for both. Any other
# StatFunc keeps receiving the original object type.
.resampling_data <- function(ind.data, StatFunc) {
  if(is.data.frame(ind.data) && (identical(StatFunc, cov) || identical(StatFunc, cor)))
    as.matrix(ind.data)
  else ind.data
}

# Runs compute(x); if it errors or returns NA, draws a new x with resample()
# and tries again, up to max.retries times, then stops with the last error.
.compute_with_retry <- function(x, compute, resample, max.retries = 100L) {
  retries <- 0L
  repeat {
    problem <- NULL
    out <- tryCatch(compute(x), error = function(e) {
      problem <<- conditionMessage(e)
      NA
    })
    if(!any(is.na(out)))
      return(list(out = out, retries = retries))
    if(is.null(problem)) problem <- "the statistic comparison returned NA"
    if(retries >= max.retries)
      stop("Could not obtain a valid resample after ", max.retries,
           " retries. Last problem: ", problem, call. = FALSE)
    retries <- retries + 1L
    x <- resample()
  }
}

# Turns the list of list(out, retries) into the data frame returned by the
# *Stat functions: a factor column X1 with the iteration number, and one
# column per element of the comparison result (named as the result, or V1, V2, ...).
.collect_comparisons <- function(results, what = "resample") {
  retries <- sum(vapply(results, function(r) r$retries, numeric(1)))
  if(retries > 0)
    warning(retries, " ", what, " samples gave an error or NA in the comparison and were redrawn.",
            call. = FALSE)
  results <- lapply(results, function(r) r$out)
  n <- length(results)
  X1 <- factor(seq_len(n), levels = seq_len(n))
  if(n == 0) return(data.frame(X1 = X1, V1 = numeric(0)))
  k <- length(results[[1]])
  nm <- names(results[[1]])
  simple <- all(vapply(results, function(x)
    is.numeric(x) && is.null(dim(x)) && length(x) == k && identical(names(x), nm), logical(1)))
  if(simple) {
    mat <- matrix(unlist(results, use.names = FALSE), nrow = n, ncol = k, byrow = TRUE)
    if(is.null(nm)) nm <- paste0("V", seq_len(k))
    out <- as.data.frame(mat, stringsAsFactors = FALSE)
    names(out) <- nm
    return(cbind(X1 = X1, out))
  }
  names(results) <- seq_len(n)
  out <- plyr::ldply(results, function(x) x)
  out[[1]] <- X1
  names(out)[1] <- "X1"
  out
}

#' R2 confidence intervals by bootstrap resampling
#'
#' Random populations are generated by  resampling 
#' the suplied data or residuals. R2 is calculated on all the
#' random population's correlation matrices, provinding a distribution based on the original data.
#'
#' @param ind.data Matrix of residuals or indiviual measurments
#' @param iterations Number of resamples to take
#' @param parallel if TRUE computations are done in parallel. Some foreach backend must be registered, like doParallel or doMC.
#' @return returns a vector with the R2 for all populations
#' @author Diogo Melo Guilherme Garcia
#' @seealso \code{\link{BootstrapRep}}, \code{\link{AlphaRep}}
#' @export
#' @examples
#' r2.dist <- BootstrapR2(iris[,1:4], 30)
#' quantile(r2.dist)
#' @keywords bootstrap
#' @keywords integration
#' @keywords repeatability
BootstrapR2 <- function (ind.data, iterations = 1000, parallel = FALSE) {
  if(is.data.frame(ind.data)) ind.data <- as.matrix(ind.data)
  it.r2 <- BootstrapStat(ind.data, iterations,
                          ComparisonFunc = function(x, y) y,
                          StatFunc = function(x) .r2(cor(x)),
                          parallel = parallel)
  return (it.r2[,2])
}
