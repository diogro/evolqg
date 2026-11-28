#' Rarefaction analysis via resampling
#'
#' Calculates the repeatability of a statistic of the data, such as
#' correlation or covariance matrix, via bootstrap resampling with
#' varying sample sizes, from 3 to the size of the original data.
#'
#' Samples of various sizes, by default without replacement (see \code{replace}), are taken from the full population, a statistic calculated
#' and compared to the full population statistic.
#'
#' A specialized plotting function displays the results in publication quality.
#' @param ind.data Matrix of residuals or individual measurments
#' @param ComparisonFunc comparison function
#' @param ... Additional arguments passed to ComparisonFunc
#' @param num.reps number of populations sampled per sample size
#' @param correlation If TRUE, correlation matrix is used, else covariance matrix. MantelCor always uses correlation matrix.
#' @param replace If true, samples are taken with replacement
#' @param parallel if TRUE computations are done in parallel. Some foreach back-end must be registered, like doParallel or doMC.
#' @details Bootstraping may be misleading with very small sample sizes. Use with caution if original sample sizes are small.
#' @return returns the mean value of comparisons from samples to original statistic, for all sample sizes.
#' @author Diogo Melo, Guilherme Garcia
#' @seealso \code{\link{BootstrapRep}}
#' @rdname Rarefaction
#' @export
#' @examples
#' ind.data <- iris[1:50,1:4]
#'
#' results.RS <- Rarefaction(ind.data, RandomSkewers, num.reps = 5)
#' #Easy parsing of results: mean comparison per sample size
#' sapply(results.RS, mean)
#'
#' \donttest{
#' results.Mantel <- Rarefaction(ind.data, MatrixCor, correlation = TRUE, num.reps = 5)
#' results.KrzCov <- Rarefaction(ind.data, KrzCor, num.reps = 5)
#' results.PCA <- Rarefaction(ind.data, PCAsimilarity, num.reps = 5)
#'}
#'
#' \dontrun{
#' #Multiple threads can be used with some foreach backend library, like doMC or doParallel
#' library(doParallel)
#' registerDoParallel(cores = 2)
#' results.KrzCov <- Rarefaction(ind.data, KrzCor, num.reps = 5, parallel = TRUE)
#'}
#' @keywords rarefaction
#' @keywords bootstap
#' @keywords repeatability
Rarefaction <- function(ind.data,
                        ComparisonFunc,
                        ...,
                        num.reps = 10,
                        correlation = FALSE,
                        replace = FALSE,
                        parallel = FALSE){
  if(correlation)  {StatFunc <- cor; c2v <- cov2cor
  } else {StatFunc <- cov; c2v <- function(x) x}
  rarefaction.list <- RarefactionStat(ind.data,
                                      StatFunc = StatFunc,
                                      ComparisonFunc = function(x, y) ComparisonFunc(c2v(x),
                                                                                     c2v(y), ...),
                                      num.reps = num.reps,
                                      replace = replace,
                                      parallel = parallel)
  return(rarefaction.list)
}

#' Non-Parametric rarefacted population samples and statistic comparison
#'
#' Calculates the repeatability of a statistic of the data, such as
#' correlation or covariance matrix, via resampling with
#' varying sample sizes, from 3 to the size of the original data.
#'
#' Samples of various sizes, by default without replacement (see \code{replace}), are taken from the full population, a statistic calculated and compared to the full population statistic.
#'
#' A specialized ploting function displays the results in publication quality.
#' @param ind.data Matrix of residuals or indiviual measurments
#' @param StatFunc Function for calculating the statistic
#' @param ComparisonFunc comparison function
#' @param ... Aditional arguments passed to ComparisonFunc
#' @param num.reps number of populations sampled per sample size
#' @param replace If true, samples are taken with replacement
#' @param parallel if TRUE computations are done in parallel. Some foreach backend must be registered, like doParallel or doMC.
#' @details Bootstraping may be misleading with very small sample sizes. Use with caution.
#'
#' All samples are drawn in the main process, in the same order whether or not \code{parallel = TRUE}; only the statistics and
#' comparisons are computed in parallel. If \code{StatFunc} raises a warning for a sample, the sample is redrawn, up to 100 times,
#' and then kept with a single summary warning. If \code{StatFunc} already warns on the full data, no samples are redrawn.
#' Data frames are converted to matrices once when \code{StatFunc} is \code{cov} or \code{cor}.
#' @return a list with one numeric vector (of length \code{num.reps}) for each sample size, from 3 to the number of individuals, named by sample size.
#' Each value is the (first element of the) comparison between the statistic of a sample and the statistic of the full data.
#' @author Diogo Melo, Guilherme Garcia
#' @seealso \code{\link{BootstrapRep}}
#' @export
#' @examples
#' ind.data <- iris[1:50,1:4]
#'
#' #Can be used to calculate any statistic via Rarefaction, not just comparisons
#' #Integration, for example:
#' results.R2 <- RarefactionStat(ind.data, cor, function(x, y) CalcEigenVar(y), num.reps = 5)
#'
#' #Easy access
#' sapply(results.R2, mean)
#'
#'\dontrun{
#' #Multiple threads can be used with some foreach backend library, like doMC or doParallel
#' library(doParallel)
#' registerDoParallel(cores = 2)
#' results.R2 <- RarefactionStat(ind.data, cor, function(x, y) CalcEigenVar(y), parallel = TRUE)
#'}
#' @keywords rarefaction
#' @keywords bootstap
#' @keywords repeatability
RarefactionStat <- function(ind.data,
                            StatFunc,
                            ComparisonFunc,
                            ...,
                            num.reps = 10,
                            replace = FALSE,
                            parallel = FALSE)
{
  if(isSymmetric(as.matrix(ind.data))) stop("input appears to be a matrix, use residuals.")
  num.ind <- dim(ind.data)[1]
  if(is.null(num.ind) || num.ind < 3) stop("at least 3 individuals are needed for rarefaction.")
  # a StatFunc that already warns on the full data would make every sample warn
  observed.warnings <- character()
  observed.stat <- withCallingHandlers(StatFunc(ind.data), warning = function(w) {
    observed.warnings <<- c(observed.warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  if(any(is.na(observed.stat)))
    stop("The statistic calculated on the observed data has missing values. Check the data for NAs.")
  max.retries <- 100L
  if(length(observed.warnings) > 0) {
    warning("StatFunc raised a warning on the observed data (", observed.warnings[1],
            "), so samples that raise warnings are not redrawn.", call. = FALSE)
    max.retries <- 0L
  }
  data <- .resampling_data(ind.data, StatFunc)
  sizes <- 3:num.ind
  # all samples are drawn here, in size and then repetition order
  indices <- unlist(lapply(sizes, function(n)
    lapply(seq_len(num.reps), function(r) sample.int(num.ind, n, replace = replace))),
    recursive = FALSE)
  worker <- function(idx) {
    retries <- 0L
    repeat {
      warned <- FALSE
      stat <- withCallingHandlers(StatFunc(data[idx, , drop = FALSE]), warning = function(w) {
        warned <<- TRUE
        invokeRestart("muffleWarning")
      })
      if(!warned || retries >= max.retries) break
      retries <- retries + 1L
      idx <- sample.int(num.ind, length(idx), replace = replace)
    }
    comparison <- ComparisonFunc(observed.stat, stat, ...)
    list(out = unname(comparison[[1]]), kept.warning = warned)
  }
  results <- .evq_lapply(indices, worker, parallel = parallel)
  kept <- sum(vapply(results, function(r) r$kept.warning, logical(1)))
  if(kept > 0 && max.retries > 0)
    warning("StatFunc raised a warning in ", kept, " samples even after ", max.retries,
            " redraws; those samples were kept.", call. = FALSE)
  values <- unlist(lapply(results, function(r) r$out), use.names = FALSE)
  if(is.null(values)) values <- numeric(0)
  split(values, factor(rep(sizes, each = num.reps), levels = sizes))
}

#' Plot Rarefaction analysis
#'
#' A specialized ploting function displays the results from Rarefaction functions in publication quality.
#'
#' @param comparison.list output from rarefaction functions can be used in ploting
#' @param y.axis Y axis lable in plot
#' @param x.axis Y axis lable in plot
#' @return ggplot2 object with rarefaction plot
#' @author Diogo Melo, Guilherme Garcia
#' @seealso \code{\link{BootstrapRep}}
#' @export
#' @importFrom ggplot2 ggplot aes geom_boxplot scale_x_continuous scale_y_continuous theme_bw .data
#'
#' @examples
#'\donttest{
#' ind.data <- iris[1:50,1:4]
#'
#' results.RS <- Rarefaction(ind.data, PCAsimilarity, num.reps = 5)
#' results.Mantel <- Rarefaction(ind.data, MatrixCor, correlation = TRUE, num.reps = 5)
#' results.KrzCov <- Rarefaction(ind.data, KrzCor, num.reps = 5)
#' results.PCA <- Rarefaction(ind.data, PCAsimilarity, num.reps = 5)
#'
#' #Plotting using ggplot2
#' a <- PlotRarefaction(results.RS, "Random Skewers")
#' b <- PlotRarefaction(results.Mantel, "Mantel")
#' c <- PlotRarefaction(results.KrzCov, "KrzCor")
#' d <- PlotRarefaction(results.PCA, "PCAsimilarity")
#'
#' if (requireNamespace("cowplot", quietly = TRUE)) {
#'   cowplot::plot_grid(a, b, c, d, labels = c("RS",
#'                                             "Mantel Correlation",
#'                                             "Krzanowski Correlation",
#'                                             "PCA Similarity"),
#'                      scale = 0.9)
#' }
#'}
#' @keywords rarefaction
#' @keywords bootstap
#' @keywords repeatability
#' @export
PlotRarefaction <- function(comparison.list, y.axis = "Statistic",
                            x.axis = "Number of sampled specimens"){
  labels <- names(comparison.list)
  if(is.null(labels)) labels <- seq_along(comparison.list)
  plot.df <- data.frame(avg.corr = as.numeric(unlist(comparison.list, use.names = FALSE)),
                        L1 = as.numeric(rep(labels, lengths(comparison.list))))
  rarefaction.plot = ggplot(plot.df, aes(.data$L1, .data$avg.corr, group = .data$L1)) +
  geom_boxplot() +
  scale_x_continuous(x.axis) +
  scale_y_continuous(y.axis) +
  theme_bw()
  return(rarefaction.plot)
}
