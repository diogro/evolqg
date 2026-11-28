#' Corrected integration value
#'
#' Calculates the Young correction for integration, using bootstrap resampling
#' Warning: CalcEigenVar is strongly preferred and should probably be used in place of this function..
#'
#' @param ind.data Matrix of individual measurments, or adjusted linear model
#' @param cv.level Coefficient of variation level chosen for integration index adjustment in linear model. Defaults to 0.06.
#' @param iterations Number of resamples to take
#' @param parallel if TRUE computations are done in parallel. Some foreach backend must be registered, like doParallel or doMC.
#' @param ... additional arguments passed to other methods
#' @details For linear models with more than one predictor, individuals are grouped by the
#' combination (interaction) of all the factors in the model. Resamples in which some group has a single
#' individual are drawn again (up to 100 times in a row).
#' @return List with adjusted integration indexes, fitted models and simulated distributions of integration indexes and mean coefficient of variation.
#' @references Young, N. M., Wagner, G. P., and Hallgrimsson, B. (2010).
#' Development and the evolvability of human limbs. Proceedings of the
#' National Academy of Sciences of the United States of America, 107(8),
#' 3400-5. doi:10.1073/pnas.0911856107
#' @author Diogo Melo, Guilherme Garcia
#' @seealso \code{\link{MeanMatrixStatistics}}, \code{\link{CalcR2}}
#' @rdname CalcR2CvCorrected
#' @export
#' @examples
#' integration.dist = CalcR2CvCorrected(iris[,1:4], iterations = 100)
#'
#' #adjusted values
#' integration.dist[[1]]
#'
#' \dontrun{
#' #ploting models
#' library(ggplot2)
#' ggplot(integration.dist$dist, aes(r2, mean_cv)) + geom_point() +
#'        geom_smooth(method = 'lm', color= 'black') + theme_bw()
#'
#' ggplot(integration.dist$dist, aes(eVals_cv, mean_cv)) + geom_point() +
#'        geom_smooth(method = 'lm', color= 'black') + theme_bw()
#'}
#' @keywords correlation
#' @keywords integration

CalcR2CvCorrected  <- function(ind.data, ...) UseMethod("CalcR2CvCorrected")

#' @rdname CalcR2CvCorrected
#' @method CalcR2CvCorrected default
#' @export
CalcR2CvCorrected.default <- function (ind.data, cv.level = 0.06, iterations = 1000, parallel = FALSE, ...) {
  cv <- function (x) return (sd(x)/mean(x))
  Stats = function(x) {
    cov.matrix = var(x)
    cor.matrix = cov2cor(cov.matrix)
    return(c(.r2(cor.matrix), cv(eigen(cov.matrix, only.values = TRUE)$values), mean (apply (x, 2, cv))))
  }
  if(is.data.frame(ind.data)) ind.data <- as.matrix(ind.data)
  it.stats <- BootstrapStat(ind.data, iterations,
                           ComparisonFunc = function(x, y) y,
                           StatFunc = Stats,
                           parallel = parallel)[,-1]
  colnames(it.stats) <- c("r2", "eVals_cv", "mean_cv")
  lm.r2 <- lm(it.stats[,1]~it.stats[,3])
  lm.eVals.cv <- lm(it.stats[,2]~it.stats[,3])
  adjusted.r2 <- lm.r2$coefficients %*% c(1, cv.level)
  adjusted.eVals.cv <- lm.eVals.cv$coefficients %*% c(1, cv.level)
  adjusted.integration  <-  c(adjusted.r2, adjusted.eVals.cv)
  names(adjusted.integration) <- c("r2", "eVals_cv")
  models <- list("r2" = lm.r2, "eVals_cv" = lm.eVals.cv)
  output <- list("adjusted.integration.index" = adjusted.integration, "models" = models, "dist" = it.stats)
  return (output)
}

#' @rdname CalcR2CvCorrected
#' @method CalcR2CvCorrected lm
#' @export
CalcR2CvCorrected.lm <- function (ind.data, cv.level = 0.06, iterations = 1000, parallel = FALSE, ...) {
    cv <- function (x) return (sd(x)/mean(x))
    model <- ind.data
    ind.data <- as.matrix(model$model[[1]])
    fac <- model$model[-1]
    df <- model$df.residual
    res <- residuals (model)
    n <- nrow (res)
    # one group per combination of all the factors in the model
    group <- if (ncol (fac) > 0) interaction (fac, drop = TRUE) else factor (rep (1, n))
    group <- as.integer (group)

    # Weighted mean coefficient of variation within groups. A group with a
    # single individual has no variance, which gives NA.
    mean.cv <- function (current.sample) {
        g <- match (group[current.sample], sort (unique (group[current.sample])))
        x <- ind.data[current.sample, , drop = FALSE]
        group.n <- tabulate (g)
        group.mean <- rowsum (x, g) / group.n
        group.sd <- sqrt (rowsum ((x - group.mean[g, , drop = FALSE])^2, g) / (group.n - 1))
        sum ((rowMeans (group.sd / group.mean) * group.n) / sum (group.n))
    }

    # Resampling happens here, in the main process, in the same order as
    # always; samples with an undefined mean cv are drawn again.
    samples <- vector ("list", iterations)
    mcv <- numeric (iterations)
    for (i in seq_len (iterations)) {
        failures <- 0L
        repeat {
            current.sample <- sample.int (n, n, replace = TRUE)
            mcv[i] <- mean.cv (current.sample)
            if (!is.na (mcv[i])) break
            failures <- failures + 1L
            if (failures >= 100L)
                stop ("Could not draw a bootstrap sample with at least two individuals in every group.",
                      call. = FALSE)
        }
        samples[[i]] <- current.sample
    }
    stats <- .evq_lapply (samples, function (current.sample) {
        sample.res <- res[current.sample, , drop = FALSE]
        corr <- cov2cor (var (sample.res) * (n - 1)/df)
        c (mean (corr[lower.tri (corr)]^2),
           cv (eigen (var (sample.res * (n - 1)/df), only.values = TRUE)$values))
    }, parallel = parallel)
    stats <- matrix (unlist (stats, use.names = FALSE), ncol = 2, byrow = TRUE)
    it.stats <- cbind (r2 = stats[, 1], evar = stats[, 2], mcv = mcv)
    colnames(it.stats) <- c("r2", "eVals_cv", "mean_cv")
    lm.r2 <- lm(it.stats[,1]~it.stats[,3])
    lm.eVals.cv <- lm(it.stats[,2]~it.stats[,3])
    adjusted.r2 <- lm.r2$coefficients %*% c(1, cv.level)
    adjusted.eVals.cv <- lm.eVals.cv$coefficients %*% c(1, cv.level)
    adjusted.integration  <-  c(adjusted.r2, adjusted.eVals.cv)
    names(adjusted.integration) <- c("r2", "eVals_cv")
    models <- list("r2" = lm.r2, "eVals_cv" = lm.eVals.cv)
    output <- list("adjusted.integration.index" = adjusted.integration, "models" = models, "dist" = it.stats)
    return (output)
}
