#' Modularity and integration analysis tool
#' 
#' Combines and compares many modularity hypothesis to a covariance matrix. Comparison values are
#' adjusted to the number of zeros in the hypothesis using a linear regression. Best hypothesis can
#' be assessed using a jack-knife procedure.
#' 
#' @param c.matrix Correlation or covariance matrix
#' @param ind.data Matrix of residuals or individual measurements
#' @param modularity.hypot Matrix of hypothesis. Each line represents a trait and each column a module.
#' if modularity.hypot[i,j] == 1, trait i is in module j.
#' @param significance Logical. Indicates if goodness of fit test should be performed. If TRUE, the output
#' has two extra columns, \code{sig} and \code{p.value}. See Details.
#' @param sample.size sample size in goodness of fit simulations via MonteCarlo 
#' @param iterations number os goodness of fit simulations
#' @param n number of jackknife samples
#' @param leave.out number of individuals to be left out of each jackknife, default is 10\%
#' @param ... currently ignored
#' @param parallel if TRUE the jackknife samples are computed in parallel. Some foreach back-end must be
#' registered, like doParallel.
#' @details When \code{significance = TRUE} two goodness of fit measures are returned for each hypothesis,
#' both based on Monte Carlo simulations of \code{iterations} samples of size \code{sample.size} drawn from the
#' hypothesised covariance matrix:
#' \itemize{
#' \item \code{sig}: the fraction of simulated gammas that are \emph{below} the observed gamma. Values near 1
#' mean the observed matrix is farther from the hypothesis than most simulated ones, that is, a poor fit.
#' \item \code{p.value}: (number of simulated gammas >= observed gamma + 1)/(\code{iterations} + 1). High
#' values mean the data are consistent with the hypothesis, low values reject it. Simulated and observed
#' matrices are both compared on the correlation scale.
#' }
#' If the best hypothesis is the null (no modules, identity matrix) a warning is issued and
#' \code{Modularity_hypothesis} is the identity matrix (each trait in its own module).
#' @note Hypothesis can be named as column names, and these will be used to make labels in the 
#' output. 
#' @return Dataframe with ranked hypothesis, ordered by the corrected gamma value
#' Jackknife will return the best hypothesis for each sample.
#' @export
#' @aliases JackKnifeMINT
#' @references 
#' Marquez, E.J. 2008. A statistical framework for testing modularity in multidimensional data. 
#' Evolution 62:2688-2708. 
#' 
#' Parsons, K.J., Marquez, E.J., Albertson, R.C. 2012. Constraint and opportunity: the genetic 
#' basis and evolution of modularity in the cichlid mandible. The American Naturalist 179:64-78.
#' 
#' http://www-personal.umich.edu/~emarquez/morph/doc/mint_man.pdf
#' @examples
#' # Creating a modular matrix:
#' modules = matrix(c(rep(c(1, 0, 0), each = 5),
#'                  rep(c(0, 1, 0), each = 5),
#'                  rep(c(0, 0, 1), each = 5)), 15)
#'                  
#' cor.hypot = CreateHypotMatrix(modules)[[4]]
#' hypot.mask = matrix(as.logical(cor.hypot), 15, 15)
#' mod.cor = matrix(NA, 15, 15)
#' mod.cor[ hypot.mask] = runif(length(mod.cor[ hypot.mask]), 0.8, 0.9) # within-modules
#' mod.cor[!hypot.mask] = runif(length(mod.cor[!hypot.mask]), 0.1, 0.2) # between-modules
#' diag(mod.cor) = 1
#' mod.cor = (mod.cor + t(mod.cor))/2 # correlation matrices should be symmetric
#' 
#' # True hypothesis and a bunch of random ones.
#' hypothetical.modules = cbind(modules, matrix(sample(c(1, 0), 4*15, replace=TRUE), 15, 4))
#' 
#' # if hypothesis columns are not named they are assigned numbers
#' colnames(hypothetical.modules) <- letters[1:7]
#'  
#' MINT(mod.cor, hypothetical.modules)
#' 
#' random_var = runif(15, 1, 10)
#' mod.data = matrix(rnorm(100*15), 100) %*% chol(sqrt(outer(random_var, random_var)) * mod.cor)
#' out_jack = JackKnifeMINT(mod.data, hypothetical.modules, n = 20)
#' 
#' library(ggplot2)
#' 
#' ggplot(out_jack, aes(rank, corrected.gamma)) + geom_point() + 
#'        geom_errorbar(aes(ymin = lower.corrected, ymax = upper.corrected))
MINT <- function(c.matrix, modularity.hypot, 
                            significance = FALSE, sample.size = NULL, iterations = 1000){
  cor.hypot.list = CombineHypot(modularity.hypot)
  if(is.list(cor.hypot.list)){
    if(is.null(colnames(modularity.hypot))) colnames(modularity.hypot) <- seq_len(dim(modularity.hypot)[2])
    gamma_df <- SingleMINT(cor.hypot.list, c.matrix, 
                               significance, sample.size, iterations)
    gamma_df = gamma_df[order(gamma_df$corrected.gamma),]
    best.mask <- attr(cor.hypot.list, "masks")[[gamma_df$.n[1]]]
    rownames(gamma_df) <- NULL
    if(!any(best.mask)){
      warning("The null hypothesis (no modules) has the best corrected gamma; ",
              "returning an identity hypothesis matrix.")
      best <- diag(nrow(modularity.hypot))
      rownames(best) <- rownames(modularity.hypot)
    } else best <- modularity.hypot[, best.mask]
    return(list(gamma_rank = gamma_df, 
                "Modularity_hypothesis" = best))
  } 
  else .mint_gamma(c.matrix, c.matrix * cor.hypot.list)
}

.mint_gamma <- function(S, S_0) sum((S - S_0)^2)

SingleMINT <- function(cor.hypot.list, c.matrix, 
                           significance = FALSE, sample.size = NULL, iterations = 1000) {
  raw.gamma = vapply(cor.hypot.list, function(x) .mint_gamma(c.matrix, c.matrix * x), numeric(1))
  num.zeros = vapply(cor.hypot.list, function(x) sum(x[lower.tri(x)] == 0), numeric(1))
  scaled.gamma = raw.gamma/max(raw.gamma)
  zeros_lm = lm(scaled.gamma ~ num.zeros)
  corrected.gamma = residuals(zeros_lm)
  gamma_df = data.frame(.id = names(cor.hypot.list), 
                        .n = seq_along(cor.hypot.list),
                        rank = rank(corrected.gamma),
                        corrected.gamma = corrected.gamma,
                        raw.gamma, stringsAsFactors = FALSE)
  if(significance){
    if(is.null(sample.size)) stop('sample.size must be non-null')
    gamma_df$sig = vapply(cor.hypot.list, function(x){
      null_dist = MonteCarloStat(c.matrix * x, sample.size, 
                                 iterations, .mint_gamma, cov)$V1 < .mint_gamma(c.matrix * x, c.matrix)
      sum(null_dist)/iterations
    }, numeric(1))
    obs.cor <- cov2cor(c.matrix)
    gamma_df$p.value = vapply(cor.hypot.list, function(x){
      S_h <- c.matrix * x
      cor.h <- cov2cor(S_h)
      root <- .mvn_root(S_h, "svd")
      null_dist <- vapply(seq_len(iterations), function(i)
        .mint_gamma(cor.h, cor(.rmvnorm(sample.size, root = root))), numeric(1))
      (sum(null_dist >= .mint_gamma(cor.h, obs.cor)) + 1)/(iterations + 1)
    }, numeric(1))
  }
  return(gamma_df)
}

#' @export
#' @rdname MINT
JackKnifeMINT <- function(ind.data, modularity.hypot, 
                                     n = 1000, leave.out = floor(dim(ind.data)[1]/10),
                                     ..., parallel = FALSE){
  if(isSymmetric(as.matrix(ind.data))) stop("input appears to be a matrix, use residuals")
  cor.hypot.list = CombineHypot(modularity.hypot)
  if(is.list(cor.hypot.list)){
    n.ind = dim(ind.data)[1]
    if(is.null(colnames(modularity.hypot))) colnames(modularity.hypot) <- seq_len(dim(modularity.hypot)[2])
    observed = SingleMINT(cor.hypot.list, cov(ind.data))
    jack.list = .evq_lapply(seq_len(n), function(i) {
      current.sample = sample(seq_len(n.ind), n.ind - leave.out)
      as.matrix(SingleMINT(cor.hypot.list, cov(ind.data[current.sample,]))[,4:5])
    }, parallel = parallel)
    out_jack = aperm(simplify2array(jack.list), c(3, 1, 2)) # iteration x hypothesis x (corrected, raw)
    jack_interval = apply(out_jack, 2:3, quantile, c(0.025, 0.975)) # jackknife confidence intervals
    out_order = apply(out_jack[,,1], 1, function(x) rank(x) == rank(observed$corrected.gamma))   # jackkkife ordering
    observed$rank.support = rowSums(out_order)/n
    observed$upper.raw = jack_interval[2,,2]
    observed$upper.corrected = jack_interval[2,,1]
    observed$lower.raw = jack_interval[1,,2]
    observed$lower.corrected = jack_interval[1,,1]
    observed = observed[order(observed$corrected.gamma),]
    rownames(observed) <- NULL
    return(observed)
  }
  else stop("use at least two hypothesis vectors, preferably more")
}
