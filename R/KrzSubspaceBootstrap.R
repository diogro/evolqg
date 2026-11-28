#' Quasi-Bayesian Krzanowski subspace comparison
#'
#' Calculates the usual Krzanowski subspace comparison using a posterior samples
#' for a set of phenotypic covariance matrices. Then, this observed comparison 
#' is contrasted to the subspace comparison across a permutation of the original
#' data. Residuals, which are used to calculate the observed P-matrices, are 
#' permuted across groups: in each replicate the pooled residuals are shuffled
#' once and split into groups with the original sample sizes, so each residual
#' is used exactly once. The null P-matrices are computed with the same degrees
#' of freedom as the observed ones. This process is repeated, creating a null distribution
#' of  subspace comparisons under the hypothesis that all P-matrices come from the
#' same population. This method is a modification on the fully Bayesian method 
#' proposed in Aguirre et. al 2013 and improved in Morrisey et al 2019.
#' 
#' @param x list of linear models from which P-matrices should be calculated
#' @param rep number of bootstrap samples to be made
#' @param MCMCsamples number of MCMCsamples for each P-matrix posterior distribution.
#' @param parallel if TRUE computations are done in parallel. Some foreach backend must be registered, like doParallel or doMC.
#' @param k number of dimensions to be retained in calculating the subspace.
#' Default is 3. It must not be larger than the number of traits.
#' @return A list with the observed and randomized eigenvalue distributions for the posterior Krz Subspace comparisons.
#' @references 
#' Aguirre, J. D., E. Hine, K. McGuigan, and M. W. Blows. 2013. “Comparing G: multivariate analysis of genetic variation in multiple populations.” Heredity 112 (February): 21–29.
#'
#' Morrissey, Michael B., Sandra Hangartner, and Keyne Monro. 2019. “A Note on Simulating Null Distributions for G Matrix Comparisons.” Evolution; International Journal of Organic Evolution 73 (12): 2512–17.
#' @seealso \code{\link{KrzSubspaceDataFrame}}, \code{\link{PlotKrzSubspace}}
#' @export
#' @examples 
#' 
#' \donttest{
#' library(plyr)
#' data(ratones)
#' 
#' model_formula = paste("cbind(", paste(names(ratones)[13:20], collapse = ", "), ") ~ SEX")
#' lm_models = dlply(ratones, .(LIN), function(df) lm(as.formula(model_formula), data = df))
#' # small rep and MCMCsamples to reduce run time, actual values should be larger
#' krz_comparsion = KrzSubspaceBootstrap(lm_models, rep = 5, MCMCsamples = 50)
#' krz_df = KrzSubspaceDataFrame(krz_comparsion)
#' PlotKrzSubspace(krz_df)
#' }
KrzSubspaceBootstrap = function(x, rep = 1, MCMCsamples = 1000, parallel = FALSE, k = 3){
  p <- ncol(residuals(x[[1]]))
  .check_krz_k(k, p)
  Ps_list <- lapply(x, function(m) unclass(BayesianCalculateMatrix(m, samples = MCMCsamples)$Ps))
  MCMC.H.val <- .krz_H_val(Ps_list, k)

  res_list <- lapply(x, residuals)
  n_list <- vapply(res_list, nrow, numeric(1))
  df_list <- vapply(x, function(m) m$df.residual, numeric(1))
  pooled <- do.call(rbind, res_list)
  group <- rep(seq_along(n_list), n_list)

  # One permutation of the pooled residuals per replicate, drawn in the master
  # process so that parallel and serial runs use the same permutations.
  perms <- lapply(seq_len(rep), function(i) split(sample(nrow(pooled)), group))

  randomKrz <- function(perm) {
    random_Ps <- lapply(seq_along(perm), function(g) {
      # BayesianCalculateMatrix only needs the residuals and their df; use the
      # df of the observed model so null and observed P-matrices are comparable
      fake_lm <- structure(list(residuals = pooled[perm[[g]], , drop = FALSE],
                                df.residual = df_list[[g]]), class = "lm")
      unclass(BayesianCalculateMatrix(fake_lm, samples = MCMCsamples)$Ps)
    })
    .krz_H_val(random_Ps, k)
  }

  rand <- .evq_lapply(perms, randomKrz, parallel = parallel)
  MCMC.H.val.random <- do.call(rbind, unname(rand))
  rownames(MCMC.H.val.random) <- seq_len(nrow(MCMC.H.val.random))
  list(observed = MCMC.H.val, random = MCMC.H.val.random)
}

# Ps_list: list (populations) of arrays samples x p x p. Returns a samples x p
# matrix with the variance of the shared-subspace eigenvectors (of the average H)
# in each posterior sample.
.krz_H_val <- function(Ps_list, k) {
  n.samples <- dim(Ps_list[[1]])[1]
  Hs <- lapply(seq_len(n.samples), function(s)
    .krz_H(lapply(Ps_list, function(Ps) Ps[s, , ]), k))
  avgH <- Reduce("+", Hs) / length(Hs)
  avgH.vec <- eigen(avgH, symmetric = TRUE)$vectors
  out <- t(vapply(Hs, function(H) colSums(avgH.vec * (H %*% avgH.vec)), numeric(ncol(avgH.vec))))
  colnames(out) <- seq_len(ncol(out))
  out
}

#' Extract confidence intervals from KrzSubspaceBootstrap
#' 
#' Returns the HPD intervals and their midpoints from the object produced by the KrzSubspaceBootstrap() function. Mainly used for ploting using PlotKrzSubspace. 
#' See example in the KrzSubspaceBootstrap function.
#' 
#' @param x output from KrzSubspaceBootstrap function.
#' @param n number of eigenvalues to use
#' @param prob Posterior probability interval. Default is 95\%.
#' 
#' @return Data frame with the posterior HPD interval (\code{lower}, \code{upper}) for the eigenvalues of the H matrix in the KrzSubspace comparison, for the observed and random (null) distributions. The \code{mean} column is the midpoint of the HPD interval, not the posterior mean.
#' 
#' @seealso \code{\link{KrzSubspaceBootstrap}}, \code{\link{PlotKrzSubspace}}
#' @export
KrzSubspaceDataFrame <- function(x, n = ncol(observed), prob = 0.95){
  observed <- x$observed
  random <- x$random
  HPD.H.val <- cbind(.hpd_interval(observed, prob = prob),
                     .hpd_interval(random, prob = prob))

  dat.observed <- as.data.frame(HPD.H.val[seq_len(n), 1:2, drop = FALSE])
  dat.random <- as.data.frame(HPD.H.val[seq_len(n), 3:4, drop = FALSE])
  dat.observed$class <- 'observed'
  dat.random$class <- 'random'
  dat.observed$PC <- dat.random$PC <- paste('PC', seq_len(n), sep = '')
  dat <- rbind(dat.observed, dat.random)
  dat$mean <-rowMeans(cbind(dat$lower, dat$upper))

  orderlist <- paste('PC', seq_len(n), sep = '')
  dat$PC <- factor(dat$PC, orderlist)

  return(dat)
}

#' Plot KrzSubspace boostrap comparison
#' 
#' Shows the null and observed distribution of eigenvalues from the Krzanowski subspace comparison 
#' @param x output from KrzSubspaceDataFrame() function.
#' @return ggplot2 object with the observed vs. random eigenvalues mean and posterior confidence intervals 
#' @importFrom ggplot2 ggplot geom_point geom_pointrange position_dodge aes xlab ylab .data
#' @export
PlotKrzSubspace <- function(x){
  plot <- ggplot(x, aes(.data$PC, .data$mean, color = .data$class,
                        group = interaction(.data$PC, .data$class))) +
    geom_point(position = position_dodge(width = 0.3)) +
    geom_pointrange(position = position_dodge(width = 0.3),
                    aes(ymin = .data$lower, ymax = .data$upper)) +
    xlab('eigenvectors') + ylab('eigenvalues')
  return(plot)
}
