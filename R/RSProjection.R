#' Random Skewers projection
#'
#' Uses Bayesian posterior samples of a set of covariance matrices to identify
#' directions of the morphospace in which these matrices differ in their amount of genetic variance.
#'
#' @param cov.matrix.array Array with dimensions traits x traits x populations x MCMCsamples. At least two populations are needed.
#' @param p significance threshold for comparison of variation in each random direction
#' @param num.vectors number of random vectors
#' @param rs_proj output from RSProjection
#' @param ncols number of columns in plot
#' @return projection of all matrices in all random vectors
#' @return set of random vectors and confidence intervals for the projections
#' @return eigen decomposition of the random vectors in directions with significant differences of variations
#' @export
#' @references Aguirre, J. D., E. Hine, K. McGuigan, and M. W. Blows. "Comparing G: multivariate analysis of genetic variation in multiple populations." Heredity 112, no. 1 (2014): 21-29.
#' @examples
#' # small MCMCsample to reduce run time, acctual sample should be larger
#' data(dentus)
#' cov.matrices = dlply(dentus, .(species), function(x) lm(as.matrix(x[,1:4])~1)) |>
#'                laply(function(x) BayesianCalculateMatrix(x, samples = 50)$Ps)
#' cov.matrices = aperm(cov.matrices, c(3, 4, 1, 2))
#' \donttest{
#' rs_proj = RSProjection(cov.matrices, p = 0.8)
#' PlotRSprojection(rs_proj, cov.matrices, ncols = 5)
#' }
#'
RSProjection <- function(cov.matrix.array, p = 0.95, num.vectors = 1000){
  if (dim(cov.matrix.array)[[1]] != dim(cov.matrix.array)[[2]]){
    stop("Covariance matrix array must be of order n x n x m x MCMCsamplele")
  }
  if (is.na(dim(cov.matrix.array)[4])) {
    stop("There are no MCMCsampleles")
  }
  n <- dim(cov.matrix.array)[[1]]
  m <- dim(cov.matrix.array)[[3]]
  MCMCsample <- dim(cov.matrix.array)[[4]]
  if (m < 2) stop("At least two populations are needed for comparison.")

  rand.vec = t(apply(matrix(rnorm(n*num.vectors), n, num.vectors), 2, Normalize))
  rownames(rand.vec) <- seq_len(num.vectors)

  # Projections of all matrices on the random vectors, diag(R M R') without
  # building the num.vectors x num.vectors product
  cov.projection <- apply(cov.matrix.array, 3:4, function(mat) rowSums((rand.vec %*% mat) * rand.vec))
  cov.projection <- aperm(cov.projection, c(3, 2, 1))
  colnames(cov.projection) <- dimnames(cov.matrix.array)[[3]]

  #setting up an index for HPD comparisons
  prs <- cbind(rep(1:m, each = m), 1:m)
  prs.comp <- prs[prs[, 1] < prs[, 2], , drop = FALSE]

  #for a given random vector, examine if the HPDintervals of any pair of  G matrices overlap
  n.pairs <- nrow(prs.comp)
  proj.score <- vapply(seq_len(num.vectors), function(v) {
    HPD <- .hpd_interval(matrix(cov.projection[, , v], MCMCsample, m), prob = p)
    as.numeric(HPD[prs.comp[,1], 1] > HPD[prs.comp[,2],2] |
               HPD[prs.comp[,2], 1] > HPD[prs.comp[,1],2])
  }, numeric(n.pairs))
  proj.score <- matrix(proj.score, num.vectors, n.pairs, byrow = TRUE)

  #collate the random vectors and the outcome of their projection on the G matrices
  vec.score <- cbind(rand.vec, proj.score)
  colnames(vec.score) <- c(1:n,
                           paste0(dimnames(cov.matrix.array)[[3]][prs.comp[, 1]],
                                  ".vs.",
                                  dimnames(cov.matrix.array)[[3]][prs.comp[, 2]]))

  sig.vec <- subset(vec.score, rowSums(vec.score[, (n+1):(n+((m^2 - m)/2)), drop = FALSE]) > 0)

  if(dim(sig.vec)[1] <= 1) {warning("There were <= 1 significant vectors, try a larger num.vectors or lower p");
    eig.R <- NA}
  else{
    #eigen decomposition of the matrix of significant random vectors
    eig.R <- eigen(cov(sig.vec[,1:n]))
    rownames(eig.R$vectors) <- dimnames(cov.matrix.array)[[1]]
    colnames(eig.R$vectors) <- c(paste("e", 1:n, sep = ""))
  }

  out = list(cov.projection = cov.projection, vec.score = vec.score, eig.R = eig.R)
  class(out) = "rs_proj"
  return(out)
}

#' @export
#' @rdname RSProjection
#' @importFrom ggplot2 theme ylab element_text facet_wrap geom_point geom_errorbar .data
PlotRSprojection <- function (rs_proj, cov.matrix.array, p = 0.95, ncols = 5){

  n <- dim(cov.matrix.array)[[1]]
  if (!is.list(rs_proj$eig.R))
    stop("rs_proj has no eigen decomposition (too few significant vectors); try a larger num.vectors or lower p in RSProjection.")

  evecs = t(rs_proj$eig.R$vectors)

  R_eigen_projection <- apply(cov.matrix.array, 3:4, function(mat) rowSums((evecs %*% mat) * evecs))
  R_eigen_projection <- aperm(R_eigen_projection, c(3, 2, 1))
  m <- dim(R_eigen_projection)[2]
  pops <- dimnames(cov.matrix.array)[[3]]
  if (is.null(pops)) pops <- as.character(seq_len(m))

  HPD.int <- lapply(seq_len(n), function(i)
    .hpd_interval(matrix(R_eigen_projection[, , i], ncol = m), prob = p))
  dat <- data.frame(Population = factor(rep(pops, each = n), levels = pops),
                    trait = rep(seq_len(n), times = m),
                    lower = unlist(lapply(seq_len(m), function(j) vapply(HPD.int, function(h) h[j, 1], 0))),
                    upper = unlist(lapply(seq_len(m), function(j) vapply(HPD.int, function(h) h[j, 2], 0))))
  dat$trait = paste0("PC", dat$trait)
  order.list = paste0("PC", 1:n)
  dat$trait = factor(dat$trait, levels = order.list)
  dat$mean = rowMeans(cbind(dat$upper, dat$lower))
  plot = ggplot(dat, aes(.data$Population, .data$mean)) + geom_point() +
    geom_errorbar(aes(ymin = .data$lower, ymax = .data$upper)) + theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ylab('Genetic Variance') + facet_wrap(~ trait, ncol = ncols, scales = "free_y")
  return(plot)
}
