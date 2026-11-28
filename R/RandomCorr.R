#' Random correlation matrix
#'
#' Internal function for generating random correlation matrices.
#' Use RandomMatrix() instead.
#' @param num.traits Number of traits in random matrix
#' @param ke Parameter for correlation matrix generation. Involves check for positive defitness
#' @return Random Matrix
#' @author Diogo Melo Edgar Zanella
#' @importFrom stats rbeta
#' @keywords randommatrices
RandCorr <- function(num.traits, ke = 10^-3){
  if(num.traits == 1) return(matrix(1, 1, 1))
  # the partial-correlation loop runs in C++ (src/fast_utils.cpp), drawing
  # the same uniforms in the same order as the R loop it replaced
  random.corr = .rand_corr_partial_cpp(as.integer(num.traits), ke)
  perm = sample(seq_len(num.traits))
  random.corr = (random.corr[perm,])[,perm]
  return (random.corr)
}

# Random LKJ(eta) correlation matrix by the vine method of Lewandowski, Kurowicka and
# Joe (2009): partial correlations at level m are Beta(beta, beta) on (-1, 1), with
# beta = eta + (k - 1 - m) / 2, and are converted to raw correlations recursively.
# eta = 1 is uniform over correlation matrices; larger eta concentrates towards the identity.
# Ported from the .rlkj function of the stanimal package (same author).
RandLKJ = function(num.traits, shape){
  # vine-method loop in C++ (src/fast_utils.cpp), same rbeta draws as the R version
  .rand_lkj_cpp(as.integer(num.traits), shape)
}
