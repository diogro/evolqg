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
  random.corr = array(0., dim = c(num.traits, num.traits))
  b = array(0., dim = c(num.traits, num.traits))
  b[lower.tri(b, diag = TRUE)] = 1.
  random.corr[2:num.traits, 1] = -1 + 2*round(runif(num.traits-1)*10^8)/10^8
  b[2:num.traits, 1] = random.corr[2:num.traits, 1]
  for (i in 2:num.traits)
    b[i, 2:i] = sqrt(1 - random.corr[i, 1]^2)
  for (i in seq_len(num.traits)[-(1:2)]){
    for (j in 2:(i-1)){
      b1 = b[j, 1:(j-1)]%*%b[i, 1:(j-1)]
      b2 = b[j, j]*b[i, j]
      z = b1 + b2
      y = b1 - b2
      if (b2 < ke){
        random.corr[i, j] = b1
      } else
        random.corr[i, j] = y + (z - y)*round(runif(1)*10^8)/10^8
      cosinv = c((random.corr[i, j] - b1)/b2)
      if (is.finite(cosinv)){
        if (cosinv > 1)
          b[i, (j+1):num.traits] = 0
        else if (cosinv < -1){
          b[i, j] = -b[i, j]
          b[i, (j+1):num.traits] = 0
        } else {
          b[i, j] = b[i, j]*cosinv
          sinTheta = sqrt(1 - cosinv^2)

          b[i, (j+1):num.traits] = b[i, (j+1):num.traits]* c(sinTheta)
        }
      }
    }
  }
  random.corr = random.corr + t(random.corr) + diag(rep(1, num.traits))
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
  k <- num.traits
  if (k == 1) return(matrix(1, 1, 1))
  P <- matrix(0, k, k)
  random.corr <- diag(k)
  for (m in seq_len(k - 1)){
    beta <- shape + (k - 1 - m) / 2
    for (i in (m + 1):k){
      P[m, i] <- 2 * rbeta(1, beta, beta) - 1  # partial correlation on (-1, 1)
      p <- P[m, i]
      if (m > 1)                              # partial -> raw correlation
        for (l in rev(seq_len(m - 1)))
          p <- p * sqrt((1 - P[l, i]^2) * (1 - P[l, m]^2)) + P[l, i] * P[l, m]
      random.corr[m, i] <- random.corr[i, m] <- p
    }
  }
  random.corr
}
