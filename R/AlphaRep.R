#' Alpha repeatability
#'
#' Calculates the matrix repeatability using the equation in Cheverud 1996
#' Quantitative genetic analysis of cranial morphology in the cotton-top
#' (Saguinus oedipus) and saddle-back (S. fuscicollis) tamarins. Journal of Evolutionary Biology 9, 5-42.
#'
#'
#' @param cor.matrix Correlation matrix
#' @param sample.size Sample size used in matrix estimation
#' @param method How the expected sampling variance of the correlations is
#'   estimated. One of \code{"cheverud"} (default), \code{"squared"} or
#'   \code{"exact"}. See Details.
#' @details The repeatability is \code{(var(r) - e) / var(r)}, where \code{r} are the
#'   off-diagonal correlations and \code{e} is the variance in \code{r} expected from
#'   sampling error alone with \code{N = sample.size}. The three methods differ in how \code{e} is estimated:
#'   \describe{
#'   \item{\code{"cheverud"}}{\code{(1 - mean(r)^2) / (N - 2)}. The
#'   equation used in Cheverud (1996), and the default so previous results are kept.}
#'   \item{\code{"squared"}}{\code{(1 - mean(r^2)) / (N - 2)}. Averages the
#'   squared correlations instead of squaring the average correlation,
#'   which does not cancel out correlations of different signs.}
#'   \item{\code{"exact"}}{\code{mean((1 - r^2)^2) / (N - 1)}. Averages
#'   the usual large-sample sampling variance of a correlation coefficient,
#'   \code{(1 - r^2)^2 / (N - 1)}, over all correlations.}
#'   }
#' @return Alpha repeatability for correlation matrix
#' @keywords repeatability
#' @seealso \code{\link{MonteCarloStat}}, \code{\link{BootstrapRep}}
#' @author Diogo Melo, Guilherme Garcia
#' @references Cheverud 1996 Quantitative genetic analysis of cranial morphology in the cotton-top
#'   (Saguinus oedipus) and saddle-back (S. fuscicollis) tamarins. Journal of Evolutionary Biology 9, 5-42.
#' @export
#' @examples
#' #For single matrices
#' cor.matrix <- RandomMatrix(10)
#' AlphaRep(cor.matrix, 10)
#' AlphaRep(cor.matrix, 100)
#' AlphaRep(cor.matrix, 100, method = "exact")
#'
#' #For many matrices
#' mat.list <- RandomMatrix(10, 100)
#' sample.sizes <- floor(runif(100, 20, 50))
#' unlist(Map(AlphaRep, mat.list, sample.sizes))

AlphaRep <- function (cor.matrix, sample.size,
                      method = c("cheverud", "squared", "exact")){
  method <- match.arg(method)
  if(!.is_cor_matrix(cor.matrix))
    stop("Matrices do not appear to be correlation matrices.")
  vec <- cor.matrix[lower.tri(cor.matrix)]
  var.erro <- switch(method,
                     cheverud = (1 - mean(vec)^2)/(sample.size - 2),
                     squared  = (1 - mean(vec^2))/(sample.size - 2),
                     exact    = mean((1 - vec^2)^2)/(sample.size - 1))
  var.vec <- var(vec)
  return((var.vec - var.erro)/var.vec)
}
