#' Parametric per trait repeatabilities
#'
#' Estimates the variance in the sample not due to measurement error
#'
#' @note Requires at least two observations per individual. The among-individual
#'   variance is computed with the coefficient n0 of Lessells and Boag (1987), so
#'   unbalanced designs (different numbers of measurements per individual) are handled.
#'   For exactly two measurements per individual n0 is 2.
#' @param ID identity of individuals
#' @param ind.data individual measurements (matrix, data.frame or a vector for a single trait)
#' @return vector of repeatabilities
#' @export
#' @references Lessels, C. M., and Boag, P. T. (1987).
#'    Unrepeatable repeatabilities: a common mistake.
#'    The Auk, 2(January), 116-121.
#' @author Guilherme Garcia
#' @examples
#' num.ind = length(iris[,1])
#' ID = rep(1:num.ind, 2)
#' ind.data = rbind(iris[,1:4], iris[,1:4]+array(rnorm(num.ind*4, 0, 0.1), dim(iris[,1:4])))
#' CalcRepeatability(ID, ind.data)
CalcRepeatability <- function (ID, ind.data)
{
  ind.data <- as.matrix(ind.data)
  ID <- as.factor(ID)
  one.column <- function (vec){
    model <- anova (lm (vec ~ ID))
    msq <- model$'Mean Sq' ## 1 between, 2 within
    # Lessells & Boag (1987) n0: average number of observations per
    # individual, which is 2 for balanced pairs and corrects unbalanced designs
    n.i <- tabulate (droplevels (ID[!is.na (vec)]))
    N <- sum (n.i)
    n0 <- (N - sum (n.i^2)/N)/(length (n.i) - 1)
    s2a <- (msq[1] - msq[2])/n0
    s2a / (s2a + msq[2])
  }
  output <- vapply (seq_len (ncol (ind.data)), function (j) one.column (ind.data[, j]), numeric (1))
  names (output) <- colnames (ind.data)
  return (output)
}
