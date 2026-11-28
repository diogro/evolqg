#' Test single modularity hypothesis using Mantel correlation
#'
#' Calculates the correlation and Mantel significance test between a hypothetical binary modularity
#' matrix and a correlation matrix. Also gives mean correlation within- and between-modules. 
#' This function is usually only called by TestModularity.
#' 
#' CalcAVG can be used when a significance test is not required.
#'
#' @param cor.hypothesis Hypothetical correlation matrix, with 1s within-modules and 0s between modules.
#' @param cor.matrix Observed empirical correlation matrix.
#' @param permutations Number of permutations used in significance calculation.
#' @param MHI Indicates if Modularity Hypothesis Index should be calculated instead of AVG Ratio.
#' @param parallel if TRUE computations are done in parallel. Some foreach back-end must be registered, like doParallel or doMC.
#' @param ... additional arguments passed to MantelCor
#' @param landmark.dim Used if permutations should be performed maintaining landmark structure in geometric morphometric data. Either 2 for 2d data or 3 for 3d data. Default is NULL for non geometric morphometric data.
#' @param withinLandmark Logical. If TRUE within-landmark correlation are used in calculation of correlation. Only used if landmark.dim is passed, default is FALSE.
#' @param perm.matrix Optional matrix of permutations (one per row, as built internally for the list method) used for the Mantel test instead of drawing new ones. Ignored if \code{landmark.dim} is given.
#' @return Returns a vector with the matrix correlation, significance via Mantel, within- and between module correlation.
#' @export
#' @importFrom Matrix bdiag
#' @rdname MantelModTest
#' @references Porto, Arthur, Felipe B. Oliveira, Leila T. Shirai, Valderes Conto, and Gabriel Marroig. 2009. "The Evolution of Modularity in the Mammalian Skull I: Morphological Integration Patterns and Magnitudes." Evolutionary Biology 36 (1): 118-35. doi:10.1007/s11692-008-9038-3.
#' 
#' Modularity and Morphometrics: Error Rates in Hypothesis Testing Guilherme Garcia, Felipe Bandoni de Oliveira, Gabriel Marroig bioRxiv 030874; doi: http://dx.doi.org/10.1101/030874
#' @author Diogo Melo, Guilherme Garcia
#' @seealso \code{\link{MantelCor}},\code{\link{CalcAVG}},\code{\link{TestModularity}}
#' @examples
#' # Create a single modularity hypothesis:
#' hypot = rep(c(1, 0), each = 6)
#' cor.hypot = CreateHypotMatrix(hypot)
#' 
#' # First with an unstructured matrix:
#' un.cor = RandomMatrix(12)
#' MantelModTest(cor.hypot, un.cor)
#' 
#' # Now with a modular matrix:
#' hypot.mask = matrix(as.logical(cor.hypot), 12, 12)
#' mod.cor = matrix(NA, 12, 12)
#' mod.cor[ hypot.mask] = runif(length(mod.cor[ hypot.mask]), 0.8, 0.9) # within-modules
#' mod.cor[!hypot.mask] = runif(length(mod.cor[!hypot.mask]), 0.3, 0.4) # between-modules
#' diag(mod.cor) = 1
#' mod.cor = (mod.cor + t(mod.cor))/2 # correlation matrices should be symmetric
#' 
#' MantelModTest(cor.hypot, mod.cor)
#' @details The list method uses the same set of permutations for every hypothesis (except when \code{landmark.dim} is given, where \code{\link{MantelCor}} is called for each one). Hypotheses with no within-module pairs (all-zero, or modules with a single trait) or with all trait pairs within modules give an undefined correlation, and a warning is issued.
#' @keywords matrixcomparison
#' @keywords matrixcorrelation
#' @keywords manteltest
MantelModTest <- function (cor.hypothesis, cor.matrix, ...) UseMethod("MantelModTest")

#' @rdname MantelModTest
#' @method MantelModTest default
#' @export
MantelModTest.default <- function (cor.hypothesis, cor.matrix,
                                   permutations = 1000, MHI = FALSE, ..., 
                                   landmark.dim = NULL, withinLandmark = FALSE,
                                   perm.matrix = NULL) {
  if(!all(cor.hypothesis == 0 | cor.hypothesis == 1)) stop("modularity hypothesis matrix should be binary")
  pairs <- cor.hypothesis[lower.tri(cor.hypothesis)]
  if(all(pairs == 0))
    warning("Hypothesis has no within-module trait pairs (all zeros or only single-trait modules); the correlation is undefined.")
  else if(all(pairs == 1))
    warning("Hypothesis places all traits in a single module (all ones); the correlation is undefined.")
  if(is.null(landmark.dim)){
    mantel.output <- suppressWarnings(
      .mantel_test(cor.matrix, cor.hypothesis, permutations, perm.matrix))
    mantel.output <- c("rho" = mantel.output$statistic, "Probability" = mantel.output$signif)
  } else{
    mantel.output <- MantelCor(cor.matrix, cor.hypothesis, permutations = permutations,
                               MHI, landmark.dim = landmark.dim, withinLandmark = withinLandmark, 
                               mod = TRUE, ...)
  }
  output = c(mantel.output, 
             CalcAVG(cor.hypothesis, cor.matrix, MHI, landmark.dim))
  return (output)
}

#' @rdname MantelModTest
#' @method MantelModTest list
#' @export
MantelModTest.list <- function (cor.hypothesis, cor.matrix,
                                permutations = 1000, MHI = FALSE, 
                                landmark.dim = NULL, withinLandmark = FALSE,
                                ..., 
                                parallel = FALSE){
  perm.matrix <- NULL
  if(is.null(landmark.dim)) # one set of permutations shared by all hypotheses
    perm.matrix <- .mantel_permutations(nrow(cor.matrix), permutations)
  output <- SingleComparisonMap(cor.hypothesis, cor.matrix,
                                function(x, y) MantelModTest(x, y,
                                                             permutations, MHI, ...,
                                                             landmark.dim = landmark.dim,
                                                             withinLandmark = withinLandmark,
                                                             perm.matrix = perm.matrix),
                                parallel = parallel)
  return(output)
}
