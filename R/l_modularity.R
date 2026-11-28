#' Defunct functions in evolqg
#'
#' These functions have been removed from evolqg and now signal an error
#' when called.
#'
#' \code{LModularity} was removed in evolqg 1.0. Modularity maximisation
#' is unreliable as a way of finding trait modules
#' (see \url{https://skewed.de/tiago/blog/modularity-harmful}). Use
#' \code{igraph::cluster_optimal} directly if you need it, or test explicit
#' hypotheses with \code{\link{TestModularity}} or \code{\link{MINT}}.
#'
#' @name evolqg-defunct
#' @aliases LModularity
#' @param ... ignored.
#' @return No return value; always signals an error.
#' @seealso \code{\link{TestModularity}}, \code{\link{MINT}}, \code{\link{Defunct}}
#' @keywords internal
NULL

#' @rdname evolqg-defunct
#' @export
LModularity <- function(...){
  .Defunct(msg = paste0("LModularity() was removed in evolqg 1.0 because modularity maximisation is ",
                        "unreliable. Use igraph::cluster_optimal() directly, or test explicit ",
                        "hypotheses with TestModularity() or MINT()."))
}
