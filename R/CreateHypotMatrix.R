#' Creates binary correlation matrices
#' 
#' Takes a binary vector or column matrix and generates list of binary correlation matrices representing
#' the partition in the vectors. 
#' @param modularity.hypot Matrix of hypothesis. Each line represents a trait and each column a module.
#' if modularity.hypot[i,j] == 1, trait i is in module j.
#' @return binary matrix or list of binary matrices. If a matrix is passed, all the vectors are combined in the 
#' last binary matrix (total hypothesis of full integration hypothesis).
#' @export
#' @examples
#' rand.hypots <- matrix(sample(c(1, 0), 30, replace=TRUE), 10, 3)
#' CreateHypotMatrix(rand.hypots) 
CreateHypotMatrix <- function(modularity.hypot){
  if(is.null(dim(modularity.hypot))) return(outer(modularity.hypot, modularity.hypot))
  num.hyp <- dim (modularity.hypot) [2]
  num.traits <- dim (modularity.hypot) [1]
  m.hyp.list <- lapply(seq_len(num.hyp), function(j) {
    x <- modularity.hypot[, j]
    outer(x, x)
  })
  names(m.hyp.list) <- seq_len(num.hyp) # as plyr::alply did; the last element stays unnamed
  m.hyp.list[[num.hyp+1]] <- matrix(as.integer (as.logical (Reduce ("+", m.hyp.list[1:num.hyp]))),
                                    num.traits, num.traits, byrow=TRUE)
  return(m.hyp.list[1:(num.hyp+1)])
}

# Maximum number of hypotheses CombineHypot accepts (2^n - 1 combinations are built).
.max_combine_hypot <- 20

# Builds the binary correlation matrix for every non-redundant combination of the
# hypothesis columns, plus the null (identity) hypothesis. Combinations that give the
# same matrix are kept once (the first one found, in increasing binary order). The
# column mask of each kept matrix is stored in attr(, "masks").
CombineHypot <- function(modularity.hypot){
  n.hypots = dim(modularity.hypot)[2]  
  if(is.null(n.hypots)) { # if single hypothesis
    cor.hypot = CreateHypotMatrix(modularity.hypot)
    diag(cor.hypot) <- 1
    return(cor.hypot)
  }
  if(n.hypots > .max_combine_hypot)
    stop("Too many hypotheses (", n.hypots, "): all 2^n combinations are tested, so at most ",
         .max_combine_hypot, " are supported.")
  if(n.hypots > 15)
    warning("Combining ", n.hypots, " hypotheses requires building 2^", n.hypots,
            " combinations, which becomes infeasible beyond about 15-20 hypotheses.")
  if(is.null(colnames(modularity.hypot))) colnames(modularity.hypot) <- seq_len(n.hypots)
  n.traits <- dim(modularity.hypot)[1]
  counter = 2^n.hypots - 1
  bits <- 2^((n.hypots - 1):0)
  masks <- lapply(seq_len(counter), function(i) bitwAnd(i, bits) > 0)
  # build the combined matrices and keep the first occurrence of each distinct,
  # non-empty one (in C++, src/fast_utils.cpp)
  combined <- .combine_hypot_cpp(matrix(as.numeric(modularity.hypot), n.traits),
                                 do.call(rbind, masks))
  keep <- combined$keep
  matrices <- combined$matrices
  # single-hypothesis combinations were built with outer(), which kept the
  # trait names as dimnames; keep doing that
  trait.names <- rownames(modularity.hypot)
  if(!is.null(trait.names)) {
    single <- vapply(masks[keep], sum, numeric(1)) == 1
    matrices[single] <- lapply(matrices[single], function(m) {
      dimnames(m) <- list(trait.names, trait.names); m })
  }
  hypot_list <- c(list(null = diag(n.traits)), matrices)
  names(hypot_list) <- c("null", vapply(masks[keep], function(mask)
    paste(colnames(modularity.hypot)[mask], collapse = "_"), character(1)))
  attr(hypot_list, "masks") <- c(list(rep(FALSE, n.hypots)), masks[keep])
  hypot_list
}

#' Create binary hypothesis
#' 
#' Takes a vetor describing a trait partition and returns a binary matrix of the partitions where each line represents a trait and each column a module. In the output matrix, if modularity.hypot[i,j] == 1, trait i is in module j.
#' @param x vector of trait partition. Each partition receive the same symbol.
#' @return Matrix of hypothesis. Each line represents a trait and each column a module.
#' if modularity.hypot[i,j] == 1, trait i is in module j.
#' @export
#' @examples
#' x = sample(c(1, 2, 3), 10, replace = TRUE)
#' Partition2HypotMatrix(x) 
Partition2HypotMatrix <- function(x){
  if(anyNA(x)) stop("x must not contain NA.")
  groups <- unique(x)
  out <- do.call(cbind, lapply(groups, function(i) as.numeric(x == i)))
  if(is.character(groups)) colnames(out) <- groups
  out
}
