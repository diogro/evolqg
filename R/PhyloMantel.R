#' Mantel test with phylogenetic permutations
#' 
#' Performs a matrix correlation with significance given by a phylogenetic Mantel Test. 
#' Pairs of rows and columns are permuted with probability proportional to their phylogenetic distance.
#' 
#' @param tree phylogenetic tree. Tip labels must match row and column names of the input matrices; matrices are reordered to match \code{tree$tip.label}.
#' @param matrix.1 pair-wise comparison/distance matrix
#' @param matrix.2 pair-wise comparison/distance matrix
#' @param ... additional parameters, currently ignored
#' @param permutations Number of permutations used in significance calculation
#' @param ComparisonFunc comparison function taking two matrices, used for both the observed value and the permutations. Default is the correlation between the lower triangles.
#' @param k determines the influence of the phylogeny. 1 is strong influence, and larger values converge to a traditional mantel test.
#' @return returns a vector with the comparison value and the proportion of permutations with a comparison at least as large as the observed one, computed as (number of permuted values >= observed + 1)/(permutations + 1). 
#' @note Triangular input (e.g. the output of \code{MatrixCor} on a list) is symmetrised before permuting.
#' @export
#' @note This method should only be used when there is no option other than representing data as pair-wise. It suffers from low power, and alternatives should be used when available.
#' @references Harmon, L. J., & Glor, R. E. (2010). Poor statistical performance of the Mantel test in phylogenetic comparative analyses. Evolution, 64(7), 2173-2178.
#' @references Lapointe, F. J., & Garland, Jr, T. (2001). A generalized permutation model for the analysis of cross-species data. Journal of Classification, 18(1), 109-127.
#' @author Diogo Melo, adapted from Harmon & Glor 2010
#' @examples 
#' data(dentus)
#' data(dentus.tree)
#' tree = dentus.tree
#' cor.matrices = dlply(dentus, .(species), function(x) cor(x[-5]))
#' comparisons = MatrixCor(cor.matrices)
#' 
#' sp.means = dlply(dentus, .(species), function(x) colMeans(x[-5]))
#' mh.dist = MultiMahalanobis(means = sp.means, cov.matrix = PhyloW(dentus.tree, cor.matrices)$'6')
#' PhyloMantel(dentus.tree, comparisons, mh.dist, k = 10000, permutations = 100)
#' 
#' #similar to MantelCor for large k:
#' \dontrun{
#' PhyloMantel(dentus.tree, comparisons, mh.dist, k = 10000)
#' MantelCor(comparisons, mh.dist)
#' }
PhyloMantel <- function(tree, matrix.1, matrix.2, ..., 
                        permutations=1000, 
                        ComparisonFunc = function(x, y) cor(x[lower.tri(x)], y[lower.tri(y)]), 
                        k = 1){ 
  matrix.1 <- .align_to_tips(.symmetrize_triangular(matrix.1), tree$tip.label, "matrix.1")
  matrix.2 <- .align_to_tips(.symmetrize_triangular(matrix.2), tree$tip.label, "matrix.2")
  corr <- ComparisonFunc(matrix.1, matrix.2) 
  prob.matrix <- phyloProb(tree, k)
  nullstats <- vapply(seq_len(permutations),
                      function(i) ComparisonFunc(matrix.1, permPhylo(prob.matrix, matrix.2)),
                      numeric(1))
  prob <- (sum(nullstats >= corr) + 1)/(permutations + 1)
  c("rho" = corr, "Probability" = prob)
}

# Reorders rows and columns of a matrix to follow the tip labels
.align_to_tips <- function(m, tips, label = "matrix"){
  m <- as.matrix(m)
  if(is.null(rownames(m)) || is.null(colnames(m))){
    if(nrow(m) != length(tips) || ncol(m) != length(tips))
      stop(label, " must have one row and column per tip.")
    rownames(m) <- colnames(m) <- tips
    return(m)
  }
  if(!all(tips %in% rownames(m)) || !all(tips %in% colnames(m)))
    stop("All tip labels must be in the row and column names of ", label, ".")
  m[tips, tips]
}

# Fills the missing triangle of triangular matrices (e.g. MatrixCor output)
.symmetrize_triangular <- function(m){
  m <- as.matrix(m)
  if(isSymmetric(unname(m))) return(m)
  up <- m[upper.tri(m)]
  low <- m[lower.tri(m)]
  empty <- function(x) all(is.na(x) | x == 0)
  if(empty(up)) m[upper.tri(m)] <- t(m)[upper.tri(m)]
  else if(empty(low)) m[lower.tri(m)] <- t(m)[lower.tri(m)]
  else stop("matrices must be symmetric or triangular.")
  m
}

# Calculates probabilities for phylogenetic permutations; from Lapointe and Garland 2001
#' @importFrom ape cophenetic.phylo
phyloProb<-function(tree, k) {
  phylo_dist <- cophenetic.phylo(tree)
  scaled_phylo_dist <- phylo_dist/max(phylo_dist)
  s <- k-scaled_phylo_dist
  s/rowSums(s)
}

# Permutes species according to phylogentic tree; returns tip names in permuted order.
# p is the probability matrix from phyloProb.
phyloPermute<-function(p) {
  tt<-rownames(p)
  nsp<-length(tt)
  order<-sample(seq_len(nsp), replace=FALSE)
  ttnew<-character(nsp)
  cpm<-p
  for(j in order[-nsp]) {
    cpm<-cpm/rowSums(cpm)
    rr<-which(rownames(cpm)==tt[j])
    pp<-cpm[rr,]
    s2<-sample(names(pp), size=1, replace=TRUE, prob=pp)
    slot<-which(tt==s2)
    rc<-which(colnames(cpm)==s2)
    ttnew[slot]<-tt[j]
    cpm<-cpm[-rr,-rc]
  }
  ttnew[which(ttnew=="")]<-tt[order[nsp]]
  
  ttnew
}	

permPhylo <- function(p, matrix.1) 
{ 
  s <- phyloPermute(p) 
  trans <- match(s, rownames(matrix.1))
  matrix.1[trans,trans] 
}
