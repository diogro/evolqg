#'Calculates ancestral states of some statistic
#'
#'Calculates weighted average of covariances matrices along a phylogeny, returning a withing-group covariance matrice for each node.
#'
#'@param tree phylogenetic tree
#'@param tip.data list of tip nodes covariance matrices
#'@param tip.sample.size vector of tip nodes sample sizes. If named, names must match the tip labels (order does not matter); if unnamed, it must be in the same order as \code{tip.data}.
#'@return list with calculated within-group matrices, using labels or numbers from tree
#'@note Internal node names are taken from \code{tree$node.label}. If these are missing, empty, repeated
#'(e.g. bootstrap support values) or repeat a tip label, they cannot identify nodes, so a warning is issued
#'and all internal nodes are named by their node numbers instead.
#'@export
#'@importFrom ape reorder.phylo
#'@import plyr
#'@examples
#'library(ape)
#'data(dentus)
#'data(dentus.tree)
#'tree <- dentus.tree
#'mat.list <- dlply(dentus, 'species', function(x) cov(x[,1:4]))
#'sample.sizes <- runif(length(tree$tip.label), 15, 20)
#'PhyloW(tree, mat.list, sample.sizes)

PhyloW <- function(tree, tip.data, tip.sample.size = NULL){
  node.names <- .phylo_node_names(tree)
  if(!all(tree$tip.label %in% names(tip.data))) stop("All tip labels must be in stat list.")
  if(is.null(tip.sample.size)){
    tip.sample.size <- rep(1, length(tip.data))
    names(tip.sample.size) <- names(tip.data)
  } else if(is.null(names(tip.sample.size))){
    if(length(tip.sample.size) != length(tip.data))
      stop("tip.sample.size must have one value per element of tip.data.")
    names(tip.sample.size) <- names(tip.data)
  } else if(!all(tree$tip.label %in% names(tip.sample.size))){
    stop("All tip labels must be in names(tip.sample.size).")
  }
  tip.sample.size <- as.list(tip.sample.size)
  ancestral.stats <- tip.data
  node.order <- unique(reorder(tree, "postorder")$edge[,1])
  for (node in node.order){
      descendants.list <- node.names[tree$edge[which(tree$edge[,1]==node),2]]
      ponderados <- Map('*', ancestral.stats[descendants.list], tip.sample.size[descendants.list])
      tip.sample.size[[node.names[node]]] <- Reduce("+", tip.sample.size[descendants.list])
      ancestral.stats[[node.names[node]]] <- Reduce("+", ponderados)/
                                              tip.sample.size[[node.names[node]]]
  }
  return(ancestral.stats)
}

# Names for all nodes of a tree (tips first, then internal nodes, indexed by node
# number). Internal nodes use tree$node.label when it identifies them uniquely and
# their node numbers otherwise (with a warning when labels exist but are unusable).
.phylo_node_names <- function(tree){
  n.tips <- length(tree$tip.label)
  if(anyDuplicated(tree$tip.label))
    stop("Tip labels must be unique.")
  internal <- as.character(n.tips + seq_len(tree$Nnode))
  if(!is.null(tree$node.label)){
    labels <- as.character(tree$node.label)
    if(any(is.na(labels) | labels == "") || anyDuplicated(labels) ||
       any(labels %in% tree$tip.label)){
      warning(structure(class = c("evolqg_node_labels", "warning", "condition"),
        list(message = paste0("tree$node.label has missing, empty or repeated labels (or labels equal ",
                              "to tip labels); using node numbers to name internal nodes."),
             call = NULL)))
    } else internal <- labels
  }
  c(tree$tip.label, internal)
}
