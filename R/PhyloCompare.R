#'Compares sister groups
#'
#'Calculates the comparison of some statistic between sister groups along a phylogeny
#'
#'@param tree phylogenetic tree
#'@param node.data list of node data
#'@param ComparisonFunc comparison function, default is PCAsimilarity
#'@param ... Additional arguments passed to ComparisonFunc
#'@return list with a data.frame of calculated comparisons for each node, using labels or numbers from tree; and a list of comparisons for plotting using phytools (see examples)
#'@note Phylogeny must be fully resolved (binary); an error is raised otherwise. \code{node.data} must be named
#'by node, as returned by \code{\link{PhyloW}}; nodes are looked up by name, not by position.
#'@author Diogo Melo
#'@export
#'@importFrom ape reorder.phylo
#'@import plyr
#'@examples
#'library(ape)
#'data(bird.orders)
#'tree <- bird.orders
#'mat.list <- RandomMatrix(5, length(tree$tip.label))
#'names(mat.list) <- tree$tip.label
#'sample.sizes <- runif(length(tree$tip.label), 15, 20)
#'phylo.state <- PhyloW(tree, mat.list, sample.sizes)
#'
#'phylo.comparisons <- PhyloCompare(tree, phylo.state)
#'
#'# plotting results on a phylogeny:
#'\dontrun{
#'library(phytools)
#'plotBranchbyTrait(tree, phylo.comparisons[[2]])
#'}
PhyloCompare <- function(tree, node.data, ComparisonFunc = PCAsimilarity, ...){
  if(!ape::is.binary(tree))
    stop("PhyloCompare requires a fully resolved (binary) tree; use ape::multi2di() to resolve polytomies.")
  node.names <- .phylo_node_names(tree)
  missing.nodes <- setdiff(node.names, names(node.data))
  if(length(missing.nodes) > 0)
    stop("node.data is missing nodes: ", paste(missing.nodes, collapse = ", "))
  node.order <- unique(reorder(tree, "postorder")$edge[,1])
  phylo.comparisons <- list()
  tree.comparisons <- list()
  for (node in node.order){
    current.nodes <- tree$edge[which(tree$edge[,1]==node),2]
    current <- ComparisonFunc(node.data[[node.names[current.nodes[1]]]],
                              node.data[[node.names[current.nodes[2]]]], ...)
    phylo.comparisons[[node.names[node]]] <- current
    tree.comparisons[[paste(node, current.nodes[1], sep = ',')]] <- current[1]
    tree.comparisons[[paste(node, current.nodes[2], sep = ',')]] <- current[1]
  }
  tree.order <- paste(tree$edge[,1], tree$edge[,2], sep = ',')
  phylo.comparisons.df <- ldply(phylo.comparisons, .id = 'node')
  return(list(df = phylo.comparisons.df, tree = tree.comparisons[tree.order]))
}
