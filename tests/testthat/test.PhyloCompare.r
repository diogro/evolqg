test_that("PhyloCompare gives correct results", 
{
  TREE <- "iris((setosa:1,versicolor:1):1,virginica:2);"
  tree.iris <- ape::read.tree(text = TREE)
  iris.cov.list <- dlply(iris, "Species", function(x) cov(x[,1:4]))
  cov.matrices <- PhyloW(tree.iris, iris.cov.list)
  RS.iris = PhyloCompare(tree.iris, cov.matrices)
  mantel.iris = PhyloCompare(tree.iris, llply(cov.matrices, cov2cor), function(x, y) cor(x[lower.tri(x)],
                                                                                         y[lower.tri(y)]))
  krz.iris = PhyloCompare(tree.iris, cov.matrices, KrzCor)
  dist.iris = PhyloCompare(tree.iris, cov.matrices, MatrixDistance, "Over")
  
  for (res in list(RS.iris, mantel.iris, krz.iris, dist.iris)) {
    expect_type(res, "list")
    expect_s3_class(res[[1]], "data.frame")
    expect_type(res[[2]], "list")
    internal.check = res[[1]][,2] == llply(unique(res[[1]]$node), function(node) unique(unlist(res[[2]][(grep(paste0(node, ','), names(res[[2]])))])))
    expect_true(all(internal.check))
  }
})

test_that("PhyloCompare compares the children of each node (hand-checked)", {
  tree <- ape::read.tree(text = "((a:1,b:1)n1:1,(c:1,d:1)n2:1)root;")
  mats <- lapply(c(c = 10, a = 1, d = 20, b = 5), function(x) matrix(x)) # not in tip order
  node.data <- PhyloW(tree, mats)
  # n1 = (1 + 5)/2 = 3, n2 = (10 + 20)/2 = 15
  res <- PhyloCompare(tree, node.data, function(x, y) x[1, 1] - y[1, 1])
  df <- res$df
  expect_equal(as.character(df$node), c("n2", "n1", "root")) # postorder
  expect_equal(df$V1, c(10 - 20, 1 - 5, 3 - 15))
  # each edge carries the comparison of its parent node
  expect_equal(names(res$tree), c("5,6", "6,1", "6,2", "5,7", "7,3", "7,4"))
  expect_equal(unlist(res$tree, use.names = FALSE), c(-12, -4, -4, -12, -10, -10))
})

test_that("PhyloCompare works with the same data when node data order changes", {
  tree <- ape::read.tree(text = "((a:1,b:1)n1:1,(c:1,d:1)n2:1)root;")
  mats <- lapply(c(a = 1, b = 5, c = 10, d = 20), function(x) matrix(x))
  nd <- PhyloW(tree, mats)
  f <- function(x, y) x[1, 1] - y[1, 1]
  expect_equal(PhyloCompare(tree, nd, f), PhyloCompare(tree, rev(nd), f))
})

test_that("PhyloCompare requires a binary tree", {
  tree <- ape::read.tree(text = "((a:1,b:1,c:1)n1:1,d:2)root;")
  mats <- lapply(c(a = 1, b = 5, c = 10, d = 20), function(x) matrix(x))
  expect_error(PhyloCompare(tree, PhyloW(tree, mats), function(x, y) 1), "binary")
})

test_that("PhyloCompare works with repeated node labels", {
  tree <- ape::read.tree(text = "((a:1,b:1)100:1,(c:1,d:1)100:1)100;")
  mats <- lapply(c(a = 1, b = 5, c = 10, d = 20), function(x) matrix(x))
  suppressWarnings(nd <- PhyloW(tree, mats))
  expect_warning(res <- PhyloCompare(tree, nd, function(x, y) x[1, 1] - y[1, 1]), "node numbers")
  expect_equal(as.character(res$df$node), c("7", "6", "5"))
  expect_equal(res$df$V1, c(-10, -4, -12))
})
